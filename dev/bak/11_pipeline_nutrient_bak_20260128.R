library(episcope)
load_config("episcope_grn.yaml")
validate_config()
normalize_config_paths()
motif_init <- init_motif_db(db)
motif_db   <- motif_init$motif_db
tf_list_all <- motif_init$tf_list
tf_list <- sort(unique(tf_list_all))
motif_ids_subset <- NULL

# ──────────────────────────────────────────────────────────────────────────────
# Turn modules ON/OFF
# ──────────────────────────────────────────────────────────────────────────────
do_load_footprints_preprocess    <- TRUE
do_tf_binding_sites_prediction   <- TRUE
do_tf_to_target_genes_prediction <- TRUE
do_build_grn                     <- FALSE
do_differential_grn              <- FALSE
do_topic_benchmark               <- FALSE
do_topic_benchmark_sig           <- FALSE
do_topic_subplot_single          <- FALSE
do_topic_subplot_bulk            <- FALSE
do_tf_centric_subplot_single     <- FALSE
do_tf_centric_subplot_bulk       <- FALSE
rerun_benchmark                  <- FALSE
reuse_edge_topics                <- TRUE
topic_benchmark_tag              <- "per_link"


# Step 0. Load footprint data and preprocess ------------------------------
if (do_load_footprints_preprocess == TRUE) {
  fp_cache_dir <- file.path(base_dir, "cache")
  fp_manifest <- load_footprints(
    root_dir = fp_root_dir,
    db_name = db,
    out_dir = file.path(fp_cache_dir, paste0("fp_", db)),
    motif_ids = motif_ids_subset,
    n_motifs = if (is.null(motif_ids_subset)) 100L else NULL
  )

  # readr::write_csv(fp_manifest, file.path(fp_cache_dir, sprintf("fp_%s_manifest.csv", db)))

  if (db == "HOCOMOCOv13" && !isTRUE(attr(fp_manifest, "from_cache"))) {
    # (Optional, when using HOCOMOCO database)
    fp_manifest <- fp_manifest_trim(fp_manifest) # renames files on disk by default
    # Overwrite every annotation CSV referenced by the manifest:
    summary_tbl <- fp_manifest_trim_annots(fp_manifest, n_workers = 18, verbose = TRUE)

    # Inspect what changed:
    dplyr::count(summary_tbl, status)
    sum(summary_tbl$n_fixed, na.rm = TRUE)
  }

  # options(future.globals.maxSize = 32 * 1024^3)
  # Align the peaks based on the peak similarity
  fp_aligned <- align_footprints(fp_manifest, mid_slop = 10L, round_digits = 1L, score_match_pct = 0.8, cache_dir = fp_cache_dir, cache_tag = db, output_mode = "distinct")

  # Example: motif clustering (PFM -> PPM) using JASPAR-like alignment
  source("R/utils_step1_motif_clustering.R")
  fp_motif_clust <- run_fp_motif_clustering(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, alpha = 0.5, min_peak_support = 50L, sim_method = "jaccard", cores = 30L, weight_by_group_size = TRUE, verbose = TRUE)

  length(unique(fp_aligned$id_map$peak_ID))
  length(unique(fp_aligned$id_map$fp_peak_bak))

  # Load inputs needed to build grn_set.
  sample_metadata <- readxl::read_excel(file.path(base_dir, "sample_metadata.xlsx"), na = "NA")
  strict_metadata <- sample_metadata |> dplyr::filter(!is.na(strict_match_rna))
  lenient_metadata <- sample_metadata |> dplyr::filter(run_grn)
  # strict_metadata_AsPC1 <- strict_metadata |> dplyr::filter(cell == "AsPC1")
  # strict_metadata_HPAFII <- strict_metadata |> dplyr::filter(cell == "HPAFII")
  # strict_metadata_Panc1 <- strict_metadata |> dplyr::filter(cell == "Panc1")

  # sample_metadata$id <- sample_metadata$ID

  # Load RNA data and ATAC data
  atac_data <- readr::read_tsv(file.path(base_dir, "All_ATAC.master_table.narrowpeaks.10mil.txt"))
  atac_out <- load_atac(atac_data, sort_peaks = TRUE)
  atac_score <- atac_out$score
  atac_overlap <- atac_out$overlap

  # rna <- readr::read_csv("inst/extdata/HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")
  # rna <- readr::read_csv("/data/homes/yl814/episcope_test/GSE87218_ATAC/GSE85243_RNA_group_averages_mapped_to_SRR.csv")
  # rna <- readr::read_csv("/data/homes/yl814/episcope_test/nutrient_stress/HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")

  rna <- readr::read_csv(file.path(base_dir, "HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv"))
  rna <- clean_hgnc(rna) # clean "HGNC" column
  # Filter rna expression genes/tfs needs to reach the threshold in at least 1 sample (group size)
  rna <- filter_rna_expr(rna, tf_list, hgnc_col = "HGNC", gene_min = threshold_gene_expr, tf_min = threshold_tf_expr, min_samples = 1L)
  # make two RNA tables with columns already named by ATAC ids ──

  # strict_rna
  smap <- dplyr::transmute(strict_metadata, old = strict_match_rna, new = id)
  strict_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
  nm <- names(strict_rna); nm[match(smap$old, nm)] <- smap$new
  strict_rna <- strict_rna |> `names<-`(nm) |> dplyr::as_tibble()

  # lenient_rna
  smap <- dplyr::transmute(lenient_metadata, old = broad_match_rna, new = id)
  lenient_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
  nm <- names(lenient_rna); nm[match(smap$old, nm)] <- smap$new
  lenient_rna <- lenient_rna |> `names<-`(nm) |> dplyr::as_tibble()

  # Build RNA and ATAC data list (used by per-motif workflows).
  rna_atac_build_args <- list(
    atac_score    = atac_score,
    atac_overlap  = atac_overlap,
    rna           = strict_rna,
    metadata      = strict_metadata,
    tf_list       = tf_list,
    motif_db      = motif_db, # tibble
    label_col     = "strict_match_rna",
    expected_n    = 23
  )

  # Build strict condition matching grn_set.
  grn_set <- build_grn_set(
    fp_score      = fp_aligned$fp_score,
    fp_bound      = fp_aligned$fp_bound,
    fp_annotation = fp_aligned$fp_annotation,
    atac_score    = atac_score,
    atac_overlap  = atac_overlap,
    rna           = strict_rna,
    metadata      = strict_metadata,
    tf_list       = tf_list,
    motif_db      = motif_db,
    label_col     = "strict_match_rna",
    expected_n    = 23
  )


  grn_set <- grn_add_rna_expressed(grn_set, label_col = "strict_match_rna", threshold_gene_expr = threshold_gene_expr)
  grn_set <- grn_add_fp_score_condition(grn_set, label_col = "strict_match_rna")
  grn_set <- grn_add_fp_bound_condition(grn_set, label_col = "strict_match_rna", threshold_fp_score = threshold_fp_score)
  grn_set <- grn_filter_fp_bound_condition(grn_set, min_bound = 1L, use_parallel = TRUE)
  # assert_fp_alignment(grn_set)

  step1_out_dir <- file.path(base_dir, "predict_tf_binding_sites")
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID")
  write_grn_outputs(grn_set, out_dir = step1_out_dir, db = db, qn_base_dir = base_dir)

}


# Step 1. Predict TF binding sites ----------------------------------------
if (do_tf_binding_sites_prediction == TRUE) {
  if (!exists("grn_set") || !is.list(grn_set)) {
    cli::cli_abort("`grn_set` not found. Run Step 0 before Step 1.")
  }

  grn_set <- grn_add_rna_condition(grn_set, label_col = "strict_match_rna")
  grn_set <- grn_add_fp_tfs(grn_set)
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID")
  step1_out_dir <- file.path(base_dir, "predict_tf_binding_sites")

  grn_set <- grn_add_fp_tf_corr(
    grn_set,
    method = "pearson",
    cores = 20L,
    chunk_size = 5000L,
    min_non_na = 5L
  )
  grn_set <- grn_filter_fp_tf_corr(
    grn_set,
    method = "pearson",
    r_thr = threshold_fp_tf_corr_r,
    p_thr = threshold_fp_tf_corr_p,
    set_active = TRUE,
    output_bed = file.path(step1_out_dir, sprintf("fp_predicted_tfbs_%s", db)),
    output_bed_condition = file.path(step1_out_dir, sprintf("fp_predicted_tfbs_%s_by_condition", db)),
    label_col = "strict_match_rna",
    verbose = FALSE
  )
  write_grn_tf_corr_outputs(grn_set, out_dir = step1_out_dir, db = db)
}

# Step 2. Connect TF-occupied enhancers to target genes -------------------
if (do_tf_to_target_genes_prediction == TRUE) {


  # atac_peak to target gene corr

  # GeneHancer based
  gh_std   <- load_genehancer_panc(file.path("inst","extdata","GeneHancer_v5.24_elite_panc.csv")) # GeneHancer_v5.24_elite_panc.csv

  # 30kb/25kb/20kb based
  # gene_annot_ref_hg38 <- episcope::episcope_build_gene_annot("hg38")
  # gene_annot_ref_mm10 <- episcope::episcope_build_gene_annot("mm10")

  # gh_std <- episcope_make_windowed_gh(
  #   peaks      = grn_set$atac_score,
  #   genes      = unique(c(grn_set$rna$HGNC, grn_set$rna$ensembl_gene_id)),
  #   flank_bp   = 30000,
  #   mode       = "TSS",
  #   gene_annot = gene_annot_ref_hg38,
  #   id_col     = "HGNC" # or "ensembl_gene_id"
  # )

  options(future.globals.maxSize = 64 * 1024^3)
  fp_res_full_pearson <- correlate_fp_to_genes(
    grn_set          = grn_set,
    gh_tbl           = gh_std,
    gene_mode        = "both",
    fp_score_tbl     = grn_set$fp_score_condition_qn,
    rna_tbl          = grn_set$rna_condition,
    fdr              = threshold_fp_gene_corr_p,
    r_abs_min        = threshold_fp_gene_corr_abs_r,
    method           = "pearson",
    workers          = 20,
    cache_dir        = file.path(base_dir, "cache", "fp_gene_corr"),
    cache_tag        = "nutrient_genehancer_pearson",
    cache_chunk_size = 5000L,
    cache_verbose    = TRUE
  )

  fp_res_full_spearman <- correlate_fp_to_genes(
    grn_set          = grn_set,
    gh_tbl           = gh_std,
    gene_mode        = "both",
    fp_score_tbl     = grn_set$fp_score_condition_qn,
    rna_tbl          = grn_set$rna_condition,
    fdr              = threshold_fp_gene_corr_p,
    r_abs_min        = threshold_fp_gene_corr_abs_r,
    method           = "spearman",
    workers          = 30,
    cache_dir        = file.path(base_dir, "cache", "fp_gene_corr"),
    cache_tag        = "nutrient_genehancer_spearman",
    cache_chunk_size = 5000L,
    cache_verbose    = TRUE
  )

  step2_out_dir <- file.path(base_dir, "connect_tf_target_genes")
  dir.create(step2_out_dir, recursive = TRUE, showWarnings = FALSE)

  fp_gene_corr_pearson_filt <- filter_fp_gene_corr_by_tf_annotation(
    fp_gene_corr_full = fp_res_full_pearson$fp_gene_corr_full,
    fp_annotation = grn_set$fp_annotation,
    r_thr = threshold_fp_gene_corr_abs_r,
    p_adj_thr = threshold_fp_gene_corr_p
  )

  fp_gene_corr_spearman_filt <- filter_fp_gene_corr_by_tf_annotation(
    fp_gene_corr_full = fp_res_full_spearman$fp_gene_corr_full,
    fp_annotation = grn_set$fp_annotation,
    r_thr = threshold_fp_gene_corr_abs_r,
    p_adj_thr = threshold_fp_gene_corr_p
  )

  fp_links_combined <- combine_fp_gene_corr_methods(
    fp_pearson  = fp_gene_corr_pearson_filt,
    fp_spearman = fp_gene_corr_spearman_filt,
    rna_tbl     = grn_set$rna_condition,
    rna_method  = "pearson",
    rna_cores   = 20
  )

  fp_links_filtered <- filter_links_by_fp_rna_criteria(
    fp_links_combined,
    fp_p_adj_thr = 0.01,
    fp_r_thr = 0.3,
    require_pos_rna = TRUE
  )

  status_res <- build_link_status_matrix(
    links = fp_links_filtered$links,
    fp_bound = grn_set$fp_bound_condition,
    rna_expressed = grn_set$rna_expressed,
    tf_col = "TF",
    gene_col = "gene_key",
    peak_col = "peak_ID",
    out_file = file.path(step2_out_dir, sprintf("tf_gene_link_status_matrix_%s.csv", db)),
    chunk_size = 50000L,
    return_keep = TRUE,
    filter_any = TRUE,
    verbose = TRUE
  )

  fp_links_kept <- fp_links_filtered$links[status_res$keep, , drop = FALSE]
  readr::write_csv(fp_links_kept, file.path(step2_out_dir, sprintf("tf_gene_links_filtered_%s.csv", db)))

  fp_gene_corr_use <- fp_res_full_pearson$fp_gene_corr_full |>
    dplyr::semi_join(
      fp_links_kept |> dplyr::select(peak_ID, gene_key),
      by = c("fp_peak" = "peak_ID", "gene_key" = "gene_key")
    )

  if (is.null(grn_set$fp_variance) || is.null(grn_set$rna_variance)) {
    hv_variance <- precompute_hvf_hvg_variance(grn_set, cores = 20, chunk_size = 50000L)
    grn_set$fp_variance <- hv_variance$fp_variance
    grn_set$rna_variance <- hv_variance$rna_variance
  }

  basal_links_step2 <- make_basal_links(
    fp_gene_corr_kept = fp_gene_corr_use,
    fp_annotation     = grn_set$fp_annotation,
    out_dir           = file.path(step2_out_dir, "basal_links_tmp"),
    prefix            = "step2",
    rna_tbl           = grn_set$rna_condition,
    rna_method        = "pearson",
    rna_cores         = 20,
    fp_variance       = grn_set$fp_variance,
    rna_variance      = grn_set$rna_variance
  )

  grn_set_cond <- prepare_grn_set_for_light_by_condition(
    grn_set,
    label_col = "strict_match_rna"
  )

  light_by_condition(
    ds = grn_set_cond,
    basal_links = basal_links_step2,
    out_dir = step2_out_dir,
    prefix = "step2",
    label_col = "strict_match_rna",
    link_score_threshold = 0,
    fp_score_threshold = 0,
    tf_expr_threshold = 0,
    fp_bound_tbl = grn_set$fp_bound_condition,
    rna_expressed_tbl = grn_set$rna_expressed,
    require_fp_bound = TRUE,
    require_gene_expr = TRUE,
    gene_expr_threshold = 1L,
    filter_active = FALSE,
    use_parallel = TRUE,
    workers = 8,
    fp_variance_tbl = grn_set$fp_variance,
    rna_variance_tbl = grn_set$rna_variance
  )



}

# Step 3. diff networks, and topic analysis -------------------------------
step2_out_dir <- file.path(base_dir, "example_tf_target_genes")

# Step 3/4 from Step 2 per-condition tables (per-cell vs 10_FBS)
step2_specs <- build_cellwise_contrasts_from_index(
  index_csv = file.path(step2_out_dir, "step2_per_condition_index.csv"),
  out_dir = step2_out_dir,
  prefix = "step2",
  ctrl_tag = "10_FBS",
  clean_names = FALSE
)
delta_link <- 1
regulated_genes <- 1.5 # 1.5 2

run_links_deltas_driver(
  specs       = step2_specs,
  clean_names = FALSE,
  parallel    = TRUE,
  restrict_to_active_both = FALSE,
  edge_change_min = delta_link,
  keep_all_cols = TRUE
)

step2_delta_csvs <- list.files(step2_out_dir, "_delta_links.csv", full.names = TRUE)
step2_de_gene_log2_abs_min <- if (regulated_genes == 1.5) 0.585 else if (regulated_genes == 2) 1 else NA_real_

step2_bulk <- episcope::filter_links_deltas_bulk(
  step2_delta_csvs,
  gene_expr_min = threshold_gene_expr,
  tf_expr_min = threshold_tf_expr,
  fp_min = threshold_fp_score,
  link_min = threshold_link_score,
  abs_delta_min = delta_link,
  apply_de_gene = TRUE,
  de_gene_log2_abs_min = step2_de_gene_log2_abs_min,
  enforce_link_expr_sign = TRUE,
  expr_dir_col = "log2FC_gene_expr",
  workers = 20
)

do_step3_topic_analysis <- TRUE
if (isTRUE(do_step3_topic_analysis)) {
  delta_links_csv <- list.files(step2_out_dir, "_delta_links\\.csv$", full.names = TRUE)
  if (!length(delta_links_csv)) cli::cli_abort("No *_delta_links.csv files found in {step2_out_dir}")

  edges_all <- load_delta_links_many(delta_links_csv, keep_original = TRUE)
  topic_root <- file.path(base_dir, "topic_models")

  motif_path <- file.path("inst", "extdata", "genome", "JASPAR2024.txt")
  motif_info <- build_tf_cluster_map_from_motif(motif_path)

  k_grid_default <- c(2:15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100)
  k_single_map <- list(HPAFII = 12L, Panc1 = 13L, AsPC1 = 10L)
  celllines <- c("AsPC1", "HPAFII", "Panc1")

  run_vae_ctf_multivi(
    edges_all = edges_all,
    out_root = topic_root,
    celllines = celllines,
    tf_cluster_map = motif_info$tf_cluster_map,
    tf_exclude = motif_info$tf_exclude,
    k_grid_default = k_grid_default,
    k_single_map = k_single_map
  )

  run_vae_doc_topic_heatmaps(topic_root)
  run_vae_topic_delta_network_plots(topic_root, step2_out_dir = step2_out_dir)
  run_vae_topic_delta_network_pathway(topic_root)
}



do_heatmap_prep <- TRUE
if (isTRUE(do_heatmap_prep)) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Heatmap prep requires data.table.")
  }
  cond_link_files <- list.files(
    step2_out_dir,
    pattern = "^step2_cond-.*_tf_gene_links\\.csv$",
    full.names = TRUE
  )
  if (!length(cond_link_files)) {
    stop("No step2_cond-*_tf_gene_links.csv files found for heatmap prep.")
  }

  heatmap_minmax_scale <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok)) {
      return(rep(NA_real_, length(x)))
    }
    rng <- range(x[ok], na.rm = TRUE)
    if (!is.finite(rng[1]) || rng[1] == rng[2]) {
      out <- rep(0, length(x))
      out[!ok] <- NA_real_
      return(out)
    }
    out <- (x - rng[1]) / (rng[2] - rng[1])
    out[!ok] <- NA_real_
    out
  }

  if (exists("markers_epi", inherits = TRUE)) {
    heatmap_markers_epi <- get("markers_epi", inherits = TRUE)
  } else {
    markers_dir <- file.path(
      base_dir, "plot_lineage_plasticity_related_subnetworks", "markers"
    )
    heatmap_markers_epi <- readr::read_tsv(
      file.path(markers_dir, "epithelial.markers_nutrient_stress_all_lines.txt"),
      show_col_types = FALSE
    ) |>
      dplyr::pull(HGNC)
  }
  if (exists("markers_mes", inherits = TRUE)) {
    heatmap_markers_mes <- get("markers_mes", inherits = TRUE)
  } else {
    markers_dir <- file.path(
      base_dir, "plot_lineage_plasticity_related_subnetworks", "markers"
    )
    heatmap_markers_mes <- readr::read_tsv(
      file.path(markers_dir, "mesenchymal.markers_nutrient_stress_all_lines.txt"),
      show_col_types = FALSE
    ) |>
      dplyr::pull(HGNC)
  }

  heatmap_sig_abs_log2fc <- 1
  if (!exists("step2_delta_csvs", inherits = TRUE)) {
    step2_delta_csvs <- list.files(step2_out_dir, "_delta_links.csv", full.names = TRUE)
  }
  read_sig_genes <- function(p) {
    hdr <- names(data.table::fread(p, nrows = 0L, showProgress = FALSE))
    if (!("gene_key" %in% hdr)) return(NULL)
    log2_col <- if ("log2FC_gene_expr" %in% hdr) {
      "log2FC_gene_expr"
    } else if ("log2fc_gene_expr" %in% hdr) {
      "log2fc_gene_expr"
    } else {
      NA_character_
    }
    if (is.na(log2_col)) return(NULL)
    dt <- data.table::fread(
      p,
      select = c("gene_key", log2_col),
      showProgress = FALSE
    )
    data.table::setnames(dt, log2_col, "log2fc_gene_expr")
    dt <- dt[is.finite(log2fc_gene_expr) & abs(log2fc_gene_expr) > heatmap_sig_abs_log2fc]
    if (!nrow(dt)) return(NULL)
    comp <- sub("_delta_links\\.csv$", "", basename(p))
    cond1 <- sub("_vs_.*$", "", comp)
    cell_line <- sub("_.*$", "", cond1)
    dt[, cell_line := cell_line]
    dt
  }
  sig_gene_tbl <- dplyr::bind_rows(Filter(Negate(is.null), lapply(step2_delta_csvs, read_sig_genes)))
  heatmap_sig_genes_all <- unique(sig_gene_tbl$gene_key)
  heatmap_sig_genes_by_cell <- lapply(
    split(sig_gene_tbl$gene_key, sig_gene_tbl$cell_line),
    unique
  )

  read_one_cond <- function(p) {
    dt <- data.table::fread(
      p,
      select = c("TF", "gene_key", "peak_ID", "fp_score", "tf_expr", "gene_expr"),
      showProgress = FALSE
    )
    cond <- sub("^step2_cond-", "", basename(p))
    cond <- sub("_tf_gene_links\\.csv$", "", cond)
    data.table::setnames(dt, c("TF", "peak_ID"), c("tf", "peak_id"))
    dt[, condition := cond]
    dt[, cell_line := sub("_.*$", "", cond)]
    tibble::as_tibble(dt)
  }

  heatmap_links_raw <- dplyr::bind_rows(lapply(cond_link_files, read_one_cond))

  if (!length(heatmap_sig_genes_all)) {
    gene_expr_tbl <- heatmap_links_raw |>
      dplyr::distinct(.data$cell_line, .data$condition, .data$gene_key, .data$gene_expr)
    ctrl_tbl <- gene_expr_tbl |>
      dplyr::filter(.data$condition == paste0(.data$cell_line, "_10_FBS")) |>
      dplyr::select(.data$cell_line, .data$gene_key, ctrl_expr = .data$gene_expr)
    gene_fc_tbl <- gene_expr_tbl |>
      dplyr::left_join(ctrl_tbl, by = c("cell_line", "gene_key")) |>
      dplyr::mutate(
        log2fc_gene_expr = log2((.data$gene_expr + 1) / (.data$ctrl_expr + 1))
      ) |>
      dplyr::filter(is.finite(.data$log2fc_gene_expr)) |>
      dplyr::filter(abs(.data$log2fc_gene_expr) > heatmap_sig_abs_log2fc)
    heatmap_sig_genes_all <- unique(gene_fc_tbl$gene_key)
    heatmap_sig_genes_by_cell <- lapply(
      split(gene_fc_tbl$gene_key, gene_fc_tbl$cell_line),
      unique
    )
  }

  heatmap_links_nonagg <- heatmap_links_raw |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(
      tf_scaled = heatmap_minmax_scale(.data$tf_expr),
      gene_scaled = heatmap_minmax_scale(.data$gene_expr),
      fp_scaled = heatmap_minmax_scale(.data$fp_score)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      hm_combo_mult = .data$tf_scaled * .data$gene_scaled * .data$fp_scaled,
      hm_combo_sum = .data$tf_scaled + .data$gene_scaled + .data$fp_scaled,
      hm_fp = .data$fp_score,
      link_id = paste(.data$tf, .data$gene_key, .data$peak_id, sep = "|"),
      epi = as.integer(.data$tf %in% heatmap_markers_epi |
        .data$gene_key %in% heatmap_markers_epi),
      mes = as.integer(.data$tf %in% heatmap_markers_mes |
        .data$gene_key %in% heatmap_markers_mes)
    )

  heatmap_links_nonagg_fp_distinct <- heatmap_links_nonagg |>
    dplyr::group_by(.data$condition, .data$peak_id) |>
    dplyr::summarise(
      fp_score = max(.data$hm_fp, na.rm = TRUE),
      epi = max(.data$epi, na.rm = TRUE),
      mes = max(.data$mes, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      fp_score = dplyr::if_else(is.finite(.data$fp_score), .data$fp_score, NA_real_),
      link_id = .data$peak_id
    )

  heatmap_links_agg <- heatmap_links_raw |>
    dplyr::group_by(.data$condition, .data$tf, .data$gene_key) |>
    dplyr::summarise(
      fp_score = sum(.data$fp_score, na.rm = TRUE),
      tf_expr = dplyr::first(.data$tf_expr),
      gene_expr = dplyr::first(.data$gene_expr),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$condition) |>
    dplyr::mutate(
      tf_scaled = heatmap_minmax_scale(.data$tf_expr),
      gene_scaled = heatmap_minmax_scale(.data$gene_expr),
      fp_scaled = heatmap_minmax_scale(.data$fp_score)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      expr_prod = .data$tf_expr * .data$gene_expr,
      expr_sum = .data$tf_expr + .data$gene_expr,
      hm_combo_mult = .data$tf_scaled * .data$gene_scaled * .data$fp_scaled,
      hm_combo_sum = .data$tf_scaled + .data$gene_scaled + .data$fp_scaled,
      hm_fp = .data$fp_score,
      hm_expr_mult = .data$expr_prod,
      hm_expr_sum = .data$expr_sum,
      link_id = paste(.data$tf, .data$gene_key, sep = "|"),
      epi = as.integer(.data$tf %in% heatmap_markers_epi |
        .data$gene_key %in% heatmap_markers_epi),
      mes = as.integer(.data$tf %in% heatmap_markers_mes |
        .data$gene_key %in% heatmap_markers_mes)
    )

  heatmap_nonagg_wide <- list(
    combo_mult = heatmap_links_nonagg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, peak_id, condition, value = hm_combo_mult) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    combo_sum = heatmap_links_nonagg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, peak_id, condition, value = hm_combo_sum) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    fp_only = heatmap_links_nonagg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, peak_id, condition, value = hm_fp) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    fp_only_distinct = heatmap_links_nonagg_fp_distinct |>
      dplyr::select(link_id, epi, mes, peak_id, condition, value = fp_score) |>
      tidyr::pivot_wider(names_from = condition, values_from = value)
  )

  heatmap_agg_wide <- list(
    combo_mult = heatmap_links_agg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, condition, value = hm_combo_mult) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    combo_sum = heatmap_links_agg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, condition, value = hm_combo_sum) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    fp_only = heatmap_links_agg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, condition, value = hm_fp) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    expr_mult = heatmap_links_agg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, condition, value = hm_expr_mult) |>
      tidyr::pivot_wider(names_from = condition, values_from = value),
    expr_sum = heatmap_links_agg |>
      dplyr::select(link_id, epi, mes, tf, gene_key, condition, value = hm_expr_sum) |>
      tidyr::pivot_wider(names_from = condition, values_from = value)
  )
}

do_heatmap_plot <- TRUE
if (isTRUE(do_heatmap_plot)) {
  source("R/utils_ggheat.R")
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Heatmap plotting requires matrixStats.")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("Heatmap plotting requires openxlsx2.")
  }

  plots_dir <- file.path(base_dir, "connect_tf_target_genes", "cluster_heatmap_per_link")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  heatmap_debug <- TRUE
  heatmap_log_path <- file.path(plots_dir, "heatmap_debug.log")
  heatmap_log <- function(...) {
    txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
    message(txt)
    cat(txt, "\n", file = heatmap_log_path, append = TRUE)
  }

  heatmap_hv_sizes <- list(
    list(tag = "hv60000", n = 60000L)
  )
  heatmap_scale_mode <- "row"
  heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100)
  heatmap_breaks <- seq(-3, 3, length.out = length(heatmap_colors) + 1)

  make_col_ann <- function(cols) {
    col_ann <- data.frame(
      CellLine = sub("^(.*?)_.*$", "\\1", cols),
      Stress = sub("^.*?_(.*)$", "\\1", cols),
      row.names = cols,
      check.names = FALSE
    ) |>
      tibble::rownames_to_column("orig") |>
      dplyr::mutate(
        Stress = dplyr::if_else(
          .data$Stress == "10_FBS",
          .data$Stress,
          .data$Stress |>
            stringr::str_remove("_\\d+$") |>
            stringr::str_remove("^([0-9]+(?:\\.[0-9]+)?(?:uM)?)_")
        )
      ) |>
      tibble::column_to_rownames("orig")
    col_ann
  }

  make_row_ann <- function(df) {
    ann <- data.frame(
      Epi = factor(ifelse(df$epi == 1L, "Yes", "No"), levels = c("No", "Yes")),
      Mes = factor(ifelse(df$mes == 1L, "Yes", "No"), levels = c("No", "Yes")),
      row.names = df$link_id,
      check.names = FALSE
    )
    ann
  }

  ann_colors <- list(
    CellLine = c(HPAFII = "#6893c6", Panc1 = "#e74b5b", AsPC1 = "#eca72e"),
    Stress = c(
      `10_FBS` = "grey",
      BCAA = "#3a78af",
      FBS = "#ef812f",
      Glc = "#308a4e",
      Gln.Arg = "#414b8c",
      Lys = "#e89c84",
      Met.Cys = "#ca2f2d",
      Trp = "#916ab6"
    ),
    Epi = c(No = "#d9d9d9", Yes = "#1b9e77"),
    Mes = c(No = "#d9d9d9", Yes = "#d95f02")
  )

  hclust_edges_with_labels <- function(hc) {
    if (is.null(hc) || !inherits(hc, "hclust")) return(NULL)
    n <- length(hc$labels)
    map_node <- function(x) ifelse(x < 0, -x, n + x)
    left <- hc$merge[, 1]
    right <- hc$merge[, 2]
    left_id <- map_node(left)
    right_id <- map_node(right)
    left_label <- rep(NA_character_, length(left))
    right_label <- rep(NA_character_, length(right))
    left_neg <- which(left < 0)
    right_neg <- which(right < 0)
    if (length(left_neg)) {
      left_label[left_neg] <- hc$labels[-left[left_neg]]
    }
    if (length(right_neg)) {
      right_label[right_neg] <- hc$labels[-right[right_neg]]
    }
    left_pos <- which(left > 0)
    right_pos <- which(right > 0)
    if (length(left_pos)) {
      left_label[left_pos] <- paste0("node_", n + left[left_pos])
    }
    if (length(right_pos)) {
      right_label[right_pos] <- paste0("node_", n + right[right_pos])
    }
    list(
      edges = data.frame(
        node_id = n + seq_len(nrow(hc$merge)),
        left_id = left_id,
        right_id = right_id,
        height = hc$height,
        left_label = left_label,
        right_label = right_label,
        left_is_leaf = left < 0,
        right_is_leaf = right < 0,
        stringsAsFactors = FALSE
      ),
      leaves = data.frame(
        leaf_id = seq_len(n),
        leaf_label = hc$labels,
        stringsAsFactors = FALSE
      )
    )
  }

  run_one_heatmap <- function(wide_tbl,
                              meta_cols,
                              label,
                              hv_n = NA_integer_,
                              hv_tag = "hv_all",
                              cut_rows = NA_integer_) {
    cond_cols <- setdiff(names(wide_tbl), meta_cols)
    if (!length(cond_cols)) {
      warning("No condition columns for heatmap: ", label)
      return(invisible(NULL))
    }

    mat <- as.matrix(wide_tbl[, cond_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- wide_tbl$link_id
    row_ann <- make_row_ann(wide_tbl)
    col_ann <- make_col_ann(colnames(mat))
    if (isTRUE(heatmap_debug)) {
      nf <- !is.finite(mat)
      heatmap_log(
        "[heatmap] label=", label, " hv=", hv_tag,
        " rows=", nrow(mat), " cols=", ncol(mat),
        " nonfinite=", sum(nf),
        " rows_nonfinite=", sum(rowSums(nf) > 0),
        " cols_nonfinite=", sum(colSums(nf) > 0)
      )
    }

    keep_idx <- seq_len(nrow(mat))
    if (is.finite(hv_n) && nrow(mat) > hv_n) {
      vars <- matrixStats::rowVars(mat, na.rm = TRUE)
      o <- order(vars, decreasing = TRUE)
      keep_idx <- o[seq_len(hv_n)]
      mat <- mat[keep_idx, , drop = FALSE]
      row_ann <- row_ann[rownames(mat), , drop = FALSE]
      if (isTRUE(heatmap_debug)) {
        nf2 <- !is.finite(mat)
        heatmap_log(
          "[heatmap] label=", label, " hv=", hv_tag,
          " after_hv rows=", nrow(mat), " cols=", ncol(mat),
          " nonfinite=", sum(nf2),
          " rows_nonfinite=", sum(rowSums(nf2) > 0),
          " cols_nonfinite=", sum(colSums(nf2) > 0)
        )
      }
    }

    scale_arg <- heatmap_scale_mode
    if (identical(heatmap_scale_mode, "row")) {
      row_sds <- matrixStats::rowSds(mat, na.rm = TRUE)
      row_means <- rowMeans(mat, na.rm = TRUE)
      const_rows <- !is.finite(row_sds) | row_sds == 0
      row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
      mat <- sweep(mat, 1, row_means, "-")
      mat <- sweep(mat, 1, row_sds, "/")
      mat[!is.finite(mat)] <- 0
      if (any(const_rows) && isTRUE(heatmap_debug)) {
        heatmap_log(
          "[heatmap] label=", label, " hv=", hv_tag,
          " const_rows=", sum(const_rows)
        )
      }
      scale_arg <- "none"
    }

    cluster_rows_arg <- TRUE
    cut_rows_use <- cut_rows
    if (nrow(mat) > 65536L) {
      if (!is.na(cut_rows)) {
        message(
          "Heatmap rows exceed 65536 (", nrow(mat),
          "); skipping cutree run for: ", label, "_", hv_tag
        )
        return(invisible(NULL))
      }
      cluster_rows_arg <- FALSE
      cut_rows_use <- NA_integer_
      message(
        "Heatmap rows exceed 65536 (", nrow(mat),
        "); disabling row clustering for: ", label, "_", hv_tag
      )
    }

    suffix <- paste0("_", hv_tag)
    if (!is.na(cut_rows)) {
      suffix <- paste0(suffix, sprintf("_k%d", cut_rows))
    }
    pdf_out <- file.path(plots_dir, sprintf("heatmap_%s%s.pdf", label, suffix))
    xlsx_out <- file.path(plots_dir, sprintf("heatmap_%s%s.xlsx", label, suffix))

    hm <- tryCatch(
      ggheat(
        mat,
        annotation_col = col_ann,
        annotation_row = row_ann,
        annotation_row_side = "right",
        annotation_colors = ann_colors,
        scale = scale_arg,
        color = heatmap_colors,
        breaks = heatmap_breaks,
        cluster_rows = cluster_rows_arg,
        treeheight_row = if (cluster_rows_arg) 50 else 0,
        clustering_method = "ward.D2",
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        cutree_rows = cut_rows_use,
        border_color = NA,
        show_rownames = FALSE,
        filename = pdf_out,
        silent = TRUE,
        width = 12,
        height = 14
      ),
      error = function(e) {
        if (isTRUE(heatmap_debug)) {
          heatmap_log("[heatmap] ERROR label=", label, " hv=", hv_tag, " msg=", conditionMessage(e))
        }
        stop(e)
      }
    )
    message("Saved heatmap: ", pdf_out)

    if (is.null(hm$matrix_ordered)) {
      warning("ggheat did not return matrix_ordered for: ", label)
      return(invisible(NULL))
    }

    meta_keep <- wide_tbl[match(rownames(hm$matrix_ordered), wide_tbl$link_id), meta_cols, drop = FALSE]
    out_mat <- as.data.frame(hm$matrix_ordered)
    out_df <- cbind(meta_keep, out_mat)

    wb <- openxlsx2::wb_workbook()$
      add_worksheet("Heatmap")$
      add_data("Heatmap", out_df, rowNames = TRUE)$
      freeze_pane("Heatmap", firstRow = TRUE, firstCol = TRUE)

    row_tree <- hclust_edges_with_labels(hm$tree_row)
    if (!is.null(row_tree)) {
      row_order <- data.frame(
        order = seq_along(hm$tree_row$order),
        index = hm$tree_row$order,
        label = hm$tree_row$labels[hm$tree_row$order]
      )
      wb$add_worksheet("row_tree_edges")$
        add_data("row_tree_edges", row_tree$edges, rowNames = FALSE)
      wb$add_worksheet("row_tree_leaves")$
        add_data("row_tree_leaves", row_tree$leaves, rowNames = FALSE)
      wb$add_worksheet("row_order")$
        add_data("row_order", row_order, rowNames = FALSE)
    }

    col_tree <- hclust_edges_with_labels(hm$tree_col)
    if (!is.null(col_tree)) {
      col_order <- data.frame(
        order = seq_along(hm$tree_col$order),
        index = hm$tree_col$order,
        label = hm$tree_col$labels[hm$tree_col$order]
      )
      wb$add_worksheet("col_tree_edges")$
        add_data("col_tree_edges", col_tree$edges, rowNames = FALSE)
      wb$add_worksheet("col_tree_leaves")$
        add_data("col_tree_leaves", col_tree$leaves, rowNames = FALSE)
      wb$add_worksheet("col_order")$
        add_data("col_order", col_order, rowNames = FALSE)
    }

    wb$save(xlsx_out)
    message("Saved heatmap matrix: ", xlsx_out)
  }

  heatmap_variants <- list(
    list(
      label = "agg_combo_mult",
      tbl = heatmap_agg_wide$combo_mult,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key")
    ),
    list(
      label = "agg_combo_sum",
      tbl = heatmap_agg_wide$combo_sum,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key")
    ),
    list(
      label = "agg_fp_only",
      tbl = heatmap_agg_wide$fp_only,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key")
    ),
    list(
      label = "agg_expr_mult",
      tbl = heatmap_agg_wide$expr_mult,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key")
    ),
    list(
      label = "agg_expr_sum",
      tbl = heatmap_agg_wide$expr_sum,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key")
    ),
    list(
      label = "nonagg_combo_mult",
      tbl = heatmap_nonagg_wide$combo_mult,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key", "peak_id")
    ),
    list(
      label = "nonagg_combo_sum",
      tbl = heatmap_nonagg_wide$combo_sum,
      meta_cols = c("link_id", "epi", "mes", "tf", "gene_key", "peak_id")
    ),
    list(
      label = "nonagg_fp_only_distinct",
      tbl = heatmap_nonagg_wide$fp_only_distinct,
      meta_cols = c("link_id", "epi", "mes", "peak_id")
    )
  )

  sig_genes_all <- if (exists("heatmap_sig_genes_all", inherits = TRUE)) {
    heatmap_sig_genes_all
  } else {
    character(0)
  }
  sig_genes_by_cell <- if (exists("heatmap_sig_genes_by_cell", inherits = TRUE)) {
    heatmap_sig_genes_by_cell
  } else {
    list()
  }
  sig_peaks_all <- if (exists("heatmap_links_nonagg", inherits = TRUE) && length(sig_genes_all)) {
    unique(heatmap_links_nonagg$peak_id[heatmap_links_nonagg$gene_key %in% sig_genes_all])
  } else {
    character(0)
  }
  sig_peaks_by_cell <- if (exists("heatmap_links_nonagg", inherits = TRUE) && length(sig_genes_by_cell)) {
    stats::setNames(
      lapply(names(sig_genes_by_cell), function(cl) {
        genes <- unique(sig_genes_by_cell[[cl]])
        if (!length(genes)) return(character(0))
        unique(heatmap_links_nonagg$peak_id[
          heatmap_links_nonagg$cell_line == cl &
            heatmap_links_nonagg$gene_key %in% genes
        ])
      }),
      names(sig_genes_by_cell)
    )
  } else {
    list()
  }

  filter_sig_wide <- function(tbl, sig_genes, sig_peaks = NULL) {
    if (!nrow(tbl)) return(tbl)
    if ("gene_key" %in% names(tbl) && length(sig_genes)) {
      tbl <- tbl[tbl$gene_key %in% sig_genes, , drop = FALSE]
    }
    if ("gene_key" %in% names(tbl) && !length(sig_genes)) {
      return(tbl[0, , drop = FALSE])
    }
    if (!is.null(sig_peaks) && "peak_id" %in% names(tbl) && length(sig_peaks)) {
      tbl <- tbl[tbl$peak_id %in% sig_peaks, , drop = FALSE]
    }
    if (!is.null(sig_peaks) && "peak_id" %in% names(tbl) && !length(sig_peaks)) {
      return(tbl[0, , drop = FALSE])
    }
    tbl
  }

  subset_wide_by_cell <- function(tbl, meta_cols, cell_line) {
    data_cols <- setdiff(names(tbl), meta_cols)
    keep_cols <- data_cols[grepl(paste0("^", cell_line, "_"), data_cols)]
    if (!length(keep_cols)) return(NULL)
    tbl[, c(meta_cols, keep_cols), drop = FALSE]
  }

  base_meta_cols <- c("link_id", "epi", "mes", "tf", "gene_key", "peak_id")
  cond_cols <- setdiff(names(heatmap_agg_wide$combo_mult), base_meta_cols)
  cell_lines <- sort(unique(sub("_.*$", "", cond_cols)))

  for (v in heatmap_variants) {
    if (is.null(v$tbl) || !nrow(v$tbl)) next
    tbl_all <- filter_sig_wide(v$tbl, sig_genes_all, sig_peaks_all)
    if (!nrow(tbl_all)) next
    for (hv in heatmap_hv_sizes) {
      run_one_heatmap(
        tbl_all,
        v$meta_cols,
        v$label,
        hv_n = hv$n,
        hv_tag = hv$tag,
        cut_rows = NA_integer_
      )
    }

    for (cell_line in cell_lines) {
      sig_genes_cell <- sig_genes_by_cell[[cell_line]]
      if (is.null(sig_genes_cell) || !length(sig_genes_cell)) next
      sig_peaks_cell <- sig_peaks_by_cell[[cell_line]]
      tbl_cell <- subset_wide_by_cell(v$tbl, v$meta_cols, cell_line)
      if (is.null(tbl_cell) || !nrow(tbl_cell)) next
      tbl_cell <- filter_sig_wide(tbl_cell, sig_genes_cell, sig_peaks_cell)
      if (!nrow(tbl_cell)) next
      for (hv in heatmap_hv_sizes) {
        run_one_heatmap(
          tbl_cell,
          v$meta_cols,
          paste0(v$label, "_", cell_line),
          hv_n = hv$n,
          hv_tag = hv$tag,
          cut_rows = NA_integer_
        )
      }
    }
  }
}

do_pathway_signatures <- TRUE
if (isTRUE(do_pathway_signatures)) {
  # pathway_dir is configured near the top of this script
  pathway_files <- list.files(
    pathway_dir,
    pattern = "significantly_enriched_pathways\\.txt$",
    full.names = TRUE
  )
  parse_overlap <- function(x) {
    x <- as.character(x)
    if (!length(x)) return(list(n = NA_integer_, total = NA_integer_, raw = NA_character_))
    pieces <- strsplit(x, "/", fixed = TRUE)[[1]]
    if (length(pieces) == 2L) {
      return(list(
        n = suppressWarnings(as.integer(pieces[[1]])),
        total = suppressWarnings(as.integer(pieces[[2]])),
        raw = x
      ))
    }
    list(n = NA_integer_, total = NA_integer_, raw = x)
  }
  parse_file <- function(f) {
    base <- basename(f)
    m <- regexec(
      "^nutrient_DE_by_stress_(.+)_strict_signatures_(Up|Down)\\.significantly_enriched_pathways\\.txt$",
      base
    )
    hit <- regmatches(base, m)[[1]]
    if (!length(hit)) return(NULL)
    comp <- hit[[2]]
    direction <- hit[[3]]
    tbl <- readr::read_tsv(f, show_col_types = FALSE)
    if (!nrow(tbl)) return(list(comp = comp, direction = direction, pathways = list()))
    tbl <- tbl |>
      dplyr::mutate(
        Adjusted.P.value = as.numeric(.data$Adjusted.P.value)
      ) |>
      dplyr::filter(is.finite(.data$Adjusted.P.value), .data$Adjusted.P.value < 0.01)
    if (!nrow(tbl)) return(list(comp = comp, direction = direction, pathways = list()))
    pathways <- lapply(seq_len(nrow(tbl)), function(i) {
      row <- tbl[i, , drop = FALSE]
      ov <- parse_overlap(row$Overlap)
      genes <- strsplit(as.character(row$Genes), ";", fixed = TRUE)[[1]]
      list(
        genes = genes[genes != ""],
        adj_p = as.numeric(row$Adjusted.P.value),
        overlap = ov$raw,
        overlap_n = ov$n,
        overlap_total = ov$total,
        db = as.character(row$datatreatment_anchor)
      )
    })
    names(pathways) <- tbl$Term
    list(comp = comp, direction = direction, pathways = pathways)
  }
  pathway_signatures <- list()
  for (f in pathway_files) {
    parsed <- parse_file(f)
    if (is.null(parsed)) next
    if (is.null(pathway_signatures[[parsed$comp]])) {
      pathway_signatures[[parsed$comp]] <- list()
    }
    pathway_signatures[[parsed$comp]][[parsed$direction]] <- parsed$pathways
  }

  if (length(pathway_signatures)) {
    count_dir <- function(x, dir) if (!is.null(x[[dir]])) length(x[[dir]]) else 0L
    pathway_signatures_summary <- data.frame(
      comparison = names(pathway_signatures),
      n_up = vapply(pathway_signatures, count_dir, integer(1), dir = "Up"),
      n_down = vapply(pathway_signatures, count_dir, integer(1), dir = "Down"),
      stringsAsFactors = FALSE
    )
    print(pathway_signatures_summary)
    example_comp <- names(pathway_signatures)[1]
    example_dirs <- names(pathway_signatures[[example_comp]])
    example_dir <- if (length(example_dirs)) example_dirs[1] else NA_character_
    if (!is.na(example_dir)) {
      message(
        "[pathway_signatures] example: ",
        example_comp,
        " ",
        example_dir,
        " n=",
        length(pathway_signatures[[example_comp]][[example_dir]])
      )
    }
  } else {
    message("[pathway_signatures] no pathway files matched filters.")
  }
}

do_topic_benchmark <- TRUE
# Topic benchmarking (delta/fc + per-condition)
if (isTRUE(do_topic_benchmark)) {
  source("R/utils_grn_lda_nmf.R")
  if (!requireNamespace("data.table", quietly = TRUE)) { stop("Topic benchmarking requires data.table.") }
  topic_runs <- list(
    list(tag = "per_link_aggregated", aggregate_peaks = TRUE),
    list(tag = "per_link_non_aggregated", aggregate_peaks = FALSE)
  )

  markers_dir <- file.path(base_dir, "plot_lineage_plasticity_related_subnetworks", "markers")
  markers_epi <- readr::read_tsv(
    file.path(markers_dir, "epithelial.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  markers_mes <- readr::read_tsv(
    file.path(markers_dir, "mesenchymal.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  if (!length(markers_epi) || !length(markers_mes)) {
    stop("markers_epi/markers_mes are empty; check marker files.")
  }

  marker_mode <- "union" # or "exclusive", then exclusive_method and exclusive_top_n will have effect
  exclusive_method <- "mean_score"
  exclusive_top_n <- 1L
  compute_gene_only <- FALSE
  merge_overlap_thresh <- 0.6
  combo_fp_weight <- 1
  combo_expr_weight <- 1
  scale_quantile <- 0.9
  scale_log1p_if_skewed <- TRUE
  scale_skew_ratio <- 10
  neg_weight_mode <- "minmax" # "minmax" or "drop"

  edges_all_tidy <- load_delta_links_all_tidy(step2_delta_csvs)
  filtered_only <- step2_bulk$filtered_paths[grepl("_delta_links_filtered\\.csv$", step2_bulk$filtered_paths)]
  union_keys <- read_union_edge_keys_from_filtered(filtered_only)
  edges_filtered_tidy <- filter_edges_all_tidy_by_union(edges_all_tidy, union_keys)
  edges_filtered_tidy <- edges_filtered_tidy |>
    dplyr::mutate(
      fp_gene_case = dplyr::coalesce(fp_bound_case, 0L) > 0L &
        dplyr::coalesce(gene_expr_flag_case, 0L) >= 1L,
      fp_gene_ctrl = dplyr::coalesce(fp_bound_ctrl, 0L) > 0L &
        dplyr::coalesce(gene_expr_flag_ctrl, 0L) >= 1L,
      tf_expr_either = dplyr::coalesce(tf_expr_flag_case, 0L) >= 1L |
        dplyr::coalesce(tf_expr_flag_ctrl, 0L) >= 1L,
      pass_topic_gate = tf_expr_either & (fp_gene_case | fp_gene_ctrl)
    ) |>
    dplyr::filter(pass_topic_gate)

  log2fc_eps <- 1e-6

  quantile_scale <- function(x,
                             q = scale_quantile,
                             log1p_if_skewed = scale_log1p_if_skewed,
                             skew_ratio = scale_skew_ratio) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    if (isTRUE(log1p_if_skewed)) {
      x_pos <- x[is.finite(x) & x >= 0]
      if (length(x_pos)) {
        med <- stats::median(x_pos, na.rm = TRUE)
        qhi <- stats::quantile(x_pos, probs = q, na.rm = TRUE, names = FALSE)
        if (is.finite(med) && med > 0 && is.finite(qhi) && qhi / med > skew_ratio) {
          x <- log1p(pmax(x, 0))
        }
      }
    }
    s <- stats::quantile(x, probs = q, na.rm = TRUE, names = FALSE)
    s <- if (is.finite(s) && s > 0) s else 1
    pmin(1, x / s)
  }

  minmax_scale <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) {
      return(rep(0, length(x)))
    }
    (x - rng[1]) / diff(rng)
  }

  apply_neg_weight_mode <- function(tbl, col) {
    if (!nrow(tbl)) return(tbl)
    if (neg_weight_mode == "drop") {
      tbl <- tbl |>
        dplyr::filter(is.finite(.data[[col]]), .data[[col]] > 0)
    } else if (neg_weight_mode == "minmax") {
      tbl[[col]] <- minmax_scale(tbl[[col]])
      tbl <- tbl |> dplyr::filter(is.finite(.data[[col]]))
    } else {
      stop("Unknown neg_weight_mode: ", neg_weight_mode)
    }
    tbl
  }

  cond_edges <- NULL
  cond_link_files <- list.files(
    step2_out_dir,
    pattern = "^step2_cond-.*_tf_gene_links\\.csv$",
    full.names = TRUE
  )
  if (length(cond_link_files)) {
    read_one_cond <- function(p) {
      dt <- data.table::fread(
        p,
        select = c("TF", "gene_key", "peak_ID", "link_score", "fp_score", "tf_expr", "gene_expr", "active_link"),
        showProgress = FALSE
      )
      cond <- sub("^step2_cond-", "", basename(p))
      cond <- sub("_tf_gene_links\\.csv$", "", cond)
      data.table::setnames(dt, c("TF", "peak_ID"), c("tf", "peak_id"))
      dt[, comparison_id := cond]
      dt[, delta_link_score := link_score]
      dt[, delta_fp_bed_score := fp_score]
      dt[, delta_tf_expr := tf_expr]
      dt[, delta_gene_expr := gene_expr]
      dt <- dt[active_link == TRUE]
      dt[, .(comparison_id, tf, gene_key, peak_id,
             delta_link_score, delta_fp_bed_score, delta_tf_expr, delta_gene_expr)]
    }
    cond_edges <- data.table::rbindlist(lapply(cond_link_files, read_one_cond), use.names = TRUE, fill = TRUE)
    cond_edges <- tibble::as_tibble(cond_edges)
  }

  for (run_cfg in topic_runs) {
    topic_benchmark_dir <- file.path(step2_out_dir, sprintf("topic_benchmark_%s", run_cfg$tag))
    doc_term_dir <- file.path(topic_benchmark_dir, "doc_terms")
    fit_dir <- file.path(topic_benchmark_dir, "fits")
    dir.create(topic_benchmark_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(doc_term_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

    aggregate_peaks <- isTRUE(run_cfg$aggregate_peaks)

    cond_edges_use <- cond_edges
    if (aggregate_peaks && !is.null(cond_edges) && nrow(cond_edges)) {
      cond_edges_use <- aggregate_edges_by_tf_gene(
        cond_edges,
        sum_cols = c("delta_fp_bed_score"),
        first_cols = c("delta_tf_expr", "delta_gene_expr")
      )
    }

    edges_filtered_use <- edges_filtered_tidy
    if (aggregate_peaks) {
      edges_filtered_use <- aggregate_edges_by_tf_gene(
        edges_filtered_tidy,
        sum_cols = c("delta_fp_bed_score", "fp_bed_score_case", "fp_bed_score_ctrl", "delta_link_score"),
        first_cols = c("delta_tf_expr", "delta_gene_expr",
                       "tf_expr_case", "tf_expr_ctrl",
                       "gene_expr_case", "gene_expr_ctrl",
                       "log2FC_tf_expr", "log2FC_gene_expr")
      )
    }

    edges_fc_log2_use <- compute_log2fc_from_case_ctrl(edges_filtered_use, eps = log2fc_eps)

    build_docs_for_label <- function(label) {
    if (label == "cond_combo") {
      if (is.null(cond_edges_use) || !nrow(cond_edges_use)) return(NULL)
      dt <- cond_edges_use |>
        dplyr::mutate(
          expr_prod = delta_tf_expr * delta_gene_expr,
          expr_scaled = quantile_scale(expr_prod),
          fp_scaled = quantile_scale(delta_fp_bed_score),
          delta_link_score = combo_fp_weight * fp_scaled + combo_expr_weight * expr_scaled
        )
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = 0,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "cond_fp") {
      if (is.null(cond_edges_use) || !nrow(cond_edges_use)) return(NULL)
      return(build_pooled_doc_term(
        cond_edges_use |> dplyr::mutate(delta_link_score = delta_fp_bed_score),
        which = "abs",
        min_abs_delta = 0,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "cond_expr") {
      if (is.null(cond_edges_use) || !nrow(cond_edges_use)) return(NULL)
      return(build_pooled_doc_term(
        cond_edges_use |>
          dplyr::mutate(delta_link_score = delta_tf_expr * delta_gene_expr),
        which = "abs",
        min_abs_delta = 0,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "diff_delta_combo") {
      dt <- edges_filtered_use |>
        dplyr::mutate(
          expr_prod = delta_tf_expr * delta_gene_expr,
          expr_scaled = quantile_scale(abs(expr_prod)),
          fp_scaled = quantile_scale(abs(delta_fp_bed_score)),
          delta_link_score = combo_fp_weight * fp_scaled + combo_expr_weight * expr_scaled
        )
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = delta_link,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "diff_delta_fp") {
      dt <- edges_filtered_use |>
        dplyr::mutate(delta_link_score = delta_fp_bed_score * sign(delta_gene_expr))
      dt <- apply_neg_weight_mode(dt, "delta_link_score")
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = if (neg_weight_mode == "minmax") 0 else delta_link,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "diff_delta_expr") {
      dt <- edges_filtered_use |>
        dplyr::mutate(delta_link_score = delta_tf_expr * delta_gene_expr)
      dt <- apply_neg_weight_mode(dt, "delta_link_score")
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = if (neg_weight_mode == "minmax") 0 else delta_link,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "diff_fc_combo") {
      dt <- edges_fc_log2_use |>
        dplyr::mutate(
          expr_prod = log2fc_tf_expr * log2fc_gene_expr,
          fp_signed = log2fc_fp_bed_score * sign(log2fc_gene_expr),
          delta_link_score = expr_prod + fp_signed
        )
      dt <- apply_neg_weight_mode(dt, "delta_link_score")
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = 0,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "diff_fc_fp") {
      dt <- edges_fc_log2_use |>
        dplyr::mutate(delta_link_score = log2fc_fp_bed_score * sign(log2fc_gene_expr))
      dt <- apply_neg_weight_mode(dt, "delta_link_score")
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = 0,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    if (label == "diff_fc_expr") {
      dt <- edges_fc_log2_use |>
        dplyr::mutate(delta_link_score = log2fc_tf_expr * log2fc_gene_expr)
      dt <- apply_neg_weight_mode(dt, "delta_link_score")
      return(build_pooled_doc_term(
        dt,
        which = "abs",
        min_abs_delta = 0,
        top_terms_per_doc = 1000,
        min_df = 2
      ))
    }
    NULL
  }

  edges_for_label <- function(label) {
    if (grepl("^cond_", label)) return(cond_edges_use)
    if (grepl("^diff_fc_", label)) return(edges_fc_log2_use)
    edges_filtered_use
  }

  variants <- list(
    list(label = "cond_combo"),
    list(label = "cond_fp"),
    list(label = "cond_expr"),
    list(label = "diff_delta_combo"),
    list(label = "diff_delta_fp"),
    list(label = "diff_delta_expr"),
    list(label = "diff_fc_combo"),
    list(label = "diff_fc_fp"),
    list(label = "diff_fc_expr")
  )
  if (is.null(cond_edges_use) || !nrow(cond_edges_use)) {
    variants <- variants[!grepl("^cond_", vapply(variants, `[[`, character(1), "label"))]
  }

  run_step2_louvain <- FALSE
  run_step2_lda <- TRUE
  run_step2_nmf <- FALSE
  bench_cores <- as.integer(Sys.getenv("TOPIC_BENCH_CORES", unset = parallel::detectCores(logical = TRUE) - 1L))
  bench_cores <- max(1L, bench_cores)

  lda_K <- c(20, 40, 60, 80)
  lda_methods <- c("VEM") # c("VEM", "Gibbs")
  gamma_cuts <- c(0.7) # c(0.2, 0.5)
  coherence_top_n <- 20L
  all_topics_path <- file.path(topic_benchmark_dir, "topic_benchmark_all_topics.csv")
  lda_res_all <- list()
  nmf_res_all <- list()
  louvain_res_all <- list()

  for (v in variants) {
    docs <- build_docs_for_label(v$label)
    if (is.null(docs) || !nrow(docs)) next
    saveRDS(docs, file.path(doc_term_dir, sprintf("doc_term_%s.rds", v$label)))

    edges_use <- edges_for_label(v$label)
    if (is.null(edges_use) || !nrow(edges_use)) next

    if (isTRUE(run_step2_louvain)) {
      weight_edges <- edges_use |>
        dplyr::mutate(louvain_weight = abs(delta_link_score))

      min_weight <- NULL
      if (!grepl("_combo$", v$label)) {
        if (grepl("^diff_fc_", v$label)) {
          min_weight <- 1
        } else if (grepl("^cond_", v$label)) {
          min_weight <- 0
        } else {
          min_weight <- delta_link
        }
      }

      edge_topics <- assign_edge_topics_louvain(
        edges_tbl = weight_edges,
        weight_col = "louvain_weight",
        min_abs_weight = min_weight,
        seed = 1L
      )
      if (nrow(edge_topics)) {
        run_id <- sprintf("louvain_%s", v$label)
        topics_path <- file.path(topic_benchmark_dir, sprintf("%s_topics.csv", run_id))
        if (isTRUE(rerun_benchmark) || !file.exists(topics_path)) {
          edge_topics_path <- file.path(fit_dir, sprintf("%s_edge_topics.csv", run_id))
          readr::write_csv(edge_topics, edge_topics_path)

          group_label <- if (grepl("^cond_", v$label)) "condition" else "comparison"
          by_group <- topic_marker_overlap_from_edge_topics_by_group(
            edge_topics = edge_topics,
            group_col = "comparison_id",
            markers_epi = markers_epi,
            markers_mes = markers_mes,
            top_n_overlap_genes = 50,
            marker_mode = marker_mode,
            exclusive_method = exclusive_method,
            exclusive_top_n = exclusive_top_n,
            compute_gene_only = compute_gene_only,
            verbose = FALSE
          )
          by_group_stats_path <- file.path(
            topic_benchmark_dir,
            sprintf("%s_topics_by_%s.csv", run_id, group_label)
          )
          by_group_best_path <- file.path(
            topic_benchmark_dir,
            sprintf("%s_best_by_%s.csv", run_id, group_label)
          )
          if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
          if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

          out <- topic_marker_overlap_from_edge_topics(
            edge_topics = edge_topics,
            markers_epi = markers_epi,
            markers_mes = markers_mes,
            universe_genes = unique(edges_use$gene_key),
            top_n_overlap_genes = 50,
            marker_mode = marker_mode,
            exclusive_method = exclusive_method,
            exclusive_top_n = exclusive_top_n,
            compute_gene_only = compute_gene_only,
            verbose = FALSE
          )

          stats_tbl <- out$stats |>
            dplyr::mutate(
              purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                  overlap_epi_n / n_entities_unique, NA_real_),
              purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                  overlap_mes_n / n_entities_unique, NA_real_),
              purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
              label = v$label,
              method = "louvain",
              K = NA_integer_,
              gamma_cutoff = NA_real_
            )

          readr::write_csv(stats_tbl, topics_path)

          best_epi <- out$best_epi
          best_mes <- out$best_mes
          louvain_res_all[[length(louvain_res_all) + 1L]] <- tibble::tibble(
            label = v$label,
            method = "louvain",
            K = NA_integer_,
            gamma_cutoff = NA_real_,
            best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
            best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
            best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
              best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
            } else {
              NA_real_
            },
            best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
            best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
            best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
              best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
            } else {
              NA_real_
            },
            topics_csv = topics_path,
            edge_topics_csv = edge_topics_path,
            topics_by_group_csv = by_group_stats_path,
            best_by_group_csv = by_group_best_path
          )
        }
      }
    }

    if (isTRUE(run_step2_lda)) {
      run_one_lda <- function(row) {
        run_id <- sprintf("lda_%s_K%d_%s_g%02d", v$label, row$K, row$method, as.integer(round(row$gamma * 100)))
        topics_path <- file.path(topic_benchmark_dir, sprintf("%s_topics.csv", run_id))
        edge_topics_path <- file.path(fit_dir, sprintf("%s_edge_topics.csv", run_id))
        if (isTRUE(reuse_edge_topics) && file.exists(edge_topics_path) &&
            (isTRUE(rerun_benchmark) || !file.exists(topics_path))) {
          edge_topics <- readr::read_csv(edge_topics_path, show_col_types = FALSE)
          if (nrow(edge_topics)) {
            group_label <- if (grepl("^cond_", v$label)) "condition" else "comparison"
            by_group <- topic_marker_overlap_from_edge_topics_by_group(
              edge_topics = edge_topics,
              group_col = "comparison_id",
              markers_epi = markers_epi,
              markers_mes = markers_mes,
              top_n_overlap_genes = 50,
              marker_mode = marker_mode,
              exclusive_method = exclusive_method,
              exclusive_top_n = exclusive_top_n,
              compute_gene_only = compute_gene_only,
              verbose = FALSE
            )
            by_group_stats_path <- file.path(
              topic_benchmark_dir,
              sprintf("%s_topics_by_%s.csv", run_id, group_label)
            )
            by_group_best_path <- file.path(
              topic_benchmark_dir,
              sprintf("%s_best_by_%s.csv", run_id, group_label)
            )
            if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
            if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

            out <- topic_marker_overlap_from_edge_topics(
              edge_topics = edge_topics,
              markers_epi = markers_epi,
              markers_mes = markers_mes,
              universe_genes = unique(edges_use$gene_key),
              top_n_overlap_genes = 50,
              marker_mode = marker_mode,
              exclusive_method = exclusive_method,
              exclusive_top_n = exclusive_top_n,
              compute_gene_only = compute_gene_only,
              verbose = FALSE
            )

            stats_tbl <- out$stats |>
              dplyr::mutate(
                purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                    overlap_epi_n / n_entities_unique, NA_real_),
                purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                    overlap_mes_n / n_entities_unique, NA_real_),
                purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
                label = v$label,
                method = row$method,
                K = row$K,
                gamma_cutoff = row$gamma,
                topic_coherence = NA_real_,
                n_terms_used = NA_integer_
              )

            readr::write_csv(stats_tbl, topics_path)

            best_epi <- out$best_epi
            best_mes <- out$best_mes
            return(tibble::tibble(
              label = v$label,
              method = row$method,
              K = row$K,
              gamma_cutoff = row$gamma,
              topic_coherence_mean = NA_real_,
              topic_coherence_median = NA_real_,
              best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
              best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
              best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
                best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
              } else {
                NA_real_
              },
              best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
              best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
              best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
                best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
              } else {
                NA_real_
              },
              topics_csv = topics_path,
              edge_topics_csv = edge_topics_path,
              topics_by_group_csv = by_group_stats_path,
              best_by_group_csv = by_group_best_path
            ))
          }
        }
        if (!isTRUE(rerun_benchmark) && file.exists(topics_path)) return(NULL)

        lda_fit <- fit_pooled_lda(
          docs,
          K = row$K,
          method = row$method,
          seed = 1,
          gamma_cutoff = row$gamma
        )

        edge_topics <- assign_edge_topics_lda(
          edges_tbl = edges_use,
          lda_fit = lda_fit,
          top_n = 1L,
          aggregate_peak = aggregate_peaks
        )
        if (!nrow(edge_topics)) return(NULL)

        readr::write_csv(edge_topics, edge_topics_path)

        group_label <- if (grepl("^cond_", v$label)) "condition" else "comparison"
        by_group <- topic_marker_overlap_from_edge_topics_by_group(
          edge_topics = edge_topics,
          group_col = "comparison_id",
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )
        by_group_stats_path <- file.path(
          topic_benchmark_dir,
          sprintf("%s_topics_by_%s.csv", run_id, group_label)
        )
        by_group_best_path <- file.path(
          topic_benchmark_dir,
          sprintf("%s_best_by_%s.csv", run_id, group_label)
        )
        if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
        if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

        out <- topic_marker_overlap_from_edge_topics(
          edge_topics = edge_topics,
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          universe_genes = unique(edges_use$gene_key),
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )

        coh_tbl <- compute_lda_topic_coherence(
          doc_term_dt = docs,
          lda_fit = lda_fit,
          top_n = coherence_top_n
        )

        stats_tbl <- out$stats |>
          dplyr::mutate(
            purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_epi_n / n_entities_unique, NA_real_),
            purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_mes_n / n_entities_unique, NA_real_),
            purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
            label = v$label,
            method = row$method,
            K = row$K,
            gamma_cutoff = row$gamma
          )

        if (nrow(coh_tbl)) {
          stats_tbl <- stats_tbl |>
            dplyr::left_join(coh_tbl, by = "topic")
        } else {
          stats_tbl$topic_coherence <- NA_real_
          stats_tbl$n_terms_used <- NA_integer_
        }

        readr::write_csv(stats_tbl, topics_path)
        readr::write_csv(lda_fit$doc_topics, file.path(fit_dir, sprintf("%s_doc_topics.csv", run_id)))
        readr::write_csv(lda_fit$doc_topic_long, file.path(fit_dir, sprintf("%s_doc_topic_long.csv", run_id)))

        best_epi <- out$best_epi
        best_mes <- out$best_mes
        coherence_mean <- if (nrow(coh_tbl)) mean(coh_tbl$topic_coherence, na.rm = TRUE) else NA_real_
        coherence_median <- if (nrow(coh_tbl)) stats::median(coh_tbl$topic_coherence, na.rm = TRUE) else NA_real_
        tibble::tibble(
          label = v$label,
          method = row$method,
          K = row$K,
          gamma_cutoff = row$gamma,
          topic_coherence_mean = coherence_mean,
          topic_coherence_median = coherence_median,
          best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
          best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
          best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
            best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
          } else {
            NA_real_
          },
          best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
          best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
          best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
            best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
          } else {
            NA_real_
          },
          topics_csv = topics_path,
          edge_topics_csv = edge_topics_path,
          topics_by_group_csv = by_group_stats_path,
          best_by_group_csv = by_group_best_path
        )
      }

      lda_grid <- expand.grid(
        gamma = gamma_cuts,
        K = lda_K,
        method = lda_methods,
        stringsAsFactors = FALSE
      )
      lda_list <- split(lda_grid, seq_len(nrow(lda_grid)))
      lda_res <- if (.Platform$OS.type != "windows" && bench_cores > 1L) {
        parallel::mclapply(lda_list, run_one_lda, mc.cores = bench_cores, mc.preschedule = FALSE)
      } else {
        lapply(lda_list, run_one_lda)
      }
      lda_res_all <- c(lda_res_all, Filter(Negate(is.null), lda_res))
    }

    if (isTRUE(run_step2_nmf)) {
      nmf_K <- c(20, 30)
      nmf_methods <- c("lee")
      nmf_gamma <- c(0.3, 0.5)

      run_one_nmf <- function(row) {
        run_id <- sprintf("nmf_%s_K%d_%s_g%02d", v$label, row$K, row$method, as.integer(round(row$gamma * 100)))
        topics_path <- file.path(topic_benchmark_dir, sprintf("%s_topics.csv", run_id))
        edge_topics_path <- file.path(fit_dir, sprintf("%s_edge_topics.csv", run_id))
        if (isTRUE(reuse_edge_topics) && file.exists(edge_topics_path) &&
            (isTRUE(rerun_benchmark) || !file.exists(topics_path))) {
          edge_topics <- readr::read_csv(edge_topics_path, show_col_types = FALSE)
          if (nrow(edge_topics)) {
            group_label <- if (grepl("^cond_", v$label)) "condition" else "comparison"
            by_group <- topic_marker_overlap_from_edge_topics_by_group(
              edge_topics = edge_topics,
              group_col = "comparison_id",
              markers_epi = markers_epi,
              markers_mes = markers_mes,
              top_n_overlap_genes = 50,
              marker_mode = marker_mode,
              exclusive_method = exclusive_method,
              exclusive_top_n = exclusive_top_n,
              compute_gene_only = compute_gene_only,
              verbose = FALSE
            )
            by_group_stats_path <- file.path(
              topic_benchmark_dir,
              sprintf("%s_topics_by_%s.csv", run_id, group_label)
            )
            by_group_best_path <- file.path(
              topic_benchmark_dir,
              sprintf("%s_best_by_%s.csv", run_id, group_label)
            )
            if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
            if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

            out <- topic_marker_overlap_from_edge_topics(
              edge_topics = edge_topics,
              markers_epi = markers_epi,
              markers_mes = markers_mes,
              universe_genes = unique(edges_use$gene_key),
              top_n_overlap_genes = 50,
              marker_mode = marker_mode,
              exclusive_method = exclusive_method,
              exclusive_top_n = exclusive_top_n,
              compute_gene_only = compute_gene_only,
              verbose = FALSE
            )

            stats_tbl <- out$stats |>
              dplyr::mutate(
                purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                    overlap_epi_n / n_entities_unique, NA_real_),
                purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                    overlap_mes_n / n_entities_unique, NA_real_),
                purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
                label = v$label,
                method = row$method,
                K = row$K,
                gamma_cutoff = row$gamma
              )

            readr::write_csv(stats_tbl, topics_path)

            best_epi <- out$best_epi
            best_mes <- out$best_mes
            return(tibble::tibble(
              label = v$label,
              method = row$method,
              K = row$K,
              gamma_cutoff = row$gamma,
              best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
              best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
              best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
                best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
              } else {
                NA_real_
              },
              best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
              best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
              best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
                best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
              } else {
                NA_real_
              },
              topics_csv = topics_path,
              edge_topics_csv = edge_topics_path,
              topics_by_group_csv = by_group_stats_path,
              best_by_group_csv = by_group_best_path
            ))
          }
        }
        if (!isTRUE(rerun_benchmark) && file.exists(topics_path)) return(NULL)

        docs_run <- docs |>
          dplyr::group_by(doc_id) |>
          dplyr::slice_max(order_by = w, n = 500, with_ties = FALSE) |>
          dplyr::ungroup()

        nmf_fit <- fit_pooled_nmf(
          docs_run,
          K = row$K,
          method = row$method,
          seed = 1,
          nrun = 1,
          gamma_cutoff = row$gamma,
          nmf_stop = "none",
          nmf_options = "-p1"
        )

        edge_topics <- assign_edge_topics_nmf(
          edges_tbl = edges_use,
          nmf_fit = nmf_fit,
          top_n = 1L,
          aggregate_peak = aggregate_peaks
        )
        if (!nrow(edge_topics)) return(NULL)

        readr::write_csv(edge_topics, edge_topics_path)

        group_label <- if (grepl("^cond_", v$label)) "condition" else "comparison"
        by_group <- topic_marker_overlap_from_edge_topics_by_group(
          edge_topics = edge_topics,
          group_col = "comparison_id",
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )
        by_group_stats_path <- file.path(
          topic_benchmark_dir,
          sprintf("%s_topics_by_%s.csv", run_id, group_label)
        )
        by_group_best_path <- file.path(
          topic_benchmark_dir,
          sprintf("%s_best_by_%s.csv", run_id, group_label)
        )
        if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
        if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

        out <- topic_marker_overlap_from_edge_topics(
          edge_topics = edge_topics,
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          universe_genes = unique(edges_use$gene_key),
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )

        stats_tbl <- out$stats |>
          dplyr::mutate(
            purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_epi_n / n_entities_unique, NA_real_),
            purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_mes_n / n_entities_unique, NA_real_),
            purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
            label = v$label,
            method = row$method,
            K = row$K,
            gamma_cutoff = row$gamma
          )

        readr::write_csv(stats_tbl, topics_path)
        readr::write_csv(nmf_fit$doc_topics, file.path(fit_dir, sprintf("%s_doc_topics.csv", run_id)))
        readr::write_csv(nmf_fit$doc_topic_long, file.path(fit_dir, sprintf("%s_doc_topic_long.csv", run_id)))

        best_epi <- out$best_epi
        best_mes <- out$best_mes
        tibble::tibble(
          label = v$label,
          method = row$method,
          K = row$K,
          gamma_cutoff = row$gamma,
          best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
          best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
          best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
            best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
          } else {
            NA_real_
          },
          best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
          best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
          best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
            best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
          } else {
            NA_real_
          },
          topics_csv = topics_path,
          edge_topics_csv = edge_topics_path,
          topics_by_group_csv = by_group_stats_path,
          best_by_group_csv = by_group_best_path
        )
      }

      nmf_grid <- expand.grid(
        gamma = nmf_gamma,
        K = nmf_K,
        method = nmf_methods,
        stringsAsFactors = FALSE
      )
      nmf_list <- split(nmf_grid, seq_len(nrow(nmf_grid)))
      nmf_res <- if (.Platform$OS.type != "windows" && bench_cores > 1L) {
        parallel::mclapply(nmf_list, run_one_nmf, mc.cores = bench_cores, mc.preschedule = FALSE)
      } else {
        lapply(nmf_list, run_one_nmf)
      }
      nmf_res_all <- c(nmf_res_all, Filter(Negate(is.null), nmf_res))
    }

    rm(docs)
    gc()
  }

  if (isTRUE(run_step2_louvain)) {
    louvain_sum <- dplyr::bind_rows(louvain_res_all)
    if (nrow(louvain_sum)) {
      readr::write_csv(louvain_sum, file.path(topic_benchmark_dir, "louvain_all_runs_summary.csv"))
    }
  }
  if (isTRUE(run_step2_lda)) {
    lda_sum <- dplyr::bind_rows(lda_res_all)
    if (nrow(lda_sum)) {
      readr::write_csv(lda_sum, file.path(topic_benchmark_dir, "lda_all_runs_summary.csv"))
    }
  }
  if (isTRUE(run_step2_nmf)) {
    nmf_sum <- dplyr::bind_rows(nmf_res_all)
    if (nrow(nmf_sum)) {
      readr::write_csv(nmf_sum, file.path(topic_benchmark_dir, "nmf_all_runs_summary.csv"))
    }
  }

  all_topics_tbl <- summarize_topic_benchmark_edge_topics(
    topic_benchmark_dir = topic_benchmark_dir,
    fit_dir = fit_dir,
    markers_epi = markers_epi,
    markers_mes = markers_mes,
    out_file = all_topics_path,
    top_n_overlap_genes = 50L,
    marker_mode = marker_mode,
    exclusive_method = exclusive_method,
    exclusive_top_n = exclusive_top_n,
    compute_gene_only = compute_gene_only
  )

  write_topic_scatter <- function(tbl, out_dir) {
    if (!nrow(tbl)) return(invisible(NULL))
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Plotting requires ggplot2.")
    }

    plot_tbl <- tbl |>
      dplyr::filter(
        is.finite(.data$epi_marker_coverage),
        is.finite(.data$mes_marker_coverage),
        is.finite(.data$purity_epi),
        is.finite(.data$purity_mes)
      ) |>
      dplyr::mutate(
        label_pretty = gsub("_", " ", gsub("^diff_", "diff ", gsub("^cond_", "cond ", .data$label))),
        K_plot_label = ifelse(is.na(.data$K), "K = NA", paste0("K = ", .data$K))
      )
    if (!nrow(plot_tbl)) return(invisible(NULL))
    k_levels <- sort(unique(plot_tbl$K[is.finite(plot_tbl$K)]))
    k_labels <- paste0("K = ", k_levels)
    if (any(is.na(plot_tbl$K))) k_labels <- c(k_labels, "K = NA")
    plot_tbl$K_plot_label <- factor(plot_tbl$K_plot_label, levels = k_labels)
    plot_tbl$label_pretty <- factor(plot_tbl$label_pretty, levels = unique(plot_tbl$label_pretty))

    base_theme <- ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "grey40"),
        strip.text.y = ggplot2::element_text(angle = 0)
      )

    max_epi_n <- max(plot_tbl$epi_marker_coverage, na.rm = TRUE)
    max_mes_n <- max(plot_tbl$mes_marker_coverage, na.rm = TRUE)
    if (!is.finite(max_epi_n) || max_epi_n < 0) max_epi_n <- 0
    if (!is.finite(max_mes_n) || max_mes_n < 0) max_mes_n <- 0

    p_epi <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes(x = .data$epi_marker_coverage, y = .data$purity_epi)
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 1.0) +
      ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
      ggplot2::scale_x_continuous(
        limits = c(0, max_epi_n + 0.01),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = "EPI marker coverage",
        y = "Epithelial purity"
      ) +
      base_theme

    p_mes <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes(x = .data$mes_marker_coverage, y = .data$purity_mes)
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 1.0) +
      ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
      ggplot2::scale_x_continuous(
        limits = c(0, max_mes_n + 0.01),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = "MES marker coverage",
        y = "Mesenchymal purity"
      ) +
      base_theme

    ggplot2::ggsave(
      filename = file.path(out_dir, "topic_scatter_epi.pdf"),
      plot = p_epi,
      width = 16,
      height = 9,
      units = "in"
    )
    ggplot2::ggsave(
      filename = file.path(out_dir, "topic_scatter_mes.pdf"),
      plot = p_mes,
      width = 16,
      height = 9,
      units = "in"
    )

    invisible(plot_tbl)
  }

  write_topic_scatter(all_topics_tbl, topic_benchmark_dir)

  merge_dir <- file.path(topic_benchmark_dir, "merged")
  fit_dir_merged <- file.path(merge_dir, "fits")
  dir.create(merge_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir_merged, recursive = TRUE, showWarnings = FALSE)

  parse_run_meta <- function(run_id) {
    if (grepl("^lda_", run_id)) {
      m <- regexec("^lda_(.+)_K(\\d+)_([^_]+)_g(\\d+)$", run_id)
      hit <- regmatches(run_id, m)[[1]]
      if (length(hit)) {
        return(list(label = hit[2], method = hit[4], K = as.integer(hit[3]),
                    gamma_cutoff = as.numeric(hit[5]) / 100))
      }
    }
    if (grepl("^nmf_", run_id)) {
      m <- regexec("^nmf_(.+)_K(\\d+)_([^_]+)_g(\\d+)$", run_id)
      hit <- regmatches(run_id, m)[[1]]
      if (length(hit)) {
        return(list(label = hit[2], method = hit[4], K = as.integer(hit[3]),
                    gamma_cutoff = as.numeric(hit[5]) / 100))
      }
    }
    if (grepl("^louvain_", run_id)) {
      label <- sub("^louvain_", "", run_id)
      return(list(label = label, method = "louvain", K = NA_integer_, gamma_cutoff = NA_real_))
    }
    list(label = run_id, method = NA_character_, K = NA_integer_, gamma_cutoff = NA_real_)
  }

  edge_files <- list.files(fit_dir, pattern = "_edge_topics\\.csv$", full.names = TRUE)
  if (length(edge_files)) {
    for (edge_path in edge_files) {
      run_id <- sub("_edge_topics\\.csv$", "", basename(edge_path))
      merged_edge_path <- file.path(fit_dir_merged, sprintf("%s_edge_topics.csv", run_id))
      topics_path <- file.path(merge_dir, sprintf("%s_topics.csv", run_id))

      edge_topics_merged <- NULL
      if (isTRUE(rerun_benchmark) || !file.exists(merged_edge_path)) {
        edge_topics <- readr::read_csv(edge_path, show_col_types = FALSE)
        if (!nrow(edge_topics)) next
        merge_res <- merge_edge_topics_by_overlap(
          edge_topics = edge_topics,
          overlap_thresh = merge_overlap_thresh
        )
        edge_topics_merged <- merge_res$edge_topics
        readr::write_csv(edge_topics_merged, merged_edge_path)
        readr::write_csv(
          merge_res$mapping,
          file.path(fit_dir_merged, sprintf("%s_topic_merge_map.csv", run_id))
        )
        readr::write_csv(
          merge_res$topic_genes,
          file.path(fit_dir_merged, sprintf("%s_topic_merge_genes.csv", run_id))
        )
      } else {
        edge_topics_merged <- readr::read_csv(merged_edge_path, show_col_types = FALSE)
      }
      if (is.null(edge_topics_merged) || !nrow(edge_topics_merged)) next

      if (isTRUE(rerun_benchmark) || !file.exists(topics_path)) {
        meta <- parse_run_meta(run_id)
        group_label <- if (grepl("^cond_", meta$label)) "condition" else "comparison"
        by_group_stats_path <- file.path(merge_dir, sprintf("%s_topics_by_%s.csv", run_id, group_label))
        by_group_best_path <- file.path(merge_dir, sprintf("%s_best_by_%s.csv", run_id, group_label))
        by_group <- topic_marker_overlap_from_edge_topics_by_group(
          edge_topics = edge_topics_merged,
          group_col = "comparison_id",
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )
        if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
        if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

        out <- topic_marker_overlap_from_edge_topics(
          edge_topics = edge_topics_merged,
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          universe_genes = unique(edge_topics_merged$gene_key),
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )

        stats_tbl <- out$stats |>
          dplyr::mutate(
            purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_epi_n / n_entities_unique, NA_real_),
            purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_mes_n / n_entities_unique, NA_real_),
            purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
            label = meta$label,
            method = meta$method,
            K = meta$K,
            gamma_cutoff = meta$gamma_cutoff,
            run_id = run_id
          )

        readr::write_csv(stats_tbl, topics_path)
      }
    }
  }

  merged_all_topics_tbl <- summarize_topic_benchmark_edge_topics(
    topic_benchmark_dir = merge_dir,
    fit_dir = fit_dir_merged,
    markers_epi = markers_epi,
    markers_mes = markers_mes,
    out_file = file.path(merge_dir, "topic_benchmark_all_topics.csv"),
    top_n_overlap_genes = 50L,
    marker_mode = marker_mode,
    exclusive_method = exclusive_method,
    exclusive_top_n = exclusive_top_n,
    compute_gene_only = compute_gene_only
  )

  write_topic_scatter(merged_all_topics_tbl, merge_dir)
}
}

do_topic_benchmark_sig <- TRUE
source("R/utils_grn_lda_nmf.R")
# Topic benchmarking (experimental directional, sig comparisons)
if (isTRUE(do_topic_benchmark_sig)) {
  source("R/utils_grn_lda_nmf.R")
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Topic benchmarking requires data.table.")
  }

  topic_runs_sig <- list(
    list(tag = "per_link_aggregated_sig_comparison_vem", aggregate_peaks = TRUE, lda_method = "VEM"),
    list(tag = "per_link_aggregated_sig_comparison_gibbs", aggregate_peaks = TRUE, lda_method = "Gibbs"),
    list(tag = "per_link_nonaggregated_sig_comparison_vem", aggregate_peaks = FALSE, lda_method = "VEM"),
    list(tag = "per_link_nonaggregated_sig_comparison_gibbs", aggregate_peaks = FALSE, lda_method = "Gibbs")
  )

  markers_dir <- file.path(base_dir, "plot_lineage_plasticity_related_subnetworks", "markers")
  markers_epi <- readr::read_tsv(
    file.path(markers_dir, "epithelial.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  markers_mes <- readr::read_tsv(
    file.path(markers_dir, "mesenchymal.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  if (!length(markers_epi) || !length(markers_mes)) {
    stop("markers_epi/markers_mes are empty; check marker files.")
  }

  marker_mode <- "union"
  exclusive_method <- "mean_score"
  exclusive_top_n <- 1L
  compute_gene_only <- FALSE

  lda_K <- c(20, 40, 60, 80)
  lda_methods <- c("VEM", "Gibbs")
  gamma_cuts <- c(0.7)
  bench_cores <- 20L
  topic_benchmark_sig_mode <- "per_cellline" # "all" or "per_cellline"
  run_cfg_parallel <- TRUE
  run_cfg_workers <- min(length(topic_runs_sig), max(1L, floor(bench_cores / 2L)))
  bench_cores_run <- if (isTRUE(run_cfg_parallel) &&
                         .Platform$OS.type != "windows" &&
                         run_cfg_workers > 1L) {
    max(1L, floor(bench_cores / run_cfg_workers))
  } else {
    bench_cores
  }

  top_terms_per_doc <- 500L
  min_df <- 2L
  log2fc_eps <- 1e-6
  target_log2fc_abs_min <- 1

  minmax_scale_vec <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) {
      return(rep(0, length(x)))
    }
    (x - rng[1]) / diff(rng)
  }

  prep_directional_combo <- function(edges_tbl,
                                     target_col,
                                     tf_col,
                                     fp_col) {
    empty_tbl <- tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      w = numeric(0)
    )
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(list(multi = empty_tbl, sum = empty_tbl))
    dt <- data.table::as.data.table(edges_tbl)
    need <- c(target_col, tf_col, fp_col)
    miss <- setdiff(need, names(dt))
    if (length(miss)) return(list(multi = empty_tbl, sum = empty_tbl))
    if (!("peak_id" %in% names(dt))) dt[, peak_id := NA_character_]

    dt <- dt[is.finite(get(target_col)) & get(target_col) != 0]
    if (!nrow(dt)) return(list(multi = empty_tbl, sum = empty_tbl))

    dt[, sign_target := sign(get(target_col))]
    dt <- dt[sign_target != 0]
    dt[, target_adj := abs(get(target_col))]
    dt[, tf_adj := get(tf_col) * sign_target]
    dt[, fp_adj := get(fp_col) * sign_target]

    dt[, tf_scaled := minmax_scale_vec(tf_adj), by = sign_target]
    dt[, fp_scaled := minmax_scale_vec(fp_adj), by = sign_target]
    dt[, target_scaled := minmax_scale_vec(target_adj), by = sign_target]

    dt[, w_multi := tf_scaled * target_scaled * fp_scaled]
    dt[, w_sum := tf_scaled + target_scaled + fp_scaled]

    dt_multi <- dt[is.finite(w_multi) & w_multi > 0,
                   .(comparison_id, tf, gene_key, peak_id, w = w_multi)]
    dt_sum <- dt[is.finite(w_sum) & w_sum > 0,
                 .(comparison_id, tf, gene_key, peak_id, w = w_sum)]

    list(multi = tibble::as_tibble(dt_multi), sum = tibble::as_tibble(dt_sum))
  }

  prep_fp_only <- function(edges_tbl, target_col, fp_col) {
    empty_tbl <- tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      w = numeric(0)
    )
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(empty_tbl)
    dt <- data.table::as.data.table(edges_tbl)
    need <- c(target_col, fp_col)
    miss <- setdiff(need, names(dt))
    if (length(miss)) return(empty_tbl)
    if (!("peak_id" %in% names(dt))) dt[, peak_id := NA_character_]

    dt <- dt[is.finite(get(target_col)) & get(target_col) != 0]
    if (!nrow(dt)) return(empty_tbl)
    dt[, w := abs(get(fp_col))]
    dt <- dt[is.finite(w) & w > 0]
    tibble::as_tibble(dt[, .(comparison_id, tf, gene_key, peak_id, w)])
  }

  prep_expr_only <- function(edges_tbl,
                             target_col,
                             tf_col,
                             sig_col,
                             sig_abs_min,
                             sum_mode = FALSE) {
    empty_tbl <- tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      w = numeric(0)
    )
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(empty_tbl)
    dt <- data.table::as.data.table(edges_tbl)
    need <- c(target_col, tf_col, sig_col)
    miss <- setdiff(need, names(dt))
    if (length(miss)) return(empty_tbl)

    dt <- dt[is.finite(get(sig_col)) & abs(get(sig_col)) > sig_abs_min]
    if (!nrow(dt)) return(empty_tbl)
    dt <- unique(dt, by = c("comparison_id", "tf", "gene_key"))
    dt[, peak_id := NA_character_]

    dt <- dt[is.finite(get(target_col)) & get(target_col) != 0]
    if (!nrow(dt)) return(empty_tbl)
    dt[, sign_target := sign(get(target_col))]
    dt <- dt[sign_target != 0]
    if (!nrow(dt)) return(empty_tbl)

    if (isTRUE(sum_mode)) {
      dt[, w_raw := get(tf_col) + get(target_col)]
      dt[sign_target < 0, w_raw := -w_raw]
    } else {
      dt[, w_raw := get(tf_col) * get(target_col)]
    }
    dt[, w := minmax_scale_vec(w_raw), by = sign_target]
    dt <- dt[is.finite(w) & w > 0]
    tibble::as_tibble(dt[, .(comparison_id, tf, gene_key, peak_id, w)])
  }

  build_doc_term_from_edges <- function(edges_tbl,
                                        top_terms_per_doc = 500L,
                                        min_df = 2L) {
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(data.table::data.table())
    dt <- data.table::as.data.table(edges_tbl)
    if (!("w" %in% names(dt))) return(data.table::data.table())
    dt <- dt[is.finite(w) & w > 0]
    if (!nrow(dt)) return(data.table::data.table())
    dt[, doc_id := paste(comparison_id, tf, sep = "::")]
    dt <- dt[, .(w = sum(w)), by = .(doc_id, comparison_id, tf, gene_key)]

    if (!is.null(top_terms_per_doc) && is.finite(top_terms_per_doc) && top_terms_per_doc > 0) {
      data.table::setorder(dt, doc_id, -w)
      dt <- dt[, head(.SD, top_terms_per_doc), by = doc_id]
    }

    df_tbl <- unique(dt[, .(doc_id, gene_key)])
    term_df <- df_tbl[, .N, by = gene_key]
    keep_terms <- term_df[N >= min_df, gene_key]
    dt <- dt[gene_key %in% keep_terms]
    dt[]
  }

  edges_all_tidy <- load_delta_links_all_tidy(step2_delta_csvs)
  edges_all_tidy <- edges_all_tidy |>
    dplyr::mutate(
      fp_gene_case = dplyr::coalesce(fp_bound_case, 0L) > 0L &
        dplyr::coalesce(gene_expr_flag_case, 0L) >= 1L,
      fp_gene_ctrl = dplyr::coalesce(fp_bound_ctrl, 0L) > 0L &
        dplyr::coalesce(gene_expr_flag_ctrl, 0L) >= 1L,
      tf_expr_either = dplyr::coalesce(tf_expr_flag_case, 0L) >= 1L |
        dplyr::coalesce(tf_expr_flag_ctrl, 0L) >= 1L,
      pass_topic_gate = tf_expr_either & (fp_gene_case | fp_gene_ctrl)
    ) |>
    dplyr::filter(pass_topic_gate)

  if (!("log2FC_tf_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2FC_tf_expr <- NA_real_
  if (!("log2FC_gene_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2FC_gene_expr <- NA_real_
  if (!("log2fc_tf_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2fc_tf_expr <- NA_real_
  if (!("log2fc_gene_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2fc_gene_expr <- NA_real_
  edges_all_tidy <- edges_all_tidy |>
    dplyr::mutate(
      log2fc_tf_expr = dplyr::coalesce(.data$log2FC_tf_expr, .data$log2fc_tf_expr),
      log2fc_gene_expr = dplyr::coalesce(.data$log2FC_gene_expr, .data$log2fc_gene_expr),
      log2fc_fp_bed_score = log2((.data$fp_bed_score_case + log2fc_eps) /
                                   (.data$fp_bed_score_ctrl + log2fc_eps))
    )

  sum_cols_sig <- intersect(
    c("delta_fp_bed_score", "log2fc_fp_bed_score"),
    names(edges_all_tidy)
  )
  first_cols_sig <- intersect(
    c("comparison_id", "tf", "gene_key",
      "delta_tf_expr", "delta_gene_expr",
      "log2fc_tf_expr", "log2fc_gene_expr"),
    names(edges_all_tidy)
  )

  run_sig_one <- function(run_cfg, edges_tbl = edges_all_tidy, run_suffix = NULL) {
    lda_res_all <- list()
    run_tag <- if (is.null(run_suffix)) run_cfg$tag else sprintf("%s_%s", run_cfg$tag, run_suffix)
    topic_benchmark_dir <- file.path(step2_out_dir, sprintf("topic_benchmark_%s", run_tag))
    doc_term_dir <- file.path(topic_benchmark_dir, "doc_terms")
    fit_dir <- file.path(topic_benchmark_dir, "fits")
    dir.create(topic_benchmark_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(doc_term_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

    edges_run <- if (isTRUE(run_cfg$aggregate_peaks)) {
      edges_agg <- aggregate_edges_by_tf_gene(
        edges_tbl,
        sum_cols = sum_cols_sig,
        first_cols = first_cols_sig
      )
      edges_agg$peak_id <- NA_character_
      edges_agg
    } else {
      edges_tbl
    }

    combo_delta <- prep_directional_combo(
      edges_run,
      target_col = "delta_gene_expr",
      tf_col = "delta_tf_expr",
      fp_col = "delta_fp_bed_score"
    )
    combo_fc <- prep_directional_combo(
      edges_run,
      target_col = "log2fc_gene_expr",
      tf_col = "log2fc_tf_expr",
      fp_col = "log2fc_fp_bed_score"
    )

    variants <- list(
      list(label = "diff_delta_combo_multi", edges = combo_delta$multi),
      list(label = "diff_delta_combo_sum", edges = combo_delta$sum),
      list(label = "diff_delta_fp", edges = prep_fp_only(edges_run, "delta_gene_expr", "delta_fp_bed_score")),
      list(label = "diff_fc_combo_multi", edges = combo_fc$multi),
      list(label = "diff_fc_combo_sum", edges = combo_fc$sum),
      list(label = "diff_fc_fp", edges = prep_fp_only(edges_run, "log2fc_gene_expr", "log2fc_fp_bed_score"))
    )

    if (isTRUE(run_cfg$aggregate_peaks)) {
      variants <- c(
        variants,
        list(
          list(label = "diff_delta_expr_multi",
               edges = prep_expr_only(edges_run, "delta_gene_expr", "delta_tf_expr",
                                      "log2fc_gene_expr", target_log2fc_abs_min, sum_mode = FALSE)),
          list(label = "diff_delta_expr_sum",
               edges = prep_expr_only(edges_run, "delta_gene_expr", "delta_tf_expr",
                                      "log2fc_gene_expr", target_log2fc_abs_min, sum_mode = TRUE)),
          list(label = "diff_fc_expr_multi",
               edges = prep_expr_only(edges_run, "log2fc_gene_expr", "log2fc_tf_expr",
                                      "log2fc_gene_expr", target_log2fc_abs_min, sum_mode = FALSE)),
          list(label = "diff_fc_expr_sum",
               edges = prep_expr_only(edges_run, "log2fc_gene_expr", "log2fc_tf_expr",
                                      "log2fc_gene_expr", target_log2fc_abs_min, sum_mode = TRUE))
        )
      )
    }

    for (v in variants) {
      if (is.null(v$edges) || !nrow(v$edges)) next

      docs <- build_doc_term_from_edges(
        v$edges,
        top_terms_per_doc = top_terms_per_doc,
        min_df = min_df
      )
      if (!nrow(docs)) next

      doc_term_path <- file.path(doc_term_dir, sprintf("doc_term_%s.rds", v$label))
      saveRDS(docs, doc_term_path)

      lda_methods_run <- run_cfg$lda_method
      lda_grid <- expand.grid(
        gamma = gamma_cuts,
        K = lda_K,
        method = lda_methods_run,
        stringsAsFactors = FALSE
      )
      lda_list <- split(lda_grid, seq_len(nrow(lda_grid)))

      run_one_lda <- function(row) {
        run_id <- sprintf("lda_%s_K%d_%s_g%d",
                          v$label, row$K, row$method, as.integer(row$gamma * 100))
        edge_topics_path <- file.path(fit_dir, sprintf("%s_edge_topics.csv", run_id))
        topics_path <- file.path(topic_benchmark_dir, sprintf("%s_topics.csv", run_id))
        by_group_stats_path <- file.path(topic_benchmark_dir, sprintf("%s_topics_by_comparison.csv", run_id))
        by_group_best_path <- file.path(topic_benchmark_dir, sprintf("%s_best_by_comparison.csv", run_id))

        if (!isTRUE(rerun_benchmark) &&
            isTRUE(reuse_edge_topics) &&
            file.exists(edge_topics_path) &&
            file.exists(topics_path)) {
          return(NULL)
        }

        lda_fit <- fit_pooled_lda(
          docs,
          K = row$K,
          method = row$method,
          seed = 1,
          gamma_cutoff = row$gamma
        )

        edge_topics <- assign_edge_topics_lda(
          edges_tbl = v$edges,
          lda_fit = lda_fit,
          top_n = 1L,
          aggregate_peak = run_cfg$aggregate_peaks
        )
        if (!nrow(edge_topics)) return(NULL)

        readr::write_csv(edge_topics, edge_topics_path)

        by_group <- topic_marker_overlap_from_edge_topics_by_group(
          edge_topics = edge_topics,
          group_col = "comparison_id",
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )
        if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
        if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

        out <- topic_marker_overlap_from_edge_topics(
          edge_topics = edge_topics,
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          universe_genes = unique(edge_topics$gene_key),
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )

        stats_tbl <- out$stats |>
          dplyr::mutate(
            purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_epi_n / n_entities_unique, NA_real_),
            purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_mes_n / n_entities_unique, NA_real_),
            purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
            label = v$label,
            method = row$method,
            K = row$K,
            gamma_cutoff = row$gamma
          )

        readr::write_csv(stats_tbl, topics_path)

        best_epi <- out$best_epi
        best_mes <- out$best_mes
        tibble::tibble(
          label = v$label,
          method = row$method,
          K = row$K,
          gamma_cutoff = row$gamma,
          best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
          best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
          best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
            best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
          } else {
            NA_real_
          },
          best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
          best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
          best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
            best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
          } else {
            NA_real_
          },
          topics_csv = topics_path,
          edge_topics_csv = edge_topics_path,
          topics_by_group_csv = by_group_stats_path,
          best_by_group_csv = by_group_best_path
        )
      }

      lda_res <- if (.Platform$OS.type != "windows" && bench_cores_run > 1L) {
        parallel::mclapply(lda_list, run_one_lda, mc.cores = bench_cores_run, mc.preschedule = FALSE)
      } else {
        lapply(lda_list, run_one_lda)
      }
      lda_res_all <- c(lda_res_all, Filter(Negate(is.null), lda_res))
      rm(docs)
      gc()
    }

    lda_sum <- dplyr::bind_rows(lda_res_all)
    if (nrow(lda_sum)) {
      readr::write_csv(lda_sum, file.path(topic_benchmark_dir, "lda_all_runs_summary.csv"))
    }

    all_topics_tbl <- summarize_topic_benchmark_edge_topics(
      topic_benchmark_dir = topic_benchmark_dir,
      fit_dir = fit_dir,
      markers_epi = markers_epi,
      markers_mes = markers_mes,
      out_file = file.path(topic_benchmark_dir, "topic_benchmark_all_topics.csv"),
      top_n_overlap_genes = 50L,
      marker_mode = marker_mode,
      exclusive_method = exclusive_method,
      exclusive_top_n = exclusive_top_n,
      compute_gene_only = compute_gene_only
    )

    write_topic_scatter_sig <- function(tbl, out_dir) {
      if (!nrow(tbl)) return(invisible(NULL))
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Plotting requires ggplot2.")
      }

      plot_tbl <- tbl |>
        dplyr::filter(
          is.finite(.data$epi_marker_coverage),
          is.finite(.data$mes_marker_coverage),
          is.finite(.data$purity_epi),
          is.finite(.data$purity_mes)
        ) |>
        dplyr::mutate(
          label_pretty = gsub("_", " ", gsub("^diff_", "diff ", gsub("^cond_", "cond ", .data$label))),
          K_plot_label = ifelse(is.na(.data$K), "K = NA", paste0("K = ", .data$K))
        )
      if (!nrow(plot_tbl)) return(invisible(NULL))

      k_levels <- sort(unique(plot_tbl$K[is.finite(plot_tbl$K)]))
      k_labels <- paste0("K = ", k_levels)
      if (any(is.na(plot_tbl$K))) k_labels <- c(k_labels, "K = NA")
      plot_tbl$K_plot_label <- factor(plot_tbl$K_plot_label, levels = k_labels)
      plot_tbl$label_pretty <- factor(plot_tbl$label_pretty, levels = unique(plot_tbl$label_pretty))

      base_theme <- ggplot2::theme_bw() +
        ggplot2::theme(
          text = ggplot2::element_text(face = "bold"),
          axis.text = ggplot2::element_text(face = "bold"),
          axis.title = ggplot2::element_text(face = "bold"),
          strip.text = ggplot2::element_text(face = "bold"),
          strip.background = ggplot2::element_rect(fill = "grey90", color = "grey40"),
          strip.text.y = ggplot2::element_text(angle = 0)
        )

      max_epi_n <- max(plot_tbl$epi_marker_coverage, na.rm = TRUE)
      max_mes_n <- max(plot_tbl$mes_marker_coverage, na.rm = TRUE)
      if (!is.finite(max_epi_n) || max_epi_n < 0) max_epi_n <- 0
      if (!is.finite(max_mes_n) || max_mes_n < 0) max_mes_n <- 0

      p_epi <- ggplot2::ggplot(
        plot_tbl,
        ggplot2::aes(x = .data$epi_marker_coverage, y = .data$purity_epi)
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.0) +
        ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
        ggplot2::scale_x_continuous(
          limits = c(0, max_epi_n + 0.01),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          x = "EPI marker coverage",
          y = "Epithelial purity"
        ) +
        base_theme

      p_mes <- ggplot2::ggplot(
        plot_tbl,
        ggplot2::aes(x = .data$mes_marker_coverage, y = .data$purity_mes)
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.0) +
        ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
        ggplot2::scale_x_continuous(
          limits = c(0, max_mes_n + 0.01),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          x = "MES marker coverage",
          y = "Mesenchymal purity"
        ) +
        base_theme

      ggplot2::ggsave(
        filename = file.path(out_dir, "topic_scatter_epi.pdf"),
        plot = p_epi,
        width = 16,
        height = 9,
        units = "in"
      )
      ggplot2::ggsave(
        filename = file.path(out_dir, "topic_scatter_mes.pdf"),
        plot = p_mes,
        width = 16,
        height = 9,
        units = "in"
      )

      invisible(plot_tbl)
    }

    write_topic_scatter_sig(all_topics_tbl, topic_benchmark_dir)
  }

  if (identical(topic_benchmark_sig_mode, "per_cellline")) {
    edges_all_tidy <- edges_all_tidy |>
      dplyr::mutate(cell_line = sub("_.*$", "", .data$comparison_id))
    edges_by_cell <- split(edges_all_tidy, edges_all_tidy$cell_line)
    for (cell_line in names(edges_by_cell)) {
      edges_cell <- edges_by_cell[[cell_line]]
      if (!nrow(edges_cell)) next
      if (isTRUE(run_cfg_parallel) && .Platform$OS.type != "windows" && run_cfg_workers > 1L) {
        parallel::mclapply(
          topic_runs_sig,
          function(cfg) run_sig_one(cfg, edges_tbl = edges_cell, run_suffix = cell_line),
          mc.cores = run_cfg_workers,
          mc.preschedule = FALSE
        )
      } else {
        lapply(
          topic_runs_sig,
          function(cfg) run_sig_one(cfg, edges_tbl = edges_cell, run_suffix = cell_line)
        )
      }
    }
  } else if (isTRUE(run_cfg_parallel) && .Platform$OS.type != "windows" && run_cfg_workers > 1L) {
    parallel::mclapply(topic_runs_sig, run_sig_one, mc.cores = run_cfg_workers, mc.preschedule = FALSE)
  } else {
    lapply(topic_runs_sig, run_sig_one)
  }
}

do_topic_benchmark_both_fp_rna <- TRUE
source("R/utils_grn_lda_nmf.R")
if (isTRUE(do_topic_benchmark_both_fp_rna)) {
  source("R/utils_grn_lda_nmf.R")
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Topic benchmarking requires data.table.")
  }

  topic_runs_fp_rna <- list(
    list(tag = "both_fp_rna_aggregated_sig_comparison_vem", aggregate_peaks = TRUE),
    list(tag = "both_fp_rna_nonaggregated_sig_comparison_vem", aggregate_peaks = FALSE)
  )

  markers_dir <- file.path(base_dir, "plot_lineage_plasticity_related_subnetworks", "markers")
  markers_epi <- readr::read_tsv(
    file.path(markers_dir, "epithelial.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  markers_mes <- readr::read_tsv(
    file.path(markers_dir, "mesenchymal.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  if (!length(markers_epi) || !length(markers_mes)) {
    stop("markers_epi/markers_mes are empty; check marker files.")
  }

  marker_mode <- "union"
  exclusive_method <- "mean_score"
  exclusive_top_n <- 1L
  compute_gene_only <- FALSE

  lda_K <- c(20, 40, 60, 80)
  gamma_cuts <- c(0.7)
  bench_cores <- 20L
  run_cfg_parallel <- TRUE
  run_cfg_workers <- min(length(topic_runs_fp_rna), max(1L, floor(bench_cores / 2L)))
  bench_cores_run <- if (isTRUE(run_cfg_parallel) &&
                         .Platform$OS.type != "windows" &&
                         run_cfg_workers > 1L) {
    max(1L, floor(bench_cores / run_cfg_workers))
  } else {
    bench_cores
  }

  top_terms_per_doc <- 500L
  min_df <- 2L
  log2fc_eps <- 1e-6

  minmax_scale_vec <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= 0) {
      return(rep(0, length(x)))
    }
    (x - rng[1]) / diff(rng)
  }

  prep_fp_only <- function(edges_tbl, target_col, fp_col) {
    empty_tbl <- tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      w = numeric(0)
    )
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(empty_tbl)
    dt <- data.table::as.data.table(edges_tbl)
    need <- c(target_col, fp_col)
    miss <- setdiff(need, names(dt))
    if (length(miss)) return(empty_tbl)
    if (!("peak_id" %in% names(dt))) dt[, peak_id := NA_character_]

    dt <- dt[is.finite(get(target_col)) & get(target_col) != 0]
    if (!nrow(dt)) return(empty_tbl)
    dt[, w := abs(get(fp_col))]
    dt <- dt[is.finite(w) & w > 0]
    tibble::as_tibble(dt[, .(comparison_id, tf, gene_key, peak_id, w)])
  }

  prep_fp_rna_edges <- function(edges_tbl, target_col, fp_col) {
    empty_tbl <- tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      w = numeric(0)
    )
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(empty_tbl)
    dt <- data.table::as.data.table(edges_tbl)
    need <- c(target_col, fp_col)
    miss <- setdiff(need, names(dt))
    if (length(miss)) return(empty_tbl)
    if (!("peak_id" %in% names(dt))) dt[, peak_id := NA_character_]

    dt <- dt[is.finite(get(target_col)) & get(target_col) != 0]
    if (!nrow(dt)) return(empty_tbl)

    dt_fp <- dt[is.finite(get(fp_col)) & get(fp_col) != 0]
    dt_fp[, w := abs(get(fp_col))]
    dt_fp <- dt_fp[is.finite(w) & w > 0]
    dt_fp[, gene_key := paste(gene_key, "FP", sep = "|")]

    dt_rna <- dt[is.finite(get(target_col)) & get(target_col) != 0]
    dt_rna[, w := abs(get(target_col))]
    dt_rna <- dt_rna[is.finite(w) & w > 0]
    dt_rna[, gene_key := paste(gene_key, "RNA", sep = "|")]

    out <- data.table::rbindlist(
      list(
        dt_fp[, .(comparison_id, tf, gene_key, peak_id, w)],
        dt_rna[, .(comparison_id, tf, gene_key, peak_id, w)]
      ),
      use.names = TRUE,
      fill = TRUE
    )
    if (!nrow(out)) return(empty_tbl)
    tibble::as_tibble(out)
  }

  build_doc_term_from_edges <- function(edges_tbl,
                                        top_terms_per_doc = 500L,
                                        min_df = 2L) {
    if (is.null(edges_tbl) || !nrow(edges_tbl)) return(data.table::data.table())
    dt <- data.table::as.data.table(edges_tbl)
    if (!("w" %in% names(dt))) return(data.table::data.table())
    dt <- dt[is.finite(w) & w > 0]
    if (!nrow(dt)) return(data.table::data.table())
    dt[, doc_id := paste(comparison_id, tf, sep = "::")]
    dt <- dt[, .(w = sum(w)), by = .(doc_id, comparison_id, tf, gene_key)]

    if (!is.null(top_terms_per_doc) && is.finite(top_terms_per_doc) && top_terms_per_doc > 0) {
      data.table::setorder(dt, doc_id, -w)
      dt <- dt[, head(.SD, top_terms_per_doc), by = doc_id]
    }

    df_tbl <- unique(dt[, .(doc_id, gene_key)])
    term_df <- df_tbl[, .N, by = gene_key]
    keep_terms <- term_df[N >= min_df, gene_key]
    dt <- dt[gene_key %in% keep_terms]
    dt[]
  }

  edges_all_tidy <- load_delta_links_all_tidy(step2_delta_csvs)
  edges_all_tidy <- edges_all_tidy |>
    dplyr::mutate(
      fp_gene_case = dplyr::coalesce(fp_bound_case, 0L) > 0L &
        dplyr::coalesce(gene_expr_flag_case, 0L) >= 1L,
      fp_gene_ctrl = dplyr::coalesce(fp_bound_ctrl, 0L) > 0L &
        dplyr::coalesce(gene_expr_flag_ctrl, 0L) >= 1L,
      tf_expr_either = dplyr::coalesce(tf_expr_flag_case, 0L) >= 1L |
        dplyr::coalesce(tf_expr_flag_ctrl, 0L) >= 1L,
      pass_topic_gate = tf_expr_either & (fp_gene_case | fp_gene_ctrl)
    ) |>
    dplyr::filter(pass_topic_gate)

  if (!("log2FC_tf_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2FC_tf_expr <- NA_real_
  if (!("log2FC_gene_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2FC_gene_expr <- NA_real_
  if (!("log2fc_tf_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2fc_tf_expr <- NA_real_
  if (!("log2fc_gene_expr" %in% names(edges_all_tidy))) edges_all_tidy$log2fc_gene_expr <- NA_real_
  edges_all_tidy <- edges_all_tidy |>
    dplyr::mutate(
      log2fc_tf_expr = dplyr::coalesce(.data$log2FC_tf_expr, .data$log2fc_tf_expr),
      log2fc_gene_expr = dplyr::coalesce(.data$log2FC_gene_expr, .data$log2fc_gene_expr),
      log2fc_fp_bed_score = log2((.data$fp_bed_score_case + log2fc_eps) /
                                   (.data$fp_bed_score_ctrl + log2fc_eps))
    )

  sum_cols_fp_rna <- intersect(
    c("delta_fp_bed_score", "log2fc_fp_bed_score"),
    names(edges_all_tidy)
  )
  first_cols_fp_rna <- intersect(
    c("comparison_id", "tf", "gene_key",
      "delta_tf_expr", "delta_gene_expr",
      "log2fc_tf_expr", "log2fc_gene_expr"),
    names(edges_all_tidy)
  )

  run_fp_rna_one <- function(run_cfg, edges_tbl = edges_all_tidy, run_suffix = NULL) {
    lda_res_all <- list()
    run_tag <- if (is.null(run_suffix)) run_cfg$tag else sprintf("%s_%s", run_cfg$tag, run_suffix)
    topic_benchmark_dir <- file.path(step2_out_dir, sprintf("topic_benchmark_%s", run_tag))
    doc_term_dir <- file.path(topic_benchmark_dir, "doc_terms")
    fit_dir <- file.path(topic_benchmark_dir, "fits")
    dir.create(topic_benchmark_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(doc_term_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

    edges_run <- if (isTRUE(run_cfg$aggregate_peaks)) {
      edges_agg <- aggregate_edges_by_tf_gene(
        edges_tbl,
        sum_cols = sum_cols_fp_rna,
        first_cols = first_cols_fp_rna
      )
      edges_agg$peak_id <- NA_character_
      edges_agg
    } else {
      edges_tbl
    }

    variants <- list(
      list(label = "diff_delta_fp",
           edges = prep_fp_only(edges_run, "delta_gene_expr", "delta_fp_bed_score")),
      list(label = "diff_fc_fp",
           edges = prep_fp_only(edges_run, "log2fc_gene_expr", "log2fc_fp_bed_score")),
      list(label = "diff_delta_fp_rna",
           edges = prep_fp_rna_edges(edges_run, "delta_gene_expr", "delta_fp_bed_score")),
      list(label = "diff_fc_fp_rna",
           edges = prep_fp_rna_edges(edges_run, "log2fc_gene_expr", "log2fc_fp_bed_score"))
    )

    for (v in variants) {
      if (is.null(v$edges) || !nrow(v$edges)) next

      docs <- build_doc_term_from_edges(
        v$edges,
        top_terms_per_doc = top_terms_per_doc,
        min_df = min_df
      )
      if (!nrow(docs)) next

      doc_term_path <- file.path(doc_term_dir, sprintf("doc_term_%s.rds", v$label))
      saveRDS(docs, doc_term_path)

      lda_grid <- expand.grid(
        gamma = gamma_cuts,
        K = lda_K,
        method = "VEM",
        stringsAsFactors = FALSE
      )
      lda_list <- split(lda_grid, seq_len(nrow(lda_grid)))

      run_one_lda <- function(row) {
        run_id <- sprintf("lda_%s_K%d_%s_g%d",
                          v$label, row$K, row$method, as.integer(row$gamma * 100))
        edge_topics_path <- file.path(fit_dir, sprintf("%s_edge_topics.csv", run_id))
        topics_path <- file.path(topic_benchmark_dir, sprintf("%s_topics.csv", run_id))
        by_group_stats_path <- file.path(topic_benchmark_dir, sprintf("%s_topics_by_comparison.csv", run_id))
        by_group_best_path <- file.path(topic_benchmark_dir, sprintf("%s_best_by_comparison.csv", run_id))

        if (!isTRUE(rerun_benchmark) &&
            isTRUE(reuse_edge_topics) &&
            file.exists(edge_topics_path) &&
            file.exists(topics_path)) {
          return(NULL)
        }

        lda_fit <- fit_pooled_lda(
          docs,
          K = row$K,
          method = row$method,
          seed = 1,
          gamma_cutoff = row$gamma
        )

        edge_topics <- assign_edge_topics_lda(
          edges_tbl = v$edges,
          lda_fit = lda_fit,
          top_n = 1L,
          aggregate_peak = run_cfg$aggregate_peaks
        )
        if (!nrow(edge_topics)) return(NULL)

        edge_topics <- edge_topics |>
          dplyr::mutate(
            gene_key_full = .data$gene_key,
            gene_key_base = sub("\\|[^|]+$", "", .data$gene_key),
            signal_type = dplyr::if_else(
              grepl("\\|", .data$gene_key),
              sub("^.*\\|", "", .data$gene_key),
              NA_character_
            )
          )
        if (all(is.na(edge_topics$signal_type)) && grepl("_fp", v$label)) {
          edge_topics$signal_type <- "FP"
        }
        edge_topics_out <- edge_topics |>
          dplyr::mutate(
            gene_key = .data$gene_key_base,
            edge_id = paste(.data$comparison_id,
                            .data$tf,
                            .data$gene_key_base,
                            dplyr::coalesce(.data$signal_type, "GENE"),
                            sep = "|")
          )

        readr::write_csv(edge_topics_out, edge_topics_path)

        edge_topics_mark <- edge_topics_out |>
          dplyr::mutate(gene_key = .data$gene_key_base)

        by_group <- topic_marker_overlap_from_edge_topics_by_group(
          edge_topics = edge_topics_mark,
          group_col = "comparison_id",
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )
        if (nrow(by_group$stats)) readr::write_csv(by_group$stats, by_group_stats_path)
        if (nrow(by_group$best)) readr::write_csv(by_group$best, by_group_best_path)

        out <- topic_marker_overlap_from_edge_topics(
          edge_topics = edge_topics_mark,
          markers_epi = markers_epi,
          markers_mes = markers_mes,
          universe_genes = unique(edge_topics_mark$gene_key),
          top_n_overlap_genes = 50,
          marker_mode = marker_mode,
          exclusive_method = exclusive_method,
          exclusive_top_n = exclusive_top_n,
          compute_gene_only = compute_gene_only,
          verbose = FALSE
        )

        stats_tbl <- out$stats |>
          dplyr::mutate(
            purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_epi_n / n_entities_unique, NA_real_),
            purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                                overlap_mes_n / n_entities_unique, NA_real_),
            purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE),
            label = v$label,
            method = row$method,
            K = row$K,
            gamma_cutoff = row$gamma
          )

        readr::write_csv(stats_tbl, topics_path)

        best_epi <- out$best_epi
        best_mes <- out$best_mes
        tibble::tibble(
          label = v$label,
          method = row$method,
          K = row$K,
          gamma_cutoff = row$gamma,
          best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
          best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
          best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_genes[[1]]) && best_epi$n_genes[[1]] > 0) {
            best_epi$overlap_epi_n[[1]] / best_epi$n_genes[[1]]
          } else {
            NA_real_
          },
          best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
          best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
          best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_genes[[1]]) && best_mes$n_genes[[1]] > 0) {
            best_mes$overlap_mes_n[[1]] / best_mes$n_genes[[1]]
          } else {
            NA_real_
          },
          topics_csv = topics_path,
          edge_topics_csv = edge_topics_path,
          topics_by_group_csv = by_group_stats_path,
          best_by_group_csv = by_group_best_path
        )
      }

      lda_res <- if (.Platform$OS.type != "windows" && bench_cores_run > 1L) {
        parallel::mclapply(lda_list, run_one_lda, mc.cores = bench_cores_run, mc.preschedule = FALSE)
      } else {
        lapply(lda_list, run_one_lda)
      }
      lda_res_all <- c(lda_res_all, Filter(Negate(is.null), lda_res))
      rm(docs)
      gc()
    }

    lda_sum <- dplyr::bind_rows(lda_res_all)
    if (nrow(lda_sum)) {
      readr::write_csv(lda_sum, file.path(topic_benchmark_dir, "lda_all_runs_summary.csv"))
    }

    all_topics_tbl <- summarize_topic_benchmark_edge_topics(
      topic_benchmark_dir = topic_benchmark_dir,
      fit_dir = fit_dir,
      markers_epi = markers_epi,
      markers_mes = markers_mes,
      out_file = file.path(topic_benchmark_dir, "topic_benchmark_all_topics.csv"),
      top_n_overlap_genes = 50L,
      marker_mode = marker_mode,
      exclusive_method = exclusive_method,
      exclusive_top_n = exclusive_top_n,
      compute_gene_only = compute_gene_only
    )

    write_topic_scatter_sig <- function(tbl, out_dir) {
      if (!nrow(tbl)) return(invisible(NULL))
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Plotting requires ggplot2.")
      }

      plot_tbl <- tbl |>
        dplyr::filter(
          is.finite(.data$epi_marker_coverage),
          is.finite(.data$mes_marker_coverage),
          is.finite(.data$purity_epi),
          is.finite(.data$purity_mes)
        ) |>
        dplyr::mutate(
          label_pretty = gsub("_", " ", gsub("^diff_", "diff ", gsub("^cond_", "cond ", .data$label))),
          K_plot_label = ifelse(is.na(.data$K), "K = NA", paste0("K = ", .data$K))
        )
      if (!nrow(plot_tbl)) return(invisible(NULL))

      k_levels <- sort(unique(plot_tbl$K[is.finite(plot_tbl$K)]))
      k_labels <- paste0("K = ", k_levels)
      if (any(is.na(plot_tbl$K))) k_labels <- c(k_labels, "K = NA")
      plot_tbl$K_plot_label <- factor(plot_tbl$K_plot_label, levels = k_labels)
      plot_tbl$label_pretty <- factor(plot_tbl$label_pretty, levels = unique(plot_tbl$label_pretty))

      base_theme <- ggplot2::theme_bw() +
        ggplot2::theme(
          text = ggplot2::element_text(face = "bold"),
          axis.text = ggplot2::element_text(face = "bold"),
          axis.title = ggplot2::element_text(face = "bold"),
          strip.text = ggplot2::element_text(face = "bold"),
          strip.background = ggplot2::element_rect(fill = "grey90", color = "grey40"),
          strip.text.y = ggplot2::element_text(angle = 0)
        )

      max_epi <- max(plot_tbl$epi_marker_coverage, na.rm = TRUE)
      max_mes <- max(plot_tbl$mes_marker_coverage, na.rm = TRUE)
      if (!is.finite(max_epi) || max_epi < 0) max_epi <- 0
      if (!is.finite(max_mes) || max_mes < 0) max_mes <- 0

      p_epi <- ggplot2::ggplot(
        plot_tbl,
        ggplot2::aes(x = .data$epi_marker_coverage, y = .data$purity_epi)
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.0) +
        ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
        ggplot2::scale_x_continuous(
          limits = c(0, max_epi + 0.01),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          x = "EPI marker coverage",
          y = "Epithelial purity"
        ) +
        base_theme

      p_mes <- ggplot2::ggplot(
        plot_tbl,
        ggplot2::aes(x = .data$mes_marker_coverage, y = .data$purity_mes)
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.0) +
        ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
        ggplot2::scale_x_continuous(
          limits = c(0, max_mes + 0.01),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          x = "MES marker coverage",
          y = "Mesenchymal purity"
        ) +
        base_theme

      ggplot2::ggsave(
        filename = file.path(out_dir, "topic_scatter_epi.pdf"),
        plot = p_epi,
        width = 16,
        height = 9,
        units = "in"
      )
      ggplot2::ggsave(
        filename = file.path(out_dir, "topic_scatter_mes.pdf"),
        plot = p_mes,
        width = 16,
        height = 9,
        units = "in"
      )
    }

    write_topic_scatter_sig(all_topics_tbl, topic_benchmark_dir)
  }

  edges_all_tidy <- edges_all_tidy |>
    dplyr::mutate(cell_line = sub("_.*$", "", .data$comparison_id))
  edges_by_cell <- split(edges_all_tidy, edges_all_tidy$cell_line)
  for (cell_line in names(edges_by_cell)) {
    edges_cell <- edges_by_cell[[cell_line]]
    if (!nrow(edges_cell)) next
    if (isTRUE(run_cfg_parallel) && .Platform$OS.type != "windows" && run_cfg_workers > 1L) {
      parallel::mclapply(
        topic_runs_fp_rna,
        function(cfg) run_fp_rna_one(cfg, edges_tbl = edges_cell, run_suffix = cell_line),
        mc.cores = run_cfg_workers,
        mc.preschedule = FALSE
      )
    } else {
      lapply(
        topic_runs_fp_rna,
        function(cfg) run_fp_rna_one(cfg, edges_tbl = edges_cell, run_suffix = cell_line)
      )
    }
  }
}

# Merge existing topic benchmark results (overlap + hclust)
do_topic_benchmark_merge <- TRUE
source("R/utils_grn_lda_nmf.R")
if (isTRUE(do_topic_benchmark_merge)) {
  source("R/utils_grn_lda_nmf.R")
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Topic merge requires data.table.")
  }

  merge_overlap_thresh <- 0.6
  merge_hclust_overlap_thresh <- merge_overlap_thresh
  merge_hclust_cut_height <- 1 - merge_hclust_overlap_thresh
  merge_dry_run <- FALSE
  merge_force_rewrite_topics <- TRUE

  topic_benchmark_merge_dirs <- list.dirs(step2_out_dir, full.names = TRUE, recursive = FALSE)
  topic_benchmark_merge_dirs <- topic_benchmark_merge_dirs[
    grepl("^topic_benchmark_", basename(topic_benchmark_merge_dirs)) &
      !grepl("^merged_", basename(topic_benchmark_merge_dirs))
  ]
  merge_allowlist <- c(
    "topic_benchmark_per_link_aggregated_sig_comparison_gibbs_Panc1",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_gibbs_Panc1",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_vem_Panc1",
    "topic_benchmark_per_link_aggregated_sig_comparison_vem_Panc1",
    "topic_benchmark_per_link_aggregated_sig_comparison_gibbs_HPAFII",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_gibbs_HPAFII",
    "topic_benchmark_per_link_aggregated_sig_comparison_gibbs",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_gibbs",
    "topic_benchmark_per_link_aggregated_sig_comparison_vem_HPAFII",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_vem_HPAFII",
    "topic_benchmark_per_link_aggregated_sig_comparison_gibbs_AsPC1",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_gibbs_AsPC1",
    "topic_benchmark_per_link_aggregated_sig_comparison_vem_AsPC1",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_vem_AsPC1",
    "topic_benchmark_per_link_aggregated_sig_comparison_vem",
    "topic_benchmark_per_link_nonaggregated_sig_comparison_vem",
    "topic_benchmark_per_link_nonaggregated_sig_comparison",
    "topic_benchmark_per_link_aggregated_sig_comparison"
  )
  topic_benchmark_merge_dirs <- topic_benchmark_merge_dirs[
    basename(topic_benchmark_merge_dirs) %in% merge_allowlist
  ]
  if (isTRUE(merge_dry_run)) {
    message("[merge dry run] folders to process:")
    if (length(topic_benchmark_merge_dirs)) {
      message(paste0(" - ", topic_benchmark_merge_dirs, collapse = "\n"))
    } else {
      message(" (none)")
    }
  }
  dry_run_shown <- FALSE

  markers_dir <- file.path(base_dir, "plot_lineage_plasticity_related_subnetworks", "markers")
  markers_epi <- readr::read_tsv(
    file.path(markers_dir, "epithelial.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  markers_mes <- readr::read_tsv(
    file.path(markers_dir, "mesenchymal.markers_nutrient_stress_all_lines.txt"),
    show_col_types = FALSE
  ) |>
    dplyr::pull(HGNC)
  if (!length(markers_epi) || !length(markers_mes)) {
    stop("markers_epi/markers_mes are empty; check marker files.")
  }

  merge_edge_topics_hclust <- function(edge_topics, overlap_thresh = 0.7) {
    dt <- data.table::as.data.table(edge_topics)
    if (!all(c("topic", "gene_key") %in% names(dt))) {
      stop("edge_topics must contain topic and gene_key for hclust merge.")
    }
    dt[, topic_raw := topic]
    dt[, topic_int := suppressWarnings(as.integer(topic))]
    use_int <- !all(is.na(dt$topic_int))
    if (use_int) {
      dt[, topic := topic_int]
    } else {
      dt[, topic := as.character(topic_raw)]
    }
    dt <- dt[!is.na(topic) & !is.na(gene_key)]
    topic_gene <- unique(dt[, .(topic, gene_key)])
    topics <- sort(unique(topic_gene$topic))
    if (length(topics) <= 1L) {
      mapping <- data.table::data.table(topic = topics, topic_merged = topics)
      dt[, topic_merged := topic]
      return(list(
        edge_topics = tibble::as_tibble(dt),
        mapping = tibble::as_tibble(mapping),
        topic_genes = tibble::as_tibble(topic_gene)
      ))
    }
    gene_sets <- split(topic_gene$gene_key, topic_gene$topic)
    n <- length(topics)
    dist_mat <- matrix(0, n, n, dimnames = list(topics, topics))
    for (i in seq_len(n - 1L)) {
      genes_i <- gene_sets[[as.character(topics[i])]]
      for (j in (i + 1L):n) {
        genes_j <- gene_sets[[as.character(topics[j])]]
        inter <- length(intersect(genes_i, genes_j))
        denom1 <- length(genes_i)
        denom2 <- length(genes_j)
        ratio1 <- if (denom1 > 0L) inter / denom1 else 0
        ratio2 <- if (denom2 > 0L) inter / denom2 else 0
        sim <- max(ratio1, ratio2)
        dist_mat[i, j] <- 1 - sim
        dist_mat[j, i] <- dist_mat[i, j]
      }
    }
    hc <- stats::hclust(stats::as.dist(dist_mat), method = "average")
    groups <- stats::cutree(hc, h = merge_hclust_cut_height)
    groups <- as.integer(groups)
    group_ids <- match(groups, sort(unique(groups)))
    mapping <- data.table::data.table(topic = topics, topic_merged = group_ids)
    if (use_int) {
      mapping[, topic := suppressWarnings(as.integer(topic))]
      mapping[, topic_merged := suppressWarnings(as.integer(topic_merged))]
    } else {
      mapping[, topic := as.character(topic)]
      mapping[, topic_merged := as.character(topic_merged)]
    }
    dt[, topic_raw := topic]
    dt[mapping, topic_merged := i.topic_merged, on = "topic"]
    dt[is.na(topic_merged), topic_merged := topic_raw]
    dt[, topic := topic_merged]
    topic_genes_merged <- dt[, .(gene_key = unique(gene_key)), by = topic]
    list(
      edge_topics = tibble::as_tibble(dt),
      mapping = tibble::as_tibble(mapping),
      topic_genes = tibble::as_tibble(topic_genes_merged)
    )
  }

  write_topic_scatter_merge <- function(tbl, out_dir) {
    if (!nrow(tbl)) return(invisible(NULL))
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Plotting requires ggplot2.")
    }
    plot_tbl <- tbl |>
      dplyr::filter(
        is.finite(.data$epi_marker_coverage),
        is.finite(.data$mes_marker_coverage),
        is.finite(.data$purity_epi),
        is.finite(.data$purity_mes)
      ) |>
      dplyr::mutate(
        label_pretty = gsub("_", " ", gsub("^diff_", "diff ", gsub("^cond_", "cond ", .data$label))),
        K_plot_label = ifelse(is.na(.data$K), "K = NA", paste0("K = ", .data$K))
      )
    if (!nrow(plot_tbl)) return(invisible(NULL))

    k_levels <- sort(unique(plot_tbl$K[is.finite(plot_tbl$K)]))
    k_labels <- paste0("K = ", k_levels)
    if (any(is.na(plot_tbl$K))) k_labels <- c(k_labels, "K = NA")
    plot_tbl$K_plot_label <- factor(plot_tbl$K_plot_label, levels = k_labels)
    plot_tbl$label_pretty <- factor(plot_tbl$label_pretty, levels = unique(plot_tbl$label_pretty))

    base_theme <- ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "grey40"),
        strip.text.y = ggplot2::element_text(angle = 0)
      )

    max_epi <- max(plot_tbl$epi_marker_coverage, na.rm = TRUE)
    max_mes <- max(plot_tbl$mes_marker_coverage, na.rm = TRUE)
    if (!is.finite(max_epi) || max_epi < 0) max_epi <- 0
    if (!is.finite(max_mes) || max_mes < 0) max_mes <- 0

    p_epi <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes(x = .data$epi_marker_coverage, y = .data$purity_epi)
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 1.0) +
      ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
      ggplot2::scale_x_continuous(
        limits = c(0, max_epi + 0.01),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = "EPI marker coverage",
        y = "Epithelial purity"
      ) +
      base_theme

    p_mes <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes(x = .data$mes_marker_coverage, y = .data$purity_mes)
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 1.0) +
      ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
      ggplot2::scale_x_continuous(
        limits = c(0, max_mes + 0.01),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = "MES marker coverage",
        y = "Mesenchymal purity"
      ) +
      base_theme

    ggplot2::ggsave(
      filename = file.path(out_dir, "topic_scatter_epi.pdf"),
      plot = p_epi,
      width = 16,
      height = 9,
      units = "in"
    )
    ggplot2::ggsave(
      filename = file.path(out_dir, "topic_scatter_mes.pdf"),
      plot = p_mes,
      width = 16,
      height = 9,
      units = "in"
    )
  }

  write_merged_run_topics <- function(stats_tbl, merge_dir) {
    if (!nrow(stats_tbl)) return(invisible(NULL))
    run_ids <- unique(stats_tbl$run_id)
    for (run_id in run_ids) {
      run_tbl <- stats_tbl |>
        dplyr::filter(.data$run_id == run_id)
      if (!nrow(run_tbl)) next

      map_file <- file.path(merge_dir, sprintf("%s_topic_mapping.csv", run_id))
      if (file.exists(map_file)) {
        map_tbl <- readr::read_csv(map_file, show_col_types = FALSE) |>
          dplyr::filter(!is.na(.data$topic), !is.na(.data$topic_merged))
        if (nrow(map_tbl)) {
          map_tbl <- map_tbl |>
            dplyr::mutate(
              topic_chr = as.character(.data$topic),
              topic_merged_chr = as.character(.data$topic_merged)
            )
          map_tbl <- if ("topic_merged_label" %in% names(map_tbl)) {
            map_tbl |>
              dplyr::mutate(topic_merged_label = as.character(.data$topic_merged_label))
          } else {
            map_tbl |>
              dplyr::group_by(.data$topic_merged_chr) |>
              dplyr::summarise(
                topic_merged_label = paste(sort(unique(.data$topic_chr)), collapse = ","),
                .groups = "drop"
              ) |>
              dplyr::right_join(
                dplyr::distinct(map_tbl, .data$topic_chr, .data$topic_merged_chr),
                by = "topic_merged_chr"
              )
          }

          label_tbl <- map_tbl |>
            dplyr::distinct(.data$topic_chr, .data$topic_merged_label)

          run_tbl <- run_tbl |>
            dplyr::mutate(topic_chr = as.character(.data$topic)) |>
            dplyr::left_join(label_tbl, by = "topic_chr") |>
            dplyr::mutate(topic = .data$topic_merged_label) |>
            dplyr::select(-dplyr::any_of(c("topic_chr", "topic_merged_label")))
        }
      }

      run_tbl <- run_tbl |>
        dplyr::mutate(topic = as.character(.data$topic)) |>
        dplyr::filter(!is.na(.data$topic), .data$topic != "") |>
        dplyr::group_by(.data$topic) |>
        dplyr::summarise(dplyr::across(dplyr::everything(), ~ dplyr::first(.x)), .groups = "drop")

      out_path <- file.path(merge_dir, sprintf("%s_topics.csv", run_id))
      if (isTRUE(merge_force_rewrite_topics) && file.exists(out_path)) {
        file.remove(out_path)
      }
      readr::write_csv(run_tbl, out_path)
      message(
        sprintf("[merge topics] %s rows=%d -> %s", run_id, nrow(run_tbl), out_path)
      )
    }
    invisible(TRUE)
  }

  merge_one_dir <- function(bench_dir) {
    fit_dir <- file.path(bench_dir, "fits")
    edge_files <- list.files(fit_dir, pattern = "_edge_topics\\.csv$", full.names = TRUE)
    if (!length(edge_files)) return(invisible(NULL))

    merge_dirs <- list(
      overlap = file.path(bench_dir, sprintf("merge_overlap%02d", as.integer(merge_overlap_thresh * 100))),
      hclust = file.path(bench_dir, sprintf("merge_hclust%02d", as.integer(merge_hclust_overlap_thresh * 100)))
    )
    merge_fit_dirs <- list(
      overlap = file.path(merge_dirs$overlap, "fits"),
      hclust = file.path(merge_dirs$hclust, "fits")
    )
    lapply(merge_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
    lapply(merge_fit_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

    for (f in edge_files) {
      base <- basename(f)
      edge_topics <- readr::read_csv(f, show_col_types = FALSE)
      if (!nrow(edge_topics)) next

      overlap_res <- merge_edge_topics_by_overlap(edge_topics, overlap_thresh = merge_overlap_thresh)
      readr::write_csv(overlap_res$edge_topics, file.path(merge_fit_dirs$overlap, base))
      overlap_map <- overlap_res$mapping |>
        dplyr::mutate(
          topic = as.integer(.data$topic),
          topic_merged = as.integer(.data$topic_merged)
        )
      overlap_label <- overlap_map |>
        dplyr::group_by(.data$topic_merged) |>
        dplyr::summarise(
          topic_merged_label = paste(sort(unique(.data$topic)), collapse = ","),
          .groups = "drop"
        )
      overlap_map <- overlap_map |>
        dplyr::left_join(overlap_label, by = "topic_merged")
      readr::write_csv(
        overlap_map,
        file.path(merge_dirs$overlap, sub("_edge_topics\\.csv$", "_topic_mapping.csv", base))
      )
      readr::write_csv(
        overlap_res$topic_genes,
        file.path(merge_dirs$overlap, sub("_edge_topics\\.csv$", "_topic_genes.csv", base))
      )

      hclust_res <- merge_edge_topics_hclust(edge_topics, overlap_thresh = merge_hclust_overlap_thresh)
      readr::write_csv(hclust_res$edge_topics, file.path(merge_fit_dirs$hclust, base))
      hclust_map <- hclust_res$mapping |>
        dplyr::mutate(
          topic = as.integer(.data$topic),
          topic_merged = as.integer(.data$topic_merged)
        )
      hclust_label <- hclust_map |>
        dplyr::group_by(.data$topic_merged) |>
        dplyr::summarise(
          topic_merged_label = paste(sort(unique(.data$topic)), collapse = ","),
          .groups = "drop"
        )
      hclust_map <- hclust_map |>
        dplyr::left_join(hclust_label, by = "topic_merged")
      readr::write_csv(
        hclust_map,
        file.path(merge_dirs$hclust, sub("_edge_topics\\.csv$", "_topic_mapping.csv", base))
      )
      readr::write_csv(
        hclust_res$topic_genes,
        file.path(merge_dirs$hclust, sub("_edge_topics\\.csv$", "_topic_genes.csv", base))
      )
    }

    merged_overlap_tbl <- summarize_topic_benchmark_edge_topics(
      topic_benchmark_dir = merge_dirs$overlap,
      fit_dir = merge_fit_dirs$overlap,
      markers_epi = markers_epi,
      markers_mes = markers_mes,
      out_file = file.path(merge_dirs$overlap, "topic_benchmark_all_topics.csv"),
      top_n_overlap_genes = 50L,
      marker_mode = "union",
      exclusive_method = "mean_score",
      exclusive_top_n = 1L,
      compute_gene_only = FALSE
    )
    write_topic_scatter_merge(merged_overlap_tbl, merge_dirs$overlap)
    write_merged_run_topics(merged_overlap_tbl, merge_dirs$overlap)

    merged_hclust_tbl <- summarize_topic_benchmark_edge_topics(
      topic_benchmark_dir = merge_dirs$hclust,
      fit_dir = merge_fit_dirs$hclust,
      markers_epi = markers_epi,
      markers_mes = markers_mes,
      out_file = file.path(merge_dirs$hclust, "topic_benchmark_all_topics.csv"),
      top_n_overlap_genes = 50L,
      marker_mode = "union",
      exclusive_method = "mean_score",
      exclusive_top_n = 1L,
      compute_gene_only = FALSE
    )
    write_topic_scatter_merge(merged_hclust_tbl, merge_dirs$hclust)
    write_merged_run_topics(merged_hclust_tbl, merge_dirs$hclust)
    invisible(TRUE)
  }

  if (isTRUE(merge_dry_run)) {
    for (bench_dir in topic_benchmark_merge_dirs) {
      fit_dir <- file.path(bench_dir, "fits")
      edge_files <- list.files(fit_dir, pattern = "_edge_topics\\.csv$", full.names = TRUE)
      if (!length(edge_files)) next

      merge_dirs <- list(
        overlap = file.path(bench_dir, sprintf("merge_overlap%02d", as.integer(merge_overlap_thresh * 100))),
        hclust = file.path(bench_dir, sprintf("merge_hclust%02d", as.integer(merge_hclust_overlap_thresh * 100)))
      )
      merge_fit_dirs <- list(
        overlap = file.path(merge_dirs$overlap, "fits"),
        hclust = file.path(merge_dirs$hclust, "fits")
      )
      if (!dry_run_shown) {
        dry_run_shown <- TRUE
        message("[merge dry run] example bench_dir: ", bench_dir)
        message("[merge dry run] overlap dir: ", merge_dirs$overlap)
        message("[merge dry run] hclust dir: ", merge_dirs$hclust)
        for (f in edge_files) {
          base <- basename(f)
          message("[merge dry run] edge_topics: ", f)
          message("  -> overlap merged: ", file.path(merge_fit_dirs$overlap, base))
          message("  -> overlap mapping: ",
                  file.path(merge_dirs$overlap, sub("_edge_topics\\.csv$", "_topic_mapping.csv", base)))
          message("  -> overlap genes: ",
                  file.path(merge_dirs$overlap, sub("_edge_topics\\.csv$", "_topic_genes.csv", base)))
          message("  -> hclust merged: ", file.path(merge_fit_dirs$hclust, base))
          message("  -> hclust mapping: ",
                  file.path(merge_dirs$hclust, sub("_edge_topics\\.csv$", "_topic_mapping.csv", base)))
          message("  -> hclust genes: ",
                  file.path(merge_dirs$hclust, sub("_edge_topics\\.csv$", "_topic_genes.csv", base)))
        }
        message("[merge dry run] overlap summary: ",
                file.path(merge_dirs$overlap, "topic_benchmark_all_topics.csv"))
        message("[merge dry run] overlap plots: ",
                file.path(merge_dirs$overlap, "topic_scatter_epi.pdf"),
                " ; ",
                file.path(merge_dirs$overlap, "topic_scatter_mes.pdf"))
        message("[merge dry run] hclust summary: ",
                file.path(merge_dirs$hclust, "topic_benchmark_all_topics.csv"))
        message("[merge dry run] hclust plots: ",
                file.path(merge_dirs$hclust, "topic_scatter_epi.pdf"),
                " ; ",
                file.path(merge_dirs$hclust, "topic_scatter_mes.pdf"))
      }
      break
    }
  } else {
    merge_use_parallel <- TRUE
    merge_workers <- parallel::detectCores(logical = TRUE)
    if (!is.numeric(merge_workers) || merge_workers < 1L) merge_workers <- 1L
    use_mc <- isTRUE(merge_use_parallel) && .Platform$OS.type != "windows" && merge_workers > 1L
    if (use_mc) {
      parallel::mclapply(
        topic_benchmark_merge_dirs,
        merge_one_dir,
        mc.cores = merge_workers,
        mc.preschedule = FALSE
      )
    } else {
      lapply(topic_benchmark_merge_dirs, merge_one_dir)
    }
  }
}

do_pathway_topic_benchmark <- TRUE
if (isTRUE(do_pathway_topic_benchmark)) {
  if (!exists("pathway_signatures") || !length(pathway_signatures)) {
    stop("pathway_signatures is empty; run do_pathway_signatures first.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Pathway benchmarking requires data.table.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Pathway benchmarking requires ggplot2.")
  }

  tmp_dir <- file.path(step2_out_dir, "tmp_r")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TMPDIR = tmp_dir, TEMP = tmp_dir, TMP = tmp_dir)

  pathway_top_n <- 20L
  pathway_min_overlap <- 1L
  pathway_sig_abs_log2fc_gene_expr <- 1
  pathway_merge_overlap_thresh <- 0.6
  pathway_merge_hclust_overlap_thresh <- 0.6
  pathway_scatter_per_pathway <- TRUE
  pathway_benchmark_log <- file.path(step2_out_dir, "pathway_benchmark_processed_dirs.txt")

  parse_pathway_comp <- function(comp) {
    parts <- strsplit(comp, "_", fixed = TRUE)[[1]]
    if (length(parts) < 3L) {
      return(list(stress_type = NA_character_, cell_line = NA_character_, stress = NA_character_))
    }
    list(
      stress_type = parts[[1]],
      cell_line = parts[[2]],
      stress = paste(parts[-c(1, 2)], collapse = "_")
    )
  }

  to_pathway_comp <- function(comp_id, known_comps) {
    comp_id <- as.character(comp_id)
    if (comp_id %in% known_comps) return(comp_id)
    case_part <- sub("_vs_.*$", "", comp_id)
    parts <- strsplit(case_part, "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2L) {
      cell <- parts[[1]]
      stress <- paste(parts[-1], collapse = "_")
      stress_parts <- strsplit(stress, "_", fixed = TRUE)[[1]]
      stress_type <- if (length(stress_parts)) tail(stress_parts, 1) else NA_character_
      key <- paste(stress_type, tolower(cell), stress, sep = "_")
      if (key %in% known_comps) return(key)
    }
    NA_character_
  }

  # Build per-comparison significant target-gene set from Step2 delta CSVs
  # (used to restrict topic genes before computing pathway coverage/purity).
  build_sig_gene_tbl_cellline <- function(step2_dir,
                                          cell_line,
                                          abs_log2fc_min = 1) {
    if (!is.character(step2_dir) || !nzchar(step2_dir) || !dir.exists(step2_dir)) return(NULL)
    if (!is.character(cell_line) || !nzchar(cell_line)) return(NULL)
    if (!is.numeric(abs_log2fc_min) || length(abs_log2fc_min) != 1L) return(NULL)

    pat <- sprintf("^%s_.+_vs_%s_10_FBS_delta_links\\.csv$", cell_line, cell_line)
    delta_files <- list.files(step2_dir, pattern = pat, full.names = TRUE)
    if (!length(delta_files)) return(NULL)

    pick_col <- function(nms, candidates) {
      hit <- intersect(candidates, nms)
      if (!length(hit)) return(NA_character_)
      hit[[1]]
    }

    out <- lapply(delta_files, function(f) {
      comp_id <- sub("_delta_links\\.csv$", "", basename(f))
      tbl <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
      if (is.null(tbl) || !nrow(tbl)) return(NULL)

      gene_col <- pick_col(names(tbl), c("gene_key", "gene", "target", "Name.Target"))
      l2fc_col <- pick_col(names(tbl), c("log2FC_gene_expr", "log2FC_gene", "log2FC_target", "log2fc_gene_expr"))
      if (is.na(gene_col) || is.na(l2fc_col)) return(NULL)

      tbl |>
        dplyr::select(
          gene_key = dplyr::all_of(gene_col),
          log2FC_gene_expr = dplyr::all_of(l2fc_col)
        ) |>
        dplyr::filter(is.finite(.data$log2FC_gene_expr)) |>
        dplyr::filter(abs(.data$log2FC_gene_expr) > abs_log2fc_min) |>
        dplyr::transmute(comparison_id = comp_id, gene_key = as.character(.data$gene_key)) |>
        dplyr::distinct()
    })

    out <- dplyr::bind_rows(Filter(Negate(is.null), out))
    if (!nrow(out)) return(NULL)
    out
  }

  build_pathway_df <- function(sig_list) {
    out <- list()
    idx <- 1L
    for (comp in names(sig_list)) {
      comp_entry <- sig_list[[comp]]
      if (!length(comp_entry)) next
      for (dir in names(comp_entry)) {
        paths <- comp_entry[[dir]]
        if (!length(paths)) next
        path_names <- names(paths)
        tbl <- data.frame(
          comparison_id = rep(comp, length(paths)),
          direction = rep(stringr::str_to_title(dir), length(paths)),
          pathway = path_names,
          adj_p = vapply(paths, function(x) x$adj_p, numeric(1)),
          overlap_n = vapply(paths, function(x) x$overlap_n, integer(1)),
          overlap_total = vapply(paths, function(x) x$overlap_total, integer(1)),
          db = vapply(paths, function(x) x$db, character(1)),
          stringsAsFactors = FALSE
        )
        tbl$genes <- lapply(paths, function(x) unique(as.character(x$genes)))
        tbl$pathway_size <- vapply(tbl$genes, length, integer(1))
        out[[idx]] <- tibble::as_tibble(tbl)
        idx <- idx + 1L
      }
    }
    dplyr::bind_rows(out)
  }

  pathway_df <- build_pathway_df(pathway_signatures)
  if (!nrow(pathway_df)) {
    stop("No pathway signature rows found after parsing.")
  }
  # Treat Up/Down equally by collapsing to a single pathway gene set per comparison.
  pathway_df <- pathway_df |>
    dplyr::group_by(.data$comparison_id, .data$pathway, .data$db) |>
    dplyr::summarise(
      direction = "All",
      adj_p = suppressWarnings(min(.data$adj_p, na.rm = TRUE)),
      overlap_n = suppressWarnings(max(.data$overlap_n, na.rm = TRUE)),
      overlap_total = suppressWarnings(max(.data$overlap_total, na.rm = TRUE)),
      genes = list(unique(unlist(.data$genes))),
      pathway_size = length(unique(unlist(.data$genes))),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      adj_p = ifelse(is.finite(.data$adj_p), .data$adj_p, NA_real_)
    ) |>
    # Keep only significant pathways, then cap to top 10 by adjusted p per comparison.
    dplyr::filter(is.finite(.data$adj_p), .data$adj_p < 0.01) |>
    dplyr::group_by(.data$comparison_id) |>
    dplyr::slice_min(order_by = .data$adj_p, n = 10, with_ties = FALSE) |>
    dplyr::ungroup()
  parsed_pathway_comp <- lapply(pathway_df$comparison_id, parse_pathway_comp)
  pathway_df$stress_type <- vapply(parsed_pathway_comp, function(x) x$stress_type, character(1))
  pathway_df$cell_line <- vapply(parsed_pathway_comp, function(x) x$cell_line, character(1))
  pathway_df$stress <- vapply(parsed_pathway_comp, function(x) x$stress, character(1))

  compute_pathway_overlap <- function(edge_topics,
                                      pathway_df_use,
                                      known_comps,
                                      top_n = 20L,
                                      min_overlap = 1L) {
    dt <- data.table::as.data.table(edge_topics)
    req <- c("comparison_id", "topic", "gene_key")
    miss <- setdiff(req, names(dt))
    if (length(miss)) return(NULL)

    dt <- dt[!is.na(comparison_id) & !is.na(topic) & !is.na(gene_key)]
    if (!nrow(dt)) return(NULL)
    dt[, comparison_id := as.character(comparison_id)]
    dt[, topic := as.integer(topic)]
    dt[, gene_key := as.character(gene_key)]

    comp_map <- data.table::data.table(comparison_id = unique(dt$comparison_id))
    comp_map[, pathway_comp := vapply(comparison_id, to_pathway_comp, character(1), known_comps = known_comps)]
    dt <- dt[comp_map, on = "comparison_id"]
    dt <- dt[!is.na(pathway_comp)]
    if (!nrow(dt)) return(NULL)

    topic_gene <- unique(dt[, .(pathway_comp, comparison_id, topic, gene_key)])
    topic_sizes <- topic_gene[, .(topic_size = data.table::uniqueN(gene_key)),
                              by = .(pathway_comp, comparison_id, topic)]

    path_dt <- data.table::as.data.table(pathway_df_use)
    data.table::setnames(path_dt, "comparison_id", "pathway_comp")
    data.table::setkey(path_dt, pathway_comp)

    comps <- intersect(unique(topic_gene$pathway_comp), unique(path_dt$pathway_comp))
    if (!length(comps)) return(NULL)

    res_list <- vector("list", length(comps))
    idx <- 1L
    for (comp in comps) {
      tg <- topic_gene[pathway_comp == comp]
      if (!nrow(tg)) next
      comp_paths <- path_dt[comp]
      if (!nrow(comp_paths)) next

      genes_by_topic <- split(tg$gene_key, tg$topic)
      size_map_tbl <- topic_sizes[pathway_comp == comp]
      size_map <- setNames(size_map_tbl$topic_size, size_map_tbl$topic)

      res_comp <- lapply(seq_len(nrow(comp_paths)), function(i) {
        genes <- comp_paths$genes[[i]]
        if (is.null(genes) || !length(genes)) return(NULL)
        genes <- unique(genes)
        res_topics <- lapply(names(genes_by_topic), function(t) {
          g_t <- genes_by_topic[[t]]
          ov <- length(intersect(genes, g_t))
          if (ov < min_overlap) return(NULL)
          cov <- ov / length(genes)
          t_size <- size_map[[t]]
          pur <- if (!is.null(t_size) && t_size > 0) ov / t_size else NA_real_
          data.table::data.table(
            comparison_id = tg$comparison_id[1],
            pathway_comp = comp,
            direction = comp_paths$direction[[i]],
            pathway = comp_paths$pathway[[i]],
            overlap_n = ov,
            pathway_size = length(genes),
            coverage = cov,
            purity = pur,
            topic = as.integer(t),
            topic_size = t_size,
            adj_p = comp_paths$adj_p[[i]],
            db = comp_paths$db[[i]]
          )
        })
        data.table::rbindlist(res_topics, use.names = TRUE, fill = TRUE)
      })

      res_comp <- data.table::rbindlist(res_comp, use.names = TRUE, fill = TRUE)
      if (!nrow(res_comp)) next
      if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
        best_paths <- res_comp[, .(
          max_cov = max(coverage, na.rm = TRUE),
          max_ov = max(overlap_n, na.rm = TRUE)
        ), by = .(pathway_comp, pathway)]
        data.table::setorder(best_paths, pathway_comp, -max_cov, -max_ov)
        keep_paths <- best_paths[, head(.SD, top_n), by = .(pathway_comp)]
        res_comp <- res_comp[keep_paths, on = .(pathway_comp, pathway), nomatch = 0]
      }
      res_list[[idx]] <- res_comp
      idx <- idx + 1L
    }

    data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  }

  parse_run_meta <- function(fname) {
    base <- sub("_edge_topics\\.csv$", "", basename(fname))
    if (grepl("^lda_", base)) {
      m <- regexec("^lda_(.+)_K(\\d+)_([^_]+)_g(\\d+)$", base)
      hit <- regmatches(base, m)[[1]]
      if (length(hit)) {
        return(list(label = hit[2], method = hit[4], K = as.integer(hit[3]),
                    gamma_cutoff = as.numeric(hit[5]) / 100, run_id = base))
      }
    }
    if (grepl("^nmf_", base)) {
      m <- regexec("^nmf_(.+)_K(\\d+)_([^_]+)_g(\\d+)$", base)
      hit <- regmatches(base, m)[[1]]
      if (length(hit)) {
        return(list(label = hit[2], method = hit[4], K = as.integer(hit[3]),
                    gamma_cutoff = as.numeric(hit[5]) / 100, run_id = base))
      }
    }
    if (grepl("^louvain_", base)) {
      label <- sub("^louvain_", "", base)
      return(list(label = label, method = "louvain", K = NA_integer_,
                  gamma_cutoff = NA_real_, run_id = base))
    }
    list(label = base, method = NA_character_, K = NA_integer_, gamma_cutoff = NA_real_, run_id = base)
  }

  pathway_stress_colors <- c(
    BCAA = "#3a78af",
    FBS = "#ef812f",
    Glc = "#308a4e",
    `Gln.Arg` = "#414b8c",
    Lys = "#e89c84",
    `Met.Cys` = "#ca2f2d",
    Trp = "#916ab6"
  )

  # Scatter summary across ALL pathways (Up + Down treated equally).
  write_pathway_scatter <- function(tbl, out_dir, per_pathway = FALSE) {
    plot_tbl <- tbl |>
      dplyr::filter(is.finite(.data$coverage), is.finite(.data$purity))
    if (!nrow(plot_tbl)) return(invisible(NULL))

    plot_tbl <- plot_tbl |>
      dplyr::mutate(
        label_pretty = gsub("_", " ", gsub("^diff_", "diff ", gsub("^cond_", "cond ", .data$label))),
        K_plot_label = ifelse(is.na(.data$K), "K = NA", paste0("K = ", .data$K))
      )
    if (!"stress_type" %in% names(plot_tbl)) {
      if ("stress" %in% names(plot_tbl)) {
        plot_tbl$stress_type <- plot_tbl$stress
      } else if ("comparison_id" %in% names(plot_tbl)) {
        parsed <- lapply(plot_tbl$comparison_id, parse_pathway_comp)
        plot_tbl$stress_type <- vapply(parsed, function(x) x$stress_type, character(1))
      } else if ("pathway_comp" %in% names(plot_tbl)) {
        parsed <- lapply(plot_tbl$pathway_comp, parse_pathway_comp)
        plot_tbl$stress_type <- vapply(parsed, function(x) x$stress_type, character(1))
      } else {
        plot_tbl$stress_type <- NA_character_
      }
    }
    plot_tbl$stress_type <- factor(plot_tbl$stress_type, levels = names(pathway_stress_colors))
    plot_tbl$label_pretty <- factor(plot_tbl$label_pretty, levels = unique(plot_tbl$label_pretty))
    k_levels <- sort(unique(plot_tbl$K[is.finite(plot_tbl$K)]))
    k_labels <- paste0("K = ", k_levels)
    if (any(is.na(plot_tbl$K))) k_labels <- c(k_labels, "K = NA")
    plot_tbl$K_plot_label <- factor(plot_tbl$K_plot_label, levels = k_labels)

    max_cov <- max(plot_tbl$coverage, na.rm = TRUE)
    if (!is.finite(max_cov) || max_cov < 0) max_cov <- 0

    base_theme <- ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "grey40"),
        strip.text.y = ggplot2::element_text(angle = 0)
      )

    p <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes(x = .data$coverage, y = .data$purity, color = .data$stress_type, shape = .data$stress_type)
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 1.0) +
      ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
      ggplot2::scale_x_continuous(
        limits = c(0, max_cov + 0.01),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_color_manual(values = pathway_stress_colors, na.value = "grey60") +
      ggplot2::scale_shape_manual(values = seq_along(pathway_stress_colors), na.value = 4) +
      ggplot2::labs(
        x = "Pathway coverage",
        y = "Topic purity",
        color = "Stress type",
        shape = "Stress type"
      ) +
      base_theme

    out_file <- file.path(out_dir, "pathway_scatter.pdf")
    old_up <- file.path(out_dir, "pathway_scatter_up.pdf")
    old_down <- file.path(out_dir, "pathway_scatter_down.pdf")
    if (file.exists(old_up)) file.remove(old_up)
    if (file.exists(old_down)) file.remove(old_down)
    ggplot2::ggsave(out_file, plot = p, width = 16, height = 9, units = "in")

    if (!isTRUE(per_pathway)) return(invisible(NULL))

    make_safe <- function(x) {
      x <- gsub("[^A-Za-z0-9_.-]", "_", x)
      substr(x, 1, 180)
    }
    per_dir <- file.path(out_dir, "pathway_scatter_per_pathway")
    dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)

    group_tbl <- plot_tbl |>
      dplyr::mutate(
        group_id = paste(.data$comparison_id, .data$pathway, sep = "__")
      )

    split_groups <- split(group_tbl, group_tbl$group_id)
    for (gid in names(split_groups)) {
      gtbl <- split_groups[[gid]]
      if (!nrow(gtbl)) next
      gtbl <- gtbl |>
        dplyr::mutate(
          label_pretty = factor(.data$label_pretty, levels = levels(plot_tbl$label_pretty)),
          K_plot_label = factor(.data$K_plot_label, levels = levels(plot_tbl$K_plot_label))
        )

      max_cov_g <- max(gtbl$coverage, na.rm = TRUE)
      if (!is.finite(max_cov_g) || max_cov_g < 0) max_cov_g <- 0

      pg <- ggplot2::ggplot(
        gtbl,
        ggplot2::aes(x = .data$coverage, y = .data$purity, color = .data$stress_type, shape = .data$stress_type)
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.0) +
        ggplot2::facet_grid(.data$K_plot_label ~ .data$label_pretty) +
        ggplot2::scale_x_continuous(
          limits = c(0, max_cov_g + 0.01),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_color_manual(values = pathway_stress_colors, na.value = "grey60") +
        ggplot2::scale_shape_manual(values = seq_along(pathway_stress_colors), na.value = 4) +
        ggplot2::labs(
          x = "Pathway coverage",
          y = "Topic purity",
          color = "Stress type",
          shape = "Stress type",
          title = gid
        ) +
        base_theme

      out_path <- file.path(per_dir, paste0(make_safe(gid), ".pdf"))
      ggplot2::ggsave(out_path, plot = pg, width = 16, height = 9, units = "in")
    }
  }

  log_pathway_dir <- function(log_file, dir_path, n_files) {
    if (!nzchar(log_file)) return(invisible(NULL))
    if (!file.exists(log_file)) {
      writeLines("timestamp\tdir_path\tn_edge_files", log_file)
    }
    entry <- sprintf("%s\t%s\t%d", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), dir_path, n_files)
    write(entry, file = log_file, append = TRUE)
  }

  run_pathway_benchmark_dir <- function(bench_dir, fit_dir, pathway_df_use, known_comps, log_file = "") {
    edge_files <- list.files(fit_dir, pattern = "_edge_topics\\.csv$", full.names = TRUE)
    if (!length(edge_files)) return(invisible(NULL))
    log_pathway_dir(log_file, fit_dir, length(edge_files))

    path_use_parallel <- TRUE
    path_workers <- parallel::detectCores(logical = TRUE)
    if (!is.numeric(path_workers) || path_workers < 1L) path_workers <- 1L
    use_mc <- isTRUE(path_use_parallel) && .Platform$OS.type != "windows" && path_workers > 1L

    worker <- function(f) {
      meta <- parse_run_meta(f)
      edge_topics <- readr::read_csv(f, show_col_types = FALSE)
      if (!nrow(edge_topics)) return(NULL)
      out <- compute_pathway_overlap(edge_topics, pathway_df_use,
                                     known_comps = known_comps,
                                     top_n = pathway_top_n,
                                     min_overlap = pathway_min_overlap)
      if (is.null(out) || !nrow(out)) return(NULL)
      out <- tibble::as_tibble(out) |>
        dplyr::mutate(
          label = meta$label,
          method = meta$method,
          K = meta$K,
          gamma_cutoff = meta$gamma_cutoff,
          run_id = meta$run_id
        )
      out
    }

    stats_list <- if (use_mc) {
      parallel::mclapply(edge_files, worker, mc.cores = path_workers, mc.preschedule = FALSE)
    } else {
      lapply(edge_files, worker)
    }

    stats_tbl <- dplyr::bind_rows(Filter(Negate(is.null), stats_list))
    if (!nrow(stats_tbl)) return(invisible(NULL))

    parsed_comp <- lapply(stats_tbl$pathway_comp, parse_pathway_comp)
    stats_tbl$stress_type <- vapply(parsed_comp, function(x) x$stress_type, character(1))
    stats_tbl$cell_line <- vapply(parsed_comp, function(x) x$cell_line, character(1))
    stats_tbl$stress <- vapply(parsed_comp, function(x) x$stress, character(1))

    readr::write_csv(stats_tbl, file.path(bench_dir, "pathway_benchmark_all_paths.csv"))
    write_pathway_scatter(stats_tbl, bench_dir, per_pathway = pathway_scatter_per_pathway)
  }

  pathway_benchmark_dirs <- list.dirs(step2_out_dir, recursive = FALSE, full.names = TRUE)
  pathway_benchmark_dirs <- pathway_benchmark_dirs[
    grepl("^topic_benchmark_per_link_.*_sig_comparison", basename(pathway_benchmark_dirs)) &
      grepl("_(AsPC1|HPAFII|Panc1)$", basename(pathway_benchmark_dirs))
  ]

  for (bench_dir in pathway_benchmark_dirs) {
    cell_line <- sub(".*_(AsPC1|HPAFII|Panc1)$", "\\1", basename(bench_dir))
    cell_line_key <- tolower(cell_line)
    pathway_df_use <- pathway_df |>
      dplyr::filter(.data$cell_line == cell_line_key)
    if (!nrow(pathway_df_use)) next
    known_comps <- unique(pathway_df_use$comparison_id)

    fit_dir <- file.path(bench_dir, "fits")
    if (dir.exists(fit_dir)) {
      run_pathway_benchmark_dir(bench_dir, fit_dir, pathway_df_use, known_comps, log_file = pathway_benchmark_log)
    }

    overlap_dir <- file.path(bench_dir, sprintf("merge_overlap%02d", as.integer(pathway_merge_overlap_thresh * 100)))
    overlap_fit <- file.path(overlap_dir, "fits")
    if (dir.exists(overlap_fit)) {
      run_pathway_benchmark_dir(overlap_dir, overlap_fit, pathway_df_use, known_comps, log_file = pathway_benchmark_log)
    }

    hclust_dir <- file.path(bench_dir, sprintf("merge_hclust%02d", as.integer(pathway_merge_hclust_overlap_thresh * 100)))
    hclust_fit <- file.path(hclust_dir, "fits")
    if (dir.exists(hclust_fit)) {
      run_pathway_benchmark_dir(hclust_dir, hclust_fit, pathway_df_use, known_comps, log_file = pathway_benchmark_log)
    }
  }

  do_pathway_topic_benchmark_single <- FALSE
  if (isTRUE(do_pathway_topic_benchmark_single)) {
    pathway_single_comp <- "AsPC1_0.4_FBS_vs_AsPC1_10_FBS"
    pathway_single_pathway <- "ATF4 human tf ARCHS4 coexpression"
    pathway_single_cell_line <- "AsPC1"
    pathway_single_edge_file <- file.path(
      step2_out_dir,
      "topic_benchmark_per_link_aggregated_sig_comparison_vem_AsPC1",
      "fits",
      "lda_diff_delta_combo_multi_K20_VEM_g70_edge_topics.csv"
    )
    if (!file.exists(pathway_single_edge_file)) {
      stop("Single-pathway test edge file not found: ", pathway_single_edge_file)
    }
    bench_dir <- dirname(pathway_single_edge_file) |> dirname()
    pathway_df_use <- pathway_df |>
      dplyr::filter(.data$cell_line == tolower(pathway_single_cell_line))
    known_comps <- unique(pathway_df_use$comparison_id)
    meta <- parse_run_meta(pathway_single_edge_file)
    edge_topics <- readr::read_csv(pathway_single_edge_file, show_col_types = FALSE)
    stats_tbl <- compute_pathway_overlap(edge_topics, pathway_df_use,
                                         known_comps = known_comps,
                                         top_n = pathway_top_n,
                                         min_overlap = pathway_min_overlap)
    if (!is.null(stats_tbl) && nrow(stats_tbl)) {
      stats_tbl <- tibble::as_tibble(stats_tbl) |>
        dplyr::mutate(
          label = meta$label,
          method = meta$method,
          K = meta$K,
          gamma_cutoff = meta$gamma_cutoff,
          run_id = meta$run_id
        ) |>
        dplyr::filter(.data$comparison_id == pathway_single_comp,
                      .data$pathway == pathway_single_pathway)

      single_out_dir <- file.path(bench_dir, "pathway_scatter_single")
      dir.create(single_out_dir, recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(stats_tbl, file.path(single_out_dir, "pathway_benchmark_all_paths_single.csv"))
      write_pathway_scatter(stats_tbl, single_out_dir, per_pathway = TRUE)
    }
  }
}

# Run and assign back to grn_set
grn_set$fp_annotation_spearman <- make_fp_annotation_corr(grn_set, method = "spearman", cores = 25L)
grn_set$fp_annotation_kendall <- make_fp_annotation_corr(grn_set, method = "kendall",  cores = 25L)
grn_set$fp_annotation_pearson <- grn_set$fp_annotation


# match_mode <- "strict"  # strict lenient
regulated_genes <- 1.5 # 1.5 2
delta_link <- 1 # 2 1
regulation_priors <- "genehancer" # "genehancer" 300kb
run_lighting_pipeline_one <- function(
  fp_gene_corr_kept,
  lighting_folder,
  grn_set,
  regulated_genes,
  delta_link,
  tf_p_thr,
  tf_r_min,
  tf_p_col = "p_adj_tf"
) {
  options(future.globals.maxSize = max(getOption("future.globals.maxSize", 500 * 1024^2), 32 * 1024^3))
  dir.create(lighting_folder, recursive = TRUE, showWarnings = FALSE)

  # Step 3. Build basal GRN & identify active regulatory edges per condition
  basal <- make_basal_links(
    fp_gene_corr_kept = fp_gene_corr_kept,
    fp_annotation     = grn_set$fp_annotation,
    out_dir           = lighting_folder,
    prefix            = "lighting",
    fp_variance       = grn_set$fp_variance,
    rna_variance      = grn_set$rna_variance
  )

  if (!all(c("r_tf", tf_p_col) %in% names(basal))) {
    stop(sprintf("Basal TF filtering requires columns r_tf and %s", tf_p_col))
  }
  basal <- basal |>
    dplyr::filter(!is.na(.data$r_tf), !is.na(.data[[tf_p_col]])) |>
    dplyr::filter(.data$r_tf >= tf_r_min) |>
    dplyr::filter(.data[[tf_p_col]] <= tf_p_thr)

  light_by_condition(
    ds = grn_set,
    basal_links = basal,
    out_dir = lighting_folder,
    prefix = "lighting",
    label_col = "strict_match_rna",
    link_score_threshold = link_score_threshold,
    fp_score_threshold = fp_score_threshold,
    tf_expr_threshold = threshold_tf_expr,
    fp_bound_tbl = grn_set$fp_bound,
    rna_expressed_tbl = grn_set$rna_expressed,
    require_fp_bound = TRUE,
    require_gene_expr = TRUE,
    gene_expr_threshold = 1L,
    filter_active = FALSE,
    use_parallel = TRUE,
    workers = 8,
    fp_variance_tbl = grn_set$fp_variance,
    rna_variance_tbl = grn_set$rna_variance
  )

  # Step 4. Perform differential GRN analysis & identify master TFs
  specs <- build_cellwise_contrasts_from_index(
    index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
    out_dir = lighting_folder,
    prefix = "lighting",
    ctrl_tag = "10_FBS",
    clean_names = FALSE
  )

  run_links_deltas_driver(
    specs       = specs,
    clean_names = FALSE,
    parallel    = TRUE,
    restrict_to_active_both = TRUE,
    edge_change_min = delta_link,
    keep_all_cols = TRUE
  )

  delta_csvs <- list.files(lighting_folder, "_delta_links.csv", full.names = TRUE)
  de_gene_log2_abs_min <- if (regulated_genes == 1.5) 0.585 else if (regulated_genes == 2) 1 else NA_real_

  bulk <- episcope::filter_links_deltas_bulk(
    delta_csvs,
    gene_expr_min = threshold_gene_expr,
    tf_expr_min = threshold_tf_expr,
    fp_min = threshold_fp_score,
    link_min = threshold_link_score,
    abs_delta_min = delta_link,
    apply_de_gene = TRUE,
    de_gene_log2_abs_min = de_gene_log2_abs_min,
    enforce_link_expr_sign = TRUE,
    expr_dir_col = "log2FC_gene_expr",
    workers = 20
  )

  # LDA across all filtered CSVs produced above
  # filtered_csvs <- bulk$filtered_paths

  # lda_out <- episcope::annotate_links_with_lda_topic_bulk(
  #   filtered_csvs,
  #   K = 20, # topics
  #   which = "abs", # |delta_link_score|
  #   min_df = 2, # drop ultra-rare terms
  #   gamma_cutoff = 0.2, # keep all topics with gamma >= 0.2 per TF
  #   parallel = TRUE, workers = 20, plan = "multisession", seed = 1
  # )

  # Summarize all annotated CSVs produced by LDA
  # annotated_csvs <- lda_out$assigned_paths

  # sum_out <- episcope::summarize_lda_annotations_bulk(
  #   annotated_csvs,
  #   edge_filter_min = delta_link, # keep edges used to call delta links
  #   parallel = TRUE, plan = "multisession", workers = 20,
  #   use_tf_r_weight = TRUE # optional: weight deltas by TF r column if present
  # )

  # Generate all bubble & hub PDFs
  # summary_csvs <- sum_out$summary_paths

  # episcope::plot_from_summary_bulk(
  #   summary_csvs,
  #   top_tf_per_topic = 20,
  #   parallel = TRUE,
  #   workers = 20,
  #   plan = "multisession",
  #   color_sigma = 2
  # )

  invisible(list(
    lighting_folder = lighting_folder,
    bulk = bulk
  ))
}

lighting_folders <- character(0)
last_run <- NULL

for (i in seq_len(nrow(filter_grid))) {
  cfg <- filter_grid[i, , drop = FALSE]

  fp_p_thr <- cfg$fp_rna_p
  fp_r_thr <- cfg$fp_rna_r
  rna_p_thr <- cfg$fp_rna_p
  rna_r_thr <- cfg$fp_rna_r

  use_atac <- isTRUE(cfg$use_atac)
  atac_p_thr <- if (use_atac) cfg$atac_p else NULL
  atac_r_thr <- if (use_atac) cfg$atac_r else NULL
  atac_same_dir <- if (use_atac) isTRUE(cfg$atac_same_dir) else FALSE

  tf_p_thr <- cfg$tf_p
  tf_r_thr <- cfg$tf_r

  fp_gene_corr_kept_rna_filtered <- filter_fp_rna_atac_corr(
    tbl = fp_gene_corr_joined,
    fp_p_thr = fp_p_thr,
    fp_r_abs_min = fp_r_thr,
    rna_p_thr = rna_p_thr,
    rna_r_abs_min = rna_r_thr,
    fp_p_col = fp_p_col_use,
    rna_p_col = rna_p_col_use,
    require_same_dir_fp_rna = TRUE,
    use_atac = use_atac,
    atac_p_thr = atac_p_thr,
    atac_r_abs_min = atac_r_thr,
    atac_p_col = atac_p_col_use,
    require_same_dir_atac_rna = atac_same_dir
  )

  # Encode the thresholds used into the output folder name.
  # Example suffix:
  #   FP[p_adj<=0.05|r>=0.5]_RNA[p_adj<=0.05|r>=0.5]_ATAC[none]
  # or
  #   ..._ATAC[p_adj<=0.01|r>=0.7|sameDirFALSE]
  p_label <- function(p_col) {
    if (grepl("adj", p_col, ignore.case = TRUE)) "pAdj" else "p"
  }

  suffix_fp <- sprintf("FP_%s%.2g_r%.1f", p_label(fp_p_col_use), fp_p_thr, fp_r_thr)
  suffix_rna <- sprintf("RNA_%s%.2g_r%.1f", p_label(rna_p_col_use), rna_p_thr, rna_r_thr)
  suffix_tf <- {
    if (!is.finite(tf_p_thr) || !is.finite(tf_r_thr)) stop("TF thresholds must be finite")
    sprintf("TF_%s%.2g_r%.1f", p_label(tf_p_col_use), tf_p_thr, tf_r_thr)
  }
  suffix_atac <- if (!use_atac) {
    "ATAC_none"
  } else {
    # guard against accidental NAs
    if (!is.finite(atac_p_thr) || !is.finite(atac_r_thr)) {
      stop("ATAC thresholds must be finite when use_atac=TRUE")
    }
    sprintf("ATAC_%s%.2g_r%.1f_sameDir%s", p_label(atac_p_col_use), atac_p_thr, atac_r_thr, if (atac_same_dir) "TRUE" else "FALSE")
  }

  lighting_folder <- file.path(
    base_dir,
    paste(
      "lighting",
      "fp_tf_corr_FDR", threshold_fp_tf_corr_p,
      regulation_priors, db,
      "regulated_genes", regulated_genes,
      "delta_link", delta_link,
      suffix_fp, suffix_rna, suffix_tf, suffix_atac,
      sep = "_"
    )
  )
  print(lighting_folder)
  lighting_folders <- c(lighting_folders, lighting_folder)

  last_run <- run_lighting_pipeline_one(
    fp_gene_corr_kept = fp_gene_corr_kept_rna_filtered,
    lighting_folder = lighting_folder,
    grn_set = grn_set,
    regulated_genes = regulated_genes,
    delta_link = delta_link,
    tf_p_thr = tf_p_thr,
    tf_r_min = tf_r_thr,
    tf_p_col = tf_p_col_use
  )
}

# Keep a final lighting_folder value for any downstream interactive code blocks.
if (length(lighting_folders) > 0) {
  lighting_folder <- lighting_folders[[length(lighting_folders)]]
}

# Expose the last run's outputs for downstream interactive analysis blocks.
if (!is.null(last_run)) {
  bulk <- last_run$bulk
  filtered_csvs <- bulk$filtered_paths
  lda_out <- last_run$lda_out
  sum_out <- last_run$sum_out
}

# Topic models ---------------------------------------------------------------
# Moved to dev/08_benchmark_topic_methods.R (source that script to run benchmarks)


# TODO TF HITS score summary plot (waterfall)

# Step 5. Generate interactive Topic & TF regulatory hub subnetwor --------

# ---- Step 5A. Single Δ-topic from a known pair ----

comp_csv <- file.path(lighting_folder, "AsPC1_10_Gln.Arg_vs_AsPC1_10_FBS_delta_links_filtered_lda_K20.csv")
summary_csv <- file.path(lighting_folder, "AsPC1_10_Gln.Arg_vs_AsPC1_10_FBS_delta_links_filtered_lda_K20_summary.csv")

# Render Topic 7 with visual rules (Δ = stress − control)
episcope::render_link_network_delta_topic_simple(
  comp_csv = comp_csv,
  summary_csv = summary_csv,
  topic_id = 7,
  edge_filter_min = 1,
  gene_fc_thresh = 1.5,
  de_reference = "str_over_ctrl",
  top_n_tfs_per_topic = 20,
  out_html = file.path(
    dirname(comp_csv),
    sprintf(
      "%s_topic-%d_subnetwork_delta.html",
      tools::file_path_sans_ext(basename(comp_csv)), 7
    )
  ),
  verbose = TRUE
)


# ---- Step 5B. Bulk Δ-topic across many summaries (auto-pairs) ----
summary_csvs <- list.files(
  lighting_folder,
  pattern = "_filtered_lda_K20_summary\\.csv$",
  full.names = TRUE
)

for (sum_path in summary_csvs) {
  # Try to find the matching annotated comparison file
  cand1 <- sub("_summary\\.csv$", ".csv", sum_path)
  cand2 <- sub("_filtered_lda_K\\d+_summary\\.csv$", "_delta_links.csv", sum_path)
  comp_csv <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else NA_character_
  if (!isTRUE(file.exists(comp_csv))) {
    cli::cli_inform("Skip (no comp CSV found for): {sum_path}")
    next
  }

  S <- readr::read_csv(sum_path, show_col_types = FALSE)
  topic_col <- if ("topic" %in% names(S)) "topic" else if ("main_topic" %in% names(S)) "main_topic" else NA_character_
  if (!isTRUE(topic_col %in% names(S))) next
  if ("topic_rank" %in% names(S)) S <- S[order(S$topic_rank), , drop = FALSE]
  topics <- head(unique(as.integer(S[[topic_col]])), 20)

  for (t in topics) {
    episcope::render_link_network_delta_topic_simple(
      comp_csv = comp_csv,
      summary_csv = sum_path,
      topic_id = t,
      edge_filter_min = 1,
      gene_fc_thresh = 1.5,
      de_reference = "str_over_ctrl",
      top_n_tfs_per_topic = 20,
      out_html = NULL,
      verbose = TRUE
    )
  }
}

# 5C TF centric plots
render_tf_hub_delta_network(
  comp_csv = file.path(
    lighting_folder,
    "BG_24h_culture_5d_vs_culture_6d_delta_links_filtered_lda_K20.csv"
  ),
  input_tf = "TBX21",
  edge_filter_min = 1,
  edge_filter_on = "either",
  gene_fc_thresh = 1.5,
  de_reference = "str_over_ctrl",
  # ring_tf_direct_only = TRUE,   # only direct TFs form the ring
  motif_db = "JASPAR2024",
  verbose = TRUE
)

comp_csvs <- list.files(lighting_folder, pattern = "_delta_links_filtered_lda_K20\\.csv$", full.names = TRUE)
tf_list <- c(
  "TBX21", "TCF3", "ZNF610", "SOX4", "ZNF85", "SP4", "REL", "MTF1", "STAT4", "JUNB",
  "MITF", "IRF1", "IRF8", "NRF1", "NFKB1", "NFKB2", "USF2", "STAT2", "STAT5A", "EGR2", "KLF6"
)

for (f in comp_csvs) {
  for (tf in tf_list) {
    res <- try(
      {
        render_tf_hub_delta_network(
          comp_csv = f,
          input_tf = tf,
          edge_filter_min = 1,
          edge_filter_on = "either",
          gene_fc_thresh = 1.5,
          de_reference = "str_over_ctrl",
          motif_db = "JASPAR2024",
          verbose = TRUE
        )
        TRUE
      },
      silent = TRUE
    )
  }
}
