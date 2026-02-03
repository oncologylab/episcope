library(episcope)

source("R/utils_config.R")
source("R/utils_logging.R")

source("R/utils_step1_footprints.R")
source("R/utils_step1_align_footprints.R")
source("R/utils_step1_grn_preprocess.R")
source("R/utils_step1_pipeline_helpers.R")
source("R/utils_step1_motif_clustering.R")
source("R/utils_step1_qc.R")

source("R/utils_step2_connect_tf_genes.R")
source("R/utils_step2_grn_pipeline.R")
source("R/utils_step2_qc.R")

source("R/utils_step3_grn_filter.R")
source("R/utils_step3_grn_diff.R")
source("R/utils_step3_topic_warplda.R")

load_config("dev/config/pdac_nutrient_stress_strict_jaspar2024.yaml")

expected_n <- if (exists("expected_n")) expected_n else NULL

# ──────────────────────────────────────────────────────────────────────────────
# Turn modules ON/OFF
# ──────────────────────────────────────────────────────────────────────────────
do_load_footprints_preprocess    <- TRUE
do_tf_binding_sites_prediction   <- TRUE
do_tf_to_target_genes_prediction <- TRUE
verbose <- TRUE


# Predict TF binding sites -------------------------------------------------
# Load footprint data and preprocess --------------------------------------
if (do_load_footprints_preprocess == TRUE) {
  # Inputs (from YAML): fp_root_dir, base_dir, db, ref_genome, thresholds, etc.
  fp_cache_dir <- file.path(base_dir, "cache")
  fp_manifest <- load_footprints(
    root_dir = fp_root_dir,
    db_name = db,
    out_dir = file.path(fp_cache_dir, paste0("fp_", db))
  )

  # readr::write_csv(fp_manifest, file.path(fp_cache_dir, sprintf("fp_%s_manifest.csv", db)))

  # options(future.globals.maxSize = 32 * 1024^3)
  # Collapse overlapping footprints and align peaks across samples
  fp_aligned <- align_footprints(fp_manifest, mid_slop = 10L, round_digits = 1L, score_match_pct = 0.8, cache_dir = fp_cache_dir, cache_tag = db, output_mode = "distinct")
  plot_fp_merge_summary(fp_aligned, out_dir = file.path(base_dir, "predict_tf_binding_sites"), db = db, verbose = TRUE)

  # Assign TF motifs to consensus motif clusters
  run_fp_motif_clustering_pre_if_needed(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, motif_db = motif_db)
  fp_motif_clust_data   <- run_fp_motif_clustering(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, motif_db = motif_db, mode = "data", target_clusters = 220, qc_mode = "fast", save_motif_db = TRUE)
  fp_motif_clust_hybrid <- run_fp_motif_clustering(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, motif_db = motif_db, mode = "hybrid", target_clusters = 165, qc_mode = "fast", save_motif_db = TRUE)

  # Build combined multi-omic data object
  omics_data <- load_multiomic_data(
    fp_aligned = fp_aligned,
    label_col = "strict_match_rna",
    expected_n = expected_n,
    tf_list = tf_list,
    motif_db = motif_db,
    threshold_gene_expr = threshold_gene_expr,
    threshold_fp_score = threshold_fp_score,
    use_parallel = TRUE,
    verbose = verbose
  )

  step1_out_dir <- file.path(base_dir, "predict_tf_binding_sites")
  write_grn_outputs(omics_data, out_dir = step1_out_dir, db = db, qn_base_dir = base_dir)
  plot_fp_norm_bound_qc(
    omics_data = omics_data,
    out_dir = step1_out_dir,
    db = db,
    threshold_fp_score = threshold_fp_score,
    max_points = 100000L,
    verbose = TRUE
  )
  plot_gene_expr_qc(
    omics_data = omics_data,
    out_dir = step1_out_dir,
    db = db,
    threshold_gene_expr = threshold_gene_expr,
    verbose = TRUE
  )
}


# Predict TF binding sites -------------------------------------------------
if (do_tf_binding_sites_prediction == TRUE) {
  if (!exists("omics_data") || !is.list(omics_data)) {
    .log_abort("`omics_data` not found. Run footprint preprocessing before TF binding site prediction.")
  }

  # Correlate TF expression vs footprint scores (canonical + all modes)
  step1_out_dir <- file.path(base_dir, "predict_tf_binding_sites")

  omics_data <- correlate_tf_to_fp(
    omics_data = omics_data,
    mode = "canonical",
    out_dir = step1_out_dir,
    label_col = "strict_match_rna",
    r_thr = threshold_fp_tf_corr_r,
    p_thr = threshold_fp_tf_corr_p,
    db = db,
    cores_pearson = 20L,
    cores_spearman = 36L,
    chunk_size = 5000L,
    min_non_na = 5L,
    qc = TRUE,
    write_bed = FALSE
  )
}


# Connect TFs to target genes ---------------------------------------------
# Connect TF-occupied enhancers to target genes ---------------------------
if (do_tf_to_target_genes_prediction == TRUE) {
  step2_out_dir <- file.path(base_dir, "connect_tf_target_genes")
  dir.create(step2_out_dir, recursive = TRUE, showWarnings = FALSE)

  # Link TFBS to candidate target genes (GeneHancer / window / loops)
  link_mode <- if (exists("tf_target_link_mode")) tf_target_link_mode else "genehancer"
  link_mode <- match.arg(link_mode, c("genehancer", "window", "loops"))

  gh_std <- NULL
  atac_gene_pairs <- NULL

  if (link_mode == "genehancer") {
    gh_path <- if (exists("genehancer_path")) genehancer_path else file.path("inst", "extdata", "GeneHancer_v5.24_elite_panc.csv")
    gh_std <- load_genehancer_panc(gh_path)
  } else if (link_mode == "window") {
    flank_bp <- if (exists("link_window_bp")) link_window_bp else 30000L
    gene_annot_ref <- if (exists("gene_annot_ref")) gene_annot_ref else episcope::episcope_build_gene_annot(ref_genome)
    id_col <- if (exists("link_gene_id_col")) link_gene_id_col else "HGNC"
    gh_std <- episcope_make_windowed_gh(
      peaks      = omics_data$atac_score,
      genes      = unique(c(omics_data$rna$HGNC, omics_data$rna$ensembl_gene_id)),
      flank_bp   = flank_bp,
      mode       = "TSS",
      gene_annot = gene_annot_ref,
      id_col     = id_col
    )
  } else if (link_mode == "loops") {
    loop_tbl <- if (exists("loop_tbl")) loop_tbl else NULL
    if (is.null(loop_tbl) && exists("loop_path") && file.exists(loop_path)) {
      loop_tbl <- if (grepl("\\.tsv$", loop_path)) {
        readr::read_tsv(loop_path, show_col_types = FALSE)
      } else {
        readr::read_csv(loop_path, show_col_types = FALSE)
      }
    }
    if (is.null(loop_tbl)) .log_abort("Loop mode requires `loop_tbl` or `loop_path`.")

    if (!"atac_peak" %in% names(loop_tbl)) {
      if (all(c("chrom", "start", "end") %in% names(loop_tbl))) {
        loop_tbl$atac_peak <- paste0(loop_tbl$chrom, ":", loop_tbl$start, "-", loop_tbl$end)
      } else if ("peak_ID" %in% names(loop_tbl)) {
        loop_tbl$atac_peak <- loop_tbl$peak_ID
      }
    }
    gene_col <- intersect(c("gene_key", "connected_gene", "gene", "target_gene", "HGNC", "ensembl_gene_id"), names(loop_tbl))[1]
    if (is.na(gene_col) || !"atac_peak" %in% names(loop_tbl)) {
      .log_abort("Loop table must contain atac_peak and a gene column (gene_key/connected_gene/gene/target_gene/HGNC/ensembl_gene_id).")
    }
    atac_gene_pairs <- loop_tbl |>
      dplyr::transmute(atac_peak = .data$atac_peak, gene_key = .data[[gene_col]]) |>
      dplyr::distinct()
  }

  # Correlate TF->gene and FP->gene (Spearman then Pearson)
  options(future.globals.maxSize = 64 * 1024^3)
  fp_res_full_pearson <- correlate_fp_to_genes(
    grn_set          = omics_data,
    atac_gene_corr_kept = atac_gene_pairs,
    gh_tbl           = gh_std,
    gene_mode        = "both",
    fp_score_tbl     = omics_data$fp_score_condition_qn,
    rna_tbl          = omics_data$rna_condition,
    fdr              = threshold_fp_gene_corr_p,
    r_abs_min        = threshold_fp_gene_corr_abs_r,
    method           = "pearson",
    workers          = 20,
    cache_dir        = file.path(base_dir, "cache", "fp_gene_corr"),
    cache_tag        = paste0("nutrient_", link_mode, "_pearson"),
    cache_chunk_size = 5000L,
    cache_verbose    = TRUE
  )

  fp_res_full_spearman <- correlate_fp_to_genes(
    grn_set          = omics_data,
    atac_gene_corr_kept = atac_gene_pairs,
    gh_tbl           = gh_std,
    gene_mode        = "both",
    fp_score_tbl     = omics_data$fp_score_condition_qn,
    rna_tbl          = omics_data$rna_condition,
    fdr              = threshold_fp_gene_corr_p,
    r_abs_min        = threshold_fp_gene_corr_abs_r,
    method           = "spearman",
    workers          = 30,
    cache_dir        = file.path(base_dir, "cache", "fp_gene_corr"),
    cache_tag        = paste0("nutrient_", link_mode, "_spearman"),
    cache_chunk_size = 5000L,
    cache_verbose    = TRUE
  )

  fp_gene_corr_pearson_filt <- filter_fp_gene_corr_by_tf_annotation(
    fp_gene_corr_full = fp_res_full_pearson$fp_gene_corr_full,
    fp_annotation = omics_data$fp_annotation,
    r_thr = threshold_fp_gene_corr_abs_r,
    p_adj_thr = threshold_fp_gene_corr_p
  )

  fp_gene_corr_spearman_filt <- filter_fp_gene_corr_by_tf_annotation(
    fp_gene_corr_full = fp_res_full_spearman$fp_gene_corr_full,
    fp_annotation = omics_data$fp_annotation,
    r_thr = threshold_fp_gene_corr_abs_r,
    p_adj_thr = threshold_fp_gene_corr_p
  )

  fp_links_combined <- combine_fp_gene_corr_methods(
    fp_pearson  = fp_gene_corr_pearson_filt,
    fp_spearman = fp_gene_corr_spearman_filt,
    rna_tbl     = omics_data$rna_condition,
    rna_method  = "pearson",
    rna_cores   = 20
  )

  fp_links_filtered <- filter_links_by_fp_rna_criteria(
    fp_links_combined,
    fp_p_adj_thr = threshold_fp_gene_corr_p,
    fp_r_thr = threshold_fp_gene_corr_abs_r,
    require_pos_rna = TRUE
  )

  fp_score_thr <- if (exists("fp_score_threshold")) fp_score_threshold else threshold_fp_score
  atac_score_thr <- if (exists("atac_score_threshold")) atac_score_threshold else 0
  require_atac_score <- isTRUE(atac_score_thr > 0)
  atac_score_tbl_use <- NULL
  if (isTRUE(require_atac_score)) {
    if (is.null(omics_data$atac_score_condition)) {
      omics_data <- grn_add_atac_score_condition(omics_data, label_col = "strict_match_rna")
    }
    atac_score_tbl_use <- omics_data$atac_score_condition
  }

  status_res <- build_link_status_matrix(
    links = fp_links_filtered$links,
    fp_bound = omics_data$fp_bound_condition,
    rna_expressed = omics_data$rna_expressed,
    tf_col = "TF",
    gene_col = "gene_key",
    peak_col = "peak_ID",
    out_file = file.path(step2_out_dir, sprintf("tf_gene_link_status_matrix_%s.csv", db)),
    chunk_size = 50000L,
    return_keep = TRUE,
    filter_any = TRUE,
    verbose = TRUE,
    fp_score_tbl = omics_data$fp_score_condition_qn,
    fp_score_threshold = fp_score_thr,
    atac_score_tbl = atac_score_tbl_use,
    atac_score_threshold = atac_score_thr,
    require_fp_bound = TRUE,
    require_gene_expr = TRUE,
    require_fp_score = TRUE,
    require_atac_score = require_atac_score,
    return_tbl = TRUE
  )

  fp_links_kept <- fp_links_filtered$links[status_res$keep, , drop = FALSE]
  readr::write_csv(fp_links_kept, file.path(step2_out_dir, sprintf("tf_gene_links_filtered_%s.csv", db)))

  build_tf_target_overview(
    links = fp_links_kept,
    link_status = status_res$status_tbl,
    out_file = file.path(step2_out_dir, sprintf("tf_target_link_overview_%s.csv", db)),
    verbose = TRUE
  )

  link_summary <- summarize_link_activity(
    link_status = status_res$status_tbl,
    out_dir = step2_out_dir,
    db = db,
    prefix = "step2",
    verbose = TRUE
  )
  plot_link_activity_qc(
    summary_total = link_summary$summary_total,
    out_dir = step2_out_dir,
    db = db,
    prefix = "step2",
    verbose = TRUE
  )

  fp_gene_corr_use <- fp_res_full_pearson$fp_gene_corr_full |>
    dplyr::semi_join(
      fp_links_kept |> dplyr::select(peak_ID, gene_key),
      by = c("fp_peak" = "peak_ID", "gene_key" = "gene_key")
    )

  if (is.null(omics_data$fp_variance) || is.null(omics_data$rna_variance)) {
    hv_variance <- precompute_hvf_hvg_variance(omics_data, cores = 20, chunk_size = 50000L)
    omics_data$fp_variance <- hv_variance$fp_variance
    omics_data$rna_variance <- hv_variance$rna_variance
  }

  basal_links_step2 <- make_basal_links(
    fp_gene_corr_kept = fp_gene_corr_use,
    fp_annotation     = omics_data$fp_annotation,
    out_dir           = file.path(step2_out_dir, "basal_links_tmp"),
    prefix            = "step2",
    rna_tbl           = omics_data$rna_condition,
    rna_method        = "pearson",
    rna_cores         = 20,
    fp_variance       = omics_data$fp_variance,
    rna_variance      = omics_data$rna_variance
  )

  omics_data_cond <- prepare_grn_set_for_light_by_condition(
    omics_data,
    label_col = "strict_match_rna"
  )

  light_by_condition(
    ds = omics_data_cond,
    basal_links = basal_links_step2,
    out_dir = step2_out_dir,
    prefix = "step2",
    label_col = "strict_match_rna",
    link_score_threshold = link_score_threshold,
    fp_score_threshold = fp_score_thr,
    tf_expr_threshold = threshold_tf_expr,
    fp_bound_tbl = omics_data$fp_bound_condition,
    rna_expressed_tbl = omics_data$rna_expressed,
    atac_score_tbl = atac_score_tbl_use,
    atac_score_threshold = atac_score_thr,
    require_atac_score = require_atac_score,
    fp_annotation_tbl = omics_data$fp_annotation,
    require_fp_bound = TRUE,
    require_gene_expr = TRUE,
    gene_expr_threshold = 1L,
    filter_active = FALSE,
    use_parallel = TRUE,
    workers = 8,
    fp_variance_tbl = omics_data$fp_variance,
    rna_variance_tbl = omics_data$rna_variance
  )

  extract_link_info_by_condition(
    out_dir = step2_out_dir,
    prefix = "step2",
    read_tables = FALSE,
    verbose = TRUE
  )
}

# Diff networks and topic analysis ----------------------------------------
step2_out_dir <- file.path(base_dir, "example_tf_target_genes")

# Per-condition tables (per-cell vs 10_FBS)
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
  if (!length(delta_links_csv)) .log_abort("No *_delta_links.csv files found in {step2_out_dir}")

  edges_all <- load_delta_links_many(delta_links_csv, keep_original = TRUE)
  topic_root <- file.path(base_dir, "topic_models_lenient")

  motif_path <- resolve_motif_db_path("JASPAR2024", ref_genome = ref_genome)
  motif_db$gene_symbol <- motif_db$HGNC
  if (is.data.frame(motif_db) && all(c("sub_cluster_name", "gene_symbol") %in% names(motif_db))) {
    motif_info <- build_tf_cluster_map_from_motif(motif_db)
  } else {
    motif_info <- build_tf_cluster_map_from_motif(motif_path)
  }

  k_grid_default <- c(2:15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100)
  k_single_map <- list(Panc1 = c(10L, 14L, 20, 35, 50, 70)) # HPAFII = 12L,  , AsPC1 = 10L
  celllines <- c("Panc1") # "AsPC1", "HPAFII",


  library(enrichR)
  library(episcope)
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
