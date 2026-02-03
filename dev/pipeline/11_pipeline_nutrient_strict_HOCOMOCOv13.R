library(episcope)

source("R/utils_config.R")
source("R/utils_logging.R")

source("R/utils_step1_footprints.R")
source("R/utils_step1_align_footprints.R")
source("R/utils_step1_grn_preprocess.R")
source("R/utils_step1_motif_clustering.R")

source("R/utils_step2_connect_tf_genes.R")
source("R/utils_step2_grn_pipeline.R")

source("R/utils_step3_grn_filter.R")
source("R/utils_step3_grn_diff.R")
source("R/utils_step3_topic_warplda.R")

load_config("dev/config/pdac_nutrient_stress_HOCOMOCOv13.yaml")



# ──────────────────────────────────────────────────────────────────────────────
# Turn modules ON/OFF
# ──────────────────────────────────────────────────────────────────────────────
do_load_footprints_preprocess    <- TRUE
do_tf_binding_sites_prediction   <- TRUE
do_tf_to_target_genes_prediction <- TRUE


# Load footprint data and preprocess --------------------------------------
if (do_load_footprints_preprocess == TRUE) {
  fp_cache_dir <- file.path(base_dir, "cache")
  fp_manifest <- load_footprints(
    root_dir = fp_root_dir,
    db_name = db,
    out_dir = file.path(fp_cache_dir, paste0("fp_", db))
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
  run_fp_motif_clustering_pre_if_needed(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, motif_db = motif_db)
  fp_motif_clust_data   <- run_fp_motif_clustering(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, motif_db = motif_db, mode = "data", target_clusters = 220, qc_mode = "none", save_motif_db = TRUE)
  fp_motif_clust_hybrid <- run_fp_motif_clustering(fp_aligned = fp_aligned, base_dir = base_dir, ref_db = db, motif_db = motif_db, mode = "hybrid", target_clusters = 220, qc_mode = "none", save_motif_db = TRUE)

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



  # grn_set <- build_grn_set(
  #   fp_score      = fp_aligned$fp_score,
  #   fp_bound      = fp_aligned$fp_bound,
  #   fp_annotation = fp_aligned$fp_annotation,
  #   atac_score    = atac_score,
  #   atac_overlap  = atac_overlap,
  #   rna           = lenient_rna,
  #   metadata      = lenient_metadata,
  #   tf_list       = tf_list,
  #   motif_db      = motif_db,
  #   label_col     = "strict_match_rna",
  #   expected_n    = 23
  # )

  grn_set <- grn_add_rna_expressed(grn_set, label_col = "strict_match_rna", threshold_gene_expr = threshold_gene_expr)
  grn_set <- grn_add_fp_score_condition(grn_set, label_col = "strict_match_rna")
  grn_set <- grn_add_fp_bound_condition(grn_set, label_col = "strict_match_rna", threshold_fp_score = threshold_fp_score)
  grn_set <- grn_filter_fp_bound_condition(grn_set, min_bound = 1L, use_parallel = TRUE)
  # assert_fp_alignment(grn_set)

  step1_out_dir <- file.path(base_dir, "predict_tf_binding_sites")
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID")
  write_grn_outputs(grn_set, out_dir = step1_out_dir, db = db, qn_base_dir = base_dir)

}


# Predict TF binding sites -------------------------------------------------
if (do_tf_binding_sites_prediction == TRUE) {
  if (!exists("grn_set") || !is.list(grn_set)) {
    .log_abort("`grn_set` not found. Run footprint preprocessing before TF binding site prediction.")
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


# Connect TF-occupied enhancers to target genes ---------------------------
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


# Diff networks and topic analysis ----------------------------------------
step2_out_dir <- file.path(base_dir, "connect_tf_target_genes")

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

  motif_path <- resolve_motif_db_path("HOCOMOCOv13", ref_genome = ref_genome)
  # motif_db$gene_symbol <- motif_db$HGNC
  if (is.data.frame(motif_db) && all(c("sub_cluster_name", "gene_symbol") %in% names(motif_db))) {
    motif_info <- build_tf_cluster_map_from_motif(motif_db)
  } else {
    motif_info <- build_tf_cluster_map_from_motif(motif_path)
  }

  k_grid_default <- c(2:15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100)
  k_single_map <- list(Panc1 = c(10L, 14L, 20, 35, 50, 70)) # HPAFII = 12L,  , AsPC1 = 10L
  celllines <- c("HPAFII") # "AsPC1", "HPAFII",


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
