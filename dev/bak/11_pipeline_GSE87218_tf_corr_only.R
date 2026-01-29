library(episcope)
load_config("episcope_grn.yaml")
source("R/utils_connect_tf_enhancers_to_target_genes.R")
source("R/utils_load_footprints_process.R")

# ──────────────────────────────────────────────────────────────────────────────
# turn heavy analyses on/off
# ──────────────────────────────────────────────────────────────────────────────
do_load_footprints_preprocess    <- TRUE
# Including
# 1. stream_write_overviews_by_motif
# 2. trim_fp_manifest
# 3. trim_manifest_annots
# 4. preprocess_align_fp_peaks
# 5. quantile_normalize_fp_unique
# 6. save_fp_aligned_normalized
# 7. process_motifs_in_parallel (filtering)


# stream_write_overviews_by_motif -> load_footprints()
# trim_fp_manifest                -> fp_manifest_trim()
# trim_manifest_annots            -> fp_manifest_trim_annots()
# preprocess_align_fp_peaks       -> align_footprints()
# quantile_normalize_fp_unique    -> qn_footprints()
# save_fp_aligned_normalized      -> save_footprints()

do_tf_binding_sites_prediction   <- TRUE
do_tf_to_target_genes_prediction <- TRUE
do_build_grn                     <- FALSE
do_differential_grn              <- FALSE
do_topic_subplot_single          <- FALSE
do_topic_subplot_bulk            <- FALSE
do_tf_centric_subplot_single     <- FALSE
do_tf_centric_subplot_bulk       <- FALSE






# fp_manifest <- load_footprints(
#   root_dir   = "Y:/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
#   db_name    = "JASPAR2024",
#   out_dir    = "Z:/episcope_test/nutrient_stress/fp_jaspar2024",
#   n_workers  = 20,
#   # n_samples  = 2,
#
#   verbose    = TRUE
# )
# readr::write_csv(fp_manifest, "Z:/episcope_test/nutrient_stress/fp_jaspar2024_manifest.csv")


# fp_manifest <- load_footprints(
#   root_dir   = "/data/homes/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
#   db_name    = "HOCOMOCOv13",
#   out_dir    = "Z:/episcope_test/nutrient_stress_test/hocomocov13",
#   n_workers  = 20,
#   n_samples  = 2,
#   verbose    = TRUE
# )
# readr::write_csv(fp_manifest, "Z:/episcope_test/nutrient_stress_test/fp_hocomocov13_manifest.csv")


# fp_manifest <- load_footprints(
#   root_dir   = "Z:/fptools_test/GSE96771_ATAC/TOBIAS_merged_peaks/",
#   db_name    = "JASPAR2024",
#   out_dir    = "Z:/episcope_test/GSE96771_ATAC",
#   n_workers  = 20,
#   n_samples  = 2,
#   verbose    = TRUE
# )
# readr::write_csv(fp_manifest, "Z:/episcope_test/GSE96771_ATAC/fp_jaspar2024_manifest.csv")

base_dir <- "/data/homes/yl814/episcope_test/GSE87218_ATAC"

# Step 1. Predict TF binding sites ----------------------------------------

if (do_load_footprints_preprocess == TRUE) {
  # fp_manifest <- load_footprints(
  #   root_dir   = "/data/homes/yl814/tbio/episcope_test/GSE87218_ATAC/TOBIAS_merged_peaks",
  #   db_name    = "JASPAR2024",
  #   out_dir    = file.path(base_dir, "fp_jaspar2024"),
  #   n_workers  = 20,
  #   # n_samples  = 2,
  #   verbose    = TRUE
  # )
  # readr::write_csv(fp_manifest, file.path(base_dir, "fp_jaspar2024_manifest.csv"))

  fp_manifest <- readr::read_csv(file.path(base_dir, "fp_jaspar2024_manifest.csv"))

  # (Optional, when using HOCOMOCO database)
  # fp_manifest <- fp_manifest_trim(fp_manifest) # renames files on disk by default
  # Overwrite every annotation CSV referenced by the manifest:
  # summary_tbl <- fp_manifest_trim_annots(fp_manifest, n_workers = 18, verbose = TRUE)

  # Inspect what changed:
  # dplyr::count(summary_tbl, status)
  # sum(summary_tbl$n_fixed, na.rm = TRUE)

  options(future.globals.maxSize = 32 * 1024^3)
  # Align the peaks based on the peak similarity
  # fp_aligned <- align_footprints(fp_manifest,
  #                                mid_slop        = 10L, # midpoint tolerance (bp)
  #                                round_digits    = 1L, # round scores before comparing vectors
  #                                score_match_pct = 0.8
  # ) # fraction of samples that must match (<=1, >=0)


  # length(unique(fp_aligned$id_map$peak_ID))
  # length(unique(fp_aligned$id_map$fp_peak_bak))

  # readr::write_csv(fp_aligned$fp_bound, file.path(base_dir, sprintf("fp_bounds_%s.csv", db)))
  # readr::write_csv(fp_aligned$fp_score, file.path(base_dir, sprintf("fp_scores_%s.csv", db)))
  # readr::write_csv(fp_aligned$fp_annotation, file.path(base_dir, sprintf("fp_annotation_%s.csv", db)))


  # fp_score_raw <- fp_aligned$fp_score
  # fp_aligned_normalized <- fp_aligned
  # fp_aligned_normalized$fp_score <- qn_footprints(fp_aligned_normalized$fp_score, id_col = "peak_ID")
  # readr::write_csv(fp_aligned_normalized$fp_score, file.path(base_dir, sprintf("fp_scores_qn_%s.csv", db)))



  # fp_aligned_normalized_manifest <- save_footprints(fp_aligned_normalized, out_dir  = paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_jaspar2024/fp_aligned_normalized_", db))
  # readr::write_csv(fp_aligned_normalized_manifest, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_jaspar2024/fp_aligned_normalized_manifest_", db, ".csv"))


  sample_metadata <- readr::read_tsv(file.path(base_dir, "GSE87218_samples_filtered.txt"), na = "NA")
  # stric_metadata   <- sample_metadata |> dplyr::filter(!is.na(strict_match_rna))
  sample_metadata$id <- sample_metadata$ID
  if (db == "jaspar2024") {
    motif_db <- readr::read_tsv("inst/extdata/genome/JASPAR2024.txt")
  } else if (db == "hocomocov13") {
    motif_db <- readr::read_tsv("inst/extdata/genome/HOCOMOCOv13.txt")
  }

  tf_list <- motif_db |>
    tidyr::separate_rows(HGNC, sep = "::") |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::distinct(HGNC) |>
    dplyr::pull(HGNC)

  # Load RNA data and ATAC data
  atac_data <- readr::read_tsv(file.path(base_dir, "All_ATAC.master_table.narrowpeaks.10mil.txt"))

  # atac_data <- readr::read_tsv("/data/homes/yl814/episcope_test/nutrient_stress/All_ATAC.master_table.narrowpeaks.10mil.txt")

  # motif_db
  # base_dir

  atac_out <- load_atac(atac_data, sort_peaks = TRUE)
  atac_score <- atac_out$score
  atac_overlap <- atac_out$overlap

  # rna <- readr::read_csv("inst/extdata/HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")
  # rna <- readr::read_csv("/data/homes/yl814/episcope_test/GSE87218_ATAC/GSE85243_RNA_group_averages_mapped_to_SRR.csv")
  # rna <- readr::read_csv("/data/homes/yl814/episcope_test/nutrient_stress/HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")

  rna <- readr::read_tsv(file.path(base_dir, "GSE85243_counts.median_ratio_normalized_group_averages_mapped_to_SRR.ALL.tsv"))
  rna <- clean_hgnc(rna) # clean "HGNC" column
  # Filter rna expression genes/tfs needs to reach the threshold in at least 1 sample (group size)
  rna <- filter_rna_expr(rna, tf_list, hgnc_col = "HGNC", gene_min = threshold_gene_expr, tf_min = threshold_tf_expr, min_samples = 1L)
  # make two RNA tables with columns already named by ATAC ids ──

  # strict_rna
  # smap <- dplyr::transmute(stric_metadata, old = strict_match_rna, new = id)
  # strict_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
  # nm <- names(strict_rna); nm[match(smap$old, nm)] <- smap$new
  # strict_rna <- strict_rna |> `names<-`(nm) |> dplyr::as_tibble()

  # Build RNA and ATAC data list
  rna_atac_build_args <- list(
    atac_score    = atac_score,
    atac_overlap  = atac_overlap,
    rna           = rna,
    metadata      = sample_metadata,
    tf_list       = tf_list,
    motif_db      = motif_db, # tibble
    label_col     = "Sample",
    expected_n    = 13
  )

  # fp_aligned_normalized_manifest <- save_footprints(fp_aligned_normalized, out_dir = file.path(base_dir, sprintf("fp_aligned_normalized_%s", db)))
  # readr::write_csv(fp_aligned_normalized_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_manifest_%s.csv", db)))

  # RNA expression and ATAC peak based, footprint filtering
  # fp_aligned_normalized_filtered_manifest <- filter_footprints(
  #   fp_manifest         = fp_aligned_normalized_manifest,
  #   out_dir             = file.path(base_dir, sprintf("fp_aligned_normalized_filtered_%s", db)),
  #   build_args          = rna_atac_build_args,
  #   n_workers           = 10,
  #   skip_existing       = FALSE,
  #   threshold_tf_expr   = threshold_tf_expr,
  #   verbose             = TRUE
  # )

  # readr::write_csv(fp_aligned_normalized_filtered_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_filtered_manifest_%s.csv", db)))
  fp_aligned_normalized_filtered_manifest <-  readr::read_csv(file.path(base_dir, sprintf("fp_aligned_normalized_filtered_manifest_%s.csv", db)))

  tfs <- tf_list[tf_list %in% unique(rna$HGNC)]

  for (tf in tfs) {
    tf_vec <- rna |> dplyr::filter(HGNC == tf) |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
    res <- tf_corr_footprints_all_tfbs(
      fp_manifest = fp_aligned_normalized_filtered_manifest,
      tf_name  = tf,
      tf_expr  = tf_vec,
      out_dir  = file.path(base_dir, "predicted_all_tfbs")
    )
  }
}


if (do_tf_binding_sites_prediction == TRUE) {

  # TF expr to footprints
  # fp_corr_manifest <- tf_corr_footprints(
  #   fp_manifest = fp_aligned_normalized_filtered_manifest,
  #   out_dir     = file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_%s", db)),
  #   build_args  = rna_atac_build_args,
  #   n_workers   = 15,
  #   cor_method  = "pearson",
  #   min_non_na  = 5L
  # )


  # readr::write_csv(fp_corr_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))
  fp_corr_manifest <- readr::read_csv(file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))

  # tf_corr_qc(fp_corr_manifest)
  # threshold_fp_tf_corr_p <- 0.3
  # combined <- tf_corr_footprints_filter(fp_corr_manifest, p_thr = threshold_fp_tf_corr_p, r_thr = threshold_fp_tf_corr_r, output_bed = file.path(base_dir, sprintf("fp_predicted_tfbs_%s", db)))

  # fp_score <- combined$fp_score
  # fp_bound <- combined$fp_bound
  # fp_annotation <- combined$fp_annotation

  # readr::write_csv(fp_score, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_score_strict_tf_filtered_corr_", db, ".csv"))
  # readr::write_csv(fp_bound, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_bound_strict_tf_filtered_corr_", db, ".csv"))
  # readr::write_csv(fp_annotation, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_annotation_strict_tf_filtered_corr_", db, ".csv"))

  # unique(fp_annotation$fp_peak)

  fp_score <- readr::read_csv(file.path(base_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db)))
  fp_bound <- readr::read_csv(file.path(base_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db)))
  fp_annotation <- readr::read_csv(file.path(base_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db)))


  # grn_set
  grn_set <- build_grn_set(
    fp_score      = fp_score,
    fp_bound      = fp_bound,
    fp_annotation = fp_annotation,
    atac_score    = atac_score,
    atac_overlap  = atac_overlap,
    rna           = rna,
    metadata      = sample_metadata,
    tf_list       = tf_list,
    motif_db      = motif_db,
    label_col     = "Sample",
    expected_n    = 13
  )
  grn_set <- filter_by_min_bound(grn_set, min_bound = 1L)
}

# (Optional correlate with all tfbs, not only canonical tfbs)
# tf_vec <- strict_rna |> dplyr::filter(HGNC == "ZEB1") |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
# res <- tf_corr_footprints_all_tfbs(
#   fp_manifest = fp_aligned_normalized_filtered_manifest,
#   tf_name  = "ZEB1",
#   tf_expr  = tf_vec,
#   out_dir  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
# )
# res$files



# Step 2. Connect TF-occupied enhancers to target genes -------------------
if (do_tf_to_target_genes_prediction == TRUE) {
  gene_annot_ref_hg38 <- episcope_build_gene_annot("hg38")
  # gene_annot_ref_mm10 <- episcope_build_gene_annot("mm10")


  # atac_peak to target gene corr

  # GeneHancer based
  # gh_std   <- load_genehancer_panc(file.path("inst","extdata","GeneHancer_v5.24_full.csv")) # GeneHancer_v5.24_elite_panc.csv

  # 30kb/25kb/20kb based
  gh_std <- episcope_make_windowed_gh(
    peaks      = grn_set$atac_score,
    genes      = unique(c(grn_set$rna$HGNC, grn_set$rna$ensembl_gene_id)),
    flank_bp   = 100000, # 100kb
    mode       = "TSS",
    gene_annot = gene_annot_ref_hg38,
    id_col     = "HGNC" # or "ensembl_gene_id"
  )

  options(future.globals.maxSize = 64 * 1024^3)

  # Run atac gene correlations (adjusted p-value < 0.05, |r| >= 0.3)
  threshold_atac_gene_corr_p <- 0.3
  # atac_res <- correlate_atac_to_genes(grn_set = grn_set, gh_tbl = gh_std, gene_mode = "both", fdr = threshold_atac_gene_corr_p, r_abs_min = threshold_atac_gene_corr_abs_r)

  # shrink to only correlated peaks
  # grn_set_filtered <- filter_grn_by_corr(grn_set, atac_res$atac_gene_corr_kept)

  # Run fp gene correlations (adjusted p-value < 0.05, |r| >= 0.3)
  threshold_fp_gene_corr_p <- 0.3
  # fp_res <- correlate_fp_to_genes(
  #   grn_set              = grn_set_filtered,
  #   atac_gene_corr_kept  = atac_res$atac_gene_corr_kept,
  #   fdr                  = threshold_fp_gene_corr_p,
  #   r_abs_min            = threshold_fp_gene_corr_abs_r,
  #   workers              = 15
  # )


  atac_res <- correlate_atac_to_genes(
    grn_set  = grn_set,
    gh_tbl   = gh_std,
    gene_mode = "both",
    fdr       = threshold_atac_gene_corr_p,
    r_abs_min = threshold_atac_gene_corr_abs_r,
    cache_dir = "/data/homes/yl814/episcope_cache/atac_gene_corr",
    cache_tag = "GSE87218_100kb",
    workers   = 15,
    cache_verbose = TRUE
  )
  fp_res <- correlate_fp_to_genes(
    grn_set             = grn_set,
    atac_gene_corr_kept = atac_res$atac_gene_corr_full,
    fdr                 = threshold_fp_gene_corr_p,
    r_abs_min           = threshold_fp_gene_corr_abs_r,
    method              = "pearson",
    workers             = 20,
    cache_dir           = "/data/homes/yl814/episcope_cache/fp_gene_corr",
    cache_tag           = "GSE87218_100kb",
    cache_chunk_size    = 5000L,
    cache_verbose       = TRUE
  )
  readr::write_csv(atac_res$atac_gene_corr_full, file.path(base_dir, sprintf("atac_gene_corr_full_100kb_%s.csv", db)))
  readr::write_csv(fp_res$fp_gene_corr_full, file.path(base_dir, sprintf("fp_gene_corr_full_100kb_%s.csv", db)))
}


if (FALSE) {
  # match_mode <- "strict"  # strict lenient
  regulated_genes <- 1.5 # 1.5 2
  delta_link <- 1 # 2 1
  regulation_priors <- "300kb" # "genehancer"
  lighting_folder <- file.path(base_dir, paste("lighting", "fp_tf_corr_FDR", threshold_fp_tf_corr_p, regulation_priors, db, "regulated_genes", regulated_genes, "delta_link", delta_link, sep = "_"))
  print(lighting_folder)


  # Step 3. Build basal GRN & identify active regulatory edges per c --------
  # basal
  basal <- make_basal_links(
    fp_gene_corr_kept = fp_res$fp_gene_corr_kept,
    fp_annotation     = grn_set$fp_annotation,
    out_dir           = lighting_folder,
    prefix            = "lighting"
  )


  light_by_condition(
    ds = grn_set,
    basal_links = basal,
    out_dir = lighting_folder,
    prefix = "lighting",
    label_col = "Sample",
    link_score_threshold = link_score_threshold,
    fp_score_threshold = fp_score_threshold,
    tf_expr_threshold = threshold_tf_expr,
    use_parallel = TRUE,
    workers = 4
  )

  # Step 4. Perform differential GRN analysis & identify master TFs ---------

  specs <- build_cellwise_contrasts_from_index(
    index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
    out_dir = lighting_folder,
    prefix = "lighting",
    ctrl_tag = "Ctrl",
    clean_names = FALSE
  )
  str(specs)


  # Adding comparisons vs. culture:
  # --- Add treatment-vs-culture matched-time contrasts -------------------------

  # map the time token in names like "BG_1h_post" -> target culture name
  time_to_culture <- c(
    "1h_post"        = "culture_1h_post",
    "4h_post"        = "culture_4h_post",
    "24h_post"       = "culture_24h_post",
    "24h_culture_5d" = "culture_6d"
  )

  # lookup: condition name -> its per-condition source path (from existing specs)
  src_lookup <- setNames(specs$cond1_source, specs$cond1_name)

  # base output folder (reuse what specs already uses)
  out_base <- dirname(specs$out_file[1])

  # pick treatments that need culture controls
  treat_rows <- specs[grepl("^(BG|LPS)_", specs$cond1_name) & specs$cond2_name == "Ctrl", , drop = FALSE]

  # derive the matched culture name for each treatment
  time_token <- sub("^[^_]+_", "", treat_rows$cond1_name) # e.g., "1h_post", "24h_culture_5d"
  matched_culture <- unname(time_to_culture[time_token])

  # keep only rows that have a valid mapping and available culture sources
  ok <- !is.na(matched_culture) & matched_culture %in% names(src_lookup) & treat_rows$cond1_name %in% names(src_lookup)

  new_specs <- tibble::tibble(
    cond1_label  = treat_rows$cond1_label[ok],
    cond2_label  = matched_culture[ok],
    cond1_source = unname(src_lookup[treat_rows$cond1_name[ok]]),
    cond2_source = unname(src_lookup[matched_culture[ok]]),
    cond1_name   = treat_rows$cond1_name[ok],
    cond2_name   = matched_culture[ok],
    out_file     = file.path(out_base, paste0(treat_rows$cond1_name[ok], "_vs_", matched_culture[ok], "_delta_links.csv"))
  )

  # append, de-duplicate if needed
  specs2 <- rbind(specs, new_specs)
  specs2 <- specs2[!duplicated(specs2[, c("cond1_name", "cond2_name")]), , drop = FALSE]

  specs <- specs2
  rm(specs2, new_specs, treat_rows, time_token, matched_culture, ok, src_lookup, out_base, time_to_culture)

  # quick check
  specs[, c("cond1_name", "cond2_name", "out_file")]


  # run_links_deltas_driver(
  #   out_dir = lighting_folder,
  #   prefix  = "lighting",
  #   index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
  #   ctrl_tag = "Ctrl",
  #   clean_names = FALSE,
  #   parallel = TRUE
  # )

  run_links_deltas_driver(
    specs       = specs,
    clean_names = FALSE,
    parallel    = TRUE
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

  bulk$filtered_paths
  bulk$filtered_dirs
  bulk$manifest_path

  # LDA across all filtered CSVs produced above
  filtered_csvs <- bulk$filtered_paths

  lda_out <- episcope::annotate_links_with_lda_topic_bulk(
    filtered_csvs,
    K = 20, # topics
    which = "abs", # |delta_link_score|
    min_df = 2, # drop ultra-rare terms
    gamma_cutoff = 0.2, # keep all topics with gamma >= 0.2 per TF
    parallel = TRUE, workers = 20, plan = "multisession", seed = 1
  )

  lda_out$assigned_paths # vector of "*_filtered_lda_K20.csv"
  lda_out$assigned_dirs # unique directories
  lda_out$manifest_path # CSV mapping filtered -> annotated

  # Summarize all annotated CSVs produced by LDA
  annotated_csvs <- lda_out$assigned_paths

  sum_out <- episcope::summarize_lda_annotations_bulk(
    annotated_csvs,
    edge_filter_min = delta_link, # keep edges used to call delta links
    parallel = TRUE, plan = "multisession", workers = 20,
    use_tf_r_weight = TRUE # optional: weight deltas by TF r column if present
  )

  sum_out$summary_paths # "*_filtered_lda_K20_summary.csv"
  sum_out$summary_dirs
  sum_out$manifest_path


  # Using outputs from summarize_lda_annotations_bulk():
  summary_csvs <- sum_out$summary_paths

  # Generate all bubble & hub PDFs (overall / activate / repress + no-cancel variants)
  episcope::plot_from_summary_bulk(
    summary_csvs,
    top_tf_per_topic = 20,
    parallel = TRUE,
    workers = 20,
    plan = "multisession",
    color_sigma = 2
  )


  # TODO TF HITS score summary plot (waterfall)

  # Step 5. Generate interactive Topic & TF regulatory hub subnetwork
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
    motif_db = "jaspar2024",
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
            motif_db = "jaspar2024",
            verbose = TRUE
          )
          TRUE
        },
        silent = TRUE
      )
    }
  }

}






