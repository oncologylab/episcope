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
do_load_footprints_preprocess    <- FALSE
do_tf_binding_sites_prediction   <- FALSE
do_tf_to_target_genes_prediction <- FALSE
verbose <- TRUE

# Load footprint data and preprocess --------------------------------------
if (do_load_footprints_preprocess == TRUE) {
  gene_symbol_col <- if (exists("gene_symbol_col")) gene_symbol_col else "HGNC"
  omics_data <- load_multiomic_data(
    config = "dev/config/pdac_nutrient_stress_strict_jaspar2024.yaml",
    genome = ref_genome,
    gene_symbol_col = gene_symbol_col,
    label_col = "strict_match_rna",
    expected_n = expected_n,
    tf_list = tf_list,
    motif_db = motif_db,
    threshold_gene_expr = threshold_gene_expr,
    threshold_fp_score = threshold_fp_score,
    do_preprocess = TRUE,
    do_motif_clustering = TRUE,
    output_mode = "distinct",
    write_outputs = TRUE,
    use_parallel = TRUE,
    verbose = verbose
  )
}

# Predict TF binding sites ------------------------------------------------
if (do_tf_binding_sites_prediction == TRUE) {
  omics_data <- correlate_tf_to_fp(
    omics_data = omics_data,
    mode = "canonical",
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
step3_root_dir <- file.path(base_dir, "diff_grn_and_regulatory_topics")

diff_res <- find_differential_links(
  config = "dev/config/pdac_nutrient_stress_strict_jaspar2024.yaml",
  compar = file.path(base_dir, "data", "episcope_comparisons.csv"),
  output_dir = step3_root_dir
)

motif_path <- file.path("inst", "extdata", "genome", "JASPAR2024.txt")
motif_db$gene_symbol <- motif_db$HGNC
if (is.data.frame(motif_db) && all(c("sub_cluster_name", "gene_symbol") %in% names(motif_db))) {
  motif_info <- build_tf_cluster_map_from_motif(motif_db)
} else {
  motif_info <- build_tf_cluster_map_from_motif(motif_path)
}

# k_grid_default <- c(2:15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100)
k_grid_default <- c(2:15, 20, 25, 30)

celllines <- c("HPAFII")
selected_k <- 10L

gene_term_modes <- c("aggregate", "unique")
include_tf_terms_opts <- c(FALSE, TRUE)
count_inputs <- c("weight", "pseudo_count_bin", "pseudo_count_log")
backends <- c("warplda", "vae")
vae_variants <- c("multivi_encoder")
library(enrichR)
for (gene_term_mode in gene_term_modes) {
  for (include_tf_terms in include_tf_terms_opts) {
    for (count_input in count_inputs) {
      for (backend in backends) {
        variant_list <- if (backend == "vae") vae_variants else "multivi_encoder"
        for (vae_variant in variant_list) {
          combo_tag <- paste0(
            "gene_", gene_term_mode,
            "_tf_", if (isTRUE(include_tf_terms)) "on" else "off",
            "_count_", count_input
          )
          combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_tag)

          topic_model_dir <- file.path(step3_root_dir, "topic_models", combo_tag)
          topic_final_dir <- file.path(step3_root_dir, "final_topics", combo_tag)

          topic_args <- list(
            pathway_source = "link_scores",
            thrP = 0.9,
            pathway_make_heatmap = FALSE,
            pathway_make_dotplot = TRUE,
            pathway_overwrite = TRUE,
            pathway_per_comparison = TRUE,
            pathway_per_comparison_dir = "per_cmpr_topic_pathway",
            pathway_split_direction = TRUE,
            run_pathway_gsea = FALSE,
            run_link_topic_scores = TRUE,
            link_topic_gate_mode = c("none", "peak_and_gene_in_set"),
            link_topic_overwrite = TRUE,
            pathway_link_scores_file = "link_topic_scores_gate_peak_and_gene_in_set.csv",
            pathway_link_scores_file_tf = "link_topic_scores_gate_peak_and_gene_in_set.csv",
            pathway_link_gene_terms_file = "topic_terms.csv",
            pathway_link_min_prob = 0,
            pathway_link_include_tf = TRUE,
            pathway_link_include_gene = TRUE,
            pathway_link_gene_min_prob = 0,
            pathway_link_tf_min_prob = 0.5,
            pathway_link_tf_max_topics = 5L,
            pathway_link_tf_top_n_per_topic = 30L
          )

          train_topic_models(
            Kgrid = k_grid_default,
            input_dir = diff_res$filtered_dir,
            output_dir = topic_model_dir,
            celllines = celllines,
            tf_cluster_map = motif_info$tf_cluster_map,
            tf_exclude = motif_info$tf_exclude,
            gene_term_mode = gene_term_mode,
            include_tf_terms = include_tf_terms,
            count_input = count_input,
            backend = backend,
            vae_variant = vae_variant,
            topic_report_args = topic_args
          )

          extract_regulatory_topics(
            k = selected_k,
            model_dir = topic_model_dir,
            output_dir = topic_final_dir,
            backend = backend,
            vae_variant = vae_variant,
            topic_report_args = topic_args
          )

          run_vae_topic_delta_network_plots(
            topic_root = topic_final_dir,
            step2_out_dir = diff_res$filtered_dir
          )
        }
      }
    }
  }
}

run_single_topic_test <- TRUE
if (isTRUE(run_single_topic_test)) {
  single_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  single_backend <- "vae"
  single_variant <- "multivi_encoder"
  single_k <- 12L
  single_model_dir <- file.path(step3_root_dir, "topic_models", single_combo_tag)
  single_final_dir <- file.path(step3_root_dir, "final_topics", single_combo_tag)

  extract_regulatory_topics(
    k = single_k,
    model_dir = single_model_dir,
    output_dir = single_final_dir,
    backend = single_backend,
    vae_variant = single_variant,
    topic_report_args = list(
      pathway_source = "link_scores",
      thrP = 0.9,
      pathway_make_heatmap = FALSE,
      pathway_make_dotplot = TRUE,
      pathway_per_comparison = TRUE,
      pathway_per_comparison_dir = "per_cmpr_topic_pathway",
      pathway_split_direction = TRUE,
      run_pathway_gsea = FALSE,
      run_link_topic_scores = TRUE,
      link_topic_gate_mode = c("none", "peak_and_gene_in_set"),
      link_topic_overwrite = TRUE,
      pathway_link_scores_file = "link_topic_scores_gate_peak_and_gene_in_set.csv",
      pathway_link_scores_file_tf = "link_topic_scores_gate_peak_and_gene_in_set.csv",
      pathway_link_gene_terms_file = "topic_terms.csv",
      pathway_link_min_prob = 0,
      pathway_link_include_tf = TRUE,
      pathway_link_include_gene = TRUE,
      pathway_link_gene_min_prob = 0,
      pathway_link_tf_min_prob = 0.5,
      pathway_link_tf_max_topics = 5L,
      pathway_link_tf_top_n_per_topic = 30L
    )
  )
}

run_gammafit_preview <- TRUE
if (isTRUE(run_gammafit_preview)) {
  preview_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  preview_backend <- "vae"
  preview_variant <- "multivi_encoder"
  preview_k <- 12L
  preview_model_dir <- file.path(step3_root_dir, "topic_models", preview_combo_tag)
  preview_final_dir <- file.path(step3_root_dir, "final_topics", preview_combo_tag)

  extract_regulatory_topics(
    k = preview_k,
    model_dir = preview_model_dir,
    output_dir = preview_final_dir,
    backend = preview_backend,
    vae_variant = preview_variant,
    topic_report_args = list(
      binarize_method = "gammafit",
      thrP = 0.9,
      pathway_make_heatmap = FALSE,
      pathway_make_dotplot = FALSE,
      run_pathway_gsea = FALSE,
      run_link_topic_scores = FALSE,
      link_topic_overwrite = TRUE
    )
  )
}

run_doc_topic_summary_preview <- TRUE
if (isTRUE(run_doc_topic_summary_preview)) {
  preview_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  preview_backend <- "vae"
  preview_variant <- "multivi_encoder"
  preview_k <- 12L
  preview_final_dir <- file.path(step3_root_dir, "final_topics", preview_combo_tag)

  out_dirs <- list.dirs(preview_final_dir, recursive = FALSE, full.names = TRUE)
  if (preview_backend == "vae") {
    patt <- paste0("_vae_joint_ctf_docs_peak_delta_fp_gene_fc_expr_", preview_variant, "_K", preview_k, "$")
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_ctf_docs_peak_delta_fp_gene_fc_expr_warplda_K", preview_k, "$"), basename(out_dirs))]
  }
  if (length(out_dirs)) {
    link_path <- file.path(out_dirs[1], "link_topic_scores_baseline.csv")
    if (file.exists(link_path)) {
      link_dt <- data.table::fread(link_path)
      plot_topic_delta_networks_from_link_scores(
        link_scores = link_dt,
        step2_out_dir = diff_res$filtered_dir,
        out_root = file.path(out_dirs[1], "doc_topic_sub_network_link_scores_baseline"),
        min_prob = 0.5,
        filter_same_direction = TRUE
      )
    } else {
      .log_warn("No link_topic_scores_baseline.csv found for preview.")
    }
  } else {
    .log_warn("No final topic directory found for preview.")
  }
}
