library(episcope)

source("R/utils_config.R")
source("R/utils_logging.R")
source("R/utils_helpers.R")

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
# Testing only
tf_binding_mode <- "all"
omics_data_rds_path <- "/data/homes/yl814/episcope_test/nutrient_stress_strict_JASPAR2024/predict_tf_binding_sites_all/01_multiomic_data_object_JASPAR2024.rds"

expected_n <- if (exists("expected_n")) expected_n else NULL
# ──────────────────────────────────────────────────────────────────────────────
# Turn modules ON/OFF
# ──────────────────────────────────────────────────────────────────────────────
do_load_footprints_preprocess    <- FALSE
do_tf_binding_sites_prediction   <- FALSE
do_tf_to_target_genes_prediction <- FALSE
do_diff_grn <- TRUE
do_topic_combo_run <- FALSE

verbose <- TRUE

# Load footprint data and preprocess --------------------------------------
if (do_load_footprints_preprocess == TRUE) {
  gene_symbol_col <- if (exists("gene_symbol_col")) gene_symbol_col else "HGNC"
  omics_data <- load_prep_multiomic_data(
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
  tf_binding_mode_use <- if (exists("tf_binding_mode")) as.character(tf_binding_mode) else "canonical"
  tf_binding_mode_use <- match.arg(tf_binding_mode_use, c("canonical", "all"))
  step1_tfbs_out_dir <- file.path(base_dir, paste0("predict_tf_binding_sites_", tf_binding_mode_use))
  omics_data <- correlate_tf_to_fp(
    omics_data = omics_data,
    mode = tf_binding_mode_use,
    out_dir = step1_tfbs_out_dir,
    label_col = "strict_match_rna",
    r_thr = threshold_fp_tf_corr_r,
    p_thr = threshold_fp_tf_corr_p,
    db = db,
    chunk_size = 5000L,
    all_mode_tf_chunk_size = if (exists("all_mode_tf_chunk_size")) as.integer(all_mode_tf_chunk_size) else 5L,
    min_non_na = 5L,
    qc = TRUE,
    write_bed = FALSE,
    use_cache = TRUE
  )
}

# Connect TFs to target genes ---------------------------------------------
if (do_tf_to_target_genes_prediction == TRUE) {
  tf_binding_mode_step2 <- if (exists("tf_binding_mode")) as.character(tf_binding_mode) else "canonical"
  tf_binding_mode_step2 <- match.arg(tf_binding_mode_step2, c("canonical", "all"))
  if (!exists("omics_data") || is.null(omics_data)) {
    omics_data_rds <- if (exists("omics_data_rds_path")) {
      omics_data_rds_path
    } else {
      file.path(
        base_dir,
        paste0("predict_tf_binding_sites_", tf_binding_mode_step2),
        sprintf("01_multiomic_data_object_%s.rds", db)
      )
    }
    if (!file.exists(omics_data_rds)) {
      .log_abort(
        paste0(
          "Step 2 requires `omics_data`. Expected RDS not found at: ",
          omics_data_rds,
          ". Set `omics_data_rds_path` to a valid file."
        )
      )
    }
    .log_inform(sprintf("Loading omics_data from %s", omics_data_rds))
    omics_data <- readRDS(omics_data_rds)
  }
  if (!is.data.frame(omics_data$rna_condition)) {
    .log_inform("`omics_data$rna_condition` missing; building from `omics_data$rna` using label_col = strict_match_rna.")
    omics_data <- grn_add_rna_condition(
      grn_set = omics_data,
      label_col = "strict_match_rna",
      verbose = TRUE
    )
  }
  if (!("tfs" %in% names(omics_data$fp_annotation) || "TF" %in% names(omics_data$fp_annotation))) {
    step1_tfbs_out_dir <- file.path(base_dir, paste0("predict_tf_binding_sites_", tf_binding_mode_step2))
    fp_ann_tf_path <- file.path(
      step1_tfbs_out_dir,
      "cache",
      sprintf("fp_annotation_strict_tf_filtered_corr_%s_%s.csv", db, tf_binding_mode_step2)
    )
    if (!file.exists(fp_ann_tf_path)) {
      .log_abort(
        paste0(
          "`omics_data$fp_annotation` has no TF column (tfs/TF), and TF-annotated file was not found at: ",
          fp_ann_tf_path
        )
      )
    }
    .log_inform(sprintf("Loading TF-annotated fp_annotation from %s", fp_ann_tf_path))
    omics_data$fp_annotation <- readr::read_csv(fp_ann_tf_path, show_col_types = FALSE)
  }
  step2_out_dir <- if (identical(tf_binding_mode_step2, "all")) {
    file.path(base_dir, "connect_tf_target_genes_all")
  } else {
    file.path(base_dir, "connect_tf_target_genes")
  }
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
  options(future.globals.maxSize = 196 * 1024^3)
  fp_gene_workers <- if (exists("fp_gene_workers")) as.integer(fp_gene_workers) else {
    if (identical(tf_binding_mode_step2, "all")) 8L else 20L
  }
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
    workers          = fp_gene_workers,
    cache_dir        = file.path(base_dir, "cache", "fp_gene_corr"),
    cache_tag        = paste0("nutrient_", tf_binding_mode_step2, "_", link_mode, "_pearson"),
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
    workers          = fp_gene_workers,
    cache_dir        = file.path(base_dir, "cache", "fp_gene_corr"),
    cache_tag        = paste0("nutrient_", tf_binding_mode_step2, "_", link_mode, "_spearman"),
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
    summary_by_tf = link_summary$summary_by_tf,
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

# Diff GRN and regulatory topics ------------------------------------------
step3_root_dir <- file.path(base_dir, "diff_grn_and_regulatory_topics")

if (do_diff_grn == TRUE) {
  .fmt_cut <- function(x) gsub("\\.", "p", as.character(signif(as.numeric(x), 4)))
  diff_modes <- c("delta", "log2fc")
  diff_res_by_mode <- list()

  for (mode_use in diff_modes) {
    fp_filter_mode <- mode_use
    mode_cutoff <- if (identical(mode_use, "delta")) fp_delta_cutoff else fp_log2fc_cutoff
    step3_out_mode <- file.path(
      base_dir,
      paste0(
        "diff_grn_and_regulatory_topics",
        "_fp_", mode_use,
        "_cutoff_", .fmt_cut(mode_cutoff),
        "_gene_log2fc_", .fmt_cut(gene_log2fc_cutoff)
      )
    )

    .log_inform(
      "Running differential links with fp_filter_mode={mode_use}, cutoff={mode_cutoff}, gene_log2fc_cutoff={gene_log2fc_cutoff} -> {step3_out_mode}"
    )
    diff_res_by_mode[[mode_use]] <- find_differential_links(
      config = NULL,
      compar = file.path(base_dir, "data", "episcope_comparisons.csv"),
      output_dir = step3_out_mode,
      overwrite_delta = FALSE,
      overwrite_filtered = FALSE,
      overwrite_tf_hubs = TRUE,
      connectivity_min_degree = 5L,
      summary_plot_format = "both"
    )
  }

  diff_res <- diff_res_by_mode[["delta"]]

  # Optional: pathway enrichment from diff_links_filtered (minimal gene_key-based run).
  do_diff_links_pathway_enrichment <- TRUE
  if (isTRUE(do_diff_links_pathway_enrichment)) {
    for (mode_use in names(diff_res_by_mode)) {
      run_diff_links_pathway_grn(
        diff_res = diff_res_by_mode[[mode_use]],
        min_genes = 5L,
        padj_cut = 0.05,
        plot_subnetwork = TRUE,
        top_n_pathways_plot = 10L,
        gene_log2fc_cutoff = gene_log2fc_cutoff,
        overwrite = TRUE,
        verbose = TRUE
      )
    }
  }

  motif_path <- file.path("inst", "extdata", "genome", "JASPAR2024.txt")
  # motif_db$gene_symbol <- motif_db$HGNC
  if (is.data.frame(motif_db) && all(c("sub_cluster_name", "gene_symbol") %in% names(motif_db))) {
    motif_info <- build_tf_cluster_map_from_motif(motif_db)
  } else {
    motif_info <- build_tf_cluster_map_from_motif(motif_path)
  }

}



if (isTRUE(do_topic_combo_run)) {
  # k_grid_default <- c(2:15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100)
  k_grid_default <- c(2:15, 20, 25, 30)

  celllines <- c("HPAFII")
  selected_k <- 10L

  gene_term_modes <- c("aggregate", "unique")
  include_tf_terms_opts <- c(FALSE, TRUE)
  count_inputs <- c("weight", "pseudo_count_bin", "pseudo_count_log")
  backends <- c("vae","warplda")
  vae_variants <- c("multivi_encoder")
  # Topic document mode: "tf_cluster" (current ctf_docs behavior) or "tf" (new tf_docs behavior).
  topic_doc_mode <- "tf_cluster"
  library(enrichR)
  library(future)
  library(future.apply)

  combo_grid <- expand.grid(
    gene_term_mode = gene_term_modes,
    include_tf_terms = include_tf_terms_opts,
    count_input = count_inputs,
    backend = backends,
    stringsAsFactors = FALSE
  )
  combo_grid$vae_variant <- ifelse(combo_grid$backend == "vae", vae_variants[1], "multivi_encoder")
  # Minimal single-combo run (set FALSE to restore full grid).
  run_single_combo_only <- TRUE
  if (isTRUE(run_single_combo_only)) {
    combo_grid <- combo_grid[
      combo_grid$gene_term_mode == "aggregate" &
        !combo_grid$include_tf_terms &
        combo_grid$count_input == "pseudo_count_log" &
        combo_grid$backend == "vae",
      ,
      drop = FALSE
    ]
  }
  thrP_modes <- c(0.8, 0.9)
  gene_prob_cutoff_modes <- c("0.3", "0.5", "max")

  # Fast rerun mode: skip training/extraction and regenerate only doc_topic_sub_network_link_* plots.
  rerun_subnetwork_only <- TRUE
  rerun_subnetwork_methods <- c("peak_and_gene", "peak_and_gene_prob")
  rerun_subnetwork_min_prob <- 0.5
  rerun_subnetwork_filter_same_direction <- FALSE
  rerun_subnetwork_min_prob_mode <- "match_link_cutoff"  # "fixed" or "match_link_cutoff"

  # Post-extraction optional modules (can be run in rerun-only mode).
  run_subnetwork_plots <- FALSE
  run_pathway_rerun <- TRUE
  run_pathway_per_comparison <- TRUE
  run_pathway_split_direction <- TRUE

  benchmark_dirname <- "benchmark_gammafit_thrP_vae_0.9_peak_and_gene_vs_peak_and_gene_prob"
  combo_error_log <- file.path(step3_root_dir, "benchmark", "topic_combo_errors.log")
  dir.create(dirname(combo_error_log), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(combo_error_log)) file.remove(combo_error_log)
  combo_failures <- vector("list", 0L)

  .append_combo_error_log <- function(line) {
    dir.create(dirname(combo_error_log), recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      {
        con <- file(combo_error_log, open = "a")
        on.exit(close(con), add = TRUE)
        writeLines(line, con = con, sep = "\n")
      },
      error = function(e) {
        .log_warn("Could not write combo error log: {conditionMessage(e)}")
        .log_warn("Original combo error: {line}")
        invisible(NULL)
      }
    )
  }

  extract_modules <- list(
    pathway = TRUE,
    doc_topic_heatmaps = TRUE,
    topic_by_comparison = TRUE,
    topic_marker_heatmap = TRUE,
    intertopic_distance = TRUE,
    ldavis = TRUE
  )
  extract_overwrite <- list(
    link_topic = TRUE,
    pathway = TRUE
  )

  .log_inform(
    "Starting topic-model combo run: {nrow(combo_grid)} combinations x {length(thrP_modes)} thrP modes x {length(gene_prob_cutoff_modes)} gene_prob cutoff modes."
  )
  for (i in seq_len(nrow(combo_grid))) {
    gene_term_mode <- combo_grid$gene_term_mode[i]
    include_tf_terms <- combo_grid$include_tf_terms[i]
    count_input <- combo_grid$count_input[i]
    backend <- combo_grid$backend[i]
    vae_variant <- combo_grid$vae_variant[i]
    link_fdr_p_use <- if (identical(backend, "warplda")) 0.01 else 0.1
    dataset_tag <- if (length(celllines)) as.character(celllines[[1]]) else "dataset"
    thrP_grid_this_combo <- if (identical(backend, "warplda")) 0.5 else thrP_modes

    combo_tag <- paste0(
      "gene_", gene_term_mode,
      "_tf_", if (isTRUE(include_tf_terms)) "on" else "off",
      "_count_", count_input
    )
    combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_tag)
    for (thrP_use in thrP_grid_this_combo) {
      for (cut_mode in gene_prob_cutoff_modes) {
        compact_run_dirname <- build_topic_compact_run_dirname(
          thrP_use = thrP_use,
          cut_mode = cut_mode,
          gene_term_mode = gene_term_mode,
          include_tf_terms = include_tf_terms,
          count_input = count_input,
          dataset_tag = dataset_tag,
          doc_mode = topic_doc_mode,
          backend = backend,
          vae_variant = vae_variant,
          k_use = selected_k
        )
        prob_cutoff_use <- if (identical(cut_mode, "max")) "max" else as.numeric(cut_mode)
        subnet_min_prob_use <- rerun_subnetwork_min_prob
        if (identical(rerun_subnetwork_min_prob_mode, "match_link_cutoff")) {
          subnet_min_prob_use <- if (identical(cut_mode, "max")) 0 else as.numeric(cut_mode)
        }
        .log_inform(
          "[{i}/{nrow(combo_grid)}] Running combo={combo_tag} backend={backend} (thrP={thrP_use}, link_method=gene_prob, gene_prob_cutoff={cut_mode}) -> {compact_run_dirname}"
        )

        topic_model_dir <- file.path(step3_root_dir, "topic_models", combo_tag)
        topic_final_dir <- file.path(step3_root_dir, compact_run_dirname)
        write_topic_directory_name_readme(
          topic_final_dir,
          fields = list(
            thrP = thrP_use,
            link_cutoff_mode = cut_mode,
            link_method = "gene_prob",
            gene_term_mode = gene_term_mode,
            include_tf_terms = include_tf_terms,
            count_input = count_input,
            doc_mode = topic_doc_mode,
            backend = backend,
            vae_variant = vae_variant,
            selected_k = selected_k,
            dataset_tag = dataset_tag
          )
        )

        topic_args <- make_topic_report_args_simple(
          thrP = thrP_use,
          link_prob_cutoff = prob_cutoff_use,
          link_fdr_p = link_fdr_p_use,
          modules = extract_modules,
          overwrite = extract_overwrite
        )

        combo_ok <- tryCatch(
          {
            if (!isTRUE(rerun_subnetwork_only)) {
              train_topic_models(
                Kgrid = k_grid_default,
                input_dir = diff_res$filtered_dir,
                output_dir = topic_model_dir,
                celllines = celllines,
                tf_cluster_map = motif_info$tf_cluster_map,
                doc_mode = topic_doc_mode,
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
                doc_mode = topic_doc_mode,
                flatten_single_output = TRUE,
                topic_report_args = topic_args
              )

              run_vae_doc_topic_heatmaps(
                topic_root = topic_final_dir,
                backend = backend,
                vae_variant = vae_variant,
                doc_mode = topic_doc_mode
              )
            } else {
              .log_inform(
                "[{i}/{nrow(combo_grid)}] rerun_subnetwork_only=TRUE; skipping model train/extract/heatmap and plotting methods={paste(rerun_subnetwork_methods, collapse = ',')}."
              )
            }

            if (isTRUE(run_subnetwork_plots)) {
              run_vae_topic_delta_network_plots(
                topic_root = topic_final_dir,
                step2_out_dir = diff_res$filtered_dir,
                min_prob = subnet_min_prob_use,
                filter_same_direction = rerun_subnetwork_filter_same_direction,
                methods = rerun_subnetwork_methods,
                backend = backend,
                vae_variant = vae_variant,
                doc_mode = topic_doc_mode
              )
            } else {
              .log_inform(
                "[{i}/{nrow(combo_grid)}] run_subnetwork_plots=FALSE; skipping subnetwork plotting."
              )
            }

            if (isTRUE(run_pathway_rerun)) {
              run_vae_topic_delta_network_pathway(
                topic_root = topic_final_dir,
                backend = backend,
                vae_variant = vae_variant,
                doc_mode = topic_doc_mode,
                top_n_per_topic = topic_args$top_n_per_topic,
                max_pathways = topic_args$max_pathways,
                per_comparison = isTRUE(run_pathway_per_comparison),
                split_direction = isTRUE(run_pathway_split_direction)
              )
            } else {
              .log_inform(
                "[{i}/{nrow(combo_grid)}] run_pathway_rerun=FALSE; skipping pathway rerun."
              )
            }
            TRUE
          },
          error = function(e) {
            err_msg <- conditionMessage(e)
            combo_failures[[length(combo_failures) + 1L]] <<- data.frame(
              row = i,
              combo_tag = as.character(combo_tag),
              backend = as.character(backend),
              gene_term_mode = as.character(gene_term_mode),
              include_tf_terms = as.logical(include_tf_terms),
              count_input = as.character(count_input),
              error = as.character(err_msg),
              stringsAsFactors = FALSE
            )
            stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
            .append_combo_error_log(sprintf(
              "[%s]\trow=%d\tcombo=%s\tbackend=%s\tgene_term_mode=%s\tinclude_tf_terms=%s\tcount_input=%s\tthrP=%s\tcutoff=%s\terror=%s",
              stamp, i, combo_tag, backend, gene_term_mode, include_tf_terms, count_input, thrP_use, cut_mode, err_msg
            ))
            .log_warn("Combo failed [{i}/{nrow(combo_grid)}]: {combo_tag} ({backend}, thrP={thrP_use}, cutoff={cut_mode}) -> {err_msg}")
            FALSE
          }
        )
        if (isTRUE(combo_ok)) {
          .log_inform("[{i}/{nrow(combo_grid)}] Completed combo={combo_tag} backend={backend}, thrP={thrP_use}, cutoff={cut_mode}.")
        }
      }
    }
  }

  summarize_topic_combo_failures(
    combo_failures = combo_failures,
    combo_error_log = combo_error_log
  )
} else if (isTRUE(verbose)) {
  .log_inform("Skipping topic-model combo run (do_topic_combo_run=FALSE).")
}


run_single_topic_test <- FALSE
if (isTRUE(do_topic_combo_run) && isTRUE(run_single_topic_test)) {
  # Pick one combo_grid row for an end-to-end troubleshooting run.
  # Example: 17 = aggregate/tf_off/pseudo_count_bin/vae
  single_combo_idx <- 17L
  if (!is.numeric(single_combo_idx) || single_combo_idx < 1L || single_combo_idx > nrow(combo_grid)) {
    .log_abort("single_combo_idx must be between 1 and {nrow(combo_grid)}.")
  }
  single_row <- combo_grid[single_combo_idx, , drop = FALSE]
  single_backend <- as.character(single_row$backend[[1]])
  single_variant <- as.character(single_row$vae_variant[[1]])
  single_combo_tag <- paste0(
    "gene_", as.character(single_row$gene_term_mode[[1]]),
    "_tf_", if (isTRUE(single_row$include_tf_terms[[1]])) "on" else "off",
    "_count_", as.character(single_row$count_input[[1]])
  )
  single_combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", single_combo_tag)
  single_k <- as.integer(selected_k)
  single_thrP_use <- if (identical(single_backend, "warplda")) 0.5 else 0.9
  single_cut_mode <- "0.5"
  single_dataset_tag <- if (length(celllines)) as.character(celllines[[1]]) else "dataset"
  single_model_dir <- file.path(step3_root_dir, "topic_models", single_combo_tag)
  single_final_dir <- file.path(
    step3_root_dir,
    build_topic_compact_run_dirname(
      thrP_use = single_thrP_use,
      cut_mode = single_cut_mode,
      gene_term_mode = as.character(single_row$gene_term_mode[[1]]),
      include_tf_terms = isTRUE(single_row$include_tf_terms[[1]]),
      count_input = as.character(single_row$count_input[[1]]),
      dataset_tag = single_dataset_tag,
      doc_mode = topic_doc_mode,
      backend = single_backend,
      vae_variant = single_variant,
      k_use = single_k
    )
  )
  .log_inform(
    "single_topic_test row={single_combo_idx}: combo={single_combo_tag}, backend={single_backend}, variant={single_variant}, k={single_k}"
  )
  single_link_fdr_p_use <- if (identical(single_backend, "warplda")) 0.01 else 0.1
  write_topic_directory_name_readme(
    single_final_dir,
    fields = list(
      thrP = single_thrP_use,
      link_cutoff_mode = single_cut_mode,
      link_method = "gene_prob",
      gene_term_mode = as.character(single_row$gene_term_mode[[1]]),
      include_tf_terms = isTRUE(single_row$include_tf_terms[[1]]),
      count_input = as.character(single_row$count_input[[1]]),
      doc_mode = topic_doc_mode,
      backend = single_backend,
      vae_variant = single_variant,
      selected_k = single_k,
      dataset_tag = single_dataset_tag
    )
  )

  single_topic_args <- make_topic_report_args_simple(
    thrP = single_thrP_use,
    link_prob_cutoff = 0.5,
    link_fdr_p = single_link_fdr_p_use,
    modules = extract_modules,
    overwrite = extract_overwrite
  )

  train_topic_models(
    Kgrid = k_grid_default,
    input_dir = diff_res$filtered_dir,
    output_dir = single_model_dir,
    celllines = celllines,
    tf_cluster_map = motif_info$tf_cluster_map,
    doc_mode = topic_doc_mode,
    tf_exclude = motif_info$tf_exclude,
    gene_term_mode = as.character(single_row$gene_term_mode[[1]]),
    include_tf_terms = isTRUE(single_row$include_tf_terms[[1]]),
    count_input = as.character(single_row$count_input[[1]]),
    backend = single_backend,
    vae_variant = single_variant,
    topic_report_args = single_topic_args
  )

  extract_regulatory_topics(
    k = single_k,
    model_dir = single_model_dir,
    output_dir = single_final_dir,
    backend = single_backend,
    vae_variant = single_variant,
    doc_mode = topic_doc_mode,
    flatten_single_output = TRUE,
    topic_report_args = single_topic_args
  )

  run_vae_doc_topic_heatmaps(
    topic_root = single_final_dir,
    backend = single_backend,
    vae_variant = single_variant,
    doc_mode = topic_doc_mode
  )

  if (isTRUE(run_subnetwork_plots)) {
    run_vae_topic_delta_network_plots(
      topic_root = single_final_dir,
      step2_out_dir = diff_res$filtered_dir,
      min_prob = rerun_subnetwork_min_prob,
      filter_same_direction = rerun_subnetwork_filter_same_direction,
      methods = rerun_subnetwork_methods,
      backend = single_backend,
      vae_variant = single_variant,
      doc_mode = topic_doc_mode
    )
  }

  run_single_topic_pathway_test <- isTRUE(run_pathway_rerun)
  if (isTRUE(run_single_topic_pathway_test)) {
    run_vae_topic_delta_network_pathway(
      topic_root = single_final_dir,
      backend = single_backend,
      vae_variant = single_variant,
      doc_mode = topic_doc_mode,
      top_n_per_topic = single_topic_args$top_n_per_topic,
      max_pathways = single_topic_args$max_pathways,
      per_comparison = isTRUE(run_pathway_per_comparison),
      split_direction = isTRUE(run_pathway_split_direction)
    )
  }

  # Quick debug checks for topic_links/topic_terms
  out_dirs <- list.dirs(single_final_dir, recursive = FALSE, full.names = TRUE)
  doc_tag <- if (identical(topic_doc_mode, "tf")) "tf" else "ctf"
  if (single_backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", single_variant, "_K", single_k, "$")
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  } else {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", single_k, "$")
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  }
  if (!length(out_dirs) && file.exists(file.path(single_final_dir, "topic_links.csv"))) {
    out_dirs <- single_final_dir
  }
  if (length(out_dirs)) {
    d <- out_dirs[[1]]
    tl_path <- file.path(d, "topic_links.csv")
    tt_path <- file.path(d, "topic_terms.csv")
    if (file.exists(tl_path)) {
      tl <- data.table::fread(tl_path)
      .log_inform(sprintf(
        "topic_links: total=%d, peak_pass=%d, gene_pass=%d, both_pass=%d, either_pass=%d",
        nrow(tl),
        sum(tl$peak_pass, na.rm = TRUE),
        sum(tl$gene_pass, na.rm = TRUE),
        sum(tl$peak_pass & tl$gene_pass, na.rm = TRUE),
        sum(tl$peak_pass | tl$gene_pass, na.rm = TRUE)
      ))
    } else {
      .log_warn("topic_links.csv not found in {d}.")
    }
    if (file.exists(tt_path)) {
      tt <- data.table::fread(tt_path)
      .log_inform(sprintf("topic_terms cols: %s", paste(names(tt), collapse = ", ")))
      if ("in_topic" %in% names(tt)) {
        in_topic_flag <- .as_logical_flag(tt$in_topic)
        in_topic_n <- sum(in_topic_flag, na.rm = TRUE)
        .log_inform(sprintf(
          "topic_terms: total=%d, in_topic=%d",
          nrow(tt),
          in_topic_n
        ))
        if (in_topic_n == 0) {
          if ("score" %in% names(tt)) {
            .log_inform(sprintf(
              "topic_terms scores: min=%.4g, max=%.4g",
              min(tt$score, na.rm = TRUE),
              max(tt$score, na.rm = TRUE)
            ))
          }
          if ("topic" %in% names(tt)) {
            topic_counts <- data.table::as.data.table(tt)[, .N, by = topic]
            .log_inform(sprintf(
              "topic_terms per topic: min=%d, max=%d",
              min(topic_counts$N, na.rm = TRUE),
              max(topic_counts$N, na.rm = TRUE)
            ))
          }
          gamma_path <- file.path(d, "topic_gamma_cutoffs.csv")
          if (file.exists(gamma_path)) {
            gc <- data.table::fread(gamma_path)
            if (all(c("peaks_gamma_cutoff", "gene_gamma_cutoff") %in% names(gc))) {
              .log_inform(sprintf(
                "peaks_gamma_cutoff: min=%.4g, max=%.4g",
                min(gc$peaks_gamma_cutoff, na.rm = TRUE),
                max(gc$peaks_gamma_cutoff, na.rm = TRUE)
              ))
              .log_inform(sprintf(
                "gene_gamma_cutoff: min=%.4g, max=%.4g",
                min(gc$gene_gamma_cutoff, na.rm = TRUE),
                max(gc$gene_gamma_cutoff, na.rm = TRUE)
              ))
            }
          }
        }
      } else {
        .log_warn("topic_terms.csv missing in_topic column in {d}.")
      }
    } else {
      .log_warn("topic_terms.csv not found in {d}.")
    }
  } else {
    .log_warn("No final topic directory found for single_topic_test.")
  }
}

run_gammafit_preview <- FALSE
if (isTRUE(do_topic_combo_run) && isTRUE(run_gammafit_preview)) {
  preview_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  preview_backend <- "vae"
  preview_variant <- "multivi_encoder"
  preview_k <- selected_k
  preview_thrP_use <- if (identical(preview_backend, "warplda")) 0.5 else 0.9
  preview_cut_mode <- "0.5"
  preview_link_fdr_p_use <- if (identical(preview_backend, "warplda")) 0.01 else 0.1
  preview_model_dir <- file.path(step3_root_dir, "topic_models", preview_combo_tag)
  preview_final_dir <- file.path(
    step3_root_dir,
    build_topic_compact_run_dirname(
      thrP_use = preview_thrP_use,
      cut_mode = preview_cut_mode,
      gene_term_mode = "aggregate",
      include_tf_terms = FALSE,
      count_input = "pseudo_count_bin",
      dataset_tag = if (length(celllines)) as.character(celllines[[1]]) else "dataset",
      doc_mode = topic_doc_mode,
      backend = preview_backend,
      vae_variant = preview_variant,
      k_use = preview_k
    )
  )

  extract_regulatory_topics(
    k = preview_k,
    model_dir = preview_model_dir,
    output_dir = preview_final_dir,
    backend = preview_backend,
    vae_variant = preview_variant,
    doc_mode = topic_doc_mode,
    flatten_single_output = TRUE,
    topic_report_args = list(
      binarize_method = "gammafit",
      thrP = preview_thrP_use,
      pathway_make_heatmap = FALSE,
      pathway_make_dotplot = FALSE,
      pathway_per_comparison = FALSE,
      run_pathway_gsea = FALSE,
      run_link_topic_scores = TRUE,
      link_topic_overwrite = TRUE,
      link_topic_method = "gene_prob",
      link_topic_prob_cutoff = 0.5,
      link_topic_fdr_q = 0.5,
      link_topic_fdr_p = preview_link_fdr_p_use
    )
  )

  out_dirs <- list.dirs(preview_final_dir, recursive = FALSE, full.names = TRUE)
  doc_tag <- if (identical(topic_doc_mode, "tf")) "tf" else "ctf"
  if (preview_backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", preview_variant, "_K", preview_k, "$")
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", preview_k, "$"), basename(out_dirs))]
  }
  if (!length(out_dirs) && file.exists(file.path(preview_final_dir, "topic_links.csv"))) {
    out_dirs <- preview_final_dir
  }
  if (length(out_dirs)) {
    d <- out_dirs[[1]]
    .log_inform("Gammafit preview outputs:")
    .log_inform("  {file.path(d, 'topic_terms_and_cutoffs_summary.pdf')}")
    .log_inform("  {file.path(d, 'topic_terms.csv')}")
    .log_inform("  {file.path(d, 'topic_links.csv')}")
  } else {
    .log_warn("No final topic directory found for gammafit preview.")
  }
}

run_doc_topic_summary_preview <- FALSE
if (isTRUE(do_topic_combo_run) && isTRUE(run_doc_topic_summary_preview)) {
  preview_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  preview_backend <- "vae"
  preview_variant <- "multivi_encoder"
  preview_k <- 12L
  preview_final_dir <- file.path(
    step3_root_dir,
    build_topic_compact_run_dirname(
      thrP_use = 0.9,
      cut_mode = "0.5",
      gene_term_mode = "aggregate",
      include_tf_terms = FALSE,
      count_input = "pseudo_count_bin",
      dataset_tag = if (length(celllines)) as.character(celllines[[1]]) else "dataset",
      doc_mode = topic_doc_mode,
      backend = preview_backend,
      vae_variant = preview_variant,
      k_use = preview_k
    )
  )

  out_dirs <- list.dirs(preview_final_dir, recursive = FALSE, full.names = TRUE)
  doc_tag <- if (identical(topic_doc_mode, "tf")) "tf" else "ctf"
  if (preview_backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", preview_variant, "_K", preview_k, "$")
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", preview_k, "$"), basename(out_dirs))]
  }
  if (!length(out_dirs) && file.exists(file.path(preview_final_dir, "topic_links.csv"))) {
    out_dirs <- preview_final_dir
  }
  if (length(out_dirs)) {
    link_path <- file.path(out_dirs[1], "topic_links.csv")
    if (file.exists(link_path)) {
      link_dt <- data.table::fread(link_path)
      for (method in c("peak_and_gene", "peak_and_gene_prob")) {
        link_scores <- .topic_links_to_link_scores(link_dt, method = method)
        if (!nrow(link_scores)) next
        plot_topic_delta_networks_from_link_scores(
          link_scores = link_scores,
          step2_out_dir = diff_res$filtered_dir,
          out_root = file.path(out_dirs[1], switch(
            as.character(method),
            peak_and_gene = "subnet_peak_gene",
            peak_and_gene_prob = "subnet_peak_gene_prob",
            paste0("subnet_", .short_link_method_tag(method))
          )),
          min_prob = 0.5,
          filter_same_direction = FALSE
        )
      }
    } else {
      .log_warn("No topic_links.csv found for preview.")
    }
  } else {
    .log_warn("No final topic directory found for preview.")
  }
}

# Benchmark: all peak-gene topic concordance scatter (all methods)
run_topic_benchmark <- FALSE
if (isTRUE(do_topic_combo_run) && isTRUE(run_topic_benchmark)) {
  .log_abort(
    "run_topic_benchmark currently expects legacy final_topics_gamma_* directories. Update benchmark utils to compact folder naming before enabling this block."
  )
  source("dev/benchmark/utils_step3_topic_benchmark.R")
  combo_grid <- expand.grid(
    gene_term_mode = gene_term_modes,
    include_tf_terms = include_tf_terms_opts,
    count_input = count_inputs,
    backend = backends,
    stringsAsFactors = FALSE
  )
  for (cut_mode in gene_prob_cutoff_modes) {
    final_topics_dirname <- paste0("final_topics_gamma_thrP_vae_0.9_peak_and_gene_vs_peak_and_gene_prob_", cut_mode)
    benchmark_dir <- file.path(
      base_dir,
      "diff_grn_and_regulatory_topics",
      paste0(benchmark_dirname, "_", cut_mode)
    )
    dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)
    plot_peak_gene_concordance_all_methods(
      combo_grid = combo_grid,
      step3_root_dir = step3_root_dir,
      final_topics_subdir = final_topics_dirname,
      k = 10L,
      vae_variant = "multivi_encoder",
      out_file = file.path(benchmark_dir, "peak_gene_concordance_all_methods.pdf"),
      point_alpha = 0.05,
      point_size = 0.4
    )
    plot_shared_topic_counts_all_methods(
      combo_grid = combo_grid,
      step3_root_dir = step3_root_dir,
      final_topics_subdir = final_topics_dirname,
      k = 10L,
      vae_variant = "multivi_encoder",
      out_file = file.path(benchmark_dir, "shared_topic_counts_all_methods.pdf")
    )
    plot_pathway_logp_hist_all_methods(
      combo_grid = combo_grid,
      step3_root_dir = step3_root_dir,
      final_topics_subdir = final_topics_dirname,
      k = 10L,
      vae_variant = "multivi_encoder",
      out_file = file.path(benchmark_dir, "pathway_logp_hist_all_methods.pdf"),
      facets_ncol = 6L
    )
    plot_pass_state_counts_all_methods(
      combo_grid = combo_grid,
      step3_root_dir = step3_root_dir,
      final_topics_subdir = final_topics_dirname,
      k = 10L,
      vae_variant = "multivi_encoder",
      out_file = file.path(benchmark_dir, "pass_state_counts_all_methods.pdf"),
      facets_ncol = 6L
    )
  }
}
