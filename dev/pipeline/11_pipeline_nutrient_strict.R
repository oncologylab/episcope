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

# Diff GRN and regulatory topics ------------------------------------------
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
combo_error_log <- file.path(step3_root_dir, "benchmark", "topic_combo_errors.log")
dir.create(dirname(combo_error_log), recursive = TRUE, showWarnings = FALSE)
if (file.exists(combo_error_log)) file.remove(combo_error_log)
combo_failures <- vector("list", 0L)

.append_combo_error_log <- function(line) {
  con <- file(combo_error_log, open = "a")
  on.exit(close(con), add = TRUE)
  writeLines(line, con = con, sep = "\n")
}

.log_inform("Starting topic-model combo run: {nrow(combo_grid)} combinations.")
for (i in seq_len(nrow(combo_grid))) {
  gene_term_mode <- combo_grid$gene_term_mode[i]
  include_tf_terms <- combo_grid$include_tf_terms[i]
  count_input <- combo_grid$count_input[i]
  backend <- combo_grid$backend[i]
  vae_variant <- combo_grid$vae_variant[i]
  thrP_use <- if (identical(backend, "warplda")) 0.5 else 0.8
  link_fdr_p_use <- if (identical(backend, "warplda")) 0.01 else 0.1

  combo_tag <- paste0(
    "gene_", gene_term_mode,
    "_tf_", if (isTRUE(include_tf_terms)) "on" else "off",
    "_count_", count_input
  )
  combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_tag)
  .log_inform(
    "[{i}/{nrow(combo_grid)}] Running combo={combo_tag} backend={backend} (thrP={thrP_use}, link_method=link_score_prob, link_prob_cutoff=0.3)."
  )

  topic_model_dir <- file.path(step3_root_dir, "topic_models", combo_tag)
  topic_final_dir <- file.path(step3_root_dir, "final_topics", combo_tag)

  topic_args <- list(
    pathway_source = "link_scores",
    thrP = thrP_use,
    pathway_make_heatmap = FALSE,
    pathway_make_dotplot = TRUE,
    pathway_overwrite = TRUE,
    pathway_per_comparison = FALSE,
    pathway_per_comparison_dir = "per_cmpr_topic_pathway",
    pathway_split_direction = TRUE,
    run_pathway_gsea = FALSE,
    run_link_topic_scores = TRUE,
    link_topic_gate_mode = "none",
    link_topic_overwrite = TRUE,
    link_topic_method = "link_score_prob",
    link_topic_prob_cutoff = 0.3,
    link_topic_fdr_q = 0.5,
    link_topic_fdr_p = link_fdr_p_use,
    pathway_link_scores_file = "topic_links.csv",
    pathway_link_scores_file_tf = "topic_links.csv",
    pathway_link_gene_terms_file = NULL,
    pathway_link_min_prob = 0,
    pathway_link_include_tf = TRUE,
    pathway_link_include_gene = TRUE,
    pathway_link_gene_min_prob = 0,
    pathway_link_tf_min_prob = 0,
    pathway_link_tf_max_topics = Inf,
    pathway_link_tf_top_n_per_topic = NA_integer_,
    top_n_per_topic = Inf,
    max_pathways = Inf
  )

  combo_ok <- tryCatch(
    {
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
        topic_report_args = topic_args
      )

      run_vae_doc_topic_heatmaps(
        topic_root = topic_final_dir,
        backend = backend,
        vae_variant = vae_variant,
        doc_mode = topic_doc_mode
      )

      run_vae_topic_delta_network_plots(
        topic_root = topic_final_dir,
        step2_out_dir = diff_res$filtered_dir
        ,backend = backend,
        vae_variant = vae_variant,
        doc_mode = topic_doc_mode
      )
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
        "[%s]\trow=%d\tcombo=%s\tbackend=%s\tgene_term_mode=%s\tinclude_tf_terms=%s\tcount_input=%s\terror=%s",
        stamp, i, combo_tag, backend, gene_term_mode, include_tf_terms, count_input, err_msg
      ))
      .log_warn("Combo failed [{i}/{nrow(combo_grid)}]: {combo_tag} ({backend}) -> {err_msg}")
      FALSE
    }
  )
  if (isTRUE(combo_ok)) {
    .log_inform("[{i}/{nrow(combo_grid)}] Completed combo={combo_tag} backend={backend}.")
  }
  # Pathway enrichment already runs inside extract_regulatory_topics(topic_report_args).
  # Keep this disabled to avoid duplicate/repeated pathway runs.
  # run_vae_topic_delta_network_pathway(
  #   topic_root = topic_final_dir,
  #   backend = backend,
  #   vae_variant = vae_variant,
  #   top_n_per_topic = topic_args$top_n_per_topic,
  #   max_pathways = topic_args$max_pathways
  # )
}

if (!length(combo_failures)) {
  .log_inform("All {nrow(combo_grid)} combinations completed successfully.")
} else {
  fail_df <- do.call(rbind, combo_failures)
  .log_warn(
    "Completed with {nrow(fail_df)} failed combination(s). Error log: {combo_error_log}"
  )
  for (j in seq_len(nrow(fail_df))) {
    .log_warn(sprintf(
      "FAILED row=%d | combo=%s | backend=%s | gene_term_mode=%s | include_tf_terms=%s | count_input=%s | error=%s",
      fail_df$row[j],
      fail_df$combo_tag[j],
      fail_df$backend[j],
      fail_df$gene_term_mode[j],
      fail_df$include_tf_terms[j],
      fail_df$count_input[j],
      fail_df$error[j]
    ))
  }
}


run_single_topic_test <- FALSE
if (isTRUE(run_single_topic_test)) {
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
  single_model_dir <- file.path(step3_root_dir, "topic_models", single_combo_tag)
  single_final_dir <- file.path(step3_root_dir, "final_topics", single_combo_tag)
  .log_inform(
    "single_topic_test row={single_combo_idx}: combo={single_combo_tag}, backend={single_backend}, variant={single_variant}, k={single_k}"
  )
  single_thrP_use <- if (identical(single_backend, "warplda")) 0.5 else 0.8
  single_link_fdr_p_use <- if (identical(single_backend, "warplda")) 0.01 else 0.1

  single_topic_args <- list(
    pathway_source = "link_scores",
    thrP = single_thrP_use,
    pathway_make_heatmap = FALSE,
    pathway_make_dotplot = TRUE,
    pathway_overwrite = TRUE,
    pathway_per_comparison = FALSE,
    pathway_per_comparison_dir = "per_cmpr_topic_pathway",
    pathway_split_direction = TRUE,
    run_pathway_gsea = FALSE,
    run_link_topic_scores = TRUE,
    link_topic_gate_mode = "none",
    link_topic_overwrite = TRUE,
    link_topic_method = "link_score_prob",
    link_topic_prob_cutoff = 0.3,
    link_topic_fdr_q = 0.5,
    link_topic_fdr_p = single_link_fdr_p_use,
    pathway_link_scores_file = "topic_links.csv",
    pathway_link_scores_file_tf = "topic_links.csv",
    pathway_link_gene_terms_file = NULL,
    pathway_link_min_prob = 0,
    pathway_link_include_tf = TRUE,
    pathway_link_include_gene = TRUE,
    pathway_link_gene_min_prob = 0,
    pathway_link_tf_min_prob = 0,
    pathway_link_tf_max_topics = Inf,
    pathway_link_tf_top_n_per_topic = NA_integer_,
    top_n_per_topic = Inf,
    max_pathways = Inf
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
    topic_report_args = single_topic_args
  )

  run_vae_doc_topic_heatmaps(
    topic_root = single_final_dir,
    backend = single_backend,
    vae_variant = single_variant,
    doc_mode = topic_doc_mode
  )

  run_vae_topic_delta_network_plots(
    topic_root = single_final_dir,
    step2_out_dir = diff_res$filtered_dir,
    backend = single_backend,
    vae_variant = single_variant,
    doc_mode = topic_doc_mode
  )

  run_single_topic_pathway_test <- FALSE
  if (isTRUE(run_single_topic_pathway_test)) {
    run_vae_topic_delta_network_pathway(
      topic_root = single_final_dir,
      backend = single_backend,
      vae_variant = single_variant,
      doc_mode = topic_doc_mode,
      top_n_per_topic = single_topic_args$top_n_per_topic,
      max_pathways = single_topic_args$max_pathways,
      per_comparison = FALSE
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
if (isTRUE(run_gammafit_preview)) {
  preview_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  preview_backend <- "vae"
  preview_variant <- "multivi_encoder"
  preview_k <- selected_k
  preview_thrP_use <- if (identical(preview_backend, "warplda")) 0.5 else 0.8
  preview_link_fdr_p_use <- if (identical(preview_backend, "warplda")) 0.01 else 0.1
  preview_model_dir <- file.path(step3_root_dir, "topic_models", preview_combo_tag)
  preview_final_dir <- file.path(step3_root_dir, "final_topics", preview_combo_tag)

  extract_regulatory_topics(
    k = preview_k,
    model_dir = preview_model_dir,
    output_dir = preview_final_dir,
    backend = preview_backend,
    vae_variant = preview_variant,
    doc_mode = topic_doc_mode,
    topic_report_args = list(
      binarize_method = "gammafit",
      thrP = preview_thrP_use,
      pathway_make_heatmap = FALSE,
      pathway_make_dotplot = FALSE,
      pathway_per_comparison = FALSE,
      run_pathway_gsea = FALSE,
      run_link_topic_scores = TRUE,
      link_topic_overwrite = TRUE,
      link_topic_method = "link_score_prob",
      link_topic_prob_cutoff = 0.3,
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
if (isTRUE(run_doc_topic_summary_preview)) {
  preview_combo_tag <- "gene_aggregate_tf_off_count_pseudo_count_bin"
  preview_backend <- "vae"
  preview_variant <- "multivi_encoder"
  preview_k <- 12L
  preview_final_dir <- file.path(step3_root_dir, "final_topics", preview_combo_tag)

  out_dirs <- list.dirs(preview_final_dir, recursive = FALSE, full.names = TRUE)
  doc_tag <- if (identical(topic_doc_mode, "tf")) "tf" else "ctf"
  if (preview_backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", preview_variant, "_K", preview_k, "$")
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", preview_k, "$"), basename(out_dirs))]
  }
  if (length(out_dirs)) {
    link_path <- file.path(out_dirs[1], "topic_links.csv")
    if (file.exists(link_path)) {
      link_dt <- data.table::fread(link_path)
      link_scores <- .topic_links_to_link_scores(link_dt, method = "peak_and_gene")
      if (nrow(link_scores)) {
        plot_topic_delta_networks_from_link_scores(
          link_scores = link_scores,
          step2_out_dir = diff_res$filtered_dir,
          out_root = file.path(out_dirs[1], "doc_topic_sub_network_link_peak_and_gene"),
          min_prob = 0.5,
          filter_same_direction = TRUE
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
run_topic_benchmark <- TRUE
if (isTRUE(run_topic_benchmark)) {
  source("dev/benchmark/utils_step3_topic_benchmark.R")
  benchmark_dir <- file.path(base_dir, "diff_grn_and_regulatory_topics", "benchmark")
  dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)
  combo_grid <- expand.grid(
    gene_term_mode = gene_term_modes,
    include_tf_terms = include_tf_terms_opts,
    count_input = count_inputs,
    backend = backends,
    stringsAsFactors = FALSE
  )
  plot_peak_gene_concordance_all_methods(
    combo_grid = combo_grid,
    step3_root_dir = step3_root_dir,
    k = 10L,
    vae_variant = "multivi_encoder",
    out_file = file.path(benchmark_dir, "peak_gene_concordance_all_methods.pdf"),
    point_alpha = 0.05,
    point_size = 0.4
  )
  plot_shared_topic_counts_all_methods(
    combo_grid = combo_grid,
    step3_root_dir = step3_root_dir,
    k = 10L,
    vae_variant = "multivi_encoder",
    out_file = file.path(benchmark_dir, "shared_topic_counts_all_methods.pdf")
  )
  plot_pathway_logp_hist_all_methods(
    combo_grid = combo_grid,
    step3_root_dir = step3_root_dir,
    k = 10L,
    vae_variant = "multivi_encoder",
    out_file = file.path(benchmark_dir, "pathway_logp_hist_all_methods.pdf"),
    facets_ncol = 6L
  )
  plot_pass_state_counts_all_methods(
    combo_grid = combo_grid,
    step3_root_dir = step3_root_dir,
    k = 10L,
    vae_variant = "multivi_encoder",
    out_file = file.path(benchmark_dir, "pass_state_counts_all_methods.pdf"),
    facets_ncol = 6L
  )
}
