library(enrichR)
library(episcope)

base_dir <- "/data/homes/yl814/episcope_test/nutrient_stress"
utils_path <- file.path("R", "utils_step3_topic_warplda.R")
if (!file.exists(utils_path)) cli::cli_abort("Missing utils file: {utils_path}")
utils_path <- normalizePath(utils_path, winslash = "/", mustWork = TRUE)
source(utils_path)

step2_out_dir <- file.path(base_dir, "connect_tf_target_genes")
delta_links_csv <- list.files(step2_out_dir, "_delta_links\\.csv$", full.names = TRUE)
if (!length(delta_links_csv)) cli::cli_abort("No *_delta_links.csv files found in {step2_out_dir}")

edges_all <- load_delta_links_many(delta_links_csv, keep_original = TRUE)
topic_root <- file.path(base_dir, "benchmark_topic_vae_model_5_rna_fp_combined") # model_selection/specified_k
K_grid_default <- c(2:15, 20, 25, 35, 40, 45, 50, 60, 70, 80, 90, 100)

motif_path <- resolve_motif_db_path("JASPAR2024", ref_genome = ref_genome)
motif_path <- "inst/extdata/genome/JASPAR2024.txt"  # fallback only
# Prefer JASPAR auto-CTF db with sub_cluster info if available
motif_db_jaspar_auto <- file.path(base_dir, "tf_motif_clustering", "motif_db_JASPAR2024_data_K220.csv")
if (file.exists(motif_db_jaspar_auto)) {
  motif_info <- build_tf_cluster_map_from_motif(motif_db_jaspar_auto)
} else {
  motif_info <- build_tf_cluster_map_from_motif(motif_path)
}
tf_cluster_map <- motif_info$tf_cluster_map
tf_exclude <- motif_info$tf_exclude

common_args <- list(
  abs_log2fc_fp_min = 0,
  # abs_log2fc_fp_min = 0.5,

  abs_delta_fp_min = 1, # set to NA_real_ for not filtering by abs_delta_fp_min
  # abs_delta_fp_min = NA_real_,

  abs_log2fc_gene_min = 1,
  require_fp_bound_either = TRUE,
  require_tf_expr_either = TRUE,
  require_gene_expr_either = TRUE,
  direction_consistency = "aligned",
  top_terms_per_doc = Inf,
  min_df = 2,
  count_method = "bin",
  count_scale = 50,
  K_grid = K_grid_default,
  iterations = 2000,
  beta = 0.1,
  seed = 123,
  binarize_method = "gammafit",
  thrP = 0.9, # gamma‑fit cutoff percentile for in_set
  top_n_terms = 500, # only used in binarize_topics() when binarize_method = "topn" to define in_set(TRUE/FALSE)
  in_set_min_terms = 1,
  pathway_use_all_terms = FALSE,
  pathway_make_heatmap = FALSE,
  top_n_per_topic = 100L,
  max_pathways = 1000L,
  pathway_tf_link_mode = "theta",
  pathway_tf_top_n_docs = 50L,
  pathway_tf_min_theta = NA_real_,
  run_pathway_gsea = FALSE,
  gsea_species = "Homo sapiens",
  gsea_nperm = 1000,
  gsea_peak_gene_agg = "max",
  pathway_source = "topic_terms",
  pathway_link_scores_file = NULL,
  pathway_link_min_prob = 0,
  pathway_link_include_tf = TRUE,
  pathway_link_include_gene = TRUE,
  pathway_per_comparison = FALSE,
  pathway_per_comparison_dir = "per_cmpr_topic_pathway",
  pathway_split_direction = TRUE,
  run_link_topic_scores = FALSE,
  link_topic_gate_mode = "none",
  link_topic_top_k = 3L,
  link_topic_min_prob = 0,
  link_topic_include_tf = FALSE,
  link_topic_chunk_size = 5000L,
  link_topic_n_cores = 1L,
  link_topic_overwrite = FALSE
)

# VAE pathway rerun defaults (used by run_vae_pathway_rerun too)
vae_pathway_source <- "link_scores"
vae_pathway_link_scores_file <- NULL
vae_pathway_link_scores_file_tf <- "link_topic_scores_gate_peak_and_gene_in_set.csv"
vae_pathway_link_gene_terms_file <- "topic_terms.csv"
vae_pathway_link_min_prob <- 0
vae_pathway_link_include_tf <- TRUE
vae_pathway_link_include_gene <- TRUE
vae_pathway_link_gene_min_prob <- 0
vae_pathway_link_tf_min_prob <- 0.5
vae_pathway_link_tf_max_topics <- 5L
vae_pathway_link_tf_top_n_per_topic <- 30L
vae_pathway_per_comparison <- TRUE
vae_pathway_per_comparison_dir <- "per_cmpr_topic_pathway"
vae_pathway_split_direction <- TRUE

celllines <- c("AsPC1", "HPAFII", "Panc1")
get_cellline <- function(x) sub("_.*$", "", x)
edges_all[, cellline := get_cellline(comparison_id)]

opt_defs <- list(
  list(label = "opt1_peak_delta_fp", suffix = "peak_delta_fp", direction_by = "fp"),
  list(label = "opt2_peak_fc_fp", suffix = "peak_fc_fp", direction_by = "fp"),
  list(label = "opt3_gene_fc_expr", suffix = "gene_fc_expr", direction_by = "gene")
)

doc_modes <- list(
  list(name = "tf_docs", doc_mode = "tf", direction_by = "gene", distinct_terms = FALSE, tf_cluster_map = NULL),
  list(name = "cmpr_docs", doc_mode = "comparison", direction_by = NULL, distinct_terms = TRUE, tf_cluster_map = NULL),
  list(name = "ctf_docs", doc_mode = "tf_cluster", direction_by = "gene", distinct_terms = TRUE, tf_cluster_map = tf_cluster_map, tf_exclude = tf_exclude)
)

.format_k_label <- function(k_grid) {
  k_grid <- unique(as.integer(k_grid))
  k_grid <- k_grid[is.finite(k_grid)]
  if (!length(k_grid)) return("Knone")
  if (length(k_grid) == 1L) return(paste0("K", k_grid))
  paste0("K", paste(k_grid, collapse = "-"))
}

run_one <- function(edges_in, label_prefix, doc_mode_cfg, opt_cfg, k_grid = NULL) {
  if (!is.null(doc_mode_cfg$tf_exclude) && length(doc_mode_cfg$tf_exclude)) {
    dt <- data.table::as.data.table(edges_in)
    dt <- dt[!toupper(tf) %in% doc_mode_cfg$tf_exclude]
    edges_in <- dt
  }
  k_label <- if (!is.null(k_grid) && length(k_grid)) paste0("_", .format_k_label(k_grid)) else ""
  dir_name <- sprintf("%s_%s_%s%s", label_prefix, doc_mode_cfg$name, opt_cfg$suffix, k_label)
  dir_path <- file.path(topic_root, dir_name)
  direction_by <- if (!is.null(doc_mode_cfg$direction_by)) doc_mode_cfg$direction_by else opt_cfg$direction_by
  local_common_args <- common_args
  if (!is.null(k_grid) && length(k_grid)) {
    local_common_args$K_grid <- as.integer(k_grid)
  }
  topic_label_cleaner <- function(label) {
    parts <- data.table::tstrsplit(label, "::", fixed = TRUE)
    base <- parts[[1]]
    dir <- if (length(parts) >= 2L) parts[[2]] else NA_character_
    base <- sub("^[^_]+_", "", base)
    base <- sub("_vs_[^_]+_10_FBS", "", base)
    out <- base
    has_dir <- !is.na(dir) & nzchar(dir)
    out[has_dir] <- paste(base[has_dir], dir[has_dir], sep = "::")
    out
  }
  args <- c(list(edges_all = edges_in), local_common_args, list(
    doc_mode = doc_mode_cfg$doc_mode,
    direction_by = direction_by,
    distinct_terms = doc_mode_cfg$distinct_terms,
    out_dir = dir_path,
    option_label = opt_cfg$label,
    topic_by_comparison_label_cleaner = topic_label_cleaner
  ))
  if (!is.null(doc_mode_cfg$tf_cluster_map)) {
    args$tf_cluster_map <- doc_mode_cfg$tf_cluster_map
  }
  do.call(run_tfdocs_warplda_one_option, args)
}

error_log <- file.path(topic_root, "run_errors.log")
run_job <- function(job) {
  label <- sprintf("%s_%s_%s", job$label_prefix, job$doc_mode_cfg$name, job$opt_cfg$suffix)
  if (!is.null(job$k_grid) && length(job$k_grid)) {
    label <- sprintf("%s_%s", label, .format_k_label(job$k_grid))
  }
  tryCatch(
    {
      do.call(run_one, job)
      NULL
    },
    error = function(e) {
      msg <- sprintf("[%s] %s: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), label, conditionMessage(e))
      cat(msg, file = error_log, append = TRUE)
      msg
    }
  )
}

run_full <- FALSE
if (run_full) {
  jobs <- list()
  job_i <- 0L
  for (doc_mode_cfg in doc_modes) {
    for (opt_cfg in opt_defs) {
      job_i <- job_i + 1L
      jobs[[job_i]] <- list(edges_in = edges_all, label_prefix = "All", doc_mode_cfg = doc_mode_cfg, opt_cfg = opt_cfg)
    }
    for (cell in celllines) {
      edges_sub <- edges_all[cellline == cell]
      if (!nrow(edges_sub)) next
      for (opt_cfg in opt_defs) {
        job_i <- job_i + 1L
        jobs[[job_i]] <- list(edges_in = edges_sub, label_prefix = cell, doc_mode_cfg = doc_mode_cfg, opt_cfg = opt_cfg)
      }
    }
  }

  mc_cores <- min(10L, length(jobs))
  if (.Platform$OS.type == "windows" || mc_cores < 2L) {
    lapply(jobs, run_job)
  } else {
    parallel::mclapply(jobs, run_job, mc.cores = mc_cores)
  }
}

run_custom <- FALSE
if (run_custom) {
  custom_doc_mode <- list(
    name = "ctf_docs",
    doc_mode = "tf_cluster",
    direction_by = "gene",
    distinct_terms = TRUE,
    tf_cluster_map = tf_cluster_map,
    tf_exclude = tf_exclude
  )
  custom_opts <- list(
    list(label = "opt1_peak_delta_fp", suffix = "peak_delta_fp", direction_by = "fp"),
    list(label = "opt3_gene_fc_expr", suffix = "gene_fc_expr", direction_by = "gene")
  )
  custom_default_k_grid <- c(20L)

  # from benchmark_topic_warp_lda_2_links_abs_log2fc_fp_min_0.5_abs_log2fc_gene_min_1_model_selection
  # custom_k_grid <- list(
  #   AsPC1 = list(
  #     peak_delta_fp = c(15L, 50L, 80L),
  #     gene_fc_expr = c(15L, 50L, 80L)
  #   ),
  #   HPAFII = list(
  #     peak_delta_fp = c(15L, 50L, 70L),
  #     gene_fc_expr = c(15L, 20L, 50L)
  #   ),
  #   Panc1 = list(
  #     peak_delta_fp = c(13L, 20L, 80L),
  #     gene_fc_expr = c(13L, 20L, 80L)
  #   )
  # )

  # from benchmark_topic_warp_lda_3_links_abs_delta_fp_min_1_abs_log2fc_gene_min_1_model_selection
  custom_k_grid <- list(
    AsPC1 = list(
      peak_delta_fp = c(15L, 50L, 80L),
      gene_fc_expr = c(15L, 50L, 80L)
    ),
    HPAFII = list(
      peak_delta_fp = c(15L, 50L, 90L),
      gene_fc_expr = c(15L, 20L, 50L)
    ),
    Panc1 = list(
      peak_delta_fp = c(13L, 20L, 80L),
      gene_fc_expr = c(13L, 20L, 80L)
    )
  )

  jobs <- list()
  job_i <- 0L
  for (cell in celllines) {
    edges_sub <- edges_all[cellline == cell]
    if (!nrow(edges_sub)) next
    for (opt_cfg in custom_opts) {
      k_override <- custom_default_k_grid
      if (cell %in% names(custom_k_grid)) {
        cell_map <- custom_k_grid[[cell]]
        if (opt_cfg$suffix %in% names(cell_map)) {
          k_override <- cell_map[[opt_cfg$suffix]]
        }
      }
      k_override <- unique(as.integer(k_override))
      k_override <- k_override[is.finite(k_override)]
      for (k_val in k_override) {
        job_i <- job_i + 1L
        jobs[[job_i]] <- list(
          edges_in = edges_sub,
          label_prefix = cell,
          doc_mode_cfg = custom_doc_mode,
          opt_cfg = opt_cfg,
          k_grid = c(k_val)
        )
      }
    }
  }

  mc_cores <- min(20L, length(jobs))
  if (.Platform$OS.type == "windows" || mc_cores < 2L) {
    lapply(jobs, run_job)
  } else {
    parallel::mclapply(jobs, run_job, mc.cores = mc_cores)
  }
}

run_mira_vae <- FALSE
run_mira_vae_extra_only <- TRUE
run_mira_vae_auto_ctf_hpa_only <- TRUE
if (isTRUE(run_mira_vae_extra_only)) {
  run_mira_vae <- TRUE
}
if (run_mira_vae) {
  vae_script <- file.path("dev", "logistic_normal_vae_topics.py")
  if (!file.exists(vae_script)) cli::cli_abort("Missing VAE script: {vae_script}")
  vae_python <- Sys.getenv("VAE_PYTHON", unset = "")
  if (!nzchar(vae_python)) {
    default_py <- "/data/homes/yl814/miniconda3/bin/python"
    if (file.exists(default_py)) vae_python <- default_py
  }
  if (!nzchar(vae_python)) vae_python <- Sys.which("python")
  if (!nzchar(vae_python)) vae_python <- Sys.which("python3")
  if (!nzchar(vae_python)) vae_python <- Sys.which("py")
  if (!nzchar(vae_python)) {
    cli::cli_abort("Python not found on PATH. Set VAE_PYTHON or vae_python to the full path of python.")
  }

  vae_k_grid_default <- K_grid_default
  vae_k_single <- c(20L, 50L, 80L)
  vae_epochs <- 200L
  vae_batch_size <- 64L
  vae_hidden <- 128L
  vae_lr <- 1e-3
  vae_seed <- 123L
  vae_device <- "cpu"
  vae_variants <- c("vae_mlp", "moetm_encoder_decoder", "multivi_encoder")
  vae_reuse_if_exists <- TRUE
  # Extra K runs for specific cellline/docmode/variant.
  # Set run_mira_vae_extra_only <- TRUE above to skip full grids and only run these.
  vae_run_extra_k_only <- isTRUE(run_mira_vae_extra_only)
  vae_extra_k_map <- list(
    HPAFII = list(
      ctf = list(
        moetm_encoder_decoder = c(15L, 40L),
        multivi_encoder = c(12L),
        vae_mlp = c(13L)
      )
    ),
    Panc1 = list(
      ctf = list(
        multivi_encoder = c(11L, 13L)
      )
    ),
    AsPC1 = list(
      ctf = list(
        multivi_encoder = c(10L)
      )
    )
  )

  auto_ctf_data_path <- file.path(base_dir, "tf_motif_clustering", "motif_db_JASPAR2024_data_K220.csv")
  auto_ctf_hybrid_path <- file.path(base_dir, "tf_motif_clustering", "motif_db_JASPAR2024_hybrid_K165.csv")
  auto_ctf_data_map <- NULL
  auto_ctf_hybrid_map <- NULL
  if (file.exists(auto_ctf_data_path)) {
    auto_ctf_data_map <- build_tf_cluster_map_from_motif(auto_ctf_data_path)$tf_cluster_map
  } else {
    cli::cli_inform("auto_ctf_data file not found: {auto_ctf_data_path}")
  }
  if (file.exists(auto_ctf_hybrid_path)) {
    auto_ctf_hybrid_map <- build_tf_cluster_map_from_motif(auto_ctf_hybrid_path)$tf_cluster_map
  } else {
    cli::cli_inform("auto_ctf_hybrid file not found: {auto_ctf_hybrid_path}")
  }

  if (isTRUE(run_mira_vae_auto_ctf_hpa_only)) {
    vae_variants <- c("multivi_encoder")
    vae_run_extra_k_only <- TRUE
    vae_extra_k_map <- list(
      HPAFII = list(
        auto_ctf_data = list(multivi_encoder = c(12L)),
        auto_ctf_hybrid = list(multivi_encoder = c(12L))
      )
    )
  }

  # Link -> topic membership (baseline + optional hard confirmation gate)
  vae_link_topic_scores <- TRUE
  vae_link_topic_gate_mode <- c("none", "peak_and_gene_in_set")
  vae_link_topic_top_k <- 3L
  vae_link_topic_min_prob <- 0
  vae_link_topic_include_tf <- FALSE
  vae_link_topic_chunk_size <- 5000L
  vae_link_topic_n_cores <- 1L
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.finite(cores)) {
    vae_link_topic_n_cores <- max(1L, min(20L, cores))
  }
  vae_link_topic_overwrite <- FALSE

  vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
  dir.create(vae_root, recursive = TRUE, showWarnings = FALSE)

  run_vae_topic_report <- function(doc_term, edges_docs, out_dir, option_label, direction_by, vae_variant, doc_mode_tag, k_grid = NULL) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!nrow(doc_term)) cli::cli_abort("doc_term is empty for VAE; relax filters.")

    data.table::fwrite(doc_term, file.path(out_dir, "doc_term.csv"))
    .save_all(out_dir, "doc_term", doc_term)

    dtm_obj <- build_sparse_dtm(doc_term)
    dtm <- dtm_obj$dtm
    .save_all(out_dir, "dtm", dtm)
    .save_all(out_dir, "dtm_index", list(doc_index = dtm_obj$doc_index, term_index = dtm_obj$term_index))

    if (is.null(k_grid) || !length(k_grid)) k_grid <- vae_k_grid_default
    k_grid_txt <- paste(as.integer(k_grid), collapse = ",")
    metrics_path <- file.path(out_dir, "model_metrics.csv")
    reuse_ok <- isTRUE(vae_reuse_if_exists) &&
      file.exists(metrics_path) &&
      dir.exists(file.path(out_dir, "vae_models"))

    if (!reuse_ok) {
      py_args <- c(
        vae_script,
        "--doc-term", file.path(out_dir, "doc_term.csv"),
        "--out-dir", out_dir,
        "--k-grid", k_grid_txt,
        "--epochs", as.character(vae_epochs),
        "--batch-size", as.character(vae_batch_size),
        "--hidden", as.character(vae_hidden),
        "--lr", as.character(vae_lr),
        "--seed", as.character(vae_seed),
        "--device", vae_device,
        "--variant", vae_variant
      )
      py_out <- tryCatch(
        system2(vae_python, py_args, stdout = TRUE, stderr = TRUE),
        error = function(e) {
          cli::cli_abort(c("Failed to run python for VAE.", i = conditionMessage(e)))
        }
      )
      writeLines(py_out, file.path(out_dir, "vae_train.log"))
      status <- attr(py_out, "status")
      if (!is.null(status) && is.finite(status) && status != 0) {
        cli::cli_abort(c("VAE python exited with non-zero status.", i = paste0("status=", status), i = paste(py_out, collapse = "\n")))
      }
    } else {
      cli::cli_inform("Reusing existing VAE outputs; skipping training for {out_dir}")
    }

    if (!file.exists(metrics_path)) cli::cli_abort("VAE did not produce model_metrics.csv in {out_dir}")
    metrics_tbl <- readr::read_csv(metrics_path, show_col_types = FALSE)
    .save_all(out_dir, "model_metrics", metrics_tbl)

    title_prefix <- paste0("Deep VAE ", vae_variant, " ", basename(out_dir))
    sel <- plot_model_selection_cistopic(metrics_tbl, file.path(out_dir, "model_selection.pdf"), title_prefix = title_prefix)
    .save_all(out_dir, "model_selection", sel)

    m <- sel$table
    idx <- sel$idx
    pick_K <- function(ix, label) {
      if (!is.finite(ix)) return(NULL)
      m$K[ix]
    }
    chosen_K <- pick_K(idx$perplexity, "perplexity")
    if (is.null(chosen_K)) chosen_K <- pick_K(idx$maximum, "maximum")
    if (is.null(chosen_K)) chosen_K <- m$K[1]

    theta_path <- file.path(out_dir, "vae_models", sprintf("theta_K%d.csv", chosen_K))
    phi_path <- file.path(out_dir, "vae_models", sprintf("phi_K%d.csv", chosen_K))
    if (!file.exists(theta_path) || !file.exists(phi_path)) {
      cli::cli_abort("Missing theta/phi outputs for K={chosen_K} in {out_dir}")
    }
    theta_df <- readr::read_csv(theta_path, show_col_types = FALSE)
    phi_df <- readr::read_csv(phi_path, show_col_types = FALSE)

    theta <- as.matrix(theta_df[, -1, drop = FALSE])
    rownames(theta) <- theta_df[[1]]
    phi <- as.matrix(phi_df[, -1, drop = FALSE])
    rownames(phi) <- phi_df[[1]]

    topic_base <- list(K = chosen_K, theta = theta, phi = phi, metrics = metrics_tbl[metrics_tbl$K == chosen_K, ])

    run_tfdocs_report_from_topic_base(
      topic_base = topic_base,
      dtm = dtm,
      edges_docs = edges_docs,
      out_dir = out_dir,
      option_label = option_label,
      direction_by = direction_by,
      binarize_method = common_args$binarize_method,
      thrP = common_args$thrP,
      top_n_terms = common_args$top_n_terms,
      in_set_min_terms = common_args$in_set_min_terms,
      pathway_use_all_terms = common_args$pathway_use_all_terms,
      pathway_make_heatmap = common_args$pathway_make_heatmap,
      top_n_per_topic = common_args$top_n_per_topic,
      max_pathways = common_args$max_pathways,
      pathway_tf_link_mode = common_args$pathway_tf_link_mode,
      pathway_tf_top_n_docs = common_args$pathway_tf_top_n_docs,
      pathway_tf_min_theta = common_args$pathway_tf_min_theta,
      run_pathway_gsea = common_args$run_pathway_gsea,
      gsea_species = common_args$gsea_species,
      gsea_nperm = common_args$gsea_nperm,
      gsea_peak_gene_agg = common_args$gsea_peak_gene_agg,
      pathway_source = vae_pathway_source,
      pathway_link_scores_file = vae_pathway_link_scores_file,
      pathway_link_scores_file_tf = if (!is.null(vae_pathway_link_scores_file_tf)) {
        file.path(out_dir, vae_pathway_link_scores_file_tf)
      } else {
        NULL
      },
      pathway_link_gene_terms_file = if (!is.null(vae_pathway_link_gene_terms_file)) {
        file.path(out_dir, vae_pathway_link_gene_terms_file)
      } else {
        NULL
      },
      pathway_link_min_prob = vae_pathway_link_min_prob,
      pathway_link_include_tf = vae_pathway_link_include_tf,
      pathway_link_include_gene = vae_pathway_link_include_gene,
      pathway_link_gene_min_prob = vae_pathway_link_gene_min_prob,
      pathway_link_tf_min_prob = vae_pathway_link_tf_min_prob,
      pathway_link_tf_max_topics = vae_pathway_link_tf_max_topics,
      pathway_link_tf_top_n_per_topic = vae_pathway_link_tf_top_n_per_topic,
      pathway_per_comparison = vae_pathway_per_comparison,
      pathway_per_comparison_dir = vae_pathway_per_comparison_dir,
      pathway_split_direction = vae_pathway_split_direction,
      run_link_topic_scores = vae_link_topic_scores,
      link_topic_gate_mode = vae_link_topic_gate_mode,
      link_topic_top_k = vae_link_topic_top_k,
      link_topic_min_prob = vae_link_topic_min_prob,
      link_topic_include_tf = vae_link_topic_include_tf,
      link_topic_chunk_size = vae_link_topic_chunk_size,
      link_topic_n_cores = vae_link_topic_n_cores,
      link_topic_overwrite = vae_link_topic_overwrite,
      title_prefix = title_prefix
    )

    invisible(TRUE)
  }

  run_vae_cell <- function(cell) {
    edges_sub <- edges_all[cellline == cell]
    if (!nrow(edges_sub)) next

    edges_filt <- filter_edges_for_tf_topics(
      edges_sub,
      abs_log2fc_fp_min = common_args$abs_log2fc_fp_min,
      abs_delta_fp_min = common_args$abs_delta_fp_min,
      abs_log2fc_gene_min = common_args$abs_log2fc_gene_min,
      require_fp_bound_either = common_args$require_fp_bound_either,
      require_tf_expr_either = common_args$require_tf_expr_either,
      require_gene_expr_either = common_args$require_gene_expr_either,
      direction_consistency = common_args$direction_consistency
    )
    if (!nrow(edges_filt)) next

    if (!is.null(tf_exclude) && length(tf_exclude)) {
      dt <- data.table::as.data.table(edges_filt)
      dt <- dt[!toupper(tf) %in% tf_exclude]
      edges_filt <- dt
    }
    if (!nrow(edges_filt)) next

    doc_modes <- list(
      list(tag = "tf", doc_mode = "tf", tf_cluster_map = NULL),
      list(tag = "ctf", doc_mode = "tf_cluster", tf_cluster_map = tf_cluster_map)
    )
    if (!is.null(auto_ctf_data_map)) {
      doc_modes <- append(doc_modes, list(
        list(tag = "auto_ctf_data", doc_mode = "tf_cluster", tf_cluster_map = auto_ctf_data_map)
      ))
    }
    if (!is.null(auto_ctf_hybrid_map)) {
      doc_modes <- append(doc_modes, list(
        list(tag = "auto_ctf_hybrid", doc_mode = "tf_cluster", tf_cluster_map = auto_ctf_hybrid_map)
      ))
    }
    if (isTRUE(run_mira_vae_auto_ctf_hpa_only)) {
      doc_modes <- Filter(function(x) x$tag %in% c("auto_ctf_data", "auto_ctf_hybrid"), doc_modes)
    }

    for (mode in doc_modes) {
      edges_docs_joint <- add_tf_docs(
        edges_filt,
        doc_mode = mode$doc_mode,
        direction_by = "gene",
        tf_cluster_map = mode$tf_cluster_map
      )
      doc_term_joint <- build_doc_term_joint(
        edges_docs_joint,
        weight_type_peak = "delta_fp",
        weight_type_gene = "fc_mag_gene",
        top_terms_per_doc = common_args$top_terms_per_doc,
        min_df = common_args$min_df,
        count_method = common_args$count_method,
        count_scale = common_args$count_scale,
        distinct_terms = TRUE,
        balance_mode = "min",
        prefix_terms = TRUE
      )

    if (nrow(doc_term_joint)) {
      for (vae_variant in vae_variants) {
        extra_k <- NULL
        if (cell %in% names(vae_extra_k_map) &&
            mode$tag %in% names(vae_extra_k_map[[cell]]) &&
            vae_variant %in% names(vae_extra_k_map[[cell]][[mode$tag]])) {
          extra_k <- vae_extra_k_map[[cell]][[mode$tag]][[vae_variant]]
        }
        extra_k <- unique(as.integer(extra_k))
        extra_k <- extra_k[is.finite(extra_k)]

        if (!isTRUE(vae_run_extra_k_only)) {
          out_joint_grid <- file.path(
            vae_root,
            paste0(cell, "_vae_joint_", mode$tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_Kgrid_default")
          )
          run_vae_topic_report(
            doc_term_joint,
            edges_docs_joint,
            out_joint_grid,
            "joint",
            "gene",
            vae_variant,
            mode$tag,
            k_grid = vae_k_grid_default
          )
          for (k_val in vae_k_single) {
            out_joint <- file.path(
              vae_root,
              paste0(cell, "_vae_joint_", mode$tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", k_val)
            )
            run_vae_topic_report(
              doc_term_joint,
              edges_docs_joint,
              out_joint,
              "joint",
              "gene",
              vae_variant,
              mode$tag,
              k_grid = c(k_val)
            )
          }
        }

        if (length(extra_k)) {
          for (k_val in extra_k) {
            out_joint <- file.path(
              vae_root,
              paste0(cell, "_vae_joint_", mode$tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", k_val)
            )
            run_vae_topic_report(
              doc_term_joint,
              edges_docs_joint,
              out_joint,
              "joint",
              "gene",
              vae_variant,
              mode$tag,
              k_grid = c(k_val)
            )
          }
        } else if (isTRUE(vae_run_extra_k_only)) {
          cli::cli_inform("Skipping {cell} {mode$tag} {vae_variant}: no extra K configured.")
        }
      }
    } else {
      cli::cli_inform("Skipping VAE joint {mode$tag}: no doc_term for {cell}")
    }
  }
  }

  vae_cells <- celllines
  if (isTRUE(run_mira_vae_auto_ctf_hpa_only)) {
    vae_cells <- "HPAFII"
  }
  mc_cores <- min(20L, length(vae_cells))
  if (.Platform$OS.type == "windows" || mc_cores < 2L) {
    lapply(vae_cells, run_vae_cell)
  } else {
    parallel::mclapply(vae_cells, run_vae_cell, mc.cores = mc_cores)
  }
}

# Extra-only postprocess (no retraining): heatmaps + subnetworks + subnet pathways
run_mira_vae_extra_postprocess <- TRUE
if (isTRUE(run_mira_vae_extra_postprocess)) {
  vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
  if (!exists("vae_extra_k_map")) {
    vae_extra_k_map <- list(
      HPAFII = list(
        ctf = list(
          moetm_encoder_decoder = c(15L, 40L),
          multivi_encoder = c(12L),
          vae_mlp = c(13L)
        )
      ),
      Panc1 = list(
        ctf = list(
          multivi_encoder = c(11L, 13L)
        )
      ),
      AsPC1 = list(
        ctf = list(
          multivi_encoder = c(10L)
        )
      )
    )
  }
  .list_extra_vae_dirs <- function(root_dir, extra_map) {
    out <- character(0)
    for (cell in names(extra_map)) {
      cell_map <- extra_map[[cell]]
      for (tag in names(cell_map)) {
        tag_map <- cell_map[[tag]]
        for (variant in names(tag_map)) {
          ks <- unique(as.integer(tag_map[[variant]]))
          ks <- ks[is.finite(ks)]
          if (!length(ks)) next
          for (k_val in ks) {
            out <- c(out, file.path(
              root_dir,
              paste0(cell, "_vae_joint_", tag, "_docs_peak_delta_fp_gene_fc_expr_", variant, "_K", k_val)
            ))
          }
        }
      }
    }
    unique(out)
  }
  extra_dirs <- .list_extra_vae_dirs(vae_root, vae_extra_k_map)
  extra_dirs <- extra_dirs[dir.exists(extra_dirs)]
  if (!length(extra_dirs)) {
    cli::cli_inform("Extra postprocess: no matching output dirs found under {vae_root}.")
  } else {
    vae_out_dirs_override <- extra_dirs
    vae_target_model_dirs_override <- extra_dirs
    run_vae_doc_topic_heatmaps_rerun <- TRUE
    run_vae_topic_delta_network_plots <- TRUE
    run_vae_topic_delta_network_pathway <- TRUE
    cli::cli_inform("Extra postprocess: {length(extra_dirs)} folder(s) queued.")
  }
}

run_vae_pathway_rerun <- FALSE
if (run_vae_pathway_rerun) {
  vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
  all_dirs <- list.dirs(vae_root, recursive = TRUE, full.names = TRUE)
  out_dirs <- all_dirs[file.exists(file.path(all_dirs, "link_topic_scores_baseline.csv"))]
  if (!length(out_dirs)) {
    out_dirs <- all_dirs[file.exists(file.path(all_dirs, "link_topic_scores_gate_peak_and_gene_in_set.csv"))]
    cli::cli_inform("VAE pathway rerun: no baseline link-score files found; falling back to gate file scan.")
  }
  if (!length(out_dirs)) {
    cli::cli_inform("VAE pathway rerun: no output dirs with link_topic_scores_* files under {vae_root}.")
    cli::cli_inform("VAE pathway rerun: set base_dir or check that link_topic_scores files exist.")
  } else {
    cli::cli_inform("VAE pathway rerun: {length(out_dirs)} target folder(s) under {vae_root}.")
  }

  rerun_one <- function(d) {
    rerun_pathway_from_link_scores(
      out_dir = d,
      link_scores_file = file.path(d, "link_topic_scores_baseline.csv"),
      tf_link_scores_file = if (!is.null(vae_pathway_link_scores_file_tf)) {
        file.path(d, vae_pathway_link_scores_file_tf)
      } else {
        NULL
      },
      gene_terms_file = if (!is.null(vae_pathway_link_gene_terms_file)) {
        file.path(d, vae_pathway_link_gene_terms_file)
      } else {
        NULL
      },
      allow_missing = TRUE,
      include_tf = vae_pathway_link_include_tf,
      include_gene = vae_pathway_link_include_gene,
      min_prob = vae_pathway_link_min_prob,
      gene_min_prob = vae_pathway_link_gene_min_prob,
      tf_min_prob = vae_pathway_link_tf_min_prob,
      tf_max_topics = vae_pathway_link_tf_max_topics,
      tf_top_n_per_topic = vae_pathway_link_tf_top_n_per_topic,
      top_n_per_topic = common_args$top_n_per_topic,
      max_pathways = common_args$max_pathways,
      per_comparison = vae_pathway_per_comparison,
      per_comparison_dir = vae_pathway_per_comparison_dir,
      split_direction = vae_pathway_split_direction,
      make_heatmap = FALSE
    )
  }

  mc_cores <- min(20L, length(out_dirs))
  if (.Platform$OS.type == "windows" || mc_cores < 2L) {
    lapply(out_dirs, rerun_one)
  } else {
    parallel::mclapply(out_dirs, rerun_one, mc.cores = mc_cores)
  }
}

run_vae_pathway_dotplot_replot <- FALSE
if (run_vae_pathway_dotplot_replot) {
  .assert_pkg("data.table")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_inform("Skipping VAE dotplot replot: ggplot2 not installed.")
  } else {
    vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
    dot_csvs <- list.files(vae_root, pattern = "_dotplot\\.csv$", recursive = TRUE, full.names = TRUE)
    if (!length(dot_csvs)) {
      cli::cli_inform("VAE dotplot replot: no *_dotplot.csv files under {vae_root}.")
    }
    replot_one <- function(csv_path) {
      dt <- data.table::fread(csv_path)
      if (!nrow(dt)) return(invisible(NULL))
      if (!("topic" %in% names(dt))) return(invisible(NULL))

      dt[, topic_num := as.integer(gsub("^Topic", "", topic))]
      dt[, score_val := ifelse(is.finite(combined_score), combined_score,
                        ifelse(is.finite(odds_ratio), odds_ratio, logp))]
      dt[, score_val := pmax(score_val, 0)]
      dt[, size_val := sqrt(score_val)]

      path_order <- dt[, .(
        min_topic = min(topic_num, na.rm = TRUE),
        max_score = max(score_val, na.rm = TRUE)
      ), by = pathway]
      path_order <- path_order[order(min_topic, -max_score)]
      dt[, pathway := factor(pathway, levels = rev(path_order$pathway))]
      dt[, topic := factor(topic, levels = paste0("Topic", sort(unique(topic_num))))]

      n_topics <- length(levels(dt$topic))
      n_paths <- length(levels(dt$pathway))
      p <- ggplot2::ggplot(
        dt,
        ggplot2::aes(x = topic, y = pathway, size = size_val, color = logp)
      ) +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::scale_color_gradient(low = "#2c7bb6", high = "#d7191c") +
        ggplot2::scale_size(range = if (n_paths > 80) c(0.6, 5) else c(1, 8)) +
        ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
        ggplot2::labs(
          x = "Topic",
          y = NULL,
          color = expression(-log[10]~"(adj.P)"),
          size = "Score"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
          axis.text.y = ggplot2::element_text(size = 8, face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
        )

      pdf_path <- sub("_dotplot\\.csv$", "_dotplot.pdf", csv_path)
      ggplot2::ggsave(
        pdf_path,
        p,
        width = max(8, n_topics * 0.6),
        height = max(6, min(160, n_paths * 0.25)),
        limitsize = FALSE
      )
    }

    mc_cores <- min(20L, length(dot_csvs))
    if (.Platform$OS.type == "windows" || mc_cores < 2L) {
      lapply(dot_csvs, replot_one)
    } else {
      parallel::mclapply(dot_csvs, replot_one, mc.cores = mc_cores)
    }
  }
}

run_vae_doc_topic_heatmaps_rerun <- FALSE
if (run_vae_doc_topic_heatmaps_rerun) {
  if (!exists("plot_tf_topic_heatmaps_from_link_scores")) {
    source(utils_path)
  }
  vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
  if (exists("vae_out_dirs_override") && length(vae_out_dirs_override)) {
    out_dirs <- vae_out_dirs_override
  } else {
    out_dirs <- list.dirs(vae_root, recursive = TRUE, full.names = TRUE)
    out_dirs <- out_dirs[file.exists(file.path(out_dirs, "link_topic_scores_baseline.csv"))]
  }

  rerun_one <- function(d) {
    if (!exists("plot_tf_topic_heatmaps_from_link_scores")) {
      source(utils_path)
    }
    if (!exists("plot_tf_topic_heatmaps_from_link_scores")) {
      stop("plot_tf_topic_heatmaps_from_link_scores not found after source(", utils_path, ").")
    }
    base_file <- file.path(d, "link_topic_scores_baseline.csv")
    if (file.exists(base_file)) {
      out_base <- file.path(d, "doc_topic_heatmaps_link_scores_baseline")
      plot_tf_topic_heatmaps_from_link_scores(
        link_scores = data.table::fread(base_file),
        out_dir = out_base,
        title_prefix = paste("Link scores baseline", basename(d), sep = " | "),
        value_col = "prob",
        min_value = 0,
        per_comparison = TRUE,
        split_direction = TRUE
      )
    }
    gate_file <- file.path(d, "link_topic_scores_gate_peak_and_gene_in_set.csv")
    if (file.exists(gate_file)) {
      out_gate <- file.path(d, "doc_topic_heatmaps_link_scores_gate_peak_and_gene_in_set")
      plot_tf_topic_heatmaps_from_link_scores(
        link_scores = data.table::fread(gate_file),
        out_dir = out_gate,
        title_prefix = paste("Link scores gate peak+gene in_set", basename(d), sep = " | "),
        value_col = "prob",
        min_value = 0,
        per_comparison = TRUE,
        split_direction = TRUE
      )
    }
  }

  mc_cores <- min(20L, length(out_dirs))
  if (.Platform$OS.type == "windows" || mc_cores < 2L) {
    lapply(out_dirs, rerun_one)
  } else {
    parallel::mclapply(out_dirs, rerun_one, mc.cores = mc_cores)
  }
}

run_vae_topic_delta_network_plots <- FALSE
if (run_vae_topic_delta_network_plots) {
  source(utils_path)
  source(file.path(dirname(utils_path), "utils_plot_tf_network_delta.R"))

  .assert_pkg("data.table")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    cli::cli_inform("Skipping topic delta network plots: htmlwidgets not installed.")
  } else {
    step2_out_dir <- file.path(base_dir, "connect_tf_target_genes")
    vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
    if (exists("vae_target_model_dirs_override") && length(vae_target_model_dirs_override)) {
      target_model_dirs <- vae_target_model_dirs_override
    } else {
      target_model_dirs <- list.dirs(vae_root, recursive = FALSE, full.names = TRUE)
      target_model_dirs <- target_model_dirs[grepl("(_ctf_docs_|_auto_ctf_)", basename(target_model_dirs))]
    }
    target_model_dirs <- sort(target_model_dirs)
    if (!length(target_model_dirs)) {
      cli::cli_abort("No *_ctf_docs_* folders found under {vae_root}")
    }

    for (target_model_dir in target_model_dirs) {
      cli::cli_inform("Delta network plots: processing {basename(target_model_dir)}")
      topic_terms_file <- file.path(target_model_dir, "topic_terms.csv")
      if (!file.exists(topic_terms_file)) {
        cli::cli_inform("topic_terms.csv not found in {target_model_dir}; skipping.")
        next
      }
      topic_terms <- data.table::fread(topic_terms_file)
      if (!("topic_num" %in% names(topic_terms))) {
        if ("topic" %in% names(topic_terms)) {
          topic_terms[, topic_num := as.integer(gsub("^Topic", "", topic))]
        } else {
          cli::cli_abort("topic_terms missing topic_num/topic.")
        }
      }
      topic_terms <- topic_terms[isTRUE(in_set) & grepl("^GENE:", term_id)]
      topic_terms[, gene_key := sub("^GENE:", "", term_id)]
      gene_terms_by_topic <- split(topic_terms$gene_key, topic_terms$topic_num)

      link_sources <- list(
        baseline = file.path(target_model_dir, "link_topic_scores_baseline.csv"),
        gate = file.path(target_model_dir, "link_topic_scores_gate_peak_and_gene_in_set.csv")
      )
      link_sources <- link_sources[file.exists(unlist(link_sources))]
      if (!length(link_sources)) {
        cli::cli_inform("No link_topic_scores_*.csv found in {target_model_dir}; skipping.")
        next
      }

    is_ctrl_tag <- function(tag) {
      ctrl_patterns <- c("(^|_)ctrl($|_)", "control", "baseline", "basal", "_10_fbs($|_)", "10fbs")
      grepl(paste(ctrl_patterns, collapse = "|"), tolower(tag), perl = TRUE)
    }
    detect_mapping <- function(df_names) {
      sc <- df_names[grepl("^link_score_", df_names)]
      if (length(sc) < 2L) stop("Could not find two 'link_score_*' columns.")
      sc <- sc[seq_len(min(2L, length(sc)))]
      suf <- sub("^link_score_", "", sc)
      ctrl_idx <- if (is_ctrl_tag(suf[1]) && !is_ctrl_tag(suf[2])) 1L else
        if (!is_ctrl_tag(suf[1]) && is_ctrl_tag(suf[2])) 2L else 2L
      str_idx <- if (ctrl_idx == 1L) 2L else 1L
      list(
        score_ctrl_col = sc[ctrl_idx],
        score_str_col  = sc[str_idx],
        sign_ctrl_col  = paste0("link_sign_", suf[ctrl_idx]),
        sign_str_col   = paste0("link_sign_", suf[str_idx]),
        tf_expr_ctrl_col   = paste0("tf_expr_",    suf[ctrl_idx]),
        tf_expr_str_col    = paste0("tf_expr_",    suf[str_idx]),
        gene_expr_ctrl_col = paste0("gene_expr_",  suf[ctrl_idx]),
        gene_expr_str_col  = paste0("gene_expr_",  suf[str_idx])
      )
    }

    delta_cache <- new.env(parent = emptyenv())
    fetch_delta_links <- function(comp) {
      if (exists(comp, envir = delta_cache, inherits = FALSE)) {
        return(get(comp, envir = delta_cache))
      }
      cand <- file.path(step2_out_dir, paste0(comp, "_delta_links_filtered.csv"))
      if (!file.exists(cand)) {
        cand <- file.path(step2_out_dir, paste0(comp, "_delta_links.csv"))
      }
      if (!file.exists(cand)) {
        cli::cli_inform("Missing delta links for comparison: {comp}")
        assign(comp, NULL, envir = delta_cache)
        return(NULL)
      }
      df <- readr::read_csv(cand, show_col_types = FALSE)
      if ("active_any" %in% names(df)) {
        df <- df[which(isTRUE(df$active_any) | df$active_any %in% c(1, "1", "TRUE", "true")), , drop = FALSE]
      }
      if (all(c("log2FC_tf_expr", "log2FC_gene_expr") %in% names(df))) {
        same_dir <- is.finite(df$log2FC_tf_expr) & is.finite(df$log2FC_gene_expr) &
          (df$log2FC_tf_expr * df$log2FC_gene_expr > 0)
        df <- df[which(same_dir), , drop = FALSE]
      }
      assign(comp, df, envir = delta_cache)
      df
    }

    for (src_name in names(link_sources)) {
      link_path <- link_sources[[src_name]]
      link_dt <- data.table::fread(link_path)
      if (!nrow(link_dt)) next

      if (!("topic_num" %in% names(link_dt))) {
        if ("topic" %in% names(link_dt)) {
          link_dt[, topic_num := as.integer(gsub("^Topic", "", topic))]
        } else {
          cli::cli_abort("link_scores missing topic_num/topic in {link_path}")
        }
      }
      if (!("prob" %in% names(link_dt))) link_dt[, prob := 1]

      doc_info <- .parse_doc_id(link_dt$doc_id)
      link_dt <- cbind(link_dt, doc_info)
      link_dt[, direction_group := .direction_group(direction)]
      link_dt <- link_dt[!is.na(comparison_id) & nzchar(comparison_id)]
      if (!nrow(link_dt)) next

      link_min_prob <- vae_pathway_link_tf_min_prob
      if (is.finite(vae_pathway_link_gene_min_prob)) {
        link_min_prob <- max(link_min_prob, vae_pathway_link_gene_min_prob, na.rm = TRUE)
      }

      out_root <- file.path(
        target_model_dir,
        if (src_name == "baseline") "doc_topic_sub_network_link_scores_baseline"
        else "doc_topic_sub_network_link_topic_scores_gate_peak_and_gene_in_set"
      )
      dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

      comps <- unique(link_dt$comparison_id)
      topics_all <- sort(unique(link_dt$topic_num))
      for (cmp in comps) {
        cmp_delta <- fetch_delta_links(cmp)
        if (is.null(cmp_delta) || !nrow(cmp_delta)) next
        ns <- names(cmp_delta)
        mp <- detect_mapping(ns)

        gene_col <- if ("gene_key" %in% ns) "gene_key" else if ("gene" %in% ns) "gene" else "target_gene"
        if (!gene_col %in% ns) {
          cli::cli_inform("No gene column found for {cmp}; skipping.")
          next
        }
        tf_col <- if ("tf" %in% ns) "tf" else "TF"
        if (!tf_col %in% ns) {
          cli::cli_inform("No TF column found for {cmp}; skipping.")
          next
        }

        for (dir_lab in unique(link_dt[comparison_id == cmp, direction_group])) {
          if (!nzchar(dir_lab)) next
          for (topic_id in topics_all) {
            link_sub <- link_dt[
              comparison_id == cmp & direction_group == dir_lab & topic_num == topic_id &
                prob >= link_min_prob & !is.na(tf) & nzchar(tf) & !is.na(gene_key) & nzchar(gene_key)
            ]
            if (!nrow(link_sub)) next
            link_pairs <- unique(link_sub[, .(tf, gene_key, prob, score)])
            if (!nrow(link_pairs)) next

            tf_list <- sort(unique(link_pairs$tf))
            gene_list <- sort(unique(link_pairs$gene_key))
            if (!length(tf_list) || !length(gene_list)) next

            cmp_dt <- data.table::as.data.table(cmp_delta)
            pair_dt <- data.table::data.table(tf = link_pairs$tf, gene = link_pairs$gene_key)
            data.table::setnames(pair_dt, c("tf", "gene"), c(tf_col, gene_col))
            sub_links <- data.table::merge.data.table(
              cmp_dt,
              pair_dt,
              by = c(tf_col, gene_col),
              all = FALSE
            )
            if (!nrow(sub_links)) next

            out_dir <- file.path(out_root, .safe_filename(cmp), .safe_filename(dir_lab))
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            input_file <- file.path(out_dir, paste0("Topic", topic_id, "_inputs.csv"))
            pairs_file <- file.path(out_dir, paste0("Topic", topic_id, "_link_pairs.csv"))
            tf_vals <- sort(unique(tf_list))
            gene_vals <- sort(unique(gene_list))
            sub_tf_vals <- sort(unique(sub_links[[tf_col]]))
            sub_gene_vals <- sort(unique(sub_links[[gene_col]]))
            input_dt <- data.frame(
              list_type = c(
                rep("tf_list", length(tf_vals)),
                rep("gene_list", length(gene_vals)),
                rep("sub_links_tf", length(sub_tf_vals)),
                rep("sub_links_gene", length(sub_gene_vals))
              ),
              value = c(tf_vals, gene_vals, sub_tf_vals, sub_gene_vals),
              stringsAsFactors = FALSE
            )
            data.table::fwrite(input_dt, input_file)
            data.table::fwrite(link_pairs, pairs_file)

            plot_title <- paste(cmp, dir_lab, paste0("Topic", topic_id), sep = " | ")
            if (!"size_by" %in% names(formals(plot_tf_network_delta))) {
              cli::cli_abort("plot_tf_network_delta does not support size_by; re-source utils_plot_tf_network_delta.R.")
            }
            args <- list(
              data = sub_links,
              plot_title = plot_title,
              layout_algo = "fr",
              physics = TRUE,
              add_direct = TRUE,
              edge_filter_min = 0,
              min_delta_abs = 0,
              keep_top_edges_per_tf = 6000,
              peak_mode = "show_all",
              show_peaks = FALSE,
              gene_fc_thresh = 1.5,
              de_reference = "str_over_ctrl",
              motif_db = "jaspar2024",
              score_ctrl_col = if (mp$score_ctrl_col %in% ns) mp$score_ctrl_col else NULL,
              score_str_col = if (mp$score_str_col %in% ns) mp$score_str_col else NULL,
              sign_ctrl_col = if (mp$sign_ctrl_col %in% ns) mp$sign_ctrl_col else NULL,
              sign_str_col = if (mp$sign_str_col %in% ns) mp$sign_str_col else NULL,
              tf_expr_ctrl_col = if (mp$tf_expr_ctrl_col %in% ns) mp$tf_expr_ctrl_col else NULL,
              tf_expr_str_col = if (mp$tf_expr_str_col %in% ns) mp$tf_expr_str_col else NULL,
              gene_expr_ctrl_col = if (mp$gene_expr_ctrl_col %in% ns) mp$gene_expr_ctrl_col else NULL,
              gene_expr_str_col = if (mp$gene_expr_str_col %in% ns) mp$gene_expr_str_col else NULL,
              size_by = "expr_max"
            )
            w <- try(do.call(plot_tf_network_delta, args), silent = TRUE)
            if (inherits(w, "try-error")) {
              err_file <- file.path(out_dir, paste0("Topic", topic_id, "_plot_error.txt"))
              err_msg <- conditionMessage(attr(w, "condition"))
              cat(sprintf("[%s] Plot failed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), err_msg),
                  file = err_file, append = TRUE)
              next
            }

            out_html <- file.path(out_dir, paste0("Topic", topic_id, ".html"))
            htmlwidgets::saveWidget(w, out_html, selfcontained = TRUE)
            if (exists(".set_html_title")) {
              .set_html_title(out_html, plot_title)
            }
            node_sizes <- attr(w, "node_sizes")
            if (!is.null(node_sizes)) {
              data.table::fwrite(
                node_sizes,
                file.path(out_dir, paste0("Topic", topic_id, "_node_sizes.csv"))
              )
            }
          }
        }
      }
    }
    }
  }
}

run_vae_topic_delta_network_pathway <- FALSE
if (run_vae_topic_delta_network_pathway) {
  .assert_pkg("data.table")
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    cli::cli_inform("Skipping subnet pathway rerun: enrichR not installed.")
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_inform("Skipping subnet pathway rerun: ggplot2 not installed.")
  } else {
    dbs <- c(
      "GO_Biological_Process_2023",
      "GO_Cellular_Component_2023",
      "GO_Molecular_Function_2023",
      "Reactome_2022",
      "WikiPathways_2024_Human"
    )
    db_short <- c(
      GO_Biological_Process_2023 = "GO:BP",
      GO_Cellular_Component_2023 = "GO:CC",
      GO_Molecular_Function_2023 = "GO:MF",
      Reactome_2022 = "Reactome",
      WikiPathways_2024_Human = "WikiPathways"
    )
    padj_cut <- 0.05
    min_genes <- 5L
    top_n_per_topic <- 100L
    max_pathways <- common_args$max_pathways

    tryCatch(
      {
        enrichR::setEnrichrSite("Enrichr")
      },
      error = function(e) {
        cli::cli_inform("Failed to set Enrichr site: {conditionMessage(e)}")
      }
    )

    vae_root <- file.path(base_dir, "benchmark_topic_vae_models")
    if (exists("vae_target_model_dirs_override") && length(vae_target_model_dirs_override)) {
      target_model_dirs <- vae_target_model_dirs_override
    } else {
      target_model_dirs <- list.dirs(vae_root, recursive = FALSE, full.names = TRUE)
      target_model_dirs <- target_model_dirs[grepl("(_ctf_docs_|_auto_ctf_)", basename(target_model_dirs))]
    }
    target_model_dirs <- sort(target_model_dirs)
    if (!length(target_model_dirs)) {
      cli::cli_abort("No *_ctf_docs_* folders found under {vae_root}")
    }

    .collect_enrichr <- function(gene_sets, log_fun) {
      if (!length(gene_sets)) return(data.table::data.table())
      res_list <- vector("list", length(gene_sets))
      names(res_list) <- names(gene_sets)
      for (nm in names(gene_sets)) {
        genes <- unique(gene_sets[[nm]])
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (length(genes) < as.integer(min_genes)) {
          log_fun(sprintf("Topic %s skipped: gene count < %d", nm, min_genes))
          next
        }
        log_fun(sprintf("Topic %s gene count: %d", nm, length(genes)))
        enr <- tryCatch(
          enrichR::enrichr(genes, dbs),
          error = function(e) {
            log_fun(sprintf("Topic %s enrichr error: %s", nm, conditionMessage(e)))
            NULL
          }
        )
        if (is.null(enr)) next
        rows <- lapply(names(enr), function(db) {
          df <- enr[[db]]
          if (is.null(df) || !nrow(df)) return(NULL)
          if (!("Adjusted.P.value" %in% names(df)) || !("Term" %in% names(df))) return(NULL)
          df <- df[is.finite(df$Adjusted.P.value) & df$Adjusted.P.value <= padj_cut, , drop = FALSE]
          if (!nrow(df)) return(NULL)
          db_label <- if (db %in% names(db_short)) db_short[[db]] else db
          term_clean <- gsub("\\s*\\([^)]*\\)$", "", df$Term)
          combined_score <- if ("Combined.Score" %in% names(df)) df$Combined.Score else NA_real_
          odds_ratio <- if ("Odds.Ratio" %in% names(df)) df$Odds.Ratio else NA_real_
          data.table::data.table(
            topic = as.integer(nm),
            pathway = paste(db_label, term_clean, sep = ": "),
            padj = as.numeric(df$Adjusted.P.value),
            pval = if ("P.value" %in% names(df)) as.numeric(df$P.value) else NA_real_,
            overlap = if ("Overlap" %in% names(df)) as.character(df$Overlap) else NA_character_,
            genes = if ("Genes" %in% names(df)) as.character(df$Genes) else NA_character_,
            logp = -log10(df$Adjusted.P.value),
            combined_score = as.numeric(combined_score),
            odds_ratio = as.numeric(odds_ratio)
          )
        })
        res_list[[nm]] <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
        log_fun(sprintf("Topic %s enriched pathways: %d (padj<=%s)", nm, nrow(res_list[[nm]]), padj_cut))
      }
      data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
    }

    .write_dotplot <- function(res_dt, dot_prefix, plot_title, log_fun) {
      plot_dt <- data.table::copy(res_dt)
      plot_dt[, topic_num := as.integer(topic)]
      plot_dt[, topic := paste0("Topic", topic_num)]
      plot_dt[, score_val := ifelse(is.finite(combined_score), combined_score,
                             ifelse(is.finite(odds_ratio), odds_ratio, logp))]
      plot_dt[, score_val := pmax(score_val, 0)]
      plot_dt[, size_val := sqrt(score_val)]
      path_order <- plot_dt[, .(
        min_topic = min(topic_num, na.rm = TRUE),
        max_score = max(score_val, na.rm = TRUE)
      ), by = pathway]
      path_order <- path_order[order(min_topic, -max_score)]
      plot_dt[, pathway := factor(pathway, levels = rev(path_order$pathway))]
      plot_dt[, topic := factor(topic, levels = paste0("Topic", sort(unique(topic_num))))]

      n_topics_plot <- length(unique(plot_dt$topic))
      n_paths_plot <- length(unique(plot_dt$pathway))
      size_range <- if (n_paths_plot > 80) c(0.6, 5) else c(1, 8)

      dot_path <- paste0(dot_prefix, ".pdf")
      dot_csv <- paste0(dot_prefix, ".csv")
      ord_dt <- plot_dt[order(
        match(pathway, levels(pathway)),
        match(topic, levels(topic))
      )]
      data.table::fwrite(ord_dt, dot_csv)

      max_path_chars <- suppressWarnings(max(nchar(as.character(plot_dt$pathway)), na.rm = TRUE))
      if (!is.finite(max_path_chars)) max_path_chars <- 0
      width <- max(8, n_topics_plot * 0.6 + max_path_chars * 0.12)
      height <- max(6, min(160, n_paths_plot * 0.25))
      p <- ggplot2::ggplot(
        plot_dt,
        ggplot2::aes(x = topic, y = pathway, size = size_val, color = logp)
      ) +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::scale_color_gradient(low = "#2c7bb6", high = "#d7191c") +
        ggplot2::scale_size(range = size_range) +
        ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
        ggplot2::labs(
          x = "Topic",
          y = NULL,
          color = expression(-log[10]~"(adj.P)"),
          size = "Score",
          title = plot_title
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
          axis.text.y = ggplot2::element_text(size = 8, face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
        )
      ggplot2::ggsave(dot_path, p, width = width, height = height, limitsize = FALSE)
      log_fun(sprintf("Dot plot CSV saved: %s", dot_csv))
      log_fun(sprintf("Dot plot PDF saved: %s", dot_path))
      invisible(TRUE)
    }

    for (target_model_dir in target_model_dirs) {
      cli::cli_inform("Delta network pathway: processing {basename(target_model_dir)}")
      input_roots <- list(
        baseline = file.path(target_model_dir, "doc_topic_sub_network_link_scores_baseline"),
        gate = file.path(target_model_dir, "doc_topic_sub_network_link_topic_scores_gate_peak_and_gene_in_set")
      )
      input_roots <- input_roots[dir.exists(unlist(input_roots))]
      if (!length(input_roots)) {
        cli::cli_inform("No doc_topic_sub_network_link_scores_* folders found in {target_model_dir}; skipping.")
        next
      }

      for (src_name in names(input_roots)) {
        in_root <- input_roots[[src_name]]
        input_files <- list.files(in_root, pattern = "Topic[0-9]+_inputs\\.csv$", recursive = TRUE, full.names = TRUE)
        if (!length(input_files)) next

        out_root <- file.path(
          target_model_dir,
          if (src_name == "baseline") "doc_topic_sub_network_pathway_link_scores_baseline"
          else "doc_topic_sub_network_pathway_link_topic_scores_gate_peak_and_gene_in_set"
        )
        dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
        log_path <- file.path(out_root, "subnet_pathway_debug.txt")
        log_msg <- function(msg) {
          stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
          cat(sprintf("[%s] %s\n", stamp, msg), file = log_path, append = TRUE)
        }

        file_info <- lapply(input_files, function(f) {
          rel <- sub(paste0("^", normalizePath(in_root, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(f, winslash = "/", mustWork = FALSE))
          parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
          if (length(parts) < 3L) return(NULL)
          topic_num <- as.integer(gsub("^Topic", "", sub("_inputs\\.csv$", "", parts[length(parts)])))
          list(
            file = f,
            comparison_id = parts[1],
            direction_group = parts[2],
            topic_num = topic_num
          )
        })
        file_info <- Filter(Negate(is.null), file_info)
        if (!length(file_info)) next
        file_dt <- data.table::rbindlist(file_info)

        for (cmp in unique(file_dt$comparison_id)) {
          for (dir_lab in unique(file_dt[comparison_id == cmp, direction_group])) {
            group_dt <- file_dt[comparison_id == cmp & direction_group == dir_lab]
            if (!nrow(group_dt)) next
            log_local <- function(msg) log_msg(sprintf("%s | %s: %s", cmp, dir_lab, msg))

            gene_sets <- list()
            for (i in seq_len(nrow(group_dt))) {
              row <- group_dt[i]
              dt <- data.table::fread(row$file)
              if (!nrow(dt)) next
              genes <- dt[list_type %in% c("gene_list", "sub_links_gene"), value]
              tfs <- dt[list_type %in% c("tf_list", "sub_links_tf"), value]
              genes <- unique(genes[!is.na(genes) & nzchar(genes)])
              tfs <- unique(tfs[!is.na(tfs) & nzchar(tfs)])
              if (isTRUE(vae_pathway_link_include_gene) && length(genes)) {
                gs <- genes
              } else {
                gs <- character(0)
              }
              if (isTRUE(vae_pathway_link_include_tf) && length(tfs)) {
                gs <- c(gs, tfs)
              }
              gene_sets[[as.character(row$topic_num)]] <- unique(gs)
            }
            if (!length(gene_sets)) {
              log_local("No gene sets; skipping.")
              next
            }

            res_sub <- .collect_enrichr(gene_sets, log_fun = log_local)
            if (!nrow(res_sub)) {
              log_local("No enriched terms at padj_cut; skipping.")
              next
            }
            res_sub <- res_sub[is.finite(logp) & logp > 0]
            if (!nrow(res_sub)) next
            res_sub <- res_sub[order(-logp), .SD[seq_len(min(.N, as.integer(top_n_per_topic)))], by = topic]
            if (nrow(res_sub)) {
              path_rank <- res_sub[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway]
              if (nrow(path_rank) > as.integer(max_pathways)) {
                keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway]
                res_sub <- res_sub[pathway %in% keep]
              }
            }
            res_sub[, topic := as.integer(topic)]
            data.table::setorder(res_sub, topic)

            out_dir <- file.path(out_root, .safe_filename(cmp))
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            plot_title <- paste(cmp, dir_lab, "Pathway enrichment (subnet inputs)", sep = " | ")
            prefix <- file.path(out_dir, paste0(.safe_filename(cmp), "_", .safe_filename(dir_lab), "_dotplot"))
            .write_dotplot(res_sub, dot_prefix = prefix, plot_title = plot_title, log_fun = log_local)
          }
        }
      }
    }
  }
}

run_integration <- FALSE
if (run_integration) {
  integration_dir <- file.path(topic_root, "integration")
  dir.create(integration_dir, showWarnings = FALSE, recursive = TRUE)

  integration_jobs <- list(
    list(cell = "HPAFII", k = 50L, gene_suffix = "gene_fc_expr", fp_suffix = "peak_delta_fp")
  )

  .theta_path_for_k <- function(run_dir, k) {
    rds_dir <- file.path(run_dir, "rds")
    if (!dir.exists(rds_dir)) return(NULL)
    preferred <- c(
      sprintf("theta_perplexity_K%d.rds", k),
      sprintf("theta_maximum_K%d.rds", k),
      sprintf("theta_derivative_K%d.rds", k),
      sprintf("theta_fallback_K%d.rds", k)
    )
    for (nm in preferred) {
      p <- file.path(rds_dir, nm)
      if (file.exists(p)) return(p)
    }
    other <- list.files(rds_dir, pattern = paste0("^theta_.*_K", k, "\\.rds$"), full.names = TRUE)
    if (length(other)) return(other[[1]])
    NULL
  }

  .normalize_rows <- function(mat) {
    mat <- as.matrix(mat)
    mat[!is.finite(mat) | mat < 0] <- 0
    rs <- rowSums(mat)
    rs[!is.finite(rs) | rs <= 0] <- 1
    mat / rs
  }

  .transform_log1p <- function(mat) {
    log1p(mat * 1000)
  }

  .transform_clr <- function(mat) {
    eps <- 1e-12
    logx <- log(mat + eps)
    logx - rowMeans(logx)
  }

  .transform_ilr <- function(mat) {
    clr <- .transform_clr(mat)
    k <- ncol(clr)
    if (k < 2L) return(clr)
    H <- stats::contr.helmert(k)
    Hn <- apply(H, 2, function(col) col / sqrt(sum(col^2)))
    clr %*% Hn
  }

  .mi_total <- function(x, y) {
    x <- .normalize_rows(x)
    y <- .normalize_rows(y)
    eps <- 1e-12
    x <- pmax(x, eps)
    y <- pmax(y, eps)
    x <- x / rowSums(x)
    y <- y / rowSums(y)
    x_m <- colMeans(x)
    y_m <- colMeans(y)
    joint <- t(x) %*% y / nrow(x)
    marg <- outer(x_m, y_m)
    sum(joint * log(joint / marg))
  }

  .mi_pointwise <- function(x, y) {
    x <- .normalize_rows(x)
    y <- .normalize_rows(y)
    eps <- 1e-12
    x <- pmax(x, eps)
    y <- pmax(y, eps)
    x <- x / rowSums(x)
    y <- y / rowSums(y)
    x_m <- colMeans(x)
    y_m <- colMeans(y)
    marg <- outer(x_m, y_m)
    vapply(seq_len(nrow(x)), function(i) {
      joint <- outer(x[i, ], y[i, ])
      sum(joint * log(joint / marg))
    }, numeric(1))
  }

  .relative_stats <- function(x, y) {
    x_norm <- sqrt(rowSums(x * x))
    y_norm <- sqrt(rowSums(y * y))
    denom <- x_norm + y_norm
    denom[!is.finite(denom) | denom <= 0] <- NA_real_
    eps <- 1e-12
    list(
      log2_ratio = log2((x_norm + eps) / (y_norm + eps)),
      rna_frac = x_norm / denom,
      fp_frac = y_norm / denom,
      rna_norm = x_norm,
      fp_norm = y_norm
    )
  }

  .get_named <- function(x, nm) {
    if (is.null(x) || is.null(names(x))) return(NA_real_)
    if (!nm %in% names(x)) return(NA_real_)
    as.numeric(x[[nm]])
  }

  .mi_grade <- function(mi_total, mi_summary) {
    if (!is.finite(mi_total)) return(NA_character_)
    p95 <- .get_named(mi_summary, "null_p95")
    p50 <- .get_named(mi_summary, "null_p50")
    if (!is.finite(p95) || !is.finite(p50)) return("unknown")
    if (mi_total >= p95) return("high")
    if (mi_total >= p50) return("moderate")
    "low"
  }

  .plot_topic_knn_pca <- function(theta_rna, theta_fp, out_file, title_prefix = NULL) {
    x_rna <- t(as.matrix(theta_rna))
    x_fp <- t(as.matrix(theta_fp))
    if (!nrow(x_rna) || !nrow(x_fp)) return(invisible(NULL))

    row_scale <- function(mat) {
      m <- rowMeans(mat)
      s <- apply(mat, 1, stats::sd)
      s[!is.finite(s) | s == 0] <- 1
      mat <- mat - m
      mat / s
    }

    x_rna <- row_scale(x_rna)
    x_fp <- row_scale(x_fp)
    x_all <- rbind(x_rna, x_fp)
    pca <- stats::prcomp(x_all, center = FALSE, scale. = FALSE)
    coords <- pca$x[, 1:2, drop = FALSE]
    n_rna <- nrow(x_rna)
    n_fp <- nrow(x_fp)
    idx_rna <- seq_len(n_rna)
    idx_fp <- n_rna + seq_len(n_fp)

    dist_all <- as.matrix(stats::dist(coords))
    dist_rf <- dist_all[idx_rna, idx_fp, drop = FALSE]
    nn_fp <- apply(dist_rf, 1, which.min)

    grDevices::pdf(out_file, width = 7, height = 6)
    tryCatch(
      {
        graphics::plot(coords[, 1], coords[, 2], type = "n", xlab = "PC1", ylab = "PC2",
                       main = title_prefix)
        graphics::segments(
          coords[idx_rna, 1], coords[idx_rna, 2],
          coords[idx_fp[nn_fp], 1], coords[idx_fp[nn_fp], 2],
          col = "#bdbdbd"
        )
        graphics::points(coords[idx_rna, 1], coords[idx_rna, 2], pch = 16, col = "#1f78b4")
        graphics::points(coords[idx_fp, 1], coords[idx_fp, 2], pch = 17, col = "#e31a1c")
        graphics::text(coords[idx_rna, 1], coords[idx_rna, 2], labels = idx_rna, pos = 3, cex = 0.8)
        graphics::text(coords[idx_fp, 1], coords[idx_fp, 2], labels = seq_len(n_fp), pos = 3, cex = 0.8)
        graphics::legend(
          "topright",
          legend = c("RNA topics", "FP topics"),
          pch = c(16, 17),
          col = c("#1f78b4", "#e31a1c"),
          bty = "n"
        )
        graphics::mtext("Lines connect each RNA topic to its nearest FP topic (KNN=1)", side = 1, line = 3, cex = 0.8)
      },
      finally = grDevices::dev.off()
    )
    invisible(TRUE)
  }

  for (job in integration_jobs) {
    run_rna <- file.path(topic_root, sprintf("%s_ctf_docs_%s_K%d", job$cell, job$gene_suffix, job$k))
    run_fp <- file.path(topic_root, sprintf("%s_ctf_docs_%s_K%d", job$cell, job$fp_suffix, job$k))
    theta_rna_path <- .theta_path_for_k(run_rna, job$k)
    theta_fp_path <- .theta_path_for_k(run_fp, job$k)
    if (is.null(theta_rna_path) || is.null(theta_fp_path)) {
      cli::cli_inform("Skipping integration: missing theta for {job$cell} K{job$k}.")
      next
    }
    theta_rna <- readRDS(theta_rna_path)
    theta_fp <- readRDS(theta_fp_path)
    if (is.null(rownames(theta_rna)) || is.null(rownames(theta_fp))) {
      cli::cli_inform("Skipping integration: theta missing rownames for {job$cell} K{job$k}.")
      next
    }
    docs <- intersect(rownames(theta_rna), rownames(theta_fp))
    if (length(docs) < 2L) {
      cli::cli_inform("Skipping integration: <2 shared docs for {job$cell} K{job$k}.")
      next
    }
    theta_rna <- theta_rna[docs, , drop = FALSE]
    theta_fp <- theta_fp[docs, , drop = FALSE]

    mi_total <- .mi_total(theta_rna, theta_fp)
    writeLines(as.character(mi_total), file.path(integration_dir, sprintf("%s_K%d_mi_total.txt", job$cell, job$k)))

    set.seed(123)
    n_null <- 100L
    mi_null_raw <- vapply(seq_len(n_null), function(i) {
      theta_fp_shuf <- theta_fp[sample(nrow(theta_fp)), , drop = FALSE]
      .mi_total(theta_rna, theta_fp_shuf)
    }, numeric(1))
    mi_null <- mi_null_raw[is.finite(mi_null_raw)]
    null_n_total <- length(mi_null_raw)
    null_n_finite <- length(mi_null)
    mi_summary <- c(
      observed = mi_total,
      null_mean = if (null_n_finite) mean(mi_null) else NA_real_,
      null_sd = if (null_n_finite > 1) stats::sd(mi_null) else NA_real_,
      null_p05 = if (null_n_finite) as.numeric(stats::quantile(mi_null, 0.05)) else NA_real_,
      null_p50 = if (null_n_finite) as.numeric(stats::quantile(mi_null, 0.5)) else NA_real_,
      null_p95 = if (null_n_finite) as.numeric(stats::quantile(mi_null, 0.95)) else NA_real_,
      empirical_p = if (null_n_finite) mean(mi_null >= mi_total) else NA_real_,
      null_n_total = null_n_total,
      null_n_finite = null_n_finite
    )
    mi_z <- NA_real_
    null_sd <- .get_named(mi_summary, "null_sd")
    null_mean <- .get_named(mi_summary, "null_mean")
    if (is.finite(null_sd) && null_sd > 0 && is.finite(null_mean)) {
      mi_z <- (mi_total - null_mean) / null_sd
    }
    mi_grade <- .mi_grade(mi_total, mi_summary)
    writeLines(
      paste(names(mi_summary), signif(mi_summary, 4), sep = ": "),
      file.path(integration_dir, sprintf("%s_K%d_mi_total_summary.txt", job$cell, job$k))
    )
    grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_mi_total_null.pdf", job$cell, job$k)), width = 7, height = 5)
    tryCatch(
      {
        null_p50 <- .get_named(mi_summary, "null_p50")
        null_p95 <- .get_named(mi_summary, "null_p95")
        xlim <- range(c(mi_null, mi_total, null_p50, null_p95), finite = TRUE)
        if (!all(is.finite(xlim))) xlim <- NULL
        if (length(mi_null) >= 2L) {
          graphics::hist(mi_null, breaks = 30, col = "#b3cde3", border = "white",
                         main = sprintf("%s K%d MI null", job$cell, job$k),
                         xlab = "Mutual information (null)",
                         xlim = xlim)
          graphics::abline(v = mi_total, col = "#e31a1c", lwd = 2)
        } else {
          graphics::plot.new()
          graphics::title(main = sprintf("%s K%d MI null", job$cell, job$k))
          graphics::text(0.5, 0.5, "Not enough finite null values to plot")
        }
        if (is.finite(null_p50)) {
          graphics::abline(v = null_p50, col = "#1f78b4", lwd = 2, lty = 2)
        }
        if (is.finite(null_p95)) {
          graphics::abline(v = null_p95, col = "#33a02c", lwd = 2, lty = 2)
        }
        graphics::legend(
          "topright",
          legend = c("observed", "null p50", "null p95"),
          col = c("#e31a1c", "#1f78b4", "#33a02c"),
          lwd = 2,
          lty = c(1, 2, 2),
          bty = "n"
        )
        graphics::mtext(
          sprintf("observed=%.3f, p=%.3f, grade=%s", mi_total, mi_summary[["empirical_p"]], mi_grade),
          side = 3,
          line = 0.2
        )
        graphics::mtext(
          "Observed outside null => strong coupling; inside null => weak coupling",
          side = 3,
          line = 1.2,
          cex = 0.8
        )
      },
      finally = grDevices::dev.off()
    )

    pointwise <- .mi_pointwise(theta_rna, theta_fp)
    data.table::fwrite(
      data.table::data.table(doc_id = docs, pointwise_mi = pointwise),
      file.path(integration_dir, sprintf("%s_K%d_mi_pointwise.csv", job$cell, job$k))
    )
    grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_mi_pointwise_hist.pdf", job$cell, job$k)), width = 7, height = 5)
    tryCatch(
      {
        graphics::hist(pointwise, breaks = 30, col = "#ccebc5", border = "white",
                       main = sprintf("%s K%d pointwise MI", job$cell, job$k),
                       xlab = "Pointwise MI")
        pw_med <- stats::median(pointwise, na.rm = TRUE)
        pw_p90 <- stats::quantile(pointwise, 0.9, na.rm = TRUE)
        if (is.finite(pw_med)) graphics::abline(v = pw_med, col = "#1f78b4", lwd = 2)
        if (is.finite(pw_p90)) graphics::abline(v = pw_p90, col = "#33a02c", lwd = 2, lty = 2)
        graphics::legend(
          "topright",
          legend = c("median", "90%"),
          col = c("#1f78b4", "#33a02c"),
          lwd = 2,
          lty = c(1, 2),
          bty = "n"
        )
        graphics::mtext("Right tail/high values = strong doc-level alignment", side = 3, line = 0.2, cex = 0.8)
      },
      finally = grDevices::dev.off()
    )
    ord <- order(pointwise, decreasing = TRUE)
    top_n <- min(20, length(ord))
    top_df <- data.table::data.table(doc_id = docs[ord][seq_len(top_n)], pointwise_mi = pointwise[ord][seq_len(top_n)])
    data.table::fwrite(top_df, file.path(integration_dir, sprintf("%s_K%d_mi_pointwise_top20.csv", job$cell, job$k)))
    grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_mi_pointwise_top20.pdf", job$cell, job$k)),
                   width = 9, height = max(6, top_n * 0.3))
    tryCatch(
      {
        graphics::par(mar = c(5, 16, 3, 2))
        graphics::barplot(
          rev(top_df$pointwise_mi),
          names.arg = rev(top_df$doc_id),
          las = 1,
          horiz = TRUE,
          col = "#8dd3c7",
          main = sprintf("%s K%d top pointwise MI", job$cell, job$k),
          xlab = "Pointwise MI"
        )
        graphics::mtext("Higher = stronger cross-modality alignment", side = 3, line = 0.2, cex = 0.8)
      },
      finally = grDevices::dev.off()
    )

    corr <- stats::cor(theta_rna, theta_fp, use = "pairwise.complete.obs")
    corr_dt <- data.table::as.data.table(corr, keep.rownames = "rna_topic")
    data.table::fwrite(
      corr_dt,
      file.path(integration_dir, sprintf("%s_K%d_topic_cross_correlation.csv", job$cell, job$k))
    )
    .plot_topic_knn_pca(
      theta_rna,
      theta_fp,
      out_file = file.path(integration_dir, sprintf("%s_K%d_topic_knn_pca.pdf", job$cell, job$k)),
      title_prefix = sprintf("%s K%d topic PCA (RNA vs FP)", job$cell, job$k)
    )
    best_idx <- apply(corr, 1, which.max)
    best_corr <- corr[cbind(seq_len(nrow(corr)), best_idx)]
    best_df <- data.table::data.table(
      rna_topic = rownames(corr),
      fp_topic = colnames(corr)[best_idx],
      corr = as.numeric(best_corr)
    )
    data.table::fwrite(
      best_df,
      file.path(integration_dir, sprintf("%s_K%d_topic_cross_correlation_best.csv", job$cell, job$k))
    )
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(
        corr,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        border_color = NA,
        main = sprintf("%s K%d topic cross-correlation\nHigher positive = aligned topics", job$cell, job$k),
        filename = file.path(integration_dir, sprintf("%s_K%d_topic_cross_correlation.pdf", job$cell, job$k))
      )
    }
    grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_topic_cross_correlation_best.pdf", job$cell, job$k)),
                   width = 7, height = max(5, nrow(best_df) * 0.25))
    tryCatch(
      {
        graphics::par(mar = c(5, 7, 3, 2))
        graphics::barplot(
          best_df$corr,
          names.arg = best_df$rna_topic,
          las = 2,
          col = "#80b1d3",
          main = sprintf("%s K%d best topic correlations", job$cell, job$k),
          ylab = "Correlation"
        )
        graphics::abline(h = 0, col = "grey60")
        graphics::mtext("Higher = stronger topic pairing across modes", side = 3, line = 0.2, cex = 0.8)
      },
      finally = grDevices::dev.off()
    )

    transforms <- list(
      log1p = .transform_log1p,
      clr = .transform_clr,
      ilr = .transform_ilr
    )
    for (nm in names(transforms)) {
      tr_fn <- transforms[[nm]]
      x_t <- tr_fn(theta_rna)
      y_t <- tr_fn(theta_fp)
      joint <- cbind(x_t, y_t)
      joint_path <- file.path(integration_dir, sprintf("%s_K%d_joint_%s.csv", job$cell, job$k, nm))
      data.table::fwrite(
        data.table::data.table(doc_id = docs, joint),
        joint_path
      )
      pca <- stats::prcomp(joint, rank. = 2, center = TRUE, scale. = TRUE)
      pca_df <- data.table::data.table(doc_id = docs, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
      data.table::fwrite(
        pca_df,
        file.path(integration_dir, sprintf("%s_K%d_joint_%s_pca.csv", job$cell, job$k, nm))
      )
      rel_stats <- .relative_stats(x_t, y_t)
      rel <- rel_stats$log2_ratio
      data.table::fwrite(
        data.table::data.table(
          doc_id = docs,
          log2_ratio = rel_stats$log2_ratio,
          rna_frac = rel_stats$rna_frac,
          fp_frac = rel_stats$fp_frac,
          rna_norm = rel_stats$rna_norm,
          fp_norm = rel_stats$fp_norm
        ),
        file.path(integration_dir, sprintf("%s_K%d_joint_%s_relative_norms.csv", job$cell, job$k, nm))
      )
      grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_joint_%s_pca.pdf", job$cell, job$k, nm)), width = 7, height = 6)
      tryCatch(
        {
          cols <- grDevices::colorRampPalette(c("#2c7bb6", "#f7f7f7", "#d7191c"))(100)
          rel_rank <- rank(rel, ties.method = "average")
          rel_col <- cols[as.integer((rel_rank - 1) / max(rel_rank) * 99) + 1]
          rel_min <- min(rel, na.rm = TRUE)
          rel_max <- max(rel, na.rm = TRUE)
          graphics::plot(
            pca_df$PC1, pca_df$PC2,
            col = rel_col, pch = 16, cex = 0.7,
            xlab = "PC1", ylab = "PC2",
            main = sprintf("%s K%d joint PCA (%s)", job$cell, job$k, nm)
          )
          graphics::mtext(
            "Color = log2(norm RNA / norm FP); red RNA-dominant, blue FP-dominant",
            side = 1,
            line = 4,
            cex = 0.8
          )
          graphics::mtext(
            "Separation by color => modality dominance; mixed colors => balanced",
            side = 1,
            line = 5,
            cex = 0.8
          )
          if (is.finite(rel_min) && is.finite(rel_max) && rel_min != rel_max) {
            q <- stats::quantile(rel, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
            q_idx <- floor((q - rel_min) / (rel_max - rel_min) * 99) + 1
            q_idx[!is.finite(q_idx)] <- 1
            graphics::legend(
              "topright",
              legend = sprintf("q%d=%.2f", c(10, 50, 90), q),
              col = cols[q_idx],
              pch = 16,
              bty = "n",
              title = "log2_ratio"
            )
          }
        },
        finally = grDevices::dev.off()
      )
      grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_joint_%s_relative_norms.pdf", job$cell, job$k, nm)), width = 7, height = 5)
      tryCatch(
        {
          graphics::hist(rel, breaks = 30, col = "#fdb462", border = "white",
                         main = sprintf("%s K%d relative norms (%s)", job$cell, job$k, nm),
                         xlab = "log2(norm RNA / norm FP)")
          graphics::abline(v = 0, col = "#1f78b4", lwd = 2)
          rel_med <- stats::median(rel, na.rm = TRUE)
          if (is.finite(rel_med)) {
            graphics::mtext(sprintf("median=%.2f", rel_med), side = 1, line = 4, cex = 0.8)
          }
          graphics::mtext(
            "Closer to 0 = balanced (often preferred); larger |value| = dominance",
            side = 1,
            line = 5,
            cex = 0.8
          )
        },
        finally = grDevices::dev.off()
      )
      grDevices::pdf(file.path(integration_dir, sprintf("%s_K%d_joint_%s_relative_balance.pdf", job$cell, job$k, nm)), width = 7, height = 5)
      tryCatch(
        {
          graphics::hist(rel_stats$rna_frac, breaks = 30, col = "#b2df8a", border = "white",
                         main = sprintf("%s K%d RNA fraction (%s)", job$cell, job$k, nm),
                         xlab = "rna_frac = ||RNA|| / (||RNA|| + ||FP||)")
          graphics::abline(v = 0.5, col = "#1f78b4", lwd = 2, lty = 2)
          frac_med <- stats::median(rel_stats$rna_frac, na.rm = TRUE)
          if (is.finite(frac_med)) {
            graphics::mtext(sprintf("median=%.2f", frac_med), side = 1, line = 4, cex = 0.8)
          }
          graphics::mtext(
            "Closer to 0.5 = balanced (often preferred); farther = dominance",
            side = 1,
            line = 5,
            cex = 0.8
          )
        },
        finally = grDevices::dev.off()
      )
    }

    summary_path <- file.path(integration_dir, sprintf("%s_K%d_integration_summary.txt", job$cell, job$k))
    summary_lines <- c(
      sprintf("MI observed: %.4f", mi_total),
      sprintf("MI null mean: %.4f", .get_named(mi_summary, "null_mean")),
      sprintf("MI null 95%%: %.4f", .get_named(mi_summary, "null_p95")),
      sprintf("MI empirical p: %.4f", .get_named(mi_summary, "empirical_p")),
      sprintf("MI null n (total/finite): %s/%s", null_n_total, null_n_finite),
      sprintf("MI z-score: %.3f", mi_z),
      sprintf("MI grade vs null: %s", mi_grade),
      sprintf("Pointwise MI median: %.4f", stats::median(pointwise)),
      sprintf("Pointwise MI 90%%: %.4f", stats::quantile(pointwise, 0.9)),
      sprintf("Best topic correlation median: %.4f", stats::median(best_df$corr)),
      sprintf("RNA fraction median (log1p): %.4f", stats::median(.relative_stats(.transform_log1p(theta_rna), .transform_log1p(theta_fp))$rna_frac, na.rm = TRUE))
    )
    writeLines(summary_lines, summary_path)

    summary_pdf <- file.path(integration_dir, sprintf("%s_K%d_integration_summary.pdf", job$cell, job$k))
    grDevices::pdf(summary_pdf, width = 7, height = 5)
    tryCatch(
      {
        graphics::plot.new()
        graphics::title(main = sprintf("%s K%d integration summary", job$cell, job$k))
        y_pos <- seq(0.9, 0.1, length.out = length(summary_lines))
        for (i in seq_along(summary_lines)) {
          graphics::text(0, y_pos[i], summary_lines[i], adj = 0, cex = 0.9)
        }
      },
      finally = grDevices::dev.off()
    )
  }
}

# Integration A - Joint model (single WarpLDA on genes+peaks)
run_joint_warplda <- FALSE
if (run_joint_warplda) {
  joint_formals <- names(formals(run_tfdocs_warplda_one_option))
  missing_joint <- setdiff(c("joint_weight_fp", "joint_weight_gene", "joint_balance"), joint_formals)
  if (length(missing_joint)) {
    cli::cli_abort(
      "Loaded run_tfdocs_warplda_one_option is missing joint args: {paste(missing_joint, collapse = ', ')}. Re-source {utils_path}."
    )
  }

  joint_doc_mode <- list(
    name = "ctf_docs",
    doc_mode = "tf_cluster",
    direction_by = "gene",
    distinct_terms = TRUE,
    tf_cluster_map = tf_cluster_map,
    tf_exclude = tf_exclude
  )

  .short_k_label <- function(k_grid, default_grid = NULL) {
    k_grid <- unique(as.integer(k_grid))
    k_grid <- k_grid[is.finite(k_grid)]
    if (!length(k_grid)) return("Knone")
    if (length(k_grid) == 1L) return(paste0("K", k_grid))
    if (!is.null(default_grid)) {
      d <- unique(as.integer(default_grid))
      d <- d[is.finite(d)]
      if (length(d) && setequal(k_grid, d)) return("Kgrid_default")
    }
    if (length(k_grid) > 6L) return(paste0("Kgrid_n", length(k_grid)))
    paste0("K", paste(k_grid, collapse = "-"))
  }

  joint_jobs <- list(
    list(
      cell = "HPAFII",
      k_grid = K_grid_default,
      fp_weight = "delta_fp",
      gene_weight = "fc_mag_gene",
      balance = "min",
      direction_by = "gene"
    ),
    list(
      cell = "HPAFII",
      k_grid = c(11L),
      fp_weight = "delta_fp",
      gene_weight = "fc_mag_gene",
      balance = "min",
      direction_by = "gene"
    ),
    list(
      cell = "HPAFII",
      k_grid = c(15L),
      fp_weight = "delta_fp",
      gene_weight = "fc_mag_gene",
      balance = "min",
      direction_by = "gene"
    ),
    list(
      cell = "HPAFII",
      k_grid = c(70L),
      fp_weight = "delta_fp",
      gene_weight = "fc_mag_gene",
      balance = "min",
      direction_by = "gene"
    )
  )

  joint_label_cleaner <- function(label) {
    parts <- data.table::tstrsplit(label, "::", fixed = TRUE)
    base <- parts[[1]]
    dir <- if (length(parts) >= 2L) parts[[2]] else NA_character_
    base <- sub("^[^_]+_", "", base)
    base <- sub("_vs_[^_]+_10_FBS", "", base)
    out <- base
    has_dir <- !is.na(dir) & nzchar(dir)
    out[has_dir] <- paste(base[has_dir], dir[has_dir], sep = "::")
    out
  }

  run_joint_job <- function(job) {
    tryCatch(
      {
        edges_sub <- edges_all[cellline == job$cell]
        if (!nrow(edges_sub)) return(NULL)
        if (!is.null(joint_doc_mode$tf_exclude) && length(joint_doc_mode$tf_exclude)) {
          dt <- data.table::as.data.table(edges_sub)
          dt <- dt[!toupper(tf) %in% joint_doc_mode$tf_exclude]
          edges_sub <- dt
        }
        fp_suffix <- if (job$fp_weight == "fc_mag_fp") "peak_fc_fp" else "peak_delta_fp"
        gene_suffix <- "gene_fc_expr"
        k_label <- .short_k_label(job$k_grid, K_grid_default)
        dir_name <- sprintf(
          "%s_%s_joint_%s_%s_%s",
          job$cell,
          joint_doc_mode$name,
          fp_suffix,
          gene_suffix,
          k_label
        )
        out_dir <- file.path(topic_root, dir_name)

        local_common_args <- common_args
        if (!is.null(job$k_grid) && length(job$k_grid)) {
          local_common_args$K_grid <- as.integer(job$k_grid)
        }

        direction_by <- if (!is.null(job$direction_by)) job$direction_by else joint_doc_mode$direction_by
        args <- c(list(edges_all = edges_sub), local_common_args, list(
          out_dir = out_dir,
          option_label = "joint",
          doc_mode = joint_doc_mode$doc_mode,
          direction_by = direction_by,
          distinct_terms = joint_doc_mode$distinct_terms,
          joint_weight_fp = job$fp_weight,
          joint_weight_gene = job$gene_weight,
          joint_balance = job$balance,
          topic_by_comparison_label_cleaner = joint_label_cleaner
        ))
        if (!is.null(joint_doc_mode$tf_cluster_map)) {
          args$tf_cluster_map <- joint_doc_mode$tf_cluster_map
        }
        do.call(run_tfdocs_warplda_one_option, args)
      },
      error = function(e) {
        cli::cli_inform("Joint job failed: {conditionMessage(e)}")
        NULL
      }
    )
  }

  if (length(joint_jobs)) {
    mc_cores <- min(20L, length(joint_jobs))
    parallel::mclapply(joint_jobs, run_joint_job, mc.cores = mc_cores)
  }
}

# Integration B - Joint model with modality-split docs (RNA vs FP)
run_joint_warplda_split_docs <- FALSE
if (run_joint_warplda_split_docs) {
  joint_doc_mode <- list(
    name = "ctf_docs",
    doc_mode = "tf_cluster",
    direction_by = "gene",
    distinct_terms = TRUE,
    tf_cluster_map = tf_cluster_map,
    tf_exclude = tf_exclude
  )

  .short_k_label <- function(k_grid, default_grid = NULL) {
    k_grid <- unique(as.integer(k_grid))
    k_grid <- k_grid[is.finite(k_grid)]
    if (!length(k_grid)) return("Knone")
    if (length(k_grid) == 1L) return(paste0("K", k_grid))
    if (!is.null(default_grid)) {
      d <- unique(as.integer(default_grid))
      d <- d[is.finite(d)]
      if (length(d) && setequal(k_grid, d)) return("Kgrid_default")
    }
    if (length(k_grid) > 6L) return(paste0("Kgrid_n", length(k_grid)))
    paste0("K", paste(k_grid, collapse = "-"))
  }

  split_job <- list(
    cell = "HPAFII",
    doc_label = "ctf_split_docs",
    direction_by = "gene",
    k_grid = K_grid_default,
    fp_weight = "delta_fp",
    gene_weight = "fc_mag_gene",
    balance = "none"
  )

  split_label_cleaner <- function(label) {
    parts <- data.table::tstrsplit(label, "::", fixed = TRUE)
    base <- parts[[1]]
    dir <- if (length(parts) >= 2L) parts[[2]] else NA_character_
    base <- sub("^[^_]+_", "", base)
    base <- sub("_vs_[^_]+_10_FBS", "", base)
    out <- base
    has_dir <- !is.na(dir) & nzchar(dir)
    out[has_dir] <- paste(base[has_dir], dir[has_dir], sep = "::")
    out
  }

  split_jobs <- list(split_job)

  run_split_job <- function(job) {
    tryCatch(
      {
        edges_sub <- edges_all[cellline == job$cell]
        if (!nrow(edges_sub)) return(NULL)
        if (!is.null(joint_doc_mode$tf_exclude) && length(joint_doc_mode$tf_exclude)) {
          dt <- data.table::as.data.table(edges_sub)
          dt <- dt[!toupper(tf) %in% joint_doc_mode$tf_exclude]
          edges_sub <- dt
        }
        edges_gene <- data.table::as.data.table(edges_sub)
        edges_gene[, doc_modality := "GENE"]
        edges_peak <- data.table::as.data.table(edges_sub)
        edges_peak[, doc_modality := "PEAK"]
        edges_split <- data.table::rbindlist(list(edges_gene, edges_peak), use.names = TRUE, fill = TRUE)

        fp_suffix <- if (job$fp_weight == "fc_mag_fp") "peak_fc_fp" else "peak_delta_fp"
        gene_suffix <- "gene_fc_expr"
        k_label <- .short_k_label(job$k_grid, K_grid_default)
        dir_name <- sprintf(
          "%s_%s_joint_%s_%s_%s",
          job$cell,
          job$doc_label,
          fp_suffix,
          gene_suffix,
          k_label
        )
        out_dir <- file.path(topic_root, dir_name)

        local_common_args <- common_args
        if (!is.null(job$k_grid) && length(job$k_grid)) {
          local_common_args$K_grid <- as.integer(job$k_grid)
        }
        local_common_args$min_df <- 1L

        direction_by <- if (!is.null(job$direction_by)) job$direction_by else joint_doc_mode$direction_by
        args <- c(list(edges_all = edges_split), local_common_args, list(
          out_dir = out_dir,
          option_label = "joint",
          doc_mode = joint_doc_mode$doc_mode,
          direction_by = direction_by,
          distinct_terms = joint_doc_mode$distinct_terms,
          joint_weight_fp = job$fp_weight,
          joint_weight_gene = job$gene_weight,
          joint_balance = job$balance,
          topic_by_comparison_label_cleaner = split_label_cleaner
        ))
        if (!is.null(joint_doc_mode$tf_cluster_map)) {
          args$tf_cluster_map <- joint_doc_mode$tf_cluster_map
        }
        do.call(run_tfdocs_warplda_one_option, args)
      },
      error = function(e) {
        cli::cli_inform("Split-doc job failed: {conditionMessage(e)}")
        NULL
      }
    )
  }

  if (length(split_jobs)) {
    mc_cores <- min(20L, length(split_jobs))
    parallel::mclapply(split_jobs, run_split_job, mc.cores = mc_cores)
  }
}
