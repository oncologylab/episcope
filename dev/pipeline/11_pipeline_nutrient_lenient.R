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

load_config("dev/config/pdac_nutrient_stress_lenient_jaspar2024.yaml")

# lenient -----------------------------------------------------------------
# Requires: omics_data_lenient (built from lenient metadata) and omics_data (for variance)
basal_links_dir <- if (exists("basal_links_dir")) basal_links_dir else base_dir
basal_links_step2 <- readr::read_csv(file.path(basal_links_dir, "connect_tf_target_genes", "basal_links_tmp/step2_overall_tf_gene_links.csv"))

step2_out_dir <- file.path(base_dir, "connect_tf_target_genes_lenient")

omics_data_lenient <- grn_add_rna_expressed(omics_data_lenient, label_col = "cell_stress_type", threshold_gene_expr = threshold_gene_expr)
omics_data_lenient <- grn_add_fp_score_condition(omics_data_lenient, label_col = "cell_stress_type")
omics_data_lenient <- grn_add_fp_bound_condition(omics_data_lenient, label_col = "cell_stress_type", threshold_fp_score = threshold_fp_score)
if (is.null(omics_data_lenient$fp_variance) || is.null(omics_data_lenient$rna_variance)) {
  hv_variance <- precompute_hvf_hvg_variance(omics_data, cores = 20, chunk_size = 50000L)
  omics_data_lenient$fp_variance <- hv_variance$fp_variance
  omics_data_lenient$rna_variance <- hv_variance$rna_variance
}

omics_data_cond <- prepare_grn_set_for_light_by_condition(
  omics_data_lenient,
  label_col = "cell_stress_type"
)

light_by_condition(
  ds = omics_data_cond,
  basal_links = basal_links_step2,
  out_dir = step2_out_dir,
  prefix = "step2",
  label_col = "cell_stress_type",
  link_score_threshold = 0,
  fp_score_threshold = 0,
  tf_expr_threshold = 0,
  fp_bound_tbl = omics_data_lenient$fp_bound_condition,
  rna_expressed_tbl = omics_data_lenient$rna_expressed,
  require_fp_bound = TRUE,
  require_gene_expr = TRUE,
  gene_expr_threshold = 1L,
  filter_active = FALSE,
  use_parallel = TRUE,
  workers = 8,
  fp_variance_tbl = omics_data_lenient$fp_variance,
  rna_variance_tbl = omics_data_lenient$rna_variance
)
step3_out_dir <- file.path(base_dir, "diff_grn_lenient")

step3_specs <- build_cellwise_contrasts_from_index(
  index_csv = file.path(step2_out_dir, "step2_per_condition_index.csv"),
  out_dir = step2_out_dir,
  prefix = "step2",
  ctrl_tag = "Ctrl",
  clean_names = FALSE
)

delta_link <- 1
regulated_genes <- 1.5 # 1.5 2

run_links_deltas_driver(
  specs       = step3_specs,
  clean_names = FALSE,
  parallel    = TRUE,
  restrict_to_active_both = FALSE,
  edge_change_min = delta_link,
  keep_all_cols = TRUE
)

step3_delta_csvs <- list.files(step2_out_dir, "_delta_links.csv", full.names = TRUE)
step3_de_gene_log2_abs_min <- if (regulated_genes == 1.5) 0.585 else if (regulated_genes == 2) 1 else NA_real_

step3_bulk <- episcope::filter_links_deltas_bulk(
  step3_delta_csvs,
  gene_expr_min = threshold_gene_expr,
  tf_expr_min = threshold_tf_expr,
  fp_min = threshold_fp_score,
  link_min = threshold_link_score,
  abs_delta_min = delta_link,
  apply_de_gene = TRUE,
  de_gene_log2_abs_min = step3_de_gene_log2_abs_min,
  enforce_link_expr_sign = TRUE,
  expr_dir_col = "log2FC_gene_expr",
  workers = 20
)
