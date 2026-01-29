library(episcope)
load_config("episcope_grn.yaml")
source("R/utils_connect_tf_enhancers_to_target_genes.R")
source("R/utils_load_footprints_process.R")

# 1. stream_write_overviews_by_motif
# 2. trim_fp_manifest
# 3. trim_manifest_annots
# 4. preprocess_align_fp_peaks
# 5. quantile_normalize_fp_unique
# 6. save_fp_aligned_normalized
# 7. process_motifs_in_parallel (filtering)

# stream_write_overviews_by_motif → load_footprints()
# trim_fp_manifest             → fp_manifest_trim()
# trim_manifest_annots         → fp_manifest_trim_annots()
# preprocess_align_fp_peaks    → align_footprints()
# quantile_normalize_fp_unique → qn_footprints()
# save_fp_aligned_normalized   → save_footprints()

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

base_dir <- "/data/homes/yl814/episcope_test/nutrient_stress"


# Step 1. Predict TF binding sites ----------------------------------------
# fp_manifest <- load_footprints(
#   root_dir   = "/data/homes/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
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

# options(future.globals.maxSize = 32 * 1024^3)

# Align the peaks based on the peak similarity
# fp_aligned <- align_footprints(fp_manifest,
#   mid_slop        = 10L, # midpoint tolerance (bp)
#   round_digits    = 1L, # round scores before comparing vectors
#   score_match_pct = 0.8
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

# fp_aligned_normalized_manifest <- save_footprints(fp_aligned_normalized, out_dir = file.path(base_dir, sprintf("fp_aligned_normalized_%s", db)))
# readr::write_csv(fp_aligned_normalized_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_manifest_%s.csv", db)))

sample_metadata <- readxl::read_excel(file.path(base_dir, "sample_metadata.xlsx"), na = "NA")
strict_metadata   <- sample_metadata |> dplyr::filter(!is.na(strict_match_rna))
strict_metadata_AsPC1 <- strict_metadata |> dplyr::filter(cell == "AsPC1")
strict_metadata_HPAFII <- strict_metadata |> dplyr::filter(cell == "HPAFII")
strict_metadata_Panc1 <- strict_metadata |> dplyr::filter(cell == "Panc1")

# sample_metadata$id <- sample_metadata$ID
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

# Build RNA and ATAC data list
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
fp_aligned_normalized_filtered_manifest <- readr::read_csv(file.path(base_dir, sprintf("fp_aligned_normalized_filtered_manifest_%s.csv", db)))

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
# fp_corr_manifest <- readr::read_csv(file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))

# tf_corr_qc(fp_corr_manifest)

# (Optional correlate with all tfbs, not only canonical tfbs)
# tfs <- tf_list[tf_list %in% unique(strict_rna$HGNC)]

# for (tf in tfs) {
#   tf_vec <- strict_rna |> dplyr::filter(HGNC == tf) |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
#   res <- tf_corr_footprints_all_tfbs(
#     fp_manifest = fp_aligned_normalized_filtered_manifest,
#     tf_name  = tf,
#     tf_expr  = tf_vec,
#     out_dir  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
#   )
# }

# tf_vec <- strict_rna |> dplyr::filter(HGNC == "FOXA2") |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
# res <- tf_corr_footprints_all_tfbs(
#   fp_manifest = fp_aligned_normalized_filtered_manifest,
#   tf_name  = "FOXA2",
#   tf_expr  = tf_vec,
#   out_dir  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
# )
# res$files
# HNF1A
# HNF4G
# FOXA1
# HNF4A
# tfs <- c("IRF1", "RARG", "SOX9", "KLF5", "FOXA2")




# path to file
# f <- "~/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs/HNF4A_overview.txt"

# read + count unique TFBS coordinates
# x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
# n_unique_tfbs <- nrow(unique(x[c("TFBS_chr", "TFBS_start", "TFBS_end")]))
# nrow(unique(fp_score[c("peak_ID")]))


# do not do any filtering for benchmarking
# combined <- tf_corr_footprints_filter(fp_corr_manifest, p_thr = 0, r_thr = 0, output_bed = file.path(base_dir, sprintf("fp_predicted_tfbs_%s", db)))

# fp_score <- combined$fp_score
# fp_bound <- combined$fp_bound
# fp_annotation <- combined$fp_annotation

# readr::write_csv(fp_score, file.path(base_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db)))
# readr::write_csv(fp_bound, file.path(base_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db)))
# readr::write_csv(fp_annotation, file.path(base_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db)))

# length(unique(fp_annotation$fp_peak))

fp_score <- readr::read_csv(file.path(base_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db)))
fp_bound <- readr::read_csv(file.path(base_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db)))
fp_annotation <- readr::read_csv(file.path(base_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db)))



# strict grn_set to build the basal network
grn_set <- build_grn_set(
  fp_score      = fp_score,
  fp_bound      = fp_bound,
  fp_annotation = fp_annotation,
  atac_score    = atac_score,
  atac_overlap  = atac_overlap,
  rna           = strict_rna,
  metadata      = strict_metadata,
  tf_list       = tf_list,
  motif_db      = motif_db,
  label_col     = "strict_match_rna",
  expected_n    = 23
)
grn_set_HPAFII <- build_grn_set(
  fp_score      = fp_score,
  fp_bound      = fp_bound,
  fp_annotation = fp_annotation,
  atac_score    = atac_score,
  atac_overlap  = atac_overlap,
  rna           = strict_rna,
  metadata      = strict_metadata_HPAFII,
  tf_list       = tf_list,
  motif_db      = motif_db,
  label_col     = "strict_match_rna",
  expected_n    = 23
)
grn_set_AsPC1 <- build_grn_set(
  fp_score      = fp_score,
  fp_bound      = fp_bound,
  fp_annotation = fp_annotation,
  atac_score    = atac_score,
  atac_overlap  = atac_overlap,
  rna           = strict_rna,
  metadata      = strict_metadata_AsPC1,
  tf_list       = tf_list,
  motif_db      = motif_db,
  label_col     = "strict_match_rna",
  expected_n    = 23
)
grn_set_Panc1 <- build_grn_set(
  fp_score      = fp_score,
  fp_bound      = fp_bound,
  fp_annotation = fp_annotation,
  atac_score    = atac_score,
  atac_overlap  = atac_overlap,
  rna           = strict_rna,
  metadata      = strict_metadata_Panc1,
  tf_list       = tf_list,
  motif_db      = motif_db,
  label_col     = "strict_match_rna",
  expected_n    = 23
)

grn_set_cellline_list <- list(grn_set_HPAFII=grn_set_HPAFII, grn_set_AsPC1=grn_set_AsPC1, grn_set_Panc1=grn_set_Panc1)

# grn_set <- filter_by_min_bound(grn_set, min_bound = 1L)


# Step 2. Connect TF-occupied enhancers to target genes -------------------

# gene_annot_ref_hg38 <- episcope::episcope_build_gene_annot("hg38")
# gene_annot_ref_mm10 <- episcope_build_gene_annot("mm10")

# atac_peak to target gene corr

# GeneHancer based
gh_std   <- load_genehancer_panc(file.path("inst","extdata","GeneHancer_v5.24_elite_panc.csv")) # GeneHancer_v5.24_elite_panc.csv

# 30kb/25kb/20kb based
# gh_std <- episcope_make_windowed_gh(
#   peaks      = grn_set$atac_score,
#   genes      = unique(c(grn_set$rna$HGNC, grn_set$rna$ensembl_gene_id)),
#   flank_bp   = 30000,
#   mode       = "TSS",
#   gene_annot = gene_annot_ref_hg38,
#   id_col     = "HGNC" # or "ensembl_gene_id"
# )

options(future.globals.maxSize = 64 * 1024^3)

# For HiC links
# for (gh_std_name in names(gh_std_hic_list)) {
#   gh_std <- gh_std_hic_list[[gh_std_name]]
#
#   print(gh_std_name)
#   atac_res <- correlate_atac_to_genes(
#     grn_set  = grn_set,
#     gh_tbl   = gh_std,
#     gene_mode = "both",
#     fdr       = threshold_atac_gene_corr_p,
#     r_abs_min = threshold_atac_gene_corr_abs_r,
#     cache_dir = "/data/homes/yl814/episcope_cache/atac_gene_corr",
#     cache_tag = "nutrient_hic",
#     workers   = 30,
#     cache_verbose = TRUE
#   )
#   fp_res_full <- correlate_fp_to_genes(
#     grn_set             = grn_set,
#     atac_gene_corr_kept = atac_res$atac_gene_corr_full,
#     fdr                 = threshold_fp_gene_corr_p,
#     r_abs_min           = threshold_fp_gene_corr_abs_r,
#     method              = "pearson",
#     workers             = 20,
#     cache_dir           = "/data/homes/yl814/episcope_cache/fp_gene_corr",
#     cache_tag           = "nutrient_hic",
#     cache_chunk_size    = 5000L,
#     cache_verbose       = TRUE
#   )
#   readr::write_csv(atac_res$atac_gene_corr_full, file.path(base_dir, sprintf("atac_gene_corr_full_%s.csv", gh_std_name)))
#   readr::write_csv(fp_res_full$fp_gene_corr_full, file.path(base_dir, sprintf("fp_gene_corr_full_%s.csv", gh_std_name)))
# }


# For Each cell line
for (grn_set_name in names(grn_set_cellline_list)) {
  print(grn_set_name)
  grn_set_current <- grn_set_cellline_list[[grn_set_name]]
  atac_res <- correlate_atac_to_genes(
    grn_set  = grn_set_current,
    gh_tbl   = gh_std,
    gene_mode = "both",
    fdr       = threshold_atac_gene_corr_p,
    r_abs_min = threshold_atac_gene_corr_abs_r,
    cache_dir = "/data/homes/yl814/episcope_cache/atac_gene_corr",
    cache_tag = sprintf("nutrient_%s",grn_set_name),
    workers   = 30,
    cache_verbose = TRUE
  )
  fp_res_full <- correlate_fp_to_genes(
    grn_set             = grn_set_current,
    atac_gene_corr_kept = atac_res$atac_gene_corr_full,
    fdr                 = threshold_fp_gene_corr_p,
    r_abs_min           = threshold_fp_gene_corr_abs_r,
    method              = "pearson",
    workers             = 20,
    cache_dir           = "/data/homes/yl814/episcope_cache/fp_gene_corr",
    cache_tag           = sprintf("nutrient_%s",grn_set_name),
    cache_chunk_size    = 5000L,
    cache_verbose       = TRUE
  )
  readr::write_csv(atac_res$atac_gene_corr_full, file.path(base_dir, sprintf("atac_gene_corr_full_%s.csv", grn_set_name)))
  readr::write_csv(fp_res_full$fp_gene_corr_full, file.path(base_dir, sprintf("fp_gene_corr_full_%s.csv", grn_set_name)))
}

