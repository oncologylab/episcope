library(episcope)
load_episcope_config("episcope_grn.yaml")

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



fp_manifest <- load_footprints(
  root_dir   = "Y:/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
  db_name    = "JASPAR2024",
  out_dir    = "Z:/episcope_test/nutrient_stress/fp_jaspar2024",
  n_workers  = 20,
  # n_samples  = 2,

  verbose    = TRUE
)
readr::write_csv(fp_manifest, "Z:/episcope_test/nutrient_stress/fp_jaspar2024_manifest.csv")


fp_manifest <- load_footprints(
  root_dir   = "/data/homes/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
  db_name    = "HOCOMOCOv13",
  out_dir    = "Z:/episcope_test/nutrient_stress_test/hocomocov13",
  n_workers  = 20,
  n_samples  = 2,
  verbose    = TRUE
)
readr::write_csv(fp_manifest, "Z:/episcope_test/nutrient_stress_test/fp_hocomocov13_manifest.csv")


fp_manifest <- load_footprints(
  root_dir   = "Z:/fptools_test/GSE96771_ATAC/TOBIAS_merged_peaks/",
  db_name    = "JASPAR2024",
  out_dir    = "Z:/episcope_test/GSE96771_ATAC",
  n_workers  = 20,
  n_samples  = 2,
  verbose    = TRUE
)
readr::write_csv(fp_manifest, "Z:/episcope_test/GSE96771_ATAC/fp_jaspar2024_manifest.csv")




# (Optional)
fp_manifest <- fp_manifest_trim(fp_manifest)  # renames files on disk by default
# Overwrite every annotation CSV referenced by the manifest:
summary_tbl <- fp_manifest_trim_annots(fp_manifest, n_workers = 18, verbose = TRUE)

# Inspect what changed:
dplyr::count(summary_tbl, status)
sum(summary_tbl$n_fixed, na.rm = TRUE)

options(future.globals.maxSize = 16 * 1024^3)
# Align the peaks based on the peak similarity
fp_aligned <- align_footprints(fp_manifest,
                                        mid_slop        = 10L, # midpoint tolerance (bp)
                                        round_digits    = 1L,  # round scores before comparing vectors
                                        score_match_pct = 0.8) # fraction of samples that must match (<=1, >=0)



length(unique(fp_aligned$id_map$peak_ID))
length(unique(fp_aligned$id_map$fp_peak_bak))
readr::write_csv(fp_aligned$fp_bound, paste0("inst/extdata/fp_bounds_", db, ".csv"))
readr::write_csv(fp_aligned$fp_score, paste0("inst/extdata/fp_scores_", db, ".csv"))
readr::write_csv(fp_aligned$fp_annotation, paste0("inst/extdata/fp_annotation_", db, ".csv"))


fp_score_raw <- fp_aligned$fp_score
fp_aligned_normalized <- fp_aligned
fp_aligned_normalized$fp_score <- qn_footprints(fp_aligned_normalized$fp_score, id_col = "peak_ID")
readr::write_csv(fp_aligned_normalized$fp_score, paste0("inst/extdata/fp_scores_qn_", db, ".csv"))

fp_aligned_normalized_manifest <- save_footprints(fp_aligned_normalized, out_dir  = paste0("inst/extdata/fp_aligned_normalized_", db))
readr::write_csv(fp_aligned_normalized_manifest, paste0("inst/extdata/fp_aligned_normalized_manifest_", db, ".csv"))




# Load RNA data and ATAC data
atac_data    <- readr::read_tsv("inst/extdata/All_ATAC.master_table.narrowpeaks.10mil.txt")

# atac_data    <- readr::read_tsv(file.path(base_dir, "CD34_BM_ATAC.master_table.narrowpeaks.10mil.txt"))


atac_out     <- make_atac_tibbles(atac_data, sort_peaks = TRUE)
atac_score   <- atac_out$score
atac_overlap <- atac_out$overlap

# rna <- readr::read_csv("inst/extdata/HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")
rna <- readr::read_csv("inst/extdata/HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")

rna <- clean_hgnc(rna) # clean "HGNC" column
# Filter rna expression genes/tfs needs to reach the threshold in at least 1 sample (group size)
rna <- filter_rna_expr(rna, tf_list, hgnc_col = "HGNC", gene_min = threshold_gene_expr, tf_min = threshold_tf_expr, min_samples = 1L)
# make two RNA tables with columns already named by ATAC ids ──

# strict_rna
smap <- dplyr::transmute(stric_metadata, old = strict_match_rna, new = id)
strict_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
nm <- names(strict_rna); nm[match(smap$old, nm)] <- smap$new
strict_rna <- strict_rna |> `names<-`(nm) |> dplyr::as_tibble()


# Build RNA and ATAC data list

strict_build_args <- list(
  atac_score    = atac_score,
  atac_overlap  = atac_overlap,
  rna           = strict_rna,
  metadata      = stric_metadata,
  tf_list       = tf_list,
  motif_db      = motif_db,
  label_col     = "strict_match_rna",
  expected_n    = 23
)



# RNA expression and ATAC peak based, footprint filtering
fp_aligned_normalized_filtered_manifest <- filter_footprints(
  fp_manifest         = fp_aligned_normalized_manifest,  # from stream_write_overviews_by_motif() returns
  out_dir             = file.path("inst","extdata",paste0("fp_aligned_normalized_filtered_", db)),
  build_args          = strict_build_args,
  n_workers           = 18,
  skip_existing       = FALSE,
  threshold_tf_expr   = threshold_tf_expr,
  verbose             = TRUE
)
readr::write_csv(fp_aligned_normalized_filtered_manifest, paste0("inst/extdata/fp_aligned_normalized_filtered_manifest_", db, ".csv"))



#
fp_corr_manifest <- process_motifs_tf_corr_in_parallel(
  fp_manifest = fp_aligned_normalized_filtered_manifest,
  out_dir     = paste0("inst/extdata/fp_aligned_normalized_filtered_corr_", db),
  build_args  = strict_build_args,
  n_workers   = 10,
  cor_method  = "pearson",
  min_non_na  = 5L
)
readr::write_csv(fp_corr_manifest, paste0("inst/extdata/fp_aligned_normalized_filtered_corr_manifest_", db, ".csv"))

