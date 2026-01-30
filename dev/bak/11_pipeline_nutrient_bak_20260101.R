library(episcope)
load_config("episcope_grn.yaml")

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

base_dir <- "/data/homes/yl814/episcope_test/nutrient_stress"


# Step 1. Predict TF binding sites ----------------------------------------
fp_manifest <- load_footprints(
  root_dir   = "/data/homes/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
  db_name    = "JASPAR2024",
  out_dir    = file.path(base_dir, "fp_jaspar2024"),
  n_workers  = 20,
  # n_samples  = 2,
  verbose    = TRUE
)
# readr::write_csv(fp_manifest, file.path(base_dir, "fp_jaspar2024_manifest.csv"))

fp_manifest <- readr::read_csv(file.path(base_dir, "cache", "fp_jaspar2024_manifest.csv"))



if (db == "hocomocov13") {
  # (Optional, when using HOCOMOCO database)
  fp_manifest <- fp_manifest_trim(fp_manifest) # renames files on disk by default
  # Overwrite every annotation CSV referenced by the manifest:
  summary_tbl <- fp_manifest_trim_annots(fp_manifest, n_workers = 18, verbose = TRUE)

  # Inspect what changed:
  dplyr::count(summary_tbl, status)
  sum(summary_tbl$n_fixed, na.rm = TRUE)
}



options(future.globals.maxSize = 32 * 1024^3)

# Align the peaks based on the peak similarity
fp_aligned <- align_footprints(fp_manifest,
  mid_slop        = 10L, # midpoint tolerance (bp)
  round_digits    = 1L, # round scores before comparing vectors
  score_match_pct = 0.8
) # fraction of samples that must match (<=1, >=0)


length(unique(fp_aligned$id_map$peak_ID))
length(unique(fp_aligned$id_map$fp_peak_bak))

# readr::write_csv(fp_aligned$fp_bound, file.path(base_dir, sprintf("cache/fp_bounds_%s.csv", db)))
# readr::write_csv(fp_aligned$fp_score, file.path(base_dir, sprintf("cache/fp_scores_%s.csv", db)))
# readr::write_csv(fp_aligned$fp_annotation, file.path(base_dir, sprintf("cache/fp_annotation_%s.csv", db)))
fp_aligned <- list()
fp_aligned$fp_bound <- readr::read_csv(file.path(base_dir, sprintf("cache/fp_bounds_%s.csv", db)))
fp_aligned$fp_score <- readr::read_csv(file.path(base_dir, sprintf("cache/fp_scores_%s.csv", db)))
fp_aligned$fp_annotation <- readr::read_csv(file.path(base_dir, sprintf("cache/fp_annotation_%s.csv", db)))


fp_score_raw <- fp_aligned$fp_score
fp_aligned_normalized <- fp_aligned
fp_aligned_normalized$fp_score <- qn_footprints(fp_aligned_normalized$fp_score, id_col = "peak_ID")
readr::write_csv(fp_aligned_normalized$fp_score, file.path(base_dir, sprintf("fp_scores_qn_%s.csv", db)))

fp_aligned_normalized_manifest <- save_footprints(fp_aligned_normalized, out_dir = file.path(base_dir, sprintf("fp_aligned_normalized_%s", db)))
readr::write_csv(fp_aligned_normalized_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_manifest_%s.csv", db)))

sample_metadata <- readxl::read_excel(file.path(base_dir, "sample_metadata.xlsx"), na = "NA")
strict_metadata <- sample_metadata |> dplyr::filter(!is.na(strict_match_rna))
lenient_metadata <- sample_metadata |> dplyr::filter(run_grn)
# strict_metadata_AsPC1 <- strict_metadata |> dplyr::filter(cell == "AsPC1")
# strict_metadata_HPAFII <- strict_metadata |> dplyr::filter(cell == "HPAFII")
# strict_metadata_Panc1 <- strict_metadata |> dplyr::filter(cell == "Panc1")

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

# lenient_rna
smap <- dplyr::transmute(lenient_metadata, old = broad_match_rna, new = id)
lenient_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
nm <- names(lenient_rna); nm[match(smap$old, nm)] <- smap$new
lenient_rna <- lenient_rna |> `names<-`(nm) |> dplyr::as_tibble()

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
fp_aligned_normalized_filtered_manifest <- filter_footprints(
  fp_manifest         = fp_aligned_normalized_manifest,
  out_dir             = file.path(base_dir, sprintf("fp_aligned_normalized_filtered_%s", db)),
  build_args          = rna_atac_build_args,
  n_workers           = 10,
  skip_existing       = FALSE,
  threshold_tf_expr   = threshold_tf_expr,
  verbose             = TRUE
)

readr::write_csv(fp_aligned_normalized_filtered_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_filtered_manifest_%s.csv", db)))
# fp_aligned_normalized_filtered_manifest <- readr::read_csv(file.path("/data/homes/yl814/episcope_test/nutrient_stress", sprintf("cache/fp_aligned_normalized_filtered_manifest_%s.csv", db)))
old_root <- "/data/homes/yl814/episcope_test/nutrient_stress/"
new_root <- "/data/homes/yl814/episcope_test/nutrient_stress/cache/"

fix_path <- function(x) {
  ifelse(is.na(x), x, sub(old_root, new_root, x, fixed = TRUE))
}

fp_aligned_normalized_filtered_manifest$score <- fix_path(fp_aligned_normalized_filtered_manifest$score)
fp_aligned_normalized_filtered_manifest$bound <- fix_path(fp_aligned_normalized_filtered_manifest$bound)
fp_aligned_normalized_filtered_manifest$annot <- fix_path(fp_aligned_normalized_filtered_manifest$annot)


# TF expr to footprints
fp_corr_manifest <- tf_corr_footprints(
  fp_manifest = fp_aligned_normalized_filtered_manifest,
  out_dir     = file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_%s", db)),
  build_args  = rna_atac_build_args,
  n_workers   = 15,
  cor_method  = "kendall",
  min_non_na  = 5L
)


readr::write_csv(fp_corr_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))
# fp_corr_manifest <- readr::read_csv(file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))

tf_corr_qc(fp_corr_manifest)

# (Optional correlate with all tfbs, not only canonical tfbs)
tfs <- tf_list[tf_list %in% unique(strict_rna$HGNC)]

for (tf in tfs) {
  tf_vec <- strict_rna |> dplyr::filter(HGNC == tf) |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
  res <- tf_corr_footprints_all_tfbs(
    fp_manifest = fp_aligned_normalized_filtered_manifest,
    tf_name  = tf,
    tf_expr  = tf_vec,
    out_dir  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
  )
}

c("HNF4A", "MAFF", "ZEB1")

tf_vec <- strict_rna |> dplyr::filter(HGNC == "HNF4A") |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
res <- tf_corr_footprints_all_tfbs(
  fp_manifest = fp_aligned_normalized_filtered_manifest,
  tf_name  = "HNF4A",
  tf_expr  = tf_vec,
  cor_method = "spearman", # c("pearson","spearman","kendall"),
  out_dir  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs_spearman"
)
res$files

# HNF1A
# HNF4G
# FOXA1
# HNF4A
tfs <- c("IRF1", "RARG", "SOX9", "KLF5", "FOXA2")




# path to file
f <- "~/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs/HNF4A_overview.txt"

# read + count unique TFBS coordinates
x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
n_unique_tfbs <- nrow(unique(x[c("TFBS_chr", "TFBS_start", "TFBS_end")]))
nrow(unique(fp_score[c("peak_ID")]))


# do not do any filtering for benchmarking
combined <- tf_corr_footprints_filter(fp_corr_manifest, p_thr = 0, r_thr = 0, output_bed = file.path(base_dir, sprintf("fp_predicted_tfbs_%s", db)))

fp_score <- combined$fp_score
fp_bound <- combined$fp_bound
fp_annotation <- combined$fp_annotation

readr::write_csv(fp_score, file.path(base_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db)))
readr::write_csv(fp_bound, file.path(base_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db)))
readr::write_csv(fp_annotation, file.path(base_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db)))

length(unique(fp_annotation$fp_peak))

fp_score <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db)))
fp_bound <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db)))
fp_annotation <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db)))



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

# Minimal helper: per-condition RNA expression flags (1/0) using a threshold.
make_rna_expressed <- function(rna_tbl, threshold_gene_expr = 10) {
  stopifnot(is.data.frame(rna_tbl))
  stopifnot(is.numeric(threshold_gene_expr), length(threshold_gene_expr) == 1L)
  sample_cols <- setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC"))
  rna_tbl |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(sample_cols),
        ~ as.integer(.x >= threshold_gene_expr)
      )
    )
}

grn_set$rna_expressed <- make_rna_expressed(grn_set$rna, threshold_gene_expr = 10)

# Minimal helper: refine FP bound by FP score threshold and ATAC overlap.
refine_fp_bound <- function(fp_bound_tbl, fp_score_tbl, atac_overlap_tbl, fp_annotation_tbl, threshold_fp_score) {
  stopifnot(is.data.frame(fp_bound_tbl), is.data.frame(fp_score_tbl), is.data.frame(atac_overlap_tbl))
  stopifnot(is.data.frame(fp_annotation_tbl))
  stopifnot(is.numeric(threshold_fp_score), length(threshold_fp_score) == 1L)

  id_bound <- names(fp_bound_tbl)[1]
  id_score <- names(fp_score_tbl)[1]
  id_atac <- names(atac_overlap_tbl)[1]

  sample_cols <- intersect(intersect(names(fp_bound_tbl), names(fp_score_tbl)), names(atac_overlap_tbl))
  sample_cols <- setdiff(sample_cols, c(id_bound, id_score, id_atac))
  if (!length(sample_cols)) stop("No shared sample columns found for fp_bound/fp_score/atac_overlap.")

  idx_score <- match(fp_bound_tbl[[id_bound]], fp_score_tbl[[id_score]])
  if (anyNA(idx_score)) stop("fp_score rows do not align with fp_bound by ID.")

  score_tbl <- fp_score_tbl[idx_score, sample_cols, drop = FALSE]

  fp_to_atac <- fp_annotation_tbl |>
    dplyr::distinct(.data$fp_peak, .data$atac_peak)
  idx_atac <- match(fp_bound_tbl[[id_bound]], fp_to_atac$fp_peak)
  if (anyNA(idx_atac)) stop("fp_bound peaks missing in fp_annotation (fp_peak).")
  atac_peaks <- fp_to_atac$atac_peak[idx_atac]
  idx_atac_tbl <- match(atac_peaks, atac_overlap_tbl[[id_atac]])
  if (anyNA(idx_atac_tbl)) stop("Mapped atac_peak missing in atac_overlap.")
  atac_tbl <- atac_overlap_tbl[idx_atac_tbl, sample_cols, drop = FALSE]

  out <- fp_bound_tbl
  for (col in sample_cols) {
    out[[col]] <- as.integer(
      out[[col]] > 0 &
        score_tbl[[col]] >= threshold_fp_score &
        atac_tbl[[col]] > 0
    )
  }
  out
}

grn_set$fp_bound <- refine_fp_bound(
  fp_bound_tbl = grn_set$fp_bound,
  fp_score_tbl = grn_set$fp_score,
  atac_overlap_tbl = grn_set$atac_overlap,
  fp_annotation_tbl = grn_set$fp_annotation,
  threshold_fp_score = threshold_fp_score
)

sum(tapply(rowSums(grn_set$fp_bound[,-1, drop = FALSE]) > 0, grn_set$fp_bound$peak_ID, any))


fp_score <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_scores_qn_%s.csv", "jaspar2024")))
fp_bound <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_bounds_%s.csv", "jaspar2024")))
fp_annotation <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_annotation_%s.csv", "jaspar2024")))


grn_set_lenient <- build_grn_set(
  fp_score      = fp_score,
  fp_bound      = fp_bound,
  fp_annotation = fp_annotation,
  atac_score    = atac_score,
  atac_overlap  = atac_overlap,
  rna           = lenient_rna,
  metadata      = lenient_metadata,
  tf_list       = tf_list,
  motif_db      = motif_db,
  label_col     = "cell_stress_type",
  expected_n    = 23
)
grn_set_lenient <- filter_by_min_bound(grn_set_lenient, min_bound = 1L)

# grn_set_HPAFII <- build_grn_set(
#   fp_score      = fp_score,
#   fp_bound      = fp_bound,
#   fp_annotation = fp_annotation,
#   atac_score    = atac_score,
#   atac_overlap  = atac_overlap,
#   rna           = strict_rna,
#   metadata      = strict_metadata_HPAFII,
#   tf_list       = tf_list,
#   motif_db      = motif_db,
#   label_col     = "strict_match_rna",
#   expected_n    = 23
# )
# grn_set_AsPC1 <- build_grn_set(
#   fp_score      = fp_score,
#   fp_bound      = fp_bound,
#   fp_annotation = fp_annotation,
#   atac_score    = atac_score,
#   atac_overlap  = atac_overlap,
#   rna           = strict_rna,
#   metadata      = strict_metadata_AsPC1,
#   tf_list       = tf_list,
#   motif_db      = motif_db,
#   label_col     = "strict_match_rna",
#   expected_n    = 23
# )
# grn_set_Panc1 <- build_grn_set(
#   fp_score      = fp_score,
#   fp_bound      = fp_bound,
#   fp_annotation = fp_annotation,
#   atac_score    = atac_score,
#   atac_overlap  = atac_overlap,
#   rna           = strict_rna,
#   metadata      = strict_metadata_Panc1,
#   tf_list       = tf_list,
#   motif_db      = motif_db,
#   label_col     = "strict_match_rna",
#   expected_n    = 23
# )

grn_set_cellline_list <- list(grn_set_HPAFII=grn_set_HPAFII, grn_set_AsPC1=grn_set_AsPC1, grn_set_Panc1=grn_set_Panc1)

grn_set <- filter_by_min_bound(grn_set, min_bound = 1L)


# Step 2. Connect TF-occupied enhancers to target genes -------------------

gene_annot_ref_hg38 <- episcope::episcope_build_gene_annot("hg38")
# gene_annot_ref_mm10 <- episcope_build_gene_annot("mm10")

# atac_peak to target gene corr

# GeneHancer based
gh_std   <- load_genehancer_panc(file.path("inst","extdata","GeneHancer_v5.24_elite_panc.csv")) # GeneHancer_v5.24_elite_panc.csv

# 30kb/25kb/20kb based
gh_std <- episcope_make_windowed_gh(
  peaks      = grn_set$atac_score,
  genes      = unique(c(grn_set$rna$HGNC, grn_set$rna$ensembl_gene_id)),
  flank_bp   = 30000,
  mode       = "TSS",
  gene_annot = gene_annot_ref_hg38,
  id_col     = "HGNC" # or "ensembl_gene_id"
)

options(future.globals.maxSize = 64 * 1024^3)

# For HiC links
for (gh_std_name in names(gh_std_hic_list)) {
  gh_std <- gh_std_hic_list[[gh_std_name]]

  print(gh_std_name)
  atac_res <- correlate_atac_to_genes(
    grn_set  = grn_set,
    gh_tbl   = gh_std,
    gene_mode = "both",
    fdr       = threshold_atac_gene_corr_p,
    r_abs_min = threshold_atac_gene_corr_abs_r,
    cache_dir = "/data/homes/yl814/episcope_cache/atac_gene_corr",
    cache_tag = "nutrient_hic",
    workers   = 30,
    cache_verbose = TRUE
  )
  fp_res_full <- correlate_fp_to_genes(
    grn_set             = grn_set,
    atac_gene_corr_kept = atac_res$atac_gene_corr_full,
    fdr                 = threshold_fp_gene_corr_p,
    r_abs_min           = threshold_fp_gene_corr_abs_r,
    method              = "pearson",
    workers             = 20,
    cache_dir           = "/data/homes/yl814/episcope_cache/fp_gene_corr",
    cache_tag           = "nutrient_hic",
    cache_chunk_size    = 5000L,
    cache_verbose       = TRUE
  )
  readr::write_csv(atac_res$atac_gene_corr_full, file.path(base_dir, sprintf("atac_gene_corr_full_%s.csv", gh_std_name)))
  readr::write_csv(fp_res_full$fp_gene_corr_full, file.path(base_dir, sprintf("fp_gene_corr_full_%s.csv", gh_std_name)))
}


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


# Run atac gene correlations (adjusted p-value < 0.05, |r| >= 0.3)
threshold_atac_gene_corr_p <- 0.05
# atac_res <- correlate_atac_to_genes(grn_set = grn_set, gh_tbl = gh_std, gene_mode = "both", fdr = threshold_atac_gene_corr_p, r_abs_min = threshold_atac_gene_corr_abs_r)
# atac_res <- correlate_atac_to_genes(
#   grn_set  = grn_set,
#   gh_tbl   = gh_std,
#   gene_mode = "both",
#   fdr       = threshold_atac_gene_corr_p,
#   r_abs_min = threshold_atac_gene_corr_abs_r,
#   cache_dir = "/data/homes/yl814/episcope_cache/atac_gene_corr",
#   cache_tag = "nutrient_100kb",
#   workers   = 30,
#   cache_verbose = TRUE
# )
atac_res <- correlate_atac_to_genes(
  grn_set  = grn_set,
  gh_tbl   = gh_std,
  gene_mode = "both",
  fdr       = threshold_atac_gene_corr_p,
  r_abs_min = threshold_atac_gene_corr_abs_r,
  cor_method = "spearman", # c("pearson", "spearman")
  cache_dir = file.path(base_dir, "cache", "atac_gene_corr"),
  cache_tag = "nutrient_genehancer_spearman",
  workers   = 30,
  cache_verbose = TRUE
)

# atac_res <- correlate_atac_to_genes(grn_set = grn_set, gh_tbl = gh_std, gene_mode = "both", fdr = threshold_atac_gene_corr_p, r_abs_min = threshold_atac_gene_corr_abs_r, cor_method ="spearman")

# shrink to only correlated peaks
# grn_set_filtered <- filter_grn_by_corr(grn_set, atac_res$atac_gene_corr_kept)
fp_res_full_pearson <- correlate_fp_to_genes(
  grn_set             = grn_set,
  atac_gene_corr_kept = atac_res$atac_gene_corr_full,
  fdr                 = threshold_fp_gene_corr_p,
  r_abs_min           = threshold_fp_gene_corr_abs_r,
  method              = "pearson",   # c("pearson", "spearman", "kendall")
  workers             = 20,
  cache_dir           = file.path(base_dir, "cache", "fp_gene_corr"),
  cache_tag           = "nutrient_genehancer",
  cache_chunk_size    = 5000L,
  cache_verbose       = TRUE
)

fp_res_full_spearman <- correlate_fp_to_genes(
  grn_set             = grn_set,
  atac_gene_corr_kept = atac_res$atac_gene_corr_full,
  fdr                 = threshold_fp_gene_corr_p,
  r_abs_min           = threshold_fp_gene_corr_abs_r,
  method              = "spearman",   # c("pearson", "spearman", "kendall")
  workers             = 20,
  cache_dir           = file.path(base_dir, "cache", "fp_gene_corr"),
  cache_tag           = "nutrient_genehancer_spearman",
  cache_chunk_size    = 5000L,
  cache_verbose       = TRUE
)

fp_res_full_kendall <- correlate_fp_to_genes(
  grn_set             = grn_set,
  atac_gene_corr_kept = atac_res$atac_gene_corr_full,
  fdr                 = threshold_fp_gene_corr_p,
  r_abs_min           = threshold_fp_gene_corr_abs_r,
  method              = "kendall",   # c("pearson", "spearman", "kendall")
  workers             = 20,
  cache_dir           = file.path(base_dir, "cache", "fp_gene_corr"),
  cache_tag           = "nutrient_genehancer_kendall",
  cache_chunk_size    = 5000L,
  cache_verbose       = TRUE
)


# Minimal helper: recompute corr_fp_tf_* for fp_annotation
make_fp_annotation_corr <- function(grn_set,
                                    method = c("spearman", "kendall"),
                                    cores = 10L,
                                    chunk_size = 5000L) {
  method <- match.arg(method)

  ann <- grn_set$fp_annotation
  fp  <- grn_set$fp_score
  rna <- grn_set$rna

  sample_cols <- intersect(
    setdiff(names(fp), "peak_ID"),
    setdiff(names(rna), c("ensembl_gene_id", "HGNC"))
  )
  if (length(sample_cols) < 3L) stop("Not enough shared samples between fp_score and rna.")

  fp_m <- as.matrix(fp[, sample_cols, drop = FALSE])
  rownames(fp_m) <- fp$peak_ID

  rna_m <- rna |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::select(HGNC, dplyr::all_of(sample_cols)) |>
    dplyr::distinct(HGNC, .keep_all = TRUE) |>
    as.data.frame()
  rownames(rna_m) <- rna_m$HGNC
  rna_m <- as.matrix(rna_m[, sample_cols, drop = FALSE])

  fp_id <- match(ann$fp_peak, rownames(fp_m))
  tf_id <- match(ann$tfs, rownames(rna_m))

  pairs_u <- unique(data.frame(fp_id = fp_id, tf_id = tf_id))
  pairs_u <- pairs_u[!is.na(pairs_u$fp_id) & !is.na(pairs_u$tf_id), , drop = FALSE]
  if (!nrow(pairs_u)) {
    out <- ann
    out$corr_fp_tf_r <- NA_real_
    out$corr_fp_tf_p <- NA_real_
    out$corr_fp_tf_p_adj <- NA_real_
    return(out)
  }

  idx <- split(seq_len(nrow(pairs_u)),
               ceiling(seq_len(nrow(pairs_u)) / chunk_size))

  worker <- function(ii) {
    sub <- pairs_u[ii, , drop = FALSE]
    r <- numeric(nrow(sub))
    p <- numeric(nrow(sub))

    for (i in seq_len(nrow(sub))) {
      x <- fp_m[sub$fp_id[i], ]
      y <- rna_m[sub$tf_id[i], ]
      ok <- is.finite(x) & is.finite(y)
      n  <- sum(ok)
      if (n < 3L) {
        r[i] <- NA_real_
        p[i] <- NA_real_
        next
      }
      ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
      r[i] <- as.numeric(unname(ct$estimate))
      p[i] <- as.numeric(ct$p.value)
    }
    data.frame(fp_id = sub$fp_id, tf_id = sub$tf_id, r = r, p = p)
  }

  use_mc <- (.Platform$OS.type != "windows") && cores > 1L
  res_list <- if (use_mc) parallel::mclapply(idx, worker, mc.cores = cores) else lapply(idx, worker)
  corr_tbl <- do.call(rbind, res_list)

  key_u <- paste(corr_tbl$fp_id, corr_tbl$tf_id, sep = "_")
  key_a <- paste(fp_id, tf_id, sep = "_")
  m <- match(key_a, key_u)

  out <- ann
  out$corr_fp_tf_r <- corr_tbl$r[m]
  out$corr_fp_tf_p <- corr_tbl$p[m]
  out$corr_fp_tf_p_adj <- stats::p.adjust(out$corr_fp_tf_p, method = "BH")
  out
}

 .row_stats_mean_var <- function(m) {
  m <- as.matrix(m)
  ok <- is.finite(m)
  m[!ok] <- 0
  n <- rowSums(ok)
  sum_x <- rowSums(m)
  sum_x2 <- rowSums(m^2)
  mean_x <- sum_x / n
  var_x <- (sum_x2 - (sum_x^2) / n) / pmax(1, n - 1)
  var_x[n < 2L] <- NA_real_
  list(mean = mean_x, var = var_x)
}

.row_stats_mean_var_chunked <- function(m, cores = 1L, chunk_size = 50000L) {
  m <- as.matrix(m)
  if (cores <= 1L || nrow(m) <= chunk_size) return(.row_stats_mean_var(m))

  idx <- split(seq_len(nrow(m)), ceiling(seq_len(nrow(m)) / chunk_size))
  worker <- function(ii) .row_stats_mean_var(m[ii, , drop = FALSE])
  use_mc <- (.Platform$OS.type != "windows") && cores > 1L
  res_list <- if (use_mc) parallel::mclapply(idx, worker, mc.cores = cores) else lapply(idx, worker)

  mean_x <- numeric(nrow(m))
  var_x <- numeric(nrow(m))
  for (k in seq_along(idx)) {
    i <- idx[[k]]
    mean_x[i] <- res_list[[k]]$mean
    var_x[i] <- res_list[[k]]$var
  }
  list(mean = mean_x, var = var_x)
}

compute_hv_variance_tbl <- function(
  tbl,
  id_cols,
  sample_cols = NULL,
  cores = 1L,
  chunk_size = 50000L
) {
  stopifnot(is.data.frame(tbl))
  if (is.null(sample_cols)) {
    drop_cols <- unique(c(id_cols, "HGNC", "ensembl_gene_id"))
    sample_cols <- setdiff(names(tbl), drop_cols)
  }
  if (!length(sample_cols)) stop("No sample columns found for variance computation.")

  m <- as.matrix(tbl[, sample_cols, drop = FALSE])
  stats_raw <- .row_stats_mean_var_chunked(m, cores = cores, chunk_size = chunk_size)

  out <- tibble::tibble(
    !!!rlang::set_names(lapply(id_cols, function(nm) tbl[[nm]]), id_cols),
    mean_raw = stats_raw$mean,
    var_raw = stats_raw$var
  )

  rsd <- sqrt(stats_raw$var) / stats_raw$mean
  rsd[!is.finite(rsd)] <- NA_real_
  out$rsd <- rsd

  out
}

precompute_hvf_hvg_variance <- function(
  grn_set,
  cores = 1L,
  chunk_size = 50000L
) {
  if (!is.list(grn_set) || !all(c("fp_score", "rna") %in% names(grn_set))) {
    stop("grn_set must contain $fp_score and $rna")
  }
  list(
    fp_variance = compute_hv_variance_tbl(
      tbl = grn_set$fp_score,
      id_cols = "peak_ID",
      cores = cores,
      chunk_size = chunk_size
    ),
    rna_variance = compute_hv_variance_tbl(
      tbl = grn_set$rna,
      id_cols = c("ensembl_gene_id", "HGNC"),
      cores = cores,
      chunk_size = chunk_size
    )
  )
}

# Run and assign back to grn_set
grn_set$fp_annotation_spearman <- make_fp_annotation_corr(grn_set, method = "spearman", cores = 25L)
grn_set$fp_annotation_kendall <- make_fp_annotation_corr(grn_set, method = "kendall",  cores = 25L)
grn_set$fp_annotation_pearson <- grn_set$fp_annotation

hv_variance <- precompute_hvf_hvg_variance(grn_set, cores = 20, chunk_size = 50000L)
grn_set$fp_variance <- hv_variance$fp_variance
grn_set$rna_variance <- hv_variance$rna_variance



# Basal network comparison across methods and FP thresholds.
fp_r_vals <- c(0.3, 0.5)
fp_p_vals <- c(0.1, 0.05, 0.01, 0.001)
basal_r_thr <- 0.3
basal_p_thr <- 0.05

fp_corr_full_by_method <- list(
  pearson  = fp_res_full_pearson$fp_gene_corr_full,
  spearman = fp_res_full_spearman$fp_gene_corr_full,
  kendall  = fp_res_full_kendall$fp_gene_corr_full
)
fp_annot_by_method <- list(
  pearson  = grn_set$fp_annotation_pearson,
  spearman = grn_set$fp_annotation_spearman,
  kendall  = grn_set$fp_annotation_kendall
)

for (method in names(fp_corr_full_by_method)) {
  fp_gene_corr_full <- fp_corr_full_by_method[[method]]
  fp_annotation_use <- fp_annot_by_method[[method]]
  if (is.null(fp_gene_corr_full) || is.null(fp_annotation_use)) next

  out_dir <- file.path(base_dir, "lighting_method_comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (fp_r_thr in fp_r_vals) {
    for (fp_p_thr in fp_p_vals) {
      fp_gene_corr_kept <- fp_gene_corr_full |>
        dplyr::filter(!is.na(.data$r_fp), !is.na(.data$p_adj_fp)) |>
        dplyr::filter(.data$r_fp >= fp_r_thr, .data$p_adj_fp < fp_p_thr)
      # dplyr::filter(abs(.data$r_fp) >= fp_r_thr, .data$p_adj_fp < fp_p_thr)

      tmp_dir <- file.path(
        tempdir(),
        sprintf("basal_tmp_method_%s_fp_r%.1f_p%.2g", method, fp_r_thr, fp_p_thr)
      )
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

      basal <- make_basal_links(
        fp_gene_corr_kept = fp_gene_corr_kept,
        fp_annotation     = fp_annotation_use,
        out_dir           = tmp_dir,
        prefix            = "lighting",
        rna_tbl           = grn_set$rna,
        rna_method        = "pearson",
        rna_cores         = max(1L, parallel::detectCores(logical = TRUE) - 2L)
      )

      basal <- basal |> dplyr::filter(.data$p_adj_tf < basal_p_thr, .data$r_tf > basal_r_thr)
      out_file <- file.path(
        out_dir,
        sprintf("basal_filtered_method_%s_fp_r%.1f_p%.2g.csv", method, fp_r_thr, fp_p_thr)
      )
      readr::write_csv(basal, out_file)
    }
  }
}

# RNA-only method comparison (RNA Pearson) using TF->gene RNA correlations.
out_dir <- file.path(base_dir, "lighting_method_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fp_gene_corr_full <- fp_res_full_pearson$fp_gene_corr_full
fp_annotation_use <- grn_set$fp_annotation_pearson
if (!is.null(fp_gene_corr_full) && !is.null(fp_annotation_use)) {
  for (fp_r_thr in fp_r_vals) {
    for (fp_p_thr in fp_p_vals) {
      tmp_dir <- file.path(
        tempdir(),
        sprintf("basal_tmp_method_rna_pearson_r%.1f_p%.2g", fp_r_thr, fp_p_thr)
      )
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

      basal <- make_basal_links(
        fp_gene_corr_kept = fp_gene_corr_full,
        fp_annotation     = fp_annotation_use,
        out_dir           = tmp_dir,
        prefix            = "lighting",
        rna_tbl           = grn_set$rna,
        rna_method        = "pearson",
        rna_cores         = 20
      )

      basal <- basal |>
        dplyr::filter(.data$p_adj_tf < basal_p_thr, .data$r_tf > basal_r_thr) |>
        dplyr::filter(!is.na(.data$r_rna_gene), !is.na(.data$p_rna_adj_gene)) |>
        dplyr::filter(abs(.data$r_rna_gene) >= fp_r_thr, .data$p_rna_adj_gene < fp_p_thr)

      out_file <- file.path(
        out_dir,
        sprintf("basal_filtered_method_rna_pearson_rna_r%.1f_p%.2g.csv", fp_r_thr, fp_p_thr)
      )
      readr::write_csv(basal, out_file)
    }
  }
}

.detect_link_cols <- function(tbl) {
  tf_col <- intersect(c("TF", "tfs", "tf"), names(tbl))[1]
  gene_col <- intersect(c("gene_key", "gene", "target", "Name.Target"), names(tbl))[1]
  peak_col <- intersect(c("peak_ID", "fp_peak", "peak_id", "peak"), names(tbl))[1]
  if (anyNA(c(tf_col, gene_col, peak_col))) {
    stop("Link columns not found. Have: ", paste(names(tbl), collapse = ", "))
  }
  list(tf = tf_col, gene = gene_col, peak = peak_col)
}

.standardize_links <- function(tbl) {
  cols <- .detect_link_cols(tbl)
  tibble::tibble(
    TF = tbl[[cols$tf]],
    gene_key = tbl[[cols$gene]],
    peak_ID = tbl[[cols$peak]]
  ) |>
    dplyr::distinct()
}

# Union links across methods (r=0.3, p=0.1) and backfill missing method stats.
union_fp_r <- 0.3
union_fp_p <- 0.1
union_files <- list(
  pearson = file.path(base_dir, "lighting_method_comparison", sprintf("basal_filtered_method_pearson_fp_r%.1f_p%.2g.csv", union_fp_r, union_fp_p)),
  spearman = file.path(base_dir, "lighting_method_comparison", sprintf("basal_filtered_method_spearman_fp_r%.1f_p%.2g.csv", union_fp_r, union_fp_p)),
  rna_pearson = file.path(base_dir, "lighting_method_comparison", sprintf("basal_filtered_method_rna_pearson_rna_r%.1f_p%.2g.csv", union_fp_r, union_fp_p))
)

if (all(vapply(union_files, file.exists, logical(1)))) {
  union_tbls <- lapply(union_files, function(f) readr::read_csv(f, show_col_types = FALSE))
  union_links <- dplyr::bind_rows(lapply(union_tbls, .standardize_links)) |>
    dplyr::distinct()

  basal_full_by_method <- list(
    pearson = make_basal_links(
      fp_gene_corr_kept = fp_res_full_pearson$fp_gene_corr_full,
      fp_annotation     = grn_set$fp_annotation_pearson,
      out_dir           = file.path(tempdir(), "basal_full_pearson"),
      prefix            = "lighting",
      rna_tbl           = grn_set$rna,
      rna_method        = "pearson",
      rna_cores         = max(1L, parallel::detectCores(logical = TRUE) - 2L)
    ),
    spearman = make_basal_links(
      fp_gene_corr_kept = fp_res_full_spearman$fp_gene_corr_full,
      fp_annotation     = grn_set$fp_annotation_spearman,
      out_dir           = file.path(tempdir(), "basal_full_spearman"),
      prefix            = "lighting",
      rna_tbl           = grn_set$rna,
      rna_method        = "pearson",
      rna_cores         = max(1L, parallel::detectCores(logical = TRUE) - 2L)
    ),
    rna_pearson = make_basal_links(
      fp_gene_corr_kept = fp_res_full_pearson$fp_gene_corr_full,
      fp_annotation     = grn_set$fp_annotation_pearson,
      out_dir           = file.path(tempdir(), "basal_full_rna_pearson"),
      prefix            = "lighting",
      rna_tbl           = grn_set$rna,
      rna_method        = "pearson",
      rna_cores         = 20
    )
  )

  for (method_key in names(union_files)) {
    full_tbl <- basal_full_by_method[[method_key]]
    cols <- .detect_link_cols(full_tbl)
    full_tbl <- full_tbl |>
      dplyr::distinct(.data[[cols$tf]], .data[[cols$gene]], .data[[cols$peak]], .keep_all = TRUE)
    out_tbl <- union_links |>
      dplyr::left_join(full_tbl, by = c(TF = cols$tf, gene_key = cols$gene, peak_ID = cols$peak))

    drop_cols <- setdiff(c(cols$tf, cols$gene, cols$peak), c("TF", "gene_key", "peak_ID"))
    if (length(drop_cols)) out_tbl <- out_tbl |> dplyr::select(-dplyr::all_of(drop_cols))

    out_tbl <- out_tbl |>
      dplyr::left_join(
        grn_set$fp_variance |>
          dplyr::select(peak_ID, rsd) |>
          dplyr::distinct(peak_ID, .keep_all = TRUE),
        by = "peak_ID"
      )

    out_file <- sub("\\.csv$", "_union.csv", union_files[[method_key]])
    readr::write_csv(out_tbl, out_file)
  }
} else {
  message("Union backfill skipped: missing one or more method files.")
}

print(grn_set$sample_metadata_used,n=23)
#cy414 10_FBS          AsPC1
#cy333 10_FBS          HPAFII
#cy423 10_FBS          Panc1

# fp_gene_corr_kept_rna <- add_tf_rna_corr_cols(
#   fp_tbl  = fp_res_full$fp_gene_corr_kept,
#   rna_tbl = grn_set$rna,
#   method  = "spearman", # c("pearson", "spearman")
#   cores   = max(1L, parallel::detectCores(logical = TRUE) - 10L)
# )
# fp_gene_corr_kept_rna <- fp_gene_corr_kept_rna |>
#   dplyr::left_join(
#     atac_res$atac_gene_corr_full,
#     by = c("atac_peak" = "atac_peak", "gene_key" = "gene_key")
#   )


# Cosine-based "shape match" between (fp_peak vs gene) and (TF vs gene) across samples.
# Adds: cos_fp_tf_shape, cos_fp_tf_shape_p, cos_fp_tf_shape_n
.cosine_shape_stat_one <- function(fp_vec, tf_vec, gene_vec) {
  ok <- is.finite(fp_vec) & is.finite(tf_vec) & is.finite(gene_vec)
  fp_vec <- fp_vec[ok]
  tf_vec <- tf_vec[ok]
  gene_vec <- gene_vec[ok]
  n <- length(gene_vec)
  if (n < 4L) return(list(cos = NA_real_, p = NA_real_, n = n))

  y_fp <- gene_vec[order(fp_vec)]
  y_tf <- gene_vec[order(tf_vec)]

  y_fp_sd <- stats::sd(y_fp)
  y_tf_sd <- stats::sd(y_tf)
  if (!is.finite(y_fp_sd) || !is.finite(y_tf_sd) || y_fp_sd == 0 || y_tf_sd == 0) {
    return(list(cos = NA_real_, p = NA_real_, n = n))
  }
  z_fp <- (y_fp - mean(y_fp)) / y_fp_sd
  z_tf <- (y_tf - mean(y_tf)) / y_tf_sd

  cos <- sum(z_fp * z_tf) / sqrt(sum(z_fp^2) * sum(z_tf^2))
  if (!is.finite(cos)) return(list(cos = NA_real_, p = NA_real_, n = n))
  cos <- max(-1, min(1, cos))

  denom <- max(1e-12, 1 - cos^2)
  t <- cos * sqrt((n - 2) / denom)
  p <- 2 * stats::pt(-abs(t), df = n - 2)

  list(cos = cos, p = p, n = n)
}

add_fp_tf_shape_cosine_cols <- function(
  tbl,
  grn_set,
  sample_ids = NULL,
  fp_peak_col = "fp_peak",
  tf_col = "tfs",
  gene_col = "gene_key"
) {
  stopifnot(is.data.frame(tbl))
  if (all(c("cos_fp_tf_shape", "cos_fp_tf_shape_p", "cos_fp_tf_shape_n") %in% names(tbl))) return(tbl)
  if (!all(c(fp_peak_col, tf_col, gene_col) %in% names(tbl))) {
    stop("tbl must contain columns: ", paste(c(fp_peak_col, tf_col, gene_col), collapse = ", "))
  }
  if (!is.list(grn_set) || !all(c("fp_score", "rna") %in% names(grn_set))) {
    stop("grn_set must contain $fp_score and $rna")
  }

  fp_tbl <- grn_set$fp_score
  rna_tbl <- grn_set$rna
  if (!"peak_ID" %in% names(fp_tbl)) stop("grn_set$fp_score must contain column 'peak_ID'")
  if (!"HGNC" %in% names(rna_tbl)) stop("grn_set$rna must contain column 'HGNC'")
  has_ens <- "ensembl_gene_id" %in% names(rna_tbl)

  if (is.null(sample_ids)) {
    sample_ids <- intersect(setdiff(names(fp_tbl), "peak_ID"), setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC")))
  }
  if (length(sample_ids) < 4L) stop("Could not infer sample_ids (need >= 4). Provide sample_ids explicitly.")

  fp_peaks_needed <- unique(as.character(tbl[[fp_peak_col]]))
  genes_needed <- unique(as.character(tbl[[gene_col]]))
  tfs_needed <- unique(as.character(tbl[[tf_col]]))
  rna_needed <- unique(c(genes_needed, tfs_needed))

  fp_sub <- fp_tbl |>
    dplyr::filter(.data$peak_ID %in% fp_peaks_needed) |>
    dplyr::select(dplyr::any_of(c("peak_ID", sample_ids))) |>
    dplyr::distinct(.data$peak_ID, .keep_all = TRUE)

  rna_sub <- if (isTRUE(has_ens)) {
    rna_tbl |>
      dplyr::filter(.data$HGNC %in% rna_needed | .data$ensembl_gene_id %in% rna_needed) |>
      dplyr::select(dplyr::any_of(c("ensembl_gene_id", "HGNC", sample_ids)))
  } else {
    rna_tbl |>
      dplyr::filter(.data$HGNC %in% rna_needed) |>
      dplyr::select(dplyr::any_of(c("HGNC", sample_ids)))
  }

  fp_mat <- as.matrix(fp_sub[, sample_ids, drop = FALSE])
  rownames(fp_mat) <- fp_sub$peak_ID

  rna_mat_hgnc <- as.matrix(rna_sub[, sample_ids, drop = FALSE])
  rownames(rna_mat_hgnc) <- rna_sub$HGNC
  if (anyDuplicated(rownames(rna_mat_hgnc))) {
    keep <- !duplicated(rownames(rna_mat_hgnc))
    rna_mat_hgnc <- rna_mat_hgnc[keep, , drop = FALSE]
  }

  rna_mat_ens <- NULL
  if (isTRUE(has_ens) && "ensembl_gene_id" %in% names(rna_sub)) {
    ens <- as.character(rna_sub$ensembl_gene_id)
    ok_ens <- !is.na(ens) & ens != ""
    if (any(ok_ens)) {
      rna_mat_ens <- as.matrix(rna_sub[ok_ens, sample_ids, drop = FALSE])
      rownames(rna_mat_ens) <- ens[ok_ens]
      if (anyDuplicated(rownames(rna_mat_ens))) {
        keep <- !duplicated(rownames(rna_mat_ens))
        rna_mat_ens <- rna_mat_ens[keep, , drop = FALSE]
      }
    }
  }

  fp_key <- as.character(tbl[[fp_peak_col]])
  tf_key <- as.character(tbl[[tf_col]])
  gene_key <- as.character(tbl[[gene_col]])

  cos_out <- rep(NA_real_, length(fp_key))
  p_out <- rep(NA_real_, length(fp_key))
  n_out <- rep(NA_integer_, length(fp_key))

  fp_idx <- match(fp_key, rownames(fp_mat))
  tf_idx_h <- match(tf_key, rownames(rna_mat_hgnc))
  gene_idx_h <- match(gene_key, rownames(rna_mat_hgnc))
  tf_idx_e <- if (!is.null(rna_mat_ens)) match(tf_key, rownames(rna_mat_ens)) else rep(NA_integer_, length(tf_key))
  gene_idx_e <- if (!is.null(rna_mat_ens)) match(gene_key, rownames(rna_mat_ens)) else rep(NA_integer_, length(gene_key))

  for (i in seq_along(fp_key)) {
    fi <- fp_idx[[i]]
    if (is.na(fi)) next

    ti <- tf_idx_h[[i]]
    tf_from_hgnc <- TRUE
    if (is.na(ti)) {
      ti <- tf_idx_e[[i]]
      tf_from_hgnc <- FALSE
    }
    if (is.na(ti)) next

    gi <- gene_idx_h[[i]]
    gene_from_hgnc <- TRUE
    if (is.na(gi)) {
      gi <- gene_idx_e[[i]]
      gene_from_hgnc <- FALSE
    }
    if (is.na(gi)) next

    fp_vec <- fp_mat[fi, ]
    tf_vec <- if (isTRUE(tf_from_hgnc)) rna_mat_hgnc[ti, ] else rna_mat_ens[ti, ]
    gene_vec <- if (isTRUE(gene_from_hgnc)) rna_mat_hgnc[gi, ] else rna_mat_ens[gi, ]

    s <- .cosine_shape_stat_one(fp_vec = fp_vec, tf_vec = tf_vec, gene_vec = gene_vec)
    cos_out[[i]] <- s$cos
    p_out[[i]] <- s$p
    n_out[[i]] <- s$n
  }

  tbl$cos_fp_tf_shape <- cos_out
  tbl$cos_fp_tf_shape_p <- p_out
  tbl$cos_fp_tf_shape_n <- n_out
  tbl
}


# simple base model  ------------------------------------------------------
threshold_fp_gene_corr_p <- 0.1

fp_r_vals <- c(0.3, 0.5, 0.7)
fp_p_vals <- c(0.1, 0.05, 0.01)
basal_r_vals <- c(0.3, 0.5, 0.7)
basal_p_vals <- c(0.1, 0.05, 0.01)

fp_r_fixed <- 0.3
fp_p_fixed <- 0.1
basal_r_fixed <- 0.3
basal_p_fixed <- 0.05
dry_run <- FALSE
workers_light <- 6L
workers_bulk <- 6L
cores_rna <- max(1L, min(6L, parallel::detectCores(logical = TRUE) - 2L))

fp_corr_full_list <- list(
  spearman = fp_res_full_spearman$fp_gene_corr_full,
  kendall = fp_res_full_kendall$fp_gene_corr_full
)

fp_annotation_list <- list(
  spearman = grn_set$fp_annotation_spearman,
  kendall = grn_set$fp_annotation_kendall
)

fp_corr_full_list <- Filter(Negate(is.null), fp_corr_full_list)
fp_annotation_list <- Filter(Negate(is.null), fp_annotation_list)
if (!length(fp_annotation_list) && !is.null(grn_set$fp_annotation)) {
  fp_annotation_list <- list(default = grn_set$fp_annotation)
}
if (!length(fp_corr_full_list) || !length(fp_annotation_list)) {
  stop("Missing fp_annotation_spearman/kendall (or fp_annotation) or fp_res_full(_kendall) in memory.")
}

run_lighting_one <- function(fp_gene_corr_kept,
                             fp_annotation_use,
                             lighting_folder,
                             basal_p_thr,
                             basal_r_thr,
                             dry_run = FALSE,
                             workers_light = 6L,
                             workers_bulk = 6L) {
  if (isTRUE(dry_run)) {
    message("[dry_run] skip: ", lighting_folder)
    return(invisible(NULL))
  }

  basal <- make_basal_links(
    fp_gene_corr_kept = fp_gene_corr_kept,
    fp_annotation     = fp_annotation_use,
    out_dir           = lighting_folder,
    prefix            = "lighting"
  )

  basal <- basal |> dplyr::filter(p_adj_tf < basal_p_thr, r_tf > basal_r_thr)

  light_by_condition(
    ds = grn_set,
    basal_links = basal,
    out_dir = lighting_folder,
    prefix = "lighting",
    label_col = "strict_match_rna",
    link_score_threshold = link_score_threshold,
    fp_score_threshold = fp_score_threshold,
    tf_expr_threshold = threshold_tf_expr,
    use_parallel = TRUE,
    workers = workers_light
  )

  specs <- build_cellwise_contrasts_from_index(
    index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
    out_dir = lighting_folder,
    prefix = "lighting",
    ctrl_tag = "10_FBS",
    clean_names = FALSE
  )
  str(specs)

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
    workers = workers_bulk
  )

    # bulk$filtered_paths
    # bulk$filtered_dirs
    # bulk$manifest_path

    # LDA across all filtered CSVs produced above
    # filtered_csvs <- bulk$filtered_paths
    # lda_out <- episcope::annotate_links_with_lda_topic_bulk(
    #   filtered_csvs,
    #   K = 20, # topics
    #   which = "abs", # |delta_link_score|
    #   min_df = 2, # drop ultra-rare terms
    #   gamma_cutoff = 0.2, # keep all topics with gamma >= 0.2 per TF
    #   parallel = TRUE, workers = 20, plan = "multisession", seed = 1
    # )

    # Summarize all annotated CSVs produced by LDA
    # annotated_csvs <- lda_out$assigned_paths
    # sum_out <- episcope::summarize_lda_annotations_bulk(
    #   annotated_csvs,
    #   edge_filter_min = delta_link, # keep edges used to call delta links
    #   parallel = TRUE, plan = "multisession", workers = 20,
    #   use_tf_r_weight = TRUE # optional: weight deltas by TF r column if present
    # )

  # Generate all bubble & hub PDFs
  # summary_csvs <- sum_out$summary_paths
  # episcope::plot_from_summary_bulk(
  #   summary_csvs,
  #   top_tf_per_topic = 20,
  #   parallel = TRUE,
  #   workers = 20,
  #   plan = "multisession",
  #   color_sigma = 2
  # )
}

ann_order <- intersect(c("kendall", "spearman"), names(fp_annotation_list))
if (!length(ann_order)) ann_order <- names(fp_annotation_list)
for (ann_name in ann_order) {
  fp_annotation_use <- fp_annotation_list[[ann_name]]
  fp_gene_corr_full <- fp_corr_full_list[[ann_name]]
  if (is.null(fp_annotation_use)) {
    stop("fp_annotation for method '", ann_name, "' is missing.")
  }
  if (is.null(fp_gene_corr_full)) {
    stop("fp_gene_corr_full for method '", ann_name, "' is missing.")
  }
  fp_gene_corr_full_rna <- add_tf_rna_corr_cols(
    fp_tbl  = fp_gene_corr_full,
    rna_tbl = grn_set$rna,
    method  = ann_name,
    cores   = cores_rna
  )

  # match_mode <- "strict"  # strict lenient
  regulated_genes <- 1.5 # 1.5 2
  delta_link <- 1 # 2 1
  regulation_priors <- "genehancer" # "genehancer" 300kb

  # (A) Basal filter grid, fixed FP filter
  fp_gene_corr_fixed <- fp_gene_corr_full_rna |>
    dplyr::filter(r_fp * r_rna > 0) |>
    dplyr::filter(abs(r_fp) > fp_r_fixed) |>
    dplyr::filter(abs(r_rna) > 0.3) |>
    dplyr::filter(p_adj_rna < 0.05) |>
    dplyr::filter(p_adj_fp < fp_p_fixed)

  for (basal_p_thr in basal_p_vals) {
    for (basal_r_thr in basal_r_vals) {
      fp_tag <- sprintf("fp_r%.1f_p%.2g", fp_r_fixed, fp_p_fixed)
      basal_tag <- sprintf("basal_r%.1f_p%.2g", basal_r_thr, basal_p_thr)
      lighting_folder <- file.path(
        base_dir,
        paste(
          "lighting",
          "fp_tf_corr_FDR", threshold_fp_tf_corr_p,
          regulation_priors, db,
          "regulated_genes", regulated_genes,
          "delta_link", delta_link,
          paste0("fpann_", ann_name),
          paste0("method_", ann_name),
          "fpfilter_fixed",
          "basalfilter_grid",
          fp_tag, basal_tag,
          sep = "_"
        )
      )
      print(lighting_folder)
      run_lighting_one(
        fp_gene_corr_fixed,
        fp_annotation_use,
        lighting_folder,
        basal_p_thr,
        basal_r_thr,
        dry_run = dry_run,
        workers_light = workers_light,
        workers_bulk = workers_bulk
      )
    }
  }

  # (B) FP filter grid, fixed basal filter
  for (fp_p_thr in fp_p_vals) {
    for (fp_r_thr in fp_r_vals) {
      fp_gene_corr_grid <- fp_gene_corr_full_rna |>
        dplyr::filter(r_fp * r_rna > 0) |>
        dplyr::filter(abs(r_fp) > fp_r_thr) |>
        dplyr::filter(abs(r_rna) > 0.3) |>
        dplyr::filter(p_adj_rna < 0.05) |>
        dplyr::filter(p_adj_fp < fp_p_thr)

      fp_tag <- sprintf("fp_r%.1f_p%.2g", fp_r_thr, fp_p_thr)
      basal_tag <- sprintf("basal_r%.1f_p%.2g", basal_r_fixed, basal_p_fixed)
      lighting_folder <- file.path(
        base_dir,
        paste(
          "lighting",
          "fp_tf_corr_FDR", threshold_fp_tf_corr_p,
          regulation_priors, db,
          "regulated_genes", regulated_genes,
          "delta_link", delta_link,
          paste0("fpann_", ann_name),
          paste0("method_", ann_name),
          "fpfilter_grid",
          "basalfilter_fixed",
          fp_tag, basal_tag,
          sep = "_"
        )
      )
      print(lighting_folder)
      run_lighting_one(
        fp_gene_corr_grid,
        fp_annotation_use,
        lighting_folder,
        basal_p_fixed,
        basal_r_fixed,
        dry_run = dry_run,
        workers_light = workers_light,
        workers_bulk = workers_bulk
      )
    }
  }
}





# start filtering grid ----------------------------------------------------
fp_gene_corr_joined <- fp_gene_corr_kept_rna



# Flexible filtering across FP / RNA / ATAC correlation columns.
# - FP:  r_fp + (p_fp or p_adj_fp)
# - RNA: r_rna + (p_rna or p_adj_rna)
# - ATAC (optional): r_atac + (p_atac or p_adj_atac)
# - Direction: always enforce FP and RNA same direction; optionally enforce ATAC and RNA same direction.
filter_fp_rna_atac_corr <- function(
  tbl,
  fp_p_thr,
  fp_r_abs_min,
  rna_p_thr,
  rna_r_abs_min,
  fp_p_col = "p_adj_fp",
  rna_p_col = "p_adj_rna",
  require_same_dir_fp_rna = TRUE,
  use_atac = FALSE,
  atac_p_thr = NULL,
  atac_r_abs_min = NULL,
  atac_p_col = "p_adj_atac",
  require_same_dir_atac_rna = FALSE
) {
  stopifnot(is.data.frame(tbl))
  stopifnot(is.numeric(fp_p_thr), length(fp_p_thr) == 1L, is.finite(fp_p_thr))
  stopifnot(is.numeric(fp_r_abs_min), length(fp_r_abs_min) == 1L, is.finite(fp_r_abs_min), fp_r_abs_min >= 0)
  stopifnot(is.numeric(rna_p_thr), length(rna_p_thr) == 1L, is.finite(rna_p_thr))
  stopifnot(is.numeric(rna_r_abs_min), length(rna_r_abs_min) == 1L, is.finite(rna_r_abs_min), rna_r_abs_min >= 0)
  stopifnot(is.character(fp_p_col), length(fp_p_col) == 1L)
  stopifnot(is.character(rna_p_col), length(rna_p_col) == 1L)
  stopifnot(is.logical(require_same_dir_fp_rna), length(require_same_dir_fp_rna) == 1L)
  stopifnot(is.logical(use_atac), length(use_atac) == 1L)
  stopifnot(is.character(atac_p_col), length(atac_p_col) == 1L)
  stopifnot(is.logical(require_same_dir_atac_rna), length(require_same_dir_atac_rna) == 1L)

  if (!fp_p_col %in% names(tbl)) stop(sprintf("Missing column '%s' in tbl", fp_p_col))
  if (!rna_p_col %in% names(tbl)) stop(sprintf("Missing column '%s' in tbl", rna_p_col))
  if (!all(c("r_fp", "r_rna") %in% names(tbl))) stop("Missing one of required columns: r_fp, r_rna")

  out <- tbl |>
    dplyr::filter(!is.na(.data$r_fp), !is.na(.data$r_rna)) |>
    dplyr::filter(!is.na(.data[[fp_p_col]]), !is.na(.data[[rna_p_col]]))

  if (isTRUE(require_same_dir_fp_rna)) {
    out <- out |> dplyr::filter(.data$r_fp * .data$r_rna > 0)
  }

  out <- out |>
    dplyr::filter(abs(.data$r_fp) >= fp_r_abs_min) |>
    dplyr::filter(abs(.data$r_rna) >= rna_r_abs_min) |>
    dplyr::filter(.data[[fp_p_col]] <= fp_p_thr) |>
    dplyr::filter(.data[[rna_p_col]] <= rna_p_thr)

  if (isTRUE(use_atac)) {
    if (!all(c("r_atac", atac_p_col) %in% names(out))) {
      stop(sprintf("ATAC filtering requested but missing columns: r_atac and/or %s", atac_p_col))
    }
    stopifnot(is.numeric(atac_p_thr), length(atac_p_thr) == 1L, is.finite(atac_p_thr))
    stopifnot(is.numeric(atac_r_abs_min), length(atac_r_abs_min) == 1L, is.finite(atac_r_abs_min), atac_r_abs_min >= 0)

    out <- out |>
      dplyr::filter(!is.na(.data$r_atac), !is.na(.data[[atac_p_col]])) |>
      dplyr::filter(abs(.data$r_atac) >= atac_r_abs_min) |>
      dplyr::filter(.data[[atac_p_col]] <= atac_p_thr)

    if (isTRUE(require_same_dir_atac_rna)) {
      out <- out |> dplyr::filter(.data$r_atac * .data$r_rna > 0)
    }
  }

  out
}


# Grid of thresholds to run (in-memory).
# Edit the vectors below (or post-filter `filter_grid`) to add/remove combinations.
fp_p_col_use <- "p_adj_fp"   # or "p_fp"
rna_p_col_use <- "p_adj_rna" # or "p_rna"
atac_p_col_use <- "p_adj_atac" # or "p_atac"
tf_p_col_use <- "p_adj_tf" # (basal TF~FP correlation) typically "p_adj_tf"

# FP/RNA thresholds (used for both FP and RNA filters)
fp_rna_p_vals <- c(0.01, 0.05, 0.10)
fp_rna_r_vals <- c(0.3, 0.5, 0.7)

# Basal TF~FP thresholds
tf_p_vals <- c(0.01, 0.05)
tf_r_vals <- c(0.3, 0.5)

# ATAC thresholds (only used when use_atac = TRUE)
atac_p_val <- 0.05
atac_r_val <- 0.5
atac_same_dir_vals <- c(FALSE, TRUE)

grid_core <- tidyr::crossing(
  fp_rna_p = fp_rna_p_vals,
  fp_rna_r = fp_rna_r_vals,
  tf_p = tf_p_vals,
  tf_r = tf_r_vals
)

filter_grid <- dplyr::bind_rows(
  grid_core |>
    dplyr::mutate(
      use_atac = FALSE,
      atac_p = NA_real_,
      atac_r = NA_real_,
      atac_same_dir = NA
    ),
  grid_core |>
    tidyr::crossing(atac_same_dir = atac_same_dir_vals) |>
    dplyr::mutate(
      use_atac = TRUE,
      atac_p = atac_p_val,
      atac_r = atac_r_val
    )
) |>
  dplyr::arrange(.data$use_atac, .data$fp_rna_p, .data$fp_rna_r, .data$atac_p, .data$atac_r, .data$atac_same_dir, .data$tf_p, .data$tf_r)

message("[filter_grid] prepared ", nrow(filter_grid), " rows in-memory")

# optional quick sanity checks
# sum(!is.na(fp_gene_corr_kept_rna_gam$coef_p_value))
# fp_gene_corr_kept_rna_gam |> dplyr::count(is.na(coef_p_value))




# fp_gene_corr_kept_rna_gam <- fp_gene_corr_kept_rna_gam |>
#   # dplyr::filter(r_atac * r_rna > 0) |>
#   dplyr::filter(abs(r_atac) > 0.3) |>
#   dplyr::filter(p_adj_atac < threshold_fp_gene_corr_p)


# readr::write_csv(atac_res$atac_gene_corr_full, file.path(base_dir, sprintf("atac_gene_corr_full_30kb_%s.csv", db)))
# readr::write_csv(fp_res_full$fp_gene_corr_full, file.path(base_dir, sprintf("fp_gene_corr_full_30kb_%s.csv", db)))

# atac_gene_corr_full_30kb <- readr::read_csv(file.path(base_dir, sprintf("atac_gene_corr_full_30kb_%s.csv", db)))
# fp_gene_corr_full_30kb <- readr::read_csv(file.path(base_dir, sprintf("fp_gene_corr_full_30kb_%s.csv", db)))



# Run fp gene correlations (adjusted p-value < 0.05, |r| >= 0.3)
# threshold_fp_gene_corr_p <- 0.3
# fp_res <- correlate_fp_to_genes(
#   grn_set              = grn_set_filtered,
#   atac_gene_corr_kept  = atac_res$atac_gene_corr_kept,
#   fdr                  = threshold_fp_gene_corr_p,
#   r_abs_min            = threshold_fp_gene_corr_abs_r,
#   workers              = 15
# )



# match_mode <- "strict"  # strict lenient
regulated_genes <- 1.5 # 1.5 2
delta_link <- 1 # 2 1
regulation_priors <- "genehancer" # "genehancer" 300kb
run_lighting_pipeline_one <- function(
  fp_gene_corr_kept,
  lighting_folder,
  grn_set,
  regulated_genes,
  delta_link,
  tf_p_thr,
  tf_r_min,
  tf_p_col = "p_adj_tf"
) {
  options(future.globals.maxSize = max(getOption("future.globals.maxSize", 500 * 1024^2), 32 * 1024^3))
  dir.create(lighting_folder, recursive = TRUE, showWarnings = FALSE)

  # Step 3. Build basal GRN & identify active regulatory edges per condition
  basal <- make_basal_links(
    fp_gene_corr_kept = fp_gene_corr_kept,
    fp_annotation     = grn_set$fp_annotation,
    out_dir           = lighting_folder,
    prefix            = "lighting"
  )

  if (!all(c("r_tf", tf_p_col) %in% names(basal))) {
    stop(sprintf("Basal TF filtering requires columns r_tf and %s", tf_p_col))
  }
  basal <- basal |>
    dplyr::filter(!is.na(.data$r_tf), !is.na(.data[[tf_p_col]])) |>
    dplyr::filter(.data$r_tf >= tf_r_min) |>
    dplyr::filter(.data[[tf_p_col]] <= tf_p_thr)

  light_by_condition(
    ds = grn_set,
    basal_links = basal,
    out_dir = lighting_folder,
    prefix = "lighting",
    label_col = "strict_match_rna",
    link_score_threshold = link_score_threshold,
    fp_score_threshold = fp_score_threshold,
    tf_expr_threshold = threshold_tf_expr,
    use_parallel = TRUE,
    workers = 8
  )

  # Step 4. Perform differential GRN analysis & identify master TFs
  specs <- build_cellwise_contrasts_from_index(
    index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
    out_dir = lighting_folder,
    prefix = "lighting",
    ctrl_tag = "10_FBS",
    clean_names = FALSE
  )

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

  # LDA across all filtered CSVs produced above
  # filtered_csvs <- bulk$filtered_paths

  # lda_out <- episcope::annotate_links_with_lda_topic_bulk(
  #   filtered_csvs,
  #   K = 20, # topics
  #   which = "abs", # |delta_link_score|
  #   min_df = 2, # drop ultra-rare terms
  #   gamma_cutoff = 0.2, # keep all topics with gamma >= 0.2 per TF
  #   parallel = TRUE, workers = 20, plan = "multisession", seed = 1
  # )

  # Summarize all annotated CSVs produced by LDA
  # annotated_csvs <- lda_out$assigned_paths

  # sum_out <- episcope::summarize_lda_annotations_bulk(
  #   annotated_csvs,
  #   edge_filter_min = delta_link, # keep edges used to call delta links
  #   parallel = TRUE, plan = "multisession", workers = 20,
  #   use_tf_r_weight = TRUE # optional: weight deltas by TF r column if present
  # )

  # Generate all bubble & hub PDFs
  # summary_csvs <- sum_out$summary_paths

  # episcope::plot_from_summary_bulk(
  #   summary_csvs,
  #   top_tf_per_topic = 20,
  #   parallel = TRUE,
  #   workers = 20,
  #   plan = "multisession",
  #   color_sigma = 2
  # )

  invisible(list(
    lighting_folder = lighting_folder,
    bulk = bulk
  ))
}

lighting_folders <- character(0)
last_run <- NULL

for (i in seq_len(nrow(filter_grid))) {
  cfg <- filter_grid[i, , drop = FALSE]

  fp_p_thr <- cfg$fp_rna_p
  fp_r_thr <- cfg$fp_rna_r
  rna_p_thr <- cfg$fp_rna_p
  rna_r_thr <- cfg$fp_rna_r

  use_atac <- isTRUE(cfg$use_atac)
  atac_p_thr <- if (use_atac) cfg$atac_p else NULL
  atac_r_thr <- if (use_atac) cfg$atac_r else NULL
  atac_same_dir <- if (use_atac) isTRUE(cfg$atac_same_dir) else FALSE

  tf_p_thr <- cfg$tf_p
  tf_r_thr <- cfg$tf_r

  fp_gene_corr_kept_rna_filtered <- filter_fp_rna_atac_corr(
    tbl = fp_gene_corr_joined,
    fp_p_thr = fp_p_thr,
    fp_r_abs_min = fp_r_thr,
    rna_p_thr = rna_p_thr,
    rna_r_abs_min = rna_r_thr,
    fp_p_col = fp_p_col_use,
    rna_p_col = rna_p_col_use,
    require_same_dir_fp_rna = TRUE,
    use_atac = use_atac,
    atac_p_thr = atac_p_thr,
    atac_r_abs_min = atac_r_thr,
    atac_p_col = atac_p_col_use,
    require_same_dir_atac_rna = atac_same_dir
  )

  # Encode the thresholds used into the output folder name.
  # Example suffix:
  #   FP[p_adj<=0.05|r>=0.5]_RNA[p_adj<=0.05|r>=0.5]_ATAC[none]
  # or
  #   ..._ATAC[p_adj<=0.01|r>=0.7|sameDirFALSE]
  p_label <- function(p_col) {
    if (grepl("adj", p_col, ignore.case = TRUE)) "pAdj" else "p"
  }

  suffix_fp <- sprintf("FP_%s%.2g_r%.1f", p_label(fp_p_col_use), fp_p_thr, fp_r_thr)
  suffix_rna <- sprintf("RNA_%s%.2g_r%.1f", p_label(rna_p_col_use), rna_p_thr, rna_r_thr)
  suffix_tf <- {
    if (!is.finite(tf_p_thr) || !is.finite(tf_r_thr)) stop("TF thresholds must be finite")
    sprintf("TF_%s%.2g_r%.1f", p_label(tf_p_col_use), tf_p_thr, tf_r_thr)
  }
  suffix_atac <- if (!use_atac) {
    "ATAC_none"
  } else {
    # guard against accidental NAs
    if (!is.finite(atac_p_thr) || !is.finite(atac_r_thr)) {
      stop("ATAC thresholds must be finite when use_atac=TRUE")
    }
    sprintf("ATAC_%s%.2g_r%.1f_sameDir%s", p_label(atac_p_col_use), atac_p_thr, atac_r_thr, if (atac_same_dir) "TRUE" else "FALSE")
  }

  lighting_folder <- file.path(
    base_dir,
    paste(
      "lighting",
      "fp_tf_corr_FDR", threshold_fp_tf_corr_p,
      regulation_priors, db,
      "regulated_genes", regulated_genes,
      "delta_link", delta_link,
      suffix_fp, suffix_rna, suffix_tf, suffix_atac,
      sep = "_"
    )
  )
  print(lighting_folder)
  lighting_folders <- c(lighting_folders, lighting_folder)

  last_run <- run_lighting_pipeline_one(
    fp_gene_corr_kept = fp_gene_corr_kept_rna_filtered,
    lighting_folder = lighting_folder,
    grn_set = grn_set,
    regulated_genes = regulated_genes,
    delta_link = delta_link,
    tf_p_thr = tf_p_thr,
    tf_r_min = tf_r_thr,
    tf_p_col = tf_p_col_use
  )
}

# Keep a final lighting_folder value for any downstream interactive code blocks.
if (length(lighting_folders) > 0) {
  lighting_folder <- lighting_folders[[length(lighting_folders)]]
}

# Expose the last run's outputs for downstream interactive analysis blocks.
if (!is.null(last_run)) {
  bulk <- last_run$bulk
  filtered_csvs <- bulk$filtered_paths
  lda_out <- last_run$lda_out
  sum_out <- last_run$sum_out
}

# Topic models ---------------------------------------------------------------
# Moved to dev/08_benchmark_topic_methods.R (source that script to run benchmarks)


# TODO TF HITS score summary plot (waterfall)

# Step 5. Generate interactive Topic & TF regulatory hub subnetwor --------

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
