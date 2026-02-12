suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(matrixStats)
  library(RColorBrewer)
  library(cli)
})

source("R/utils_logging.R")
source("R/utils_ggheat.R")
source("R/utils_step1_footprints.R")
source("R/utils_step1_align_footprints.R")

# -----------------------------------------------------------------------------
# Benchmark footprint heatmaps:
# - raw footprint score (before normalization; sample columns)
# - aligned + quantile-normalized footprint score (condition columns)
# Uses one shared annotation palette but does not force identical clustering/order.
# -----------------------------------------------------------------------------

base_dir <- Sys.getenv(
  "FP_HEATMAP_BASE_DIR",
  unset = "/data/homes/yl814/episcope_test/nutrient_stress_strict_JASPAR2024/predict_tf_binding_sites"
)
db_use <- Sys.getenv("FP_HEATMAP_DB", unset = "JASPAR2024")
project_base_dir <- dirname(base_dir)
cache_dir <- file.path(project_base_dir, "cache")
fp_root_dir <- Sys.getenv("FP_HEATMAP_FP_ROOT_DIR", unset = "")
metadata_csv <- Sys.getenv(
  "FP_HEATMAP_METADATA_CSV",
  unset = file.path(project_base_dir, "data", "sample_metadata_strict.csv")
)
prealign_csv <- Sys.getenv(
  "FP_HEATMAP_PREALIGN_CSV",
  unset = file.path(cache_dir, sprintf("fp_scores_prealign_%s.csv", db_use))
)
manifest_csv <- file.path(cache_dir, sprintf("fp_%s_manifest.csv", db_use))
manifest_out_dir <- file.path(cache_dir, sprintf("fp_%s", db_use))
qn_csv <- file.path(base_dir, sprintf("03_fp_score_qn_%s.csv", db_use))
rds_file <- file.path(base_dir, sprintf("01_multiomic_data_object_%s.rds", db_use))
out_dir <- Sys.getenv(
  "FP_HEATMAP_OUT_DIR",
  unset = file.path(base_dir, "fp_alignment_heatmap")
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) {
  .log_abort("Could not create output directory: {.path {out_dir}}")
}
if (file.access(out_dir, mode = 2) != 0) {
  .log_abort("Output directory is not writable: {.path {out_dir}}")
}

if (!file.exists(qn_csv)) {
  .log_abort("Required input is missing: {.path {qn_csv}}")
}

if (!file.exists(prealign_csv)) {
  .log_inform("Pre-alignment score table not found. Generating: {.path {prealign_csv}}")

  fp_manifest <- NULL
  if (file.exists(manifest_csv)) {
    .log_inform("Loading existing motif manifest: {.path {manifest_csv}}")
    fp_manifest <- readr::read_csv(manifest_csv, show_col_types = FALSE)
  } else {
    if (!nzchar(fp_root_dir)) {
      .log_abort(
        "Missing pre-alignment score file and manifest. Set `FP_HEATMAP_FP_ROOT_DIR` to build motif manifest."
      )
    }
    .log_inform("Building motif manifest from root dir: {.path {fp_root_dir}}")
    fp_manifest <- load_footprints(
      root_dir = fp_root_dir,
      db_name = db_use,
      out_dir = manifest_out_dir,
      skip_existing = TRUE,
      verbose = TRUE
    )
  }

  # Generate pre-alignment score table without running full pipeline.
  suppressWarnings(
    align_footprints(
      fp_manifest,
      cache_dir = cache_dir,
      cache_tag = db_use,
      output_mode = "distinct",
      use_cache = FALSE,
      write_cache = FALSE,
      save_prealign_score = TRUE,
      prealign_score_path = prealign_csv,
      prealign_output_mode = "distinct",
      prealign_only = TRUE,
      verbose = TRUE
    )
  )

  if (!file.exists(prealign_csv)) {
    .log_abort("Failed to create pre-alignment score table: {.path {prealign_csv}}")
  }
}

.log_inform("Loading pre-alignment raw footprint score: {.path {prealign_csv}}")
raw_tbl <- readr::read_csv(prealign_csv, show_col_types = FALSE) |>
  dplyr::distinct(.data$peak_ID, .keep_all = TRUE)

.log_inform("Loading aligned/QN footprint score: {.path {qn_csv}}")
qn_tbl <- readr::read_csv(qn_csv, show_col_types = FALSE) |>
  dplyr::distinct(.data$peak_ID, .keep_all = TRUE)

if (file.exists(rds_file)) {
  .log_inform("Loading metadata from RDS: {.path {rds_file}}")
  omics_data <- readRDS(rds_file)
  if (!"sample_metadata_used" %in% names(omics_data)) {
    .log_abort("`sample_metadata_used` is not present in the RDS object.")
  }
  meta <- omics_data$sample_metadata_used
} else if (file.exists(metadata_csv)) {
  .log_inform("Loading metadata from CSV: {.path {metadata_csv}}")
  meta <- readr::read_csv(metadata_csv, show_col_types = FALSE)
} else {
  .log_abort("No metadata source found: {.path {rds_file}} or {.path {metadata_csv}}")
}
need_meta_cols <- c("id", "strict_match_rna")
if (!all(need_meta_cols %in% names(meta))) {
  miss <- setdiff(need_meta_cols, names(meta))
  .log_abort("Metadata is missing required columns: {paste(miss, collapse = ', ')}")
}

# Build sample-id -> condition-label map for raw columns.
id_to_label <- meta |>
  dplyr::filter(!is.na(.data$strict_match_rna), nzchar(.data$strict_match_rna)) |>
  dplyr::distinct(.data$id, .data$strict_match_rna)

raw_ids <- setdiff(names(raw_tbl), "peak_ID")
map_tbl <- tibble::tibble(id = raw_ids) |>
  dplyr::left_join(id_to_label, by = "id")

if (anyNA(map_tbl$strict_match_rna)) {
  n_miss <- sum(is.na(map_tbl$strict_match_rna))
  miss_ids <- map_tbl$id[is.na(map_tbl$strict_match_rna)]
  .log_warn(
    "{n_miss} raw sample column(s) have no strict_match_rna label and will be dropped: {paste(miss_ids, collapse = ', ')}"
  )
}

map_tbl <- map_tbl |>
  dplyr::filter(!is.na(.data$strict_match_rna), nzchar(.data$strict_match_rna))

if (!nrow(map_tbl)) {
  .log_abort("No raw sample columns could be mapped to strict_match_rna labels.")
}

# Keep raw at sample level (before alignment/normalization view).
raw_ids_keep <- map_tbl$id
raw_mat <- as.matrix(raw_tbl[, raw_ids_keep, drop = FALSE])
storage.mode(raw_mat) <- "numeric"
rownames(raw_mat) <- raw_tbl$peak_ID

qn_cols <- setdiff(names(qn_tbl), "peak_ID")
qn_mat <- as.matrix(qn_tbl[, qn_cols, drop = FALSE])
storage.mode(qn_mat) <- "numeric"
rownames(qn_mat) <- qn_tbl$peak_ID

# Shared rows between raw and qn.
shared_peaks <- intersect(rownames(raw_mat), rownames(qn_mat))
if (!length(shared_peaks)) {
  .log_abort("No shared peak_ID rows between raw and QN matrices.")
}

raw_aligned <- raw_mat[shared_peaks, , drop = FALSE]
qn_aligned <- qn_mat[shared_peaks, , drop = FALSE]

# Select high-variance peaks by QN matrix.
hv_sizes <- 20000L
vars <- matrixStats::rowVars(qn_aligned, na.rm = TRUE)
ord <- order(vars, decreasing = TRUE)

parse_condition_annotation <- function(labels) {
  tibble::tibble(label = labels) |>
    dplyr::mutate(
      CellLine = stringr::str_replace(.data$label, "_.*$", ""),
      StressRaw = stringr::str_replace(.data$label, "^.*?_", ""),
      Stress = dplyr::if_else(
        .data$StressRaw == "10_FBS",
        .data$StressRaw,
        .data$StressRaw |>
          stringr::str_remove("_\\d+$") |>
          stringr::str_remove("^([0-9]+(?:\\.[0-9]+)?(?:uM)?)_")
      )
    ) |>
    dplyr::select("label", "CellLine", "Stress")
}

raw_label_lookup <- stats::setNames(map_tbl$strict_match_rna, map_tbl$id)

raw_col_ann <- parse_condition_annotation(unname(raw_label_lookup[colnames(raw_aligned)])) |>
  dplyr::mutate(sample_id = colnames(raw_aligned), .before = 1) |>
  dplyr::select("sample_id", "CellLine", "Stress") |>
  tibble::column_to_rownames("sample_id")

qn_col_ann <- parse_condition_annotation(colnames(qn_aligned)) |>
  tibble::column_to_rownames("label")

all_col_ann <- dplyr::bind_rows(
  raw_col_ann |> tibble::rownames_to_column("id"),
  qn_col_ann |> tibble::rownames_to_column("id")
)

all_col_ann <- all_col_ann |>
  dplyr::mutate(
    CellLine = as.character(.data$CellLine),
    Stress = as.character(.data$Stress)
  )

cell_levels <- unique(all_col_ann$CellLine)
stress_levels <- unique(all_col_ann$Stress)

cell_pal <- c("#6893c6", "#e74b5b", "#eca72e", "#43aa8b", "#f94144", "#577590")
stress_pal <- c("grey", "#3a78af", "#ef812f", "#308a4e", "#414b8c", "#e89c84", "#ca2f2d", "#916ab6", "#f4a261", "#2a9d8f")

annotation_colors <- list(
  CellLine = stats::setNames(cell_pal[seq_along(cell_levels)], cell_levels),
  Stress = stats::setNames(stress_pal[seq_along(stress_levels)], stress_levels)
)

heat_cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100)

write_ordered_xlsx <- function(mat, out_xlsx) {
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    .log_warn("Package `openxlsx2` not installed; skipping XLSX export: {.path {out_xlsx}}")
    return(invisible(NULL))
  }
  out_df <- as.data.frame(mat)
  out_df <- cbind(peak_ID = rownames(mat), out_df)

  wb <- openxlsx2::wb_workbook()$
    add_worksheet("Heatmap")$
    add_data("Heatmap", out_df, rowNames = FALSE)$
    freeze_pane("Heatmap", firstRow = TRUE, firstCol = TRUE)

  wb$save(out_xlsx)
}

n_use <- min(length(ord), hv_sizes)
keep_idx <- ord[seq_len(n_use)]

qn_top <- qn_aligned[keep_idx, , drop = FALSE]
raw_top <- raw_aligned[rownames(qn_top), , drop = FALSE]

qn_pdf <- file.path(out_dir, sprintf("fp_after_alignment_qn_hv%d_heatmap.pdf", n_use))
raw_pdf <- file.path(out_dir, sprintf("fp_raw_before_alignment_hv%d_heatmap.pdf", n_use))
qn_xlsx <- file.path(out_dir, sprintf("fp_after_alignment_qn_hv%d_heatmap.xlsx", n_use))
raw_xlsx <- file.path(out_dir, sprintf("fp_raw_before_alignment_hv%d_heatmap.xlsx", n_use))

combined_vals <- c(as.vector(raw_top), as.vector(qn_top))
combined_vals <- combined_vals[is.finite(combined_vals)]
if (!length(combined_vals)) {
  .log_abort("No finite values available for heatmap color scale.")
}
rng <- range(combined_vals, na.rm = TRUE)
if (rng[1] == rng[2]) {
  rng <- rng + c(-1e-6, 1e-6)
}
my_breaks <- seq(rng[1], rng[2], length.out = length(heat_cols) + 1L)

.log_inform("Plotting QN heatmap (top {n_use} HV peaks): {.path {qn_pdf}}")
hm_qn <- ggheat(
  qn_top,
  annotation_col = qn_col_ann[colnames(qn_top), , drop = FALSE],
  annotation_colors = annotation_colors,
  scale = "none",
  color = heat_cols,
  breaks = my_breaks,
  clustering_method = "ward.D2",
  distfun = function(x) dist(x, method = "euclidean"),
  border_color = NA,
  show_rownames = FALSE,
  filename = qn_pdf,
  silent = TRUE,
  width = 10,
  height = 18
)

if (is.null(hm_qn$matrix_ordered)) {
  .log_abort("ggheat did not return `matrix_ordered`; please use updated R/utils_ggheat.R.")
}
qn_ordered <- as.matrix(hm_qn$matrix_ordered)

.log_inform("Plotting raw heatmap before alignment/QN (top {n_use} HV peaks): {.path {raw_pdf}}")
hm_raw <- ggheat(
  raw_top,
  annotation_col = raw_col_ann[colnames(raw_top), , drop = FALSE],
  annotation_colors = annotation_colors,
  scale = "none",
  color = heat_cols,
  breaks = my_breaks,
  clustering_method = "ward.D2",
  distfun = function(x) dist(x, method = "euclidean"),
  border_color = NA,
  show_rownames = FALSE,
  filename = raw_pdf,
  silent = TRUE,
  width = 10,
  height = 18
)

if (is.null(hm_raw$matrix_ordered)) {
  .log_abort("ggheat did not return `matrix_ordered` for raw heatmap.")
}
raw_ordered <- as.matrix(hm_raw$matrix_ordered)

write_ordered_xlsx(qn_ordered, qn_xlsx)
write_ordered_xlsx(raw_ordered, raw_xlsx)

.log_inform("Saved outputs for hv={n_use}: raw-before + after-alignment/QN heatmaps and ordered matrices.")

.log_inform("Done. Output directory: {.path {out_dir}}")
