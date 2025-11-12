# render_tf_hub_delta_network(
#   comp_csv  = "Z:/episcope/inst/extdata/lighting/AsPC1_Glc_vs_AsPC1_Ctrl_delta_links_filtered_lda_K20.csv",
#   input_tf  = "HNF4A",
#   out_html  = "AsPC1_Glc_vs_AsPC1_Ctrl_delta_links_HNF4A.html",
#   edge_filter_min = 2,
#   edge_filter_on  = "either",
#   gene_fc_thresh  = 1.5,
#   de_reference    = "str_over_ctrl",
#   motif_db = "jaspar2024",
#   ring_tf_only_direct = FALSE
# )
# --- CONFIG ---
# src_dir  <- "Y:/cy232/cutntag/humanPDAC/Nutrients_basal/Episcope/select_network"
# out_dir  <- "Z:/episcope_test_nutrient_stress/select_network_plots"

src_dir <- out_dir <- lighting_folder
layout   <- "fr"           # fixed single layout
show_peaks <- FALSE        # keep triangles hidden by default
set.seed(1L)
motif_db <- "jaspar2024"

# --- helpers (no library(); explicit namespacing only) ---
is_ctrl_tag <- function(tag) {
  # heuristics for picking the control column
  any(grepl(c("(^|_)ctrl($|_)", "control", "baseline", "basal", "_10_fbs($|_)", "10fbs"),
            tolower(tag), perl = TRUE))
}
detect_mapping <- function(df_names) {
  # choose the two link_score_* columns
  sc <- df_names[grepl("^link_score_", df_names)]
  if (length(sc) < 2L) stop("Could not find two 'link_score_*' columns.")
  # prefer exactly two; if more, keep the first two deterministically
  sc <- sc[seq_len(min(2L, length(sc)))]
  suf <- sub("^link_score_", "", sc)
  # assign ctrl/str by heuristics on the suffix
  ctrl_idx <- if (is_ctrl_tag(suf[1]) && !is_ctrl_tag(suf[2])) 1L else
    if (!is_ctrl_tag(suf[1]) && is_ctrl_tag(suf[2])) 2L else 2L # fall back: 2nd is ctrl
  str_idx  <- if (ctrl_idx == 1L) 2L else 1L
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
# pick_motif_db <- function(base_name) {
#   if (grepl("hocomoco", base_name, ignore.case = TRUE)) "hocomocov13" else "jaspar2024"
# }

# --- run ---
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(src_dir, pattern = "_K20.csv$", full.names = TRUE)
if (!length(files)) stop("No .txt or .csv files found in: ", src_dir)

for (in_path in files) {
  base_name <- sub("\\.txt$", "", basename(in_path))
  message(">>> Processing: ", base_name)

  # read
  tf_network <- readr::read_csv(in_path, show_col_types = FALSE)
  ns <- names(tf_network)

  # per-file mappings and motif_db
  mp <- detect_mapping(ns)
  # motif_db <- pick_motif_db(base_name)

  # plot
  w <- try(
    plot_tf_network_delta(
      data = tf_network,
      plot_title = base_name,
      layout_algo = layout,
      physics     = TRUE,          # uses stabilization+freeze inside the function
      add_direct  = TRUE,
      edge_filter_min = 0,
      min_delta_abs   = 0,
      keep_top_edges_per_tf = 6000,
      peak_mode   = "show_all",
      show_peaks  = show_peaks,
      gene_fc_thresh  = 1.5,
      de_reference    = "str_over_ctrl",
      motif_db        = motif_db,
      # auto mappings
      score_ctrl_col      = mp$score_ctrl_col,
      score_str_col       = mp$score_str_col,
      sign_ctrl_col       = if (mp$sign_ctrl_col %in% ns) mp$sign_ctrl_col else NULL,
      sign_str_col        = if (mp$sign_str_col  %in% ns) mp$sign_str_col  else NULL,
      tf_expr_ctrl_col    = if (mp$tf_expr_ctrl_col   %in% ns) mp$tf_expr_ctrl_col   else NULL,
      tf_expr_str_col     = if (mp$tf_expr_str_col    %in% ns) mp$tf_expr_str_col    else NULL,
      gene_expr_ctrl_col  = if (mp$gene_expr_ctrl_col %in% ns) mp$gene_expr_ctrl_col else NULL,
      gene_expr_str_col   = if (mp$gene_expr_str_col  %in% ns) mp$gene_expr_str_col  else NULL
    ),
    silent = TRUE
  )
  if (inherits(w, "try-error")) {
    message("    x Failed: ", conditionMessage(attr(w, "condition")))
    next
  }

  out_file <- file.path(out_dir, paste0(base_name, "_layout_", layout, ".html"))
  htmlwidgets::saveWidget(w, out_file, selfcontained = FALSE)
  message("    ✓ Saved: ", normalizePath(out_file, mustWork = FALSE))
}









