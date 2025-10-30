# utils_grn_link_network_tf_delta.R — TF-centric Δ subnetwork (single panel)
# Author: Yaoxiang Li (cleaned & packaged)
# Updated: 2025-10-30
#
# NOTE:
# - Public API: render_tf_hub_delta_network()
# - All other objects are internal (helpers, constants).
# - Behavior and visuals are preserved to match your original script exactly.
# - No library() calls; explicit namespacing used throughout.

#' TF-hub Δ subnetwork (single panel)
#'
#' Build and save an interactive, TF-centric **delta-links** network from a single
#' “comparison filtered clustered” CSV where per-condition measures already live side-by-side
#' (columns like `*_ctrl` and `*_str`). The visualization encodes **Δ = stress − control**
#' on edges and preserves your legacy styling and HTML structure (node shapes, colors,
#' widths, legends, and jump menus).
#'
#' Conceptually, the plot places the **root TF** at the center, arranges a ring of TFs around it,
#' and shows how regulatory links fan out to target genes through shared or unique enhancer
#' “peak” pseudonodes. TF box fill reflects signed TF RNA log2FC, gene dot size reflects the
#' maximum absolute Z-score of expression across sides, and edge hue/width convey the sign and
#' magnitude of Δ link strength. Dashed gray segments correspond to TF→peak, thin gray bridges
#' connect peak→gene, and direct TF→gene Δ edges are drawn faintly when also covered by peaks.
#'
#' @param comp_csv Path to the comparison CSV (e.g., `*_delta_links_filtered_lda_K20.csv`).
#' @param input_tf Root TF symbol to visualize at the center.
#' @param out_html Optional output HTML path (default: alongside \code{comp_csv},
#'   named `\{file\}_tfhub_delta_\{TF\}.html`).
#' @param edge_filter_min Numeric minimum absolute `link_score_*` required for an edge to be
#'   considered (filtering uses score magnitude only; “active” flags still affect styling).
#' @param edge_filter_on One of \code{"either"}, \code{"stress"}, \code{"control"}, \code{"both"},
#'   indicating on which side(s) the magnitude threshold must pass.
#' @param gene_fc_thresh Fold-change threshold used for up/down border coloring (default `1.5`).
#' @param de_reference One of \code{"str_over_ctrl"} (default) or \code{"ctrl_over_str"}; controls
#'   the sign convention for TF box fill and gene border log2FC.
#' @param ring_tf_outgoing \code{"context"} (legacy default), \code{"full"}, or \code{"full_topk"}.
#'   Controls how **ring TFs** (non-root TFs) contribute outgoing edges:
#'   \itemize{
#'     \item \strong{context}: only show edges needed to explain co-regulated genes with the root TF,
#'           plus TF↔TF context.
#'     \item \strong{full}: include all outgoing edges from ring TFs within the TF set.
#'     \item \strong{full_topk}: like \emph{full} but per-TF keep only the top-K strongest edges.
#'   }
#' @param ring_tf_topk Integer; only used if \code{ring_tf_outgoing="full_topk"} (default `250`).
#' @param ring_tf_direct_only Logical; if \code{TRUE}, the ring of TFs includes \emph{only} TFs that are
#'   directly connected to the root TF (either the root regulates them or they regulate the root).
#'   No additional TFs are pulled in merely because they co-regulate the same gene(s). When
#'   \code{FALSE} (default), behavior matches the current implementation (direct TFs plus those
#'   second-order TFs brought in by co-regulation context).
#' @param motif_db Character; \code{"hocomocov13"} (default) or \code{"jaspar2024"}—used to establish
#'   the TF symbol universe for TF vs. gene typing.
#' @param verbose Logical; print progress messages.
#'
#' @details
#' Expected columns (some are auto-detected with flexible tags): \cr
#' \code{TF}, \code{gene_key}, \code{peak_id}, \code{link_score_ctrl}, \code{link_score_str},
#' \code{link_sign_ctrl}, \code{link_sign_str}, optional activity flags
#' (\code{active_ctrl}, \code{active_str}), TF expression
#' (\code{tf_expr_ctrl}, \code{tf_expr_str}), and gene expression
#' (\code{gene_expr_ctrl}, \code{gene_expr_str}). Column names are resolved using the inferred
#' control/stress tags from the filename header (e.g., `A_vs_B`) with several common aliases.
#'
#' @return (Invisibly) the path to the written HTML file.
#'
#' @section Visual encodings:
#' \itemize{
#'   \item \strong{TF nodes}: box fill = signed TF RNA log2FC (blue→white→red);
#'         label font size ∝ max(TF RNA).
#'   \item \strong{Gene nodes}: dot size ∝ max absolute Z(expr) across sides; border color reflects
#'         up/down per \code{gene_fc_thresh} under \code{de_reference}.
#'   \item \strong{Edges}: color = red (Δ>0) or blue (Δ<0); width = fixed bins by |Δ|;
#'         TF→peak dashed gray; peak→gene neutral gray bridges; direct TF→gene Δ drawn faint.
#' }
#'
#' @seealso \code{render_link_network_delta_topic()}, \code{render_link_network_for_topic()}
#'
#' @author Yaoxiang Li and Chunling Yi
#' @export
#'
#' @examples
#' \dontrun{
#' render_tf_hub_delta_network(
#'   comp_csv  = system.file(
#'     "extdata/lighting_strict_hocomocov13_regulated_genes_1.5_delta_link_1",
#'     "AsPC1_Glc_vs_AsPC1_Ctrl_delta_links_filtered_lda_K20.csv",
#'     package = "episcope"
#'   ),
#'   input_tf  = "HSF2",
#'   edge_filter_min = 1,
#'   edge_filter_on  = "either",
#'   gene_fc_thresh  = 1.5,
#'   de_reference    = "str_over_ctrl",
#'   ring_tf_direct_only = TRUE,   # only direct TFs form the ring
#'   motif_db        = "hocomocov13",
#'   verbose         = TRUE
#' )
#' }
render_tf_hub_delta_network <- function(
    comp_csv,
    input_tf,
    out_html = NULL,
    edge_filter_min = 2,
    edge_filter_on = c("either","stress","control","both"),
    gene_fc_thresh = 1.5,
    de_reference = c("str_over_ctrl","ctrl_over_str"),
    ring_tf_outgoing = "full",   # "context", "full", "full_topk"
    ring_tf_topk = 250,
    motif_db = c("hocomocov13","jaspar2024"),
    ring_tf_only_direct = FALSE,   # <-- NEW: only direct TFs around the ring

    verbose = TRUE
){
  edge_filter_on <- match.arg(edge_filter_on)
  de_reference   <- match.arg(de_reference)
  motif_db       <- match.arg(motif_db)

  # Ensure TF symbol universe is loaded for this DB (cached once per session)
  .init_known_tf_symbols(motif_db = motif_db)

  stopifnot(file.exists(comp_csv))
  .llog("Reading comparison CSV: ", comp_csv, verbose = verbose)
  C <- readr::read_csv(comp_csv, show_col_types = FALSE)

  base <- basename(comp_csv)
  hdr  <- sub("_delta_links.*$", "", base)
  bits <- strsplit(hdr, "_vs_", fixed = TRUE)[[1]]
  stress_tag <- bits[[1]] %||% "stress"
  ctrl_tag   <- bits[[2]] %||% "control"

  border_up_for <- if (de_reference == "ctrl_over_str") "control" else "stress"

  core <- .build_nodes_edges_from_root_tf(
    comp_df = C,
    input_tf = input_tf,
    edge_filter_min = edge_filter_min,
    edge_filter_on  = edge_filter_on,
    gene_fc_thresh  = gene_fc_thresh,
    de_reference    = de_reference,
    cond1_tag = ctrl_tag,
    cond2_tag = stress_tag,
    border_up_for = border_up_for,
    ring_tf_outgoing = ring_tf_outgoing,
    ring_tf_topk = ring_tf_topk,
    ring_tf_only_direct = ring_tf_only_direct,  # <-- NEW: pass through

    verbose = verbose
  )

  nodes  <- core$nodes
  coords <- .layout_root_tf_circles(
    nodes,
    child_tfs         = core$direct_tfs,
    common_sig        = core$common_sig,
    common_nonsig     = core$common_nonsig,
    uniq_sig          = core$uniq_sig,
    uniq_nonsig       = core$uniq_nonsig,
    ed_sub            = core$ed_sub,
    root_uniq_sig     = core$root_uniq_sig,
    root_uniq_nonsig  = core$root_uniq_nonsig,
    layout_mode       = "ring",
    kk_maxiter        = 350
  )

  built <- .tfhub_build_delta_panel(core, coords, stress_tag = stress_tag, ctrl_tag = ctrl_tag)
  nd  <- built$nodes
  edv <- built$edges

  main_hdr <- sprintf(
    "%s | TF hub Δ-links: %s  (Δ = %s − %s; TF→PEAK; bridges)",
    hdr, input_tf, stress_tag, ctrl_tag
  )

  widget <- visNetwork::visNetwork(nd, edv, width = "100%", height = "700px", main = main_hdr) |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                           nodesIdSelection = TRUE) |>
    visNetwork::visEdges(scaling = list(min = 1, max = 10), smooth = FALSE) |>
    visNetwork::visNodes(borderWidth = 1, font = list(bold = TRUE),
                         margin = list(top=6,right=10,bottom=6,left=10)) |>
    visNetwork::visPhysics(enabled = FALSE) |>
    visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) |>
    .attach_vis_jumpers() |>
    .tf_on_top()

  page <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title(main_hdr),
      htmltools::tags$style(htmltools::HTML("
        body { font-family: Helvetica, Arial, sans-serif; margin: 16px; }
        .grn-legend {
          position: relative !important;
          float: none !important;
          clear: both;
          width: min(1200px, 96vw);
          margin: 12px auto 6px;
          box-sizing: border-box;
          display: flex;
          gap: 18px;
          align-items: center;
          justify-content: center;
          flex-wrap: wrap;
          font-size: 12px;
          color: #555;
          white-space: normal !important;
          word-break: normal !important;
          overflow: visible;
        }
        .grn-legend .grad  { display:inline-block; width:40px; height:14px; background:linear-gradient(90deg,#4575b4,#FFFFFF,#d73027); border:1px solid #BFBFBF; margin-right:6px;}
        .grn-legend .dot   { display:inline-block; width:12px; height:12px; background:#D9D9D9; border-radius:50%; border:1px solid #BFBFBF; margin-right:6px;}
        .grn-legend .plus  { display:inline-block; width:24px; height:3px; background:#d73027; margin-right:6px;}
        .grn-legend .minus { display:inline-block; width:24px; height:3px; background:#4575b4; margin-right:6px;}
      "))
    ),
    htmltools::tags$body(
      htmltools::h3(main_hdr),
      widget,
      htmltools::div(
        class="grn-legend",
        htmltools::HTML(paste0(
          "<span>Legend:</span>",
          "<span><span class='grad'></span>TF (box) fill by signed log2FC(TF RNA): blue (−) → white (0) → red (+); TF label size ∝ max(TF RNA)</span>",
          "<span><span class='dot'></span>Gene (dot) size ∝ max |Z(expr)| across sides</span>",
          "<span>Triangle = peak pseudonode (shared enhancer)</span>",
          "<span><em>Dashed</em> edges = TF→peak; bridges (peak→gene) are thin solid gray; direct Δ TF→gene are solid but slightly faded</span>",
          "<span>Edge color: red if Δ&gt;0, blue if Δ&lt;0 (Δ = ", stress_tag, " − ", ctrl_tag, ")</span>",
          "<span>Edge width binned by |Δ link_score|; arrow tip: pointed for link_sign = + (activation), blunt bar for link_sign = − (repression)</span>",
          "<span>Reserved pie slice shows genes uniquely regulated by the center TF</span>"
        ))
      )
    )
  )

  if (is.null(out_html)) {
    file_base <- tools::file_path_sans_ext(basename(comp_csv))
    out_html <- file.path(
      dirname(comp_csv),
      sprintf("%s_tfhub_delta_%s.html", file_base, .safe(input_tf))
    )
  }
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  htmltools::save_html(
    htmltools::as.tags(page),
    file   = out_html,
    libdir = paste0(tools::file_path_sans_ext(basename(out_html)), "_files")
  )
  .llog("✓ wrote: ", normalizePath(out_html, mustWork = FALSE), verbose = verbose)
  invisible(out_html)
}

# ─────────────────────────────────────────────────────────────────────────────
# INTERNALS (helpers, constants) — do not export
# ─────────────────────────────────────────────────────────────────────────────

# Colors / knobs
.col_edge_pos   <- "#D86B5E"
.col_edge_neg   <- "#6F9ED4"
.col_edge_zero  <- "#C8C8C8"
.col_gene_fill  <- "#D9D9D9"
.col_border_none <- "#BFBFBF"
.col_border_up   <- "#D86B5E"
.col_border_down <- "#6F9ED4"

.zero_edge_alpha  <- 0.6
.zero_edge_width  <- 0.5
.zero_edge_dashes <- FALSE
.edge_alpha_range <- c(0.10, 0.85)
.edge_width_range <- c(0.6, 5)
.edge_width_gamma <- 0.75

.tf_font_range <- c(22, 44)
.tf_font_gamma <- 1
.tf_label_halo_w <- 0.25
.tf_label_halo_c <- "#BFBFBF"

.gene_size_range <- c(1, 5)
.gene_size_gamma <- 0.9

.layout_gene_min_sep_px <- 30
.peak_node_px <- 6

`%||%` <- function(a, b) if (is.null(a)) b else a
.safe <- function(x) gsub("[^A-Za-z0-9_.-]+","_", x)
.llog <- function(..., verbose = TRUE){ if (isTRUE(verbose)) message("[tfhubΔ] ", paste0(..., collapse = "")) }

robust_z <- function(x){
  x <- as.numeric(x)
  m <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(madv) || madv == 0) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
    return((x - m)/(sx + 1e-12))
  }
  (x - m)/(madv + 1e-12)
}
clamp <- function(x, lo, hi) pmax(pmin(x, hi), lo)
.to_rgba <- function(hex, alpha){
  hex <- gsub("^#", "", hex)
  r <- strtoi(substr(hex,1,2),16L); g <- strtoi(substr(hex,3,4),16L); b <- strtoi(substr(hex,5,6),16L)
  sprintf("rgba(%d,%d,%d,%.3f)", r,g,b, pmax(pmin(alpha,1),0))
}


# Known TF universe (cached per session, package-safe)
.tf_cache <- new.env(parent = emptyenv())

.get_tf_syms <- function() {
  get0("KNOWN_TF_SYMBOLS", envir = .tf_cache, inherits = FALSE)
}

.set_tf_syms <- function(x, motif_db) {
  assign("KNOWN_TF_SYMBOLS", x, envir = .tf_cache)
  assign("KNOWN_TF_DB", motif_db, envir = .tf_cache)
  invisible(TRUE)
}

.get_motif_cluster_tsv <- function(motif_db = c("hocomocov13","jaspar2024")){
  motif_db <- match.arg(motif_db)
  if (motif_db == "jaspar2024") {
    p <- system.file("extdata","genome","JASPAR2024.txt", package = "episcope")
  } else {
    p <- system.file("extdata","genome","HOCOMOCO_v13_motif_cluster_definition_with_sub_cluster.txt",
                     package = "episcope")
  }
  if (!nzchar(p) || !file.exists(p)) {
    stop("Motif cluster definition TSV not found for motif_db='", motif_db,
         "'. Ensure the file is installed in inst/extdata/genome/.")
  }
  p
}

.load_known_tf_symbols <- function(path){
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (!"HGNC" %in% names(df)) stop("HGNC column not found in: ", path)
  syms <- df$HGNC
  syms <- unlist(strsplit(syms, "::", fixed = TRUE), use.names = FALSE)
  syms <- trimws(as.character(syms))
  syms <- syms[nzchar(syms)]
  unique(syms)
}

.init_known_tf_symbols <- function(motif_db = c("hocomocov13","jaspar2024")){
  motif_db <- match.arg(motif_db)
  cur_db   <- get0("KNOWN_TF_DB", envir = .tf_cache, inherits = FALSE)
  cur_set  <- .get_tf_syms()
  if (is.null(cur_set) || is.null(cur_db) || !identical(cur_db, motif_db)) {
    path <- .get_motif_cluster_tsv(motif_db)
    .set_tf_syms(.load_known_tf_symbols(path), motif_db)
  }
  invisible(TRUE)
}

.is_known_tf <- function(x){
  x <- trimws(as.character(x))
  if (is.null(.get_tf_syms())) .init_known_tf_symbols("hocomocov13")
  x %in% (.get_tf_syms() %||% character(0))
}



# Column picking / coercions
.pick_col <- function(df, candidates, label, required = TRUE){
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    if (isTRUE(required)) stop("Missing required column for ", label,
                               " (one of: ", paste(candidates, collapse=", "), ").")
    return(NULL)
  }
  hit[[1]]
}
.detect_cols <- function(df, cond1_tag = NULL, cond2_tag = NULL){
  .tag_variants <- function(x){
    if (is.null(x) || !nzchar(x)) return(character(0))
    x0 <- as.character(x)
    x1 <- gsub("[^A-Za-z0-9_.-]+", "_", x0)
    x2 <- gsub("[^A-Za-z0-9]+", "_", x0)
    unique(c(x0, x1, x2))
  }
  v1 <- .tag_variants(cond1_tag)
  v2 <- .tag_variants(cond2_tag)
  add_dyn <- function(fixed, bases, tags){
    dyn <- as.vector(outer(bases, tags, function(b, t) paste0(b, "_", t)))
    unique(c(fixed, dyn))
  }
  list(
    tf     = .pick_col(df, c("TF","tf"), "TF"),
    gene   = .pick_col(df, c("gene_key","gene","target"), "gene_key"),
    peak   = .pick_col(df, c("peak_id","peak_ID","peak","peak_key","peakID"), "peak_id"),

    score_ctrl = .pick_col(df, add_dyn(c("edge_score_ctrl","link_score_ctrl","score_ctrl"),
                                       c("edge_score","link_score","score"), v1), "edge_score_ctrl"),
    score_str  = .pick_col(df, add_dyn(c("edge_score_str","link_score_str","score_str"),
                                       c("edge_score","link_score","score"), v2), "edge_score_str"),
    sign_ctrl  = .pick_col(df, add_dyn(c("link_sign_ctrl","edge_sign_ctrl","sign_ctrl"),
                                       c("link_sign","edge_sign","sign"), v1), "link_sign_ctrl"),
    sign_str   = .pick_col(df, add_dyn(c("link_sign_str","edge_sign_str","sign_str"),
                                       c("link_sign","edge_sign","sign"), v2), "link_sign_str"),

    active_ctrl = .pick_col(df, add_dyn(c("active_ctrl","is_active_ctrl"),
                                        c("active","is_active"), v1), "active_ctrl", required = FALSE),
    active_str  = .pick_col(df, add_dyn(c("active_str","is_active_str"),
                                        c("active","is_active"), v2), "active_str", required = FALSE),

    r_tf_ctrl = .pick_col(df, add_dyn(c("r_tf_ctrl","r_tf_cond2","r_tf_A","r_tf_1"),
                                      c("r_tf"), v1), "r_tf_ctrl", required = FALSE),
    r_tf_str  = .pick_col(df, add_dyn(c("r_tf_str","r_tf_cond1","r_tf_B","r_tf_2"),
                                      c("r_tf"), v2), "r_tf_str", required = FALSE),

    tf_expr_ctrl = .pick_col(df, add_dyn(c("tf_expr_ctrl","ctrl_tf_expr","tf_expr_cond2","tf_expr_A","tf_expr_1"),
                                         c("tf_expr"), v1), "tf_expr_ctrl", required = FALSE),
    tf_expr_str  = .pick_col(df, add_dyn(c("tf_expr_str","str_tf_expr","tf_expr_cond1","tf_expr_B","tf_expr_2"),
                                         c("tf_expr"), v2), "tf_expr_str", required = FALSE),
    gene_expr_ctrl = .pick_col(df, add_dyn(c("gene_expr_ctrl","expr_ctrl","ctrl_gene_expr","gene_expr_1","gene_expr_A"),
                                           c("gene_expr","expr"), v1), "gene_expr_ctrl", required = FALSE),
    gene_expr_str  = .pick_col(df, add_dyn(c("gene_expr_str","expr_str","str_gene_expr","gene_expr_2","gene_expr_B"),
                                           c("gene_expr","expr"), v2), "gene_expr_str", required = FALSE)
  )
}
.as_sign <- function(v){
  if (is.numeric(v)) {
    s <- ifelse(v > 0, 1L, ifelse(v < 0, -1L, NA_integer_))
    s[!is.finite(v)] <- NA_integer_
    return(s)
  }
  vv <- tolower(trimws(as.character(v)))
  dplyr::case_when(
    vv %in% c("+","pos","positive","act","activation","1","up","plus") ~  1L,
    vv %in% c("-","neg","negative","rep","repression","-1","down","minus") ~ -1L,
    TRUE ~ NA_integer_
  )
}
.as_bool <- function(v){
  if (is.logical(v)) return(v)
  if (is.numeric(v)) return(v != 0)
  vv <- tolower(trimws(as.character(v)))
  ifelse(vv %in% c("1","true","t","yes","y"), TRUE,
         ifelse(vv %in% c("0","false","f","no","n",""), FALSE, NA))
}

.compute_pass_flags <- function(ed, min_abs, col_score, col_active = NULL){
  sc <- suppressWarnings(as.numeric(ed[[col_score]]))
  ok_mag <- is.finite(sc) & (abs(sc) >= min_abs)
  if (!is.null(col_active) && col_active %in% names(ed)) {
    act <- .as_bool(ed[[col_active]]); act[is.na(act)] <- FALSE
  } else act <- rep(TRUE, nrow(ed))
  ok_mag & act
}
.compute_caps <- function(v_pos, v_neg, q = 0.995){
  pos <- stats::quantile(v_pos[v_pos > 0], q, na.rm = TRUE, names = FALSE)
  neg <- stats::quantile(v_neg[v_neg < 0], 1 - q, na.rm = TRUE, names = FALSE)
  list(
    pos = ifelse(is.finite(pos), pos, suppressWarnings(max(v_pos[v_pos > 0], na.rm = TRUE))),
    neg = ifelse(is.finite(neg), neg, suppressWarnings(min(v_neg[v_neg < 0], na.rm = TRUE)))
  )
}

# Aggregate per-peak→gene (TF-independent)
.aggregate_peak_gene <- function(ed){
  ed |>
    dplyr::group_by(peak_id, gene_key) |>
    dplyr::summarise(
      score_ctrl = dplyr::first(stats::na.omit(score_ctrl)),
      score_str  = dplyr::first(stats::na.omit(score_str)),
      sign_ctrl  = dplyr::first(stats::na.omit(sign_ctrl)),
      sign_str   = dplyr::first(stats::na.omit(sign_str)),
      pass_ctrl  = any(pass_ctrl, na.rm = TRUE),
      pass_str   = any(pass_str,  na.rm = TRUE),
      .groups = "drop"
    )
}

# Styling helpers
.fixed_delta_width_bin <- function(mag_abs){
  dplyr::case_when(
    is.finite(mag_abs) & mag_abs >  8 ~ 4,
    is.finite(mag_abs) & mag_abs >  4 ~ 3,
    is.finite(mag_abs) & mag_abs >  2 ~ 2,
    is.finite(mag_abs) & mag_abs >  0 ~ 1,
    TRUE ~ 1
  )
}
.fan_multi_edges <- function(edges){
  if (!nrow(edges)) return(edges)
  ee <- edges |>
    dplyr::group_by(from, to) |>
    dplyr::mutate(.k = dplyr::row_number(), .n = dplyr::n()) |>
    dplyr::ungroup()
  ee[, c("id","from","to","width","color","dashes","title")]
}
.force_edge_cols <- function(df){
  if (!nrow(df)) {
    return(tibble::tibble(
      id=character(0), from=character(0), to=character(0),
      width=numeric(0), color=character(0),
      dashes=logical(0), title=character(0),
      arrows=I(list())
    ))
  }
  if (!"dashes" %in% names(df)) df$dashes <- FALSE
  if (!"arrows" %in% names(df)) df$arrows <- ""

  if (is.list(df$dashes)) {
    df$dashes <- vapply(df$dashes, function(x) TRUE, logical(1))
  } else {
    df$dashes <- as.logical(df$dashes); df$dashes[is.na(df$dashes)] <- FALSE
  }

  normalize_one_arrow <- function(val){
    if (is.character(val) && length(val) == 1) {
      vv <- trimws(val)
      if (identical(vv, "to"))     return(list(to   = list(enabled = TRUE,  type = "arrow")))
      if (identical(vv, "from"))   return(list(from = list(enabled = TRUE,  type = "arrow")))
      if (identical(vv, "to;bar")) return(list(to   = list(enabled = TRUE,  type = "bar")))
      return(list(to = list(enabled = FALSE)))
    }
    if (is.list(val)) return(val)
    list(to = list(enabled = FALSE))
  }
  if (!is.list(df$arrows)) {
    df$arrows <- I(lapply(df$arrows, normalize_one_arrow))
  } else {
    if (length(df$arrows) == 1L && nrow(df) > 1L) {
      df$arrows <- I(replicate(nrow(df), df$arrows[[1]], simplify = FALSE))
    } else {
      df$arrows <- I(lapply(df$arrows, normalize_one_arrow))
    }
  }
  df[, c("id","from","to","width","color","dashes","title","arrows")]
}

.style_edges_one_side <- function(ed, side = c("ctrl","str"),
                                  caps, alpha_range = .edge_alpha_range,
                                  width_range = .edge_width_range, width_gamma = .edge_width_gamma,
                                  zero_col = .col_edge_zero, zero_alpha = .zero_edge_alpha,
                                  zero_width = .zero_edge_width, zero_dashes = .zero_edge_dashes,
                                  log_max_joint = NULL){
  side <- match.arg(side)
  s <- if (side == "ctrl") ed$score_ctrl else ed$score_str
  signv <- if (side == "ctrl") ed$sign_ctrl else ed$sign_str
  pass  <- if (side == "ctrl") ed$pass_ctrl else ed$pass_str

  sval <- signv * abs(s); sval[!is.finite(sval)] <- 0
  cap  <- ifelse(sval >= 0, caps$pos, abs(caps$neg))
  mag_cap <- pmin(abs(s), cap) / pmax(cap, 1e-9)
  alpha <- alpha_range[1] + (alpha_range[2] - alpha_range[1]) * mag_cap
  alpha[!pass] <- zero_alpha

  base_hex <- ifelse(signv >= 0, .col_edge_pos, .col_edge_neg)
  base_hex[!pass] <- zero_col
  col <- mapply(.to_rgba, base_hex, alpha, USE.NAMES = FALSE)

  w <- rep(zero_width, nrow(ed))
  idx <- which(pass & is.finite(s) & abs(s) > 0)
  if (length(idx)) {
    mag_log <- log2(1 + abs(s[idx]))
    max_log <- if (is.null(log_max_joint)) max(mag_log, na.rm = TRUE) else log_max_joint
    if (is.finite(max_log) && max_log > 0) {
      sc <- (mag_log / max_log)^width_gamma
      w[idx] <- width_range[1] + sc * (width_range[2] - width_range[1])
    } else w[idx] <- mean(width_range)
  }

  tibble::tibble(
    id    = paste(ed$TF, ed$gene_key, ed$peak_id, side, sep = "|"),
    from  = ed$TF,
    to    = ed$gene_key,
    width = w,
    color = col,
    dashes = ifelse(pass, FALSE, zero_dashes),
    title = sprintf("%s → %s<br>peak: %s<br>%s: score=%s sign=%s",
                    ed$TF, ed$gene_key, ed$peak_id, side,
                    ifelse(is.finite(s), sprintf("%.3f", s), "NA"),
                    ifelse(is.finite(signv), ifelse(signv>0, "+", "−"), "NA"))
  )
}
.style_edges_delta_topic <- function(ed, alpha_range = .edge_alpha_range,
                                     width_range = .edge_width_range, width_gamma = 1,
                                     cap_abs = NULL){
  delta <- suppressWarnings(as.numeric(ed$score_str) - as.numeric(ed$score_ctrl))
  delta[!is.finite(delta)] <- 0
  mag <- abs(delta)

  cap <- if (is.null(cap_abs)) {
    x <- mag[mag > 0]
    q <- stats::quantile(x, 0.995, na.rm = TRUE, names = FALSE)
    if (!is.finite(q) || q <= 0) q <- suppressWarnings(max(x, na.rm = TRUE))
    if (!is.finite(q) || q <= 0) 1 else q
  } else cap_abs

  s <- pmin(mag, cap) / pmax(cap, 1e-9)
  if (!isTRUE(all.equal(width_gamma, 1))) s <- s^width_gamma
  alpha <- alpha_range[1] + s * (alpha_range[2] - alpha_range[1])

  width <- .fixed_delta_width_bin(mag)
  base_hex <- ifelse(delta >= 0, .col_edge_pos, .col_edge_neg)
  col      <- mapply(.to_rgba, base_hex, alpha, USE.NAMES = FALSE)

  link_sign <- ifelse(is.na(ed$sign_ctrl), ed$sign_str, ed$sign_ctrl)
  link_sign <- ifelse(link_sign >= 0, 1L, ifelse(link_sign < 0, -1L, NA_integer_))

  out <- tibble::tibble(
    id    = paste(ed$TF, ed$gene_key, ed$peak_id, "delta", sep = "|"),
    from  = ed$TF,
    to    = ed$gene_key,
    width = width,
    color = col,
    dashes = FALSE,
    .link_sign = link_sign,
    .delta     = delta,
    title = sprintf(
      "%s → %s<br>peak: %s<br>Δ score (str − ctrl) = %.3f<br>ctrl: %.3f (sign %s), str: %.3f (sign %s)",
      ed$TF, ed$gene_key, ed$peak_id,
      delta,
      ed$score_ctrl, ifelse(is.na(ed$sign_ctrl), "NA", ifelse(ed$sign_ctrl>0, "+", "−")),
      ed$score_str,  ifelse(is.na(ed$sign_str),  "NA", ifelse(ed$sign_str >0, "+", "−"))
    )
  )
  fanned <- .fan_multi_edges(out[, c("id","from","to","width","color","dashes","title")])
  sign_map <- stats::setNames(out$.link_sign, out$id)
  fanned$arrows <- I(lapply(fanned$id, function(k){
    sg <- sign_map[[k]]
    if (is.na(sg) || sg >= 0) list(to = list(enabled = TRUE, type = "arrow"))
    else                      list(to = list(enabled = TRUE, type = "bar"))
  }))
  fanned
}

# Pseudo-peak nodes & unique-peak nodes
.add_peak_pseudonodes <- function(nodes, coords, ed, peak_prefix = "PEAK:",
                                  w_gene = 0.60, r_peak = NULL, tf_ang_map = NULL,
                                  jitter_sd_px = 6){
  tf_per_peak <- ed |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(n_tf = dplyr::n_distinct(TF), .groups = "drop")
  multi_peaks <- tf_per_peak$peak_id[tf_per_peak$n_tf > 1]
  if (!length(multi_peaks)) {
    return(list(nodes = nodes, coords = coords,
                peak_ids = character(0),
                bridge_pairs = tibble::tibble()))
  }

  peak_nodes <- tibble::tibble(
    id = paste0(peak_prefix, multi_peaks),
    type = "peak",
    tf_expr_ctrl = NA_real_, tf_expr_str = NA_real_,
    gene_expr_ctrl = NA_real_, gene_expr_str = NA_real_,
    n_active_ctrl = 0L, n_active_str = 0L,
    l2fc = NA_real_, de_flag_border = "none"
  )
  nodes2 <- dplyr::bind_rows(nodes, peak_nodes) |>
    dplyr::distinct(id, .keep_all = TRUE)

  assoc <- ed |>
    dplyr::filter(peak_id %in% multi_peaks) |>
    dplyr::distinct(peak_id, TF, gene_key)

  if (isTRUE(is.numeric(r_peak)) && length(r_peak) == 1 && is.finite(r_peak) &&
      !is.null(tf_ang_map) && length(tf_ang_map)) {
    ang_pk <- assoc |>
      dplyr::group_by(peak_id) |>
      dplyr::summarise(
        a = {
          aa <- unname(tf_ang_map[intersect(TF, names(tf_ang_map))])
          if (!length(aa) || !any(is.finite(aa))) stats::runif(1, 0, 2*pi) else Arg(mean(exp(1i * aa)))
        },
        .groups = "drop"
      )
    pos_pk <- ang_pk |>
      dplyr::mutate(
        id = paste0(peak_prefix, peak_id),
        x  = r_peak * cos(a) + stats::rnorm(dplyr::n(), 0, jitter_sd_px),
        y  = r_peak * sin(a) + stats::rnorm(dplyr::n(), 0, jitter_sd_px)
      ) |>
      dplyr::select(id, x, y)
  } else {
    pos_tf   <- coords |>
      dplyr::select(id, x, y) |>
      dplyr::rename(TF = id)
    pos_gene <- coords |>
      dplyr::select(id, x, y) |>
      dplyr::rename(gene_key = id)

    agg_tf <- assoc |>
      dplyr::left_join(pos_tf, by = "TF") |>
      dplyr::group_by(peak_id) |>
      dplyr::summarise(x_tf = mean(x, na.rm = TRUE),
                       y_tf = mean(y, na.rm = TRUE), .groups = "drop")
    agg_g  <- assoc |>
      dplyr::left_join(pos_gene, by = "gene_key") |>
      dplyr::group_by(peak_id) |>
      dplyr::summarise(x_g = mean(x, na.rm = TRUE),
                       y_g = mean(y, na.rm = TRUE), .groups = "drop")

    pos_pk <- dplyr::full_join(agg_tf, agg_g, by = "peak_id") |>
      dplyr::mutate(
        x = dplyr::case_when(is.finite(x_tf) & is.finite(x_g) ~ (1 - w_gene) * x_tf + w_gene * x_g,
                             is.finite(x_g) ~ x_g,
                             TRUE ~ x_tf),
        y = dplyr::case_when(is.finite(y_tf) & is.finite(y_g) ~ (1 - w_gene) * y_tf + w_gene * y_g,
                             is.finite(y_g) ~ y_g,
                             TRUE ~ y_tf),
        id = paste0(peak_prefix, peak_id)
      ) |>
      dplyr::select(id, x, y)
  }

  coords2 <- dplyr::bind_rows(coords, pos_pk) |>
    dplyr::distinct(id, .keep_all = TRUE)

  bridge_pairs <- assoc |>
    dplyr::distinct(peak_id, gene_key) |>
    dplyr::mutate(peak_node = paste0(peak_prefix, peak_id))

  list(nodes = nodes2, coords = coords2,
       peak_ids = paste0(peak_prefix, multi_peaks),
       bridge_pairs = bridge_pairs)
}
.add_unique_peak_pseudonodes <- function(nodes, coords, ed, upeak_prefix = "UPEAK:",
                                         w_gene = 0.60, jitter_sd_px = 4){
  tf_per_peak <- ed |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(n_tf = dplyr::n_distinct(TF), .groups = "drop")
  uniq_peaks <- tf_per_peak$peak_id[tf_per_peak$n_tf == 1]
  if (!length(uniq_peaks)) {
    return(list(nodes = nodes, coords = coords,
                upeak_ids = character(0),
                bridge_pairs = tibble::tibble()))
  }

  assoc <- ed |>
    dplyr::filter(peak_id %in% uniq_peaks) |>
    dplyr::distinct(peak_id, TF, gene_key)

  pos_tf   <- coords |>
    dplyr::select(id, x, y) |>
    dplyr::rename(TF = id)
  pos_gene <- coords |>
    dplyr::select(id, x, y) |>
    dplyr::rename(gene_key = id)

  agg_tf <- assoc |>
    dplyr::left_join(pos_tf,   by = "TF") |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(x_tf = mean(x, na.rm = TRUE),
                     y_tf = mean(y, na.rm = TRUE), .groups = "drop")
  agg_g  <- assoc |>
    dplyr::left_join(pos_gene, by = "gene_key") |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(x_g = mean(x, na.rm = TRUE),
                     y_g = mean(y, na.rm = TRUE), .groups = "drop")

  pos_up <- dplyr::full_join(agg_tf, agg_g, by = "peak_id") |>
    dplyr::mutate(
      x = dplyr::case_when(is.finite(x_tf) & is.finite(x_g) ~ (1 - w_gene) * x_tf + w_gene * x_g,
                           is.finite(x_g) ~ x_g,
                           TRUE ~ x_tf),
      y = dplyr::case_when(is.finite(y_tf) & is.finite(y_g) ~ (1 - w_gene) * y_tf + w_gene * y_g,
                           is.finite(y_g) ~ y_g,
                           TRUE ~ y_tf),
      id = paste0(upeak_prefix, peak_id)
    ) |>
    dplyr::transmute(id, x = x + stats::rnorm(dplyr::n(), 0, jitter_sd_px),
                     y = y + stats::rnorm(dplyr::n(), 0, jitter_sd_px))

  u_nodes <- tibble::tibble(
    id = pos_up$id,
    type = "peak",
    tf_expr_ctrl = NA_real_, tf_expr_str = NA_real_,
    gene_expr_ctrl = NA_real_, gene_expr_str = NA_real_,
    n_active_ctrl = 0L, n_active_str = 0L,
    l2fc = NA_real_, de_flag_border = "none"
  )
  nodes2  <- dplyr::bind_rows(nodes, u_nodes) |>
    dplyr::distinct(id, .keep_all = TRUE)
  coords2 <- dplyr::bind_rows(coords, pos_up)  |>
    dplyr::distinct(id, .keep_all = TRUE)

  bridge_pairs <- assoc |>
    dplyr::distinct(peak_id, gene_key) |>
    dplyr::mutate(peak_node = paste0(upeak_prefix, peak_id))

  list(nodes = nodes2, coords = coords2,
       upeak_ids = pos_up$id,
       bridge_pairs = bridge_pairs)
}

.neutral_bridge_edges <- function(bridge_pairs, alpha = 0.60, w = 0.8){
  if (!nrow(bridge_pairs)) return(tibble::tibble(
    id=character(0), from=character(0), to=character(0), width=numeric(0),
    color=character(0), dashes=logical(0), title=character(0),
    arrows=I(list())
  ))
  col <- .to_rgba(.col_edge_zero, alpha)
  tibble::tibble(
    id    = paste("bridge", bridge_pairs$peak_node, bridge_pairs$gene_key, sep = "|"),
    from  = bridge_pairs$peak_node,
    to    = bridge_pairs$gene_key,
    width = w,
    color = col,
    dashes = FALSE,
    title  = sprintf("peak %s → %s", gsub("^PEAK:","", bridge_pairs$peak_node), bridge_pairs$gene_key),
    arrows = I(replicate(nrow(bridge_pairs), list(to = list(enabled = FALSE)), simplify = FALSE))
  )
}

.style_tf_peak_edges_rbin <- function(ed, side = c("ctrl","str"), gray_hex = .col_edge_zero){
  side <- match.arg(side)
  rcol <- if (side == "ctrl") "r_tf_ctrl" else "r_tf_str"

  df <- ed |>
    dplyr::group_by(TF, peak_id) |>
    dplyr::summarise(r = suppressWarnings(max(.data[[rcol]], na.rm = TRUE)), .groups = "drop")
  df$r[!is.finite(df$r)] <- -Inf

  width <- dplyr::case_when(
    is.finite(df$r) & df$r > 0.9 ~ 4,
    is.finite(df$r) & df$r > 0.7 ~ 3,
    is.finite(df$r) & df$r > 0.5 ~ 2,
    is.finite(df$r) & df$r > 0.3 ~ 1,
    TRUE ~ 1
  )
  tibble::tibble(
    id    = paste("tfpeak", df$TF, df$peak_id, side, sep="|"),
    from  = df$TF,
    to    = paste0("PEAK:", df$peak_id),
    width = width,
    color = .to_rgba(gray_hex, 0.85),
    dashes = FALSE,
    title  = sprintf("%s → peak %s<br>r_tf_%s = %s",
                     df$TF, df$peak_id, side,
                     ifelse(is.finite(df$r), sprintf('%.3f', df$r), 'NA')),
    arrows = I(replicate(nrow(df), list(to = list(enabled = FALSE)), simplify = FALSE))
  )
}
.style_tf_upeak_edges_rbin <- function(ed, side = c("ctrl","str"), gray_hex = .col_edge_zero){
  side <- match.arg(side)
  rcol <- if (side == "ctrl") "r_tf_ctrl" else "r_tf_str"

  tf_per_peak <- ed |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(n_tf = dplyr::n_distinct(TF), .groups = "drop")
  uniq_peaks <- tf_per_peak$peak_id[tf_per_peak$n_tf == 1]

  df <- ed |>
    dplyr::filter(peak_id %in% uniq_peaks) |>
    dplyr::group_by(TF, peak_id) |>
    dplyr::summarise(r = suppressWarnings(max(.data[[rcol]], na.rm = TRUE)),
                     .groups = "drop")
  df$r[!is.finite(df$r)] <- -Inf

  width <- dplyr::case_when(
    is.finite(df$r) & df$r > 0.9 ~ 4,
    is.finite(df$r) & df$r > 0.7 ~ 3,
    is.finite(df$r) & df$r > 0.5 ~ 2,
    is.finite(df$r) & df$r > 0.3 ~ 1,
    TRUE ~ 1
  )
  tibble::tibble(
    id    = paste("tfupeak", df$TF, df$peak_id, side, sep="|"),
    from  = df$TF,
    to    = paste0("UPEAK:", df$peak_id),
    width = width,
    color = .to_rgba(gray_hex, 0.85),
    dashes = FALSE,
    title  = sprintf("%s → peak %s (unique)<br>r_tf_%s = %s",
                     df$TF, df$peak_id, side,
                     ifelse(is.finite(df$r), sprintf('%.3f', df$r), 'NA')),
    arrows = I(replicate(nrow(df), list(to = list(enabled = FALSE)), simplify = FALSE))
  )
}

.style_peak_gene_one_side <- function(pg, side = c("ctrl","str"),
                                      caps,
                                      alpha_range = .edge_alpha_range,
                                      width_range = .edge_width_range,
                                      width_gamma = .edge_width_gamma,
                                      zero_col = .col_edge_zero, zero_alpha = .zero_edge_alpha,
                                      zero_width = .zero_edge_width, zero_dashes = FALSE,
                                      log_max_joint = NULL){
  side <- match.arg(side)
  s  <- if (side == "ctrl") pg$score_ctrl else pg$score_str
  sg <- if (side == "ctrl") pg$sign_ctrl  else pg$sign_str
  ps <- if (side == "ctrl") pg$pass_ctrl  else pg$pass_str

  sval <- sg * abs(s); sval[!is.finite(sval)] <- 0
  cap  <- ifelse(sval >= 0, caps$pos, abs(caps$neg))
  mag_cap <- pmin(abs(s), cap) / pmax(cap, 1e-9)

  alpha <- alpha_range[1] + (alpha_range[2] - alpha_range[1]) * mag_cap
  alpha[!ps] <- zero_alpha
  base_hex <- ifelse(sg >= 0, .col_edge_pos, .col_edge_neg)
  base_hex[!ps] <- zero_col
  col <- mapply(.to_rgba, base_hex, alpha, USE.NAMES = FALSE)

  w <- rep(zero_width, nrow(pg))
  idx <- which(ps & is.finite(s) & abs(s) > 0)
  if (length(idx)) {
    mag_log <- log2(1 + abs(s[idx]))
    max_log <- if (is.null(log_max_joint)) max(mag_log, na.rm = TRUE) else log_max_joint
    if (is.finite(max_log) && max_log > 0) {
      sc <- (mag_log / max_log)^width_gamma
      w[idx] <- width_range[1] + sc * (width_range[2] - width_range[1])
    } else w[idx] <- mean(width_range)
  }
  tibble::tibble(
    id    = paste("pk2g", pg$peak_id, pg$gene_key, side, sep="|"),
    from  = paste0("PEAK:", pg$peak_id),
    to    = pg$gene_key,
    width = w,
    color = col,
    dashes = zero_dashes,
    title = sprintf("peak %s → %s<br>%s: score=%s sign=%s",
                    pg$peak_id, pg$gene_key, side,
                    ifelse(is.finite(s), sprintf("%.3f", s), "NA"),
                    ifelse(is.finite(sg), ifelse(sg>0, "+", "−"), "NA")),
    arrows = I(replicate(nrow(pg), list(to = list(enabled = FALSE)), simplify = FALSE))
  )
}
.style_edges_delta_pg <- function(pg, alpha_range = .edge_alpha_range,
                                  width_range = .edge_width_range, width_gamma = 1,
                                  cap_abs = NULL){
  delta <- suppressWarnings(as.numeric(pg$score_str) - as.numeric(pg$score_ctrl))
  delta[!is.finite(delta)] <- 0
  mag <- abs(delta)

  cap <- if (is.null(cap_abs)) {
    x <- mag[mag > 0]
    q <- stats::quantile(x, 0.995, na.rm = TRUE, names = FALSE)
    if (!is.finite(q) || q <= 0) q <- suppressWarnings(max(x, na.rm = TRUE))
    if (!is.finite(q) || q <= 0) 1 else q
  } else cap_abs

  s <- pmin(mag, cap) / pmax(cap, 1e-9)
  if (!isTRUE(all.equal(width_gamma, 1))) s <- s^width_gamma
  alpha <- alpha_range[1] + s * (alpha_range[2] - alpha_range[1])

  width <- .fixed_delta_width_bin(mag)
  base_hex <- ifelse(delta >= 0, .col_edge_pos, .col_edge_neg)
  col      <- mapply(.to_rgba, base_hex, alpha, USE.NAMES = FALSE)

  link_sign <- ifelse(is.na(pg$sign_ctrl), pg$sign_str, pg$sign_ctrl)
  link_sign <- ifelse(link_sign >= 0, 1L, ifelse(link_sign < 0, -1L, NA_integer_))

  tibble::tibble(
    id    = paste("pk2g", pg$peak_id, pg$gene_key, "delta", sep="|"),
    from  = paste0("PEAK:", pg$peak_id),
    to    = pg$gene_key,
    width = width,
    color = col,
    dashes = FALSE,
    arrows = ifelse(is.na(link_sign) | link_sign >= 0, "to", "to;bar"),
    title = sprintf(
      "peak %s → %s<br>Δ score (str − ctrl) = %.3f<br>ctrl: %.3f (sign %s), str: %.3f (sign %s)",
      pg$peak_id, pg$gene_key,
      delta,
      pg$score_ctrl, ifelse(is.na(pg$sign_ctrl), "NA", ifelse(pg$sign_ctrl>0, "+", "−")),
      pg$score_str,  ifelse(is.na(pg$sign_str),  "NA", ifelse(pg$sign_str >0, "+", "−"))
    )
  )
}
# Reuse the same styling as shared peak→gene Δ edges
.style_edges_delta_upg <- function(pg, alpha_range = .edge_alpha_range,
                                   width_range = .edge_width_range, width_gamma = 1,
                                   cap_abs = NULL){
  .style_edges_delta_pg(pg, alpha_range = alpha_range,
                        width_range = width_range, width_gamma = width_gamma,
                        cap_abs = cap_abs)
}
.drop_direct_edges_covered_by_peaks <- function(dir_edges, edsub, peak_nodes){
  if (!nrow(dir_edges) || !length(peak_nodes)) return(dir_edges)
  bare <- unique(gsub("^(PEAK|UPEAK):", "", peak_nodes))
  covered_tf_gene <- edsub |>
    dplyr::filter(peak_id %in% bare) |>
    dplyr::distinct(TF, gene_key) |>
    dplyr::rename(from = TF, to = gene_key)
  dplyr::anti_join(dir_edges, covered_tf_gene, by = c("from","to"))
}

# Layout & scene construction
.bump_peak_nodes_radially <- function(coords, nodes, push_px = 22){
  if (!nrow(coords) || !nrow(nodes)) return(coords)
  pk_ids <- nodes$id[nodes$type == "peak"]
  if (!length(pk_ids)) return(coords)
  i <- match(pk_ids, coords$id); i <- i[is.finite(i)]
  if (!length(i)) return(coords)
  x <- coords$x[i]; y <- coords$y[i]
  r <- sqrt(x*x + y*y); th <- atan2(y, x)
  r2 <- r + push_px
  coords$x[i] <- r2 * cos(th)
  coords$y[i] <- r2 * sin(th)
  coords
}

# (Large) layout routine — preserved from your script with identical math/ordering.
.layout_root_tf_circles <- function(
    nodes,
    child_tfs,
    common_sig, common_nonsig,
    uniq_sig, uniq_nonsig,
    ed_sub,
    root_uniq_sig    = character(0),
    root_uniq_nonsig = character(0),
    root_wedge_center_deg = 270,
    root_wedge_span_deg   = 40,
    tf_wedge_guard_deg    = 12,
    gene_min_sep_px = .layout_gene_min_sep_px,
    ring_gap_px     = 34,
    tf_sep_mult     = 1.8,
    peak_sep_mult   = 1.0,
    base_r2_min     = 160,
    inner_ring_stretch = 5,
    outer_ring_stretch = 5,
    max_iter_nudge  = 4,
    layout_mode   = c("ring","kk_relax"),
    kk_maxiter    = 300,
    tf_min_r_frac = 0.75,
    tf_max_r_frac = 0.96
){
  layout_mode <- match.arg(layout_mode)

  .circmean <- function(a) if (!length(a)) NA_real_ else Arg(mean(exp(1i * a)))
  .even_angles <- function(n, start = 0) if (n <= 0) numeric(0) else (start + (0:(n-1)) * 2*pi / n) %% (2*pi)
  .even_angles_excluding <- function(n, hole_center, hole_span){
    if (n <= 0) return(numeric(0))
    hole_span <- max(0, min(2*pi, hole_span))
    if (hole_span >= 2*pi - 1e-6) return(.even_angles(n, 0))
    start <- (hole_center + hole_span/2) %% (2*pi)
    L <- 2*pi - hole_span
    (start + (0:(n-1)) * L / n) %% (2*pi)
  }
  .even_angles_in_span <- function(n, center, span){
    if (n <= 0) return(numeric(0))
    span <- max(0, min(2*pi, span))
    if (span <= 0) return(numeric(0))
    start <- (center - span/2) %% (2*pi)
    (start + ((seq_len(n) - 0.5) * (span/n))) %% (2*pi)
  }
  .project_into_wedge <- function(a, center, span, margin = 0.02){
    d <- atan2(sin(a - center), cos(a - center))
    d <- pmin(pmax(d, -span/2 + margin), span/2 - margin)
    (center + d) %% (2*pi)
  }
  .clamp_outside_wedge <- function(a, c0, span, margin = 0.02){
    d <- atan2(sin(a - c0), cos(a - c0))
    inside <- abs(d) <= span/2
    ifelse(inside, (c0 + sign(d) * (span/2 + margin)) %% (2*pi), a)
  }
  .min_r_for <- function(n, sep_px){ if (!is.finite(n) || n <= 1) return(sep_px * 0.6); (n * sep_px) / (2*pi) }
  .resolve_collisions <- function(pos, min_sep_px = gene_min_sep_px, max_iter = max_iter_nudge){
    if (!nrow(pos)) return(pos)
    for (it in seq_len(max_iter)) {
      moved <- FALSE
      for (i in seq_len(nrow(pos)-1L)) {
        xi <- pos$x[i]; yi <- pos$y[i]
        for (j in (i+1L):nrow(pos)) {
          dx <- pos$x[j] - xi; dy <- pos$y[j] - yi; d <- sqrt(dx*dx + dy*dy)
          if (!is.finite(d) || d >= min_sep_px) next
          if (d < 1e-6) { ang <- stats::runif(1, 0, 2*pi); dx <- cos(ang); dy <- sin(ang); d <- 1e-6 }
          push <- 0.5 * (min_sep_px - d); ux <- dx/d; uy <- dy/d
          pos$x[i] <- pos$x[i] - ux*push; pos$y[i] <- pos$y[i] - uy*push
          pos$x[j] <- pos$x[j] + ux*push; pos$y[j] <- pos$y[j] + uy*push
          moved <- TRUE
        }
      }
      if (!moved) break
    }
    pos
  }
  .geom_median <- function(xx, yy, tol = 1e-6, maxit = 200){
    if (!length(xx)) return(c(0,0))
    x <- mean(xx, na.rm = TRUE); y <- mean(yy, na.rm = TRUE)
    for (k in seq_len(maxit)){
      dx <- x - xx; dy <- y - yy; d  <- sqrt(dx*dx + dy*dy)
      if (any(d < tol, na.rm = TRUE)) { i <- which.min(d); return(c(xx[i], yy[i])) }
      w <- 1 / pmax(d, tol)
      x_new <- sum(w * xx)/sum(w); y_new <- sum(w * yy)/sum(w)
      if (sqrt((x_new - x)^2 + (y_new - y)^2) < tol) return(c(x_new, y_new))
      x <- x_new; y <- y_new
    }
    c(x, y)
  }
  .api_deg_to_math_rad <- function(deg_api) ((deg_api + 180) %% 360) * pi/180
  .deg_to_rad          <- function(deg) (deg %% 360) * pi/180

  tf_all_non_root <- nodes$id[nodes$type == "TF" & !nodes$is_root]
  root_id <- nodes$id[nodes$is_root]

  tf_with_outgoing <- intersect(tf_all_non_root, unique(ed_sub$TF))
  tf_only_target   <- setdiff(tf_all_non_root, tf_with_outgoing)
  tf_only_target   <- intersect(tf_only_target, unique(ed_sub$gene_key))

  pk_multi <- ed_sub |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(n_tf = dplyr::n_distinct(TF), .groups="drop") |>
    dplyr::filter(n_tf > 1) |>
    dplyr::pull(peak_id)

  tf_pool_for_rings <- setdiff(tf_all_non_root, tf_only_target)
  tf_with_shared <- ed_sub |>
    dplyr::filter(peak_id %in% pk_multi) |>
    dplyr::pull(TF) |>
    unique() |>
    intersect(tf_pool_for_rings)
  tf_no_shared <- setdiff(tf_pool_for_rings, tf_with_shared)

  root_uniq <- unique(c(root_uniq_sig, root_uniq_nonsig))
  root_uniq <- intersect(root_uniq, nodes$id[nodes$type != "TF"])

  genes_all <- union(setdiff(nodes$id[nodes$type != "TF"], root_uniq),
                     tf_only_target)

  regs_by_gene <- ed_sub |>
    dplyr::mutate(peak_node = ifelse(peak_id %in% pk_multi,
                                     paste0("PEAK:", peak_id), NA_character_)) |>
    dplyr::group_by(gene_key) |>
    dplyr::summarise(regs = list(unique(c(TF, stats::na.omit(peak_node)))), .groups="drop") |>
    dplyr::filter(gene_key %in% genes_all)

  nreg <- vapply(regs_by_gene$regs, length, 1L)
  genes_coreg  <- regs_by_gene$gene_key[nreg >= 2]
  genes_unique <- setdiff(genes_all, genes_coreg)

  wedge_c <- .api_deg_to_math_rad(root_wedge_center_deg)
  wedge_t <- .deg_to_rad(root_wedge_span_deg + 2*tf_wedge_guard_deg)
  wedge_g <- .deg_to_rad(root_wedge_span_deg)

  th_shared <- .even_angles_excluding(length(tf_with_shared), wedge_c, wedge_t)
  th_solo   <- .even_angles_excluding(length(tf_no_shared),   wedge_c, wedge_t)

  tf_sep   <- max(gene_min_sep_px, gene_min_sep_px * tf_sep_mult)
  pk_sep   <- gene_min_sep_px * peak_sep_mult
  g_sep    <- gene_min_sep_px

  n_r2 <- length(tf_with_shared)
  n_r3 <- length(pk_multi)
  n_r4 <- length(genes_coreg)
  n_r5 <- length(tf_no_shared)
  n_r6 <- length(unique(c(genes_unique, root_uniq)))

  r2_min_need <- .min_r_for(n_r2, tf_sep)
  r3_min_need <- .min_r_for(n_r3, pk_sep)
  r4_min_need <- .min_r_for(n_r4, g_sep)
  r5_min_need <- .min_r_for(n_r5, tf_sep)
  r6_min_need <- .min_r_for(n_r6, g_sep)

  r2 <- max(base_r2_min, r2_min_need)
  r3 <- max(r2 + ring_gap_px * inner_ring_stretch, r3_min_need)
  r4 <- max(r3 + ring_gap_px * inner_ring_stretch, r4_min_need)
  r5 <- max(r4 + ring_gap_px * outer_ring_stretch, r5_min_need)
  r6 <- max(r5 + ring_gap_px * outer_ring_stretch, r6_min_need)

  pos_root <- tibble::tibble(id = root_id, x = 0, y = 0)
  pos_tf_sh <- tibble::tibble(
    id = tf_with_shared,
    x  = r2 * cos(th_shared) + stats::rnorm(length(tf_with_shared),0,3),
    y  = r2 * sin(th_shared) + stats::rnorm(length(tf_with_shared),0,3)
  )
  pos_tf_so <- tibble::tibble(
    id = tf_no_shared,
    x  = r5 * cos(th_solo) + stats::rnorm(length(tf_no_shared),0,3),
    y  = r5 * sin(th_solo) + stats::rnorm(length(tf_no_shared),0,3)
  )

  tf_pos <- dplyr::bind_rows(pos_tf_sh, pos_tf_so, pos_root)
  tf_xy  <- split(tf_pos[,c("x","y")], tf_pos$id)

  .coreg_xy <- function(gid){
    regs <- regs_by_gene$regs[[match(gid, regs_by_gene$gene_key)]]
    tf_regs <- intersect(regs, names(tf_xy))
    if (length(tf_regs) == 0) {
      th <- (wedge_c + pi) %% (2*pi)
      return(c(r4 * cos(th), r4 * sin(th)))
    }
    if (length(tf_regs) == 2) {
      x1 <- tf_xy[[tf_regs[1]]]$x; y1 <- tf_xy[[tf_regs[1]]]$y
      x2 <- tf_xy[[tf_regs[2]]]$x; y2 <- tf_xy[[tf_regs[2]]]$y
      mx <- (x1 + x2)/2; my <- (y1 + y2)/2
      ang <- atan2(my, mx); ang <- .clamp_outside_wedge(ang, wedge_c, wedge_g, margin = 0.02)
      r   <- sqrt(mx*mx + my*my); r <- min(max(r, r3 + 0.25*(r4 - r3)), r4)
      return(c(r * cos(ang), r * sin(ang)))
    }
    xx <- vapply(tf_regs, function(k) tf_xy[[k]]$x, numeric(1))
    yy <- vapply(tf_regs, function(k) tf_xy[[k]]$y, numeric(1))
    gm <- .geom_median(xx, yy)
    ang <- atan2(gm[2], gm[1]); ang <- .clamp_outside_wedge(ang, wedge_c, wedge_g, margin = 0.02)
    r   <- sqrt(sum(gm^2));      r   <- min(max(r, r3 + 0.25*(r4 - r3)), r4)
    c(r * cos(ang), r * sin(ang))
  }
  if (length(genes_coreg)) {
    gm_mat <- t(vapply(genes_coreg, .coreg_xy, numeric(2)))
    pos_coreg <- tibble::tibble(
      id  = genes_coreg,
      x   = gm_mat[,1] + stats::rnorm(nrow(gm_mat), 0, 3.5),
      y   = gm_mat[,2] + stats::rnorm(nrow(gm_mat), 0, 3.5)
    )
  } else pos_coreg <- tibble::tibble(id = character(0), x = numeric(0), y = numeric(0))

  ang_tf_map <- with(dplyr::bind_rows(pos_tf_sh, pos_tf_so), stats::setNames(atan2(y, x), id))
  .gene_anchor <- function(gid){
    regs <- regs_by_gene$regs[[match(gid, regs_by_gene$gene_key)]]
    aa <- unname(ang_tf_map[intersect(regs, names(ang_tf_map))]); aa <- aa[is.finite(aa)]
    if (!length(aa)) stats::runif(1,0,2*pi) else {
      d <- Arg(mean(exp(1i * aa))); d
    }
  }
  ang_gene_unique <- stats::setNames(vapply(genes_unique, .gene_anchor, 0.0), genes_unique)
  ang_gene_unique <- .clamp_outside_wedge(ang_gene_unique, wedge_c, wedge_g)

  pos_unique <- tibble::tibble(
    id = names(ang_gene_unique),
    x  = r6 * cos(unname(ang_gene_unique)) + stats::rnorm(length(ang_gene_unique),0,4),
    y  = r6 * sin(unname(ang_gene_unique)) + stats::rnorm(length(ang_gene_unique),0,4)
  )

  th_root <- {
    nru <- length(root_uniq)
    if (nru == 0) numeric(0) else {
      inner <- max(0, (root_wedge_span_deg - 8))
      .even_angles_in_span(nru, wedge_c, .deg_to_rad(inner))
    }
  }
  pos_rootuniq <- tibble::tibble(
    id = root_uniq,
    x  = r6 * cos(th_root) + stats::rnorm(length(th_root),0,4),
    y  = r6 * sin(th_root) + stats::rnorm(length(th_root),0,4)
  )

  peak_regs <- ed_sub |>
    dplyr::filter(peak_id %in% pk_multi) |>
    dplyr::group_by(peak_id) |>
    dplyr::summarise(tfs = list(unique(TF)), .groups = "drop")

  pos_peaks <- tibble::tibble(id = character(0), x = numeric(0), y = numeric(0))
  if (nrow(peak_regs)) {
    r_lo <- r2 + 0.35 * (r3 - r2)
    r_hi <- r4 - 0.35 * (r4 - r3)
    if (!is.finite(r_lo) || !is.finite(r_hi) || r_hi <= r_lo) {
      r_lo <- r3 * 0.9; r_hi <- r4 * 0.95
    }
    tf_xy <- tf_xy
    pk_xy <- lapply(seq_len(nrow(peak_regs)), function(i){
      tfi <- intersect(peak_regs$tfs[[i]], names(tf_xy))
      if (!length(tfi)) {
        th <- stats::runif(1, 0, 2*pi)
        return(c(r3 * cos(th), r3 * sin(th)))
      }
      if (length(tfi) == 2) {
        x1 <- tf_xy[[tfi[1]]]$x; y1 <- tf_xy[[tfi[1]]]$y
        x2 <- tf_xy[[tfi[2]]]$x; y2 <- tf_xy[[tfi[2]]]$y
        mx <- (x1 + x2)/2; my <- (y1 + y2)/2
        ang <- atan2(my, mx); ang <- .clamp_outside_wedge(ang, wedge_c, wedge_g, margin = 0.02)
        r   <- sqrt(mx*mx + my*my); r <- min(max(r, r_lo), r_hi)
        return(c(r * cos(ang), r * sin(ang)))
      } else {
        xx <- vapply(tfi, function(k) tf_xy[[k]]$x, numeric(1))
        yy <- vapply(tfi, function(k) tf_xy[[k]]$y, numeric(1))
        gm <- .geom_median(xx, yy)
        ang <- atan2(gm[2], gm[1]); ang <- .clamp_outside_wedge(ang, wedge_c, wedge_g, margin = 0.02)
        r   <- sqrt(sum(gm^2));  r <- min(max(r, r_lo), r_hi)
        return(c(r * cos(ang), r * sin(ang)))
      }
    })
    pk_xy <- do.call(rbind, pk_xy)
    pos_peaks <- tibble::tibble(
      id = paste0("PEAK:", peak_regs$peak_id),
      x  = pk_xy[,1] + stats::rnorm(nrow(peak_regs), 0, 3.0),
      y  = pk_xy[,2] + stats::rnorm(nrow(peak_regs), 0, 3.0)
    )
  }

  coords <- dplyr::bind_rows(pos_root, pos_tf_sh, pos_tf_so, pos_coreg, pos_unique, pos_rootuniq, pos_peaks) |>
    dplyr::distinct(id, .keep_all = TRUE)
  coords <- .resolve_collisions(coords, min_sep_px = gene_min_sep_px, max_iter = max_iter_nudge)

  if (length(root_uniq)) {
    i  <- match(root_uniq, coords$id)
    th <- atan2(coords$y[i], coords$x[i])
    arc_step <- if (is.finite(r6) && r6 > 0) gene_min_sep_px / r6 else 0.02
    th <- th + stats::rnorm(length(i), mean = 0, sd = arc_step * 0.35)
    th <- .project_into_wedge(th, wedge_c, wedge_g, margin = 0.02)
    coords$x[i] <- r6 * cos(th)
    coords$y[i] <- r6 * sin(th)
  }
  coords <- .resolve_collisions(coords, min_sep_px = gene_min_sep_px, max_iter = max_iter_nudge)
  coords
}

.tfhub_build_delta_panel <- function(core, coords, stress_tag = "stress", ctrl_tag = "control"){
  nodes <- core$nodes
  edsub <- core$ed_sub

  tf_ang_map <- with(coords[coords$id %in% nodes$id[nodes$type=="TF"], ], stats::setNames(atan2(y, x), id))
  r_tf_med <- stats::median(sqrt(coords$x[nodes$type=="TF" & !nodes$is_root]^2 +
                                   coords$y[nodes$type=="TF" & !nodes$is_root]^2), na.rm = TRUE)
  r_gene_q <- stats::quantile(sqrt(coords$x[nodes$type!="TF"]^2 + coords$y[nodes$type!="TF"]^2),
                              probs = 0.25, na.rm = TRUE)
  r_peak_ring <- r_tf_med + 0.72 * (r_gene_q - r_tf_med) + 12

  pk  <- .add_peak_pseudonodes(nodes, coords, edsub, r_peak = r_peak_ring, tf_ang_map = tf_ang_map)
  nodes  <- pk$nodes
  coords <- pk$coords

  upk <- .add_unique_peak_pseudonodes(nodes, coords, edsub, w_gene = 0.08, jitter_sd_px = 2)
  nodes  <- upk$nodes
  coords <- upk$coords
  coords <- .bump_peak_nodes_radially(coords, nodes, push_px = 22)

  tf_idx  <- nodes$type == "TF"
  tf_l2fc <- nodes$l2fc[tf_idx]
  max_abs <- suppressWarnings(max(abs(tf_l2fc), na.rm = TRUE)); if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  tf_signed <- rep(NA_real_, nrow(nodes)); tf_signed[tf_idx] <- tf_l2fc / max_abs
  ramp_div <- function(u_signed){
    u_signed <- tidyr::replace_na(u_signed, 0)
    w <- 0.5 + 0.5 * clamp(u_signed, -1, 1)
    grDevices::rgb(grDevices::colorRamp(c(.col_edge_neg, "#FFFFFF", .col_edge_pos))(w)/255)
  }
  gene_z_joint <- pmax(abs(nodes$node_z_ctrl), abs(nodes$node_z_str), na.rm = TRUE)
  gene_size <- {
    x <- pmax(0, as.numeric(gene_z_joint))
    if (!any(is.finite(x))) rep(mean(.gene_size_range), nrow(nodes)) else {
      s <- (x / max(x, na.rm = TRUE))^.gene_size_gamma
      .gene_size_range[1] + s * (.gene_size_range[2] - .gene_size_range[1])
    }
  }
  thr <- log2(1.5)
  border_col <- ifelse(!is.finite(nodes$l2fc) | abs(nodes$l2fc) < thr, .col_border_none,
                       ifelse(nodes$l2fc > 0, .col_border_up, .col_border_down))
  fc_abs <- 2^abs(nodes$l2fc)
  bw_gene <- dplyr::case_when(
    is.finite(fc_abs) & fc_abs >  8 ~ 4,
    is.finite(fc_abs) & fc_abs >  4 ~ 3,
    is.finite(fc_abs) & fc_abs >  2 ~ 2,
    is.finite(fc_abs) & fc_abs >  0 ~ 1,
    TRUE ~ 1
  )
  border_width <- ifelse(nodes$type == "TF", 1, bw_gene)
  border_width[nodes$type == "peak"] <- 0.5

  pos <- coords[match(nodes$id, coords$id), c("x","y")]
  pos$x[!is.finite(pos$x)] <- 0; pos$y[!is.finite(pos$y)] <- 0
  tf_font <- {
    x <- pmax(nodes$tf_expr_ctrl, nodes$tf_expr_str, na.rm = TRUE)
    x <- pmax(0, as.numeric(x))
    if (!any(is.finite(x)) || max(x, na.rm = TRUE) == 0) rep(mean(.tf_font_range), nrow(nodes)) else {
      s <- (x / max(x, na.rm = TRUE))^.tf_font_gamma
      .tf_font_range[1] + s * (.tf_font_range[2] - .tf_font_range[1])
    }
  }
  nd <- data.frame(
    id    = nodes$id,
    label = ifelse(nodes$type == "TF", nodes$id, ""),
    x = pos$x, y = pos$y,
    color.background = ifelse(nodes$type == "TF", ramp_div(tf_signed), .col_gene_fill),
    color.border     = border_col,
    borderWidth      = border_width,
    value = ifelse(nodes$type == "TF", 1,
                   ifelse(nodes$type == "peak", NA, gene_size)),
    size  = ifelse(nodes$type == "peak", .peak_node_px, NA),
    title = ifelse(
      nodes$type == "TF",
      sprintf("TF %s%s%s<br>log2FC(TF RNA): %.3f<br>#active edges: ctrl=%s, str=%s",
              nodes$id,
              ifelse(is.finite(nodes$tf_expr_ctrl), sprintf("<br>expr_ctrl=%.3f", nodes$tf_expr_ctrl), ""),
              ifelse(is.finite(nodes$tf_expr_str ), sprintf("<br>expr_str =%.3f", nodes$tf_expr_str ), ""),
              ifelse(is.finite(nodes$l2fc), nodes$l2fc, NA),
              nodes$n_active_ctrl, nodes$n_active_str),
      ifelse(
        nodes$type == "peak",
        sprintf("Peak %s", gsub("^(UPEAK|PEAK):","", nodes$id)),
        sprintf("Gene %s%s%s",
                nodes$id,
                ifelse(is.finite(nodes$gene_expr_ctrl), sprintf("<br>expr_ctrl=%.3f", nodes$gene_expr_ctrl), ""),
                ifelse(is.finite(nodes$gene_expr_str ), sprintf("<br>expr_str =%.3f", nodes$gene_expr_str ), "")  # <-- add "" here
        )
      )
    ),
    shape = ifelse(nodes$type == "TF", "box",
                   ifelse(nodes$type == "peak", "triangle", "dot")),
    physics = FALSE, fixed = FALSE,
    "font.color"       = "#2b2b2b",
    "font.size"        = ifelse(nodes$type == "TF", tf_font, 0),
    "font.face"        = "arial",
    "font.strokeWidth" = ifelse(nodes$type == "TF", .tf_label_halo_w, 0),
    "font.strokeColor" = ifelse(nodes$type == "TF", .tf_label_halo_c, NA),
    stringsAsFactors = FALSE
  )

  e_tfpk <- .style_tf_peak_edges_rbin(edsub, "ctrl") |>
    dplyr::filter(to %in% pk$peak_ids)
  e_tfu  <- .style_tf_upeak_edges_rbin(edsub, "ctrl") |>
    dplyr::filter(to %in% upk$upeak_ids)

  pg_all   <- .aggregate_peak_gene(edsub)
  pg_share <- pg_all |>
    dplyr::filter(paste0("PEAK:",  .data$peak_id) %in% pk$peak_ids)
  pg_uniq  <- pg_all |>
    dplyr::filter(paste0("UPEAK:", .data$peak_id) %in% upk$upeak_ids)

  delta_all <- abs(suppressWarnings(as.numeric(pg_all$score_str) - as.numeric(pg_all$score_ctrl)))
  cap_joint <- stats::quantile(delta_all[delta_all > 0], 0.995, na.rm = TRUE, names = FALSE)
  if (!is.finite(cap_joint) || cap_joint <= 0) cap_joint <- suppressWarnings(max(delta_all, na.rm = TRUE))
  if (!is.finite(cap_joint) || cap_joint <= 0) cap_joint <- 1

  e_pg_delta  <- .style_edges_delta_pg (pg_share, cap_abs = cap_joint)
  e_upg_delta <- .style_edges_delta_upg(pg_uniq,  cap_abs = cap_joint)

  covered_nodes <- c(pk$peak_ids, upk$upeak_ids)
  e_delta_direct <- .style_edges_delta_topic(edsub, cap_abs = cap_joint) |>
    dplyr::mutate(id = paste("dir", .data$id, sep="|"),
                  dashes = FALSE,
                  color  = gsub(",\\s*[0-9.]+\\)$", ",0.35)", .data$color))

  edv <- dplyr::bind_rows(
    .force_edge_cols(e_tfpk),
    .force_edge_cols(e_tfu),
    .force_edge_cols(e_pg_delta),
    .force_edge_cols(e_upg_delta),
    .force_edge_cols(e_delta_direct)
  )
  list(nodes = nd, edges = edv)
}

# Build core subgraph (selection, stats, styling prep)
.build_nodes_edges_from_root_tf <- function(
    comp_df,
    input_tf,
    edge_filter_min = 2,
    edge_filter_on = c("either","stress","control","both"),
    gene_fc_thresh = 1.5,
    de_reference = c("str_over_ctrl","ctrl_over_str"),
    cond1_tag = NULL,
    cond2_tag = NULL,
    border_up_for = c("stress","control"),
    ring_tf_outgoing = c("context","full","full_topk"),
    ring_tf_topk = 150,
    ring_tf_only_direct = FALSE,
    verbose = TRUE
){
  edge_filter_on <- match.arg(edge_filter_on)
  de_reference   <- match.arg(de_reference)
  border_up_for  <- match.arg(border_up_for)
  ring_tf_outgoing <- match.arg(ring_tf_outgoing)

  cols <- .detect_cols(comp_df, cond1_tag = cond1_tag, cond2_tag = cond2_tag)
  tf_col   <- cols$tf;    gene_col <- cols$gene;    peak_col <- cols$peak
  sc_ctrl  <- cols$score_ctrl;  sc_str <- cols$score_str
  sg_ctrl  <- cols$sign_ctrl;   sg_str <- cols$sign_str
  ac_ctrl  <- cols$active_ctrl; ac_str <- cols$active_str
  tfc      <- cols$tf_expr_ctrl; tfs   <- cols$tf_expr_str
  gxc      <- cols$gene_expr_ctrl; gxs  <- cols$gene_expr_str
  rtc      <- cols$r_tf_ctrl; rts <- cols$r_tf_str

  .llog(
    "mapped cols: score_ctrl='", sc_ctrl, "', score_str='", sc_str,
    "', sign_ctrl='", sg_ctrl, "', sign_str='", sg_str,
    "', active_ctrl='", ac_ctrl %||% "NULL",
    "', active_str='",  ac_str  %||% "NULL", "'",
    verbose = verbose
  )

  ed0 <- comp_df |>
    dplyr::transmute(
      TF       = .data[[tf_col]],
      gene_key = .data[[gene_col]],
      peak_id  = .data[[peak_col]],
      score_ctrl = suppressWarnings(as.numeric(.data[[sc_ctrl]])),
      score_str  = suppressWarnings(as.numeric(.data[[sc_str]])),
      sign_ctrl  = .as_sign(.data[[sg_ctrl]]),
      sign_str   = .as_sign(.data[[sg_str]]),
      active_ctrl = if (!is.null(ac_ctrl)) .as_bool(.data[[ac_ctrl]]) else NA,
      active_str  = if (!is.null(ac_str))  .as_bool(.data[[ac_str]])  else NA,
      tf_expr_ctrl   = if (!is.null(tfc)) suppressWarnings(as.numeric(.data[[tfc]])) else NA_real_,
      tf_expr_str    = if (!is.null(tfs)) suppressWarnings(as.numeric(.data[[tfs]])) else NA_real_,
      gene_expr_ctrl = if (!is.null(gxc)) suppressWarnings(as.numeric(.data[[gxc]])) else NA_real_,
      gene_expr_str  = if (!is.null(gxs)) suppressWarnings(as.numeric(.data[[gxs]])) else NA_real_,
      r_tf_ctrl = if (!is.null(rtc)) suppressWarnings(as.numeric(.data[[rtc]])) else NA_real_,
      r_tf_str  = if (!is.null(rts)) suppressWarnings(as.numeric(.data[[rts]])) else NA_real_
    ) |>
    dplyr::distinct(TF, gene_key, peak_id, .keep_all = TRUE) |>
    dplyr::filter(!is.na(TF), !is.na(gene_key), !is.na(peak_id))

  ed0$TF       <- trimws(as.character(ed0$TF))
  ed0$gene_key <- trimws(as.character(ed0$gene_key))
  input_tf     <- trimws(as.character(input_tf))

  .llog(">> Building TF-hub for '", input_tf, "': total edges=", nrow(ed0), verbose = verbose)

  ed0 <- dplyr::mutate(
    ed0,
    pass_ctrl_score = .compute_pass_flags(dplyr::cur_data(), edge_filter_min, "score_ctrl", col_active = NULL),
    pass_str_score  = .compute_pass_flags(dplyr::cur_data(), edge_filter_min, "score_str",  col_active = NULL),
    pass_ctrl = .compute_pass_flags(dplyr::cur_data(), edge_filter_min, "score_ctrl", "active_ctrl"),
    pass_str  = .compute_pass_flags(dplyr::cur_data(), edge_filter_min, "score_str",  "active_str")
  )

  keep <- switch(edge_filter_on,
                 control = ed0$pass_ctrl_score,
                 stress  = ed0$pass_str_score,
                 both    = ed0$pass_ctrl_score & ed0$pass_str_score,
                 either  = ed0$pass_ctrl_score | ed0$pass_str_score)
  keep[is.na(keep)] <- FALSE
  ed <- ed0[keep, , drop = FALSE]

  .llog("After filters (by SCORE only): edges=", nrow(ed), verbose = verbose)

  known_tfs <- .get_tf_syms() %||% character(0)

  child_tfs <- ed |>
    dplyr::filter(TF == input_tf) |>
    dplyr::pull(gene_key) |>
    unique() |>
    intersect(known_tfs)

  parent_tfs <- ed |>
    dplyr::filter(gene_key == input_tf) |>
    dplyr::pull(TF) |>
    unique() |>
    intersect(known_tfs)

  direct_tfs <- union(child_tfs, parent_tfs)
  tf_ring_core <- direct_tfs

  ext_parents <- ed |>
    dplyr::filter(gene_key %in% direct_tfs) |>
    dplyr::pull(TF) |>
    unique()
  direct_tfs <- union(direct_tfs, intersect(ext_parents, known_tfs))

  tf_in_plot <- unique(c(if (.is_known_tf(input_tf)) input_tf else character(0), direct_tfs))

  # Optionally disable second-order TFs (only direct TFs around the ring)
  second_order_tfs <- if (isTRUE(ring_tf_only_direct)) {
    character(0)
  } else {
    ed |>
      dplyr::filter(TF %in% direct_tfs) |>
      dplyr::pull(gene_key) |>
      unique() |>
      intersect(known_tfs)
  }
  tf_in_plot <- union(tf_in_plot, second_order_tfs)

  root_targets <- ed |>
    dplyr::filter(TF == input_tf) |>
    dplyr::pull(gene_key) |>
    unique()

  common_genes <- character(0); unique_genes <- character(0)

  if (length(direct_tfs)) {
    tf_all_in_plot <- unique(c(if (.is_known_tf(input_tf)) input_tf else character(0),
                               direct_tfs, second_order_tfs))
    tfs_for_context <- if (ring_tf_outgoing == "context") tf_ring_core else setdiff(tf_all_in_plot, input_tf)

    conn_targets <- ed |>
      dplyr::filter(TF %in% tfs_for_context) |>
      dplyr::pull(gene_key) |>
      unique()

    common_genes <- intersect(root_targets, conn_targets)
    unique_genes <- setdiff(conn_targets, root_targets)

    if (ring_tf_outgoing == "context") {
      ed_sub <- ed |>
        dplyr::filter(
          TF == input_tf |
            (TF %in% parent_tfs & gene_key == input_tf) |
            (TF %in% direct_tfs & gene_key %in% c(common_genes, unique_genes)) |
            (TF %in% tf_all_in_plot & gene_key %in% tf_all_in_plot)
        )
    } else {
      ed_full <- ed |>
        dplyr::filter(TF %in% tf_all_in_plot)

      if (ring_tf_outgoing == "full_topk" && is.finite(ring_tf_topk) && ring_tf_topk > 0) {
        ed_full <- ed_full |>
          dplyr::mutate(smax = pmax(abs(score_ctrl), abs(score_str), na.rm = TRUE)) |>
          dplyr::group_by(TF, gene_key) |>
          dplyr::slice_max(order_by = smax, n = 1, with_ties = FALSE) |>
          dplyr::ungroup() |>
          dplyr::group_by(TF) |>
          dplyr::slice_max(order_by = smax, n = ring_tf_topk, with_ties = TRUE) |>
          dplyr::ungroup()
      }

      ed_sub <- dplyr::bind_rows(
        ed_full,
        ed |>
          dplyr::filter(TF %in% parent_tfs, gene_key == input_tf),
        ed |>
          dplyr::filter(TF %in% tf_all_in_plot, gene_key %in% tf_all_in_plot)
      )
    }

  } else {
    ed_sub <- ed |>
      dplyr::filter(TF == input_tf)
    tfs_for_context <- character(0)
  }

  ed_sub <- ed_sub |>
    dplyr::distinct(TF, gene_key, peak_id, .keep_all = TRUE)

  is_gene_like <- function(x) !.is_known_tf(x)
  common_genes <- common_genes[is_gene_like(common_genes)]
  unique_genes <- unique_genes[is_gene_like(unique_genes)]
  root_uniq    <- setdiff(root_targets, c(common_genes, unique_genes))
  root_uniq    <- root_uniq[is_gene_like(root_uniq)]

  any_pass_by_gene <- function(df){
    if (!nrow(df)) return(tibble::tibble(gene_key = character(0), any_pass = logical(0)))
    df |>
      dplyr::mutate(any_pass = pass_ctrl | pass_str) |>
      dplyr::group_by(gene_key) |>
      dplyr::summarise(any_pass = any(any_pass, na.rm = TRUE), .groups = "drop")
  }

  common_df <- any_pass_by_gene(ed_sub |>
                                  dplyr::filter(gene_key %in% common_genes))
  common_sig    <- common_df$gene_key[common_df$any_pass]
  common_nonsig <- setdiff(common_genes, common_sig)

  uniq_df <- any_pass_by_gene(
    ed_sub |>
      dplyr::filter(gene_key %in% unique_genes,
                    TF %in% tfs_for_context)
  )
  uniq_sig    <- uniq_df$gene_key[uniq_df$any_pass]
  uniq_nonsig <- setdiff(unique_genes, uniq_sig)

  runiq_df <- any_pass_by_gene(ed_sub |>
                                 dplyr::filter(gene_key %in% root_uniq,
                                               TF == input_tf))
  root_uniq_sig    <- runiq_df$gene_key[runiq_df$any_pass]
  root_uniq_nonsig <- setdiff(root_uniq, root_uniq_sig)

  node_ids   <- unique(c(input_tf, direct_tfs, common_genes, unique_genes, root_uniq, second_order_tfs))
  node_types <- ifelse(.is_known_tf(node_ids), "TF", "gene")
  nodes <- tibble::tibble(id = node_ids, type = node_types, is_root = node_ids == input_tf)

  tf_ctrl <- ed0 |>
    dplyr::group_by(TF) |>
    dplyr::summarise(tf_expr_ctrl = stats::median(tf_expr_ctrl, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(id = TF)
  tf_str  <- ed0 |>
    dplyr::group_by(TF) |>
    dplyr::summarise(tf_expr_str = stats::median(tf_expr_str, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(id = TF)
  g_ctrl  <- ed0 |>
    dplyr::group_by(gene_key) |>
    dplyr::summarise(gene_expr_ctrl = stats::median(gene_expr_ctrl, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(id = gene_key)
  g_str   <- ed0 |>
    dplyr::group_by(gene_key) |>
    dplyr::summarise(gene_expr_str = stats::median(gene_expr_str, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(id = gene_key)

  nodes <- nodes |>
    dplyr::left_join(tf_ctrl, by = "id") |>
    dplyr::left_join(tf_str,  by = "id") |>
    dplyr::left_join(g_ctrl,  by = "id") |>
    dplyr::left_join(g_str,   by = "id") |>
    dplyr::mutate(
      node_raw_ctrl = dplyr::if_else(type == "TF", dplyr::coalesce(tf_expr_ctrl, 0), dplyr::coalesce(gene_expr_ctrl, 0)),
      node_raw_str  = dplyr::if_else(type == "TF", dplyr::coalesce(tf_expr_str,  0), dplyr::coalesce(gene_expr_str,  0))
    ) |>
    dplyr::mutate(
      node_z_ctrl = robust_z(node_raw_ctrl),
      node_z_str  = robust_z(node_raw_str)
    )

  deg_ctrl <- ed_sub |>
    dplyr::filter(pass_ctrl) |>
    dplyr::count(TF, name = "n_active_ctrl") |>
    dplyr::rename(id = TF)
  deg_str  <- ed_sub |>
    dplyr::filter(pass_str ) |>
    dplyr::count(TF, name = "n_active_str" ) |>
    dplyr::rename(id = TF)

  nodes <- nodes |>
    dplyr::left_join(deg_ctrl, by = "id") |>
    dplyr::left_join(deg_str,  by = "id") |>
    dplyr::mutate(
      n_active_ctrl = tidyr::replace_na(n_active_ctrl, 0L),
      n_active_str  = tidyr::replace_na(n_active_str,  0L)
    )

  .scale_font <- function(x, range = .tf_font_range, gamma = .tf_font_gamma){
    x <- pmax(0, as.numeric(x))
    if (!any(is.finite(x)) || max(x, na.rm = TRUE) == 0) return(rep(mean(range), length(x)))
    s <- (x / max(x, na.rm = TRUE))^gamma
    range[1] + s * (range[2] - range[1])
  }
  nodes$tf_font_ctrl <- NA_real_
  nodes$tf_font_str  <- NA_real_
  tf_idx <- nodes$type == "TF"
  if (any(tf_idx)) {
    nodes$tf_font_ctrl[tf_idx] <- .scale_font(nodes$tf_expr_ctrl[tf_idx])
    nodes$tf_font_str [tf_idx] <- .scale_font(nodes$tf_expr_str [tf_idx])
  }

  eps <- 1e-9; thr <- log2(gene_fc_thresh)
  l2fc_gene <- if (de_reference == "ctrl_over_str") {
    log2((nodes$gene_expr_ctrl + eps)/(nodes$gene_expr_str + eps))
  } else {
    log2((nodes$gene_expr_str  + eps)/(nodes$gene_expr_ctrl + eps))
  }
  l2fc_tf <- if (de_reference == "ctrl_over_str") {
    log2((nodes$tf_expr_ctrl + eps)/(nodes$tf_expr_str + eps))
  } else {
    log2((nodes$tf_expr_str  + eps)/(nodes$tf_expr_ctrl + eps))
  }
  nodes$l2fc <- ifelse(nodes$type == "TF", l2fc_tf, l2fc_gene)

  nodes$de_flag_border <- if (border_up_for == "stress") {
    dplyr::case_when(
      is.finite(nodes$l2fc) & nodes$l2fc <= -thr ~ "up",
      is.finite(nodes$l2fc) & nodes$l2fc >=  thr ~ "down",
      TRUE ~ "none"
    )
  } else {
    dplyr::case_when(
      is.finite(nodes$l2fc) & nodes$l2fc >=  thr ~ "up",
      is.finite(nodes$l2fc) & nodes$l2fc <= -thr ~ "down",
      TRUE ~ "none"
    )
  }

  signed_ctrl <- ed_sub$score_ctrl * ifelse(is.na(ed_sub$sign_ctrl), 0, ed_sub$sign_ctrl)
  signed_str  <- ed_sub$score_str  * ifelse(is.na(ed_sub$sign_str),  0, ed_sub$sign_str)
  caps <- .compute_caps(
    v_pos = c(signed_ctrl[signed_ctrl > 0], signed_str[signed_str > 0]),
    v_neg = c(signed_ctrl[signed_ctrl < 0], signed_str[signed_str < 0]),
    q = 0.995
  )
  v_for_log <- c(abs(ed_sub$score_ctrl[ed_sub$pass_ctrl]), abs(ed_sub$score_str[ed_sub$pass_str]))
  v_for_log <- v_for_log[is.finite(v_for_log) & v_for_log > 0]
  log_max_joint <- if (length(v_for_log)) max(log2(1 + v_for_log)) else 1

  edges_ctrl <- .style_edges_one_side(ed_sub, "ctrl", caps, log_max_joint = log_max_joint) |>
    .fan_multi_edges()
  edges_str  <- .style_edges_one_side(ed_sub, "str",  caps, log_max_joint = log_max_joint) |>
    .fan_multi_edges()

  tf_set <- nodes$id[nodes$type == "TF"]
  edges_ctrl$arrows <- ifelse((edges_ctrl$from %in% tf_set) & (edges_ctrl$to %in% tf_set), "to", "")
  edges_str$arrows  <- ifelse((edges_str$from  %in% tf_set) & (edges_str$to  %in% tf_set),  "to", "")

  list(
    nodes = nodes,
    edges_ctrl = edges_ctrl,
    edges_str  = edges_str,
    caps = caps,
    child_tfs = child_tfs,
    parent_tfs = parent_tfs,
    direct_tfs = direct_tfs,
    common_sig = common_sig,
    common_nonsig = common_nonsig,
    uniq_sig = uniq_sig,
    uniq_nonsig = uniq_nonsig,
    root_uniq_sig    = root_uniq_sig,
    root_uniq_nonsig = root_uniq_nonsig,
    ed_sub = ed_sub,
    input_tf = input_tf
  )
}

# UI helpers
.attach_vis_jumpers <- function(widget, tf_label = "Jump to TF…", gene_label = "Jump to gene…"){
  if (is.null(widget$elementId) || !nzchar(widget$elementId)) {
    widget$elementId <- paste0("vis-", paste(sample(c(letters, 0:9), 8, TRUE), collapse = ""))
  }
  nodes <- widget$x$nodes
  if (is.null(nodes)) return(widget)
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  if (is.null(nodes$id))    nodes$id    <- nodes$label
  if (is.null(nodes$label)) nodes$label <- as.character(nodes$id)

  isTF <- rep(FALSE, nrow(nodes))
  if (!is.null(nodes$kind))  isTF <- isTF | nodes$kind  %in% "TF"
  if (!is.null(nodes$group)) isTF <- isTF | nodes$group %in% "TF"
  if (!is.null(nodes$shape)) isTF <- isTF | tolower(as.character(nodes$shape)) %in% "box"
  isTF[is.na(isTF)] <- FALSE

  isPeak <- rep(FALSE, nrow(nodes))
  if (!is.null(nodes$type))  isPeak <- isPeak | tolower(as.character(nodes$type))  %in% "peak"
  if (!is.null(nodes$group)) isPeak <- isPeak | tolower(as.character(nodes$group)) %in% "peak"
  if (!is.null(nodes$shape)) isPeak <- isPeak | tolower(as.character(nodes$shape)) %in% "triangle"
  if (!is.null(nodes$id))    isPeak <- isPeak | grepl("^(?:U?PEAK:)", as.character(nodes$id))
  if (!is.null(nodes$label)) isPeak <- isPeak | grepl("^(?:U?PEAK:)", as.character(nodes$label))
  isPeak[is.na(isPeak)] <- FALSE

  tfs   <- nodes[isTF, c("id","label"), drop = FALSE]
  genes <- nodes[(!isTF & !isPeak), c("id","label"), drop = FALSE]

  if (nrow(genes)) {
    bad_id <- grepl("^(?:U?PEAK:|pk:)", as.character(genes$id))
    if (any(bad_id)) genes <- genes[!bad_id, , drop = FALSE]
  }

  ord <- function(df) order(tolower(ifelse(is.na(df$label) | df$label == "", as.character(df$id), df$label)))
  if (nrow(tfs))   tfs   <- tfs[ord(tfs),]
  if (nrow(genes)) genes <- genes[ord(genes),]

  make_options <- function(df){
    if (!nrow(df)) return(list())
    apply(df, 1, function(r){
      htmltools::tags$option(
        value = as.character(r[["id"]]),
        ifelse(is.na(r[["label"]]) || r[["label"]] == "", r[["id"]], r[["label"]])
      )
    })
  }

  sel_tf_id   <- paste0(widget$elementId, "-selTF")
  sel_gene_id <- paste0(widget$elementId, "-selGene")

  controls <- htmlwidgets::prependContent(
    widget,
    htmltools::tags$div(
      class = "vis-jumpers",
      style = "margin-bottom:6px;",
      htmltools::tags$select(
        id = sel_tf_id, style = "width:100%;margin:4px 0;display:block;",
        htmltools::tags$option(value = "", tf_label),
        make_options(tfs)
      ),
      htmltools::tags$select(
        id = sel_gene_id, style = "width:100%;margin:4px 0;display:block;",
        htmltools::tags$option(value = "", gene_label),
        make_options(genes)
      )
    )
  )
  js <- sprintf("
    function(el, x){
      var inst = this;
      function getNet(i){ return i && i.network ? i.network : (i && i.body && i.body.network ? i.body.network : null); }
      function selectAndFocus(raw){
        var id = raw;
        var net = getNet(inst); if(!net) return;
        try{
          var ds = net.body && net.body.data && net.body.data.nodes;
          if(!ds || (ds.get(id) === null || ds.get(id) === undefined)){ return; }
          try{ net.setOptions({interaction:{hover:false}}); }catch(e){}
          net.setSelection({nodes:[id], edges:[]}, {unselectAll:true});
          net.fit({ nodes:[id], animation:{ duration:600, easingFunction:'easeInOutQuad' } });
          setTimeout(function(){ try{ net.setOptions({interaction:{hover:true}}); }catch(e){} }, 50);
        } catch(e){}
      }
      var selTF   = document.getElementById(%s);
      var selGene = document.getElementById(%s);
      if(!selTF || !selGene) return;
      selTF.addEventListener('change',  function(){ if(this.value){ selGene.selectedIndex=0; selectAndFocus(this.value); } });
      selGene.addEventListener('change',function(){ if(this.value){ selTF.selectedIndex=0;  selectAndFocus(this.value); } });
    }",
                jsonlite::toJSON(sel_tf_id), jsonlite::toJSON(sel_gene_id)
  )
  htmlwidgets::onRender(controls, js)
}
.attach_vis_sync <- function(widget, role = c("ctrl","str"), pair_id){
  role <- match.arg(role)
  stopifnot(is.character(pair_id), nchar(pair_id) > 0)
  js <- sprintf("
    function(el, x){
      var gid  = %s;
      var role = %s;
      window.__visPairs = window.__visPairs || {};
      var store = window.__visPairs[gid] = window.__visPairs[gid] || {};
      function getNet(inst){ return inst && inst.network ? inst.network :
               (inst && inst.body && inst.body.network ? inst.body.network : null); }
      store[role] = this;
      function wire(){
        var A = store['ctrl'], B = store['str'];
        if(!A || !B){ setTimeout(wire,120); return; }
        var netA = getNet(A), netB = getNet(B);
        if(!netA || !netB){ setTimeout(wire,120); return; }
        store._lock = false;
        function syncView(from, to){
          if(store._lock) return; store._lock = true;
          try{ to.moveTo({ position: from.getViewPosition(), scale: from.getScale(), animation:false }); }
          finally { store._lock = false; }
        }
        function syncAllNodes(from, to){
          if(store._lock) return; store._lock = true;
          try{
            var pos = from.getPositions(); var ids = Object.keys(pos);
            for (var i=0;i<ids.length;i++){ var p = pos[ids[i]]; if(p) to.moveNode(ids[i], p.x, p.y); }
          } finally { store._lock = false; }
        }
        setTimeout(function(){
          if(role==='str'){ syncAllNodes(netB, netA); syncView(netB, netA); }
        },300);
        function wireDrag(src, dst){
          src.on('dragEnd', function(){ syncAllNodes(src, dst); syncView(src, dst); });
          src.on('zoom',    function(){ syncView(src, dst); });
        }
        wireDrag(netA, netB); wireDrag(netB, netA);
      }
      setTimeout(wire,0);
    }
  ", jsonlite::toJSON(pair_id), jsonlite::toJSON(role))
  htmlwidgets::onRender(widget, js)
}
.tf_on_top <- function(widget){
  js <- "
  function(el, x){
    function getNet(i){ return (i && i.network) ? i.network :
                        (i && i.body && i.body.network ? i.body.network : null); }
    var net = getNet(this); if(!net) return;
    function sortTFTop(){
      try{
        var body = net.body;
        if(!body || !body.nodeIndices) return;
        var ids = body.nodeIndices.slice();
        ids.sort(function(a,b){
          var A = body.nodes[a], B = body.nodes[b];
          var aTF = A && A.options && (A.options.shape === 'box' || A.options.group === 'TF');
          var bTF = B && B.options && (B.options.shape === 'box' || B.options.group === 'TF');
          return (aTF === bTF) ? 0 : (aTF ? 1 : -1);
        });
        body.nodeIndices = ids;
      } catch(e){}
    }
    net.on('initRedraw',  sortTFTop);
    net.on('afterDrawing', sortTFTop);
    setTimeout(sortTFTop, 0);
  }"
  htmlwidgets::onRender(widget, js)
}
