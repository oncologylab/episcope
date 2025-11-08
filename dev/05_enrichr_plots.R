#' Paginated Enrichr dotplots (solid dots, nutrient border, y-beeswarm repel)
#'
#' - x = -log10(FDR) (capped at 20; fixed axis 0–20; dots stay inside)
#' - one dot per nutrient condition per term (Mixed expanded)
#' - **vertical beeswarm repel** inside each y row so identical p-values don't look different
#' - border color = nutrient (your palette); inner fill = Count (optional)
#'
#' @param fill_by_count logical; if TRUE (default) inner fill encodes #Genes (Count),
#'   else fill is neutral and Count is size only.
save_enrichr_faceted_pages <- function(
    x,
    fdr_col = "Adjusted.P.value",
    extra_non_tf = 5L,
    case_insensitive = TRUE,
    wrap_width = 48L,
    drop_go_ids = TRUE,
    exclude_databases = c("ARCHS4_TFs_Coexp"),
    ncol = 4L,
    nrow = 4L,
    height_per_facet = 3.0,
    width = 16,
    fill_by_count = TRUE,
    out_pdf = "enrichr_tf_facets_paginated.pdf"
) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

  # ---------- load ----------
  df <- x
  if (is.character(x) && length(x) == 1L) {
    ext <- tolower(tools::file_ext(x))
    if (ext %in% c("tsv","txt"))       df <- readr::read_tsv(x, show_col_types = FALSE)
    else if (ext == "csv")             df <- readr::read_csv(x, show_col_types = FALSE)
    else if (ext %in% c("xlsx","xls")) df <- readxl::read_excel(x)
    else cli::cli_abort("Unsupported file type: {ext}")
  }
  if (!is.data.frame(df)) cli::cli_abort("`x` must be a data.frame or a readable file path.")

  nm <- names(df)
  nm <- gsub("\\s+", ".", nm)
  nm <- sub("^Old\\.", "Old_", nm)
  nm <- sub("^Adjusted\\.P\\.value$", "Adjusted.P.value", nm)
  nm <- sub("^P\\.value$", "P.value", nm)
  names(df) <- nm

  need <- c("Term","Overlap","Genes","regulation", fdr_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) cli::cli_abort("Missing columns: {paste(miss, collapse=', ')}")

  if ("database" %in% names(df) && length(exclude_databases)) {
    df <- df[!(df$database %in% exclude_databases), , drop = FALSE]
  }

  # ---------- parse ----------
  parts <- strsplit(as.character(df$Overlap), "/", fixed = TRUE)
  df$k <- vapply(parts, function(y) suppressWarnings(as.numeric(y[1L])), numeric(1))
  df$n <- vapply(parts, function(y) suppressWarnings(as.numeric(y[2L])), numeric(1))
  df$gene_count <- ifelse(is.na(df$Genes) | df$Genes == "",
                          0L,
                          vapply(strsplit(df$Genes, ";", fixed = TRUE), length, integer(1)))

  eps <- 1e-300
  df$FDR       <- as.numeric(df[[fdr_col]])
  df$mlog10FDR <- -log10(pmax(df$FDR, eps))
  df$tf        <- sub("_.*$", "", as.character(df$regulation))

  df$highlight <- if (case_insensitive) {
    mapply(function(term, tf) grepl(tf, term, ignore.case = TRUE), df$Term, df$tf)
  } else {
    mapply(function(term, tf) grepl(tf, term, fixed = TRUE), df$Term, df$tf)
  }

  lbl <- df$Term
  if (drop_go_ids) lbl <- gsub("\\s*\\(GO:\\d+\\)$", "", lbl)
  df$Term_label <- stringr::str_wrap(lbl, width = wrap_width)

  df <- df[order(df$regulation, df$FDR, -df$gene_count, df$Term_label), , drop = FALSE]
  keep_list <- lapply(split(df, df$regulation), function(dd) {
    hi  <- dd[dd$highlight, , drop = FALSE]
    non <- dd[!dd$highlight, , drop = FALSE]
    if (nrow(non)) non <- non[order(non$FDR, -non$gene_count, non$Term_label), , drop = FALSE]
    non <- utils::head(non, extra_non_tf)
    out <- rbind(hi, non)
    out[order(out$FDR, -out$gene_count, out$Term_label), , drop = FALSE]
  })
  d_keep <- dplyr::bind_rows(keep_list)

  # ---------- explode conditions ----------
  cond_col <- if ("conditon" %in% names(d_keep)) "conditon" else if ("condition" %in% names(d_keep)) "condition" else NA_character_
  nutrient_levels <- c("BCAA","FBS","Glc","Gln.Arg","Lys","Met.Cys","Trp","Mixed","Other")

  if (!is.na(cond_col)) {
    d_long <- tidyr::separate_rows(d_keep, !!rlang::sym(cond_col), sep = ";\\s*")
    d_long$nutrient <- vapply(d_long[[cond_col]], function(s) {
      m <- regmatches(s, regexpr("BCAA|FBS|Glc|Gln\\.Arg|Lys|Met\\.Cys|Trp", s, ignore.case = TRUE, perl = TRUE))
      if (length(m) == 1L) {
        m <- gsub("\\s+", "", m)
        m <- gsub("(?i)glc","Glc", m, perl=TRUE)
        m <- gsub("(?i)fbs","FBS", m, perl=TRUE)
        m <- gsub("(?i)bcaa","BCAA", m, perl=TRUE)
        m <- gsub("(?i)gln\\.arg","Gln.Arg", m, perl=TRUE)
        m <- gsub("(?i)lys","Lys", m, perl=TRUE)
        m <- gsub("(?i)met\\.cys","Met.Cys", m, perl=TRUE)
        m <- gsub("(?i)trp","Trp", m, perl=TRUE)
        m
      } else if (!is.na(s) && grepl(";", s)) {
        "Mixed"
      } else {
        "Other"
      }
    }, character(1))
  } else {
    d_long <- d_keep
    d_long$nutrient <- "Other"
  }
  d_long$nutrient <- factor(d_long$nutrient, levels = nutrient_levels)

  # cap x; NO horizontal offsets (we'll repel vertically)
  d_long$mlog10FDR_cap <- pmin(pmax(d_long$mlog10FDR, 0), 20)

  # facet-wise y order (top = most significant)
  d_long <- dplyr::group_by(d_long, regulation)
  d_long <- dplyr::mutate(d_long, Term_label = factor(Term_label, levels = rev(unique(Term_label))))
  d_long <- dplyr::ungroup(d_long)

  # pagination
  per_page <- as.integer(ncol) * as.integer(nrow)
  regs  <- unique(d_long$regulation)
  pages <- split(regs, ceiling(seq_along(regs) / per_page))
  page_h <- max(4, as.integer(nrow) * height_per_facet)

  # palettes
  nutrient_cols <- c(
    BCAA   = "#3a78af",
    FBS    = "#ef812f",
    Glc    = "#308a4e",
    `Gln.Arg` = "#414b8c",
    Lys    = "#e89c84",
    `Met.Cys` = "#ca2f2d",
    Trp    = "#916ab6",
    Mixed  = "#8a8a8a",
    Other  = "#666666"
  )
  plasma_cols <- viridisLite::plasma(11, direction = -1)

  if (!is.null(out_pdf)) {
    grDevices::pdf(out_pdf, width = width, height = page_h, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  res <- vector("list", length(pages))
  for (i in seq_along(pages)) {
    regs_i <- pages[[i]]
    dd <- d_long[d_long$regulation %in% regs_i, , drop = FALSE]

    aes_common <- ggplot2::aes(
      x = mlog10FDR_cap, y = Term_label,
      size = gene_count,
      color = nutrient
    )
    if (fill_by_count) {
      aes_common$fill <- dd$gene_count
    } else {
      aes_common$fill <- I("#ffffff")  # neutral inner
    }

    p <- ggplot2::ggplot(dd, aes_common) +
      # vertical beeswarm repel inside each y row (keeps x exact)
      ggbeeswarm::geom_beeswarm(
        shape = 21, stroke = 1.0, cex = NA, alpha = 0.95,
        priority = "density", groupOnX = FALSE, dodge.width = 0.5
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 20),
        breaks = seq(0, 20, by = 5),
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::scale_size_area(name = "Count", max_size = 9, breaks = c(10,20,40,60)) +
      ggplot2::scale_color_manual(name = "Condition", values = nutrient_cols, drop = FALSE)
    if (fill_by_count) {
      p <- p + ggplot2::scale_fill_gradientn(name = "Count", colours = plasma_cols)
    } else {
      p <- p + ggplot2::guides(fill = "none")
    }
    p <- p +
      ggplot2::facet_wrap(~ regulation, ncol = ncol, scales = "free_y") +
      ggplot2::labs(
        title = "Pathway terms per TF regulation set (filtered & wrapped)",
        subtitle = "x = −log10(FDR) (capped at 20); solid dots; vertical beeswarm repel; border = nutrient",
        x = expression(-log[10]~FDR), y = "Term"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(linewidth = 0.2),
        panel.grid.major.y = ggplot2::element_line(linewidth = 0.2),
        panel.spacing = grid::unit(1.0, "lines"),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(size = 9, margin = ggplot2::margin(r = 4)),
        axis.title.x = ggplot2::element_text(face = "bold"),
        axis.title.y = ggplot2::element_text(face = "bold"),
        strip.clip = "off",
        plot.margin = grid::unit(c(0.5, 0.8, 0.5, 1.2), "lines"),
        legend.key.height = grid::unit(0.6, "lines")
      )

    res[[i]] <- p
    if (!is.null(out_pdf)) print(p)
  }

  if (!is.null(out_pdf)) cli::cli_inform(c(v = paste0("Saved PDF: ", normalizePath(out_pdf))))
  invisible(list(pages = res, data = d_long))
}


df <- readxl::read_excel(file.path(
  "inst","extdata","EnrichR_top_pathways_regulated_by_Top_TFs_GeneEnhancer_regulated_genes_2_delta_link_1_unfiltered.xlsx"))

save_enrichr_faceted_pages(
  df,
  ncol = 3,     # 4 columns
  nrow = 4,     # 4 rows  → 16 facets per page
  extra_non_tf = 5,
  height_per_facet = 3.0,
  out_pdf = "enrichr_tf_facets_paginated.pdf"
)
# End of file
