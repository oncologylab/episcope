# Gene TSS annotation helpers.

.gene_tss_normalize_ref_genome <- function(ref_genome = NULL) {
  ref <- if (is.null(ref_genome) || !length(ref_genome) || is.na(ref_genome[[1L]])) {
    .cfg_get("ref_genome")
  } else {
    ref_genome[[1L]]
  }
  ref <- tolower(as.character(ref))
  if (!length(ref) || is.na(ref[[1L]]) || !nzchar(ref[[1L]])) ref <- "hg38"
  aliases <- c(grch38 = "hg38", hg19 = "hg19", grcm38 = "mm10", mm10 = "mm10", hg38 = "hg38")
  if (ref %in% names(aliases)) return(unname(aliases[[ref]]))
  ref
}

.gene_tss_chr_style <- function(chr) {
  chr <- as.character(chr)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr[chr %in% c("MT", "Mt", "mt", "M")] <- "M"
  paste0("chr", chr)
}

.gene_tss_read_table <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    .log_abort("gene_tss path must be a single non-empty file path.")
  }
  if (!file.exists(path)) .log_abort("gene_tss file not found: {path}")
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "rds")) {
    out <- readRDS(path)
  } else if (identical(ext, "parquet")) {
    if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet gene_tss files.")
    out <- arrow::read_parquet(path)
  } else {
    out <- readr::read_csv(path, show_col_types = FALSE)
  }
  if (!is.data.frame(out)) .log_abort("gene_tss file must contain a data.frame.")
  out
}

.gene_tss_from_ensdb <- function(ref_genome, genes = NULL) {
  ref <- .gene_tss_normalize_ref_genome(ref_genome)
  ensdb_pkg <- switch(
    ref,
    hg38 = "EnsDb.Hsapiens.v86",
    mm10 = "EnsDb.Mmusculus.v79",
    .log_abort("Automatic gene TSS loading is only configured for ref_genome hg38 and mm10. Provide a gene_tss table for ref_genome={ref}.")
  )
  if (!requireNamespace(ensdb_pkg, quietly = TRUE)) {
    .log_abort("Package {ensdb_pkg} is required for automatic {ref} gene TSS loading. Install it or provide gene_tss in project.yaml.")
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    .log_abort("Package AnnotationDbi is required for automatic gene TSS loading.")
  }
  ensdb <- get(ensdb_pkg, envir = asNamespace(ensdb_pkg))
  keys <- if (is.null(genes)) AnnotationDbi::keys(ensdb, keytype = "SYMBOL") else unique(as.character(genes))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  ann <- AnnotationDbi::select(
    ensdb,
    keys = keys,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND")
  )
  ann <- ann[!is.na(ann$SYMBOL) & !is.na(ann$SEQNAME) & !is.na(ann$GENESEQSTART) & !is.na(ann$GENESEQEND), , drop = FALSE]
  if (!nrow(ann)) .log_abort("No {ref} gene TSS records were resolved from {ensdb_pkg}.")
  strand_int <- suppressWarnings(as.integer(ann$SEQSTRAND))
  tss <- ifelse(is.finite(strand_int) & strand_int < 0L, as.integer(ann$GENESEQEND), as.integer(ann$GENESEQSTART))
  strand <- ifelse(is.finite(strand_int) & strand_int < 0L, "-", "+")
  out <- data.frame(
    target_gene = as.character(ann$SYMBOL),
    target_chr = .gene_tss_chr_style(ann$SEQNAME),
    target_tss = as.integer(tss),
    target_strand = strand,
    stringsAsFactors = FALSE
  )
  allowed <- paste0("chr", c(1:22, "X", "Y", "M"))
  if (identical(ref, "mm10")) allowed <- paste0("chr", c(1:19, "X", "Y", "M"))
  out <- out[out$target_chr %in% allowed & !duplicated(out$target_gene), , drop = FALSE]
  tibble::as_tibble(out)
}

#' Load gene TSS annotations
#'
#' Loads a normalized gene TSS table for Module 2. If `gene_tss` is supplied,
#' it is read and normalized. Otherwise the table is built from an installed
#' EnsDb annotation matching `ref_genome`.
#'
#' @param gene_tss Optional data frame or path to a CSV, CSV.GZ, Parquet, or RDS
#'   gene TSS table.
#' @param ref_genome Reference genome. Supported automatic values are `hg38`
#'   and `mm10`.
#' @param genes Optional gene symbols to restrict automatic EnsDb loading.
#' @param verbose Emit concise progress messages.
#' @return A tibble with `target_gene`, `target_chr`, `target_tss`, and
#'   `target_strand`.
#' @keywords internal
#' @export
load_gene_tss <- function(gene_tss = NULL, ref_genome = NULL, genes = NULL, verbose = TRUE) {
  ref <- .gene_tss_normalize_ref_genome(ref_genome)
  if (is.null(gene_tss) || (is.character(gene_tss) && length(gene_tss) == 1L && !nzchar(gene_tss))) {
    if (isTRUE(verbose)) .log_inform("Loading gene TSS annotations for ref_genome={ref} from installed EnsDb.")
    out <- .gene_tss_from_ensdb(ref, genes = genes)
  } else if (is.character(gene_tss) && length(gene_tss) == 1L) {
    if (isTRUE(verbose)) .log_inform("Loading gene TSS annotations from {gene_tss}.")
    out <- .gene_tss_read_table(gene_tss)
  } else if (is.data.frame(gene_tss)) {
    out <- gene_tss
  } else {
    .log_abort("gene_tss must be NULL, a data.frame, or a single file path.")
  }
  out <- .module2_normalize_gene_tss(out)
  out$target_chr <- .gene_tss_chr_style(out$target_chr)
  out <- out[!duplicated(out$target_gene), , drop = FALSE]
  tibble::as_tibble(out)
}

.module2_resolve_gene_tss <- function(gene_tss, cfg, multiomic_data = NULL, verbose = TRUE) {
  if (!is.null(gene_tss)) return(load_gene_tss(gene_tss, ref_genome = .module2_cfg_value(cfg, "ref_genome", NULL), verbose = FALSE))
  gene_tss_path <- .module2_cfg_value(cfg, "gene_tss", NULL)
  if (is.null(gene_tss_path)) gene_tss_path <- .module2_cfg_value(cfg, "gene_tss_path", NULL)
  genes <- NULL
  if (!is.null(multiomic_data) && is.list(multiomic_data$matrices) && is.matrix(multiomic_data$matrices$gene_expr)) {
    genes <- rownames(multiomic_data$matrices$gene_expr)
  }
  load_gene_tss(
    gene_tss = gene_tss_path,
    ref_genome = .module2_cfg_value(cfg, "ref_genome", NULL),
    genes = genes,
    verbose = verbose
  )
}
