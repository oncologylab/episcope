.gene_symbol_orgdb_package <- function(species = NULL, ref_genome = NULL) {
  species <- if (is.null(species) || !length(species) || is.na(species[[1L]]) || !nzchar(as.character(species)[[1L]])) {
    if (!is.null(ref_genome) && length(ref_genome) && !is.na(ref_genome[[1L]]) && nzchar(as.character(ref_genome)[[1L]])) {
      .pathway_species_from_ref_genome(ref_genome)
    } else {
      "human"
    }
  } else {
    species
  }
  species <- .normalize_pathway_species(species)
  if (identical(species, "mouse")) {
    "org.Mm.eg.db"
  } else {
    "org.Hs.eg.db"
  }
}

.gene_symbol_heuristic_case <- function(x, species = c("human", "mouse")) {
  species <- match.arg(species)
  x <- trimws(as.character(x))
  out <- x
  valid <- !is.na(out) & nzchar(out) & grepl("^[A-Za-z][A-Za-z0-9_.-]*$", out)
  if (identical(species, "human")) {
    out[valid] <- toupper(out[valid])
  } else {
    upper_like <- valid & grepl("^[A-Z0-9_.-]+$", out)
    out[upper_like] <- paste0(substr(out[upper_like], 1L, 1L), tolower(substr(out[upper_like], 2L, nchar(out[upper_like]))))
  }
  out[!valid] <- NA_character_
  out
}

.gene_symbol_cache <- new.env(parent = emptyenv())

.gene_symbol_load_maps <- function(species = NULL, ref_genome = NULL) {
  pkg <- .gene_symbol_orgdb_package(species = species, ref_genome = ref_genome)
  cache_key <- pkg
  if (exists(cache_key, envir = .gene_symbol_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gene_symbol_cache, inherits = FALSE))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace(pkg, quietly = TRUE)) {
    return(NULL)
  }
  db <- getExportedValue(pkg, pkg)
  keytypes <- AnnotationDbi::keytypes(db)
  symbol_keys <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  symbol_keys <- unique(as.character(symbol_keys[!is.na(symbol_keys) & nzchar(symbol_keys)]))
  symbol_exact <- stats::setNames(symbol_keys, symbol_keys)

  case_dt <- data.table::data.table(
    match_key = toupper(symbol_keys),
    canonical = symbol_keys
  )
  case_dt <- case_dt[, .(
    canonical = if (data.table::uniqueN(canonical) == 1L) canonical[[1L]] else NA_character_,
    ambiguous = data.table::uniqueN(canonical) > 1L
  ), by = match_key]
  symbol_case <- stats::setNames(case_dt$canonical, case_dt$match_key)
  symbol_case_ambiguous <- stats::setNames(case_dt$ambiguous, case_dt$match_key)

  load_select_map <- function(keytype, match_type) {
    if (!keytype %in% keytypes) return(data.table::data.table())
    keys <- AnnotationDbi::keys(db, keytype = keytype)
    keys <- unique(as.character(keys[!is.na(keys) & nzchar(keys)]))
    if (!length(keys)) return(data.table::data.table())
    sel <- suppressMessages(AnnotationDbi::select(
      db,
      keys = keys,
      keytype = keytype,
      columns = "SYMBOL"
    ))
    sel <- data.table::as.data.table(sel)
    if (!all(c(keytype, "SYMBOL") %in% names(sel))) return(data.table::data.table())
    sel <- sel[!is.na(get(keytype)) & nzchar(get(keytype)) & !is.na(SYMBOL) & nzchar(SYMBOL)]
    if (!nrow(sel)) return(data.table::data.table())
    sel[, match_key := toupper(as.character(get(keytype)))]
    sel <- unique(sel[, .(match_key, canonical = as.character(SYMBOL))])
    out <- sel[, .(
      canonical = if (data.table::uniqueN(canonical) == 1L) canonical[[1L]] else NA_character_,
      ambiguous = data.table::uniqueN(canonical) > 1L
    ), by = match_key]
    out[, match_type := match_type]
    out
  }

  alias_dt <- load_select_map("ALIAS", "alias")
  ensembl_dt <- load_select_map("ENSEMBL", "ensembl")
  entrez_dt <- load_select_map("ENTREZID", "entrez")

  maps <- list(
    package = pkg,
    species = if (identical(pkg, "org.Mm.eg.db")) "mouse" else "human",
    symbol_exact = symbol_exact,
    symbol_case = symbol_case,
    symbol_case_ambiguous = symbol_case_ambiguous,
    alias = alias_dt,
    ensembl = ensembl_dt,
    entrez = entrez_dt
  )
  assign(cache_key, maps, envir = .gene_symbol_cache)
  maps
}

#' Resolve gene symbols to formal species-specific symbols
#'
#' Resolves human or mouse gene names using Bioconductor organism annotation
#' packages when available, with a deterministic case-normalized fallback when
#' formal annotation is unavailable or a symbol has no database match.
#'
#' @param genes Character vector of gene names or identifiers.
#' @param species Species name, either `"human"` or `"mouse"`.
#' @param ref_genome Optional reference genome such as `"hg38"` or `"mm10"`.
#'   Used to infer `species` when `species` is `NULL`.
#' @param use_annotation Logical; use AnnotationDbi and organism annotation
#'   packages when available.
#' @param allow_heuristic Logical; use case-normalized fallback for unmatched
#'   symbols or when annotation packages are unavailable.
#'
#' @return A data frame with one row per input gene and columns describing the
#'   canonical symbol, matching key, match source, and ambiguity status.
#' @export
resolve_gene_symbols <- function(genes,
                                 species = NULL,
                                 ref_genome = NULL,
                                 use_annotation = TRUE,
                                 allow_heuristic = TRUE) {
  .assert_pkg("data.table")
  species <- if (is.null(species) || !length(species) || is.na(species[[1L]]) || !nzchar(as.character(species)[[1L]])) {
    if (!is.null(ref_genome) && length(ref_genome) && !is.na(ref_genome[[1L]]) && nzchar(as.character(ref_genome)[[1L]])) {
      .pathway_species_from_ref_genome(ref_genome)
    } else {
      "human"
    }
  } else {
    species
  }
  species <- .normalize_pathway_species(species)
  input <- as.character(genes)
  query_vec <- trimws(input)
  query_vec[is.na(input)] <- NA_character_
  out <- data.table::data.table(
    input_gene = input,
    query = query_vec,
    canonical_symbol = NA_character_,
    match_key = NA_character_,
    match_type = "unmatched",
    matched = FALSE,
    ambiguous = FALSE,
    annotation_source = NA_character_
  )
  valid <- !is.na(query_vec) & nzchar(query_vec)
  maps <- if (isTRUE(use_annotation)) .gene_symbol_load_maps(species = species, ref_genome = ref_genome) else NULL

  if (!is.null(maps) && any(valid)) {
    exact_idx <- valid & query_vec %in% names(maps$symbol_exact)
    if (any(exact_idx)) {
      out[exact_idx, `:=`(
        canonical_symbol = unname(maps$symbol_exact[query_vec[exact_idx]]),
        match_type = "official_exact",
        matched = TRUE,
        annotation_source = maps$package
      )]
    }

    unresolved <- valid & !out$matched
    case_key <- toupper(query_vec)
    case_idx <- unresolved & case_key %in% names(maps$symbol_case)
    if (any(case_idx)) {
      canon <- unname(maps$symbol_case[case_key[case_idx]])
      amb <- unname(maps$symbol_case_ambiguous[case_key[case_idx]])
      rows <- which(case_idx)
      out[rows, `:=`(
        canonical_symbol = canon,
        match_type = ifelse(amb, "ambiguous_official_case", "official_case"),
        matched = !amb & !is.na(canon),
        ambiguous = amb,
        annotation_source = maps$package
      )]
    }

    apply_map <- function(map_dt, type_label) {
      unresolved_inner <- valid & !out$matched & !out$ambiguous
      if (!nrow(map_dt) || !any(unresolved_inner)) return(invisible(NULL))
      idx <- match(case_key[unresolved_inner], map_dt$match_key)
      hit <- !is.na(idx)
      if (!any(hit)) return(invisible(NULL))
      rows <- which(unresolved_inner)[hit]
      map_rows <- map_dt[idx[hit]]
      out[rows, `:=`(
        canonical_symbol = map_rows$canonical,
        match_type = ifelse(map_rows$ambiguous, paste0("ambiguous_", type_label), type_label),
        matched = !map_rows$ambiguous & !is.na(map_rows$canonical),
        ambiguous = map_rows$ambiguous,
        annotation_source = maps$package
      )]
      invisible(NULL)
    }
    apply_map(maps$alias, "alias")
    apply_map(maps$ensembl, "ensembl")
    apply_map(maps$entrez, "entrez")
  }

  fallback_idx <- valid & !out$matched & !out$ambiguous & isTRUE(allow_heuristic)
  if (any(fallback_idx)) {
    fallback <- .gene_symbol_heuristic_case(query_vec[fallback_idx], species = species)
    ok <- !is.na(fallback) & nzchar(fallback)
    rows <- which(fallback_idx)
    out[rows[ok], `:=`(
      canonical_symbol = fallback[ok],
      match_type = "heuristic_case",
      matched = TRUE,
      annotation_source = if (is.null(maps)) "heuristic_no_orgdb" else "heuristic_unmatched"
    )]
  }
  out[matched == TRUE, match_key := toupper(canonical_symbol)]
  out[, query := NULL]
  out[]
}

.gene_symbol_key_table <- function(genes,
                                   species = NULL,
                                   ref_genome = NULL,
                                   use_annotation = TRUE,
                                   allow_heuristic = TRUE,
                                   input_col = "gene") {
  resolved <- resolve_gene_symbols(
    genes,
    species = species,
    ref_genome = ref_genome,
    use_annotation = use_annotation,
    allow_heuristic = allow_heuristic
  )
  out <- data.table::data.table(
    gene = resolved$input_gene,
    gene_key__ = resolved$match_key,
    gene_canonical = resolved$canonical_symbol,
    gene_match_type = resolved$match_type,
    gene_matched = resolved$matched,
    gene_ambiguous = resolved$ambiguous,
    gene_annotation_source = resolved$annotation_source
  )
  if (!identical(input_col, "gene")) data.table::setnames(out, "gene", input_col)
  out[]
}

.gene_symbol_match_summary <- function(match_type) {
  x <- as.character(match_type)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  tab <- sort(table(x), decreasing = TRUE)
  paste(paste0(names(tab), "=", as.integer(tab)), collapse = ";")
}

.write_gene_symbol_conversion_audit <- function(out_dir, tables) {
  .assert_pkg("data.table")
  if (is.null(out_dir) || !nzchar(as.character(out_dir)[[1L]]) || is.null(tables) || !length(tables)) {
    return(invisible(NULL))
  }
  details <- data.table::rbindlist(lapply(names(tables), function(source_name) {
    dt <- data.table::as.data.table(tables[[source_name]])
    if (!nrow(dt)) return(data.table::data.table())
    need <- c("gene", "gene_canonical", "gene_match_type", "gene_matched", "gene_ambiguous")
    if (!all(need %in% names(dt))) return(data.table::data.table())
    out <- unique(dt[, ..need])
    out[, source := source_name]
    out
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(details)) return(invisible(NULL))
  data.table::setcolorder(details, c("source", setdiff(names(details), "source")))
  details[, gene_matched := .as_logical_flag(gene_matched)]
  details[, gene_ambiguous := .as_logical_flag(gene_ambiguous)]
  details[, is_formal_match := gene_matched & !grepl("^heuristic", as.character(gene_match_type))]
  details[, is_heuristic_match := gene_matched & grepl("^heuristic", as.character(gene_match_type))]
  audit <- details[, .(
    total = .N,
    matched = sum(gene_matched, na.rm = TRUE),
    matched_percent = 100 * mean(gene_matched, na.rm = TRUE),
    formal_matched = sum(is_formal_match, na.rm = TRUE),
    formal_percent = 100 * mean(is_formal_match, na.rm = TRUE),
    heuristic_matched = sum(is_heuristic_match, na.rm = TRUE),
    ambiguous = sum(gene_ambiguous, na.rm = TRUE),
    unmatched = sum(!gene_matched & !gene_ambiguous, na.rm = TRUE),
    match_summary = .gene_symbol_match_summary(gene_match_type)
  ), by = source]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(audit, file.path(out_dir, "gene_symbol_conversion_audit.csv"))
  data.table::fwrite(details, file.path(out_dir, "gene_symbol_conversion_details.csv"))
  invisible(audit)
}
