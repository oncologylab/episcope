# FILE: R/utils_grn_diff.R
# ---------------------------------------------------------------------------
# General TF–gene links delta computation (episcope)
# Author: Yaoxiang Li
# Updated: 2025-11-10
#
# Summary
#   Load TF–gene link tables, harmonize columns, and compute per-link deltas
#   for link_score, fp_bed_score, tf_expr, gene_expr, active (cond1 - cond2).
#   Includes bulk/driver helpers for per-cell contrasts.
#
# Conventions
#   - roxygen2 docs; explicit namespacing
#   - cli::cli_abort()/cli_warn()/cli_inform()
#   - NO rlang::.data; prefer bare names or base subsetting
#   - Programmatic selection uses [[ ]] and base ops
#   - Deterministic ordering
# ---------------------------------------------------------------------------

# =============================
# Internal utilities (helpers)
# =============================

#' @keywords internal
#' @importFrom cli cli_inform
.log <- function(msg, verbose = TRUE) {
  if (isTRUE(verbose)) cli::cli_inform(c(v = paste0("[utils_grn_diff] ", msg)))
  invisible(NULL)
}

#' @keywords internal
.safe_label <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x, perl = TRUE)
}

#' @keywords internal
#' Derive a readable condition name from a file path or symbol
.derive_cond_name <- function(x, fallback = "cond") {
  nm <- tryCatch(deparse(substitute(x)), error = function(e) NA_character_)
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    b <- basename(x)
    nm <- sub("\\.csv$", "", b, ignore.case = TRUE)
  }
  if (!is.character(nm) || is.na(nm) || nm == "") nm <- fallback
  nm
}

#' @keywords internal
#' Assert identical key sets and uniqueness across two tables
#' @importFrom cli cli_abort
.assert_identical_link_keys <- function(t1, t2,
                                        key_cols = c("tf", "gene_key", "peak_id"),
                                        verbose = TRUE) {
  # per-table duplicates
  dup1 <- t1[duplicated(t1[key_cols]), key_cols, drop = FALSE]
  dup2 <- t2[duplicated(t2[key_cols]), key_cols, drop = FALSE]

  if (nrow(dup1)) {
    ex <- utils::head(unique(dup1), 5)
    cli::cli_abort(c(
      "Duplicate keys detected in cond1.",
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }
  if (nrow(dup2)) {
    ex <- utils::head(unique(dup2), 5)
    cli::cli_abort(c(
      "Duplicate keys detected in cond2.",
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }

  # identical key sets
  k1 <- unique(t1[key_cols])
  k2 <- unique(t2[key_cols])

  tag1 <- interaction(k1, drop = TRUE)
  tag2 <- interaction(k2, drop = TRUE)

  only_in_1 <- k1[!(tag1 %in% tag2), , drop = FALSE]
  only_in_2 <- k2[!(tag2 %in% tag1), , drop = FALSE]

  if (nrow(only_in_1)) {
    ex <- utils::head(only_in_1, 5)
    cli::cli_abort(c(
      "Some keys exist only in cond1.",
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }
  if (nrow(only_in_2)) {
    ex <- utils::head(only_in_2, 5)
    cli::cli_abort(c(
      "Some keys exist only in cond2.",
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }
  invisible(TRUE)
}

#' @keywords internal
.collapse_links_by_key <- function(tbl,
                                   key_cols = c("tf", "gene_key", "peak_id"),
                                   which_cond = "cond",
                                   method = c("max_link_score", "first"),
                                   verbose = TRUE) {
  method <- match.arg(method)
  n0 <- nrow(tbl)
  if (!n0) return(tbl)

  if (method == "first") {
    out <- tbl[!duplicated(tbl[key_cols]), , drop = FALSE]
  } else {
    # keep row with max link_score within key group (break ties by first)
    ord <- order(tbl[[key_cols[1]]], tbl[[key_cols[2]]], tbl[[key_cols[3]]],
                 -as.numeric(tbl[["link_score"]]), method = "radix", na.last = TRUE)
    tmp <- tbl[ord, , drop = FALSE]
    out <- tmp[!duplicated(tmp[key_cols]), , drop = FALSE]
    # restore stable order by keys
    o2 <- order(out[[key_cols[1]]], out[[key_cols[2]]], out[[key_cols[3]]],
                method = "radix", na.last = TRUE)
    out <- out[o2, , drop = FALSE]
  }

  n1 <- nrow(out)
  if (n1 < n0) .log(paste0("Collapsed ", (n0 - n1), " duplicate row(s) in ", which_cond,
                           " using method=", method, "."), verbose = verbose)
  out
}

#' @keywords internal
.reconcile_link_signs <- function(t1, t2,
                                  key_cols = c("tf", "gene_key", "peak_id"),
                                  policy = c("error", "prefer_cond1", "prefer_max_magnitude", "drop"),
                                  cond1_name = "cond1", cond2_name = "cond2",
                                  verbose = TRUE) {
  policy <- match.arg(policy)

  j <- merge(t1[c(key_cols, "link_sign", "link_score")],
             t2[c(key_cols, "link_sign", "link_score")],
             by = key_cols, all = FALSE, suffixes = c("_1","_2"))
  bad <- j[!is.na(j[["link_sign_1"]]) &
             !is.na(j[["link_sign_2"]]) &
             j[["link_sign_1"]] != j[["link_sign_2"]], , drop = FALSE]

  if (!nrow(bad)) {
    .log("Sanity: link_sign agrees across conditions for all non-NA keys.", verbose = verbose)
    return(list(t1 = t1, t2 = t2))
  }

  if (policy == "error") {
    ex <- utils::head(bad, 5)
    cli::cli_abort(c(
      paste0("link_sign mismatch between ", cond1_name, " and ", cond2_name, " for ", nrow(bad), " key(s)."),
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }

  if (policy == "drop") {
    tag_bad <- interaction(bad[key_cols], drop = TRUE)
    keep1 <- !(interaction(t1[key_cols], drop = TRUE) %in% tag_bad)
    keep2 <- !(interaction(t2[key_cols], drop = TRUE) %in% tag_bad)
    .log(paste0("Dropped ", sum(!keep1 | !keep2), " key(s) with link_sign mismatch."), verbose = verbose)
    return(list(t1 = t1[keep1, , drop = FALSE], t2 = t2[keep2, , drop = FALSE]))
  }

  if (policy == "prefer_cond1") {
    # Copy cond1's sign to cond2 for mismatched keys
    tag_bad <- interaction(bad[key_cols], drop = TRUE)
    tag_t2  <- interaction(t2[key_cols],  drop = TRUE)
    idx <- match(tag_t2, tag_bad)
    need <- which(!is.na(idx))
    if (length(need)) {
      t2$link_sign[need] <- bad$link_sign_1[idx[need]]
    }
    .log(paste0("Harmonized link_sign for ", nrow(bad), " key(s) to match ", cond1_name, "."), verbose = verbose)
    return(list(t1 = t1, t2 = t2))
  }

  # prefer_max_magnitude
  prefer2 <- abs(bad[["link_score_2"]]) > abs(bad[["link_score_1"]])

  # keys where cond2 wins -> set cond1 sign to cond2's
  if (any(prefer2)) {
    k12 <- bad[prefer2, c(key_cols, "link_sign_2"), drop = FALSE]
    colnames(k12)[ncol(k12)] <- "link_sign_ref"
    tag_t1 <- interaction(t1[key_cols], drop = TRUE)
    idx <- match(tag_t1, interaction(k12[key_cols], drop = TRUE))
    pos <- which(!is.na(idx))
    if (length(pos)) t1$link_sign[pos] <- k12$link_sign_ref[idx[pos]]
  }
  # others -> set cond2 sign to cond1's
  if (any(!prefer2)) {
    k21 <- bad[!prefer2, c(key_cols, "link_sign_1"), drop = FALSE]
    colnames(k21)[ncol(k21)] <- "link_sign_ref"
    tag_t2 <- interaction(t2[key_cols], drop = TRUE)
    idx <- match(tag_t2, interaction(k21[key_cols], drop = TRUE))
    pos <- which(!is.na(idx))
    if (length(pos)) t2$link_sign[pos] <- k21$link_sign_ref[idx[pos]]
  }

  .log(paste0("Harmonized link_sign for ", nrow(bad), " mismatched key(s) by larger magnitude."), verbose = verbose)
  list(t1 = t1, t2 = t2)
}

# ==============================
# Public API: Load + Harmonize
# ==============================

#' Load TF–gene links table from CSV or tibble
#'
#' Harmonizes common columns, ensures required IDs, and standardizes `link_sign`.
#' If `link_sign` is absent, it is derived from `r_gene` when available, else
#' from the sign of `link_score` (and `link_score` is abs()-ed accordingly).
#'
#' @param x A `data.frame`/`tibble` **or** a character file path to a CSV.
#' @param verbose Logical; emit progress.
#'
#' @return A tibble with standardized columns and types.
#' @export
#'
#' @importFrom readr read_csv
#' @importFrom tibble as_tibble
#' @importFrom janitor clean_names make_clean_names
#' @importFrom cli cli_warn
load_links_table <- function(x, verbose = TRUE) {
  .log("Loading links table…", verbose = verbose)

  tbl <- if (is.character(x) && length(x) == 1L) {
    .log(paste0("from CSV: ", x), verbose = verbose)
    readr::read_csv(x, show_col_types = FALSE)
  } else if (inherits(x, "data.frame")) {
    .log("from in-memory tibble/data.frame", verbose = verbose)
    x
  } else {
    cli::cli_abort("`x` must be a file path or a data.frame/tibble.")
  }

  tbl <- tibble::as_tibble(janitor::clean_names(tbl))

  # Required IDs
  req_ids <- c("tf", "gene_key", "peak_id")
  miss <- setdiff(req_ids, names(tbl))
  if (length(miss)) cli::cli_abort("Missing required columns: {paste(miss, collapse = ', ')}")

  # active_link -> active (if needed)
  if (!("active" %in% names(tbl)) && "active_link" %in% names(tbl)) {
    tbl[["active"]] <- tbl[["active_link"]]
    .log("Mapped active_link -> active", verbose = verbose)
  }
  if (!("active" %in% names(tbl))) tbl[["active"]] <- NA

  # Unify link_score
  if (!("link_score" %in% names(tbl))) {
    alts <- c("edge_weight", "link_weight", "weight")
    alt_hit <- alts[alts %in% names(tbl)]
    if (length(alt_hit)) {
      tbl[["link_score"]] <- as.numeric(tbl[[alt_hit[1L]]])
      .log(paste0("Mapped ", alt_hit[1L], " -> link_score"), verbose = verbose)
    } else {
      tbl[["link_score"]] <- NA_real_
      cli::cli_warn("No link_score/edge_weight/link_weight/weight found; created NA link_score.")
    }
  }

  # Map fp_score -> fp_bed_score if needed
  if (!("fp_bed_score" %in% names(tbl))) {
    alts <- c("fp_score", "footprint_score", "fp")
    alt_hit <- alts[alts %in% names(tbl)]
    if (length(alt_hit)) {
      tbl[["fp_bed_score"]] <- as.numeric(tbl[[alt_hit[1L]]])
      .log(paste0("Mapped ", alt_hit[1L], " -> fp_bed_score"), verbose = verbose)
    } else {
      tbl[["fp_bed_score"]] <- NA_real_
    }
  }

  # Ensure typical numeric columns exist
  for (cc in c("tf_expr", "gene_expr")) if (!(cc %in% names(tbl))) tbl[[cc]] <- NA_real_

  # Coerce numeric/logical for known fields
  num_cols <- intersect(c("link_score","fp_bed_score","tf_expr","gene_expr",
                          "r_gene","p_gene","p_adj_gene","r_tf","p_tf","p_adj_tf"),
                        names(tbl))
  for (cc in num_cols) tbl[[cc]] <- as.numeric(tbl[[cc]])
  tbl[["active"]] <- as.logical(tbl[["active"]])

  # Sanity logs
  ns <- tbl[["link_score"]]
  .log(paste0("Sanity: link_score present; pos=", sum(ns > 0, na.rm = TRUE),
              ", neg=", sum(ns < 0, na.rm = TRUE),
              ", zero=", sum(ns == 0, na.rm = TRUE),
              "; link_sign present: ", "link_sign" %in% names(tbl), "."),
       verbose = verbose)

  # Build/standardize link_sign
  if (!("link_sign" %in% names(tbl))) {
    if ("r_gene" %in% names(tbl)) {
      rg <- tbl[["r_gene"]]
      tbl[["link_sign"]] <- ifelse(is.na(rg), NA_character_, ifelse(rg < 0, "-", "+"))
      tbl[["link_score"]] <- abs(tbl[["link_score"]])
      .log("Created link_sign from r_gene sign; set link_score := abs(link_score).", verbose = verbose)
    } else {
      # fallback: from link_score sign
      pos <- sum(ns > 0, na.rm = TRUE); neg <- sum(ns < 0, na.rm = TRUE)
      if (pos > 0 && neg > 0) {
        tbl[["link_sign"]] <- ifelse(ns < 0, "-", ifelse(ns > 0, "+", "+"))
        tbl[["link_score"]] <- abs(ns)
        .log("Derived link_sign from link_score; set link_score := abs(link_score).", verbose = verbose)
      } else if (pos > 0 && neg == 0) {
        tbl[["link_sign"]] <- "+"
        .log("Derived link_sign='+' (link_score non-negative).", verbose = verbose)
      } else if (pos == 0 && neg > 0) {
        tbl[["link_sign"]] <- "-"
        tbl[["link_score"]] <- abs(ns)
        .log("Derived link_sign='-' (link_score non-positive); set link_score := abs(link_score).", verbose = verbose)
      } else {
        tbl[["link_sign"]] <- ifelse(is.na(ns), NA_character_, "+")
        .log("Derived link_sign: zeros -> '+', NAs -> NA.", verbose = verbose)
      }
    }
  } else {
    v <- as.character(tbl[["link_sign"]])
    vt <- trimws(v)
    vl <- tolower(vt)
    mapped <- ifelse(vt %in% c("+","-"), vt,
                     ifelse(vl %in% c("plus","positive","pos","p","up"), "+",
                            ifelse(vl %in% c("minus","negative","neg","n","down"), "-", NA_character_)))
    changed <- sum(!is.na(vt) & vt != mapped, na.rm = TRUE)
    invalid <- sum(is.na(mapped) & !is.na(vt))
    if (changed > 0) .log(paste0("Standardized ", changed, " link_sign value(s) to '+'/'-'."), verbose = verbose)
    if (invalid > 0) cli::cli_warn("{invalid} link_sign value(s) unrecognized; set to NA.")
    tbl[["link_sign"]] <- mapped
  }

  # Keep canonical columns (plus common stats if present)
  keep <- c("tf","gene_key","peak_id","link_score","fp_bed_score","tf_expr","gene_expr","active",
            "link_sign","r_gene","p_gene","p_adj_gene","r_tf","p_tf","p_adj_tf")
  keep <- intersect(keep, names(tbl))
  tbl[keep]
}

# ==============================
# Public API: Compare (Two / Bulk)
# ==============================

#' Compare two conditions (cond1 - cond2) for TF–gene links
#'
#' Computes deltas for `link_score`, `fp_bed_score`, `tf_expr`, `gene_expr`,
#' plus log2FC for expressions with a pseudocount. Enforces identical key sets
#' (`tf,gene_key,peak_id`) when `strict = TRUE` and reconciles `link_sign`
#' according to a chosen policy.
#'
#' @param cond1,cond2 A tibble/data.frame or CSV path.
#' @param cond1_name,cond2_name Optional condition names; default derived from inputs.
#' @param clean_names Logical; if TRUE, sanitize condition names via `janitor::make_clean_names()`.
#' @param strict Logical; if TRUE, enforce identical key sets.
#' @param dedupe Logical; if TRUE, collapse duplicate keys per condition.
#' @param dedupe_method One of `c("max_link_score","first")`.
#' @param sign_policy One of `c("error","prefer_cond1","prefer_max_magnitude","drop")`.
#' @param pseudocount Non-negative numeric added for log2FC calculations (default 1).
#' @param verbose Logical.
#'
#' @return A tibble with IDs, per-condition values, deltas, and log2FC columns.
#' @export
#'
#' @importFrom janitor make_clean_names
# BEGIN EDIT: replace the entire compare_links_two_conditions() with this version
compare_links_two_conditions <- function(cond1, cond2,
                                         cond1_name = NULL,
                                         cond2_name = NULL,
                                         clean_names = FALSE,
                                         strict = TRUE,
                                         dedupe = TRUE,
                                         dedupe_method = c("max_link_score","first"),
                                         sign_policy = c("error","prefer_cond1","prefer_max_magnitude","drop"),
                                         pseudocount = 1,
                                         verbose = TRUE) {
  dedupe_method <- match.arg(dedupe_method)
  sign_policy   <- match.arg(sign_policy)

  pc <- as.numeric(pseudocount)
  if (!is.finite(pc) || pc < 0) cli::cli_abort("`pseudocount` must be a non-negative finite number.")

  t1 <- load_links_table(cond1, verbose = verbose)
  t2 <- load_links_table(cond2, verbose = verbose)

  if (is.null(cond1_name)) cond1_name <- .derive_cond_name(cond1, "cond1")
  if (is.null(cond2_name)) cond2_name <- .derive_cond_name(cond2, "cond2")
  if (isTRUE(clean_names)) {
    cond1_name <- janitor::make_clean_names(cond1_name)
    cond2_name <- janitor::make_clean_names(cond2_name)
  }
  .log(paste0("Comparing ", cond1_name, " vs ", cond2_name, " (delta = cond1 - cond2)"), verbose = verbose)
  .log(paste0("Using pseudocount=", pc, " for log2FC_* columns."), verbose = verbose)

  key_cols <- c("tf","gene_key","peak_id")

  # Deduplicate by key
  if (isTRUE(dedupe)) {
    t1 <- .collapse_links_by_key(t1, key_cols, which_cond = cond1_name, method = dedupe_method, verbose = verbose)
    t2 <- .collapse_links_by_key(t2, key_cols, which_cond = cond2_name, method = dedupe_method, verbose = verbose)
  }

  if (isTRUE(strict)) .assert_identical_link_keys(t1, t2, key_cols = key_cols, verbose = verbose)

  # Reconcile link_sign across conditions
  harmon <- .reconcile_link_signs(t1, t2, key_cols = key_cols,
                                  policy = sign_policy,
                                  cond1_name = cond1_name, cond2_name = cond2_name,
                                  verbose = verbose)
  t1 <- harmon$t1; t2 <- harmon$t2

  # Deterministic ordering
  ord1 <- order(t1$tf, t1$gene_key, t1$peak_id, method = "radix")
  ord2 <- order(t2$tf, t2$gene_key, t2$peak_id, method = "radix")
  t1 <- t1[ord1, , drop = FALSE]
  t2 <- t2[ord2, , drop = FALSE]

  s1 <- paste0("_", cond1_name)
  s2 <- paste0("_", cond2_name)

  # --- Build per-condition frames with explicit dynamic names (robust) ---
  t1s <- data.frame(tf = t1$tf, gene_key = t1$gene_key, peak_id = t1$peak_id, check.names = FALSE)
  t1s[[paste0("link_score",   s1)]] <- t1$link_score
  t1s[[paste0("link_sign",    s1)]] <- t1$link_sign
  t1s[[paste0("fp_bed_score", s1)]] <- t1$fp_bed_score
  t1s[[paste0("tf_expr",      s1)]] <- t1$tf_expr
  t1s[[paste0("gene_expr",    s1)]] <- t1$gene_expr
  t1s[[paste0("active",       s1)]] <- t1$active
  t1s[[paste0("r_tf",         s1)]] <- t1$r_tf

  t2s <- data.frame(tf = t2$tf, gene_key = t2$gene_key, peak_id = t2$peak_id, check.names = FALSE)
  t2s[[paste0("link_score",   s2)]] <- t2$link_score
  t2s[[paste0("link_sign",    s2)]] <- t2$link_sign
  t2s[[paste0("fp_bed_score", s2)]] <- t2$fp_bed_score
  t2s[[paste0("tf_expr",      s2)]] <- t2$tf_expr
  t2s[[paste0("gene_expr",    s2)]] <- t2$gene_expr
  t2s[[paste0("active",       s2)]] <- t2$active
  t2s[[paste0("r_tf",         s2)]] <- t2$r_tf

  # Join
  jj <- merge(t1s, t2s, by = c("tf","gene_key","peak_id"), all.x = TRUE, all.y = FALSE, sort = FALSE)

  # Assert the expected columns exist (prevents silent NULL -> numeric(0))
  lk1 <- paste0("link_score", s1); lk2 <- paste0("link_score", s2)
  fp1 <- paste0("fp_bed_score", s1); fp2 <- paste0("fp_bed_score", s2)
  tf1 <- paste0("tf_expr", s1); tf2 <- paste0("tf_expr", s2)
  ge1 <- paste0("gene_expr", s1); ge2 <- paste0("gene_expr", s2)
  must_have <- c(lk1, lk2, fp1, fp2, tf1, tf2, ge1, ge2)
  missing <- setdiff(must_have, names(jj))
  if (length(missing)) {
    cli::cli_abort(c(
      "Joined table is missing expected columns.",
      i = "Missing: {paste(missing, collapse = ', ')}",
      i = "This usually means dynamic column naming failed upstream."
    ))
  }

  # Safe log2FC
  .safe_log2fc <- function(a, b, pc) {
    out <- suppressWarnings(log2((a + pc) / (b + pc)))
    out[!is.finite(out)] <- NA_real_
    out
  }

  # Deltas & log2FC
  jj[["delta_link_score"]]   <- as.numeric(jj[[lk1]]) - as.numeric(jj[[lk2]])
  jj[["delta_fp_bed_score"]] <- as.numeric(jj[[fp1]]) - as.numeric(jj[[fp2]])
  jj[["delta_tf_expr"]]      <- as.numeric(jj[[tf1]]) - as.numeric(jj[[tf2]])
  jj[["delta_gene_expr"]]    <- as.numeric(jj[[ge1]]) - as.numeric(jj[[ge2]])
  jj[["log2FC_tf_expr"]]     <- .safe_log2fc(as.numeric(jj[[tf1]]), as.numeric(jj[[tf2]]), pc)
  jj[["log2FC_gene_expr"]]   <- .safe_log2fc(as.numeric(jj[[ge1]]), as.numeric(jj[[ge2]]), pc)

  # Attach available stats from cond1 (optional)
  stat_cols <- intersect(c("r_gene","p_gene","p_adj_gene","r_tf","p_tf","p_adj_tf"),
                         intersect(names(t1), names(t2)))
  if (length(stat_cols)) {
    add <- t1[c("tf","gene_key","peak_id", stat_cols)]
    jj <- merge(jj, add, by = c("tf","gene_key","peak_id"), all.x = TRUE, sort = FALSE)
  }

  # Final column order
  left_cols  <- c("tf","gene_key","peak_id")
  per_cond   <- c(
    lk1, lk2,
    paste0("link_sign", s1), paste0("link_sign", s2),
    fp1, fp2, tf1, tf2, ge1, ge2,
    paste0("active", s1), paste0("active", s2),
    paste0("r_tf", s1), paste0("r_tf", s2)
  )
  delta_cols <- c("delta_link_score","delta_fp_bed_score","delta_tf_expr","delta_gene_expr")
  extra_cols <- c("log2FC_tf_expr","log2FC_gene_expr")
  keep <- unique(c(left_cols, per_cond, delta_cols, extra_cols, stat_cols))
  jj[keep]
}
# END EDIT


#' Bulk comparisons (optional parallel)
#'
#' Runs `compare_links_two_conditions()` across a spec table and optionally writes
#' per-contrast CSVs. Uses `future.apply` for parallelization.
#'
#' @param specs A data.frame/tibble with columns:
#'   `cond1_source`, `cond2_source`, `cond1_name`, `cond2_name`, `out_file`.
#' @inheritParams compare_links_two_conditions
#'
#' @return (invisible) list of result tibbles; also writes CSVs when `out_file` is set.
#' @export
#'
#' @importFrom readr write_csv
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession
#' @importFrom parallel detectCores
compare_links_bulk <- function(specs,
                               clean_names = FALSE,
                               parallel = TRUE,
                               plan = "multisession",
                               workers = NULL,
                               dedupe = TRUE,
                               pseudocount = 1,
                               verbose = TRUE) {
  if (!is.data.frame(specs) || nrow(specs) == 0L) cli::cli_abort("compare_links_bulk: empty specs.")

  do_one <- function(row) {
    c1 <- row$cond1_source
    c2 <- row$cond2_source
    n1 <- if (!is.null(row$cond1_name)) row$cond1_name else NULL
    n2 <- if (!is.null(row$cond2_name)) row$cond2_name else NULL
    out <- row$out_file

    res <- compare_links_two_conditions(
      c1, c2, n1, n2,
      clean_names        = clean_names,
      strict             = TRUE,
      dedupe             = dedupe,
      pseudocount        = pseudocount,
      verbose            = verbose
    )
    if (!is.null(out) && is.character(out) && nzchar(out)) {
      dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(res, out)
      .log(paste0("✓ wrote: ", basename(out)), verbose = verbose)
    }
    res
  }

  if (!isTRUE(parallel)) {
    invisible(lapply(split(specs, seq_len(nrow(specs))), do_one))
  } else {
    if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    strategy <- if (is.character(plan)) {
      ns <- asNamespace("future")
      if (exists(plan, envir = ns, mode = "function")) get(plan, envir = ns) else future::multisession
    } else if (is.function(plan)) plan else future::multisession

    .log(paste0("Launching future plan=",
                if (is.character(plan)) plan else "custom",
                ", workers=", workers), verbose = verbose)

    oplan <- future::plan()
    on.exit(future::plan(oplan), add = TRUE)
    future::plan(strategy, workers = workers)

    invisible(future.apply::future_lapply(split(specs, seq_len(nrow(specs))), do_one))
  }
}

#' # ==============================
#' # Build specs from index (cell-wise contrasts)
#' # ==============================
#'
#' #' Build per-cell “everything vs 10_FBS” contrasts from index CSV
#' #'
#' #' The index CSV must contain a `label` column (e.g., "HPAFII_10_FBS", "HPAFII_0_FBS", ...).
#' #' For each `cell` (prefix before first underscore), constructs stress vs control
#' #' pairs and emits a spec table consumable by `compare_links_bulk()`.
#' #'
#' #' @param index_csv Path to `lighting_per_condition_index.csv` (must contain `label`).
#' #' @param out_dir Directory containing per-condition links CSVs.
#' #' @param prefix File prefix (default `"lighting"`).
#' #' @param ctrl_tag Control tag appended after `<cell>_` (default `"10_FBS"`).
#' #' @param clean_names Logical; sanitize condition names in specs.
#' #' @param verbose Logical.
#' #'
#' #' @return A tibble with columns `cond1_label`, `cond2_label`, `cond1_source`,
#' #'   `cond2_source`, `cond1_name`, `cond2_name`, `out_file`.
#' #' @export
#' #'
#' #' @importFrom readr read_csv
#' #' @importFrom tibble tibble as_tibble
#' build_cellwise_contrasts_from_index <- function(index_csv,
#'                                                 out_dir,
#'                                                 prefix = "lighting",
#'                                                 ctrl_tag = "10_FBS",
#'                                                 clean_names = FALSE,
#'                                                 verbose = TRUE) {
#'   .log(paste0("Reading index: ", index_csv), verbose = verbose)
#'   idx <- readr::read_csv(index_csv, show_col_types = FALSE)
#'   if (!("label" %in% names(idx))) cli::cli_abort("Index missing required column `label`.")
#'
#'   # label -> cell (prefix before first underscore)
#'   label_chr <- as.character(idx[["label"]])
#'   cell_chr  <- sub("_.*$", "", label_chr)
#'
#'   ctrl_label <- paste0(cell_chr, "_", ctrl_tag)
#'   is_ctrl <- label_chr == ctrl_label
#'
#'   ctrl_df <- unique(data.frame(cell = cell_chr[is_ctrl],
#'                                control = label_chr[is_ctrl],
#'                                stringsAsFactors = FALSE))
#'
#'   stress_df <- unique(data.frame(cell = cell_chr[!is_ctrl],
#'                                  stress = label_chr[!is_ctrl],
#'                                  stringsAsFactors = FALSE))
#'
#'   pairs <- merge(ctrl_df, stress_df, by = "cell", all = FALSE, sort = TRUE)
#'   pairs <- pairs[order(pairs$cell, pairs$stress, method = "radix"), , drop = FALSE]
#'
#'   if (!nrow(pairs)) {
#'     cli::cli_warn("No contrasts built — ensure '<cell>_{ctrl_tag}' exists in labels.")
#'     return(tibble::tibble())
#'   }
#'
#'   mk_path <- function(lab) {
#'     file.path(out_dir, sprintf("%s_cond-%s_tf_gene_links.csv", prefix, .safe_label(lab)))
#'   }
#'   name_fun <- if (isTRUE(clean_names)) janitor::make_clean_names else identity
#'
#'   cond1_label <- pairs$stress
#'   cond2_label <- pairs$control
#'
#'   specs <- data.frame(
#'     cond1_label  = cond1_label,
#'     cond2_label  = cond2_label,
#'     cond1_source = vapply(cond1_label, mk_path, character(1)),
#'     cond2_source = vapply(cond2_label, mk_path, character(1)),
#'     cond1_name   = vapply(cond1_label, function(z) name_fun(z), character(1)),
#'     cond2_name   = vapply(cond2_label, function(z) name_fun(z), character(1)),
#'     out_file     = file.path(out_dir,
#'                              sprintf("%s_vs_%s_delta_links.csv",
#'                                      .safe_label(cond1_label), .safe_label(cond2_label))),
#'     stringsAsFactors = FALSE,
#'     check.names = FALSE
#'   )
#'
#'   # Warn about missing inputs but still return specs
#'   all_paths <- unique(c(specs$cond1_source, specs$cond2_source))
#'   for (p in all_paths) {
#'     if (!file.exists(p)) cli::cli_warn("File not found (will fail when run): {p}")
#'   }
#'
#'   tibble::as_tibble(specs)
#' }
# ===============================
# Build specs from index (cell-wise or global-control)
# ===============================

#' Build contrasts from lighting index: either per-cell or global control
#'
#' If the index contains an exact match to `ctrl_tag` (e.g., "Ctrl"), build
#' contrasts of **every other label vs that global control**.
#'
#' Otherwise, fall back to the original per-cell behavior where `label` has the
#' form `<cell>_<condition>` and the control is `<cell>_{ctrl_tag}` (e.g., "HPAFII_10_FBS").
#'
#' @param index_csv Path to `lighting_per_condition_index.csv` (must contain `label`).
#' @param out_dir Directory containing per-condition links CSVs.
#' @param prefix File prefix used by `light_by_condition()` (default `"lighting"`).
#' @param ctrl_tag Control tag. For global control, this is the **exact label** (e.g., `"Ctrl"`).
#'   For per-cell control, this is the suffix after `<cell>_` (e.g., `"10_FBS"`).
#' @param clean_names Logical; sanitize condition names in specs with `janitor::make_clean_names()`.
#' @param verbose Logical.
#'
#' @return Tibble with columns `cond1_label`, `cond2_label`, `cond1_source`,
#'   `cond2_source`, `cond1_name`, `cond2_name`, `out_file`.
#' @export
#'
#' @importFrom readr read_csv
#' @importFrom tibble tibble as_tibble
build_cellwise_contrasts_from_index <- function(index_csv,
                                                out_dir,
                                                prefix = "lighting",
                                                ctrl_tag = "10_FBS",
                                                clean_names = FALSE,
                                                verbose = TRUE) {
  .log <- function(msg, verbose = TRUE) if (isTRUE(verbose)) message("[utils_grn_diff] ", msg)
  .safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

  .log(paste0("Reading index:\n  ", index_csv), verbose = verbose)
  idx <- readr::read_csv(index_csv, show_col_types = FALSE)
  if (!("label" %in% names(idx))) cli::cli_abort("Index missing required column `label`.")

  label_chr <- as.character(idx[["label"]])

  # --- Branch 1: Global control present (exact match to ctrl_tag) ---
  if (ctrl_tag %in% label_chr) {
    ctrl_label <- ctrl_tag
    stress_labels <- setdiff(label_chr, ctrl_label)
    if (!length(stress_labels)) {
      cli::cli_warn("No contrasts built — only the global control `{ctrl_tag}` exists.")
      return(tibble::tibble())
    }

    mk_path <- function(lab) file.path(out_dir, sprintf("%s_cond-%s_tf_gene_links.csv",
                                                        prefix, .safe_label(lab)))
    name_fun <- if (isTRUE(clean_names)) janitor::make_clean_names else identity

    specs <- data.frame(
      cond1_label  = stress_labels,
      cond2_label  = rep(ctrl_label, length(stress_labels)),
      cond1_source = vapply(stress_labels, mk_path, character(1)),
      cond2_source = rep(mk_path(ctrl_label), length(stress_labels)),
      cond1_name   = vapply(stress_labels, function(z) name_fun(z), character(1)),
      cond2_name   = rep(name_fun(ctrl_label), length(stress_labels)),
      out_file     = file.path(out_dir,
                               sprintf("%s_vs_%s_delta_links.csv",
                                       .safe_label(stress_labels), .safe_label(ctrl_label))),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    # warn about missing inputs but still return specs
    all_paths <- unique(c(specs$cond1_source, specs$cond2_source))
    for (p in all_paths) if (!file.exists(p)) cli::cli_warn("File not found (will fail when run): {p}")

    return(tibble::as_tibble(specs))
  }

  # --- Branch 2: Per-cell control (<cell>_{ctrl_tag}) ---
  cell_chr  <- sub("_.*$", "", label_chr)              # prefix before first underscore
  ctrl_label <- paste0(cell_chr, "_", ctrl_tag)
  is_ctrl <- label_chr == ctrl_label

  ctrl_df <- unique(data.frame(cell = cell_chr[is_ctrl],
                               control = label_chr[is_ctrl],
                               stringsAsFactors = FALSE))
  stress_df <- unique(data.frame(cell = cell_chr[!is_ctrl],
                                 stress = label_chr[!is_ctrl],
                                 stringsAsFactors = FALSE))
  pairs <- merge(ctrl_df, stress_df, by = "cell", all = FALSE, sort = TRUE)
  pairs <- pairs[order(pairs$cell, pairs$stress, method = "radix"), , drop = FALSE]

  if (!nrow(pairs)) {
    cli::cli_warn("No contrasts built — ensure either a global control '{ctrl_tag}' exists in labels, or per-cell '<cell>_{ctrl_tag}' labels are present.")
    return(tibble::tibble())
  }

  mk_path <- function(lab) file.path(out_dir, sprintf("%s_cond-%s_tf_gene_links.csv",
                                                      prefix, .safe_label(lab)))
  name_fun <- if (isTRUE(clean_names)) janitor::make_clean_names else identity

  cond1_label <- pairs$stress
  cond2_label <- pairs$control

  specs <- data.frame(
    cond1_label  = cond1_label,
    cond2_label  = cond2_label,
    cond1_source = vapply(cond1_label, mk_path, character(1)),
    cond2_source = vapply(cond2_label, mk_path, character(1)),
    cond1_name   = vapply(cond1_label, function(z) name_fun(z), character(1)),
    cond2_name   = vapply(cond2_label, function(z) name_fun(z), character(1)),
    out_file     = file.path(out_dir,
                             sprintf("%s_vs_%s_delta_links.csv",
                                     .safe_label(cond1_label), .safe_label(cond2_label))),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  all_paths <- unique(c(specs$cond1_source, specs$cond2_source))
  for (p in all_paths) if (!file.exists(p)) cli::cli_warn("File not found (will fail when run): {p}")

  tibble::as_tibble(specs)
}

# ==============================
# Driver
# ==============================
#' Run link-delta comparisons from a prebuilt specs table
#'
#' @param specs Tibble/data.frame with columns:
#'   `cond1_source`, `cond2_source`, `cond1_name`, `cond2_name`, `out_file`.
#' @param clean_names Logical; pass-through to `compare_links_bulk()`.
#' @param parallel Logical; pass-through to `compare_links_bulk()`.
#' @param plan Future plan name or function; pass-through to `compare_links_bulk()`.
#' @param workers Integer or NULL; pass-through to `compare_links_bulk()`.
#' @param dedupe Logical; pass-through to `compare_links_bulk()`.
#' @param pseudocount Numeric; pass-through to `compare_links_bulk()`.
#' @param verbose Logical; pass-through to `compare_links_bulk()`.
#'
#' @return (invisible) list of result tibbles.
#' @export
run_links_deltas_driver <- function(specs,
                                    clean_names = FALSE,
                                    parallel    = TRUE,
                                    plan        = "multisession",
                                    workers     = NULL,
                                    dedupe      = TRUE,
                                    pseudocount = 1,
                                    verbose     = TRUE) {
  if (!is.data.frame(specs) || nrow(specs) == 0L) {
    cli::cli_abort("`specs` must be a non-empty data.frame/tibble.")
  }
  need <- c("cond1_source","cond2_source","cond1_name","cond2_name","out_file")
  miss <- setdiff(need, names(specs))
  if (length(miss)) {
    cli::cli_abort(c("`specs` missing required columns.", i = paste(miss, collapse = ", ")))
  }

  # Optional: warn if any input files are absent (do not stop)
  paths <- unique(c(specs$cond1_source, specs$cond2_source))
  for (p in paths) if (!file.exists(p)) cli::cli_warn("File not found (may fail later): {p}")

  compare_links_bulk(specs,
                     clean_names = clean_names,
                     parallel    = parallel,
                     plan        = plan,
                     workers     = workers,
                     dedupe      = dedupe,
                     pseudocount = pseudocount,
                     verbose     = verbose)
}

#' #' Driver: run all cell-wise contrasts and write CSVs
#' #'
#' #' Convenience wrapper that builds specs from an index and runs
#' #' `compare_links_bulk()` with your chosen parallel settings.
#' #'
#' #' @inheritParams compare_links_two_conditions
#' #' @param out_dir Base directory with inputs and where outputs are written.
#' #' @param prefix File prefix (default `"lighting"`).
#' #' @param index_csv Path to index CSV (default `file.path(out_dir, "lighting_per_condition_index.csv")`).
#' #' @param ctrl_tag Control tag (default `"10_FBS"`).
#' #'
#' #' @return (invisible) list of result tibbles.
#' #' @export
#' run_links_deltas_driver <- function(out_dir = "inst/extdata/lighting",
#'                                     prefix = "lighting",
#'                                     index_csv = file.path(out_dir, "lighting_per_condition_index.csv"),
#'                                     ctrl_tag = "10_FBS",
#'                                     clean_names = FALSE,
#'                                     parallel = TRUE,
#'                                     plan = "multisession",
#'                                     workers = NULL,
#'                                     dedupe = TRUE,
#'                                     pseudocount = 1,
#'                                     verbose = TRUE) {
#'
#'   specs <- build_cellwise_contrasts_from_index(index_csv, out_dir, prefix, ctrl_tag,
#'                                                clean_names = clean_names, verbose = verbose)
#'   if (!nrow(specs)) return(invisible(list()))
#'
#'   compare_links_bulk(specs,
#'                      clean_names        = clean_names,
#'                      parallel           = parallel,
#'                      plan               = plan,
#'                      workers            = workers,
#'                      dedupe             = dedupe,
#'                      pseudocount        = pseudocount,
#'                      verbose            = verbose)
#' }

# ---------------------------------------------------------------------------
# Examples (not run)
# ---------------------------------------------------------------------------
# res <- compare_links_two_conditions(
#   "inst/extdata/lighting/lighting_cond-HPAFII_0_FBS_tf_gene_links.csv",
#   "inst/extdata/lighting/lighting_cond-HPAFII_10_FBS_tf_gene_links.csv",
#   clean_names   = FALSE,
#   dedupe        = TRUE,
#   dedupe_method = "max_link_score",
#   sign_policy   = "error"
# )
#
# specs <- build_cellwise_contrasts_from_index(
#   index_csv = "inst/extdata/lighting/lighting_per_condition_index.csv",
#   out_dir   = "inst/extdata/lighting",
#   prefix    = "lighting",
#   ctrl_tag  = "10_FBS",
#   clean_names = FALSE
# )
# run_links_deltas_driver(
#   out_dir = "inst/extdata/lighting",
#   prefix  = "lighting",
#   index_csv = "inst/extdata/lighting/lighting_per_condition_index.csv",
#   ctrl_tag = "10_FBS",
#   clean_names = FALSE,
#   parallel = TRUE
# )
