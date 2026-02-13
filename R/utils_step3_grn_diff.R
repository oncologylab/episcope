# FILE: R/utils_step3_grn_diff.R
# ---------------------------------------------------------------------------
# General TF-gene links delta computation
# Author: Yaoxiang Li
# Updated: 2025-11-10
#
# Summary
#   Load TF-gene link tables, harmonize columns, and compute per-link deltas
#   for link_score, fp_bed_score, tf_expr, gene_expr, active (cond1 - cond2).
#   Includes bulk/driver helpers for per-cell contrasts.
#
# Conventions
#   - roxygen2 docs; explicit namespacing
#   - .log_abort()/cli_warn()/cli_inform()
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
  if (isTRUE(verbose)) .log_inform(c(v = paste0("[utils_grn_diff] ", msg)))
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
    .log_abort(c(
      "Duplicate keys detected in cond1.",
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }
  if (nrow(dup2)) {
    ex <- utils::head(unique(dup2), 5)
    .log_abort(c(
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
    .log_abort(c(
      "Some keys exist only in cond1.",
      i = "Examples:\n{paste(utils::capture.output(print(ex)), collapse = '\n')}"
    ))
  }
  if (nrow(only_in_2)) {
    ex <- utils::head(only_in_2, 5)
    .log_abort(c(
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
    .log_abort(c(
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

#' Load TF-gene links table from CSV or tibble
#'
#' Harmonizes common columns, ensures required IDs, and standardizes `link_sign`.
#' If `link_sign` is absent, it is derived from `r_gene` when available, else
#' from the sign of `link_score` (and `link_score` is abs()-ed accordingly).
#'
#' @param x A `data.frame`/`tibble` **or** a character file path to a CSV.
#' @param verbose Logical; emit progress.
#' @param keep_all_cols Logical; if TRUE, preserve extra columns beyond the
#'   canonical link fields (default TRUE).
#'
#' @return A tibble with standardized columns and types.
#' @export
#'
#' @importFrom readr read_csv
#' @importFrom tibble as_tibble
#' @importFrom janitor clean_names make_clean_names
#' @importFrom cli cli_warn
load_links_table <- function(x, verbose = TRUE, keep_all_cols = TRUE) {
  .log("Loading links table.", verbose = verbose)

  tbl <- if (is.character(x) && length(x) == 1L) {
    .log(paste0("from CSV: ", x), verbose = verbose)
    readr::read_csv(x, show_col_types = FALSE)
  } else if (inherits(x, "data.frame")) {
    .log("from in-memory tibble/data.frame", verbose = verbose)
    x
  } else {
    .log_abort("`x` must be a file path or a data.frame/tibble.")
  }

  tbl <- tibble::as_tibble(janitor::clean_names(tbl))

  # Required IDs
  req_ids <- c("tf", "gene_key", "peak_id")
  miss <- setdiff(req_ids, names(tbl))
  if (length(miss)) .log_abort("Missing required columns: {paste(miss, collapse = ', ')}")

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
      .log_warn("No link_score/edge_weight/link_weight/weight found; created NA link_score.")
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
    if (invalid > 0) .log_warn("{invalid} link_sign value(s) unrecognized; set to NA.")
    tbl[["link_sign"]] <- mapped
  }

  # Keep canonical columns (plus common stats if present); optionally preserve extras
  keep <- c("tf","gene_key","peak_id","link_score","fp_bed_score","tf_expr","gene_expr","active",
            "link_sign","r_gene","p_gene","p_adj_gene","r_tf","p_tf","p_adj_tf")
  keep <- intersect(keep, names(tbl))
  if (isTRUE(keep_all_cols)) {
    rest <- setdiff(names(tbl), keep)
    tbl <- tbl[c(keep, rest)]
  } else {
    tbl <- tbl[keep]
  }
  tbl
}

# ==============================
# Public API: Compare (Two / Bulk)
# ==============================

#' Compare two conditions (cond1 - cond2) for TF-gene links
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
#' @param restrict_to_active_both Logical; if TRUE, set delta/log2FC columns to NA
#'   unless a link is active in both conditions (default TRUE).
#' @param keep_all_cols Logical; if TRUE, retain extra columns from per-condition
#'   inputs and carry them through (default TRUE).
#' @param edge_change_min Numeric; minimum |delta_link_score| used to classify
#'   edges as gain/loss/shared (default 0).
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
                                         restrict_to_active_both = TRUE,
                                         edge_change_min = 0,
                                         keep_all_cols = TRUE,
                                         verbose = TRUE) {
  dedupe_method <- match.arg(dedupe_method)
  sign_policy   <- match.arg(sign_policy)

  pc <- as.numeric(pseudocount)
  if (!is.finite(pc) || pc < 0) .log_abort("`pseudocount` must be a non-negative finite number.")
  edge_min <- as.numeric(edge_change_min)
  if (!is.finite(edge_min) || edge_min < 0) {
    .log_warn("edge_change_min must be a non-negative finite number; defaulting to 0.")
    edge_min <- 0
  }

  t1 <- load_links_table(cond1, verbose = verbose, keep_all_cols = keep_all_cols)
  t2 <- load_links_table(cond2, verbose = verbose, keep_all_cols = keep_all_cols)

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
  t1s[[paste0("r_tf",         s1)]] <- if ("r_tf" %in% names(t1)) t1$r_tf else NA_real_

  t2s <- data.frame(tf = t2$tf, gene_key = t2$gene_key, peak_id = t2$peak_id, check.names = FALSE)
  t2s[[paste0("link_score",   s2)]] <- t2$link_score
  t2s[[paste0("link_sign",    s2)]] <- t2$link_sign
  t2s[[paste0("fp_bed_score", s2)]] <- t2$fp_bed_score
  t2s[[paste0("tf_expr",      s2)]] <- t2$tf_expr
  t2s[[paste0("gene_expr",    s2)]] <- t2$gene_expr
  t2s[[paste0("active",       s2)]] <- t2$active
  t2s[[paste0("r_tf",         s2)]] <- if ("r_tf" %in% names(t2)) t2$r_tf else NA_real_

  static_cols <- intersect(
    c("n_used_tf", "motif",
      "r_gene", "p_gene", "p_adj_gene", "p_tf", "p_adj_tf",
      "r_rna_gene", "p_rna_gene", "p_rna_adj_gene",
      "fp_mean_raw", "fp_var_raw", "fp_rsd",
      "gene_mean_raw", "gene_var_raw", "gene_rsd"),
    intersect(names(t1), names(t2))
  )
  base_cols <- c(key_cols,
                 "link_score", "link_sign", "fp_bed_score", "tf_expr", "gene_expr", "active", "r_tf",
                 static_cols)
  extra_cols <- setdiff(intersect(names(t1), names(t2)), base_cols)
  extra_cols <- setdiff(extra_cols, c("active_link", "fp_score", "footprint_score", "fp"))

  if (length(extra_cols)) {
    for (cc in extra_cols) {
      t1s[[paste0(cc, s1)]] <- t1[[cc]]
      t2s[[paste0(cc, s2)]] <- t2[[cc]]
    }
  }

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
    .log_abort(c(
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
  .safe_fc <- function(a, b, pc) {
    out <- suppressWarnings((a + pc) / (b + pc))
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
  jj[["fc_link_score"]]      <- .safe_fc(as.numeric(jj[[lk1]]), as.numeric(jj[[lk2]]), pc)
  jj[["fc_fp_bed_score"]]    <- .safe_fc(as.numeric(jj[[fp1]]), as.numeric(jj[[fp2]]), pc)
  jj[["fc_tf_expr"]]         <- .safe_fc(as.numeric(jj[[tf1]]), as.numeric(jj[[tf2]]), pc)
  jj[["fc_gene_expr"]]       <- .safe_fc(as.numeric(jj[[ge1]]), as.numeric(jj[[ge2]]), pc)

  rk1 <- suppressWarnings(rank(-as.numeric(jj[[lk1]]), ties.method = "average", na.last = "keep"))
  rk2 <- suppressWarnings(rank(-as.numeric(jj[[lk2]]), ties.method = "average", na.last = "keep"))
  jj[[paste0("rank_link_score", s1)]] <- rk1
  jj[[paste0("rank_link_score", s2)]] <- rk2
  jj[["rank_shift_link_score"]] <- rk1 - rk2

  dlink <- as.numeric(jj[["delta_link_score"]])
  edge_change <- ifelse(is.na(dlink), NA_character_,
                        ifelse(dlink >  edge_min, "gain",
                               ifelse(dlink < -edge_min, "loss", "shared")))
  jj[["edge_change"]] <- edge_change

  # Active gating (both conditions)
  a1_raw <- jj[[paste0("active", s1)]]
  a2_raw <- jj[[paste0("active", s2)]]
  a1 <- as.logical(a1_raw)
  a2 <- as.logical(a2_raw)
  a1[is.na(a1)] <- FALSE
  a2[is.na(a2)] <- FALSE
  jj[["active_both"]] <- a1 & a2
  jj[["active_any"]]  <- a1 | a2

  has_active <- any(!is.na(a1_raw)) && any(!is.na(a2_raw))
  if (isTRUE(restrict_to_active_both) && isTRUE(has_active)) {
    idx_off <- !jj[["active_both"]]
    if (any(idx_off)) {
      jj[["delta_link_score"]][idx_off]   <- NA_real_
      jj[["delta_fp_bed_score"]][idx_off] <- NA_real_
      jj[["delta_tf_expr"]][idx_off]      <- NA_real_
      jj[["delta_gene_expr"]][idx_off]    <- NA_real_
      jj[["log2FC_tf_expr"]][idx_off]     <- NA_real_
      jj[["log2FC_gene_expr"]][idx_off]   <- NA_real_
      jj[["fc_link_score"]][idx_off]      <- NA_real_
      jj[["fc_fp_bed_score"]][idx_off]    <- NA_real_
      jj[["fc_tf_expr"]][idx_off]         <- NA_real_
      jj[["fc_gene_expr"]][idx_off]       <- NA_real_
      jj[["rank_shift_link_score"]][idx_off] <- NA_real_
      jj[["edge_change"]][idx_off]        <- NA_character_
    }
  }

  # Attach available static stats from cond1 (optional)
  if (length(static_cols)) {
    add <- t1[c("tf","gene_key","peak_id", static_cols)]
    jj <- merge(jj, add, by = c("tf","gene_key","peak_id"), all.x = TRUE, sort = FALSE)
  }

  # Final column order
  left_cols  <- c("tf","gene_key","peak_id")
  per_cond   <- c(
    lk1, lk2,
    paste0("link_sign", s1), paste0("link_sign", s2),
    fp1, fp2, tf1, tf2, ge1, ge2,
    paste0("active", s1), paste0("active", s2),
    paste0("r_tf", s1), paste0("r_tf", s2),
    paste0("rank_link_score", s1), paste0("rank_link_score", s2)
  )
  per_cond_extra <- unlist(lapply(extra_cols, function(z) c(paste0(z, s1), paste0(z, s2))))
  delta_cols <- c("delta_link_score","delta_fp_bed_score","delta_tf_expr","delta_gene_expr")
  log2_cols <- c("log2FC_tf_expr","log2FC_gene_expr")
  fc_cols <- c("fc_link_score","fc_fp_bed_score","fc_tf_expr","fc_gene_expr")
  keep <- unique(c(left_cols, per_cond, per_cond_extra,
                   delta_cols, log2_cols, fc_cols,
                   "rank_shift_link_score", "edge_change",
                   "active_both", "active_any",
                   static_cols))
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
                               restrict_to_active_both = TRUE,
                               edge_change_min = 0,
                               keep_all_cols = TRUE,
                               verbose = TRUE) {
  if (!is.data.frame(specs) || nrow(specs) == 0L) .log_abort("compare_links_bulk: empty specs.")
  n_specs <- nrow(specs)
  if (isTRUE(verbose)) {
    .log_inform("Computing link deltas for {n_specs} comparison(s).", n_specs = n_specs)
  }
  worker_verbose <- if (isTRUE(parallel)) FALSE else verbose

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
      restrict_to_active_both = restrict_to_active_both,
      edge_change_min    = edge_change_min,
      keep_all_cols      = keep_all_cols,
      verbose            = worker_verbose
    )
    if (!is.null(out) && is.character(out) && nzchar(out)) {
      dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(res, out)
      if (isTRUE(verbose) && !isTRUE(parallel)) {
        .log(paste0("?wrote: ", basename(out)), verbose = verbose)
      }
    }
    res
  }

  if (!isTRUE(parallel)) {
    invisible(lapply(split(specs, seq_len(nrow(specs))), do_one))
  } else {
    if (is.null(workers)) {
      workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    }
    if (workers > n_specs) {
      workers <- n_specs
    }
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

    if (isTRUE(verbose) && requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(along = seq_len(n_specs))
        do_one_wrap <- function(row) {
          on.exit(p(), add = TRUE)
          do_one(row)
        }
        invisible(future.apply::future_lapply(split(specs, seq_len(n_specs)), do_one_wrap))
      })
    } else {
      invisible(future.apply::future_lapply(split(specs, seq_len(n_specs)), do_one))
    }
  }
}

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
#' @param out_dir Directory where delta outputs should be written.
#' @param input_dir Directory containing per-condition links CSVs. Defaults to
#'   `out_dir` for backward compatibility.
#' @param prefix File prefix used by `light_by_condition()` (default `"lighting"`).
#' @param input_prefix File prefix for input per-condition links (defaults to `prefix`).
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
                                                input_dir = out_dir,
                                                prefix = "lighting",
                                                input_prefix = prefix,
                                                ctrl_tag = "10_FBS",
                                                clean_names = FALSE,
                                                labels_only = FALSE,
                                                verbose = TRUE) {
  .log <- function(msg, verbose = TRUE) {
    if (isTRUE(verbose)) .log_inform(c(v = paste0("[utils_grn_diff] ", msg)))
  }
  .safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
  .mk_cond_link_path <- function(lab) {
    s_lab <- .safe_label(lab)
    new_path <- file.path(input_dir, sprintf("%s_%s_tf_gene_links.csv", input_prefix, s_lab))
    new_path_sub <- file.path(input_dir, "per_condition_link_matrices", sprintf("%s_%s_tf_gene_links.csv", input_prefix, s_lab))
    old_path <- file.path(input_dir, sprintf("%s_cond-%s_tf_gene_links.csv", input_prefix, s_lab))
    old_path_sub <- file.path(input_dir, "per_condition_link_matrices", sprintf("%s_cond-%s_tf_gene_links.csv", input_prefix, s_lab))
    if (file.exists(new_path)) return(new_path)
    if (file.exists(new_path_sub)) return(new_path_sub)
    if (file.exists(old_path_sub)) return(old_path_sub)
    if (file.exists(old_path)) return(old_path)
    new_path_sub
  }

  .log(paste0("Reading index:\n  ", index_csv), verbose = verbose)
  idx <- readr::read_csv(index_csv, show_col_types = FALSE)
  if (!("label" %in% names(idx))) .log_abort("Index missing required column `label`.")

  label_chr <- as.character(idx[["label"]])

  # --- Branch 1: Global control present (exact match to ctrl_tag) ---
  if (ctrl_tag %in% label_chr) {
    ctrl_label <- ctrl_tag
    stress_labels <- setdiff(label_chr, ctrl_label)
    if (!length(stress_labels)) {
      .log_warn("No contrasts built ?only the global control `{ctrl_tag}` exists.")
      return(tibble::tibble())
    }

    mk_path <- function(lab) .mk_cond_link_path(lab)
    name_fun <- if (isTRUE(clean_names)) janitor::make_clean_names else identity

    specs <- data.frame(
      cond1_label  = stress_labels,
      cond2_label  = rep(ctrl_label, length(stress_labels)),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    if (!isTRUE(labels_only)) {
      specs$cond1_source <- vapply(stress_labels, mk_path, character(1))
      specs$cond2_source <- rep(mk_path(ctrl_label), length(stress_labels))
      specs$cond1_name   <- vapply(stress_labels, function(z) name_fun(z), character(1))
      specs$cond2_name   <- rep(name_fun(ctrl_label), length(stress_labels))
      specs$out_file     <- file.path(out_dir,
                                      sprintf("%s_vs_%s_delta_links.csv",
                                              .safe_label(stress_labels), .safe_label(ctrl_label)))
    }

    # warn about missing inputs but still return specs
    if (!isTRUE(labels_only)) {
      all_paths <- unique(c(specs$cond1_source, specs$cond2_source))
      for (p in all_paths) if (!file.exists(p)) .log_warn("File not found (will fail when run): {p}")
    }

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
    .log_warn("No contrasts built ?ensure either a global control '{ctrl_tag}' exists in labels, or per-cell '<cell>_{ctrl_tag}' labels are present.")
    return(tibble::tibble())
  }

  mk_path <- function(lab) .mk_cond_link_path(lab)
  name_fun <- if (isTRUE(clean_names)) janitor::make_clean_names else identity

  cond1_label <- pairs$stress
  cond2_label <- pairs$control

  specs <- data.frame(
    cond1_label  = cond1_label,
    cond2_label  = cond2_label,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  if (!isTRUE(labels_only)) {
    specs$cond1_source <- vapply(cond1_label, mk_path, character(1))
    specs$cond2_source <- vapply(cond2_label, mk_path, character(1))
    specs$cond1_name   <- vapply(cond1_label, function(z) name_fun(z), character(1))
    specs$cond2_name   <- vapply(cond2_label, function(z) name_fun(z), character(1))
    specs$out_file     <- file.path(out_dir,
                                    sprintf("%s_vs_%s_delta_links.csv",
                                            .safe_label(cond1_label), .safe_label(cond2_label)))
  }

  if (!isTRUE(labels_only)) {
    all_paths <- unique(c(specs$cond1_source, specs$cond2_source))
    for (p in all_paths) if (!file.exists(p)) .log_warn("File not found (will fail when run): {p}")
  }

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
#' @param restrict_to_active_both Logical; if TRUE, set delta/log2FC columns to NA
#'   unless a link is active in both conditions (default TRUE).
#' @param keep_all_cols Logical; if TRUE, retain extra per-condition columns (default TRUE).
#' @param edge_change_min Numeric; minimum |delta_link_score| used to classify
#'   edges as gain/loss/shared (default 0).
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
                                    restrict_to_active_both = TRUE,
                                    edge_change_min = 0,
                                    keep_all_cols = TRUE,
                                    input_dir = NULL,
                                    input_prefix = "lighting",
                                    out_dir = NULL,
                                    verbose     = TRUE) {
  if (!is.data.frame(specs) || nrow(specs) == 0L) {
    .log_abort("`specs` must be a non-empty data.frame/tibble.")
  }
  if (!all(c("cond1_label", "cond2_label") %in% names(specs))) {
    .log_abort("`specs` must include cond1_label and cond2_label.")
  }

  .safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
  .mk_cond_link_path <- function(lab) {
    s_lab <- .safe_label(lab)
    new_path <- file.path(input_dir, sprintf("%s_%s_tf_gene_links.csv", input_prefix, s_lab))
    new_path_sub <- file.path(input_dir, "per_condition_link_matrices", sprintf("%s_%s_tf_gene_links.csv", input_prefix, s_lab))
    old_path <- file.path(input_dir, sprintf("%s_cond-%s_tf_gene_links.csv", input_prefix, s_lab))
    old_path_sub <- file.path(input_dir, "per_condition_link_matrices", sprintf("%s_cond-%s_tf_gene_links.csv", input_prefix, s_lab))
    if (file.exists(new_path)) return(new_path)
    if (file.exists(new_path_sub)) return(new_path_sub)
    if (file.exists(old_path_sub)) return(old_path_sub)
    if (file.exists(old_path)) return(old_path)
    new_path_sub
  }
  name_fun <- if (isTRUE(clean_names)) janitor::make_clean_names else identity

  if (!all(c("cond1_source", "cond2_source", "cond1_name", "cond2_name", "out_file") %in% names(specs))) {
    if (is.null(input_dir) || is.null(out_dir)) {
      .log_abort("`specs` is missing path columns; set `input_dir` and `out_dir` to derive them.")
    }
    specs$cond1_source <- vapply(specs$cond1_label, .mk_cond_link_path, character(1))
    specs$cond2_source <- vapply(specs$cond2_label, .mk_cond_link_path, character(1))
    specs$cond1_name <- vapply(specs$cond1_label, function(z) name_fun(z), character(1))
    specs$cond2_name <- vapply(specs$cond2_label, function(z) name_fun(z), character(1))
    specs$out_file <- file.path(
      out_dir,
      sprintf(
        "%s_vs_%s_delta_links.csv",
        .safe_label(specs$cond1_label),
        .safe_label(specs$cond2_label)
      )
    )
  }

  # Optional: warn if any input files are absent (do not stop)
  paths <- unique(c(specs$cond1_source, specs$cond2_source))
  for (p in paths) if (!file.exists(p)) .log_warn("File not found (may fail later): {p}")

  compare_links_bulk(specs,
                     clean_names = clean_names,
                     parallel    = parallel,
                     plan        = plan,
                     workers     = workers,
                     dedupe      = dedupe,
                     pseudocount = pseudocount,
                     restrict_to_active_both = restrict_to_active_both,
                     edge_change_min = edge_change_min,
                     keep_all_cols = keep_all_cols,
                     verbose     = verbose)
}

#' Find differential links (Module 3)
#'
#' Runs link-delta comparisons and filtering using a config + comparison list.
#'
#' @param config YAML config path or preloaded config (global env).
#' @param compar Optional CSV path with columns `cond1_label`, `cond2_label`.
#' @param output_dir Output directory for differential links (will create subfolders).
#' @param qc Logical; reserved for future QC outputs.
#' @param input_dir Optional override for step2 link directory (default: base_dir/connect_tf_target_genes).
#' @param input_prefix File prefix for step2 links (default "step2").
#' @param ctrl_tag Control tag for per-cell contrasts (default "10_FBS").
#' @param clean_names Logical; pass-through to comparisons builder.
#' @param parallel Logical; pass-through to delta computation.
#' @param workers Integer; pass-through to delta computation.
#' @param summary_tag Optional tag to insert into TF hub summary filenames.
#' @param write_tf_hubs_fp If TRUE, write TF hub summaries using delta FP scores.
#' @param overwrite_delta Logical; if TRUE, recompute `diff_links` even if files exist.
#' @param overwrite_filtered Logical; if TRUE, recompute `diff_links_filtered` even if files exist.
#' @param overwrite_tf_hubs Logical; if TRUE, always regenerate `master_tf_summary` outputs.
#' @param connectivity_min_degree Integer cutoff for TF inclusion in connectivity
#'   heatmaps. A TF is shown only if it connects to at least this many other TFs
#'   (undirected degree).
#' @param summary_plot_format Output format for QC differential-link summary plots:
#'   "pdf", "html", or "both".
#' @param summary_html_selfcontained Logical; passed to
#'   \code{htmlwidgets::saveWidget()} when writing HTML summaries.
#' @export
find_differential_links <- function(config,
                                    compar = NULL,
                                    output_dir = "diff_grn_and_regulatory_topics",
                                    qc = TRUE,
                                    input_dir = NULL,
                                    input_prefix = "step2",
                                    ctrl_tag = "10_FBS",
                                    clean_names = FALSE,
                                    parallel = TRUE,
                                    workers = NULL,
                                    summary_tag = NULL,
                                    write_tf_hubs_fp = TRUE,
                                    overwrite_delta = FALSE,
                                    overwrite_filtered = FALSE,
                                    overwrite_tf_hubs = TRUE,
                                    connectivity_min_degree = 1L,
                                    summary_plot_format = c("pdf", "html", "both"),
                                    summary_html_selfcontained = TRUE) {
  if (is.character(config) && length(config) == 1L && file.exists(config)) {
    load_config(config)
  }
  if (!exists("base_dir")) {
    .log_abort("`base_dir` not found. Provide a valid config or load_config() first.")
  }

  step2_out_dir <- if (!is.null(input_dir)) input_dir else file.path(base_dir, "connect_tf_target_genes")
  out_root <- if (grepl("^/", output_dir)) output_dir else file.path(base_dir, output_dir)
  delta_dir <- file.path(out_root, "diff_links")
  filtered_dir <- file.path(out_root, "diff_links_filtered")

  legacy_delta_dir <- file.path(out_root, "delta_links")
  legacy_filtered_dir <- file.path(out_root, "filtered_delta_links")
  if (!dir.exists(delta_dir) && dir.exists(legacy_delta_dir)) {
    moved <- file.rename(legacy_delta_dir, delta_dir)
    if (isTRUE(moved)) {
      .log_inform("Migrated legacy folder {legacy_delta_dir} -> {delta_dir}.")
    } else {
      .log_warn("Could not migrate legacy folder {legacy_delta_dir}; continuing with {delta_dir}.")
    }
  }
  if (!dir.exists(filtered_dir) && dir.exists(legacy_filtered_dir)) {
    moved <- file.rename(legacy_filtered_dir, filtered_dir)
    if (isTRUE(moved)) {
      .log_inform("Migrated legacy folder {legacy_filtered_dir} -> {filtered_dir}.")
    } else {
      .log_warn("Could not migrate legacy folder {legacy_filtered_dir}; continuing with {filtered_dir}.")
    }
  }

  dir.create(delta_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(compar)) {
    if (is.data.frame(compar)) {
      comp_tbl <- compar
    } else if (is.character(compar) && length(compar) == 1L) {
      comp_tbl <- readr::read_csv(compar, show_col_types = FALSE)
    } else {
      .log_abort("compar must be a data.frame/tibble or a path to a CSV.")
    }
    if (!all(c("cond1_label", "cond2_label") %in% names(comp_tbl))) {
      .log_abort("compar must include cond1_label and cond2_label.")
    }
    specs <- comp_tbl[, c("cond1_label", "cond2_label")]
  } else {
    idx_csv <- file.path(step2_out_dir, "per_condition_link_matrices", "step2_per_condition_index.csv")
    if (!file.exists(idx_csv)) {
      idx_csv <- file.path(step2_out_dir, "step2_per_condition_index.csv")
    }
    specs <- build_cellwise_contrasts_from_index(
      index_csv = idx_csv,
      out_dir = delta_dir,
      input_dir = step2_out_dir,
      prefix = "step3",
      input_prefix = input_prefix,
      ctrl_tag = ctrl_tag,
      clean_names = clean_names,
      labels_only = TRUE
    )
  }

  delta_link <- if (exists("delta_link")) delta_link else 1
  de_gene_log2_abs_min <- if (exists("de_gene_log2_abs_min")) de_gene_log2_abs_min else 0.5
  de_tf_log2_abs_min <- if (exists("de_tf_log2_abs_min")) de_tf_log2_abs_min else 0.5
  fp_delta_min <- if (exists("fp_delta_min")) fp_delta_min else 0.5

  existing_delta <- list.files(delta_dir, "_delta_links\\.csv$", full.names = TRUE)
  if (!isTRUE(overwrite_delta) && length(existing_delta) > 0L) {
    .log_inform("Using existing delta links from {delta_dir} (set overwrite_delta=TRUE to rebuild).")
  } else {
    run_links_deltas_driver(
      specs = specs,
      clean_names = clean_names,
      parallel = parallel,
      workers = workers,
      restrict_to_active_both = FALSE,
      edge_change_min = delta_link,
      keep_all_cols = TRUE,
      input_dir = step2_out_dir,
      input_prefix = input_prefix,
      out_dir = delta_dir
    )
  }

  delta_csvs <- list.files(delta_dir, "_delta_links.csv$", full.names = TRUE)
  if (!length(delta_csvs)) .log_abort("No *_delta_links.csv files found in {delta_dir}")

  existing_filtered_manifest <- file.path(filtered_dir, "cache", "filtered_links_manifest.csv")
  existing_filtered_files <- list.files(filtered_dir, "_filtered_links_(up|down)\\.csv$", full.names = TRUE)
  if (!isTRUE(overwrite_filtered) && file.exists(existing_filtered_manifest) && length(existing_filtered_files) > 0L) {
    .log_inform("Using existing filtered delta links from {filtered_dir} (set overwrite_filtered=TRUE to rebuild).")
  } else {
    filter_links_deltas_bulk(
      delta_csvs,
      gene_expr_min = threshold_gene_expr,
      tf_expr_min = threshold_tf_expr,
      fp_min = threshold_fp_score,
      link_min = threshold_link_score,
      abs_delta_min = delta_link,
      apply_de_gene = TRUE,
      de_gene_log2_abs_min = de_gene_log2_abs_min,
      de_tf_log2_abs_min = de_tf_log2_abs_min,
      fp_delta_min = fp_delta_min,
      split_direction = TRUE,
      write_combined = FALSE,
      enforce_link_expr_sign = TRUE,
      expr_dir_col = "log2FC_gene_expr",
      workers = workers,
      filtered_dir = filtered_dir
    )
  }

  if (isTRUE(qc)) {
    .plot_fp_gene_volcanoes_from_filtered(
      filtered_dir = filtered_dir,
      out_root = out_root,
      output_format = summary_plot_format,
      html_selfcontained = summary_html_selfcontained
    )
  }

  if (isTRUE(write_tf_hubs_fp)) {
    tf_hub_dir <- file.path(out_root, "master_tf_summary")
    dir.create(tf_hub_dir, recursive = TRUE, showWarnings = FALSE)
    filtered_files <- list.files(filtered_dir, "_filtered_links_(up|down)\\.csv$", full.names = TRUE)
    if (length(filtered_files)) {
      base_ids <- unique(sub("_(up|down)\\.csv$", "", basename(filtered_files)))
      for (bid in base_ids) {
        files <- filtered_files[grepl(paste0("^", gsub("([\\+\\-\\(\\)\\[\\]\\.\\^\\$\\|\\?\\*\\+])", "\\\\\\1", bid), "_(up|down)\\.csv$"), basename(filtered_files))]
        if (!length(files)) next
        df_list <- lapply(files, function(p) {
          x <- readr::read_csv(p, show_col_types = FALSE)
          if (grepl("_up\\.csv$", basename(p), ignore.case = TRUE)) {
            x$direction_group <- "Up"
          } else if (grepl("_down\\.csv$", basename(p), ignore.case = TRUE)) {
            x$direction_group <- "Down"
          }
          x
        })
        df_all <- dplyr::bind_rows(df_list)
        if (!nrow(df_all)) next
        contrast_id <- sub("_filtered_links.*$", "", bid)
        cp <- .contrast_parts(contrast_id)
        summary_df <- .summarize_delta_fp_for_tf_hubs(df_all, cond1 = cp$cond1, cond2 = cp$cond2)
        if (!nrow(summary_df)) next
        tag <- if (is.null(summary_tag) || !nzchar(summary_tag)) "fp" else summary_tag
        out_csv <- file.path(tf_hub_dir, paste0(contrast_id, "_master_tf_summary.csv"))
        out_pdf <- file.path(tf_hub_dir, paste0(contrast_id, "_master_tf_summary.pdf"))
        out_pdf_waterfall <- file.path(tf_hub_dir, paste0(contrast_id, "_tf_target_direction_waterfall.pdf"))
        out_pdf_heatmap <- file.path(tf_hub_dir, paste0(contrast_id, "_tf_tf_connectivity_heatmap.pdf"))
        if (!isTRUE(overwrite_tf_hubs) &&
            file.exists(out_csv) &&
            file.exists(out_pdf) &&
            file.exists(out_pdf_waterfall) &&
            file.exists(out_pdf_heatmap)) {
          .log_inform("Skipping existing master_tf_summary for {contrast_id} (set overwrite_tf_hubs=TRUE to rebuild).")
          next
        }
        readr::write_csv(summary_df, out_csv)
        contrast <- .contrast_from_file(out_csv)
        .plot_tf_hubs_fp(
          summary_df = summary_df,
          out_pdf = out_pdf,
          links_df = df_all,
          title_text = paste0("TF hubs (delta fp_score) - ", contrast),
          cond1_label = cp$cond1,
          cond2_label = cp$cond2,
          min_y_for_label = 3,
          size_min = 1,
          size_max = 4,
          connectivity_min_degree = connectivity_min_degree
        )
      }
    }
  }

  invisible(list(
    delta_dir = delta_dir,
    filtered_dir = filtered_dir
  ))
}

.plot_fp_gene_volcanoes_from_filtered <- function(filtered_dir,
                                                  out_root,
                                                  output_format = c("pdf", "html", "both"),
                                                  html_selfcontained = TRUE) {
  .assert_pkg("ggplot2")
  .assert_pkg("readr")
  output_format <- match.arg(output_format)
  do_pdf <- output_format %in% c("pdf", "both")
  do_html <- output_format %in% c("html", "both")
  if (isTRUE(do_html) && !requireNamespace("plotly", quietly = TRUE)) {
    if (identical(output_format, "html")) {
      .log_abort("summary_plot_format='html' requires the `plotly` package.")
    } else {
      .log_warn("plotly not installed; skipping HTML differential link summaries.")
      do_html <- FALSE
    }
  }
  if (isTRUE(do_html)) .assert_pkg("htmlwidgets")

  files <- list.files(filtered_dir, "_filtered_links_(up|down)\\.csv$", full.names = TRUE)
  if (!length(files)) return(invisible(NULL))
  base_ids <- unique(sub("_(up|down)\\.csv$", "", basename(files)))
  out_dir <- filtered_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_list <- list()
  all_x <- numeric(0)
  all_y <- numeric(0)
  all_gene_count_rows <- list()

  .pick_gene_col <- function(nm) {
    cands <- c("gene_key", "gene", "target_gene")
    hit <- cands[cands %in% nm]
    if (!length(hit)) return(NA_character_)
    hit[[1]]
  }

  .count_unique_genes_by_dir <- function(df, stage, comparison_id, direction_hint = NULL) {
    if (!is.data.frame(df) || !nrow(df)) return(data.frame())
    gene_col <- .pick_gene_col(names(df))
    if (is.na(gene_col)) return(data.frame())
    gvals <- as.character(df[[gene_col]])
    gvals[!is.na(gvals) & nzchar(gvals)]
    dir_vals <- rep(NA_character_, length(gvals))
    if (!is.null(direction_hint)) {
      dir_vals <- rep(as.character(direction_hint), length(gvals))
    } else if ("direction_group" %in% names(df)) {
      dir_vals <- as.character(df$direction_group)
    } else if ("direction" %in% names(df)) {
      dir_vals <- as.character(df$direction)
    } else if ("log2FC_gene_expr" %in% names(df)) {
      lfc <- suppressWarnings(as.numeric(df$log2FC_gene_expr))
      dir_vals <- ifelse(is.finite(lfc) & lfc >= 0, "Up", ifelse(is.finite(lfc) & lfc < 0, "Down", NA_character_))
    }
    keep <- !is.na(gvals) & nzchar(gvals) & !is.na(dir_vals) & dir_vals %in% c("Up", "Down")
    if (!any(keep)) return(data.frame())
    tmp <- data.frame(
      comparison_id = comparison_id,
      direction = dir_vals[keep],
      gene_key = gvals[keep],
      stage = stage,
      stringsAsFactors = FALSE
    )
    out <- tmp |>
      dplyr::distinct(.data$comparison_id, .data$direction, .data$stage, .data$gene_key) |>
      dplyr::count(.data$comparison_id, .data$direction, .data$stage, name = "n_gene_unique")
    as.data.frame(out)
  }

  for (bid in base_ids) {
    sel <- files[grepl(paste0("^", gsub("([\\+\\-\\(\\)\\[\\]\\.\\^\\$\\|\\?\\*\\+])", "\\\\\\1", bid), "_(up|down)\\.csv$"), basename(files))]
    if (!length(sel)) next
    contrast_id <- sub("_filtered_links$", "", bid)
    df_list <- lapply(sel, function(p) readr::read_csv(p, show_col_types = FALSE))
    df <- dplyr::bind_rows(df_list)
    if (!nrow(df)) next
    if (!all(c("delta_fp_bed_score", "log2FC_gene_expr") %in% names(df))) {
      .log_warn("Skipping volcano; missing delta_fp_bed_score or log2FC_gene_expr in {bid}.")
      next
    }
    df$log2FC_gene_expr <- suppressWarnings(as.numeric(df$log2FC_gene_expr))
    df$delta_fp_bed_score <- suppressWarnings(as.numeric(df$delta_fp_bed_score))
    df <- df[is.finite(df$log2FC_gene_expr) & is.finite(df$delta_fp_bed_score), , drop = FALSE]
    if (!nrow(df)) next
    df$direction <- ifelse(df$log2FC_gene_expr >= 0, "Up", "Down")
    tf_col <- if ("tf" %in% names(df)) "tf" else if ("TF" %in% names(df)) "TF" else NULL
    gene_col <- if ("gene_key" %in% names(df)) "gene_key" else if ("gene" %in% names(df)) "gene" else if ("target_gene" %in% names(df)) "target_gene" else NULL
    if (!is.null(tf_col) && !is.null(gene_col)) {
      df$tf_target <- paste0(as.character(df[[tf_col]]), "::", as.character(df[[gene_col]]))
    } else {
      df$tf_target <- "NA::NA"
    }
    df$hover_text <- paste0(
      df$tf_target,
      "<br>gene log2FC: ", signif(df$log2FC_gene_expr, 4),
      "<br>delta FP: ", signif(df$delta_fp_bed_score, 4),
      "<br>Direction: ", df$direction
    )
    all_x <- c(all_x, df$log2FC_gene_expr)
    all_y <- c(all_y, df$delta_fp_bed_score)
    contrast <- .contrast_from_file(bid)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC_gene_expr, y = delta_fp_bed_score, color = direction, text = hover_text)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.4) +
      ggplot2::geom_point(alpha = 0.6, size = 0.7) +
      ggplot2::scale_color_manual(values = c(Up = "#d73027", Down = "#4575b4")) +
      ggplot2::labs(
        title = contrast,
        x = "gene log2FC",
        y = "delta FP",
        color = "Direction"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title.x = ggplot2::element_text(face = "bold", size = 8),
        axis.title.y = ggplot2::element_text(face = "bold", size = 8),
        axis.text = ggplot2::element_text(size = 7),
        legend.position = "none"
    )
    plot_list[[length(plot_list) + 1L]] <- p

    # Gene-key counts by direction: before = delta links, after = filtered links.
    count_rows_this <- list()
    for (p in sel) {
      dlab <- if (grepl("_up\\.csv$", basename(p), ignore.case = TRUE)) "Up" else "Down"
      dfi <- readr::read_csv(p, show_col_types = FALSE)
      cr <- .count_unique_genes_by_dir(
        df = dfi,
        stage = "After filtering",
        comparison_id = contrast_id,
        direction_hint = dlab
      )
      if (nrow(cr)) count_rows_this[[length(count_rows_this) + 1L]] <- cr
    }
    delta_path <- file.path(out_root, "diff_links", paste0(contrast_id, "_delta_links.csv"))
    if (!file.exists(delta_path)) {
      delta_path <- file.path(out_root, "delta_links", paste0(contrast_id, "_delta_links.csv"))
    }
    if (file.exists(delta_path)) {
      delta_df <- readr::read_csv(delta_path, show_col_types = FALSE)
      crb <- .count_unique_genes_by_dir(
        df = delta_df,
        stage = "Before filtering",
        comparison_id = contrast_id,
        direction_hint = NULL
      )
      if (nrow(crb)) count_rows_this[[length(count_rows_this) + 1L]] <- crb
    }
    gene_counts_this <- if (length(count_rows_this)) dplyr::bind_rows(count_rows_this) else data.frame()
    if (nrow(gene_counts_this)) {
      gene_counts_this$direction <- factor(gene_counts_this$direction, levels = c("Up", "Down"))
      gene_counts_this$stage <- factor(gene_counts_this$stage, levels = c("Before filtering", "After filtering"))
      all_gene_count_rows[[length(all_gene_count_rows) + 1L]] <- gene_counts_this
    }

    if (isTRUE(do_pdf)) {
      out_pdf <- file.path(out_dir, paste0(contrast_id, "_differential_links_summary.pdf"))
      grDevices::pdf(out_pdf, width = 6.2, height = 4.6, onefile = TRUE)
      print(p)
      if (nrow(gene_counts_this)) {
        p_cnt <- ggplot2::ggplot(
          gene_counts_this,
          ggplot2::aes(x = .data$direction, y = .data$n_gene_unique, fill = .data$stage)
        ) +
          ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
          ggplot2::scale_fill_manual(values = c("Before filtering" = "#9ecae1", "After filtering" = "#3182bd")) +
          ggplot2::labs(
            title = paste0(.contrast_from_file(bid), " | Unique gene_key counts"),
            x = "Direction",
            y = "Unique gene_key (N)",
            fill = NULL
          ) +
          ggplot2::theme_classic(base_size = 9) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
            axis.title = ggplot2::element_text(face = "bold"),
            legend.position = "top"
          )
        print(p_cnt)
      }
      grDevices::dev.off()
    }
    if (isTRUE(do_html)) {
      out_html <- file.path(out_dir, paste0(contrast_id, "_differential_links_summary.html"))
      wp <- plotly::ggplotly(p, tooltip = "text")
      htmlwidgets::saveWidget(wp, out_html, selfcontained = isTRUE(html_selfcontained))
    }
  }
  if (!length(plot_list) || !isTRUE(do_pdf)) return(invisible(TRUE))

  xlim <- range(all_x, finite = TRUE)
  ylim <- range(all_y, finite = TRUE)
  if (!all(is.finite(xlim))) xlim <- c(-1, 1)
  if (!all(is.finite(ylim))) ylim <- c(-1, 1)

  pdf_path <- file.path(out_dir, "differential_links_summary_plots.pdf")
  grDevices::pdf(pdf_path, width = 16, height = 9.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  n_per_page <- 12L
  total <- length(plot_list)
  page_idx <- seq_len(ceiling(total / n_per_page))
  for (pg in page_idx) {
    start <- (pg - 1L) * n_per_page + 1L
    end <- min(pg * n_per_page, total)
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 3, ncol = 4)
    grid::pushViewport(grid::viewport(layout = lay))
    idx <- 0L
    for (i in start:end) {
      idx <- idx + 1L
      row <- ((idx - 1L) %/% 4L) + 1L
      col <- ((idx - 1L) %% 4L) + 1L
      grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
      grid::grid.draw(ggplot2::ggplotGrob(
        plot_list[[i]] + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
      ))
      grid::popViewport()
    }
    grid::popViewport()
  }

  if (isTRUE(do_pdf) && length(all_gene_count_rows)) {
    counts_all <- dplyr::bind_rows(all_gene_count_rows)
    if (nrow(counts_all)) {
      counts_all$comparison_direction <- paste0(counts_all$comparison_id, " | ", counts_all$direction)
      counts_all$stage <- factor(counts_all$stage, levels = c("Before filtering", "After filtering"))
      ord_tbl <- counts_all |>
        dplyr::filter(.data$stage == "Before filtering") |>
        dplyr::arrange(dplyr::desc(.data$n_gene_unique)) |>
        dplyr::distinct(.data$comparison_direction, .keep_all = FALSE)
      ord_levels <- ord_tbl$comparison_direction
      if (!length(ord_levels)) ord_levels <- unique(counts_all$comparison_direction)
      counts_all$comparison_direction <- factor(counts_all$comparison_direction, levels = ord_levels)
      p_all_counts <- ggplot2::ggplot(
        counts_all,
        ggplot2::aes(x = .data$comparison_direction, y = .data$n_gene_unique, fill = .data$stage)
      ) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
        ggplot2::scale_fill_manual(values = c("Before filtering" = "#9ecae1", "After filtering" = "#3182bd")) +
        ggplot2::labs(
          title = "Unique gene_key counts by comparison-direction",
          x = "Comparison | Direction",
          y = "Unique gene_key (N)",
          fill = NULL
        ) +
        ggplot2::theme_classic(base_size = 9) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.title = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          legend.position = "top"
        )
      ggplot2::ggsave(
        filename = file.path(out_dir, "differential_links_gene_key_counts_summary.pdf"),
        plot = p_all_counts,
        width = 14,
        height = max(5, min(14, 2.8 + 0.18 * length(unique(counts_all$comparison_direction)))),
        dpi = 200,
        limitsize = FALSE
      )
    }
  }
  invisible(TRUE)
}

.contrast_from_file <- function(f) {
  b <- basename(f)
  stem <- sub("\\.(csv|pdf|html)$", "", b, ignore.case = TRUE)
  stem <- sub("_master_tf_summary.*$", "", stem)
  stem <- sub("_differential_links_summary.*$", "", stem)
  stem <- sub("_filtered_links.*$", "", stem)
  stem <- sub("_delta_links.*$", "", stem)
  parts <- strsplit(stem, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) == 2) paste(parts[1], "vs", parts[2]) else stem
}

.contrast_parts <- function(f) {
  b <- basename(f)
  stem <- sub("\\.(csv|pdf|html)$", "", b, ignore.case = TRUE)
  stem <- sub("_master_tf_summary.*$", "", stem)
  stem <- sub("_differential_links_summary.*$", "", stem)
  stem <- sub("_filtered_links.*$", "", stem)
  stem <- sub("_delta_links.*$", "", stem)
  parts <- strsplit(stem, "_vs_", fixed = TRUE)[[1]]
  cond1 <- if (length(parts) >= 1) parts[1] else NA_character_
  cond2 <- if (length(parts) >= 2) parts[2] else NA_character_
  list(cond1 = cond1, cond2 = cond2, label = if (length(parts) == 2) paste(cond1, "vs", cond2) else stem)
}

.path_add_postfix <- function(path, postfix, new_ext = ".pdf") {
  noext <- sub("\\.[^.]+$", "", path)
  paste0(noext, ".", postfix, new_ext)
}

.summarize_delta_fp_for_tf_hubs <- function(df, cond1 = NULL, cond2 = NULL) {
  tf_col <- if ("tf" %in% names(df)) "tf" else if ("TF" %in% names(df)) "TF" else NULL
  gene_col <- if ("gene_key" %in% names(df)) "gene_key" else if ("gene" %in% names(df)) "gene" else if ("target_gene" %in% names(df)) "target_gene" else NULL
  if (is.null(tf_col) || is.null(gene_col)) return(tibble::tibble())

  delta_col <- if ("delta_fp_bed_score" %in% names(df)) "delta_fp_bed_score" else if ("delta_fp" %in% names(df)) "delta_fp" else NULL
  if (is.null(delta_col)) return(tibble::tibble())

  tf_expr_max <- NA_real_
  if (is.character(cond1) && nzchar(cond1) && is.character(cond2) && nzchar(cond2)) {
    tf_c1 <- paste0("tf_expr_", cond1)
    tf_c2 <- paste0("tf_expr_", cond2)
    if (all(c(tf_c1, tf_c2) %in% names(df))) {
      tf_expr_max <- pmax(
        suppressWarnings(as.numeric(df[[tf_c1]])),
        suppressWarnings(as.numeric(df[[tf_c2]])),
        na.rm = TRUE
      )
    }
  }
  if (all(is.na(tf_expr_max)) && "tf_expr_max" %in% names(df)) {
    tf_expr_max <- suppressWarnings(as.numeric(df$tf_expr_max))
  }
  if (all(is.na(tf_expr_max)) && all(c("tf_expr_ctrl", "tf_expr_str") %in% names(df))) {
    tf_expr_max <- pmax(
      suppressWarnings(as.numeric(df$tf_expr_ctrl)),
      suppressWarnings(as.numeric(df$tf_expr_str)),
      na.rm = TRUE
    )
  }
  tf_log2_fc <- if ("log2FC_tf_expr" %in% names(df)) {
    suppressWarnings(as.numeric(df$log2FC_tf_expr))
  } else if ("delta_tf_expr" %in% names(df)) {
    suppressWarnings(as.numeric(df$delta_tf_expr))
  } else if ("tf_log2_fc" %in% names(df)) {
    suppressWarnings(as.numeric(df$tf_log2_fc))
  } else {
    NA_real_
  }

  delta_vec <- suppressWarnings(as.numeric(df[[delta_col]]))
  reg_sign <- sign(delta_vec)
  delta_oriented <- delta_vec
  delta_abs <- abs(delta_vec)

  df0 <- tibble::tibble(
    TF = df[[tf_col]],
    gene = df[[gene_col]],
    delta = delta_vec,
    delta_oriented = delta_oriented,
    delta_abs = delta_abs,
    reg_sign = reg_sign,
    tf_expr_max = tf_expr_max,
    tf_log2_fc = tf_log2_fc
  )
  df0 <- df0[is.finite(df0$delta), , drop = FALSE]
  if (!nrow(df0)) return(tibble::tibble())

  df0 |>
    dplyr::group_by(.data$TF) |>
    dplyr::summarise(
      tf_delta_sum                  = sum(.data$delta_oriented, na.rm = TRUE),
      tf_delta_sum_activate         = sum(ifelse(.data$reg_sign > 0,  .data$delta_oriented, 0), na.rm = TRUE),
      tf_delta_sum_repress          = sum(ifelse(.data$reg_sign < 0,  .data$delta_oriented, 0), na.rm = TRUE),
      tf_delta_sum_abs              = sum(.data$delta_abs, na.rm = TRUE),
      tf_delta_sum_abs_activate     = sum(ifelse(.data$reg_sign > 0,  .data$delta_abs, 0), na.rm = TRUE),
      tf_delta_sum_abs_repress      = sum(ifelse(.data$reg_sign < 0,  .data$delta_abs, 0), na.rm = TRUE),
      tf_n_links                    = dplyr::n(),
      tf_n_links_activate           = sum(.data$reg_sign > 0, na.rm = TRUE),
      tf_n_links_repress            = sum(.data$reg_sign < 0, na.rm = TRUE),
      tf_expr_max                   = suppressWarnings(max(.data$tf_expr_max, na.rm = TRUE)),
      tf_log2_fc                    = stats::median(.data$tf_log2_fc, na.rm = TRUE),
      tf_hits_hub                   = NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      topic = 1L,
      topic_mean_abs_delta = mean(.data$tf_delta_sum_abs, na.rm = TRUE),
      topic_n_TFs = dplyr::n(),
      topic_n_links = sum(.data$tf_n_links, na.rm = TRUE),
      topic_rank = 1L
    )
}

.robust_z <- function(x) {
  x <- as.numeric(x)
  m <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(madv) || madv == 0) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
    return((x - m) / (sx + 1e-12))
  }
  (x - m) / (madv + 1e-12)
}

.clamp <- function(x, lo, hi) pmax(pmin(x, hi), lo)

.plot_tf_hubs_fp <- function(summary_df,
                             out_pdf,
                             links_df = NULL,
                             title_text = NULL,
                             cond1_label = NULL,
                             cond2_label = NULL,
                             min_y_for_label = 3,
                             size_min = 1,
                             size_max = 4,
                             color_sigma = 2,
                             base_size = 12,
                             width_in = 10,
                             height_in = 6.5,
                             dpi = 200,
                             connectivity_min_degree = 1L) {
  .assert_pkg("ggplot2")
  .assert_pkg("ggrepel")
  .assert_pkg("scales")
  .assert_pkg("pheatmap")
  .assert_pkg("RColorBrewer")

  need <- c(
    "TF","tf_expr_max","tf_log2_fc",
    "tf_n_links","tf_delta_sum","tf_delta_sum_abs"
  )
  if (!all(need %in% names(summary_df))) {
    .log_abort("TF hub summary missing required columns.")
  }

  c1 <- if (is.null(cond1_label) || !nzchar(cond1_label)) "cond1" else cond1_label
  c2 <- if (is.null(cond2_label) || !nzchar(cond2_label)) "cond2" else cond2_label

  TF <- summary_df |>
    dplyr::group_by(.data$TF) |>
    dplyr::summarise(
      tf_links = sum(.data$tf_n_links, na.rm = TRUE),
      tf_sum_delta = sum(.data$tf_delta_sum, na.rm = TRUE),
      tf_sum_abs_delta = sum(.data$tf_delta_sum_abs, na.rm = TRUE),
      tf_expr_max = max(.data$tf_expr_max, na.rm = TRUE),
      tf_log2_fc_med = stats::median(.data$tf_log2_fc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      tf_links_plus1 = tf_links + 1L,
      tf_log2_expr_max = log2(pmax(tf_expr_max, 1e-9))
    )

  TF$x_sum <- sign(TF$tf_sum_delta) * log2(abs(TF$tf_sum_delta) + 1)
  TF$log2_fc_z_c <- .clamp(.robust_z(TF$tf_log2_fc_med), -abs(color_sigma), abs(color_sigma))

  lab_df <- TF |>
    dplyr::filter(tf_links >= min_y_for_label)

  y_breaks <- 2^(0:ceiling(log(max(TF$tf_links_plus1, na.rm = TRUE), base = 2)))
  y_labels <- pmax(y_breaks - 1, 0)

  if (is.null(title_text)) {
    title_text <- "TF hubs"
  }

  caption_text <- paste0(
    "delta fp_score = fp_score(", c1, ") - fp_score(", c2, ").\n",
    "Y: log2(1 + #links). Size: log2(max TF RNA across ", c1, " & ", c2, "). Color: TF log2FC z-score."
  )

  p <- ggplot2::ggplot(TF, ggplot2::aes(x = x_sum, y = tf_links_plus1)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(size = tf_log2_expr_max, color = log2_fc_z_c),
                        shape = 16, alpha = 0.8) +
    ggrepel::geom_text_repel(
      data = lab_df, ggplot2::aes(label = TF),
      size = 1.8, fontface = "bold",
      max.overlaps = 100, box.padding = 0.15, point.padding = 0.15,
      min.segment.length = 0.1, segment.alpha = 0.4
    ) +
    ggplot2::scale_size_continuous(range = c(size_min, size_max), name = "expr max (log2)") +
    ggplot2::scale_color_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      midpoint = 0, limits = c(-abs(color_sigma), abs(color_sigma)),
      oob = scales::squish, name = "TF log2FC (z)"
    ) +
    ggplot2::labs(
      title = title_text,
      x = "sum delta fp_score per TF (log2-scaled magnitude)",
      y = "number of differential links per TF in [log2(1 + count)] scale",
      caption = caption_text
    ) +
    ggplot2::scale_y_continuous(
      trans = scales::log_trans(base = 2),
      breaks = y_breaks, labels = y_labels,
      expand = ggplot2::expansion(mult = c(0.02, 0.06))
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1)),
      color = ggplot2::guide_colorbar(order = 2)
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(face = "bold")
    )

  out_base <- sub("\\.pdf$", "", basename(out_pdf))
  out_base <- sub("_master_tf_summary$", "", out_base)
  out_pdf_waterfall <- file.path(dirname(out_pdf), paste0(out_base, "_tf_target_direction_waterfall.pdf"))
  out_pdf_heatmap <- file.path(dirname(out_pdf), paste0(out_base, "_tf_tf_connectivity_heatmap.pdf"))

  # Master summary: page 1 only.
  ggplot2::ggsave(out_pdf, p, width = width_in, height = height_in, dpi = dpi, limitsize = FALSE)

  # Additional pages use the original per-link table when available.
  if (is.data.frame(links_df) && nrow(links_df)) {
    link_dt <- data.table::as.data.table(links_df)
    tf_col <- if ("tf" %in% names(link_dt)) "tf" else if ("TF" %in% names(link_dt)) "TF" else NULL
    gene_col <- if ("gene_key" %in% names(link_dt)) "gene_key" else if ("gene" %in% names(link_dt)) "gene" else if ("target_gene" %in% names(link_dt)) "target_gene" else NULL
    peak_col <- if ("peak_id" %in% names(link_dt)) "peak_id" else NULL

    if (!is.null(tf_col) && !is.null(gene_col)) {
      # Normalize direction.
      if ("direction_group" %in% names(link_dt)) {
        link_dt[, .dir := tolower(trimws(as.character(direction_group)))]
        link_dt[grepl("^up", .dir), .dir2 := "Up"]
        link_dt[grepl("^down", .dir), .dir2 := "Down"]
      } else if ("log2FC_gene_expr" %in% names(link_dt)) {
        link_dt[, .dir2 := ifelse(suppressWarnings(as.numeric(log2FC_gene_expr)) >= 0, "Up", "Down")]
      } else if ("delta_gene_expr" %in% names(link_dt)) {
        link_dt[, .dir2 := ifelse(suppressWarnings(as.numeric(delta_gene_expr)) >= 0, "Up", "Down")]
      } else {
        link_dt[, .dir2 := "Up"]
      }
      link_dt <- link_dt[.dir2 %in% c("Up", "Down")]

      # De-duplicate links once.
      if (!is.null(peak_col)) {
        link_dt <- unique(link_dt[, .(TF = as.character(get(tf_col)), gene_key = as.character(get(gene_col)), peak_id = as.character(get(peak_col)), direction = .dir2)])
      } else {
        link_dt <- unique(link_dt[, .(TF = as.character(get(tf_col)), gene_key = as.character(get(gene_col)), direction = .dir2)])
      }
      link_dt <- link_dt[!is.na(TF) & nzchar(TF) & !is.na(gene_key) & nzchar(gene_key)]

      if (nrow(link_dt)) {
        tf_universe <- unique(link_dt$TF)
        tf_upper <- toupper(tf_universe)
        link_dt[, target_type := ifelse(toupper(gene_key) %in% tf_upper, "TF target", "Gene target")]

        # Waterfall uses unique gene counts (not link counts), with TF-level display cutoff.
        gene_dt <- unique(link_dt[, .(TF, gene_key, direction, target_type)])
        tf_gene_totals <- gene_dt[, .(n_unique_genes = data.table::uniqueN(gene_key)), by = TF]
        tf_keep <- tf_gene_totals[n_unique_genes > 10, TF]
        gene_dt_plot <- gene_dt[TF %in% tf_keep]

        # Sort TFs from largest net positive to largest net negative by unique-gene counts.
        tf_order <- gene_dt_plot[
          , .(net_n = sum(ifelse(direction == "Up", 1L, -1L), na.rm = TRUE)),
          by = TF
        ][order(net_n, decreasing = TRUE)]$TF
        bar_dt_all <- gene_dt_plot[
          , .(n = ifelse(.BY$direction == "Up", .N, -.N)),
          by = .(TF, target_type, direction)
        ]

        # Stacked waterfall data prep by TF.
        x_max_abs <- suppressWarnings(max(abs(bar_dt_all$n), na.rm = TRUE))
        if (!is.finite(x_max_abs) || x_max_abs <= 0) x_max_abs <- 1
        x_lim <- c(-x_max_abs, x_max_abs)
        tf_per_page <- 45L
        tf_chunks <- split(tf_order, ceiling(seq_along(tf_order) / tf_per_page))
        n_pages <- length(tf_chunks)

        # Standalone waterfall PDF as a single page with dynamic height and fixed x-axis.
        if (length(tf_order) > 0L) {
          n_tf_total <- length(tf_order)
          tf_axis_size_single <- max(4, min(8, 260 / max(10, n_tf_total)))
          waterfall_height <- max(height_in, min(28, 2.8 + 0.11 * n_tf_total))
          bar_dt_single <- data.table::copy(bar_dt_all)
          bar_dt_single[, TF := factor(TF, levels = rev(tf_order))]
          p2_single <- ggplot2::ggplot(bar_dt_single, ggplot2::aes(x = n, y = TF, fill = direction, alpha = target_type)) +
            ggplot2::geom_col(width = 0.82) +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
            ggplot2::scale_fill_manual(values = c("Up" = "#C21807", "Down" = "#1565C0")) +
            ggplot2::scale_alpha_manual(values = c("TF target" = 1.0, "Gene target" = 0.45)) +
            ggplot2::labs(
              title = paste0(
                sub("^TF hubs \\(delta fp_score\\) - ", "", title_text),
                " | TF target-gene waterfall"
              ),
              x = "Signed unique-gene count (Up positive, Down negative)",
              y = "TF",
              fill = "Direction",
              alpha = "Target class",
              caption = "Counts are unique genes per TF. TFs shown: >10 unique genes. Red=Up, Blue=Down; shading encodes target class."
            ) +
            ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
            ggplot2::theme_classic(base_size = base_size) +
            ggplot2::theme(
              text = ggplot2::element_text(face = "bold"),
              plot.title = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(b = 10), face = "bold"),
              axis.title.x = ggplot2::element_text(face = "bold"),
              axis.title.y = ggplot2::element_text(face = "bold"),
              axis.text.y = ggplot2::element_text(size = tf_axis_size_single, face = "bold"),
              axis.text.x = ggplot2::element_text(face = "bold"),
              legend.position = "right",
              plot.caption = ggplot2::element_text(size = max(7, base_size - 3), margin = ggplot2::margin(t = 8)),
              plot.margin = ggplot2::margin(t = 16, r = 40, b = 22, l = 14)
            )
          ggplot2::ggsave(out_pdf_waterfall, p2_single, width = width_in, height = waterfall_height, dpi = dpi, limitsize = FALSE)
        } else {
          grDevices::pdf(out_pdf_waterfall, width = width_in, height = height_in, onefile = TRUE)
          graphics::plot.new()
          graphics::title(main = paste0(title_text, " | TF target-direction waterfall"))
          graphics::text(0.5, 0.5, "No TFs pass cutoff: unique gene count > 10")
          grDevices::dev.off()
        }

        # Page 3: TF-to-TF connectivity clustered heatmap.
        tf2tf <- link_dt[toupper(gene_key) %in% tf_upper, .(TF_src = TF, TF_tgt = gene_key)]
        if (nrow(tf2tf)) {
          conn <- tf2tf[, .N, by = .(TF_src, TF_tgt)]
          tf_levels <- sort(unique(c(conn$TF_src, conn$TF_tgt)))
          mat <- matrix(0, nrow = length(tf_levels), ncol = length(tf_levels),
                        dimnames = list(tf_levels, tf_levels))
          for (ii in seq_len(nrow(conn))) {
            mat[conn$TF_src[ii], conn$TF_tgt[ii]] <- conn$N[ii]
          }

          deg_cut <- suppressWarnings(as.integer(connectivity_min_degree))
          if (!is.finite(deg_cut) || is.na(deg_cut) || deg_cut < 0L) deg_cut <- 1L
          adj_undir_keep <- (mat > 0) | (t(mat) > 0)
          diag(adj_undir_keep) <- FALSE
          deg_vec <- rowSums(adj_undir_keep, na.rm = TRUE)
          tf_keep <- names(deg_vec)[deg_vec >= deg_cut]
          if (length(tf_keep) >= 2L) {
            mat <- mat[tf_keep, tf_keep, drop = FALSE]
          } else {
            tryCatch({
              grDevices::pdf(out_pdf_heatmap, width = 10, height = 10, onefile = TRUE)
              pdf_dev <- grDevices::dev.cur()
              graphics::plot.new()
              graphics::title(main = paste0(sub("^TF hubs \\(delta fp_score\\) - ", "", title_text), " | TF-to-TF connectivity heatmap"))
              graphics::text(0.5, 0.58, "No TFs pass connectivity cutoff")
              graphics::text(0.5, 0.50, paste0("Cutoff: undirected degree >= ", deg_cut))
              graphics::text(0.5, 0.42, paste0("TFs passing cutoff: ", length(tf_keep)))
              try(grDevices::dev.off(pdf_dev), silent = TRUE)
            }, error = function(e) {
              .log_warn("Failed writing cutoff-empty standalone TF connectivity heatmap PDF: {conditionMessage(e)}")
            })
            mat <- NULL
          }

          if (!is.null(mat) && is.matrix(mat) && nrow(mat) >= 2L) {
            mat_log <- log1p(mat)
            tf_fs <- max(6.5, min(13, 320 / max(10, nrow(mat_log))))
            cell_sz <- max(2, min(12, 240 / max(12, nrow(mat_log))))
            # Build an undirected shortest-path connectivity score for indirect links.
            # Score(i,j) = 1 / (1 + number of intermediate TF nodes on shortest path i->j),
            # so direct links have score 1 and disconnected pairs have score 0.
            adj_undir <- (mat > 0) | (t(mat) > 0)
            diag(adj_undir) <- FALSE
            n_tf <- nrow(adj_undir)
            dist_mat <- matrix(Inf, n_tf, n_tf, dimnames = dimnames(adj_undir))
            for (src in seq_len(n_tf)) {
              d <- rep(Inf, n_tf)
              d[src] <- 0
              q <- src
              while (length(q)) {
                v <- q[[1]]
                q <- q[-1]
                nb <- which(adj_undir[v, ] & is.infinite(d))
                if (length(nb)) {
                  d[nb] <- d[v] + 1
                  q <- c(q, nb)
                }
              }
              dist_mat[src, ] <- d
            }
            nodes_needed <- pmax(dist_mat - 1, 0)
            conn_indirect <- matrix(0, n_tf, n_tf, dimnames = dimnames(dist_mat))
            conn_indirect[is.finite(nodes_needed)] <- 1 / (1 + nodes_needed[is.finite(nodes_needed)])
            diag(conn_indirect) <- 1

            # Third connectivity mode (one-layer composite):
            # score(i,j) = I(i->j) + I(j->i) + 0.5 * |targets(i) ∩ targets(j)|
            # where targets(.) are gene_key targets (including TF targets).
            tf_levels_now <- rownames(mat)
            edge_gene <- unique(link_dt[, .(TF = as.character(TF), gene_key = as.character(gene_key))])
            edge_gene <- edge_gene[TF %in% tf_levels_now]
            if (nrow(edge_gene)) {
              gene_levels <- sort(unique(edge_gene$gene_key))
              inc <- matrix(0L, nrow = length(tf_levels_now), ncol = length(gene_levels),
                            dimnames = list(tf_levels_now, gene_levels))
              r_idx <- match(edge_gene$TF, tf_levels_now)
              c_idx <- match(edge_gene$gene_key, gene_levels)
              keep_idx <- which(!is.na(r_idx) & !is.na(c_idx))
              if (length(keep_idx)) {
                inc[cbind(r_idx[keep_idx], c_idx[keep_idx])] <- 1L
              }
              shared_targets <- inc %*% t(inc)
            } else {
              shared_targets <- matrix(0, nrow = n_tf, ncol = n_tf, dimnames = dimnames(mat))
            }
            direct_bin <- (mat > 0) * 1L
            conn_layer1 <- direct_bin + t(direct_bin) + 0.5 * shared_targets
            conn_layer1 <- as.matrix(conn_layer1)
            diag(conn_layer1) <- 0

            # Standalone square heatmap PDF with two pages: direct and indirect.
            # Use ggplot tiles on a dedicated device to avoid no-page PDF issues.
            tryCatch({
              old_dev <- grDevices::dev.cur()
              tf_n <- nrow(mat_log)
              side_in <- max(10, min(24, 4 + 0.15 * tf_n))
              grDevices::pdf(out_pdf_heatmap, width = side_in, height = side_in, onefile = TRUE)
              pdf_dev <- grDevices::dev.cur()
              grDevices::dev.set(pdf_dev)

            .plot_mat_tile <- function(mat_in, title_in, fill_title, palette_fn, fill_transform = "none") {
              ord <- stats::hclust(stats::dist(mat_in), method = "complete")$order
              mat_ord <- mat_in[ord, ord, drop = FALSE]
              dt <- data.table::as.data.table(as.table(mat_ord))
              names(dt) <- c("TF_row", "TF_col", "value")
              dt[, value_raw := as.numeric(value)]
              if (identical(fill_transform, "log1p")) {
                dt[, value_fill := log1p(pmax(value_raw, 0))]
                max_raw <- max(dt$value_raw, na.rm = TRUE)
                raw_breaks <- pretty(c(0, max_raw), n = 5)
                raw_breaks <- raw_breaks[is.finite(raw_breaks) & raw_breaks >= 0]
                fill_breaks <- log1p(raw_breaks)
                fill_labels <- format(round(raw_breaks, 2), trim = TRUE, scientific = FALSE)
              } else {
                dt[, value_fill := value_raw]
                fill_breaks <- ggplot2::waiver()
                fill_labels <- ggplot2::waiver()
              }
              dt[, TF_row := factor(as.character(TF_row), levels = rev(rownames(mat_ord)))]
              dt[, TF_col := factor(as.character(TF_col), levels = colnames(mat_ord))]
              fs <- max(6.5, min(13, 340 / max(10, nrow(mat_ord))))
              p_hm <- ggplot2::ggplot(dt, ggplot2::aes(x = TF_col, y = TF_row, fill = value_fill)) +
                ggplot2::geom_tile() +
                ggplot2::scale_fill_gradientn(
                  colors = palette_fn(128),
                  name = fill_title,
                  breaks = fill_breaks,
                  labels = fill_labels
                ) +
                ggplot2::coord_fixed() +
                ggplot2::labs(
                  title = title_in,
                  x = "TF target",
                  y = "TF source",
                  caption = paste0("TF inclusion cutoff: undirected degree >= ", deg_cut)
                ) +
                ggplot2::theme_minimal(base_size = 11) +
                ggplot2::theme(
                  text = ggplot2::element_text(face = "bold"),
                  plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                  axis.title.x = ggplot2::element_text(face = "bold"),
                  axis.title.y = ggplot2::element_text(face = "bold"),
                  axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = fs, face = "bold"),
                  axis.text.y = ggplot2::element_text(size = fs, face = "bold"),
                  plot.caption = ggplot2::element_text(face = "bold", hjust = 0, size = max(8, fs * 0.9)),
                  panel.grid = ggplot2::element_blank()
                )
              print(p_hm)
            }

            heat_cols <- grDevices::colorRampPalette(
              rev(RColorBrewer::brewer.pal(7, "RdYlBu"))
            )(100)

            .plot_mat_tile(
              conn_indirect,
              paste0(
                sub("^TF hubs \\(delta fp_score\\) - ", "", title_text),
                " | TF-to-TF indirect connectivity (clustered)\nScore(i,j)=1/(1 + # min nodes); direct=1, disconnected=0."
              ),
              "indirect score",
              function(n) heat_cols,
              fill_transform = "none"
            )
            .plot_mat_tile(
              conn_layer1,
              paste0(
                sub("^TF hubs \\(delta fp_score\\) - ", "", title_text),
                " | TF-to-TF one-layer composite connectivity (clustered)\nScore(i,j)=I(i->j)+I(j->i)+0.5*#shared targets (color uses log1p(score))"
              ),
              "composite score",
              function(n) heat_cols,
              fill_transform = "log1p"
            )
              # Close standalone PDF immediately and return to the original device.
              try(grDevices::dev.off(pdf_dev), silent = TRUE)
              devs2 <- grDevices::dev.list()
              if (!is.null(devs2) && old_dev %in% devs2) {
                try(grDevices::dev.set(old_dev), silent = TRUE)
              }
            }, error = function(e) {
              .log_warn("Failed writing standalone TF connectivity heatmap PDF: {conditionMessage(e)}")
            })
          }
        } else {
          # Ensure standalone connectivity PDF still has one explanatory page.
          tryCatch({
            grDevices::pdf(out_pdf_heatmap, width = 10, height = 10, onefile = TRUE)
            pdf_dev <- grDevices::dev.cur()
            graphics::plot.new()
            graphics::title(main = paste0(title_text, " | TF-to-TF connectivity heatmap"))
            graphics::text(0.5, 0.5, "No TF->TF links found")
            try(grDevices::dev.off(pdf_dev), silent = TRUE)
          }, error = function(e) {
            .log_warn("Failed writing empty standalone TF connectivity heatmap PDF: {conditionMessage(e)}")
          })
        }
      }
    }
  }

  invisible(out_pdf)
}

# ---------------------------------------------------------------------------
# Examples (not run)
# ---------------------------------------------------------------------------
#' @title Run Pathway Enrichment and Optional Pathway Subnet Plots
#' @description Minimal helper for `diff_links_filtered`: runs pathway enrichment
#'   from distinct `gene_key` values per comparison-direction and optionally
#'   renders pathway subnet HTML plots using existing network plotting utilities.
#' Minimal pathway enrichment on filtered differential-link files
#'
#' For each `*_filtered_links_up.csv` / `*_filtered_links_down.csv` under
#' `diff_links_filtered`, reads distinct `gene_key`, runs Enrichr pathway
#' enrichment (same DB defaults as later topic-pathway steps), sorts by adjusted
#' p-value (small to large), and writes a CSV in the same folder.
#'
#' @param diff_res Optional list returned by [find_differential_links()].
#' @param filtered_dir Directory containing filtered link files.
#' @param dbs Enrichr databases.
#' @param min_genes Minimum distinct genes to run enrichment.
#' @param padj_cut Keep rows with adjusted p-value <= this cutoff.
#' @param plot_subnetwork If TRUE, render pathway subnet HTMLs using
#'   \code{plot_tf_network_delta()}.
#' @param top_n_pathways_plot Maximum pathways (smallest adjusted p-value) to
#'   plot per comparison-direction.
#' @param pathway_gene_overlap_thresh Overlap threshold (relative to smaller
#'   gene set) used to treat pathways as redundant before plotting.
#' @param subnetwork_dirname Output subfolder (under \code{filtered_dir}) for
#'   pathway subnetworks.
#' @param html_selfcontained Passed to \code{htmlwidgets::saveWidget()} for
#'   pathway subnet HTML.
#' @param overwrite Overwrite existing enrichment CSVs.
#' @param verbose Emit concise messages.
#'
#' @return Invisible character vector of output CSV paths.
#' @export
run_diff_links_pathway_grn <- function(diff_res = NULL,
                                       filtered_dir = NULL,
                                       dbs = c("GO_Biological_Process_2023",
                                               "GO_Cellular_Component_2023",
                                               "GO_Molecular_Function_2023",
                                               "Reactome_2022",
                                               "WikiPathways_2024_Human"),
                                       min_genes = 5L,
                                       padj_cut = 0.05,
                                       plot_subnetwork = TRUE,
                                       top_n_pathways_plot = 10L,
                                       pathway_gene_overlap_thresh = 0.8,
                                       subnetwork_dirname = "subnet_pathway",
                                       html_selfcontained = FALSE,
                                       overwrite = FALSE,
                                       verbose = TRUE) {
  if (is.null(filtered_dir) && is.list(diff_res) && !is.null(diff_res$filtered_dir)) {
    filtered_dir <- diff_res$filtered_dir
  }
  if (!is.character(filtered_dir) || !nzchar(filtered_dir) || !dir.exists(filtered_dir)) {
    .log_abort("`filtered_dir` must be an existing directory (or provide `diff_res$filtered_dir`).")
  }
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    .log_abort("`enrichR` package is required for pathway enrichment.")
  }
  if (!exists(".ensure_enrichr_ready", mode = "function")) {
    .log_abort("Missing internal helper `.ensure_enrichr_ready()`. In source()-based runs, source('R/utils_helpers.R') before this file.")
  }
  .ensure_enrichr_ready(site = "Enrichr", verbose = verbose)

  files <- list.files(filtered_dir, "_filtered_links_(up|down)\\.csv$", full.names = TRUE)
  if (!length(files)) {
    .log_warn("No *_filtered_links_(up|down).csv files found in {filtered_dir}.")
    return(invisible(character(0)))
  }

  can_plot <- isTRUE(plot_subnetwork) &&
    requireNamespace("htmlwidgets", quietly = TRUE) &&
    exists("plot_tf_network_delta", mode = "function")
  if (isTRUE(plot_subnetwork) && !isTRUE(can_plot) && isTRUE(verbose)) {
    .log_warn("Pathway subnetwork plotting skipped: need htmlwidgets and plot_tf_network_delta().")
  }
  subnetwork_root <- file.path(filtered_dir, subnetwork_dirname)
  if (isTRUE(can_plot)) dir.create(subnetwork_root, recursive = TRUE, showWarnings = FALSE)
  summary_rows <- list()

  .split_genes <- function(x) {
    if (!length(x) || is.na(x) || !nzchar(x)) return(character(0))
    g <- unlist(strsplit(as.character(x), "[,;|/]+"))
    g <- trimws(g)
    g[!is.na(g) & nzchar(g)]
  }
  .select_nonredundant_pathways <- function(tbl, overlap_thr = 0.8) {
    if (!is.data.frame(tbl) || !nrow(tbl)) return(tbl)
    if (!all(c("genes", "adjusted_p") %in% names(tbl))) return(tbl)
    n <- nrow(tbl)
    gene_sets <- lapply(tbl$genes, .split_genes)
    gene_sets <- lapply(gene_sets, unique)
    sizes <- vapply(gene_sets, length, integer(1))
    keep <- rep(TRUE, n)
    for (i in seq_len(n)) {
      if (!keep[i]) next
      gi <- gene_sets[[i]]
      if (!length(gi)) next
      for (j in seq.int(i + 1L, n)) {
        if (j > n || !keep[j]) next
        gj <- gene_sets[[j]]
        if (!length(gj)) next
        ov <- length(intersect(gi, gj))
        min_sz <- min(length(gi), length(gj))
        if (min_sz <= 0L) next
        frac_small <- ov / min_sz
        containment <- ov == min_sz
        if (!(isTRUE(containment) || frac_small >= as.numeric(overlap_thr))) next

        if (sizes[i] > sizes[j]) {
          keep[j] <- FALSE
        } else if (sizes[j] > sizes[i]) {
          keep[i] <- FALSE
          break
        } else {
          pi <- suppressWarnings(as.numeric(tbl$adjusted_p[i]))
          pj <- suppressWarnings(as.numeric(tbl$adjusted_p[j]))
          if (!is.finite(pi)) pi <- Inf
          if (!is.finite(pj)) pj <- Inf
          if (pi <= pj) {
            keep[j] <- FALSE
          } else {
            keep[i] <- FALSE
            break
          }
        }
      }
    }
    out <- tbl[keep, , drop = FALSE]
    out <- out[order(out$adjusted_p), , drop = FALSE]
    out
  }

  out_files <- character(0)
  for (csv in files) {
    file_base <- sub("\\.csv$", "", basename(csv), ignore.case = TRUE)
    dir_label <- if (grepl("_up$", file_base, ignore.case = TRUE)) "Up" else "Down"
    comparison_id <- sub("_(up|down)$", "", file_base, ignore.case = TRUE)
    comparison_id <- sub("_filtered_links$", "", comparison_id, ignore.case = TRUE)

    out_csv <- sub("\\.csv$", "_pathway_enrichment.csv", csv, ignore.case = TRUE)
    if (!isTRUE(overwrite) && file.exists(out_csv)) {
      out_files <- c(out_files, out_csv)
      res <- readr::read_csv(out_csv, show_col_types = FALSE)
    } else {
      dat <- readr::read_csv(csv, show_col_types = FALSE, col_select = c("gene_key"))
      genes <- unique(as.character(dat$gene_key))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      if (length(genes) < as.integer(min_genes)) {
        if (isTRUE(verbose)) .log_inform("Skip {basename(csv)}: distinct gene_key={length(genes)} < {as.integer(min_genes)}.")
        next
      }

      enr <- tryCatch(enrichR::enrichr(genes, dbs), error = function(e) NULL)
      if (is.null(enr)) {
        .log_warn("Enrichr failed for {basename(csv)}.")
        next
      }

      rows <- lapply(names(enr), function(db) {
        tbl <- enr[[db]]
        if (!is.data.frame(tbl) || !nrow(tbl) || !all(c("Term", "Adjusted.P.value") %in% names(tbl))) return(NULL)
        out <- data.frame(
          db = db,
          pathway = as.character(tbl$Term),
          adjusted_p = suppressWarnings(as.numeric(tbl$Adjusted.P.value)),
          p_value = if ("P.value" %in% names(tbl)) suppressWarnings(as.numeric(tbl$P.value)) else NA_real_,
          combined_score = if ("Combined.Score" %in% names(tbl)) suppressWarnings(as.numeric(tbl$Combined.Score)) else NA_real_,
          overlap = if ("Overlap" %in% names(tbl)) as.character(tbl$Overlap) else NA_character_,
          genes = if ("Genes" %in% names(tbl)) as.character(tbl$Genes) else NA_character_,
          stringsAsFactors = FALSE
        )
        out
      })
      res <- dplyr::bind_rows(rows)
      if (!nrow(res)) next
      res <- res[is.finite(res$adjusted_p) & res$adjusted_p <= as.numeric(padj_cut), , drop = FALSE]
      if (!nrow(res)) next
      res <- res[order(res$adjusted_p), , drop = FALSE]
      readr::write_csv(res, out_csv)
      out_files <- c(out_files, out_csv)
      if (isTRUE(verbose)) .log_inform("Saved pathway enrichment: {out_csv}")
    }

    if (!isTRUE(can_plot) || !nrow(res)) next

    links_full <- readr::read_csv(csv, show_col_types = FALSE)
    if (!all(c("gene_key") %in% names(links_full))) next
    if (!all(c("adjusted_p", "pathway", "genes") %in% names(res))) next
    res <- res[is.finite(res$adjusted_p), , drop = FALSE]
    if (!nrow(res)) next
    res <- res[order(res$adjusted_p), , drop = FALSE]
    res <- .select_nonredundant_pathways(
      res,
      overlap_thr = pathway_gene_overlap_thresh
    )
    if (!nrow(res)) next
    plot_rows <- res[seq_len(min(as.integer(top_n_pathways_plot), nrow(res))), , drop = FALSE]
    out_dir <- file.path(subnetwork_root, .safe_label(comparison_id), dir_label)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    for (i in seq_len(nrow(plot_rows))) {
      p_name <- as.character(plot_rows$pathway[i])
      p_genes <- .split_genes(plot_rows$genes[i])
      if (!length(p_genes)) next
      sub_links <- links_full[as.character(links_full$gene_key) %in% p_genes, , drop = FALSE]
      if (!nrow(sub_links)) next

      p_slug <- .safe_label(p_name)
      if (!nzchar(p_slug)) p_slug <- paste0("pathway_", i)
      if (nchar(p_slug) > 120L) p_slug <- substr(p_slug, 1L, 120L)
      out_html <- file.path(out_dir, paste0(p_slug, ".html"))
      if (!isTRUE(overwrite) && file.exists(out_html)) {
        next
      }

      plot_title <- paste(comparison_id, dir_label, p_name, sep = " | ")
      w <- try(
        plot_tf_network_delta(
          data = sub_links,
          plot_title = plot_title,
          layout_algo = "fr",
          physics = TRUE,
          add_direct = TRUE,
          edge_filter_min = 0,
          min_delta_abs = 0,
          keep_top_edges_per_tf = 6000,
          peak_mode = "show_all",
          show_peaks = FALSE,
          gene_fc_thresh = 1.5,
          de_reference = "str_over_ctrl",
          motif_db = "jaspar2024"
        ),
        silent = TRUE
      )
      if (inherits(w, "try-error")) next
      htmlwidgets::saveWidget(w, out_html, selfcontained = isTRUE(html_selfcontained))
      if (exists(".set_html_title")) .set_html_title(out_html, plot_title)

      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        comparison = comparison_id,
        direction = dir_label,
        pathway_rank = i,
        pathway = p_name,
        adjusted_p = as.numeric(plot_rows$adjusted_p[i]),
        n_links = nrow(sub_links),
        html = out_html,
        stringsAsFactors = FALSE
      )
    }
  }

  if (isTRUE(can_plot) && length(summary_rows)) {
    summary_tbl <- dplyr::bind_rows(summary_rows)
    readr::write_csv(summary_tbl, file.path(subnetwork_root, "pathway_sub_network_summary.csv"))
  }

  invisible(out_files)
}

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
