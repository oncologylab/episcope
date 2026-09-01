# Module 2 TF-regulon browser payload and standalone HTML writer.

.module2_regulon_join_links <- function(module2, tf = NULL, target_gene = NULL) {
  links <- query_module2_links(
    module2,
    tf = tf,
    target_gene = target_gene,
    pass_only = TRUE
  )
  links <- data.table::as.data.table(links)
  if (!nrow(links)) return(data.table::data.table())
  candidates <- data.table::as.data.table(.module2_report_read_table(
    module2,
    "module2_fp_target_candidates",
    columns = c("candidate_id", "fp_id", "target_gene")
  ))
  tf_corr <- data.table::as.data.table(.module2_report_read_table(
    module2,
    "module2_tf_target_corr",
    columns = c("tf", "target_gene", "best_r", "pass")
  ))
  fp_corr <- data.table::as.data.table(.module2_report_read_table(
    module2,
    "module2_fp_target_corr",
    columns = c("fp_id", "target_gene", "best_r", "pass")
  ))
  links[, `:=`(
    tf = as.character(tf),
    target_gene = as.character(target_gene)
  )]
  candidates[, target_gene := as.character(target_gene)]
  tf_corr[, `:=`(
    tf = as.character(tf),
    target_gene = as.character(target_gene)
  )]
  fp_corr[, target_gene := as.character(target_gene)]
  out <- candidates[links, on = "candidate_id", nomatch = 0L]
  out <- tf_corr[pass %in% TRUE][out, on = c("tf", "target_gene"), nomatch = 0L]
  data.table::setnames(out, "best_r", "rna_r", skip_absent = TRUE)
  out <- fp_corr[pass %in% TRUE][out, on = c("fp_id", "target_gene"), nomatch = 0L]
  data.table::setnames(out, "best_r", "fp_r", skip_absent = TRUE)
  if (!"link_id" %in% names(out)) out[, link_id := paste(tf, fp_id, target_gene, sep = "::")]
  out[, `:=`(
    tf = as.character(tf),
    target_gene = as.character(target_gene),
    fp_id = as.character(fp_id),
    fp_r = suppressWarnings(as.numeric(fp_r)),
    rna_r = suppressWarnings(as.numeric(rna_r))
  )]
  unique(out[is.finite(fp_r) & is.finite(rna_r)], by = "link_id")
}

.module2_regulon_condition_payload <- function(links,
                                                nodes,
                                                multiomic_data = NULL,
                                                conditions = NULL) {
  active <- NULL
  if (is.null(multiomic_data)) {
    return(list(condition_links = data.table::data.table(), node_expression = data.table::data.table()))
  }
  if (!is_multiomic_object(multiomic_data)) {
    .log_abort("`multiomic_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  }
  validate_multiomic_object(multiomic_data)
  mats <- multiomic_data$matrices
  available <- Reduce(intersect, list(
    colnames(mats$fp_score), colnames(mats$fp_bound),
    colnames(mats$gene_expr), colnames(mats$gene_on)
  ))
  if (is.null(conditions)) conditions <- available
  conditions <- intersect(as.character(conditions), available)
  if (!length(conditions)) .log_abort("No requested regulon conditions are available in `multiomic_data`.")
  condition_labels <- stats::setNames(conditions, conditions)
  sample_meta <- multiomic_data$samples
  if (is.data.frame(sample_meta) && nrow(sample_meta)) {
    key_cols <- intersect(c("name", "condition_id", "id", "sample_id"), names(sample_meta))
    label_cols <- intersect(c("name", "nutrient_expanded_label", "condition_id"), names(sample_meta))
    label_col <- label_cols[1L]
    if (length(key_cols) && length(label_col)) {
      proposed <- vapply(conditions, function(cc) {
        idx <- NA_integer_
        for (nm in key_cols) {
          hit <- match(cc, as.character(sample_meta[[nm]]))
          if (!is.na(hit)) {
            idx <- hit
            break
          }
        }
        value <- if (is.na(idx)) cc else as.character(sample_meta[[label_col]][idx])
        if (is.na(value) || !nzchar(value)) cc else value
      }, character(1L))
      duplicate_groups <- split(seq_along(proposed), proposed)
      for (idx in duplicate_groups) {
        if (length(idx) > 1L) {
          proposed[idx] <- paste0(proposed[idx], " (replicate ", seq_along(idx), ")")
        }
      }
      condition_labels[] <- proposed
    }
  }

  link_dt <- data.table::as.data.table(links)
  fp_idx <- match(link_dt$fp_id, rownames(mats$fp_score))
  tf_idx <- match(link_dt$from, toupper(rownames(mats$gene_expr)))
  gene_idx <- match(link_dt$to, toupper(rownames(mats$gene_expr)))
  condition_rows <- lapply(conditions, function(cc) {
    score <- as.numeric(mats$fp_score[fp_idx, cc])
    bound <- as.logical(mats$fp_bound[fp_idx, cc])
    tf_on <- as.logical(mats$gene_on[tf_idx, cc])
    target_on <- as.logical(mats$gene_on[gene_idx, cc])
    out <- data.table::data.table(
      condition = unname(condition_labels[[cc]]),
      from = link_dt$from,
      to = link_dt$to,
      fp_id = link_dt$fp_id,
      fp_score = score,
      active = is.finite(score) & bound %in% TRUE & tf_on %in% TRUE & target_on %in% TRUE
    )
    out[active %in% TRUE]
  })
  node_ids <- as.character(nodes$id)
  node_idx <- match(node_ids, toupper(rownames(mats$gene_expr)))
  expression_rows <- lapply(conditions, function(cc) {
    data.table::data.table(
      condition = unname(condition_labels[[cc]]),
      id = node_ids,
      expression = as.numeric(mats$gene_expr[node_idx, cc])
    )
  })
  list(
    condition_links = data.table::rbindlist(condition_rows, use.names = TRUE),
    node_expression = data.table::rbindlist(expression_rows, use.names = TRUE),
    condition_labels = condition_labels
  )
}

.module2_regulon_build <- function(module2,
                                   tfs,
                                   top_n = 100L,
                                   default_top_n = NULL,
                                   target_genes = NULL,
                                   multiomic_data = NULL,
                                   conditions = NULL,
                                   default_condition = NULL,
                                   max_linked_tfs = 50L,
                                   max_linked_tfs_per_target = 2L,
                                   supporting_module2 = NULL,
                                   default_fp_r_cutoff = 0.5,
                                   expression_pseudocount = 1,
                                   title = NULL) {
  tf <- target_gene <- fp_r <- rna_r <- target_rank <- NULL
  tfs <- unique(toupper(as.character(tfs)))
  tfs <- tfs[nzchar(tfs)]
  seed <- .module2_regulon_join_links(module2, tf = tfs)
  if (!is.null(target_genes)) {
    target_genes <- unique(toupper(trimws(as.character(target_genes))))
    target_genes <- target_genes[nzchar(target_genes)]
    seed <- seed[toupper(target_gene) %chin% target_genes]
  }
  if (!nrow(seed)) return(NULL)
  ranked <- seed[, .(
    best_fp_r = max(fp_r, na.rm = TRUE),
    best_rna_r = max(rna_r, na.rm = TRUE)
  ), by = target_gene]
  data.table::setorder(ranked, -best_fp_r, -best_rna_r, target_gene)
  if (is.finite(top_n)) ranked <- head(ranked, as.integer(top_n))
  ranked[, target_rank := seq_len(.N)]
  targets <- ranked$target_gene

  supporting_source <- if (is.null(supporting_module2)) module2 else supporting_module2
  all_links <- .module2_regulon_join_links(supporting_source, target_gene = targets)
  if (!is.null(supporting_module2)) {
    seed_links <- seed[target_gene %chin% targets]
    all_links <- data.table::rbindlist(
      list(seed_links, all_links[!tf %chin% tfs]),
      use.names = TRUE,
      fill = TRUE
    )
  }
  if (!nrow(all_links)) return(NULL)
  pair_rank <- all_links[, .(
    pair_fp_r = max(fp_r, na.rm = TRUE),
    pair_rna_r = max(rna_r, na.rm = TRUE)
  ), by = .(tf, target_gene)]
  pair_rank <- merge(pair_rank, ranked[, .(target_gene, target_rank)], by = "target_gene", sort = FALSE)
  seed_pairs <- pair_rank[tf %chin% tfs]
  supporting_pairs <- pair_rank[!tf %chin% tfs]
  data.table::setorder(supporting_pairs, target_gene, -pair_fp_r, -pair_rna_r, tf)
  keep_supporting <- character()
  if (
    nrow(supporting_pairs) &&
      is.finite(max_linked_tfs) && max_linked_tfs > 0L &&
      is.finite(max_linked_tfs_per_target) && max_linked_tfs_per_target > 0L
  ) {
    supporting_pairs <- supporting_pairs[, head(.SD, as.integer(max_linked_tfs_per_target)), by = target_gene]
    regulator_rank <- supporting_pairs[, .(
      n_targets = data.table::uniqueN(target_gene),
      best_fp_r = max(pair_fp_r, na.rm = TRUE),
      best_rna_r = max(pair_rna_r, na.rm = TRUE)
    ), by = tf]
    data.table::setorder(regulator_rank, -n_targets, -best_fp_r, -best_rna_r, tf)
    keep_supporting <- head(regulator_rank$tf, as.integer(max_linked_tfs))
  } else {
    supporting_pairs <- supporting_pairs[0]
  }
  keep_pairs <- data.table::rbindlist(list(
    seed_pairs[, .(tf, target_gene)],
    supporting_pairs[tf %chin% keep_supporting, .(tf, target_gene)]
  ))
  all_links <- keep_pairs[all_links, on = c("tf", "target_gene"), nomatch = 0L]
  all_links <- merge(all_links, ranked[, .(target_gene, target_rank)], by = "target_gene", sort = FALSE)
  links <- all_links[, .(
    from = tf,
    to = target_gene,
    fp_id,
    fp_r,
    rna_r,
    target_rank,
    is_seed_edge = tf %chin% tfs
  )]
  from_ids <- unique(links$from)
  node_ids <- unique(c(from_ids, links$to))
  nodes <- data.table::data.table(id = node_ids)
  nodes <- merge(nodes, ranked[, .(id = target_gene, target_rank)], by = "id", all.x = TRUE, sort = FALSE)
  nodes[, `:=`(
    node_type = ifelse(id %chin% from_ids, "TF", "Gene"),
    node_role = ifelse(id %chin% tfs, "seed_tf", ifelse(id %chin% targets, "selected_target", "supporting_tf")),
    is_seed_tf = id %chin% tfs,
    is_target_gene = id %chin% targets,
    perturbation_log2fc = NA_real_
  )]
  condition_payload <- .module2_regulon_condition_payload(
    links = links,
    nodes = nodes,
    multiomic_data = multiomic_data,
    conditions = conditions
  )
  if (!is.null(default_condition) && default_condition %in% names(condition_payload$condition_labels)) {
    default_condition <- unname(condition_payload$condition_labels[[default_condition]])
  }
  if (is.null(title)) {
    target_label <- if (!is.null(target_genes)) {
      sprintf("selected-%d", nrow(ranked))
    } else if (is.finite(top_n)) {
      sprintf("top-%d", as.integer(top_n))
    } else {
      "all"
    }
    title <- sprintf("%s Module2 %s target network", paste(tfs, collapse = "/"), target_label)
  }
  .module2_regulon_payload(
    nodes = nodes,
    links = links,
    condition_links = condition_payload$condition_links,
    node_expression = condition_payload$node_expression,
    seed_tfs = tfs,
    default_condition = default_condition,
    default_top_n = default_top_n,
    default_fp_r_cutoff = default_fp_r_cutoff,
    expression_pseudocount = expression_pseudocount,
    title = title
  )
}

.module2_regulon_payload <- function(nodes,
                                     links,
                                     condition_links = NULL,
                                     node_expression = NULL,
                                     seed_tfs = NULL,
                                     default_condition = NULL,
                                     default_top_n = NULL,
                                     default_fp_r_cutoff = 0.5,
                                     expression_pseudocount = 1,
                                     title = "Module 2 TF-target network") {
  id <- from <- to <- fp_id <- fp_r <- rna_r <- target_rank <- NULL
  condition <- fp_score <- active <- expression <- NULL
  node_type <- node_role <- is_seed_tf <- is_target_gene <- NULL
  perturbation_log2fc <- NULL

  node_dt <- data.table::as.data.table(nodes)
  link_dt <- data.table::as.data.table(links)
  if (!nrow(node_dt) || !nrow(link_dt)) {
    return(list(nodes = data.table::data.table(), links = data.table::data.table()))
  }
  if (!"id" %in% names(node_dt) && "node_id" %in% names(node_dt)) {
    data.table::setnames(node_dt, "node_id", "id")
  }
  if (!"fp_id" %in% names(link_dt)) {
    link_dt[, fp_id := if ("best_peak_ID" %in% names(link_dt)) as.character(best_peak_ID) else paste0("edge_", .I)]
  }
  if (!"fp_r" %in% names(link_dt)) {
    link_dt[, fp_r := if ("edge_r" %in% names(link_dt)) as.numeric(edge_r) else if ("abs_edge_score" %in% names(link_dt)) as.numeric(abs_edge_score) else 0]
  }
  if (!"rna_r" %in% names(link_dt)) link_dt[, rna_r := 1]
  if (!"target_rank" %in% names(link_dt)) link_dt[, target_rank := NA_integer_]
  if (!"is_seed_edge" %in% names(link_dt)) link_dt[, is_seed_edge := FALSE]
  link_dt[, `:=`(
    from = toupper(as.character(from)),
    to = toupper(as.character(to)),
    fp_id = as.character(fp_id),
    fp_r = suppressWarnings(as.numeric(fp_r)),
    rna_r = suppressWarnings(as.numeric(rna_r)),
    target_rank = suppressWarnings(as.integer(target_rank)),
    is_seed_edge = as.logical(is_seed_edge)
  )]
  link_dt <- unique(link_dt[, .(from, to, fp_id, fp_r, rna_r, target_rank, is_seed_edge)])
  link_dt <- link_dt[nzchar(from) & nzchar(to) & is.finite(fp_r) & is.finite(rna_r)]

  node_defaults <- list(
    node_type = "Gene",
    node_role = "selected_target",
    is_seed_tf = FALSE,
    is_target_gene = FALSE,
    target_rank = NA_integer_,
    perturbation_log2fc = NA_real_
  )
  for (nm in names(node_defaults)) {
    if (!nm %in% names(node_dt)) node_dt[, (nm) := node_defaults[[nm]]]
  }
  seed_tfs <- unique(toupper(as.character(seed_tfs)))
  seed_tfs <- seed_tfs[nzchar(seed_tfs)]
  node_dt[, `:=`(
    id = toupper(as.character(id)),
    node_type = as.character(node_type),
    node_role = as.character(node_role),
    is_seed_tf = as.logical(is_seed_tf) | id %chin% seed_tfs,
    is_target_gene = as.logical(is_target_gene),
    target_rank = suppressWarnings(as.integer(target_rank)),
    perturbation_log2fc = suppressWarnings(as.numeric(perturbation_log2fc))
  )]
  node_dt <- unique(node_dt[, .(id, node_type, node_role, is_seed_tf, is_target_gene, target_rank, perturbation_log2fc)], by = "id")

  condition_dt <- data.table::as.data.table(condition_links)
  if (nrow(condition_dt)) {
    required <- c("condition", "from", "to", "fp_id", "fp_score", "active")
    missing <- setdiff(required, names(condition_dt))
    if (length(missing)) .log_abort("Condition regulon payload is missing: {paste(missing, collapse = ', ')}")
    condition_dt[, `:=`(
      condition = as.character(condition),
      from = toupper(as.character(from)),
      to = toupper(as.character(to)),
      fp_id = as.character(fp_id),
      fp_score = suppressWarnings(as.numeric(fp_score)),
      active = as.logical(active)
    )]
    condition_dt <- unique(condition_dt[, .(condition, from, to, fp_id, fp_score, active)])
    condition_dt <- condition_dt[nzchar(condition)]
  } else {
    condition_dt <- data.table::data.table(
      condition = character(), from = character(), to = character(),
      fp_id = character(), fp_score = numeric(), active = logical()
    )
  }

  expression_dt <- data.table::as.data.table(node_expression)
  if (nrow(expression_dt)) {
    if (!"id" %in% names(expression_dt) && "node_id" %in% names(expression_dt)) {
      data.table::setnames(expression_dt, "node_id", "id")
    }
    expression_dt[, `:=`(
      condition = as.character(condition),
      id = toupper(as.character(id)),
      expression = suppressWarnings(as.numeric(expression))
    )]
    expression_dt <- unique(expression_dt[, .(condition, id, expression)])
  } else {
    expression_dt <- data.table::data.table(condition = character(), id = character(), expression = numeric())
  }
  available_conditions <- sort(unique(c(condition_dt$condition, expression_dt$condition)))
  if (!is.null(default_condition) && !default_condition %in% available_conditions) {
    replicate_match <- available_conditions[startsWith(
      available_conditions,
      paste0(default_condition, " (replicate ")
    )]
    if (length(replicate_match)) default_condition <- replicate_match[[1L]]
  }
  if (is.null(default_condition) || !default_condition %in% available_conditions) {
    default_condition <- if (length(available_conditions)) available_conditions[[1L]] else ""
  }
  max_target_rank <- suppressWarnings(max(link_dt$target_rank, na.rm = TRUE))
  if (!is.finite(max_target_rank) || max_target_rank < 1) max_target_rank <- 1L
  if (is.null(default_top_n)) {
    default_top_n <- max_target_rank
  } else {
    default_top_n <- as.numeric(default_top_n)[[1L]]
    if (!is.finite(default_top_n)) default_top_n <- max_target_rank
    default_top_n <- max(1L, min(as.integer(default_top_n), as.integer(max_target_rank)))
  }
  minimum_fp_r <- min(link_dt$fp_r, na.rm = TRUE)
  default_fp_r_cutoff <- as.numeric(default_fp_r_cutoff)[[1L]]
  if (!is.finite(default_fp_r_cutoff)) default_fp_r_cutoff <- 0.5
  default_fp_r_cutoff <- max(-1, min(1, default_fp_r_cutoff))
  list(
    title = as.character(title),
    seed_tfs = seed_tfs,
    default_condition = as.character(default_condition),
    default_top_n = as.integer(default_top_n),
    expression_pseudocount = as.numeric(expression_pseudocount)[[1L]],
    minimum_fp_r = minimum_fp_r,
    default_fp_r_cutoff = default_fp_r_cutoff,
    nodes = node_dt,
    links = link_dt,
    condition_links = condition_dt,
    node_expression = expression_dt
  )
}

.module2_regulon_write_html <- function(payload, out_html, max_edges_default = 0L) {
  if (!nrow(payload$nodes) || !nrow(payload$links)) {
    dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
    writeLines("<!doctype html><html><head><meta charset=\"utf-8\"></head><body><b>No edges to plot.</b></body></html>", out_html, useBytes = TRUE)
    return(out_html)
  }
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  json <- function(x) jsonlite::toJSON(x, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null", digits = 7)
  page_title <- .module2_report_browser_html_escape(payload$title)
  packed_payload <- list(
    links = .module2_report_browser_browser_payload_to_columnar(payload$links),
    condition_links = .module2_report_browser_browser_payload_to_columnar(payload$condition_links),
    node_expression = .module2_report_browser_browser_payload_to_columnar(payload$node_expression)
  )
  packed_base64 <- .module2_report_browser_encode_browser_json_deflate_base64(
    packed_payload,
    digits = 15L
  )
  default_fp_r_cutoff <- payload$default_fp_r_cutoff
  if (!length(default_fp_r_cutoff) || !is.finite(default_fp_r_cutoff)) {
    default_fp_r_cutoff <- min(suppressWarnings(as.numeric(payload$links$fp_r)), na.rm = TRUE)
  }
  if (!is.finite(default_fp_r_cutoff)) default_fp_r_cutoff <- -1
  minimum_fp_r <- payload$minimum_fp_r
  if (!length(minimum_fp_r) || !is.finite(minimum_fp_r)) {
    minimum_fp_r <- min(suppressWarnings(as.numeric(payload$links$fp_r)), na.rm = TRUE)
  }
  if (!is.finite(minimum_fp_r)) minimum_fp_r <- -1
  expression_pseudocount <- payload$expression_pseudocount
  if (!length(expression_pseudocount) || !is.finite(expression_pseudocount)) expression_pseudocount <- 1
  default_top_n <- payload$default_top_n
  if (!length(default_top_n) || !is.finite(default_top_n) || default_top_n < 1) default_top_n <- 1L
  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    sprintf("<title>%s</title>", page_title),
    "<style>",
    "body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif}",
    ".wrap{max-width:min(calc((100vh - 245px)*1.7778),calc(100vw - 24px));margin:0 auto;padding:12px 12px 16px}",
    ".top{border-bottom:1px solid #d6d6d0;padding-bottom:9px;margin-bottom:8px}h1{font-size:21px;line-height:1.18;margin:0 0 5px;font-weight:700}",
    ".meta,.note{font-size:11px;line-height:1.35;color:#555;font-weight:700}.meta{margin-bottom:8px}",
    ".controls{display:flex;gap:7px 10px;align-items:center;flex-wrap:wrap}.control{display:flex;gap:5px;align-items:center;font-size:11px;font-weight:700;white-space:nowrap}",
    "select,input,button{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #999;border-radius:3px;padding:5px 7px;background:#fff;color:#111}",
    "input[type=range]{width:110px;padding:0}input[type=number]{width:58px}#fpValue{width:82px}input[type=color]{width:30px;height:27px;padding:2px}select{max-width:210px}",
    "button{cursor:pointer;background:#222;color:#fff;border-color:#222}.danger{background:#9f1239;border-color:#9f1239}.secondary{background:#fff;color:#222;border-color:#777}",
    ".nodePicker{position:relative}.nodePicker>summary{list-style:none;cursor:pointer;font:700 12px Arial;border:1px solid #999;border-radius:3px;padding:6px 8px;background:#fff;min-width:180px}.nodePicker>summary::-webkit-details-marker{display:none}.nodePickerPanel{position:absolute;z-index:20;top:calc(100% + 3px);left:0;width:360px;padding:8px;border:1px solid #888;border-radius:4px;background:#fff;box-shadow:0 5px 18px rgba(0,0,0,.18)}.nodePickerPanel input[type=search],.nodePickerPanel textarea{box-sizing:border-box;width:100%;margin-bottom:6px}.nodePickerPanel textarea{min-height:68px;resize:vertical;font:12px Arial}.nodePickerActions{display:flex;gap:6px;margin-bottom:6px}.nodeChecklist{max-height:300px;overflow:auto;border:1px solid #ddd;padding:4px;background:#fafafa}.nodeCheck{display:flex;align-items:center;gap:6px;padding:3px 4px;font-size:11px;font-weight:700}.nodeCheck:hover{background:#eef2f7}.nodePasteStatus{font-size:10px;color:#555;margin:3px 0 6px}",
    ".appearancePanel{width:100%;display:flex;gap:7px 11px;align-items:center;flex-wrap:wrap;padding:4px 0 0;border:0;background:transparent}.appearancePanel .control{font-size:10px}.appearancePanel input[type=number]{width:52px}",
    ".canvas{position:relative;width:100%;aspect-ratio:16/9;max-height:calc(100vh - 245px);border:1px solid #d6d6d0;background:#fff;overflow:hidden;cursor:grab;box-shadow:0 1px 2px rgba(0,0,0,.04)}",
    ".canvas.panning{cursor:grabbing}svg{width:100%;height:100%;display:block;background:#fff}.tooltip{position:absolute;display:none;background:rgba(17,17,17,.92);color:#fff;font:700 12px Arial;padding:7px 8px;border-radius:3px;pointer-events:none;max-width:380px;line-height:1.35}",
    "</style></head><body><div class=\"wrap\"><div class=\"top\">",
    sprintf("<h1>%s</h1>", page_title),
    "<div class=\"meta\">Passed Module 2 links only. Min FP R and Top targets qualify targets from seed-TF links; supporting TFs never introduce targets. All predictions is condition-independent. A condition edge score is the sum of FP score x FP R x RNA R over active supporting footprints. Comparisons use Cond1 minus Cond2.</div>",
    "<div class=\"controls\">",
    "<label class=\"control\">Cond1 <select id=\"cond1Select\"></select><input id=\"cond1Color\" type=\"color\" value=\"#B2182B\" title=\"Cond1 comparison color\"></label><label class=\"control\">Cond2 <select id=\"cond2Select\"></select><input id=\"cond2Color\" type=\"color\" value=\"#2166AC\" title=\"Cond2 comparison color\"></label>",
    "<label class=\"control\">Min FP R &gt; <input id=\"fpRange\" type=\"range\" min=\"-1\" max=\"1\" step=\"0.01\" value=\"0.5\"><input id=\"fpValue\" type=\"number\" min=\"-1\" max=\"1\" step=\"0.01\" value=\"0.5\"></label>",
    "<label class=\"control\">Top targets after FP R <input id=\"topRange\" type=\"range\" min=\"0\" step=\"1\"><input id=\"topValue\" type=\"number\" min=\"0\" step=\"1\"></label>",
    "<label class=\"control\">Max edges <input id=\"maxEdges\" type=\"number\" min=\"0\"></label>",
    "<label class=\"control\">Layout <select id=\"layout\"><option value=\"force\">Force</option><option value=\"cose\">Compound spring</option><option value=\"breadthfirst\">Breadth-first</option><option value=\"radial\">Radial</option><option value=\"columns\">Columns</option><option value=\"hierarchy\">Hierarchy</option><option value=\"concentric\">Concentric</option><option value=\"circle\">Circle</option><option value=\"grid\">Grid</option><option value=\"spiral\">Spiral</option><option value=\"clustered\">Clustered</option></select></label>",
    "<label class=\"control\">Spacing <input id=\"spacing\" type=\"range\" min=\"0.1\" max=\"2\" step=\"0.01\" value=\"0.5\"><input id=\"spacingValue\" type=\"number\" min=\"0.1\" max=\"2\" step=\"0.01\" value=\"0.5\"></label>",
    "<label class=\"control\">Labels <input id=\"labels\" type=\"checkbox\" checked></label><label class=\"control\">Arrows <input id=\"arrows\" type=\"checkbox\" checked></label>",
    "<details id=\"nodePicker\" class=\"nodePicker\"><summary id=\"nodePickerSummary\">Select nodes</summary><div class=\"nodePickerPanel\"><input id=\"nodeSearch\" type=\"search\" placeholder=\"Search node names...\" aria-label=\"Search node names\"><div class=\"nodePickerActions\"><button id=\"nodeShowAll\" type=\"button\" class=\"secondary\">Check all</button><button id=\"nodeHideAll\" type=\"button\" class=\"secondary\">Uncheck all</button><button id=\"undo\" type=\"button\" class=\"secondary\">Undo</button></div><textarea id=\"nodePaste\" placeholder=\"Paste node names separated by commas, semicolons, spaces, tabs, or new lines\"></textarea><button id=\"nodeApply\" type=\"button\">Show pasted nodes only</button><div id=\"nodePasteStatus\" class=\"nodePasteStatus\"></div><div id=\"nodeChecklist\" class=\"nodeChecklist\"></div></div></details>",
    "<div class=\"appearancePanel\" id=\"appearance\" aria-label=\"Appearance controls\">",
    "<label class=\"control\">Seed TF <input id=\"seedColor\" type=\"color\" value=\"#C2410C\"></label><label class=\"control\">Other TF <input id=\"tfColor\" type=\"color\" value=\"#2166AC\"></label><label class=\"control\">Target <input id=\"geneColor\" type=\"color\" value=\"#9CA3A0\"></label><label class=\"control\">Neutral <input id=\"neutralColor\" type=\"color\" value=\"#E0E3E1\"></label><label class=\"control\">Single edge <input id=\"edgeColor\" type=\"color\" value=\"#6B7280\"></label>",
    "<label class=\"control\">Target dot min/max <input id=\"geneMin\" type=\"number\" min=\"2\" max=\"30\" step=\"0.5\" value=\"4\"><input id=\"geneMax\" type=\"number\" min=\"2\" max=\"30\" step=\"0.5\" value=\"10\"></label><label class=\"control\">TF box min/max <input id=\"tfMin\" type=\"number\" min=\"6\" max=\"40\" step=\"0.5\" value=\"10\"><input id=\"tfMax\" type=\"number\" min=\"6\" max=\"40\" step=\"0.5\" value=\"18\"></label>",
    "<label class=\"control\">Edge min/max <input id=\"edgeMin\" type=\"number\" min=\"0.05\" max=\"12\" step=\"0.05\" value=\"0.25\"><input id=\"edgeMax\" type=\"number\" min=\"0.05\" max=\"12\" step=\"0.05\" value=\"2\"></label><label class=\"control\">Edge opacity <input id=\"edgeOpacity\" type=\"range\" min=\"0.05\" max=\"1\" step=\"0.01\" value=\"0.58\"><input id=\"edgeOpacityValue\" type=\"number\" min=\"0.05\" max=\"1\" step=\"0.01\" value=\"0.58\"></label><label class=\"control\">Node opacity <input id=\"nodeOpacity\" type=\"range\" min=\"0.1\" max=\"1\" step=\"0.01\" value=\"0.95\"><input id=\"nodeOpacityValue\" type=\"number\" min=\"0.1\" max=\"1\" step=\"0.01\" value=\"0.95\"></label>",
    "<label class=\"control\">Arrow size <input id=\"arrowSize\" type=\"number\" min=\"3\" max=\"20\" step=\"0.5\" value=\"8\"></label><label class=\"control\">TF label <input id=\"tfLabelSize\" type=\"number\" min=\"6\" max=\"24\" step=\"0.5\" value=\"10\"></label><label class=\"control\">Target label <input id=\"geneLabelSize\" type=\"number\" min=\"6\" max=\"24\" step=\"0.5\" value=\"9\"></label><button id=\"resetAppearance\" class=\"secondary\" type=\"button\">Reset appearance</button></div>",
    "<button id=\"reset\">Reset</button><button id=\"exportSvg\">Export SVG</button>",
    "</div></div><p class=\"note\" id=\"stats\"></p>",
    "<div class=\"canvas\" id=\"canvas\"><div class=\"tooltip\" id=\"tooltip\"></div><svg id=\"networkSvg\" viewBox=\"0 0 1600 900\" role=\"img\" aria-label=\"TF target network\"><g id=\"viewLayer\"><g id=\"edgeLayer\"></g><g id=\"arrowLayer\"></g><g id=\"nodeLayer\"></g><g id=\"labelLayer\"></g></g></svg></div></div>",
    "<script>",
    sprintf("const FULL_NODES=%s;", json(payload$nodes)),
    sprintf("const REGULON_PAYLOAD_DEFLATE_BASE64=%s;", json(packed_base64)),
    "let FULL_LINKS=[],FULL_CONDITION_LINKS=[],FULL_NODE_EXPRESSION=[];",
    sprintf("const SEED_TFS=%s;", json(as.list(payload$seed_tfs))),
    sprintf("const DEFAULT_CONDITION=%s;", json(payload$default_condition)),
    sprintf("const EXPRESSION_PSEUDOCOUNT=%s;", format(expression_pseudocount, scientific = FALSE, trim = TRUE)),
    sprintf("const DEFAULT_FP_R_CUTOFF=%s;", format(default_fp_r_cutoff, scientific = FALSE, trim = TRUE)),
    sprintf("const MINIMUM_FP_R=%s;", format(minimum_fp_r, scientific = FALSE, trim = TRUE)),
    sprintf("const DEFAULT_TOP_N=%d;", as.integer(default_top_n)),
    sprintf("const DEFAULT_MAX_EDGES=%d;", as.integer(max_edges_default)),
    "const WIDTH=1600,HEIGHT=900,NS='http://www.w3.org/2000/svg';",
    "const $=id=>document.getElementById(id),svg=$('networkSvg'),viewLayer=$('viewLayer'),edgeLayer=$('edgeLayer'),arrowLayer=$('arrowLayer'),nodeLayer=$('nodeLayer'),labelLayer=$('labelLayer'),canvas=$('canvas'),tooltip=$('tooltip');",
    .module2_report_browser_illustrator_svg_js(),
    "const ctl={c1:$('cond1Select'),c2:$('cond2Select'),c1color:$('cond1Color'),c2color:$('cond2Color'),fp:$('fpRange'),fpv:$('fpValue'),top:$('topRange'),topv:$('topValue'),max:$('maxEdges'),layout:$('layout'),space:$('spacing'),spacev:$('spacingValue'),labels:$('labels'),arrows:$('arrows'),stats:$('stats'),picker:$('nodePicker'),pickerSummary:$('nodePickerSummary'),nodeSearch:$('nodeSearch'),nodePaste:$('nodePaste'),nodePasteStatus:$('nodePasteStatus'),nodeChecklist:$('nodeChecklist'),seedColor:$('seedColor'),tfColor:$('tfColor'),geneColor:$('geneColor'),neutralColor:$('neutralColor'),edgeColor:$('edgeColor'),geneMin:$('geneMin'),geneMax:$('geneMax'),tfMin:$('tfMin'),tfMax:$('tfMax'),edgeMin:$('edgeMin'),edgeMax:$('edgeMax'),edgeOpacity:$('edgeOpacity'),edgeOpacityv:$('edgeOpacityValue'),nodeOpacity:$('nodeOpacity'),nodeOpacityv:$('nodeOpacityValue'),arrowSize:$('arrowSize'),tfLabelSize:$('tfLabelSize'),geneLabelSize:$('geneLabelSize')};",
    "let state={nodes:[],edges:[],pickerIds:[],deleted:new Set(),undo:[],fixed:new Map(),selected:'',drag:null,pan:null,view:{x:0,y:0,k:1},needsFit:true,layoutFallback:'',requestedTop:null,conditionRequestedTop:null,predictionAutoAll:false,lastCond1:null};",
    "const CONDITION_CACHE=new Map(),EXPRESSION_CACHE=new Map(),TARGET_CACHE=new Map();",
    "const esc=s=>String(s).replace(/[&<>\"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;'}[c]));",
    "const num=(x,d=0)=>Number.isFinite(Number(x))?Number(x):d;const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));",
    "async function decodePayload(b64){if(!('DecompressionStream'in window))throw new Error('A current browser with DecompressionStream support is required.');const bin=atob(b64),bytes=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)bytes[i]=bin.charCodeAt(i);const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));return JSON.parse(await new Response(stream).text())}",
    "function rowsFromColumns(x){if(Array.isArray(x))return x;const cols=x&&x.columns?x.columns:[],data=x&&x.data?x.data:[],n=data.length?data[0].length:0,rows=new Array(n);for(let i=0;i<n;i++){const row={};for(let j=0;j<cols.length;j++)row[cols[j]]=data[j][i];rows[i]=row}return rows}",
    "function conditions(){return [...new Set(FULL_CONDITION_LINKS.map(x=>x.condition).concat(FULL_NODE_EXPRESSION.map(x=>x.condition)))].sort((a,b)=>a.localeCompare(b));}",
    "function linkKey(x){return x.from+'\\u001f'+x.to+'\\u001f'+x.fp_id}function condMap(c){if(CONDITION_CACHE.has(c))return CONDITION_CACHE.get(c);const m=new Map();FULL_CONDITION_LINKS.filter(x=>x.condition===c).forEach(x=>m.set(linkKey(x),x));CONDITION_CACHE.set(c,m);return m}",
    "function exprMap(c){if(EXPRESSION_CACHE.has(c))return EXPRESSION_CACHE.get(c);const m=new Map();FULL_NODE_EXPRESSION.filter(x=>x.condition===c).forEach(x=>m.set(x.id,num(x.expression,NaN)));EXPRESSION_CACHE.set(c,m);return m}",
    "function targetLimit(){return Math.max(0,Math.floor(num(ctl.topv.value,0)))}function fpCut(){return num(ctl.fpv.value,.3)}function displayMode(){return ctl.c1.value?(ctl.c2.value?'comparison':'single'):'prediction'}",
    "function eligibleTargets(cut){const key=Number(cut).toFixed(12);if(TARGET_CACHE.has(key))return TARGET_CACHE.get(key);const ranks=new Map();FULL_LINKS.forEach(l=>{if(!l.is_seed_edge||num(l.fp_r,-Infinity)<=cut)return;const r=num(l.target_rank,Infinity),old=ranks.get(l.to);if(old===undefined||r<old)ranks.set(l.to,r)});const out=[...ranks].sort((a,b)=>a[1]-b[1]||a[0].localeCompare(b[0])).map(x=>x[0]);TARGET_CACHE.set(key,out);return out}function syncTopBounds(){const eligible=eligibleTargets(fpCut()).length,autoAll=displayMode()==='prediction'&&state.predictionAutoAll,requested=autoAll?eligible:Math.max(0,Math.floor(num(state.requestedTop,targetLimit())));ctl.top.max=eligible;ctl.topv.max=eligible;const next=Math.min(requested,eligible);ctl.top.value=next;ctl.topv.value=next;return eligible}",
    "function aggregate(){const mode=displayMode(),c1=condMap(ctl.c1.value),c2=condMap(ctl.c2.value),cut=fpCut(),eligible=eligibleTargets(cut),allowed=new Set(eligible.slice(0,targetLimit())),groups=new Map();FULL_LINKS.forEach(l=>{if(!allowed.has(l.to)||num(l.fp_r,-Infinity)<=cut)return;const key=l.from+'\\u001f'+l.to;let g=groups.get(key);if(!g){g={from:l.from,to:l.to,target_rank:l.target_rank,is_seed_edge:!!l.is_seed_edge,fp_r:-Infinity,rna_r:-Infinity,score1:0,score2:0,n:0};groups.set(key,g)}g.fp_r=Math.max(g.fp_r,num(l.fp_r));g.rna_r=Math.max(g.rna_r,num(l.rna_r));g.n++;const a=c1.get(linkKey(l)),b=c2.get(linkKey(l));if(a&&a.active&&Number.isFinite(Number(a.fp_score)))g.score1+=num(a.fp_score)*num(l.fp_r)*num(l.rna_r);if(b&&b.active&&Number.isFinite(Number(b.fp_score)))g.score2+=num(b.fp_score)*num(l.fp_r)*num(l.rna_r)});let edges=[...groups.values()].map(e=>{e.predictionScore=Math.abs(e.fp_r*e.rna_r);e.delta=e.score1-e.score2;e.score=mode==='prediction'?e.predictionScore:(mode==='single'?e.score1:Math.abs(e.delta));return e}).filter(e=>mode==='prediction'||e.score>0);edges.sort((x,y)=>y.score-x.score||x.from.localeCompare(y.from)||x.to.localeCompare(y.to));const mx=Math.floor(num(ctl.max.value,0));if(mx>0)edges=edges.slice(0,mx);const availableActive=new Set(SEED_TFS);edges.forEach(e=>{availableActive.add(e.from);availableActive.add(e.to)});const availableNodes=FULL_NODES.filter(n=>availableActive.has(n.id));edges=edges.filter(e=>!state.deleted.has(e.from)&&!state.deleted.has(e.to));const linkedRegulators=new Set(edges.map(e=>e.from));let nodes=availableNodes.filter(n=>!state.deleted.has(n.id)&&(n.node_type!=='TF'||n.is_seed_tf||linkedRegulators.has(n.id))).map(n=>Object.assign({},n));return{nodes,edges,availableNodes,mode,eligibleCount:eligible.length,expr1:exprMap(ctl.c1.value),expr2:exprMap(ctl.c2.value)};}",
    "function nodeMetrics(n,g){const a=g.expr1.get(n.id),b=g.expr2.get(n.id),p=Math.max(0,num(EXPRESSION_PSEUDOCOUNT,1));n.expr1=a;n.expr2=b;n.maxExpr=g.mode==='comparison'?Math.max(Number.isFinite(a)?a:0,Number.isFinite(b)?b:0):(Number.isFinite(a)?a:0);n.log2fc=Number.isFinite(a)&&Number.isFinite(b)?Math.log2((a+p)/(b+p)):NaN;}",
    "function seeded(id){let h=2166136261;for(let i=0;i<id.length;i++){h^=id.charCodeAt(i);h=Math.imul(h,16777619)}return(h>>>0)/4294967295}",
    "function ring(items,r,a0=-Math.PI/2){items.forEach((n,i)=>{const a=a0+2*Math.PI*i/Math.max(1,items.length);n.x=800+Math.cos(a)*r;n.y=450+Math.sin(a)*r})}",
    "function columns(nodes){const tf=nodes.filter(n=>n.node_type==='TF'),genes=nodes.filter(n=>n.node_type!=='TF');[tf,genes].forEach((v,k)=>v.sort((a,b)=>num(a.target_rank,1e9)-num(b.target_rank,1e9)||a.id.localeCompare(b.id)).forEach((n,i)=>{n.x=800+(k?330:-330);n.y=70+i*760/Math.max(1,v.length-1)}))}",
    "function breadth(nodes,edges){const seeds=nodes.filter(n=>n.is_seed_tf),seen=new Set(seeds.map(n=>n.id)),levels=[seeds],adj=new Map();edges.forEach(e=>{if(!adj.has(e.from))adj.set(e.from,[]);adj.get(e.from).push(e.to)});for(let k=0;k<3;k++){const next=[];(levels[k]||[]).forEach(n=>(adj.get(n.id)||[]).forEach(id=>{if(!seen.has(id)){seen.add(id);const x=nodes.find(q=>q.id===id);if(x)next.push(x)}}));if(next.length)levels.push(next)}const rest=nodes.filter(n=>!seen.has(n.id));if(rest.length)levels.push(rest);levels.forEach((v,k)=>v.sort((a,b)=>a.id.localeCompare(b.id)).forEach((n,i)=>{n.x=180+k*(1240/Math.max(1,levels.length-1));n.y=80+i*740/Math.max(1,v.length-1)}))}",
    "function baseLayout(nodes,edges,name){if(name==='columns'||name==='hierarchy'){columns(nodes);return}if(name==='breadthfirst'){breadth(nodes,edges);return}if(name==='grid'){const q=Math.ceil(Math.sqrt(nodes.length));nodes.forEach((n,i)=>{n.x=180+(i%q)*1240/Math.max(1,q-1);n.y=100+Math.floor(i/q)*700/Math.max(1,q-1)});return}if(name==='circle'){ring(nodes,340);return}const seeds=nodes.filter(n=>n.is_seed_tf),tfs=nodes.filter(n=>n.node_type==='TF'&&!n.is_seed_tf),genes=nodes.filter(n=>n.node_type!=='TF');seeds.forEach((n,i)=>{n.x=800+(i-(seeds.length-1)/2)*90;n.y=450});ring(tfs,name==='concentric'?180:230);ring(genes,name==='spiral'?330:360,name==='clustered'?0:-Math.PI/2);if(name==='spiral')genes.forEach((n,i)=>{const a=i*2.4;n.x=800+Math.cos(a)*(90+i*7);n.y=450+Math.sin(a)*(90+i*7)});}",
    "function spring(nodes,edges,collision){baseLayout(nodes,edges,'radial');const by=new Map(nodes.map(n=>[n.id,n]));for(let it=0;it<(collision?150:90);it++){const f=new Map(nodes.map(n=>[n.id,{x:0,y:0}]));for(let i=0;i<nodes.length;i++)for(let j=i+1;j<nodes.length;j++){const a=nodes[i],b=nodes[j],dx=a.x-b.x,dy=a.y-b.y,d2=Math.max(100,dx*dx+dy*dy),d=Math.sqrt(d2),rep=(collision?15000:8500)/d2,ux=dx/d,uy=dy/d;f.get(a.id).x+=ux*rep;f.get(a.id).y+=uy*rep;f.get(b.id).x-=ux*rep;f.get(b.id).y-=uy*rep}edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;const dx=b.x-a.x,dy=b.y-a.y,d=Math.max(1,Math.hypot(dx,dy)),pull=(d-(collision?185:230))*.005;f.get(a.id).x+=dx/d*pull;f.get(a.id).y+=dy/d*pull;f.get(b.id).x-=dx/d*pull;f.get(b.id).y-=dy/d*pull});nodes.forEach(n=>{if(n.is_seed_tf)return;n.x=clamp(n.x+f.get(n.id).x,45,1555);n.y=clamp(n.y+f.get(n.id).y,45,855)})}}",
    "function layout(nodes,edges){const name=ctl.layout.value;state.layoutFallback='';if((name==='force'||name==='cose')&&nodes.length>500){baseLayout(nodes,edges,'grid');state.layoutFallback=' Large-network grid fallback is active because more than 500 nodes are displayed.'}else if(name==='force'||name==='cose')spring(nodes,edges,name==='cose');else baseLayout(nodes,edges,name);const s=clamp(num(ctl.space.value,.5),.1,2);ctl.space.value=s;ctl.spacev.value=s;nodes.forEach(n=>{n.x=800+(n.x-800)*s;n.y=450+(n.y-450)*s})}",
    "function svgEl(tag,attrs){const x=document.createElementNS(NS,tag);Object.entries(attrs||{}).forEach(([k,v])=>x.setAttribute(k,String(v)));return x}",
    "function hexRgb(x){const s=String(x||'').trim();if(/^#[0-9a-f]{6}$/i.test(s)){const q=parseInt(s.slice(1),16);return[(q>>16)&255,(q>>8)&255,q&255]}const m=s.match(/[\\d.]+/g);return m&&m.length>=3?m.slice(0,3).map(Number):[128,128,128]}function mixColor(a,b,t){const x=hexRgb(a),y=hexRgb(b);return'rgb('+x.map((v,i)=>Math.round(v+(y[i]-v)*clamp(t,0,1))).join(',')+')'}function labelContrast(fill){const c=hexRgb(fill),lum=(.2126*c[0]+.7152*c[1]+.0722*c[2])/255;return lum>.62?{fill:'#111111',stroke:'#FFFFFF'}:{fill:'#FFFFFF',stroke:'#4B5563'}}",
    "function boundedPair(a,b,da,db,lo,hi){let x=clamp(num(a.value,da),lo,hi),y=clamp(num(b.value,db),lo,hi);if(x>y)y=x;a.value=x;b.value=y;return[x,y]}function robustMax(values){const x=values.filter(Number.isFinite).filter(v=>v>=0).sort((a,b)=>a-b);return x.length?Math.max(1e-9,x[Math.min(x.length-1,Math.floor((x.length-1)*.95))]):1}function prepareNodeScale(g){g.tfExprMax=robustMax(g.nodes.filter(n=>n.node_type==='TF').map(n=>n.maxExpr));g.geneExprMax=robustMax(g.nodes.filter(n=>n.node_type!=='TF').map(n=>n.maxExpr))}",
    "function nodeRadius(n,g){const pair=n.node_type==='TF'?boundedPair(ctl.tfMin,ctl.tfMax,10,18,6,40):boundedPair(ctl.geneMin,ctl.geneMax,4,10,2,30),lo=pair[0],hi=pair[1];if(g.mode==='prediction')return n.is_seed_tf?hi:(lo+hi)/2;const limit=n.node_type==='TF'?g.tfExprMax:g.geneExprMax,t=Math.sqrt(clamp(Math.max(0,n.maxExpr)/Math.max(limit,1e-9),0,1));return lo+(hi-lo)*t}",
    "function fcColor(fc){if(!Number.isFinite(fc))return ctl.neutralColor.value;const t=clamp(Math.abs(fc)/3,0,1),end=fc>=0?ctl.c1color.value:ctl.c2color.value;return mixColor(ctl.neutralColor.value,end,t)}",
    "function roleColor(n,g){if(g.mode==='comparison')return fcColor(n.log2fc);if(n.is_seed_tf)return ctl.seedColor.value;if(n.node_type==='TF')return ctl.tfColor.value;return ctl.geneColor.value}",
    "function clear(){[edgeLayer,arrowLayer,nodeLayer,labelLayer].forEach(x=>x.replaceChildren())}",
    "function nodeBoundary(n,r,ux,uy){if(n.node_type!=='TF')return r;const hw=1.6*r,hh=.8*r,corner=Math.min(4,hw,hh),tx=Math.abs(ux)>1e-9?hw/Math.abs(ux):Infinity,ty=Math.abs(uy)>1e-9?hh/Math.abs(uy):Infinity,t=Math.min(tx,ty),px=ux*t,py=uy*t;if(Math.abs(px)<=hw-corner+1e-9||Math.abs(py)<=hh-corner+1e-9)return t;const cx=(ux<0?-1:1)*(hw-corner),cy=(uy<0?-1:1)*(hh-corner),dot=ux*cx+uy*cy,disc=Math.max(0,dot*dot-(cx*cx+cy*cy-corner*corner));return dot+Math.sqrt(disc)}",
    "function edgeEnds(a,b,ra,rb){const dx=b.x-a.x,dy=b.y-a.y,d=Math.max(1,Math.hypot(dx,dy)),ux=dx/d,uy=dy/d,sourceDistance=nodeBoundary(a,ra,ux,uy),targetDistance=Math.min(d,Math.max(0,nodeBoundary(b,rb,-ux,-uy)-1)),startDistance=Math.min(sourceDistance,Math.max(0,d-targetDistance));return{x1:a.x+ux*startDistance,y1:a.y+uy*startDistance,x2:b.x-ux*targetDistance,y2:b.y-uy*targetDistance,ux,uy}}",
    "function drawArrow(p,color,opacity,e){const s=clamp(num(ctl.arrowSize.value,8),3,20),w=s*.5,x=p.x2,y=p.y2,pts=[[x,y],[x-p.ux*s-p.uy*w,y-p.uy*s+p.ux*w],[x-p.ux*s+p.uy*w,y-p.uy*s-p.ux*w]].map(q=>q.join(',')).join(' ');ctl.arrowSize.value=s;arrowLayer.appendChild(svgEl('polygon',{points:pts,fill:color,'fill-opacity':opacity,'data-from':e.from,'data-to':e.to,'data-title':'direction'}))}",
    "function edgeTitle(e,g){let out=e.from+' -> '+e.to+'\\nFP R: '+e.fp_r.toFixed(3)+'\\nRNA R: '+e.rna_r.toFixed(3);if(g.mode==='prediction')return out+'\\nPrediction score: '+e.score.toFixed(3);if(g.mode==='single')return out+'\\n'+ctl.c1.value+' edge score: '+e.score1.toFixed(3);return out+'\\n'+ctl.c1.value+' edge score: '+e.score1.toFixed(3)+'\\n'+ctl.c2.value+' edge score: '+e.score2.toFixed(3)+'\\nDelta ('+ctl.c1.value+' - '+ctl.c2.value+'): '+e.delta.toFixed(3)}",
    "function nodeTitle(n,g){let out=n.id+'\\nRole: '+n.node_role;if(g.mode==='single')return out+'\\nExpression ('+ctl.c1.value+'): '+(Number.isFinite(n.expr1)?n.expr1.toFixed(3):'unavailable');if(g.mode==='comparison')return out+'\\nRNA log2FC ('+ctl.c1.value+' / '+ctl.c2.value+'): '+(Number.isFinite(n.log2fc)?n.log2fc.toFixed(3):'unavailable')+'\\nMax expression: '+(Number.isFinite(n.maxExpr)?n.maxExpr.toFixed(3):'unavailable');return out}",
    "function render(){const g=aggregate();g.nodes.forEach(n=>nodeMetrics(n,g));prepareNodeScale(g);layout(g.nodes,g.edges);g.nodes.forEach(n=>{const p=state.fixed.get(n.id);if(p){n.x=p.x;n.y=p.y}});state.nodes=g.nodes;state.edges=g.edges;state.pickerIds=g.availableNodes.map(n=>n.id);clear();const by=new Map(g.nodes.map(n=>[n.id,n])),edgeRange=boundedPair(ctl.edgeMin,ctl.edgeMax,.25,2,.05,12),edgeLimit=robustMax(g.edges.map(e=>Math.abs(e.score))),op=clamp(num(ctl.edgeOpacity.value,.58),.05,1),nodeOp=clamp(num(ctl.nodeOpacity.value,.95),.1,1),tfLabel=clamp(num(ctl.tfLabelSize.value,10),6,24),geneLabel=clamp(num(ctl.geneLabelSize.value,9),6,24);ctl.edgeOpacity.value=op;ctl.edgeOpacityv.value=op;ctl.nodeOpacity.value=nodeOp;ctl.nodeOpacityv.value=nodeOp;ctl.tfLabelSize.value=tfLabel;ctl.geneLabelSize.value=geneLabel;g.edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;const p=edgeEnds(a,b,nodeRadius(a,g),nodeRadius(b,g)),color=g.mode==='comparison'?(e.delta>=0?ctl.c1color.value:ctl.c2color.value):ctl.edgeColor.value,width=edgeRange[0]+(edgeRange[1]-edgeRange[0])*Math.sqrt(clamp(Math.abs(e.score)/edgeLimit,0,1)),line=svgEl('line',{x1:p.x1,y1:p.y1,x2:p.x2,y2:p.y2,stroke:color,'stroke-width':width,'stroke-opacity':op,fill:'none','data-title':edgeTitle(e,g)});edgeLayer.appendChild(line);if(ctl.arrows.checked)drawArrow(p,color,op,e)});g.nodes.forEach(n=>{const r=nodeRadius(n,g),fill=roleColor(n,g),shape=n.node_type==='TF'?'rect':'circle',el=shape==='rect'?svgEl('rect',{x:n.x-r*1.6,y:n.y-r*.8,width:r*3.2,height:r*1.6,rx:4}):svgEl('circle',{cx:n.x,cy:n.y,r});el.setAttribute('fill',fill);el.setAttribute('fill-opacity',nodeOp);el.setAttribute('stroke',n.id===state.selected?'#D7263D':'#FFFFFF');el.setAttribute('stroke-width',n.id===state.selected?'4':'1.2');el.setAttribute('data-id',n.id);el.setAttribute('data-title',nodeTitle(n,g));nodeLayer.appendChild(el);if(ctl.labels.checked){const contrast=n.node_type==='TF'?labelContrast(fill):{fill:'#111111',stroke:'#FFFFFF'},t=svgEl('text',{x:n.node_type==='TF'?n.x:n.x+r+4,y:n.y,'font-family':'Arial, Helvetica, sans-serif','font-size':n.node_type==='TF'?tfLabel:geneLabel,'font-weight':'700','text-anchor':n.node_type==='TF'?'middle':'start','dominant-baseline':'middle',fill:contrast.fill,stroke:contrast.stroke,'stroke-width':n.node_type==='TF'?1.8:2.6,'paint-order':'stroke'});t.textContent=n.id;labelLayer.appendChild(t)}});refreshSelect();const displayedSeedTargets=new Set(g.edges.filter(e=>e.is_seed_edge).map(e=>e.to)).size,prefix=g.nodes.length+' nodes, '+g.edges.length+' aggregated edges, '+displayedSeedTargets+' displayed seed-TF targets; '+g.eligibleCount+' seed-TF targets pass Min FP R. ';ctl.stats.textContent=prefix+(g.mode==='comparison'?'Cond1 color: stronger in '+ctl.c1.value+'; Cond2 color: stronger in '+ctl.c2.value+'. Edge width: absolute Cond1 - Cond2 edge-score delta. Node size: maximum expression.':g.mode==='single'?'Edge width: '+ctl.c1.value+' edge score. Node size: '+ctl.c1.value+' expression.':'All predictions: edge width reflects the embedded prediction score.')+state.layoutFallback;if(state.needsFit){state.needsFit=false;fitView()}else applyView()}",
    "function updateNodePickerSummary(){ctl.pickerSummary.textContent='Select nodes: '+state.nodes.length+' / '+state.pickerIds.length+' shown'}function snapshotVisibility(){state.undo.push([...state.deleted]);if(state.undo.length>50)state.undo.shift()}function applyShownSet(shown){snapshotVisibility();state.pickerIds.forEach(id=>{if(shown.has(id))state.deleted.delete(id);else state.deleted.add(id)});if(state.deleted.has(state.selected))state.selected='';render()}function renderNodeChecklist(){const current=new Set(state.pickerIds),q=ctl.nodeSearch.value.trim().toUpperCase(),rows=FULL_NODES.filter(n=>current.has(n.id)&&(!q||n.id.toUpperCase().includes(q))).sort((a,b)=>a.id.localeCompare(b.id));ctl.nodeChecklist.replaceChildren();rows.forEach(n=>{const label=document.createElement('label');label.className='nodeCheck';const box=document.createElement('input');box.type='checkbox';box.checked=!state.deleted.has(n.id);box.dataset.nodeId=n.id;box.addEventListener('change',()=>{snapshotVisibility();if(box.checked)state.deleted.delete(n.id);else state.deleted.add(n.id);if(state.deleted.has(state.selected))state.selected='';render()});const text=document.createElement('span');text.textContent=n.id;label.append(box,text);ctl.nodeChecklist.appendChild(label)});ctl.nodePasteStatus.textContent=rows.length+' node name'+(rows.length===1?'':'s')+' match the search'}function parsePastedNodes(raw){const current=new Set(state.pickerIds),known=new Map(FULL_NODES.filter(n=>current.has(n.id)).map(n=>[n.id.toUpperCase(),n.id])),matched=new Set(),unmatched=[];const cleaned=String(raw||'').replace(/\\bc\\s*\\(/gi,' ').replace(/[\\[\\]{}()\"'`]/g,' '),chunks=cleaned.split(/[,;|\\t\\r\\n]+/).map(x=>x.trim()).filter(Boolean);const add=value=>{const hit=known.get(value.trim().toUpperCase());if(hit){matched.add(hit);return true}return false};chunks.forEach(chunk=>{if(add(chunk))return;chunk.split(/\\s+/).filter(Boolean).forEach(token=>{if(!add(token))unmatched.push(token)})});if(!chunks.length&&cleaned.trim())add(cleaned.trim());return{matched,unmatched:[...new Set(unmatched)]}}function refreshSelect(){updateNodePickerSummary();if(ctl.picker.open)renderNodeChecklist()}",
    "function applyView(){viewLayer.setAttribute('transform','translate('+state.view.x+' '+state.view.y+') scale('+state.view.k+')')}",
    "function fitView(){const box=viewLayer.getBBox(),pad=55;if(!(box.width>0&&box.height>0)){applyView();return}const k=clamp(Math.min((WIDTH-2*pad)/box.width,(HEIGHT-2*pad)/box.height),.25,2.5);state.view={x:WIDTH/2-(box.x+box.width/2)*k,y:HEIGHT/2-(box.y+box.height/2)*k,k};applyView()}",
    "function point(ev){const r=svg.getBoundingClientRect();return{x:(ev.clientX-r.left)*WIDTH/r.width,y:(ev.clientY-r.top)*HEIGHT/r.height}}function world(ev){const p=point(ev);return{x:(p.x-state.view.x)/state.view.k,y:(p.y-state.view.y)/state.view.k}}",
    "function redraw(){render()}",
    "function exportTag(){if(!ctl.c1.value)return'prediction';return ctl.c2.value?'cond1_'+ctl.c1.value+'_vs_cond2_'+ctl.c2.value:'cond1_'+ctl.c1.value}",
    "function exportSvg(){const viewBox=craftgrnSvgViewBox(state.view,WIDTH,HEIGHT),clone=svg.cloneNode(true),cv=clone.querySelector('#viewLayer');if(cv)cv.removeAttribute('transform');craftgrnPrepareIllustratorSvg(svg,clone,viewBox,WIDTH,HEIGHT);const text=new XMLSerializer().serializeToString(clone),blob=new Blob([text],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=(document.title+'_'+exportTag()).replace(/[^A-Za-z0-9_.-]+/g,'_')+'.svg';document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(a.href),1000)}",
    "function setPair(range,value,v,integer,fallback,redraw=true){const lo=num(range.min,-Infinity),hi=num(range.max,Infinity),step=num(range.step,0);let x=clamp(num(v,fallback),lo,hi);if(integer)x=Math.round(x);else if(step>0){const base=Number.isFinite(lo)?lo:0;x=base+Math.round((x-base)/step)*step;x=Number(clamp(x,lo,hi).toFixed(12))}range.value=x;value.value=range.value;if(redraw)render();return num(range.value,x)}",
    "function syncPair(range,value,integer,fallback,afterSet=null){const set=v=>{const out=setPair(range,value,v,integer,fallback,false);if(afterSet)afterSet(out);else render();return out};range.addEventListener('input',()=>set(range.value));value.addEventListener('input',()=>{if(value.value!==''&&Number.isFinite(Number(value.value)))set(value.value)});value.addEventListener('change',()=>set(value.value))}",
    "function fillConditions(select,rows,includeNone,current,noneLabel='None'){select.replaceChildren();if(includeNone){const none=document.createElement('option');none.value='';none.textContent=noneLabel;select.appendChild(none)}rows.forEach(c=>{const option=document.createElement('option');option.value=c;option.textContent=c;select.appendChild(option)});select.value=rows.includes(current)||includeNone&&current===''?current:includeNone?'':rows[0]||''}",
    "function syncCond2(){const current=ctl.c2.value,rows=conditions().filter(c=>c!==ctl.c1.value);fillConditions(ctl.c2,rows,true,rows.includes(current)?current:'','None')}",
    "function syncConditionControls(){const enabled=!!ctl.c1.value;ctl.c1.disabled=!conditions().length;ctl.c2.disabled=!enabled||ctl.c2.options.length<=1;ctl.c1color.disabled=!enabled;ctl.c2color.disabled=!ctl.c2.value}",
    "async function init(){ctl.stats.textContent='Loading embedded network data...';const packed=await decodePayload(REGULON_PAYLOAD_DEFLATE_BASE64);FULL_LINKS=rowsFromColumns(packed.links);FULL_CONDITION_LINKS=rowsFromColumns(packed.condition_links);FULL_NODE_EXPRESSION=rowsFromColumns(packed.node_expression);const cs=conditions();if(cs.length){fillConditions(ctl.c1,cs,true,cs.includes(DEFAULT_CONDITION)?DEFAULT_CONDITION:cs[0],'All predictions');syncCond2()}else{fillConditions(ctl.c1,[],true,'','All predictions');fillConditions(ctl.c2,[],true,'','None')}syncConditionControls();const allowedMin=clamp(num(MINIMUM_FP_R,-1),-1,1),initialFp=clamp(num(DEFAULT_FP_R_CUTOFF,.5),allowedMin,1),mx=eligibleTargets(initialFp).length,initialTop=clamp(Math.floor(num(DEFAULT_TOP_N,mx)),0,mx);ctl.fp.min=allowedMin;ctl.fpv.min=allowedMin;ctl.fp.value=initialFp;ctl.fpv.value=initialFp;ctl.top.max=mx;ctl.topv.max=mx;ctl.top.value=initialTop;ctl.topv.value=initialTop;state.requestedTop=initialTop;state.conditionRequestedTop=initialTop;state.lastCond1=ctl.c1.value;ctl.max.value=DEFAULT_MAX_EDGES;syncPair(ctl.fp,ctl.fpv,false,initialFp,()=>{syncTopBounds();render()});syncPair(ctl.top,ctl.topv,true,initialTop,value=>{state.requestedTop=value;if(displayMode()==='prediction')state.predictionAutoAll=false;else state.conditionRequestedTop=value;render()});syncPair(ctl.space,ctl.spacev,false,.5);syncPair(ctl.edgeOpacity,ctl.edgeOpacityv,false,.58);syncPair(ctl.nodeOpacity,ctl.nodeOpacityv,false,.95);render();renderNodeChecklist()}",
    "ctl.c1.addEventListener('change',()=>{const wasPrediction=state.lastCond1==='',isPrediction=ctl.c1.value==='';if(!wasPrediction&&isPrediction){state.conditionRequestedTop=state.requestedTop;state.predictionAutoAll=true;state.requestedTop=eligibleTargets(fpCut()).length}else if(wasPrediction&&!isPrediction){state.predictionAutoAll=false;state.requestedTop=Math.max(0,Math.floor(num(state.conditionRequestedTop,DEFAULT_TOP_N)))}state.lastCond1=ctl.c1.value;syncCond2();syncConditionControls();syncTopBounds();render()});ctl.c2.addEventListener('change',()=>{syncConditionControls();render()});[ctl.labels,ctl.arrows,ctl.max].forEach(x=>{x.addEventListener('input',render);x.addEventListener('change',render)});[ctl.c1color,ctl.c2color,ctl.seedColor,ctl.tfColor,ctl.geneColor,ctl.neutralColor,ctl.edgeColor,ctl.geneMin,ctl.geneMax,ctl.tfMin,ctl.tfMax,ctl.edgeMin,ctl.edgeMax,ctl.arrowSize,ctl.tfLabelSize,ctl.geneLabelSize].forEach(x=>{x.addEventListener('input',render);x.addEventListener('change',render)});ctl.layout.addEventListener('change',()=>{state.fixed.clear();state.needsFit=true;render()});ctl.nodeSearch.addEventListener('input',renderNodeChecklist);ctl.picker.addEventListener('toggle',()=>{if(ctl.picker.open)renderNodeChecklist()});$('nodeShowAll').addEventListener('click',()=>applyShownSet(new Set(state.pickerIds)));$('nodeHideAll').addEventListener('click',()=>applyShownSet(new Set()));$('nodeApply').addEventListener('click',()=>{const parsed=parsePastedNodes(ctl.nodePaste.value);applyShownSet(parsed.matched);ctl.nodePasteStatus.textContent=parsed.matched.size+' matched; '+parsed.unmatched.length+' unmatched'+(parsed.unmatched.length?': '+parsed.unmatched.slice(0,20).join(', '):'')});",
    "const APPEARANCE_DEFAULTS={cond1Color:'#B2182B',cond2Color:'#2166AC',seedColor:'#C2410C',tfColor:'#2166AC',geneColor:'#9CA3A0',neutralColor:'#E0E3E1',edgeColor:'#6B7280',geneMin:4,geneMax:10,tfMin:10,tfMax:18,edgeMin:.25,edgeMax:2,arrowSize:8,tfLabelSize:10,geneLabelSize:9};function resetAppearance(){Object.entries(APPEARANCE_DEFAULTS).forEach(([id,value])=>{$(id).value=value});setPair(ctl.edgeOpacity,ctl.edgeOpacityv,.58,false,.58,false);setPair(ctl.nodeOpacity,ctl.nodeOpacityv,.95,false,.95,false);render()}",
    "$('undo').addEventListener('click',()=>{const previous=state.undo.pop();if(previous){state.deleted=new Set(previous);if(state.deleted.has(state.selected))state.selected='';render();renderNodeChecklist()}});$('resetAppearance').addEventListener('click',resetAppearance);$('reset').addEventListener('click',()=>{state.deleted.clear();state.undo=[];state.fixed.clear();state.selected='';state.view={x:0,y:0,k:1};state.needsFit=true;render();renderNodeChecklist()});$('exportSvg').addEventListener('click',exportSvg);",
    "svg.addEventListener('wheel',ev=>{ev.preventDefault();const p=point(ev),old=state.view.k,next=clamp(old*(ev.deltaY<0?1.15:1/1.15),.25,6);state.view.x=p.x-(p.x-state.view.x)*next/old;state.view.y=p.y-(p.y-state.view.y)*next/old;state.view.k=next;applyView()},{passive:false});",
    "nodeLayer.addEventListener('mousedown',ev=>{const id=ev.target.getAttribute('data-id');if(!id)return;ev.stopPropagation();state.selected=id;state.drag=id;render()});svg.addEventListener('mousedown',ev=>{if(ev.target.getAttribute&&ev.target.getAttribute('data-id'))return;const p=point(ev);state.pan={sx:p.x,sy:p.y,x:state.view.x,y:state.view.y};canvas.classList.add('panning')});svg.addEventListener('mousemove',ev=>{if(state.drag){const p=world(ev),x=clamp(p.x,30,1570),y=clamp(p.y,30,870);state.fixed.set(state.drag,{x,y});render()}else if(state.pan){const p=point(ev);state.view.x=state.pan.x+p.x-state.pan.sx;state.view.y=state.pan.y+p.y-state.pan.sy;applyView()}const title=ev.target.getAttribute?ev.target.getAttribute('data-title'):'';if(title){tooltip.innerHTML=esc(title).replace(/\\n/g,'<br>');tooltip.style.display='block';const b=canvas.getBoundingClientRect();tooltip.style.left=(ev.clientX-b.left+12)+'px';tooltip.style.top=(ev.clientY-b.top+12)+'px'}else tooltip.style.display='none'});",
    "window.addEventListener('mouseup',()=>{state.drag=null;state.pan=null;canvas.classList.remove('panning')});svg.addEventListener('mouseleave',()=>{state.drag=null;state.pan=null;tooltip.style.display='none'});init().catch(error=>{ctl.stats.textContent='Failed to load network: '+error.message;console.error(error)});",
    "</script></body></html>"
  )
  writeLines(html, out_html, useBytes = TRUE)
  out_html
}
