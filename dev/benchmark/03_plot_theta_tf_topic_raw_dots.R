#!/usr/bin/env Rscript

# Plot unaggregated condition::TF theta values by TF and topic.

args <- commandArgs(trailingOnly = TRUE)
default_run_root <- paste0(
  "/data/homes/yl814/episcope_test/",
  "nutrient_stress_strict_JASPAR2024_expanded/",
  "regulatory_topics_hpafii_condition_models_de_gene_filtered/",
  "topic_runs/",
  "run_012_lda_gene_peak_K30_shared_de_union_expr_scaled_tf_peak_scaled"
)
theta_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path(default_run_root, "topic_models", "vae_models", "theta_K30.csv")
}
output_file <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  file.path(
    default_run_root,
    "review",
    "run012_lda_K30_tf_topic_raw_theta_dotplot.pdf"
  )
}
topic_terms_file <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  run_root <- dirname(dirname(dirname(normalizePath(
    theta_file,
    winslash = "/",
    mustWork = TRUE
  ))))
  file.path(run_root, "topic_extraction", "topic_terms.csv")
}

cutoffs <- c(0.3, 0.5)
grid_columns <- 6L

required_packages <- c("data.table", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
if (!file.exists(theta_file)) {
  stop("Missing theta file: ", theta_file, call. = FALSE)
}
if (!file.exists(topic_terms_file)) {
  stop("Missing final topic-term file: ", topic_terms_file, call. = FALSE)
}

theta_table <- data.table::fread(theta_file, showProgress = FALSE)
if (ncol(theta_table) < 2L || !identical(names(theta_table)[[1L]], "doc_id")) {
  stop("theta must have doc_id followed by topic columns.", call. = FALSE)
}
if (anyDuplicated(theta_table$doc_id)) {
  stop("theta contains duplicated document IDs.", call. = FALSE)
}
if (any(!grepl("^.+::[^:]+$", theta_table$doc_id))) {
  stop("Every theta document ID must use condition::TF format.", call. = FALSE)
}

topic_names <- names(theta_table)[-1L]
topic_numbers <- suppressWarnings(as.integer(sub("^Topic", "", topic_names)))
if (anyNA(topic_numbers) || anyDuplicated(topic_numbers)) {
  stop("theta topic columns must have unique Topic<number> names.", call. = FALSE)
}
topic_order <- order(topic_numbers)
topic_names <- topic_names[topic_order]
topic_numbers <- topic_numbers[topic_order]
theta_table <- theta_table[, c("doc_id", topic_names), with = FALSE]

theta_matrix <- as.matrix(theta_table[, ..topic_names])
storage.mode(theta_matrix) <- "double"
if (any(!is.finite(theta_matrix))) {
  stop("theta contains non-finite values.", call. = FALSE)
}
if (any(theta_matrix < 0 | theta_matrix > 1)) {
  stop("theta values must be within [0, 1].", call. = FALSE)
}
row_error <- abs(rowSums(theta_matrix) - 1)
if (max(row_error) > 1e-8) {
  stop(
    "theta rows are not normalized; maximum absolute row-sum error is ",
    format(max(row_error), scientific = TRUE),
    ".",
    call. = FALSE
  )
}

document_table <- data.table::data.table(
  doc_id = as.character(theta_table$doc_id),
  condition_id = sub("::[^:]+$", "", as.character(theta_table$doc_id)),
  tf = sub("^.*::", "", as.character(theta_table$doc_id))
)
tf_levels <- sort(unique(document_table$tf), method = "radix")
condition_levels <- sort(unique(document_table$condition_id), method = "radix")
rows_per_page <- length(tf_levels)
page_height <- max(20, 3.5 + 0.135 * length(tf_levels))
tf_index <- stats::setNames(seq_along(tf_levels), tf_levels)
document_table[, tf_index := unname(tf_index[tf])]
document_table[, page := ((tf_index - 1L) %/% rows_per_page) + 1L]
n_pages <- max(document_table$page)

topic_terms <- data.table::fread(topic_terms_file, showProgress = FALSE)
required_term_columns <- c(
  "term_id", "term_group", "topic_num", "in_topic", "assignment_method"
)
missing_term_columns <- setdiff(required_term_columns, names(topic_terms))
if (length(missing_term_columns)) {
  stop(
    "Final topic-term file is missing: ",
    paste(missing_term_columns, collapse = ", "),
    call. = FALSE
  )
}
in_topic_flag <- toupper(as.character(topic_terms$in_topic)) %in% c("TRUE", "T", "1")
tf_term_assignments <- topic_terms[
  term_group == "GENE" & in_topic_flag,
  .(
    tf_key = toupper(sub("^GENE:", "", as.character(term_id))),
    topic_num = suppressWarnings(as.integer(topic_num)),
    assignment_method = as.character(assignment_method)
  )
]
if (!nrow(tf_term_assignments) ||
    any(!grepl("^gammafit_maxprob", tf_term_assignments$assignment_method))) {
  stop(
    "Final GENE-term rows must use GammaFit+MaxProb assignment.",
    call. = FALSE
  )
}
if (tf_term_assignments[, .N, by = tf_key][N > 1L, .N]) {
  stop(
    "Final GammaFit+MaxProb GENE terms assign at least one TF to multiple topics.",
    call. = FALSE
  )
}
tf_term_assignments <- unique(
  tf_term_assignments[tf_key %in% toupper(tf_levels), .(tf_key, topic_num)]
)
if (any(!tf_term_assignments$topic_num %in% topic_numbers)) {
  stop("Final TF-term assignments contain topics absent from theta.", call. = FALSE)
}

theta_long <- data.table::melt(
  theta_table,
  id.vars = "doc_id",
  measure.vars = topic_names,
  variable.name = "topic",
  value.name = "theta",
  variable.factor = FALSE
)
theta_long <- merge(
  theta_long,
  document_table,
  by = "doc_id",
  all.x = TRUE,
  sort = FALSE
)
theta_long[, topic_num := topic_numbers[match(topic, topic_names)]]
data.table::setorder(theta_long, tf_index, topic_num, condition_id)

make_page_data <- function(cutoff, page_number) {
  tf_start <- (page_number - 1L) * rows_per_page + 1L
  tf_end <- min(length(tf_levels), page_number * rows_per_page)
  page_tfs <- tf_levels[seq.int(tf_start, tf_end)]
  page_rows <- data.table::data.table(
    tf = page_tfs,
    tf_row = rev(seq_along(page_tfs))
  )
  cells <- data.table::CJ(
    topic_num = topic_numbers,
    tf = page_tfs,
    unique = TRUE
  )
  cells[, tf_key := toupper(tf)]
  cells <- merge(cells, page_rows, by = "tf", sort = FALSE)
  term_borders <- merge(
    cells,
    tf_term_assignments,
    by = c("tf_key", "topic_num"),
    all = FALSE,
    sort = FALSE
  )

  tiles <- theta_long[
    tf %in% page_tfs & theta >= cutoff,
    .(condition_id, tf, topic_num, theta)
  ]
  data.table::setorder(tiles, tf, topic_num, -theta, condition_id)
  tiles[, slot := seq_len(.N), by = .(tf, topic_num)]
  tiles[, slot_row := (slot - 1L) %/% grid_columns]
  tiles[, slot_column := (slot - 1L) %% grid_columns]
  tiles[, row_count := max(slot_row) + 1L, by = .(tf, topic_num)]
  tiles[, column_count := .N, by = .(tf, topic_num, slot_row)]
  tiles[, x := topic_num + (slot_column - (column_count - 1) / 2) * 0.16]
  tiles[, y_offset := (slot_row - (row_count - 1) / 2) * 0.285]
  tiles <- merge(tiles, page_rows, by = "tf", sort = FALSE)
  tiles[, y := tf_row + y_offset]

  list(
    cells = cells,
    tiles = tiles,
    term_borders = term_borders,
    page_rows = page_rows,
    tf_start = tf_start,
    tf_end = tf_end
  )
}

validate_cutoff <- function(cutoff) {
  expected <- sum(theta_matrix >= cutoff)
  observed <- theta_long[theta >= cutoff, .N]
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(
      "Tile-count validation failed for cutoff ",
      cutoff,
      ": expected ",
      expected,
      ", observed ",
      observed,
      ".",
      call. = FALSE
    )
  }

  expected_cells <- theta_long[
    theta >= cutoff,
    .(expected = .N),
    by = .(tf, topic_num)
  ]
  observed_cells <- data.table::rbindlist(
    lapply(seq_len(n_pages), function(page_number) {
      make_page_data(cutoff, page_number)$tiles[
        ,
        .(observed = .N),
        by = .(tf, topic_num)
      ]
    })
  )
  checked <- merge(
    expected_cells,
    observed_cells,
    by = c("tf", "topic_num"),
    all = TRUE
  )
  checked[is.na(expected), expected := 0L]
  checked[is.na(observed), observed := 0L]
  if (checked[expected != observed, .N]) {
    stop(
      "Per-cell tile-count validation failed for cutoff ",
      cutoff,
      ".",
      call. = FALSE
    )
  }
  as.integer(observed)
}

tile_counts <- vapply(cutoffs, validate_cutoff, integer(1L))
names(tile_counts) <- format(cutoffs, nsmall = 1)
border_count <- sum(vapply(seq_len(n_pages), function(page_number) {
  nrow(make_page_data(cutoffs[[1L]], page_number)$term_borders)
}, integer(1L)))
if (!identical(as.integer(border_count), as.integer(nrow(tf_term_assignments)))) {
  stop(
    "TF-term border validation failed: expected ",
    nrow(tf_term_assignments),
    ", observed ",
    border_count,
    ".",
    call. = FALSE
  )
}

plot_page <- function(cutoff, page_number) {
  page_data <- make_page_data(cutoff, page_number)
  page_rows <- page_data$page_rows
  first_tf <- page_rows$tf[[1L]]
  last_tf <- page_rows$tf[[nrow(page_rows)]]

  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = page_data$cells,
      ggplot2::aes(x = topic_num, y = tf_row),
      width = 0.96,
      height = 0.94,
      fill = "white",
      color = "#E2E7E4",
      linewidth = 0.15
    ) +
    ggplot2::geom_tile(
      data = page_data$tiles,
      ggplot2::aes(x = x, y = y, fill = theta),
      width = 0.15,
      height = 0.36,
      color = "white",
      linewidth = 0.08
    ) +
    ggplot2::geom_tile(
      data = page_data$term_borders,
      ggplot2::aes(x = topic_num, y = tf_row),
      width = 0.96,
      height = 0.94,
      fill = NA,
      color = "#D62728",
      linewidth = 0.85
    ) +
    ggplot2::scale_x_continuous(
      breaks = topic_numbers,
      labels = paste0("T", topic_numbers),
      limits = range(topic_numbers) + c(-0.55, 0.55),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = page_rows$tf_row,
      labels = page_rows$tf,
      limits = c(0.45, nrow(page_rows) + 0.55),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#F7FBFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"),
      values = scales::rescale(c(0, 0.3, 0.5, 0.75, 1)),
      limits = c(0, 1),
      breaks = c(0, 0.3, 0.5, 0.7, 1),
      oob = scales::squish,
      name = "Raw theta"
    ) +
    ggplot2::labs(
      title = sprintf(
        "Raw condition::TF theta by TF and topic | theta >= %.1f",
        cutoff
      ),
      subtitle = sprintf(
        "TFs %d-%d of %d (%s to %s) | page %d of %d",
        page_data$tf_start,
        page_data$tf_end,
        length(tf_levels),
        first_tf,
        last_tf,
        page_number,
        n_pages
      ),
      x = "Topic",
      y = "TF",
      caption = paste0(
        "One tile per passing condition::TF document; tile fill is its raw theta. ",
        "Red border: GENE:<TF> final GammaFit+MaxProb topic assignment. ",
        "No theta aggregation. ",
        length(condition_levels),
        " conditions."
      )
    ) +
    ggplot2::theme_minimal(base_size = 9, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 9,
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 10,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 9,
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
      plot.title = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 9,
        hjust = 0.5
      ),
      plot.caption = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 9,
        hjust = 0
      ),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 9
      ),
      legend.text = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = 9
      ),
      legend.key.height = grid::unit(0.38, "in"),
      plot.margin = ggplot2::margin(10, 12, 8, 10)
    )
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
grDevices::cairo_pdf(
  output_file,
  width = 16,
  height = page_height,
  family = "Helvetica",
  onefile = TRUE,
  bg = "white"
)
device_open <- TRUE
on.exit(if (isTRUE(device_open)) grDevices::dev.off(), add = TRUE)
for (cutoff in cutoffs) {
  for (page_number in seq_len(n_pages)) {
    print(plot_page(cutoff, page_number))
  }
}
grDevices::dev.off()
device_open <- FALSE

message("Wrote: ", output_file)
message(
  "Documents: ",
  nrow(theta_table),
  "; conditions: ",
  length(condition_levels),
  "; TFs: ",
  length(tf_levels),
  "; topics: ",
  length(topic_names)
)
message(
  "Passing tiles: ",
  paste0(names(tile_counts), "=", tile_counts, collapse = "; ")
)
message("TF-term topic borders: ", border_count)
