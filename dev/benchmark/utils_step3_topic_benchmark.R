## Benchmark helpers for Module 3 topic models
##
## This file provides light-weight utilities for inspecting topic-term
## probabilities and peak->gene concordance for specific terms.

.find_topic_model_subdir <- function(model_dir,
                                     backend = c("vae", "warplda"),
                                     vae_variant = "multivi_encoder",
                                     doc_mode = c("tf_cluster", "tf")) {
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  subdirs <- list.dirs(model_dir, recursive = FALSE, full.names = TRUE)
  if (backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_Kgrid$")
  } else {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_Kgrid$")
  }
  subdirs <- subdirs[grepl(patt, basename(subdirs))]
  if (!length(subdirs)) return(NULL)
  subdirs[[1]]
}

.load_topic_phi <- function(model_dir,
                            k,
                            backend = c("vae", "warplda"),
                            vae_variant = "multivi_encoder",
                            doc_mode = c("tf_cluster", "tf")) {
  backend <- match.arg(backend)
  subdir <- .find_topic_model_subdir(model_dir, backend = backend, vae_variant = vae_variant, doc_mode = doc_mode)
  if (is.null(subdir)) {
    .log_abort("No topic model directory found in {model_dir}.")
  }
  phi_path <- file.path(subdir, "vae_models", sprintf("phi_K%d.csv", as.integer(k)))
  if (!file.exists(phi_path)) {
    .log_abort("Missing phi for K={k}: {phi_path}")
  }
  if (!file.exists(phi_path)) {
    .log_abort("Missing phi file: {phi_path}")
  }
  phi_df <- data.table::fread(phi_path, data.table = FALSE)
  if (ncol(phi_df) < 2) {
    .log_abort("phi CSV has fewer than 2 columns: {phi_path}")
  }
  if (is.na(names(phi_df)[1]) || names(phi_df)[1] == "") {
    names(phi_df)[1] <- "term_id"
  }
  phi <- as.matrix(phi_df[, -1, drop = FALSE])
  rownames(phi) <- phi_df[[1]]
  list(phi = phi, phi_path = phi_path)
}

.load_edges_docs <- function(model_dir,
                             backend = c("vae", "warplda"),
                             vae_variant = "multivi_encoder",
                             doc_mode = c("tf_cluster", "tf")) {
  backend <- match.arg(backend)
  subdir <- .find_topic_model_subdir(model_dir, backend = backend, vae_variant = vae_variant, doc_mode = doc_mode)
  if (is.null(subdir)) {
    .log_abort("No topic model directory found in {model_dir}.")
  }
  edges_path <- file.path(subdir, "rds", "edges_docs.rds")
  if (!file.exists(edges_path)) {
    .log_abort("Missing edges_docs: {edges_path}")
  }
  readRDS(edges_path)
}

plot_peak_gene_topic_concordance <- function(peak_id,
                                             model_dir,
                                             k,
                                             backend = c("vae", "warplda"),
                                             vae_variant = "multivi_encoder",
                                             topic_id = NULL,
                                             mode = c("topic_means", "single_topic"),
                                             out_file = NULL) {
  backend <- match.arg(backend)
  mode <- match.arg(mode)

  phi_obj <- .load_topic_phi(model_dir, k, backend = backend, vae_variant = vae_variant, doc_mode = doc_mode)
  phi <- phi_obj$phi
  edges_docs <- .load_edges_docs(model_dir, backend = backend, vae_variant = vae_variant, doc_mode = doc_mode)

  peak_term <- if (grepl("^PEAK:", peak_id)) peak_id else paste0("PEAK:", peak_id)
  if (!peak_term %in% rownames(phi)) {
    .log_abort("Peak term not found in phi: {peak_term}")
  }

  ed <- data.table::as.data.table(edges_docs)
  if (!all(c("peak_id", "gene_key") %in% names(ed))) {
    .log_abort("edges_docs missing peak_id/gene_key columns.")
  }
  genes <- unique(ed[peak_id == gsub("^PEAK:", "", peak_term), gene_key])
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) {
    .log_abort("No connected genes found for peak {peak_id}.")
  }
  gene_terms <- paste0("GENE:", genes)
  gene_terms <- gene_terms[gene_terms %in% rownames(phi)]
  if (!length(gene_terms)) {
    .log_abort("Connected genes not found in phi term list.")
  }

  topic_cols <- colnames(phi)
  if (is.null(topic_cols) || !length(topic_cols)) {
    topic_cols <- paste0("Topic", seq_len(ncol(phi)))
    colnames(phi) <- topic_cols
  }

  if (mode == "single_topic") {
    if (is.null(topic_id)) {
      topic_id <- which.max(phi[peak_term, ])
    }
    topic_col <- if (paste0("Topic", topic_id) %in% topic_cols) paste0("Topic", topic_id) else topic_cols[[topic_id]]
    peak_prob <- as.numeric(phi[peak_term, topic_col])
    gene_prob <- as.numeric(phi[gene_terms, topic_col])
    df <- data.frame(
      gene = sub("^GENE:", "", gene_terms),
      peak_prob = peak_prob,
      gene_prob = gene_prob,
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = peak_prob, y = gene_prob)) +
      ggplot2::geom_jitter(width = 0.002, height = 0, alpha = 0.7, size = 1.5) +
      ggplot2::labs(
        title = paste0("Peak vs gene topic prob (Topic ", topic_id, ")"),
        x = "Peak topic probability",
        y = "Gene topic probability"
      ) +
      ggplot2::theme_classic(base_size = 11)
  } else {
    peak_probs <- as.numeric(phi[peak_term, ])
    gene_mat <- phi[gene_terms, , drop = FALSE]
    gene_means <- colMeans(gene_mat, na.rm = TRUE)
    df <- data.frame(
      topic = seq_along(peak_probs),
      peak_prob = peak_probs,
      gene_prob = gene_means,
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = peak_prob, y = gene_prob)) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggrepel::geom_text_repel(ggplot2::aes(label = topic), size = 3, max.overlaps = 50) +
      ggplot2::labs(
        title = "Peak vs mean gene topic probability (per topic)",
        x = "Peak topic probability",
        y = "Mean gene topic probability"
      ) +
      ggplot2::theme_classic(base_size = 11)
  }

  if (!is.null(out_file) && nzchar(out_file)) {
    ggplot2::ggsave(out_file, p, width = 5.5, height = 4.2, dpi = 200)
  }
  p
}

plot_peak_gene_concordance_all <- function(model_dir,
                                           k,
                                           backend = c("vae", "warplda"),
                                           vae_variant = "multivi_encoder",
                                           doc_mode = c("tf_cluster", "tf"),
                                           topic_terms_path = NULL,
                                           final_topics_dir = NULL,
                                           out_file = NULL,
                                           point_alpha = 0.05,
                                           point_size = 0.4,
                                           facets_ncol = 5,
                                           facets_nrow = 3,
                                           verbose = TRUE) {
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"

  if (isTRUE(verbose)) {
    message("[topic_benchmark] Loading phi for K=", k, " (", backend, "/", vae_variant, ") from ", model_dir)
  }
  phi_obj <- .load_topic_phi(model_dir, k, backend = backend, vae_variant = vae_variant, doc_mode = doc_mode)
  phi <- phi_obj$phi
  if (isTRUE(verbose)) {
    message("[topic_benchmark] raw phi dims: ", nrow(phi), " x ", ncol(phi))
  }
  term_in_rows <- any(grepl("^(PEAK:|GENE:)", rownames(phi)))
  term_in_cols <- any(grepl("^(PEAK:|GENE:)", colnames(phi)))
  if (!term_in_rows && term_in_cols) {
    phi <- t(phi)
  }
  if (isTRUE(verbose)) {
    message("[topic_benchmark] phi dims (terms x topics): ", nrow(phi), " x ", ncol(phi))
  }
  edges_docs <- .load_edges_docs(model_dir, backend = backend, vae_variant = vae_variant, doc_mode = doc_mode)
  if (isTRUE(verbose)) {
    message("[topic_benchmark] edges_docs rows: ", nrow(edges_docs))
    message("[topic_benchmark] edges_docs cols: ", paste(names(edges_docs), collapse = ", "))
  }

  if (!all(c("peak_id", "gene_key") %in% names(edges_docs))) {
    .log_abort("edges_docs missing peak_id/gene_key columns.")
  }

  peak_terms <- paste0("PEAK:", edges_docs$peak_id)
  gene_terms <- paste0("GENE:", edges_docs$gene_key)
  if (isTRUE(verbose)) {
    message("[topic_benchmark] example peak term: ", head(peak_terms, 1))
    message("[topic_benchmark] example gene term: ", head(gene_terms, 1))
  }
  keep <- peak_terms %in% rownames(phi) & gene_terms %in% rownames(phi)
  if (isTRUE(verbose)) {
    message("[topic_benchmark] pairs total: ", length(keep), "; kept: ", sum(keep))
    message("[topic_benchmark] phi term prefix counts: PEAK=",
            sum(grepl("^PEAK:", rownames(phi))), ", GENE=", sum(grepl("^GENE:", rownames(phi))))
  }
  if (!any(keep)) {
    .log_abort("No peak-gene pairs found in phi term list.")
  }

  peak_terms <- peak_terms[keep]
  gene_terms <- gene_terms[keep]

  peak_mat <- phi[peak_terms, , drop = FALSE]
  gene_mat <- phi[gene_terms, , drop = FALSE]
  topic_ids <- seq_len(ncol(phi))

  topic_terms <- NULL
  link_topic_pairs <- NULL
  link_topic_df <- NULL
  if (is.null(topic_terms_path) && !is.null(final_topics_dir) && nzchar(final_topics_dir)) {
    out_dirs <- list.dirs(final_topics_dir, recursive = FALSE, full.names = TRUE)
    if (backend == "vae") {
      patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", as.integer(k), "$")
      out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
    } else {
      patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", as.integer(k), "$")
      out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
    }
    if (length(out_dirs)) {
      topic_terms_path <- file.path(out_dirs[[1]], "topic_terms.csv")
    }
  }
  if (!is.null(topic_terms_path) && nzchar(topic_terms_path) && file.exists(topic_terms_path)) {
    topic_terms <- readr::read_csv(topic_terms_path, show_col_types = FALSE, name_repair = "minimal")
    if (!all(c("topic", "term_id", "in_topic") %in% names(topic_terms))) {
      .log_abort("topic_terms.csv must include topic, term_id, in_topic.")
    }
    topic_terms <- topic_terms[isTRUE(topic_terms$in_topic), , drop = FALSE]
    if (isTRUE(verbose)) {
      message("[topic_benchmark] Using topic_terms filter from ", topic_terms_path)
    }
  }
  if (!is.null(final_topics_dir) && nzchar(final_topics_dir) && length(out_dirs)) {
    link_path <- file.path(out_dirs[[1]], "topic_links.csv")
    if (file.exists(link_path)) {
      link_df <- data.table::fread(link_path)
      if (all(c("peak_id", "gene_key", "topic_num", "peak_pass", "gene_pass") %in% names(link_df))) {
        link_df <- link_df[isTRUE(link_df$peak_pass) & isTRUE(link_df$gene_pass), , drop = FALSE]
        link_df$peak_term <- paste0("PEAK:", link_df$peak_id)
        link_df$gene_term <- paste0("GENE:", link_df$gene_key)
        link_topic_pairs <- split(paste0(link_df$peak_term, "||", link_df$gene_term), link_df$topic_num)
        link_topic_df <- link_df[, c("topic_num", "peak_term", "gene_term")]
        if (isTRUE(verbose)) {
          message("[topic_benchmark] Using link-topic pairs from ", link_path)
        }
      }
    }
  }

  .build_df <- function(use_topic_terms = TRUE) {
    df_list <- vector("list", length(topic_ids))
    for (t in topic_ids) {
      if (isTRUE(use_topic_terms) && !is.null(topic_terms)) {
        in_terms <- topic_terms$term_id[topic_terms$topic == t]
        peaks_in <- in_terms[grepl("^PEAK:", in_terms)]
        genes_in <- in_terms[grepl("^GENE:", in_terms)]
        keep_t <- (peak_terms %in% peaks_in) & (gene_terms %in% genes_in)
        keep_t[is.na(keep_t)] <- FALSE
      } else if (!is.null(link_topic_pairs)) {
        pair_key <- paste0(peak_terms, "||", gene_terms)
        pairs_in <- link_topic_pairs[[as.character(t)]]
        keep_t <- pair_key %in% pairs_in
        keep_t[is.na(keep_t)] <- FALSE
      } else {
        keep_t <- rep(TRUE, length(peak_terms))
      }
      idx <- which(keep_t)
      if (!length(idx)) {
        df_list[[t]] <- data.frame(
          topic = integer(0),
          peak_prob = numeric(0),
          gene_prob = numeric(0)
        )
      } else {
        df_list[[t]] <- data.frame(
          topic = rep(t, length(idx)),
          peak_prob = as.numeric(peak_mat[idx, t]),
          gene_prob = as.numeric(gene_mat[idx, t])
        )
      }
    }
    do.call(rbind, df_list)
  }

  if (!is.null(link_topic_df)) {
    link_topic_df$topic_num <- as.integer(link_topic_df$topic_num)
    keep_topic <- link_topic_df$topic_num >= 1L & link_topic_df$topic_num <= ncol(phi)
    keep_peak <- link_topic_df$peak_term %in% rownames(phi)
    keep_gene <- link_topic_df$gene_term %in% rownames(phi)
    link_topic_df <- link_topic_df[keep_topic & keep_peak & keep_gene, , drop = FALSE]
    if (!nrow(link_topic_df)) {
      .log_abort("No link-topic pairs found in phi term list.")
    }
    peak_idx <- match(link_topic_df$peak_term, rownames(phi))
    gene_idx <- match(link_topic_df$gene_term, rownames(phi))
    topic_idx <- link_topic_df$topic_num
    df <- data.frame(
      topic = topic_idx,
      peak_prob = phi[cbind(peak_idx, topic_idx)],
      gene_prob = phi[cbind(gene_idx, topic_idx)]
    )
  } else {
    df <- .build_df(use_topic_terms = TRUE)
  }
  if (!nrow(df)) {
    .log_abort("No peak-gene pairs remain after link-topic filtering.")
  }

  xlim <- range(log10(df$peak_prob + 1e-6), finite = TRUE)
  ylim <- range(log10(df$gene_prob + 1e-6), finite = TRUE)
  if (!all(is.finite(xlim))) xlim <- c(-6, 0)
  if (!all(is.finite(ylim))) ylim <- c(-6, 0)

  per_page <- facets_ncol * facets_nrow
  pages <- split(topic_ids, ceiling(seq_along(topic_ids) / per_page))

  if (!is.null(out_file) && nzchar(out_file)) {
    grDevices::pdf(out_file, width = 11, height = 8.5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  for (pg in seq_along(pages)) {
    topics_pg <- pages[[pg]]
    df_pg <- df[df$topic %in% topics_pg, , drop = FALSE]
    df_all <- df
    df_all$topic <- "ALL"
    df_pg$topic <- as.character(df_pg$topic)
    df_plot <- rbind(
      df_pg,
      df_all
    )
    stats_df <- lapply(split(df_plot, df_plot$topic), function(d) {
      if (nrow(d) < 3) {
        return(data.frame(topic = unique(d$topic), label = "r=NA, p=NA"))
      }
      ct <- stats::cor.test(d$peak_prob, d$gene_prob, method = "pearson")
      data.frame(
        topic = unique(d$topic),
        label = sprintf("r=%.2f, p=%.2g", ct$estimate, ct$p.value)
      )
    })
    stats_df <- do.call(rbind, stats_df)
    df_plot$peak_prob_log10 <- log10(df_plot$peak_prob + 1e-6)
    df_plot$gene_prob_log10 <- log10(df_plot$gene_prob + 1e-6)
    stats_df$label <- stats_df$label
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = peak_prob_log10, y = gene_prob_log10)) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::facet_wrap(~ topic, ncol = facets_ncol, nrow = facets_nrow) +
      ggplot2::labs(
        title = "All peak-gene topic probabilities",
        x = "log10(PEAK topic probability + 1e-6)",
        y = "log10(GENE topic probability + 1e-6)"
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        strip.text = ggplot2::element_text(size = 8)
      ) +
      ggplot2::geom_text(
        data = stats_df,
        ggplot2::aes(x = xlim[1], y = ylim[2], label = label),
        inherit.aes = FALSE,
        hjust = 0, vjust = 1, size = 2.4, fontface = "bold"
      )
    if (!is.null(out_file) && nzchar(out_file)) {
      print(p)
    }
  }

  if (is.null(out_file) || !nzchar(out_file)) {
    return(p)
  }
  invisible(NULL)
}

.build_peak_gene_df <- function(model_dir,
                                k,
                                backend,
                                vae_variant,
                                final_topics_dir,
                                doc_mode = c("tf_cluster", "tf"),
                                verbose = TRUE) {
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  phi_obj <- .load_topic_phi(model_dir, k, backend = backend, vae_variant = vae_variant)
  phi <- phi_obj$phi

  term_in_rows <- any(grepl("^(PEAK:|GENE:)", rownames(phi)))
  term_in_cols <- any(grepl("^(PEAK:|GENE:)", colnames(phi)))
  if (!term_in_rows && term_in_cols) {
    phi <- t(phi)
  }

  edges_docs <- .load_edges_docs(model_dir, backend = backend, vae_variant = vae_variant)
  if (!all(c("peak_id", "gene_key") %in% names(edges_docs))) {
    .log_abort("edges_docs missing peak_id/gene_key columns.")
  }
  peak_terms <- paste0("PEAK:", edges_docs$peak_id)
  gene_terms <- paste0("GENE:", edges_docs$gene_key)
  keep <- peak_terms %in% rownames(phi) & gene_terms %in% rownames(phi)
  peak_terms <- peak_terms[keep]
  gene_terms <- gene_terms[keep]
  if (!length(peak_terms)) {
    .log_abort("No peak-gene pairs found in phi term list.")
  }

  out_dirs <- list.dirs(final_topics_dir, recursive = FALSE, full.names = TRUE)
  if (backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", as.integer(k), "$")
  } else {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", as.integer(k), "$")
  }
  out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
  if (!length(out_dirs)) {
    .log_abort("No final_topics subdir found for K={k} in {final_topics_dir}")
  }
  link_path <- file.path(out_dirs[[1]], "topic_links.csv")
  if (!file.exists(link_path)) {
    .log_abort("Missing topic_links.csv at {link_path}")
  }

  link_df <- data.table::fread(link_path)
  if (!all(c("peak_id", "gene_key", "topic_num", "peak_pass", "gene_pass") %in% names(link_df))) {
    .log_abort("topic_links.csv missing required columns.")
  }
  link_df <- link_df[isTRUE(link_df$peak_pass) & isTRUE(link_df$gene_pass), , drop = FALSE]
  link_df$topic_num <- as.integer(link_df$topic_num)
  link_df$peak_term <- paste0("PEAK:", link_df$peak_id)
  link_df$gene_term <- paste0("GENE:", link_df$gene_key)

  keep_topic <- link_df$topic_num >= 1L & link_df$topic_num <= ncol(phi)
  keep_peak <- link_df$peak_term %in% rownames(phi)
  keep_gene <- link_df$gene_term %in% rownames(phi)
  link_df <- link_df[keep_topic & keep_peak & keep_gene, , drop = FALSE]
  if (!nrow(link_df)) {
    .log_abort("No link-topic pairs found in phi term list.")
  }

  peak_idx <- match(link_df$peak_term, rownames(phi))
  gene_idx <- match(link_df$gene_term, rownames(phi))
  topic_idx <- link_df$topic_num

  df <- data.frame(
    topic = topic_idx,
    peak_prob = phi[cbind(peak_idx, topic_idx)],
    gene_prob = phi[cbind(gene_idx, topic_idx)]
  )

  xlim <- range(log10(df$peak_prob + 1e-6), finite = TRUE)
  ylim <- range(log10(df$gene_prob + 1e-6), finite = TRUE)
  if (!all(is.finite(xlim))) xlim <- c(-6, 0)
  if (!all(is.finite(ylim))) ylim <- c(-6, 0)

  list(df = df, xlim = xlim, ylim = ylim)
}

plot_peak_gene_concordance_all_methods <- function(combo_grid,
                                                   step3_root_dir,
                                                   final_topics_subdir = "final_topics",
                                                   k,
                                                   vae_variant = "multivi_encoder",
                                                   doc_mode = c("tf_cluster", "tf"),
                                                   out_file,
                                                   point_alpha = 0.05,
                                                   point_size = 0.4,
                                                   facets_ncol = 6,
                                                   density_bins = 140L,
                                                   methods = c("peak_and_gene", "link_score_prob"),
                                                   verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  .assert_pkg("readr")
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"

  methods <- intersect(methods, c("peak_and_gene", "link_score_prob"))
  if (!length(methods)) {
    .log_abort("methods must include peak_and_gene and/or link_score_prob.")
  }
  if (missing(out_file) || is.null(out_file) || !nzchar(out_file)) {
    .log_abort("out_file must be provided for peak-gene concordance outputs.")
  }
  out_base <- tools::file_path_sans_ext(out_file)
  out_ext <- tools::file_ext(out_file)
  if (!nzchar(out_ext)) out_ext <- "pdf"
  method_out_files <- setNames(
    vapply(methods, function(m) sprintf("%s_%s.%s", out_base, m, out_ext), character(1)),
    methods
  )
  method_out_files_density <- setNames(
    vapply(methods, function(m) sprintf("%s_rrho2_like_%s.%s", out_base, m, out_ext), character(1)),
    methods
  )
  out_file_density_all <- sprintf("%s_rrho2_like_all.%s", out_base, out_ext)

  combo_grid$combo_tag <- paste0(
    "gene_", combo_grid$gene_term_mode,
    "_tf_", ifelse(combo_grid$include_tf_terms, "on", "off"),
    "_count_", combo_grid$count_input
  )
  combo_grid$combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_grid$combo_tag)

  # Match the ordering used by histogram benchmark for direct visual comparison.
  combo_grid$backend <- factor(combo_grid$backend, levels = c("warplda", "vae"))
  combo_grid$gene_term_mode <- factor(combo_grid$gene_term_mode, levels = c("aggregate", "unique"))
  combo_grid$include_tf_terms <- factor(combo_grid$include_tf_terms, levels = c(FALSE, TRUE))
  combo_grid$count_input <- factor(combo_grid$count_input, levels = c("weight", "pseudo_count_bin", "pseudo_count_log"))
  combo_grid <- combo_grid[order(combo_grid$backend, combo_grid$gene_term_mode, combo_grid$include_tf_terms, combo_grid$count_input), , drop = FALSE]
  combo_grid$backend <- as.character(combo_grid$backend)
  combo_grid$gene_term_mode <- as.character(combo_grid$gene_term_mode)
  combo_grid$include_tf_terms <- as.logical(as.character(combo_grid$include_tf_terms))
  combo_grid$count_input <- as.character(combo_grid$count_input)

  .row_label <- function(row) {
    backend_lbl <- if (identical(row$backend, "warplda")) "LDA" else "MultiVI"
    gene_lbl <- if (identical(row$gene_term_mode, "aggregate")) "Aggregate" else "Distinct"
    tf_lbl <- if (isTRUE(row$include_tf_terms)) "TF_on" else "TF_off"
    count_lbl <- switch(
      row$count_input,
      weight = "Weight",
      pseudo_count_bin = "Bin",
      pseudo_count_log = "Log",
      row$count_input
    )
    paste(backend_lbl, gene_lbl, tf_lbl, count_lbl)
  }

  panel_levels <- vapply(seq_len(nrow(combo_grid)), function(i) .row_label(combo_grid[i, , drop = FALSE]), character(1))
  .as_flag <- function(x) {
    if (is.logical(x)) return(!is.na(x) & x)
    if (is.numeric(x) || is.integer(x)) return(is.finite(x) & x != 0)
    x <- tolower(trimws(as.character(x)))
    x %in% c("true", "t", "1", "yes", "y")
  }

  all_data <- vector("list", 0L)
  all_data_unfiltered <- vector("list", 0L)
  for (method in methods) {
    for (i in seq_len(nrow(combo_grid))) {
      row <- combo_grid[i, , drop = FALSE]
      final_dir <- file.path(step3_root_dir, final_topics_subdir, row$combo_tag)
      out_dirs <- list.dirs(final_dir, recursive = FALSE, full.names = TRUE)
      if (row$backend == "vae") {
        patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", as.integer(k), "$")
      } else {
        patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", as.integer(k), "$")
      }
      out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
      if (!length(out_dirs)) {
        if (isTRUE(verbose)) .log_warn("[topic_benchmark] Missing final_topics for {row$combo_tag}")
        next
      }
      link_path <- file.path(out_dirs[[1]], "topic_links.csv")
      if (!file.exists(link_path)) {
        if (isTRUE(verbose)) .log_warn("[topic_benchmark] Missing topic_links.csv for {row$combo_tag}")
        next
      }
      dt <- data.table::fread(link_path)
      dt_all <- data.table::copy(dt)
      req <- c("peak_id", "gene_key", "peak_score", "gene_score", "peak_pass", "gene_pass")
      if (!all(req %in% names(dt))) {
        if (isTRUE(verbose)) .log_warn("[topic_benchmark] topic_links missing required columns for {row$combo_tag}")
        next
      }
      dt_all[, peak_score := suppressWarnings(as.numeric(peak_score))]
      dt_all[, gene_score := suppressWarnings(as.numeric(gene_score))]
      dt_all <- dt_all[is.finite(peak_score) & is.finite(gene_score)]
      if (nrow(dt_all)) {
        dt_all <- dt_all[, .(
          peak_score = max(peak_score, na.rm = TRUE),
          gene_score = max(gene_score, na.rm = TRUE)
        ), by = .(peak_id, gene_key)]
        dt_all <- dt_all[is.finite(peak_score) & is.finite(gene_score)]
        if (nrow(dt_all)) {
          dt_all[, panel := .row_label(row)]
          all_data_unfiltered[[length(all_data_unfiltered) + 1L]] <- dt_all[, .(panel, peak_score, gene_score)]
        }
      }
      if (method == "peak_and_gene") {
        dt <- dt[.as_flag(peak_pass) & .as_flag(gene_pass)]
      } else {
        if (!("link_pass" %in% names(dt))) {
          if (isTRUE(verbose)) .log_warn("[topic_benchmark] topic_links missing link_pass for {row$combo_tag}")
          next
        }
        dt <- dt[.as_flag(link_pass)]
      }
      if (!nrow(dt)) next
      dt[, peak_score := suppressWarnings(as.numeric(peak_score))]
      dt[, gene_score := suppressWarnings(as.numeric(gene_score))]
      dt <- dt[is.finite(peak_score) & is.finite(gene_score)]
      if (!nrow(dt)) next
      # One point per unique peak-target link.
      dt <- dt[, .(
        peak_score = max(peak_score, na.rm = TRUE),
        gene_score = max(gene_score, na.rm = TRUE)
      ), by = .(peak_id, gene_key)]
      dt <- dt[is.finite(peak_score) & is.finite(gene_score)]
      if (!nrow(dt)) next
      dt[, method := method]
      dt[, panel := .row_label(row)]
      all_data[[length(all_data) + 1L]] <- dt[, .(method, panel, peak_score, gene_score)]
    }
  }

  if (!length(all_data)) {
    .log_abort("No peak-gene scores loaded from topic_links.csv.")
  }

  plot_df <- data.table::rbindlist(all_data, use.names = TRUE, fill = TRUE)
  plot_df <- plot_df[is.finite(peak_score) & is.finite(gene_score)]
  if (!nrow(plot_df)) {
    .log_abort("No finite peak/gene scores to plot.")
  }
  plot_df[, peak_score_log10 := suppressWarnings(log10(peak_score))]
  plot_df[, gene_score_log10 := suppressWarnings(log10(gene_score))]
  plot_df <- plot_df[is.finite(peak_score_log10) & is.finite(gene_score_log10)]
  if (!nrow(plot_df)) {
    .log_abort("No finite log10(peak_score)/log10(gene_score) values to plot.")
  }
  plot_df_all <- NULL
  if (length(all_data_unfiltered)) {
    plot_df_all <- data.table::rbindlist(all_data_unfiltered, use.names = TRUE, fill = TRUE)
    plot_df_all <- plot_df_all[is.finite(peak_score) & is.finite(gene_score)]
    if (nrow(plot_df_all)) {
      plot_df_all[, peak_score_log10 := suppressWarnings(log10(peak_score))]
      plot_df_all[, gene_score_log10 := suppressWarnings(log10(gene_score))]
      plot_df_all <- plot_df_all[is.finite(peak_score_log10) & is.finite(gene_score_log10)]
    } else {
      plot_df_all <- NULL
    }
  }

  # Shared axes across methods for direct comparability.
  xlim <- range(plot_df$peak_score_log10, finite = TRUE)
  ylim <- range(plot_df$gene_score_log10, finite = TRUE)
  if (!all(is.finite(xlim))) xlim <- c(-6, 0)
  if (!all(is.finite(ylim))) ylim <- c(-6, 0)
  if (data.table::is.data.table(plot_df_all) && nrow(plot_df_all)) {
    xlim_all <- range(plot_df_all$peak_score_log10, finite = TRUE)
    ylim_all <- range(plot_df_all$gene_score_log10, finite = TRUE)
    if (all(is.finite(xlim_all))) xlim <- range(c(xlim, xlim_all), finite = TRUE)
    if (all(is.finite(ylim_all))) ylim <- range(c(ylim, ylim_all), finite = TRUE)
  }
  # Use a shared square bin width in data units for RRHO2-like density panels.
  density_bins <- max(10L, as.integer(density_bins))
  xr <- diff(xlim)
  yr <- diff(ylim)
  max_range <- max(xr, yr)
  square_binwidth <- max_range / density_bins
  if (!is.finite(square_binwidth) || square_binwidth <= 0) {
    square_binwidth <- 0.1
  }
  n_panels <- length(panel_levels)
  nrow_facets <- ceiling(n_panels / as.integer(facets_ncol))
  width <- max(12, as.integer(facets_ncol) * 2.4)
  height <- max(8, nrow_facets * 2.4)

  for (m in methods) {
    dsub <- plot_df[method == m]
    if (!nrow(dsub)) {
      if (isTRUE(verbose)) .log_warn("[topic_benchmark] No data for method={m} in peak-gene concordance.")
      next
    }
    dsub[, panel := factor(panel, levels = panel_levels)]
    stats_dt <- dsub[, {
      x <- as.numeric(peak_score)
      y <- as.numeric(gene_score)
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]
      y <- y[ok]
      n <- length(x)
      if (n >= 3L && stats::sd(x) > 0 && stats::sd(y) > 0) {
        r <- suppressWarnings(stats::cor(x, y, method = "pearson"))
        if (is.finite(r) && abs(r) < 1) {
          tval <- r * sqrt((n - 2) / (1 - r^2))
          pval <- 2 * stats::pt(-abs(tval), df = n - 2)
          lab <- sprintf("R=%.2f, p=%s", r, formatC(pval, format = "e", digits = 2))
        } else if (is.finite(r) && abs(r) == 1) {
          lab <- sprintf("R=%.2f, p=0", r)
        } else {
          lab <- "R=NA, p=NA"
        }
      } else {
        lab <- "R=NA, p=NA"
      }
      .(label = lab)
    }, by = panel]
    stats_dt[, panel := factor(panel, levels = levels(dsub$panel))]
    x_annot <- xlim[1] + 0.02 * diff(xlim)
    y_annot <- ylim[2] - 0.02 * diff(ylim)
    stats_dt[, `:=`(x = x_annot, y = y_annot)]

    p <- ggplot2::ggplot(dsub, ggplot2::aes(x = peak_score_log10, y = gene_score_log10)) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::facet_wrap(~ panel, ncol = as.integer(facets_ncol), drop = FALSE) +
      ggplot2::geom_text(
        data = stats_dt,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 2.4,
        fontface = "bold"
      ) +
      ggplot2::labs(
        title = paste0("Peak-gene score concordance (", m, ")"),
        x = "log10(peak_score)",
        y = "log10(gene_score)"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(size = 8, face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold")
      )
    ggplot2::ggsave(method_out_files[[m]], p, width = width, height = height, limitsize = FALSE)
    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote peak-gene concordance: {method_out_files[[m]]}")

    # RRHO2-like density view on the same transformed scores and panel layout.
    rrho2_palette <- c("#4815AA", "#3783FF", "#4DE94C", "#FFEE00", "#FF8C00", "#F60000")
    p_density <- ggplot2::ggplot(dsub, ggplot2::aes(x = peak_score_log10, y = gene_score_log10)) +
      ggplot2::stat_bin_2d(binwidth = c(square_binwidth, square_binwidth), na.rm = TRUE) +
      ggplot2::scale_fill_gradientn(colours = rrho2_palette, trans = "sqrt") +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::facet_wrap(~ panel, ncol = as.integer(facets_ncol), drop = FALSE) +
      ggplot2::geom_text(
        data = stats_dt,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 2.4,
        fontface = "bold"
      ) +
      ggplot2::labs(
        title = paste0("Peak-gene RRHO2-like density (", m, ")"),
        x = "log10(peak_score)",
        y = "log10(gene_score)",
        fill = "N"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(size = 8, face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(face = "bold")
      )
    ggplot2::ggsave(method_out_files_density[[m]], p_density, width = width, height = height, limitsize = FALSE)
    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote RRHO2-like density: {method_out_files_density[[m]]}")
  }

  if (data.table::is.data.table(plot_df_all) && nrow(plot_df_all)) {
    dsub <- data.table::copy(plot_df_all)
    dsub[, panel := factor(panel, levels = panel_levels)]
    stats_dt <- dsub[, {
      x <- as.numeric(peak_score)
      y <- as.numeric(gene_score)
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]
      y <- y[ok]
      n <- length(x)
      if (n >= 3L && stats::sd(x) > 0 && stats::sd(y) > 0) {
        r <- suppressWarnings(stats::cor(x, y, method = "pearson"))
        if (is.finite(r) && abs(r) < 1) {
          tval <- r * sqrt((n - 2) / (1 - r^2))
          pval <- 2 * stats::pt(-abs(tval), df = n - 2)
          lab <- sprintf("R=%.2f, p=%s", r, formatC(pval, format = "e", digits = 2))
        } else if (is.finite(r) && abs(r) == 1) {
          lab <- sprintf("R=%.2f, p=0", r)
        } else {
          lab <- "R=NA, p=NA"
        }
      } else {
        lab <- "R=NA, p=NA"
      }
      .(label = lab)
    }, by = panel]
    stats_dt[, panel := factor(panel, levels = levels(dsub$panel))]
    x_annot <- xlim[1] + 0.02 * diff(xlim)
    y_annot <- ylim[2] - 0.02 * diff(ylim)
    stats_dt[, `:=`(x = x_annot, y = y_annot)]

    rrho2_palette <- c("#4815AA", "#3783FF", "#4DE94C", "#FFEE00", "#FF8C00", "#F60000")
    p_density_all <- ggplot2::ggplot(dsub, ggplot2::aes(x = peak_score_log10, y = gene_score_log10)) +
      ggplot2::stat_bin_2d(binwidth = c(square_binwidth, square_binwidth), na.rm = TRUE) +
      ggplot2::scale_fill_gradientn(colours = rrho2_palette, trans = "sqrt") +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::facet_wrap(~ panel, ncol = as.integer(facets_ncol), drop = FALSE) +
      ggplot2::geom_text(
        data = stats_dt,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 2.4,
        fontface = "bold"
      ) +
      ggplot2::labs(
        title = "Peak-gene RRHO2-like density (all links)",
        x = "log10(peak_score)",
        y = "log10(gene_score)",
        fill = "N"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(size = 8, face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(face = "bold")
      )
    ggplot2::ggsave(out_file_density_all, p_density_all, width = width, height = height, limitsize = FALSE)
    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote RRHO2-like density: {out_file_density_all}")
  }

  invisible(NULL)
}

.count_shared_topics <- function(df, item_col, topic_col) {
  if (!is.data.frame(df) || !all(c(item_col, topic_col) %in% names(df))) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- df[!is.na(df[[item_col]]) & df[[item_col]] != "", , drop = FALSE]
  df <- df[!is.na(df[[topic_col]]), , drop = FALSE]
  if (!nrow(df)) return(NULL)
  df <- df[, c(item_col, topic_col), drop = FALSE]
  if (!is.data.frame(df) || ncol(df) < 2L) return(NULL)
  df <- df[!duplicated(df), , drop = FALSE]
  n_topics <- as.integer(table(df[[item_col]]))
  n_topics <- n_topics[n_topics >= 1L]
  if (!length(n_topics)) {
    return(data.frame(n_topics = integer(0), n_items = integer(0), stringsAsFactors = FALSE))
  }
  tab <- table(n_topics)
  out <- data.frame(
    n_topics = as.integer(names(tab)),
    n_items = as.integer(tab),
    stringsAsFactors = FALSE
  )
  out[order(out$n_topics), , drop = FALSE]
}

plot_shared_topic_counts_all_methods <- function(combo_grid,
                                                 step3_root_dir,
                                                 final_topics_subdir = "final_topics",
                                                 k,
                                                 vae_variant = "multivi_encoder",
                                                 doc_mode = c("tf_cluster", "tf"),
                                                 out_file,
                                                 methods = c("peak_and_gene", "link_score_prob"),
                                                 verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"

  methods <- intersect(methods, c("peak_and_gene", "link_score_prob"))
  if (!length(methods)) {
    .log_abort("methods must include peak_and_gene and/or link_score_prob.")
  }
  if (missing(out_file) || is.null(out_file) || !nzchar(out_file)) {
    .log_abort("out_file must be provided for shared-topic benchmark outputs.")
  }

  out_base <- tools::file_path_sans_ext(out_file)
  out_ext <- tools::file_ext(out_file)
  if (!nzchar(out_ext)) out_ext <- "pdf"
  method_out_files <- setNames(
    vapply(methods, function(m) sprintf("%s_%s.%s", out_base, m, out_ext), character(1)),
    methods
  )

  combo_grid$combo_tag <- paste0(
    "gene_", combo_grid$gene_term_mode,
    "_tf_", ifelse(combo_grid$include_tf_terms, "on", "off"),
    "_count_", combo_grid$count_input
  )
  combo_grid$combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_grid$combo_tag)

  # Order for easier comparisons.
  combo_grid$backend <- factor(combo_grid$backend, levels = c("warplda", "vae"))
  combo_grid$gene_term_mode <- factor(combo_grid$gene_term_mode, levels = c("aggregate", "unique"))
  combo_grid$include_tf_terms <- factor(combo_grid$include_tf_terms, levels = c(FALSE, TRUE))
  combo_grid$count_input <- factor(combo_grid$count_input, levels = c("weight", "pseudo_count_bin", "pseudo_count_log"))
  combo_grid <- combo_grid[order(combo_grid$backend, combo_grid$gene_term_mode, combo_grid$include_tf_terms, combo_grid$count_input), , drop = FALSE]
  combo_grid$backend <- as.character(combo_grid$backend)
  combo_grid$gene_term_mode <- as.character(combo_grid$gene_term_mode)
  combo_grid$include_tf_terms <- as.logical(as.character(combo_grid$include_tf_terms))
  combo_grid$count_input <- as.character(combo_grid$count_input)

  for (method in methods) {
    row_payload <- vector("list", nrow(combo_grid))
    for (i in seq_len(nrow(combo_grid))) {
      row <- combo_grid[i, ]
      final_dir <- file.path(step3_root_dir, final_topics_subdir, row$combo_tag)
      out_dirs <- list.dirs(final_dir, recursive = FALSE, full.names = TRUE)
      if (row$backend == "vae") {
        patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", as.integer(k), "$")
      } else {
        patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", as.integer(k), "$")
      }
      out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
      if (!length(out_dirs)) {
        if (isTRUE(verbose)) .log_warn("[topic_benchmark] Missing final_topics for {row$combo_tag}")
        row_payload[[i]] <- list(
          label = paste("Missing:", row$combo_tag),
          links = data.frame(n_topics = integer(0), n_items = integer(0)),
          genes = data.frame(n_topics = integer(0), n_items = integer(0)),
          tfs = data.frame(n_topics = integer(0), n_items = integer(0)),
          paths = data.frame(n_topics = integer(0), n_items = integer(0))
        )
        next
      }
      final_subdir <- out_dirs[[1]]

      link_path <- file.path(final_subdir, "topic_links.csv")
      link_df <- NULL
      if (file.exists(link_path)) {
        link_df <- data.table::fread(link_path)
        if (all(c("peak_id", "gene_key", "topic_num", "peak_pass", "gene_pass") %in% names(link_df))) {
          if (method == "peak_and_gene") {
            link_df <- link_df[.as_logical_flag(peak_pass) & .as_logical_flag(gene_pass)]
          } else {
            if (!("link_pass" %in% names(link_df))) {
              if (isTRUE(verbose)) .log_warn("[topic_benchmark] topic_links missing link_pass for {row$combo_tag}")
              link_df <- NULL
            } else {
              link_df <- link_df[.as_logical_flag(link_pass)]
            }
          }
          link_df$link_id <- paste(link_df$peak_id, link_df$gene_key, sep = "::")
        } else {
          link_df <- NULL
        }
      }

      link_counts <- .count_shared_topics(link_df, "link_id", "topic_num")
      gene_counts <- .count_shared_topics(link_df, "gene_key", "topic_num")

      if (!is.null(link_df) && nrow(link_df) && "tf" %in% names(link_df)) {
        tf_tbl <- link_df[!is.na(tf) & tf != "", .(tf, topic_num)]
        tf_counts <- .count_shared_topics(tf_tbl, "tf", "topic_num")
      } else {
        tf_counts <- data.frame(n_topics = integer(0), n_items = integer(0), stringsAsFactors = FALSE)
      }

      pathway_path <- file.path(final_subdir, paste0("topic_pathway_enrichment_", method, "_dotplot.csv"))
      pathway_counts <- NULL
      if (file.exists(pathway_path)) {
        path_df <- data.table::fread(pathway_path)
        if (all(c("topic", "pathway") %in% names(path_df))) {
          pathway_counts <- .count_shared_topics(path_df, "pathway", "topic")
        }
      }

      backend_lbl <- if (identical(row$backend, "warplda")) "LDA" else "MultiVI"
      gene_lbl <- if (identical(row$gene_term_mode, "aggregate")) "Aggregate" else "Distinct"
      tf_lbl <- if (isTRUE(row$include_tf_terms)) "TF_on" else "TF_off"
      count_lbl <- switch(
        row$count_input,
        weight = "Weight",
        pseudo_count_bin = "Bin",
        pseudo_count_log = "Log",
        row$count_input
      )
      label <- paste0(
        "Shared-topic counts (", method, "): ",
        backend_lbl, " ", gene_lbl, " ", tf_lbl, " ", count_lbl
      )

      row_payload[[i]] <- list(
        label = label,
        links = if (is.null(link_counts)) data.frame(n_topics = integer(0), n_items = integer(0)) else link_counts,
        genes = if (is.null(gene_counts)) data.frame(n_topics = integer(0), n_items = integer(0)) else gene_counts,
        tfs = if (is.null(tf_counts)) data.frame(n_topics = integer(0), n_items = integer(0)) else tf_counts,
        paths = if (is.null(pathway_counts)) data.frame(n_topics = integer(0), n_items = integer(0)) else pathway_counts
      )
    }

    row_payload <- Filter(Negate(is.null), row_payload)
    if (!length(row_payload)) {
      if (isTRUE(verbose)) .log_warn("[topic_benchmark] No rows to plot for method={method}.")
      next
    }

    .max_y <- function(lst, key) {
      vals <- vapply(lst, function(x) {
        df <- x[[key]]
        if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(0)
        mx <- suppressWarnings(max(df$n_items, na.rm = TRUE))
        if (!is.finite(mx)) 0 else as.numeric(mx)
      }, numeric(1))
      max(1, max(vals, na.rm = TRUE))
    }
    y_links <- .max_y(row_payload, "links")
    y_genes <- .max_y(row_payload, "genes")
    y_tfs <- .max_y(row_payload, "tfs")
    y_paths <- .max_y(row_payload, "paths")

    .topic_levels <- function(lst, key) {
      vals <- unique(unlist(lapply(lst, function(x) {
        df <- x[[key]]
        if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(integer(0))
        as.integer(df$n_topics)
      })))
      vals <- vals[is.finite(vals) & vals >= 1L]
      if (!length(vals)) vals <- 1L
      as.character(sort(vals))
    }
    x_links <- .topic_levels(row_payload, "links")
    x_genes <- .topic_levels(row_payload, "genes")
    x_tfs <- .topic_levels(row_payload, "tfs")
    x_paths <- .topic_levels(row_payload, "paths")

    .plot_bar <- function(df, title, y_max, x_levels) {
      if (!is.data.frame(df) || !nrow(df)) {
        df <- data.frame(n_topics = integer(0), n_items = numeric(0))
      } else {
        df <- df[is.finite(df$n_topics) & !is.na(df$n_topics) & df$n_topics >= 1L, , drop = FALSE]
      }
      df$bar_group <- ifelse(as.integer(df$n_topics) == 1L, "n1", "n_other")
      df$n_topics <- factor(as.character(df$n_topics), levels = x_levels)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = n_topics, y = n_items)) +
        ggplot2::geom_col(ggplot2::aes(fill = bar_group), width = 0.8, na.rm = TRUE) +
        ggplot2::scale_fill_manual(values = c(n1 = "#e34a33", n_other = "#3182bd"), guide = "none") +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::coord_cartesian(ylim = c(0, y_max), expand = TRUE) +
        ggplot2::labs(
          title = title,
          x = "Topics shared (N)",
          y = "Count"
        ) +
        ggplot2::theme_classic(base_size = 9) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(face = "bold"),
          axis.text = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          aspect.ratio = 1
        )
      if (!nrow(df)) {
        p <- p + ggplot2::annotate("text", x = x_levels[[1]], y = y_max * 0.5, label = "No data", fontface = "bold", size = 3)
      }
      p
    }

    if (!requireNamespace("patchwork", quietly = TRUE)) {
      .log_abort("plot_shared_topic_counts_all_methods requires package {.pkg patchwork}.")
    }

    row_plots <- lapply(row_payload, function(x) {
      p_links <- .plot_bar(x$links, paste0(x$label, "\nLinks"), y_links, x_links)
      p_genes <- .plot_bar(x$genes, "Genes", y_genes, x_genes)
      p_tfs <- .plot_bar(x$tfs, "TFs", y_tfs, x_tfs)
      p_path <- .plot_bar(x$paths, "Pathways", y_paths, x_paths)
      (p_links | p_genes | p_tfs | p_path)
    })

    panel_w <- 4.0
    page_w <- panel_w * 4
    page_h <- max(4.5, length(row_plots) * panel_w + 0.8)
    full_plot <- patchwork::wrap_plots(row_plots, ncol = 1)
    ggplot2::ggsave(method_out_files[[method]], full_plot, width = page_w, height = page_h, limitsize = FALSE)

    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote {method} shared-topic benchmark: {method_out_files[[method]]}")
  }

  invisible(NULL)
}

plot_pathway_logp_hist_all_methods <- function(combo_grid,
                                               step3_root_dir,
                                               final_topics_subdir = "final_topics",
                                               k,
                                               vae_variant = "multivi_encoder",
                                               doc_mode = c("tf_cluster", "tf"),
                                               out_file,
                                               methods = c("peak_and_gene", "link_score_prob"),
                                               facets_ncol = 6,
                                               min_bins = 80L,
                                               max_bins = 260L,
                                               verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"

  methods <- intersect(methods, c("peak_and_gene", "link_score_prob"))
  if (!length(methods)) {
    .log_abort("methods must include peak_and_gene and/or link_score_prob.")
  }
  if (missing(out_file) || is.null(out_file) || !nzchar(out_file)) {
    .log_abort("out_file must be provided for pathway logp histogram outputs.")
  }

  out_base <- tools::file_path_sans_ext(out_file)
  out_ext <- tools::file_ext(out_file)
  if (!nzchar(out_ext)) out_ext <- "pdf"
  method_out_files <- setNames(
    vapply(methods, function(m) sprintf("%s_%s.%s", out_base, m, out_ext), character(1)),
    methods
  )

  combo_grid$combo_tag <- paste0(
    "gene_", combo_grid$gene_term_mode,
    "_tf_", ifelse(combo_grid$include_tf_terms, "on", "off"),
    "_count_", combo_grid$count_input
  )
  combo_grid$combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_grid$combo_tag)

  combo_grid$backend <- factor(combo_grid$backend, levels = c("warplda", "vae"))
  combo_grid$gene_term_mode <- factor(combo_grid$gene_term_mode, levels = c("aggregate", "unique"))
  combo_grid$include_tf_terms <- factor(combo_grid$include_tf_terms, levels = c(FALSE, TRUE))
  combo_grid$count_input <- factor(combo_grid$count_input, levels = c("weight", "pseudo_count_bin", "pseudo_count_log"))
  combo_grid <- combo_grid[order(combo_grid$backend, combo_grid$gene_term_mode, combo_grid$include_tf_terms, combo_grid$count_input), , drop = FALSE]

  combo_grid$backend <- as.character(combo_grid$backend)
  combo_grid$gene_term_mode <- as.character(combo_grid$gene_term_mode)
  combo_grid$include_tf_terms <- as.logical(as.character(combo_grid$include_tf_terms))
  combo_grid$count_input <- as.character(combo_grid$count_input)

  .row_label <- function(row) {
    backend_lbl <- if (identical(row$backend, "warplda")) "LDA" else "MultiVI"
    gene_lbl <- if (identical(row$gene_term_mode, "aggregate")) "Aggregate" else "Distinct"
    tf_lbl <- if (isTRUE(row$include_tf_terms)) "TF_on" else "TF_off"
    count_lbl <- switch(
      row$count_input,
      weight = "Weight",
      pseudo_count_bin = "Bin",
      pseudo_count_log = "Log",
      row$count_input
    )
    paste(backend_lbl, gene_lbl, tf_lbl, count_lbl)
  }

  # Gather all logp values first so bin width and x-range are globally shared.
  all_data <- vector("list", 0L)
  all_values <- numeric(0)
  panel_levels <- setNames(character(nrow(combo_grid)), combo_grid$combo_tag)
  for (i in seq_len(nrow(combo_grid))) {
    panel_levels[[i]] <- .row_label(combo_grid[i, , drop = FALSE])
  }
  panel_levels <- unname(panel_levels)

  for (method in methods) {
    for (i in seq_len(nrow(combo_grid))) {
      row <- combo_grid[i, , drop = FALSE]
      final_dir <- file.path(step3_root_dir, final_topics_subdir, row$combo_tag)
      out_dirs <- list.dirs(final_dir, recursive = FALSE, full.names = TRUE)
      if (row$backend == "vae") {
        patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", as.integer(k), "$")
      } else {
        patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", as.integer(k), "$")
      }
      out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
      if (!length(out_dirs)) next
      final_subdir <- out_dirs[[1]]
      path_csv <- file.path(final_subdir, paste0("topic_pathway_enrichment_", method, "_dotplot.csv"))
      if (!file.exists(path_csv)) next
      dt <- data.table::fread(path_csv)
      if (!("logp" %in% names(dt))) next
      vals <- suppressWarnings(as.numeric(dt$logp))
      vals <- vals[is.finite(vals) & vals >= 0]
      if (!length(vals)) next
      all_values <- c(all_values, vals)
      all_data[[length(all_data) + 1L]] <- data.frame(
        method = method,
        panel = .row_label(row),
        logp = vals,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(all_data) || !length(all_values)) {
    .log_warn("[topic_benchmark] No pathway logp values found for histogram benchmark.")
    return(invisible(NULL))
  }

  # Shared bin width across all methods+combos (Freedman-Diaconis fallback),
  # with capped bin count for visible bars.
  xr <- range(all_values, finite = TRUE)
  if (!all(is.finite(xr)) || diff(xr) <= 0) {
    xr <- c(0, 1)
  }
  bw <- suppressWarnings(2 * stats::IQR(all_values, na.rm = TRUE) / (length(all_values)^(1 / 3)))
  if (!is.finite(bw) || bw <= 0) bw <- diff(xr) / 30
  if (!is.finite(bw) || bw <= 0) bw <- 0.1
  n_bins <- as.integer(ceiling(diff(xr) / bw))
  if (!is.finite(n_bins) || n_bins < as.integer(min_bins)) n_bins <- as.integer(min_bins)
  if (n_bins > as.integer(max_bins)) n_bins <- as.integer(max_bins)
  bw <- diff(xr) / n_bins
  if (!is.finite(bw) || bw <= 0) bw <- 0.1
  # Shared y-axis max by counting with global breaks.
  y_max <- 1L
  for (d in all_data) {
    breaks <- seq(xr[1], xr[2] + bw, by = bw)
    h <- hist(d$logp, breaks = breaks, plot = FALSE)
    y_max <- max(y_max, max(h$counts, na.rm = TRUE))
  }

  .fit_gamma <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x) & x > 0]
    if (length(x) < 10L) {
      return(list(ok = FALSE, shape = NA_real_, rate = NA_real_, scale = NA_real_))
    }
    if (!requireNamespace("MASS", quietly = TRUE)) {
      return(list(ok = FALSE, shape = NA_real_, rate = NA_real_, scale = NA_real_))
    }
    fit <- tryCatch(
      suppressWarnings(MASS::fitdistr(x, densfun = "gamma")),
      error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$estimate)) {
      return(list(ok = FALSE, shape = NA_real_, rate = NA_real_, scale = NA_real_))
    }
    shape <- unname(as.numeric(fit$estimate[["shape"]]))
    rate <- unname(as.numeric(fit$estimate[["rate"]]))
    scale <- if (is.finite(rate) && rate > 0) 1 / rate else NA_real_
    list(ok = TRUE, shape = shape, rate = rate, scale = scale)
  }

  # Quantitative panel-level summaries for cross-panel comparison.
  panel_metrics <- lapply(all_data, function(d) {
    vals <- as.numeric(d$logp)
    vals <- vals[is.finite(vals) & vals >= 0]
    if (!length(vals)) return(NULL)
    gf <- .fit_gamma(vals)
    data.frame(
      method = as.character(d$method[[1]]),
      panel = as.character(d$panel[[1]]),
      n_pathways = length(vals),
      logp_mean = mean(vals),
      logp_median = stats::median(vals),
      logp_iqr = stats::IQR(vals),
      logp_q90 = stats::quantile(vals, 0.90, names = FALSE, type = 7),
      logp_q95 = stats::quantile(vals, 0.95, names = FALSE, type = 7),
      frac_pathways_logp_ge_2 = mean(vals >= 2),
      frac_pathways_logp_ge_3 = mean(vals >= 3),
      frac_pathways_logp_ge_5 = mean(vals >= 5),
      gamma_fit_shape = gf$shape,
      gamma_fit_rate = gf$rate,
      gamma_fit_scale = gf$scale,
      gamma_fit_mean = if (isTRUE(gf$ok)) gf$shape * gf$scale else NA_real_,
      gamma_fit_variance = if (isTRUE(gf$ok)) gf$shape * (gf$scale^2) else NA_real_,
      gamma_fit_tail_prob_logp_ge_3 = if (isTRUE(gf$ok)) 1 - stats::pgamma(3, shape = gf$shape, rate = gf$rate) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  panel_metrics <- Filter(Negate(is.null), panel_metrics)
  if (length(panel_metrics)) {
    metrics_dt <- data.table::rbindlist(lapply(panel_metrics, data.table::as.data.table), use.names = TRUE, fill = TRUE)
    metrics_dt[, method := factor(method, levels = c("peak_and_gene", "link_score_prob"))]
    metrics_dt <- metrics_dt[order(method, panel)]
    metrics_file <- sprintf("%s_panel_metrics.csv", out_base)
    readr::write_csv(as.data.frame(metrics_dt), metrics_file)
    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote pathway panel metrics: {metrics_file}")

    # Make a compact best-practice summary figure to compare methods/panels.
    metrics_plot <- as.data.frame(metrics_dt, stringsAsFactors = FALSE)
    metrics_plot$method <- as.character(metrics_plot$method)
    metrics_plot$panel_method <- paste(metrics_plot$method, metrics_plot$panel, sep = " | ")
    metrics_key <- c(
      "frac_pathways_logp_ge_5",
      "frac_pathways_logp_ge_3",
      "logp_q95",
      "logp_median",
      "gamma_fit_tail_prob_logp_ge_3"
    )
    metrics_key <- metrics_key[metrics_key %in% names(metrics_plot)]
    if (length(metrics_key)) {
      for (mcol in metrics_key) {
        vals <- as.numeric(metrics_plot[[mcol]])
        s <- stats::sd(vals, na.rm = TRUE)
        mu <- mean(vals, na.rm = TRUE)
        z <- if (is.finite(s) && s > 0) (vals - mu) / s else rep(0, length(vals))
        metrics_plot[[paste0("z_", mcol)]] <- z
      }
      z_cols <- paste0("z_", metrics_key)
      metrics_plot$composite_right_shift_score <- rowMeans(metrics_plot[, z_cols, drop = FALSE], na.rm = TRUE)
      metrics_plot <- metrics_plot[order(-metrics_plot$composite_right_shift_score, metrics_plot$panel_method), , drop = FALSE]
      metrics_plot$panel_method <- factor(metrics_plot$panel_method, levels = rev(unique(metrics_plot$panel_method)))

      long_dt <- data.table::melt(
        data.table::as.data.table(metrics_plot),
        id.vars = c("method", "panel", "panel_method", "composite_right_shift_score"),
        measure.vars = z_cols,
        variable.name = "metric",
        value.name = "z_score"
      )
      long_dt[, metric := gsub("^z_", "", metric)]
      metric_labels <- c(
        "frac_pathways_logp_ge_5" = "Frac logp >= 5",
        "frac_pathways_logp_ge_3" = "Frac logp >= 3",
        "logp_q95" = "Q95(logp)",
        "logp_median" = "Median(logp)",
        "gamma_fit_tail_prob_logp_ge_3" = "Gamma tail P(logp >= 3)"
      )
      metric_levels <- names(metric_labels)[names(metric_labels) %in% unique(long_dt$metric)]
      long_dt[, metric := factor(metric, levels = metric_levels, labels = unname(metric_labels[metric_levels]))]

      metrics_pdf <- sprintf("%s_panel_metrics.pdf", out_base)
      p_metrics <- ggplot2::ggplot(long_dt, ggplot2::aes(x = metric, y = panel_method, fill = z_score)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
          low = "#2166ac",
          mid = "white",
          high = "#b2182b",
          midpoint = 0,
          name = "Z-score"
        ) +
        ggplot2::labs(
          title = "Pathway logp panel metrics (higher is better)",
          subtitle = "Rows ranked by composite right-shift score (top = best)",
          x = "Metric",
          y = "Method | Panel"
        ) +
        ggplot2::theme_classic(base_size = 9) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(face = "bold", angle = 30, hjust = 1),
          axis.text.y = ggplot2::element_text(face = "bold", size = 7),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(face = "bold", hjust = 0.5),
          legend.title = ggplot2::element_text(face = "bold"),
          legend.text = ggplot2::element_text(face = "bold")
        )
      fig_h <- max(8, 0.25 * length(unique(long_dt$panel_method)))
      ggplot2::ggsave(metrics_pdf, p_metrics, width = 11, height = fig_h, limitsize = FALSE)
      if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote pathway panel metrics figure: {metrics_pdf}")
    }
  }

  for (method in methods) {
    method_frames <- Filter(function(x) identical(x$method[[1]], method), all_data)
    if (!length(method_frames)) {
      if (isTRUE(verbose)) .log_warn("[topic_benchmark] No data for method={method}.")
      next
    }
    dtm <- data.table::rbindlist(lapply(method_frames, data.table::as.data.table), use.names = TRUE, fill = TRUE)
    dtm <- dtm[is.finite(logp)]
    if (!nrow(dtm)) next
    dtm[, panel := factor(panel, levels = panel_levels)]
    # Build per-panel gamma curve for optional overlay.
    split_logp <- split(dtm$logp, dtm$panel)
    gamma_curve_dt <- lapply(names(split_logp), function(panel_name) {
      v <- split_logp[[panel_name]]
      v <- as.numeric(v)
      v <- v[is.finite(v) & v > 0]
      if (!length(v)) return(NULL)
      gf <- .fit_gamma(v)
      if (!isTRUE(gf$ok)) return(NULL)
      xs <- seq(xr[1], xr[2], length.out = 300L)
      dens <- stats::dgamma(xs, shape = gf$shape, rate = gf$rate)
      data.frame(panel = panel_name, x = xs, density = dens, stringsAsFactors = FALSE)
    })
    gamma_curve_dt <- Filter(Negate(is.null), gamma_curve_dt)
    if (length(gamma_curve_dt)) {
      gamma_curve_dt <- data.table::rbindlist(lapply(gamma_curve_dt, data.table::as.data.table), use.names = TRUE, fill = TRUE)
      gamma_curve_dt[, panel := factor(panel, levels = panel_levels)]
    }

    p <- ggplot2::ggplot(dtm, ggplot2::aes(x = logp)) +
      ggplot2::geom_histogram(
        binwidth = bw,
        boundary = xr[1],
        closed = "right",
        fill = "#2c7fb8",
        color = "grey95",
        linewidth = 0.1
      ) +
      ggplot2::facet_wrap(~panel, ncol = as.integer(facets_ncol), drop = FALSE) +
      ggplot2::coord_cartesian(xlim = xr, ylim = c(0, y_max), expand = FALSE) +
      ggplot2::labs(
        title = paste0("Pathway logp histogram (", method, ")"),
        x = "-log10(adjusted p-value)",
        y = "Count"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold", size = 8),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )

    n_panels <- length(panel_levels)
    nrow_facets <- ceiling(n_panels / as.integer(facets_ncol))
    width <- max(12, as.integer(facets_ncol) * 2.4)
    height <- max(8, nrow_facets * 2.4)
    ggplot2::ggsave(method_out_files[[method]], p, width = width, height = height, limitsize = FALSE)
    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote pathway logp histogram: {method_out_files[[method]]}")

    # Save overlay with gamma fit curve for quantitative visual comparison.
    density_aes <- if (utils::packageVersion("ggplot2") >= "3.4.0") {
      ggplot2::aes(y = ggplot2::after_stat(density))
    } else {
      ggplot2::aes(y = ..density..)
    }
    p_overlay <- ggplot2::ggplot(dtm, ggplot2::aes(x = logp)) +
      ggplot2::geom_histogram(
        density_aes,
        binwidth = bw,
        boundary = xr[1],
        closed = "right",
        fill = "grey90",
        color = "grey70",
        linewidth = 0.1
      ) +
      ggplot2::facet_wrap(~panel, ncol = as.integer(facets_ncol), drop = FALSE) +
      ggplot2::coord_cartesian(xlim = xr, expand = FALSE) +
      ggplot2::labs(
        title = paste0("Pathway logp with gamma fit (", method, ")"),
        x = "-log10(adjusted p-value)",
        y = "Density"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold", size = 8),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
    if (length(gamma_curve_dt)) {
      p_overlay <- p_overlay +
        ggplot2::geom_line(
          data = gamma_curve_dt,
          ggplot2::aes(x = x, y = density),
          inherit.aes = FALSE,
          color = "#d7301f",
          linewidth = 0.5
        )
    }
    overlay_file <- sprintf("%s_gammafit_%s.%s", out_base, method, out_ext)
    ggplot2::ggsave(overlay_file, p_overlay, width = width, height = height, limitsize = FALSE)
    if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote pathway gammafit overlay: {overlay_file}")
  }

  invisible(NULL)
}

plot_pass_state_counts_all_methods <- function(combo_grid,
                                               step3_root_dir,
                                               final_topics_subdir = "final_topics",
                                               k,
                                               vae_variant = "multivi_encoder",
                                               doc_mode = c("tf_cluster", "tf"),
                                               out_file,
                                               facets_ncol = 6,
                                               verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"

  if (missing(out_file) || is.null(out_file) || !nzchar(out_file)) {
    .log_abort("out_file must be provided for pass-state benchmark outputs.")
  }

  combo_grid$combo_tag <- paste0(
    "gene_", combo_grid$gene_term_mode,
    "_tf_", ifelse(combo_grid$include_tf_terms, "on", "off"),
    "_count_", combo_grid$count_input
  )
  combo_grid$combo_tag <- gsub("[^A-Za-z0-9_.-]+", "_", combo_grid$combo_tag)

  combo_grid$backend <- factor(combo_grid$backend, levels = c("warplda", "vae"))
  combo_grid$gene_term_mode <- factor(combo_grid$gene_term_mode, levels = c("aggregate", "unique"))
  combo_grid$include_tf_terms <- factor(combo_grid$include_tf_terms, levels = c(FALSE, TRUE))
  combo_grid$count_input <- factor(combo_grid$count_input, levels = c("weight", "pseudo_count_bin", "pseudo_count_log"))
  combo_grid <- combo_grid[order(combo_grid$backend, combo_grid$gene_term_mode, combo_grid$include_tf_terms, combo_grid$count_input), , drop = FALSE]

  combo_grid$backend <- as.character(combo_grid$backend)
  combo_grid$gene_term_mode <- as.character(combo_grid$gene_term_mode)
  combo_grid$include_tf_terms <- as.logical(as.character(combo_grid$include_tf_terms))
  combo_grid$count_input <- as.character(combo_grid$count_input)

  .row_label <- function(row) {
    backend_lbl <- if (identical(row$backend, "warplda")) "LDA" else "MultiVI"
    gene_lbl <- if (identical(row$gene_term_mode, "aggregate")) "Aggregate" else "Distinct"
    tf_lbl <- if (isTRUE(row$include_tf_terms)) "TF_on" else "TF_off"
    count_lbl <- switch(
      row$count_input,
      weight = "Weight",
      pseudo_count_bin = "Bin",
      pseudo_count_log = "Log",
      row$count_input
    )
    paste(backend_lbl, gene_lbl, tf_lbl, count_lbl)
  }

  method_levels <- c("link_score_prob", "peak_and_gene")
  method_labels <- c(
    "link_score_prob" = "Link score prob",
    "peak_and_gene" = "Peak and gene"
  )
  status_levels <- c("Pass", "Fail")
  status_colors <- c(
    "Pass" = "#2c7fb8",
    "Fail" = "#bdbdbd"
  )

  panel_levels <- vapply(seq_len(nrow(combo_grid)), function(i) .row_label(combo_grid[i, , drop = FALSE]), character(1))
  all_rows <- vector("list", nrow(combo_grid))

  for (i in seq_len(nrow(combo_grid))) {
    row <- combo_grid[i, , drop = FALSE]
    panel <- .row_label(row)
    final_dir <- file.path(step3_root_dir, final_topics_subdir, row$combo_tag)
    out_dirs <- list.dirs(final_dir, recursive = FALSE, full.names = TRUE)
    if (row$backend == "vae") {
      patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_K", as.integer(k), "$")
    } else {
      patt <- paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_K", as.integer(k), "$")
    }
    out_dirs <- out_dirs[grepl(patt, basename(out_dirs))]
    counts <- data.frame(
      method = rep(method_levels, each = length(status_levels)),
      status = rep(status_levels, times = length(method_levels)),
      count = 0L,
      stringsAsFactors = FALSE
    )
    if (length(out_dirs)) {
      link_path <- file.path(out_dirs[[1]], "topic_links.csv")
      if (file.exists(link_path)) {
        dt <- data.table::fread(link_path)
        req <- c("peak_pass", "gene_pass", "peak_id", "gene_key")
        if (all(req %in% names(dt))) {
          if (!("tf" %in% names(dt)) && isTRUE(verbose)) {
            .log_warn("[topic_benchmark] Missing tf in topic_links.csv for {row$combo_tag}; using peak_id+gene_key link key.")
          }
          if ("tf" %in% names(dt)) {
            dt <- dt[!is.na(tf) & tf != "" & !is.na(peak_id) & peak_id != "" & !is.na(gene_key) & gene_key != ""]
            dt[, link_id := paste(tf, peak_id, gene_key, sep = "::")]
          } else {
            dt <- dt[!is.na(peak_id) & peak_id != "" & !is.na(gene_key) & gene_key != ""]
            dt[, link_id := paste(peak_id, gene_key, sep = "::")]
          }
          if (nrow(dt)) {
            # Strictly count each tf-peak-gene link once in this plot.
            dt <- unique(dt, by = c("link_id", "peak_pass", "gene_pass"))
            if (!("link_pass" %in% names(dt))) {
              if (isTRUE(verbose)) .log_warn("[topic_benchmark] Missing link_pass in topic_links.csv for {row$combo_tag}")
              dt[, link_pass := FALSE]
            }
            link_status <- dt[, .(
              pass_link_score_prob = any(.as_logical_flag(link_pass), na.rm = TRUE),
              pass_peak_and_gene = any(.as_logical_flag(peak_pass) & .as_logical_flag(gene_pass), na.rm = TRUE)
            ), by = .(link_id)]

            n_total <- nrow(link_status)
            n_pass_prob <- sum(link_status$pass_link_score_prob, na.rm = TRUE)
            n_pass_and <- sum(link_status$pass_peak_and_gene, na.rm = TRUE)

            counts$count[counts$method == "link_score_prob" & counts$status == "Pass"] <- as.integer(n_pass_prob)
            counts$count[counts$method == "link_score_prob" & counts$status == "Fail"] <- as.integer(n_total - n_pass_prob)
            counts$count[counts$method == "peak_and_gene" & counts$status == "Pass"] <- as.integer(n_pass_and)
            counts$count[counts$method == "peak_and_gene" & counts$status == "Fail"] <- as.integer(n_total - n_pass_and)
          }
        } else if (isTRUE(verbose)) {
          .log_warn("[topic_benchmark] Missing required columns in topic_links.csv for {row$combo_tag}")
        }
      } else if (isTRUE(verbose)) {
        .log_warn("[topic_benchmark] Missing topic_links.csv for {row$combo_tag}")
      }
    } else if (isTRUE(verbose)) {
      .log_warn("[topic_benchmark] Missing final_topics for {row$combo_tag}")
    }

    counts$panel <- panel
    all_rows[[i]] <- counts[, c("panel", "method", "status", "count")]
  }

  plot_df <- data.table::rbindlist(lapply(all_rows, data.table::as.data.table), use.names = TRUE, fill = TRUE)
  plot_df[, panel := factor(panel, levels = panel_levels)]
  plot_df[, method := factor(method, levels = method_levels, labels = unname(method_labels[method_levels]))]
  plot_df[, status := factor(status, levels = status_levels)]
  y_max <- max(1, suppressWarnings(max(plot_df$count, na.rm = TRUE)))
  if (!is.finite(y_max)) y_max <- 1

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method, y = count, fill = status)) +
    ggplot2::geom_col(width = 0.8, na.rm = TRUE, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_fill_manual(values = status_colors, drop = FALSE) +
    ggplot2::facet_wrap(~panel, ncol = as.integer(facets_ncol), drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, y_max), expand = FALSE) +
    ggplot2::labs(
      title = "Link assignment pass/fail by method",
      x = "Method",
      y = "Count",
      fill = "Status",
      caption = "Pass: link assigned to at least one topic in any doc_id under the method.\nLink key: tf + peak_id + gene_key."
    ) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(face = "bold", hjust = 0)
    )

  n_panels <- length(panel_levels)
  nrow_facets <- ceiling(n_panels / as.integer(facets_ncol))
  width <- max(12, as.integer(facets_ncol) * 2.4)
  height <- max(8, nrow_facets * 2.4)
  ggplot2::ggsave(out_file, p, width = width, height = height, limitsize = FALSE)
  if (isTRUE(verbose)) .log_inform("[topic_benchmark] Wrote pass-state count benchmark: {out_file}")

  invisible(NULL)
}
