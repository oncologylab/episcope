test_that("public Module 3 topic-model defaults use native WarpLDA", {
  expect_identical(eval(formals(train_topic_models)$backend)[[1L]], "warplda")
  expect_identical(eval(formals(extract_regulatory_topics)$backend)[[1L]], "warplda")
})

test_that("WarpLDA training honors explicit topic count input", {
  root <- withr::local_tempdir()
  filtered_dir <- file.path(root, "filtered")
  out_dir <- file.path(root, "models")
  dir.create(filtered_dir, recursive = TRUE)

  links <- data.table::data.table(
    comparison_id = rep("CondA_vs_CondB", 6),
    comparison_label = rep("CondA vs CondB", 6),
    cond1_id = rep("CondA", 6),
    cond2_id = rep("CondB", 6),
    cond1_label = rep("CondA", 6),
    cond2_label = rep("CondB", 6),
    tf = rep(c("TF1", "TF2"), each = 3),
    gene_key = paste0("G", 1:6),
    peak_id = paste0("P", 1:6),
    log2fc_gene = c(1.2, 1.4, 1.6, 1.1, 1.3, 1.5),
    log2fc_fp = c(1.0, 1.1, 1.2, 1.0, 1.1, 1.2),
    delta_fp_score = c(0.8, 0.9, 1.0, 0.8, 0.9, 1.0),
    log2FC_fp_score = c(1.0, 1.1, 1.2, 1.0, 1.1, 1.2),
    log2FC_gene_expr = c(1.2, 1.4, 1.6, 1.1, 1.3, 1.5),
    tf_expr_cond1 = 20,
    tf_expr_cond2 = 20,
    gene_expr_cond1 = 20,
    gene_expr_cond2 = 20,
    fp_score_cond1 = 5,
    fp_score_cond2 = 5
  )
  up_path <- file.path(filtered_dir, "CondA_vs_CondB_filtered_links_up.csv")
  data.table::fwrite(links, up_path)
  data.table::fwrite(
    data.table::data.table(
      comparison_id = "CondA_vs_CondB",
      comparison_label = "CondA vs CondB",
      up_path = up_path,
      down_path = NA_character_
    ),
    file.path(filtered_dir, "filtered_links_manifest.csv")
  )

  train_topic_models(
    Kgrid = 2L,
    input_dir = filtered_dir,
    output_dir = out_dir,
    analysis_label = "count_input_test",
    doc_mode = "tf",
    doc_design = "comparison",
    fp_term_mode = "unique",
    gene_term_mode = "unique",
    include_tf_terms = TRUE,
    min_df = 1,
    count_method = "bin",
    count_input = "pseudo_count_log",
    backend = "warplda",
    warplda_iterations = 1L,
    local_threads = 1L,
    reuse_if_exists = FALSE,
    abs_log2fc_fp_min = 0,
    abs_log2fc_gene_min = 0,
    abs_delta_fp_min = 0,
    direction_consistency = "none",
    require_fp_bound_either = FALSE,
    require_tf_expr_either = FALSE,
    require_gene_expr_either = FALSE,
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    threshold_tf_expr = 0,
    save_full_doc_term_csv = FALSE,
    flat_output = TRUE
  )

  doc_term <- readRDS(file.path(out_dir, "rds", "doc_term.rds"))
  dtm <- readRDS(file.path(out_dir, "rds", "dtm.rds"))
  metrics <- data.table::fread(file.path(out_dir, "model_metrics.csv"))
  summary <- data.table::fread(file.path(out_dir, "topic_input_summary.csv"))

  expect_equal(sum(dtm@x), sum(doc_term$pseudo_count_log))
  expect_false(isTRUE(all.equal(sum(dtm@x), sum(doc_term$pseudo_count_bin))))
  expect_equal(unique(summary$count_input_effective), "pseudo_count_log")
  expect_equal(unique(metrics$count_input_effective), "pseudo_count_log")
})

test_that("Module 3 fresh training clears stale model artifacts", {
  out_dir <- withr::local_tempdir()
  stale_paths <- c(
    file.path(out_dir, "model_metrics.csv"),
    file.path(out_dir, "rds", "model_metrics.rds"),
    file.path(out_dir, "rds", "theta_maximum_K2.rds"),
    file.path(out_dir, "tmp_models", "fit_K2.rds"),
    file.path(out_dir, "vae_models", "theta_K2.csv"),
    file.path(out_dir, "vae_progress.tsv")
  )
  for (dir in unique(dirname(stale_paths))) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  for (path in stale_paths) {
    writeLines("stale", path)
  }

  expect_true(any(file.exists(stale_paths)))
  craftgrn:::.reset_topic_model_artifacts(out_dir, backend = "warplda", reuse_if_exists = FALSE)
  expect_false(file.exists(file.path(out_dir, "model_metrics.csv")))
  expect_false(file.exists(file.path(out_dir, "rds", "model_metrics.rds")))
  expect_false(file.exists(file.path(out_dir, "rds", "theta_maximum_K2.rds")))
  expect_false(dir.exists(file.path(out_dir, "tmp_models")))
  expect_false(dir.exists(file.path(out_dir, "vae_models")))
  expect_false(file.exists(file.path(out_dir, "vae_progress.tsv")))

  dir.create(file.path(out_dir, "tmp_models"), recursive = TRUE, showWarnings = FALSE)
  keep_path <- file.path(out_dir, "tmp_models", "fit_K2.rds")
  writeLines("keep", keep_path)
  craftgrn:::.reset_topic_model_artifacts(out_dir, backend = "warplda", reuse_if_exists = TRUE)
  expect_true(file.exists(keep_path))
})

test_that("condition doc_tf applies condition thresholds and TF self terms", {
  edges <- data.table::data.table(
    comparison_id = c("C1", "C1"),
    case_id = c("CondA", "CondA"),
    ctrl_id = c("CondB", "CondB"),
    tf = c("TF1", "TF1"),
    gene_key = c("G1", "G2"),
    peak_id = c("P1", "P2"),
    fp_score_case = c(3, 1),
    fp_score_ctrl = c(3, 3),
    gene_expr_case = c(12, 9),
    gene_expr_ctrl = c(12, 12),
    tf_expr_case = c(11, 11),
    tf_expr_ctrl = c(9, 9)
  )

  docs <- add_condition_tf_docs(edges, doc_mode = "tf")
  terms <- build_doc_term_condition_union(
    docs,
    count_method = "log",
    threshold_gene_expr = 10,
    threshold_fp_score = 2,
    threshold_tf_expr = 10,
    include_tf_terms = TRUE,
    require_tf_expr = TRUE
  )

  expect_true("CondA::TF1" %in% terms$doc_id)
  expect_false("CondB::TF1" %in% terms$doc_id)

  cond_terms <- terms[doc_id == "CondA::TF1", term_id]
  expect_true("GENE:G1" %in% cond_terms)
  expect_true("GENE:TF1" %in% cond_terms)
  expect_true("PEAK:P1" %in% cond_terms)
  expect_false("GENE:G2" %in% cond_terms)
  expect_false("PEAK:P2" %in% cond_terms)
})

test_that("WarpLDA one-option runner exposes pathway dotplot control", {
  fmls <- names(formals(run_tfdocs_warplda_one_option))
  expect_true("pathway_make_heatmap" %in% fmls)
  expect_true("pathway_make_dotplot" %in% fmls)
})

test_that("condition document terms use one condition-specific value per term", {
  edges <- data.table::data.table(
    condition_label = rep("CondA", 4),
    tf_doc = rep("TF1", 4),
    tf = rep("TF1", 4),
    gene_key = c("G1", "G1", "G2", "G2"),
    peak_id = c("P1", "P1", "P2", "P2"),
    fp_score_condition = c(10, 10, 20, 20),
    gene_expr_condition = c(10, 10, 10, 10),
    tf_expr_condition = c(10, 10, 10, 10)
  )

  terms <- build_doc_term_condition_union(
    edges,
    count_method = "log",
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    threshold_tf_expr = 0,
    include_tf_terms = TRUE,
    require_tf_expr = TRUE
  )

  expect_equal(nrow(terms[term_id == "GENE:G1"]), 1L)
  expect_equal(nrow(terms[term_id == "PEAK:P1"]), 1L)
  expect_equal(nrow(terms[term_id == "GENE:TF1"]), 1L)
  expect_equal(terms[term_id == "GENE:G1", weight], 10)
  expect_equal(terms[term_id == "PEAK:P1", weight], 10)
  expect_equal(terms[term_id == "PEAK:P2", weight], 20)
  expect_equal(terms[term_id == "GENE:TF1", weight], 10)
})

test_that("condition document terms warn on inconsistent repeated values", {
  gene_bad <- data.table::data.table(
    condition_label = c("CondA", "CondA"),
    tf_doc = c("TF1", "TF1"),
    tf = c("TF1", "TF1"),
    gene_key = c("G1", "G1"),
    peak_id = c("P1", "P2"),
    fp_score_condition = c(10, 10),
    gene_expr_condition = c(10, 11),
    tf_expr_condition = c(10, 10)
  )
  expect_warning(
    build_doc_term_condition_union(gene_bad, threshold_gene_expr = 0, threshold_fp_score = 0),
    "repeated gene document-term group"
  )

  peak_bad <- data.table::copy(gene_bad)
  peak_bad[, `:=`(
    gene_expr_condition = c(10, 10),
    peak_id = c("P1", "P1"),
    fp_score_condition = c(10, 11)
  )]
  expect_warning(
    build_doc_term_condition_union(peak_bad, threshold_gene_expr = 0, threshold_fp_score = 0),
    "repeated peak document-term group"
  )

  tf_bad <- data.table::copy(gene_bad)
  tf_bad[, `:=`(
    gene_expr_condition = c(10, 10),
    tf_expr_condition = c(10, 11)
  )]
  expect_warning(
    build_doc_term_condition_union(
      tf_bad,
      threshold_gene_expr = 0,
      threshold_fp_score = 0,
      threshold_tf_expr = 0,
      include_tf_terms = TRUE,
      require_tf_expr = TRUE
    ),
    "repeated TF self document-term group"
  )
})

test_that("comparison doc_tf uses direction-relevant condition thresholds", {
  edges <- data.table::data.table(
    comparison_id = c("C1", "C1", "C1"),
    tf = c("TF1", "TF1", "TF1"),
    gene_key = c("G_up", "G_down", "G_low"),
    peak_id = c("P_up", "P_down", "P_low"),
    log2fc_gene = c(1.5, -1.4, 1.2),
    log2fc_fp = c(1.2, -1.1, 1.1),
    delta_fp = c(0.8, -0.7, 0.6),
    fc_mag_fp = c(2.3, 2.1, 2.1),
    fc_mag_gene = c(2.8, 2.6, 2.4),
    fc_mag_tf = c(2.2, 2.2, 2.2),
    tf_expr_case = c(11, 8, 11),
    tf_expr_ctrl = c(8, 12, 8),
    gene_expr_case = c(12, 8, 9),
    gene_expr_ctrl = c(8, 13, 8),
    fp_score_case = c(3, 1, 3),
    fp_score_ctrl = c(1, 3, 1)
  )

  docs <- add_tf_docs(edges, doc_mode = "tf", direction_by = "gene")
  terms <- build_doc_term_joint(
    docs,
    weight_type_peak = "log2fc_fp",
    weight_type_gene = "log2fc_gene",
    min_df = 1,
    count_method = "log",
    distinct_terms = TRUE,
    gene_term_mode = "aggregate",
    include_tf_terms = TRUE,
    balance_mode = "min",
    threshold_gene_expr = 10,
    threshold_fp_score = 2,
    threshold_tf_expr = 10,
    require_condition_thresholds = TRUE
  )

  up_terms <- terms[doc_id == "C1::TF1::Target-Up", term_id]
  down_terms <- terms[doc_id == "C1::TF1::Target-Down", term_id]

  expect_true("GENE:G_up" %in% up_terms)
  expect_true("PEAK:P_up" %in% up_terms)
  expect_true("GENE:TF1" %in% up_terms)
  expect_false("GENE:G_low" %in% up_terms)
  expect_false("PEAK:P_low" %in% up_terms)

  expect_true("GENE:G_down" %in% down_terms)
  expect_true("PEAK:P_down" %in% down_terms)
  expect_true("GENE:TF1" %in% down_terms)
})

test_that("comparison doc_tf drops zero-direction rows before threshold selection", {
  edges <- data.table::data.table(
    comparison_id = c("C1", "C1"),
    tf = c("TF1", "TF1"),
    gene_key = c("G0", "G1"),
    peak_id = c("P0", "P1"),
    log2fc_gene = c(0, 1),
    tf_expr_case = c(10, 10),
    tf_expr_ctrl = c(10, 10),
    gene_expr_case = c(10, 10),
    gene_expr_ctrl = c(10, 10),
    fp_score_case = c(2, 2),
    fp_score_ctrl = c(2, 2)
  )

  docs <- add_tf_docs(edges, doc_mode = "tf", direction_by = "gene")

  expect_equal(nrow(docs), 1L)
  expect_equal(docs$doc_id, "C1::TF1::Target-Up")
  expect_equal(docs$gene_expr_condition, 10)
})

test_that("comparison topic documents preserve readable comparison labels", {
  edges <- data.table::data.table(
    comparison_id = c("Fib_BIRT_vs_BI", "Fib_BIRT_vs_BI"),
    comparison_label = c("Fibroblast BIRT vs BI", "Fibroblast BIRT vs BI"),
    tf = c("TF1", "TF2"),
    gene_key = c("G1", "G2"),
    peak_id = c("P1", "P2"),
    log2fc_gene = c(1, -1),
    tf_expr_cond1 = c(10, 10),
    tf_expr_cond2 = c(10, 10),
    gene_expr_cond1 = c(10, 10),
    gene_expr_cond2 = c(10, 10),
    fp_score_cond1 = c(2, 2),
    fp_score_cond2 = c(2, 2)
  )

  docs <- add_tf_docs(edges, doc_mode = "tf", direction_by = "gene")

  expect_equal(
    docs$doc_id,
    c("Fib_BIRT_vs_BI::TF1::Target-Up", "Fib_BIRT_vs_BI::TF2::Target-Down")
  )
  expect_equal(
    docs$doc_display_label,
    c("Fibroblast BIRT vs BI::Target-Up", "Fibroblast BIRT vs BI::Target-Down")
  )
  expect_equal(unique(docs$comparison_label), "Fibroblast BIRT vs BI")
})

test_that("topic comparison heatmap labels prefer comparison metadata labels", {
  dt <- data.table::data.table(
    doc_id = c("Fib_BIRT_vs_BI::TF1::Target-Up", "Fib_BIRT_vs_BI::TF2::Target-Down"),
    comparison_id = c("Fib_BIRT_vs_BI", "Fib_BIRT_vs_BI"),
    tf_doc = c("TF1", "TF2"),
    direction = c("Target-Up", "Target-Down"),
    direction_label = c("Target-Up", "Target-Down"),
    Topic1 = c(0.8, 0.2),
    Topic2 = c(0.2, 0.8)
  )
  edges_docs <- data.table::data.table(
    doc_id = c("Fib_BIRT_vs_BI::TF1::Target-Up", "Fib_BIRT_vs_BI::TF2::Target-Down"),
    comparison_id = c("Fib_BIRT_vs_BI", "Fib_BIRT_vs_BI"),
    comparison_label = c("Fibroblast BIRT vs BI", "Fibroblast BIRT vs BI"),
    doc_display_label = c("Fibroblast BIRT vs BI::Target-Up", "Fibroblast BIRT vs BI::Target-Down")
  )

  out <- craftgrn:::.resolve_topic_comparison_labels(dt, edges_docs = edges_docs)

  expect_equal(
    out$comparison_label,
    c("Fibroblast BIRT vs BI::Target-Up", "Fibroblast BIRT vs BI::Target-Down")
  )
  expect_equal(out$comparison_id, c("Fib_BIRT_vs_BI", "Fib_BIRT_vs_BI"))
})

test_that("Module 3 edge loading applies readable labels from manifest", {
  tmp <- withr::local_tempdir()
  edges <- data.table::data.table(
    comparison_id = "Fib_BIRT_vs_BI",
    comparison_label = "Fib_BIRT_vs_BI",
    tf = "TF1",
    gene_key = "G1",
    peak_id = "P1"
  )
  data.table::fwrite(edges, file.path(tmp, "Fib_BIRT_vs_BI_filtered_links_up.csv"))
  data.table::fwrite(data.table::data.table(
    comparison_id = "Fib_BIRT_vs_BI",
    comparison_label = "Fibroblast BIRT vs BI",
    up_path = file.path(tmp, "Fib_BIRT_vs_BI_filtered_links_up.csv"),
    down_path = file.path(tmp, "Fib_BIRT_vs_BI_filtered_links_down.csv")
  ), file.path(tmp, "filtered_links_manifest.csv"))

  out <- craftgrn:::.apply_module3_manifest_labels(edges, tmp)

  expect_equal(out$comparison_label[[1]], "Fibroblast BIRT vs BI")
})

test_that("condition doc_ctf wrapper keeps cluster-level documents", {
  edges <- data.table::data.table(
    comparison_id = "C1",
    case_id = "CondA",
    ctrl_id = "CondB",
    tf = "TF1",
    gene_key = "G1",
    peak_id = "P1",
    fp_score_case = 3,
    fp_score_ctrl = 3,
    gene_expr_case = 12,
    gene_expr_ctrl = 12,
    tf_expr_case = 11,
    tf_expr_ctrl = 11
  )
  cluster_map <- c(TF1 = "ClusterA")

  docs <- add_condition_tf_cluster_docs(edges, tf_cluster_map = cluster_map)

  expect_true(all(c("CondA::ClusterA", "CondB::ClusterA") %in% docs$doc_id))
  expect_false(any(grepl("CondA::TF1", docs$doc_id, fixed = TRUE)))
})

test_that("condition FP term modes build unique, aggregated, and weighted terms", {
  edges <- data.table::data.table(
    condition_label = "CondA",
    tf_doc = "TF1",
    tf = "TF1",
    gene_key = c("G1", "G1"),
    peak_id = c("P1", "P2"),
    fp_score_condition = c(2, 3),
    gene_expr_condition = c(5, 5),
    tf_expr_condition = c(10, 10)
  )

  uniq <- build_doc_term_condition_union(
    edges,
    count_method = "log",
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    fp_term_mode = "unique"
  )
  expect_true(all(c("PEAK:P1", "PEAK:P2", "GENE:G1") %in% uniq$term_id))

  aggr <- build_doc_term_condition_union(
    edges,
    count_method = "log",
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    fp_term_mode = "aggregate"
  )
  expect_true(all(c("PEAK:G1", "GENE:G1") %in% aggr$term_id))
  expect_false(any(c("PEAK:P1", "PEAK:P2") %in% aggr$term_id))
  expect_equal(aggr[term_id == "GENE:G1", weight], 5)
  expect_equal(aggr[term_id == "PEAK:G1", weight], 5)

  weighted <- build_doc_term_condition_union(
    edges,
    count_method = "log",
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    fp_term_mode = "aggregate_weight"
  )
  expect_true("GENE:G1" %in% weighted$term_id)
  expect_false(any(grepl("^PEAK:", weighted$term_id)))
  expect_equal(weighted[term_id == "GENE:G1", weight], 25)
})

test_that("condition fp_aggr sanity check ignores expected TF self term only", {
  edges <- data.table::data.table(
    condition_label = "CondA",
    tf_doc = "TF1",
    tf = "TF1",
    gene_key = c("G1", "TF1", "TF2"),
    peak_id = c("P1", "P_self", "P2"),
    fp_score_condition = c(2, 4, 3),
    gene_expr_condition = c(5, 7, 6),
    tf_expr_condition = c(10, 10, 10)
  )

  expect_warning(
    terms <- build_doc_term_condition_union(
      edges,
      count_method = "log",
      threshold_gene_expr = 0,
      threshold_fp_score = 0,
      threshold_tf_expr = 0,
      include_tf_terms = TRUE,
      require_tf_expr = TRUE,
      fp_term_mode = "aggregate"
    ),
    NA
  )
  expect_true(all(c("GENE:TF1", "GENE:TF2", "PEAK:TF1", "PEAK:TF2") %in% terms$term_id))
})

test_that("VAE variant aliases normalize to Python variants", {
  expect_equal(.normalize_vae_python_variant("vmlp"), "vae_mlp")
  expect_equal(.normalize_vae_python_variant("model_vmlp"), "vae_mlp")
  expect_equal(.normalize_vae_python_variant("mve"), "multivi_encoder")
  expect_equal(.normalize_vae_python_variant("model_mve"), "multivi_encoder")
  expect_equal(.normalize_vae_python_variant("moetm"), "moetm_encoder_decoder")
  expect_equal(.normalize_vae_python_variant("model_moetm"), "moetm_encoder_decoder")
  expect_equal(.normalize_vae_python_variant("vae_mlp"), "vae_mlp")
  expect_error(.normalize_vae_python_variant("bad_variant"), "Unsupported VAE variant")
})

test_that("VAE cache plan falls back to CSV when Python pyarrow is unavailable", {
  csv_plan <- .vae_doc_term_cache_plan(has_r_arrow = TRUE, python_has_pyarrow = FALSE, save_full_doc_term_csv = FALSE)
  expect_false(csv_plan$write_arrow)
  expect_true(csv_plan$save_full_doc_term_csv)

  arrow_plan <- .vae_doc_term_cache_plan(has_r_arrow = TRUE, python_has_pyarrow = TRUE, save_full_doc_term_csv = FALSE)
  expect_true(arrow_plan$write_arrow)
  expect_false(arrow_plan$save_full_doc_term_csv)
})

test_that("WarpLDA resume requires existing fit artifacts, not metrics alone", {
  fit_dir <- tempfile("warplda-fits-")
  dir.create(fit_dir)
  fit_k2 <- file.path(fit_dir, "fit_K2.rds")
  fit_k4 <- tempfile(fileext = ".rds")
  saveRDS(list(K = 2L, metrics = list(perplexity = 1, loglik_approx = -1, n_tokens = 10)), fit_k2)
  saveRDS(list(K = 4L, metrics = list(perplexity = 1, loglik_approx = -1, n_tokens = 10)), fit_k4)

  k_grid <- c(2L, 3L, 4L)
  fit_files_all <- file.path(fit_dir, sprintf("fit_K%d.rds", k_grid))
  existing_metrics <- data.table::data.table(
    K = c(2L, 3L, 4L),
    perplexity = 1,
    loglik = -1,
    n_tokens = 10,
    fit_file = c(fit_k2, file.path(fit_dir, "fit_K3.rds"), fit_k4)
  )

  status <- .warplda_completed_from_cache(k_grid, existing_metrics, fit_files_all)
  expect_equal(status$done_from_file, c(TRUE, FALSE, FALSE))
  expect_equal(status$done_from_metrics, c(TRUE, FALSE, TRUE))
})

test_that("gammafit topic-term scope preserves existing topic-group behavior", {
  score_mat <- matrix(
    c(
      0.9, 0.1, 0.3, 0.2,
      0.8, 0.7, 0.95, 0.1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2"))
  )

  default_terms <- binarize_topics(score_mat, method = "gammafit", min_terms = 1L)
  explicit_terms <- binarize_topics(
    score_mat,
    method = "gammafit",
    min_terms = 1L,
    gammafit_scope = "topic_term_group"
  )

  expect_equal(default_terms, explicit_terms)
  expect_equal(names(default_terms), c("topic", "term_id", "score", "in_topic"))
  expect_true(default_terms[default_terms$topic == 1L & default_terms$term_id == "PEAK:P1", "in_topic"])
  expect_true(default_terms[default_terms$topic == 2L & default_terms$term_id == "PEAK:P1", "in_topic"])
  expect_true(default_terms[default_terms$topic == 1L & default_terms$term_id == "GENE:G1", "in_topic"])
  expect_true(default_terms[default_terms$topic == 2L & default_terms$term_id == "GENE:G1", "in_topic"])
})

test_that("gammafit global term-group scope pools scores across topics", {
  score_mat <- matrix(
    c(
      0.9, 0.1, 0.3, 0.2,
      0.8, 0.7, 0.95, 0.1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2"))
  )

  cutoffs <- .gammafit_cutoffs_by_termclass(
    score_mat,
    min_terms = 1L,
    gammafit_scope = "global_term_group"
  )
  expect_equal(unique(cutoffs$peaks_gamma_cutoff), 0.9)
  expect_equal(unique(cutoffs$gene_gamma_cutoff), 0.95)

  global_terms <- binarize_topics(
    score_mat,
    method = "gammafit",
    min_terms = 1L,
    gammafit_scope = "global_term_group"
  )
  expect_true(global_terms[global_terms$topic == 1L & global_terms$term_id == "PEAK:P1", "in_topic"])
  expect_false(global_terms[global_terms$topic == 2L & global_terms$term_id == "PEAK:P1", "in_topic"])
  expect_false(global_terms[global_terms$topic == 1L & global_terms$term_id == "GENE:G1", "in_topic"])
  expect_true(global_terms[global_terms$topic == 2L & global_terms$term_id == "GENE:G1", "in_topic"])
  expect_true(all(c("term_group", "gammafit_scope", "gamma_cutoff") %in% names(global_terms)))
})

test_that("topic link gammafit scope is passed through to gene-only links", {
  edges <- data.table::data.table(
    doc_id = "D1",
    tf = "TF1",
    peak_id = "P1",
    gene_key = "G1"
  )
  score_mat <- matrix(
    c(
      0.9, 0.1, 0.3,
      0.8, 0.7, 0.95
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("PEAK:P1", "PEAK:P2", "GENE:G1"))
  )

  topic_scope <- compute_topic_links(
    edges,
    score_mat,
    fp_term_mode = "aggregate_weight",
    binarize_method = "gammafit",
    gammafit_scope = "topic_term_group",
    link_method = "gammafit",
    min_terms = 1L,
    overwrite = TRUE
  )
  global_scope <- compute_topic_links(
    edges,
    score_mat,
    fp_term_mode = "aggregate_weight",
    binarize_method = "gammafit",
    gammafit_scope = "global_term_group",
    link_method = "gammafit",
    min_terms = 1L,
    overwrite = TRUE
  )

  expect_equal(sort(topic_scope$topic_num), c(1L, 2L))
  expect_equal(global_scope$topic_num, 2L)
  expect_equal(unique(global_scope$gene_gamma_cutoff), 0.95)
})

test_that("gammafit topic link summary reports pre-filter scored rows", {
  edges <- data.table::data.table(
    doc_id = c("D1", "D1"),
    tf = c("TF1", "TF1"),
    peak_id = c("P1", "P2"),
    gene_key = c("G1", "G2")
  )
  score_mat <- matrix(
    c(
      0.9, 0.1, 0.8, 0.2,
      0.2, 0.7, 0.1, 0.9
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2"))
  )
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L),
    term_id = c("PEAK:P1", "GENE:G1"),
    in_topic = TRUE
  )
  out_dir <- tempfile("module3-gammafit-summary-")
  dir.create(out_dir, recursive = TRUE)

  compute_topic_links(
    edges,
    score_mat,
    topic_terms = topic_terms,
    binarize_method = "topn",
    link_method = "gammafit",
    pass_file = file.path(out_dir, "topic_links_pass.csv"),
    output_mode = "pass",
    overwrite = TRUE
  )
  summary <- data.table::fread(file.path(out_dir, "topic_link_summary.csv"))

  expect_equal(summary$n_scored_rows, 4)
  expect_equal(summary$n_pass_rows, 1)
})

test_that("TF document topic assignment data aligns document and term panels", {
  theta <- matrix(
    c(0.8, 0.2, 0.1, 0.9, 0.6, 0.4),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("CondA::TF1", "CondA::TF2", "CondB::TF1"),
      c("Topic1", "Topic2")
    )
  )
  phi <- matrix(
    c(0.7, 0.3, 0.2, 0.8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:TF1", "GENE:TF2"))
  )
  topic_terms <- data.table::data.table(
    topic = c("Topic1", "Topic2"),
    topic_num = c(1L, 2L),
    term_id = c("GENE:TF1", "GENE:TF2"),
    in_topic = TRUE
  )

  plot_dt <- .prepare_tf_doc_topic_assignment_heatmap_data(
    theta = theta,
    phi = phi,
    topic_terms = topic_terms,
    doc_design = "condition",
    aggregate_fun = "max"
  )

  expect_true(all(c("Aggregate", "CondA", "CondB") %in% plot_dt$page_label))
  expect_true(all(c("TF doc score", "TF term score") %in% plot_dt$panel))
  expect_equal(unique(plot_dt[page_label == "Aggregate" & tf == "TF1" & topic == "Topic1" & panel == "TF doc score", value]), 1)
  expect_true(plot_dt[page_label == "CondA" & tf == "TF1" & topic == "Topic1" & panel == "TF term score", term_pass])
  expect_equal(
    unique(plot_dt[page_label == "CondA" & panel == "TF doc score", as.character(topic)]),
    unique(plot_dt[page_label == "CondA" & panel == "TF term score", as.character(topic)])
  )
})

test_that("comparison FP term modes build aggregated peak and weighted gene terms", {
  edges <- data.table::data.table(
    comparison_id = "C1",
    tf = "TF1",
    gene_key = c("G1", "G1"),
    peak_id = c("P1", "P2"),
    log2fc_gene = c(1, 1),
    log2fc_fp = c(1, 1),
    delta_fp = c(2, 3),
    fc_mag_fp = c(2, 3),
    fc_mag_gene = c(5, 5),
    fc_mag_tf = c(10, 10),
    tf_expr_case = c(10, 10),
    tf_expr_ctrl = c(10, 10),
    gene_expr_case = c(10, 10),
    gene_expr_ctrl = c(10, 10),
    fp_score_case = c(2, 3),
    fp_score_ctrl = c(2, 3)
  )
  docs <- add_tf_docs(edges, doc_mode = "tf", direction_by = "gene")

  aggr <- build_doc_term_joint(
    docs,
    weight_type_peak = "log2fc_fp",
    weight_type_gene = "log2fc_gene",
    min_df = 1,
    count_method = "log",
    fp_term_mode = "aggregate",
    balance_mode = "min"
  )
  expect_true(all(c("PEAK:G1", "GENE:G1") %in% aggr$term_id))
  expect_false(any(c("PEAK:P1", "PEAK:P2") %in% aggr$term_id))

  weighted <- build_doc_term_joint(
    docs,
    weight_type_peak = "log2fc_fp",
    weight_type_gene = "log2fc_gene",
    min_df = 1,
    count_method = "log",
    fp_term_mode = "aggregate_weight",
    balance_mode = "min"
  )
  expect_true("GENE:G1" %in% weighted$term_id)
  expect_false(any(grepl("^PEAK:", weighted$term_id)))
  expect_equal(weighted[term_id == "GENE:G1", weight], 2)
})

test_that("doc-term cache writes preview and optional full CSV", {
  tmp <- tempfile("doc_term_cache_")
  dir.create(tmp)
  dt <- data.table::data.table(
    doc_id = paste0("D", 1:120),
    term_id = paste0("GENE:G", 1:120),
    weight = 1,
    pseudo_count = 1,
    pseudo_count_bin = 1,
    pseudo_count_log = 1
  )
  paths <- write_doc_term_cache(dt, tmp, save_full_doc_term_csv = FALSE)
  expect_true(file.exists(paths$preview))
  expect_false(file.exists(file.path(tmp, "doc_term.csv")))
  expect_equal(nrow(data.table::fread(paths$preview)), 100L)

  paths_full <- write_doc_term_cache(dt, tmp, save_full_doc_term_csv = TRUE)
  expect_true(file.exists(paths_full$csv))
})

test_that("parallel topic link scoring matches serial scoring", {
  skip_on_os("windows")
  edges <- data.table::data.table(
    doc_id = c("D1", "D1", "D2", "D2"),
    tf = c("TF1", "TF1", "TF2", "TF2"),
    peak_id = c("P1", "P2", "P1", "P2"),
    gene_key = c("G1", "G2", "G1", "G2")
  )
  score_mat <- matrix(
    c(
      0.9, 0.1, 0.8, 0.2,
      0.2, 0.7, 0.1, 0.9
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(score_mat) <- c("topic_1", "topic_2")
  colnames(score_mat) <- c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2")
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L, 2L, 2L),
    term_id = c("PEAK:P1", "GENE:G1", "PEAK:P2", "GENE:G2"),
    in_topic = TRUE
  )

  serial <- compute_topic_links(
    edges,
    score_mat,
    topic_terms = topic_terms,
    binarize_method = "topn",
    link_method = "gammafit",
    chunk_size = 1L,
    n_cores = 1L,
    overwrite = TRUE
  )
  parallel <- compute_topic_links(
    edges,
    score_mat,
    topic_terms = topic_terms,
    binarize_method = "topn",
    link_method = "gammafit",
    chunk_size = 1L,
    n_cores = 2L,
    overwrite = TRUE
  )

  key_cols <- c("doc_id", "tf", "peak_id", "gene_key", "topic_num")
  data.table::setorderv(serial, key_cols)
  data.table::setorderv(parallel, key_cols)
  expect_equal(parallel, serial, ignore_attr = TRUE)
})
