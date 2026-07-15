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

test_that("VAE training reuses complete K values and trains missing K values", {
  out_dir <- withr::local_tempdir()
  models_dir <- file.path(out_dir, "vae_models")
  dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

  data.table::fwrite(
    data.table::data.table(K = 2L, perplexity = 10, loglik = -10),
    file.path(out_dir, "model_metrics.csv")
  )
  data.table::fwrite(
    data.table::data.table(K = 2L),
    file.path(out_dir, "vae_model_manifest.csv")
  )
  data.table::fwrite(
    data.table::data.table(doc_id = c("D1", "D2"), Topic1 = c(0.7, 0.3), Topic2 = c(0.3, 0.7)),
    file.path(models_dir, "theta_K2.csv")
  )
  data.table::fwrite(
    data.table::data.table(term_id = c("T1", "T2"), Topic1 = c(0.8, 0.2), Topic2 = c(0.2, 0.8)),
    file.path(models_dir, "phi_K2.csv")
  )
  writeLines("old model", file.path(models_dir, "model_K2.pt"))

  fake_trainer <- file.path(out_dir, "fake_vae.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "val <- function(flag) args[match(flag, args) + 1L]",
    "out_dir <- val('--out-dir')",
    "ks <- as.integer(strsplit(val('--k-grid'), ',', fixed = TRUE)[[1]])",
    "dir.create(file.path(out_dir, 'vae_models'), recursive = TRUE, showWarnings = FALSE)",
    "metrics <- data.frame(K = ks, perplexity = 100 + ks, loglik = -100 - ks)",
    "write.csv(metrics, file.path(out_dir, 'model_metrics.csv'), row.names = FALSE)",
    "write.csv(data.frame(K = ks), file.path(out_dir, 'vae_model_manifest.csv'), row.names = FALSE)",
    "for (k in ks) {",
    "  theta <- data.frame(doc_id = c('D1', 'D2'))",
    "  phi <- data.frame(term_id = c('T1', 'T2'))",
    "  for (i in seq_len(k)) { theta[[paste0('Topic', i)]] <- 1 / k; phi[[paste0('Topic', i)]] <- 1 / k }",
    "  write.csv(theta, file.path(out_dir, 'vae_models', paste0('theta_K', k, '.csv')), row.names = FALSE)",
    "  write.csv(phi, file.path(out_dir, 'vae_models', paste0('phi_K', k, '.csv')), row.names = FALSE)",
    "  writeLines('model', file.path(out_dir, 'vae_models', paste0('model_K', k, '.pt')))",
    "}"
  ), fake_trainer)

  doc_term <- data.table::data.table(
    doc_id = rep(c("D1", "D2"), each = 2),
    term_id = rep(c("T1", "T2"), 2),
    pseudo_count_bin = c(1, 1, 1, 1),
    pseudo_count_log = c(1, 1, 1, 1),
    weight = c(1, 1, 1, 1)
  )

  craftgrn:::run_vae_topic_report_py(
    doc_term = doc_term,
    edges_docs = data.table::data.table(),
    out_dir = out_dir,
    option_label = "joint",
    direction_by = "gene",
    vae_script = fake_trainer,
    k_grid = c(2L, 3L),
    vae_python = file.path(R.home("bin"), "Rscript"),
    reuse_if_exists = TRUE,
    do_report = FALSE,
    count_input = "pseudo_count_bin"
  )

  metrics <- data.table::fread(file.path(out_dir, "model_metrics.csv"))
  manifest <- data.table::fread(file.path(out_dir, "vae_model_manifest.csv"))

  expect_equal(sort(metrics$K), c(2L, 3L))
  expect_equal(sort(manifest$K), c(2L, 3L))
  expect_true(file.exists(file.path(models_dir, "theta_K2.csv")))
  expect_true(file.exists(file.path(models_dir, "theta_K3.csv")))
  expect_equal(metrics[K == 2L, perplexity], 10)
  expect_equal(metrics[K == 3L, perplexity], 103)
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

test_that("Module 3 condition links are built directly from Module 2 links", {
  links <- data.table::data.table(
    tf = c("TF1", "TF1", "TF2", "TF3"),
    fp_id = c("P1", "P2", "P1", "P3"),
    target_gene = c("G1", "G2", "G1", "G3"),
    module2_link_pass = TRUE,
    tf_expression_target_r = c(0.8, 0.7, 0.6, 0.5),
    fp_target_rna_r = c(0.9, 0.8, 0.7, 0.6)
  )
  module2 <- list(links = links, module2_links = links)

  fp_score <- matrix(
    c(3, 1, 4, 3, 1, 4),
    nrow = 3,
    dimnames = list(c("P1", "P2", "P3"), c("CondA", "CondB"))
  )
  fp_bound <- fp_score >= 2
  gene_expr <- matrix(
    c(12, 8, 8, 20, 20, 20, 20, 20, 20, 20),
    nrow = 5,
    dimnames = list(c("TF1", "TF2", "TF3", "G1", "G2"), c("CondA", "CondB"))
  )
  gene_on <- gene_expr >= 10
  multiomic <- list(
    schema = "craftgrn_multiomic_v1",
    samples = data.frame(condition_id = c("CondA", "CondB")),
    features = list(
      fp = data.frame(fp_id = c("P1", "P2", "P3")),
      gene = data.frame(gene_id = c("TF1", "TF2", "TF3", "G1", "G2"))
    ),
    matrices = list(
      fp_score = fp_score,
      fp_bound = fp_bound,
      gene_expr = gene_expr,
      gene_on = gene_on
    ),
    refs = list(tf = c("TF1", "TF2", "TF3"))
  )
  class(multiomic) <- c("craftgrn_multiomic", "list")

  root <- withr::local_tempdir()
  out <- module3_prepare_condition_links(
    module2 = module2,
    multiomic_data = multiomic,
    conditions = c("CondA", "CondB"),
    output_dir = root,
    threshold_fp_score = 2,
    threshold_gene_expr = 10,
    threshold_tf_expr = 10,
    output_format = "csv",
    overwrite = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(file.path(root, "condition_links_manifest.csv")))
  expect_equal(sort(out$manifest$condition_id), c("CondA", "CondB"))
  cond_a <- data.table::fread(out$manifest[condition_id == "CondA", path])
  cond_b <- data.table::fread(out$manifest[condition_id == "CondB", path])

  expect_equal(nrow(cond_a), 1L)
  expect_equal(cond_a$doc_id, "CondA::TF1")
  expect_equal(cond_a$gene_key, "G1")
  expect_equal(cond_a$peak_id, "P1")
  expect_equal(cond_a$fp_score_condition, 3)
  expect_equal(cond_a$gene_expr_condition, 20)
  expect_equal(cond_a$tf_expr_condition, 12)

  expect_equal(nrow(cond_b), 2L)
  expect_setequal(cond_b$doc_id, c("CondB::TF1", "CondB::TF2"))
  expect_false("TF3" %in% cond_b$tf)
})

test_that("condition links resolve project-prefixed matrix condition columns", {
  links <- data.table::data.table(
    tf = "TF1",
    fp_id = "P1",
    target_gene = "G1",
    module2_link_pass = TRUE
  )
  module2 <- list(links = links, module2_links = links)
  fp_score <- matrix(3, nrow = 1, dimnames = list("P1", "GSE_TEST_CondA"))
  fp_bound <- matrix(TRUE, nrow = 1, dimnames = list("P1", "GSE_TEST_CondA"))
  gene_expr <- matrix(
    c(12, 20),
    nrow = 2,
    dimnames = list(c("TF1", "G1"), "GSE_TEST_CondA")
  )
  gene_on <- gene_expr >= 10
  multiomic <- list(
    schema = "craftgrn_multiomic_v1",
    project = list(project_id = "GSE_TEST"),
    samples = data.frame(condition_id = "CondA"),
    features = list(
      fp = data.frame(fp_id = "P1"),
      gene = data.frame(gene_id = c("TF1", "G1"))
    ),
    matrices = list(
      fp_score = fp_score,
      fp_bound = fp_bound,
      gene_expr = gene_expr,
      gene_on = gene_on
    ),
    refs = list(tf = "TF1")
  )
  class(multiomic) <- c("craftgrn_multiomic", "list")

  root <- withr::local_tempdir()
  out <- module3_prepare_condition_links(
    module2 = module2,
    multiomic_data = multiomic,
    conditions = "CondA",
    output_dir = root,
    threshold_fp_score = 2,
    threshold_gene_expr = 10,
    threshold_tf_expr = 10,
    output_format = "csv",
    overwrite = TRUE,
    verbose = FALSE
  )
  cond_a <- data.table::fread(out$manifest$path[[1L]])

  expect_equal(out$manifest$condition_id, "CondA")
  expect_equal(out$manifest$matrix_condition_id, "GSE_TEST_CondA")
  expect_equal(cond_a$condition_id, "CondA")
  expect_equal(cond_a$matrix_condition_id, "GSE_TEST_CondA")
  expect_equal(cond_a$doc_id, "CondA::TF1")
})

test_that("condition-native topic training does not require differential link columns", {
  root <- withr::local_tempdir()
  input_dir <- file.path(root, "condition_links")
  model_dir <- file.path(root, "models")
  dir.create(input_dir, recursive = TRUE)

  condition_links <- data.table::data.table(
    condition_id = rep(c("CondA", "CondB"), each = 6),
    condition_label = rep(c("CondA", "CondB"), each = 6),
    doc_id = paste0(rep(c("CondA", "CondB"), each = 6), "::", rep(c("TF1", "TF2"), each = 3, times = 2)),
    tf_doc = rep(c("TF1", "TF2"), each = 3, times = 2),
    tf = rep(c("TF1", "TF2"), each = 3, times = 2),
    gene_key = paste0("G", rep(1:6, times = 2)),
    peak_id = paste0("P", rep(1:6, times = 2)),
    fp_score_condition = c(5, 4, 3, 5, 4, 3, 4, 5, 3, 4, 5, 3),
    gene_expr_condition = c(20, 18, 16, 21, 19, 17, 18, 20, 16, 19, 21, 17),
    tf_expr_condition = 30
  )
  links_path <- file.path(input_dir, "all_condition_links.csv")
  data.table::fwrite(condition_links, links_path)
  data.table::fwrite(
    data.table::data.table(
      condition_id = c("CondA", "CondB"),
      path = links_path,
      format = "csv",
      n_links = c(6L, 6L)
    ),
    file.path(input_dir, "condition_links_manifest.csv")
  )

  train_topic_models(
    Kgrid = 2L,
    input_dir = input_dir,
    output_dir = model_dir,
    analysis_label = "condition_native",
    doc_mode = "tf",
    doc_design = "condition",
    input_source = "condition_links",
    fp_term_mode = "unique",
    gene_term_mode = "unique",
    include_tf_terms = TRUE,
    min_df = 1,
    count_method = "bin",
    backend = "warplda",
    warplda_iterations = 1L,
    local_threads = 1L,
    reuse_if_exists = FALSE,
    threshold_gene_expr = 10,
    threshold_fp_score = 2,
    threshold_tf_expr = 10,
    save_full_doc_term_csv = FALSE,
    flat_output = TRUE
  )

  summary <- data.table::fread(file.path(model_dir, "topic_input_summary.csv"))
  doc_term <- readRDS(file.path(model_dir, "rds", "doc_term.rds"))

  expect_equal(summary$input_source, "condition_links")
  expect_equal(summary$n_link_rows_after_filter, 12)
  expect_setequal(unique(sub("::.*$", "", doc_term$doc_id)), c("CondA", "CondB"))
  expect_true(file.exists(file.path(model_dir, "vae_models", "theta_K2.csv")))
})

test_that("module3_construct_docs accepts condition-link inputs", {
  root <- withr::local_tempdir()
  input_dir <- file.path(root, "condition_links")
  docs_dir <- file.path(root, "topic_documents")
  dir.create(input_dir, recursive = TRUE)

  condition_links <- data.table::data.table(
    condition_id = rep(c("CondA", "CondB"), each = 2),
    condition_label = rep(c("CondA", "CondB"), each = 2),
    doc_id = c("CondA::TF1", "CondA::TF2", "CondB::TF1", "CondB::TF2"),
    tf_doc = c("TF1", "TF2", "TF1", "TF2"),
    tf = c("TF1", "TF2", "TF1", "TF2"),
    gene_key = c("G1", "G2", "G1", "G2"),
    peak_id = c("P1", "P2", "P1", "P2"),
    fp_score_condition = c(5, 4, 5, 4),
    gene_expr_condition = c(20, 18, 20, 18),
    tf_expr_condition = c(30, 30, 30, 30)
  )
  links_path <- file.path(input_dir, "condition_links.csv")
  data.table::fwrite(condition_links, links_path)
  data.table::fwrite(
    data.table::data.table(
      condition_id = c("CondA", "CondB"),
      path = links_path,
      format = "csv",
      n_links = c(2L, 2L)
    ),
    file.path(input_dir, "condition_links_manifest.csv")
  )

  out <- module3_construct_docs(
    filtered_dir = input_dir,
    output_dir = docs_dir,
    input_source = "condition_links",
    doc_mode = "tf",
    doc_design = "condition",
    fp_term_mode = "unique",
    gene_term_mode = "unique",
    include_tf_terms = TRUE,
    count_method = "log",
    threshold_gene_expr = 10,
    threshold_fp_score = 2,
    threshold_tf_expr = 10,
    min_df = 1,
    overwrite = TRUE,
    verbose = FALSE
  )

  summary <- data.table::fread(file.path(docs_dir, "topic_input_summary.csv"))
  doc_term <- readRDS(file.path(docs_dir, "rds", "doc_term.rds"))

  expect_false(isTRUE(out$reused))
  expect_equal(summary$input_source, "condition_links")
  expect_setequal(unique(sub("::.*$", "", doc_term$doc_id)), c("CondA", "CondB"))
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

test_that("topic score methods expose specificity and legacy rowmax scores", {
  phi <- matrix(
    c(
      0.6, 0.3, 0.1,
      0.6, 0.05, 0.35
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:Common", "GENE:T1", "GENE:T2"))
  )

  rowmax <- score_terms_normtop(phi, method = "rowmax_phi")
  expect_equal(rowmax, phi / apply(phi, 1, max), tolerance = 1e-12)

  specificity <- score_terms_normtop(phi, method = "normtop_specificity")
  expect_equal(dim(specificity), dim(phi))
  expect_equal(rownames(specificity), rownames(phi))
  expect_equal(colnames(specificity), colnames(phi))
  expect_equal(specificity["Topic1", "GENE:T1"], 1, tolerance = 1e-8)
  expect_equal(specificity["Topic2", "GENE:T2"], 1, tolerance = 1e-8)
  expect_equal(specificity["Topic1", "GENE:Common"], 0, tolerance = 1e-8)
  expect_equal(specificity["Topic2", "GENE:Common"], 0, tolerance = 1e-8)
  expect_lt(specificity["Topic1", "GENE:Common"], rowmax["Topic1", "GENE:Common"])
  expect_lt(specificity["Topic2", "GENE:Common"], rowmax["Topic2", "GENE:Common"])
})

test_that("gammafit diagnostics report score method and minimum term forcing", {
  score_mat <- matrix(
    c(
      seq(0.01, 0.12, length.out = 12),
      seq(0.12, 0.01, length.out = 12)
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(score_mat) <- c("Topic1", "Topic2")
  colnames(score_mat) <- paste0("GENE:G", seq_len(12))
  topic_terms <- binarize_topics(
    score_mat,
    method = "gammafit",
    thrP = 0.9999,
    min_terms = 3L,
    gammafit_scope = "topic_term_group"
  )

  diagnostics <- .gammafit_diagnostics_by_termclass(
    score_mat,
    topic_terms = topic_terms,
    topic_score_method = "normtop_specificity",
    thrP = 0.9999,
    min_terms = 3L,
    gammafit_scope = "topic_term_group"
  )

  expect_s3_class(diagnostics, "data.table")
  expect_true(all(c(
    "topic_score_method", "gammafit_scope", "topic_num", "term_group",
    "positive_count", "zero_fraction", "gamma_shape", "gamma_rate",
    "gamma_cutoff", "selected_by_gamma", "selected_after_min_terms",
    "forced_min_terms"
  ) %in% names(diagnostics)))
  gene_rows <- diagnostics[term_group == "GENE"]
  expect_equal(nrow(gene_rows), 2L)
  expect_true(all(gene_rows$topic_score_method == "normtop_specificity"))
  expect_true(all(gene_rows$positive_count == 12L))
  expect_true(all(is.finite(gene_rows$zero_fraction)))
  expect_true(all(is.finite(gene_rows$gamma_shape)))
  expect_true(all(is.finite(gene_rows$gamma_rate)))
  expect_true(all(gene_rows$selected_after_min_terms >= 3L))
  expect_true(any(gene_rows$forced_min_terms))
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

test_that("topic-link preflight blocks excessive Cartesian scoring", {
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
  withr::local_options(list(craftgrn.topic_link.max_scored_rows = 3))

  expect_error(
    compute_topic_links(
      edges,
      score_mat,
      binarize_method = "topn",
      link_method = "gammafit",
      overwrite = TRUE
    ),
    "would materialize too many edge-topic rows"
  )
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
  coverage <- data.table::fread(file.path(out_dir, "topic_item_coverage_counts.csv"))

  expect_equal(summary$n_scored_rows, 4)
  expect_equal(summary$n_pass_rows, 1)
  expect_setequal(coverage$unit, c("Terms", "Genes", "TF-gene-doc links", "Links", "TFs"))
  expect_equal(coverage[unit == "Terms" & status == "Pass", count], 2)
  expect_equal(coverage[unit == "Terms" & status == "Fail", count], 2)
  expect_equal(coverage[unit == "Terms" & status == "Pass", fraction], 0.5)
  expect_equal(coverage[unit == "Terms" & status == "Pass", percent], 50)
  expect_equal(coverage[unit == "Genes" & status == "Pass", count], 1)
  expect_equal(coverage[unit == "Genes" & status == "Fail", count], 1)
  expect_equal(coverage[unit == "Genes" & status == "Pass", percent], 50)
  expect_equal(coverage[unit == "Links" & status == "Pass", count], 1)
  expect_equal(coverage[unit == "Links" & status == "Fail", count], 1)
  expect_equal(coverage[unit == "Links" & status == "Pass", percent], 50)
  expect_equal(coverage[unit == "TF-gene-doc links" & status == "Pass", count], 1)
  expect_equal(coverage[unit == "TF-gene-doc links" & status == "Fail", count], 1)
  expect_equal(coverage[unit == "TF-gene-doc links" & status == "Pass", percent], 50)
  expect_equal(coverage[unit == "TFs" & status == "Pass", count], 1)
  expect_equal(coverage[unit == "TFs" & status == "Fail", count], 0)
  expect_equal(coverage[unit == "TFs" & status == "Pass", percent], 100)
})

test_that("term coverage summary lines report model and term-group percentages", {
  score_mat <- matrix(
    c(
      0.9, 0.8, 0.7, 0.1,
      0.2, 0.1, 0.3, 0.4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Topic1", "Topic2"),
      c("GENE:G1", "GENE:G2", "PEAK:P1", "PEAK:P2")
    )
  )
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L, 1L),
    term_id = c("GENE:G1", "PEAK:P1", "PEAK:P2"),
    in_topic = TRUE
  )

  lines <- .topic_term_coverage_summary_lines(topic_terms, score_mat)

  expect_equal(
    lines,
    c(
      "Terms: 3 / 4 = 75.00%",
      "GENE terms: 1 / 2 = 50.00%",
      "PEAK terms: 2 / 2 = 100.00%"
    )
  )
})

test_that("assignment coverage summary table includes term groups and item coverage", {
  score_mat <- matrix(
    c(
      0.9, 0.8, 0.7, 0.1,
      0.2, 0.1, 0.3, 0.4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Topic1", "Topic2"),
      c("GENE:G1", "GENE:G2", "PEAK:P1", "PEAK:P2")
    )
  )
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L, 1L),
    term_id = c("GENE:G1", "PEAK:P1", "PEAK:P2"),
    in_topic = TRUE
  )
  item_coverage <- data.table::data.table(
    unit = c("Genes", "Genes", "TF-gene-doc links", "TF-gene-doc links", "Links", "Links", "TFs", "TFs"),
    status = rep(c("Pass", "Fail"), 4),
    count = c(4L, 6L, 5L, 5L, 12L, 8L, 3L, 1L),
    total = c(10L, 10L, 10L, 10L, 20L, 20L, 4L, 4L),
    percent = c(40, 60, 50, 50, 60, 40, 75, 25)
  )

  tbl <- .topic_assignment_coverage_summary_table(topic_terms, score_mat, item_coverage)
  tbl_no_peak_expanded <- .topic_assignment_coverage_summary_table(
    topic_terms,
    score_mat,
    item_coverage,
    show_peak_expanded_link_coverage = FALSE
  )

  expect_equal(tbl[label == "Terms", label_text], "3 / 4 = 75.00%")
  expect_equal(tbl[label == "GENE terms", label_text], "1 / 2 = 50.00%")
  expect_equal(tbl[label == "PEAK terms", label_text], "2 / 2 = 100.00%")
  expect_equal(tbl[label == "Genes", label_text], "4 / 10 = 40.00%")
  expect_equal(tbl[label == "TF-gene-doc links", label_text], "5 / 10 = 50.00%")
  expect_equal(tbl[label == "TF-peak-gene links", label_text], "12 / 20 = 60.00%")
  expect_equal(tbl[label == "TFs", label_text], "3 / 4 = 75.00%")
  expect_false("TF-peak-gene links" %in% tbl_no_peak_expanded$label)
  expect_true("TF-gene-doc links" %in% tbl_no_peak_expanded$label)
})

test_that("theta_and_terms requires both document theta and topic terms", {
  edges <- data.table::data.table(
    doc_id = c("D1", "D2"),
    tf = c("TF1", "TF2"),
    peak_id = c("P1", "P2"),
    gene_key = c("G1", "G2")
  )
  score_mat <- matrix(
    c(
      0.9, 0.8, 0.7, 0.6,
      0.1, 0.2, 0.3, 0.4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Topic1", "Topic2"),
      c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2")
    )
  )
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L, 1L, 1L),
    term_id = c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2"),
    in_topic = TRUE
  )
  theta <- matrix(
    c(0.5, 0.5, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("D1", "D2"), c("Topic1", "Topic2"))
  )
  out_dir <- tempfile("module3-theta-and-terms-")
  dir.create(out_dir, recursive = TRUE)

  res <- compute_topic_links(
    edges,
    score_mat,
    topic_terms = topic_terms,
    theta = theta,
    topic_tf_membership_cutoff = 0.3,
    binarize_method = "topn",
    link_method = "theta_and_terms",
    pass_file = file.path(out_dir, "topic_links_pass.csv"),
    output_mode = "pass",
    overwrite = TRUE
  )

  pass <- data.table::fread(file.path(out_dir, "topic_links_pass.csv"))
  coverage <- data.table::fread(file.path(out_dir, "topic_item_coverage_counts.csv"))
  expect_equal(nrow(res), 1L)
  expect_equal(nrow(pass), 1L)
  expect_equal(pass$doc_id, "D1")
  expect_equal(pass$topic_num, 1L)
  expect_equal(pass$theta, 0.5)
  expect_equal(pass$theta_cutoff, 0.3)
  expect_true(pass$theta_pass)
  expect_equal(coverage[unit == "Links" & status == "Pass", count], 1)
  expect_equal(coverage[unit == "Links" & status == "Fail", count], 1)
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

test_that("per-comparison pathway gene sets use theta documents and topic terms", {
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L, 1L, 2L),
    term_id = c("GENE:G1", "GENE:G2", "PEAK:G3", "GENE:G4"),
    in_topic = c(TRUE, TRUE, TRUE, TRUE)
  )
  theta <- matrix(
    c(
      0.8, 0.2,
      0.1, 0.9,
      0.7, 0.3
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c(
        "CmpA::TF1::Target-Up",
        "CmpA::TF2::Target-Down",
        "CmpB::TF1::Target-Up"
      ),
      c("Topic1", "Topic2")
    )
  )
  edges <- data.table::data.table(
    doc_id = c(
      "CmpA::TF1::Target-Up",
      "CmpA::TF1::Target-Up",
      "CmpA::TF1::Target-Up",
      "CmpA::TF2::Target-Down",
      "CmpB::TF1::Target-Up"
    ),
    tf = c("TF1", "TF1", "TF1", "TF2", "TF1"),
    peak_id = paste0("P", 1:5),
    gene_key = c("G1", "G3", "G_not_topic", "G4", "G2")
  )

  out <- topic_gene_sets_by_comparison_terms(
    topic_terms = topic_terms,
    edges_docs = edges,
    theta = theta,
    theta_min = 0.3,
    include_peak_terms = TRUE,
    doc_design = "comparison"
  )

  expect_setequal(out[comparison_id == "CmpA" & direction_group == "Up" & topic == 1L, gene], c("G1", "G3"))
  expect_equal(out[comparison_id == "CmpA" & direction_group == "Down" & topic == 2L, gene], "G4")
  expect_equal(out[comparison_id == "CmpB" & direction_group == "Up" & topic == 1L, gene], "G2")
  expect_false("G_not_topic" %in% out$gene)
})

test_that("overall pathway gene sets use assigned topic terms before link projection", {
  topic_terms <- data.table::data.table(
    topic = c(1L, 1L, 1L, 2L, 2L),
    term_id = c("GENE:G1", "PEAK:G3", "PEAK:P4", "GENE:G4", "GENE:G_fail"),
    in_topic = c(TRUE, TRUE, TRUE, TRUE, FALSE)
  )
  edges <- data.table::data.table(
    peak_id = c("P4", "P5"),
    gene_key = c("G5", "G_not_topic")
  )
  gene_sets <- topic_gene_sets_from_terms(
    topic_terms = topic_terms,
    edges_docs = edges,
    option_label = "joint",
    use_all_terms = FALSE,
    include_peak_terms = TRUE
  )

  expect_setequal(gene_sets[["1"]], c("G1", "G3", "G5"))
  expect_equal(gene_sets[["2"]], "G4")
  expect_false("G_fail" %in% unlist(gene_sets, use.names = FALSE))
  expect_false("G_not_topic" %in% unlist(gene_sets, use.names = FALSE))
})

test_that("overall topic pathway dotplot shows a readable expanded default", {
  expect_equal(formals(plot_topic_pathway_enrichment_heatmap)$dot_top_n_per_topic, 25L)
  expect_equal(formals(run_tfdocs_report_from_topic_base)$dot_top_n_per_topic, 25L)
})

test_that("topic-link enumeration is explicit opt-in", {
  defaults <- make_topic_report_args_simple(
    thrP = 0.9,
    link_prob_cutoff = 0.3,
    link_fdr_p = NA_real_
  )
  requested <- make_topic_report_args_simple(
    thrP = 0.9,
    link_prob_cutoff = 0.3,
    link_fdr_p = NA_real_,
    extraction_steps = "topic_links"
  )

  expect_false(defaults$run_link_topic_scores)
  expect_true(requested$run_link_topic_scores)
})

test_that("Module 3 K extraction worker planning preserves RAM headroom", {
  gib <- 1024^3
  root <- tempfile("module3-k-worker-plan-")
  dir.create(root, recursive = TRUE)

  roomy <- .module3_extraction_k_worker_plan(
    n_tasks = 8L,
    model_dir = root,
    k_max_workers = 4L,
    k_memory_gb = 16,
    k_memory_reserve_gb = 32,
    cores = 36L,
    available_bytes = 128 * gib,
    os_type = "unix"
  )
  tight <- .module3_extraction_k_worker_plan(
    n_tasks = 8L,
    model_dir = root,
    k_max_workers = 4L,
    k_memory_gb = 16,
    k_memory_reserve_gb = 32,
    cores = 36L,
    available_bytes = 40 * gib,
    os_type = "unix"
  )
  link_scoring <- .module3_extraction_k_worker_plan(
    n_tasks = 8L,
    model_dir = root,
    k_memory_gb = 16,
    cores = 36L,
    link_scoring = TRUE,
    available_bytes = 128 * gib,
    os_type = "unix"
  )

  expect_equal(roomy$workers, 4L)
  expect_equal(tight$workers, 1L)
  expect_equal(link_scoring$workers, 1L)
  expect_match(link_scoring$reason, "topic-link scoring", fixed = TRUE)
})

test_that("Module 3 K extraction workers keep worker output quiet", {
  tasks <- list(list(k = 2L), list(k = 3L))
  expect_silent(
    .module3_run_extraction_k_tasks(
      tasks = tasks,
      worker_fun = function(task) {
        message("worker message")
        task$k
      },
      workers = if (.Platform$OS.type == "unix") 2L else 1L,
      cores = 2L
    )
  )
})

test_that("module3_extract_topics supports a safe multi-K extraction", {
  testthat::skip_if_not_installed("Matrix")
  root <- tempfile("module3-multi-k-extract-")
  model_dir <- file.path(root, "models")
  output_dir <- file.path(root, "extraction")
  dir.create(file.path(model_dir, "vae_models"), recursive = TRUE)
  dir.create(file.path(model_dir, "rds"), recursive = TRUE)

  theta2 <- data.table::data.table(
    doc_id = c("CondA::TF1", "CondB::TF1"),
    Topic1 = c(0.8, 0.2),
    Topic2 = c(0.2, 0.8)
  )
  theta3 <- data.table::data.table(
    doc_id = c("CondA::TF1", "CondB::TF1"),
    Topic1 = c(0.7, 0.1),
    Topic2 = c(0.2, 0.7),
    Topic3 = c(0.1, 0.2)
  )
  phi2 <- data.table::data.table(
    term_id = c("Topic1", "Topic2"),
    `PEAK:G1` = c(0.8, 0.2),
    `PEAK:G2` = c(0.2, 0.8),
    `GENE:G1` = c(0.8, 0.2),
    `GENE:G2` = c(0.2, 0.8)
  )
  phi3 <- data.table::data.table(
    term_id = c("Topic1", "Topic2", "Topic3"),
    `PEAK:G1` = c(0.7, 0.2, 0.1),
    `PEAK:G2` = c(0.1, 0.7, 0.2),
    `GENE:G1` = c(0.7, 0.2, 0.1),
    `GENE:G2` = c(0.1, 0.7, 0.2)
  )
  data.table::fwrite(theta2, file.path(model_dir, "vae_models", "theta_K2.csv"))
  data.table::fwrite(theta3, file.path(model_dir, "vae_models", "theta_K3.csv"))
  data.table::fwrite(phi2, file.path(model_dir, "vae_models", "phi_K2.csv"))
  data.table::fwrite(phi3, file.path(model_dir, "vae_models", "phi_K3.csv"))

  dtm <- Matrix::Matrix(matrix(c(3, 1, 2, 1, 1, 2, 1, 3), nrow = 2, byrow = TRUE), sparse = TRUE)
  rownames(dtm) <- theta2$doc_id
  colnames(dtm) <- c("PEAK:G1", "PEAK:G2", "GENE:G1", "GENE:G2")
  edges <- data.table::data.table(
    doc_id = theta2$doc_id,
    tf = "TF1",
    peak_id = c("P1", "P2"),
    gene_key = c("G1", "G2")
  )
  saveRDS(dtm, file.path(model_dir, "rds", "dtm.rds"))
  saveRDS(edges, file.path(model_dir, "rds", "edges_docs.rds"))

  module3_extract_topics(
    k = c(2L, 3L),
    model_dir = model_dir,
    output_dir = output_dir,
    backend = "warplda",
    doc_mode = "tf",
    weight_label = "peak_score_gene_expr",
    k_workers = 2L,
    k_max_workers = 2L,
    k_memory_gb = 1,
    k_memory_reserve_gb = 0,
    cores = 2L,
    verbose = FALSE,
    topic_report_args = list(
      fp_term_mode = "aggregate",
      in_topic_min_terms = 1L,
      run_link_topic_scores = FALSE,
      run_gammafit_summary = FALSE,
      run_pathway_enrichment = FALSE,
      run_topic_term_heatmap = FALSE,
      run_topic_by_comparison_heatmaps = FALSE,
      run_intertopic_distance_map = FALSE
    )
  )

  expect_true(file.exists(file.path(output_dir, "K2", "topic_terms.csv")))
  expect_true(file.exists(file.path(output_dir, "K3", "topic_terms.csv")))
  expect_true(file.exists(file.path(output_dir, "K2", "topic_extraction_step_timing.csv")))
  expect_true(file.exists(file.path(output_dir, "K3", "topic_extraction_step_timing.csv")))
})

test_that("overall topic pathway output no longer writes a heatmap PDF", {
  testthat::skip_if_not_installed("pheatmap")

  out_dir <- file.path(tempdir(), paste0("topic-pathway-no-heatmap-", sample.int(1e8, 1L)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "topic_pathway_enrichment_heatmap.pdf")
  topic_terms <- data.table::data.table(
    topic = c("Topic1", "Topic1", "Topic1"),
    topic_num = c(1L, 1L, 1L),
    term_id = c("GENE:A", "GENE:B", "GENE:C"),
    score = c(0.8, 0.7, 0.6),
    in_topic = TRUE
  )

  with_mocked_bindings(
    plot_topic_pathway_enrichment_heatmap(
      topic_terms = topic_terms,
      edges_docs = NULL,
      option_label = "opt3_gene_fc_expr",
      out_file = out_file,
      dbs = "Toy_DB",
      pathway_species = "human",
      min_genes = 2L,
      top_n_per_topic = Inf,
      make_heatmap = TRUE,
      make_dotplot = FALSE,
      pathway_backend = "enrichly"
    ),
    .pathway_backend_available = function(pathway_backend) TRUE,
    .run_enrichr_cached = function(...) {
      list(Toy_DB = data.frame(
        Term = "Term one",
        Adjusted.P.value = 0.01,
        P.value = 0.001,
        Overlap = "2/10",
        Genes = "A;B",
        Combined.Score = 8,
        Odds.Ratio = 2,
        query_size = 3L,
        term_size = 10L,
        background_size = 20000L,
        check.names = FALSE
      ))
    }
  )

  expect_true(file.exists(file.path(out_dir, "topic_pathway_enrichment_topic_terms.csv")))
  expect_false(file.exists(out_file))
})

test_that("human mouse pathway best selector keeps one best species row per pathway", {
  dt <- data.table::data.table(
    topic = c(1L, 1L, 1L, 2L),
    pathway = c(
      "KEGG: T cell receptor signaling pathway",
      "KEGG: T cell receptor signaling pathway",
      "Reactome: Interferon Signaling",
      "WikiPathways: Apoptosis"
    ),
    database = c("KEGG_2021_Human", "KEGG_2019_Mouse", "Reactome_2022", "WikiPathways_2024_Mouse"),
    database_label = c("KEGG", "KEGG", "Reactome", "WikiPathways"),
    pathway_term = c(
      "T cell receptor signaling pathway",
      "T cell receptor signaling pathway",
      "Interferon Signaling",
      "Apoptosis"
    ),
    pathway_species = c("human", "mouse", "human", "mouse"),
    padj = c(0.02, 0.01, 0.03, 0.04),
    pval = c(0.002, 0.001, 0.003, 0.004),
    overlap = c("2/10", "3/11", "1/8", "2/9"),
    overlap_hits = c(2L, 3L, 1L, 2L),
    genes = c("A;B", "A;B;C", "D", "E;F"),
    logp = -log10(c(0.02, 0.01, 0.03, 0.04)),
    combined_score = c(10, 8, 4, 5),
    odds_ratio = c(2, 3, 1, 2),
    cluster_size = c(20L, 20L, 20L, 30L),
    query_size = c(20L, 20L, 20L, 30L),
    term_size = c(10L, 11L, 8L, 9L),
    background_size = 20000L
  )

  out <- .select_best_human_mouse_pathways(dt)
  tcr <- out[topic == 1L & pathway_norm_key == "kegg:t cell receptor signaling pathway"]

  expect_equal(nrow(out), 3L)
  expect_equal(nrow(tcr), 1L)
  expect_equal(tcr$selected_pathway_species, "mouse")
  expect_equal(tcr$selected_database, "KEGG_2019_Mouse")
  expect_equal(tcr$human_padj, 0.02)
  expect_equal(tcr$mouse_padj, 0.01)
  expect_equal(tcr$human_overlap_hits, 2L)
  expect_equal(tcr$mouse_overlap_hits, 3L)
})

test_that("per-comparison topic-term pathway wrapper retests overall pathway genes", {
  topic_terms <- data.table::data.table(
    topic = c("Topic1", "Topic1", "Topic1", "Topic2"),
    topic_num = c(1L, 1L, 1L, 2L),
    term_id = c("GENE:Atp2b4", "GENE:Selplg", "GENE:G4", "GENE:G5"),
    in_topic = TRUE
  )
  theta <- matrix(
    c(
      0.72, 0.28,
      0.25, 0.75
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("CmpA::TF1::Target-Up", "CmpA::TF2::Target-Down"),
      c("Topic1", "Topic2")
    )
  )
  edges <- data.table::data.table(
    doc_id = c(
      "CmpA::TF1::Target-Up",
      "CmpA::TF1::Target-Up",
      "CmpA::TF1::Target-Up",
      "CmpA::TF2::Target-Down"
    ),
    tf = c("TF1", "TF1", "TF1", "TF2"),
    peak_id = paste0("P", 1:4),
    gene_key = c("Atp2b4", "Selplg", "G_not_topic", "G5")
  )
  out_dir <- file.path(tempdir(), paste0("topic-term-pathway-", sample.int(1e8, 1L)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      topic = c(1L, 1L, 2L),
      pathway = c("Reactome: Path A", "Reactome: Path B", "Reactome: Path C"),
      padj = c(0.01, 0.02, 0.03),
      pval = c(0.001, 0.002, 0.003),
      overlap = c("3/10", "1/8", "1/9"),
      overlap_hits = c(3L, 1L, 1L),
      genes = c("ATP2B4;SELPLG;G3", "G4", "G5"),
      logp = -log10(c(0.01, 0.02, 0.03)),
      combined_score = c(10, 5, 3),
      odds_ratio = c(2, 1.5, 1.2),
      cluster_size = c(4L, 4L, 1L),
      term_size = NA,
      background_size = NA
    ),
    file.path(out_dir, "topic_pathway_enrichment_topic_terms.csv")
  )

  plot_topic_pathway_enrichment_by_comparison_terms(
    topic_terms = topic_terms,
    edges_docs = edges,
    theta = theta,
    out_dir = out_dir,
    pathway_backend = "enrichly",
    background_size = 20000L
  )

  out <- data.table::fread(file.path(out_dir, "per_comparison_topic_pathway_enrichment.csv"))
  expect_true(all(c("comparison_id", "direction_group", "topic", "pathway", "pval", "padj") %in% names(out)))
  up_path_a <- out[comparison_id == "CmpA" & direction_group == "Up" & topic == 1L & pathway == "Reactome: Path A"]
  expect_equal(up_path_a$query_size, 2L)
  expect_equal(up_path_a$overlap_hits, 2L)
  expect_equal(up_path_a$overlap_genes, "Atp2b4;Selplg")
  expect_equal(
    up_path_a$pval,
    stats::phyper(2 - 1, 10, 20000 - 10, 2, lower.tail = FALSE),
    tolerance = 1e-12
  )
  down_path_c <- out[comparison_id == "CmpA" & direction_group == "Down" & topic == 2L & pathway == "Reactome: Path C"]
  expect_equal(down_path_c$query_size, 1L)
  expect_equal(down_path_c$overlap_hits, 1L)
  expect_false(file.exists(file.path(out_dir, "per_comparison_pathway_topic_terms")))
  expect_false(file.exists(file.path(out_dir, "CmpA_Up_topic_term_dotplot.pdf")))
})

test_that("per-comparison pathway retest resolves formal human aliases", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  topic_terms <- data.table::data.table(
    topic = "Topic1",
    topic_num = 1L,
    term_id = "GENE:P53",
    in_topic = TRUE,
    score = 0.9
  )
  theta <- matrix(
    0.9,
    nrow = 1,
    dimnames = list("CmpAlias::TF1::Target-Up", "Topic1")
  )
  edges <- data.table::data.table(
    doc_id = "CmpAlias::TF1::Target-Up",
    tf = "TF1",
    peak_id = "peak1",
    gene_key = "P53"
  )
  out_dir <- file.path(tempdir(), paste0("topic-term-pathway-alias-", sample.int(1e8, 1L)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      topic = 1L,
      pathway = "Reactome: TP53 Pathway",
      padj = 0.01,
      pval = 0.001,
      overlap = "1/20",
      overlap_hits = 1L,
      genes = "TP53",
      logp = 2,
      combined_score = 10,
      odds_ratio = 2,
      cluster_size = 1L,
      term_size = NA,
      background_size = NA
    ),
    file.path(out_dir, "topic_pathway_enrichment_topic_terms.csv")
  )

  plot_topic_pathway_enrichment_by_comparison_terms(
    topic_terms = topic_terms,
    edges_docs = edges,
    theta = theta,
    out_dir = out_dir,
    pathway_backend = "enrichly",
    pathway_species = "human",
    background_size = 20000L
  )

  out <- data.table::fread(file.path(out_dir, "per_comparison_topic_pathway_enrichment.csv"))
  row <- out[comparison_id == "CmpAlias" & direction_group == "Up" & topic == 1L]
  expect_equal(row$query_size, 1L)
  expect_equal(row$overlap_hits, 1L)
  expect_equal(row$overlap_genes, "P53")
  expect_equal(row$overlap_gene_symbols, "TP53")
  expect_match(row$gene_match_summary, "alias=1", fixed = TRUE)
})

test_that("condition pathway retest writes condition-labeled outputs", {
  topic_terms <- data.table::data.table(
    topic = c("Topic1", "Topic1", "Topic2"),
    topic_num = c(1L, 1L, 2L),
    term_id = c("GENE:GeneA", "GENE:GeneB", "GENE:GeneC"),
    in_topic = TRUE,
    score = c(0.9, 0.8, 0.7)
  )
  theta <- matrix(
    c(
      0.8, 0.2,
      0.7, 0.3,
      0.1, 0.9
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("CondA::TF1", "CondA::TF2", "CondB::TF1"),
      c("Topic1", "Topic2")
    )
  )
  edges <- data.table::data.table(
    doc_id = c("CondA::TF1", "CondA::TF2", "CondB::TF1"),
    tf = c("TF1", "TF2", "TF1"),
    peak_id = paste0("P", 1:3),
    gene_key = c("GeneA", "GeneB", "GeneC")
  )
  out_dir <- file.path(tempdir(), paste0("condition-topic-pathway-", sample.int(1e8, 1L)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      topic = c(1L, 2L),
      pathway = c("Reactome: Condition Path A", "Reactome: Condition Path B"),
      padj = c(0.01, 0.02),
      pval = c(0.001, 0.002),
      overlap = c("2/10", "1/8"),
      overlap_hits = c(2L, 1L),
      genes = c("GENEA;GENEB", "GENEC"),
      logp = -log10(c(0.01, 0.02)),
      combined_score = c(10, 5),
      odds_ratio = c(2, 1.5),
      cluster_size = c(2L, 1L),
      term_size = NA,
      background_size = NA
    ),
    file.path(out_dir, "topic_pathway_enrichment_topic_terms.csv")
  )

  plot_topic_pathway_enrichment_by_condition_terms(
    topic_terms = topic_terms,
    edges_docs = edges,
    theta = theta,
    out_dir = out_dir,
    pathway_backend = "enrichly",
    pathway_species = "human",
    background_size = 20000L
  )

  out_file <- file.path(out_dir, "per_condition_topic_pathway_enrichment.csv")
  expect_true(file.exists(out_file))
  expect_true(file.exists(file.path(out_dir, "per_condition_topic_pathway_debug.txt")))
  expect_false(file.exists(file.path(out_dir, "per_comparison_topic_pathway_enrichment.csv")))
  out <- data.table::fread(out_file)
  expect_true(all(c(
    "condition_id", "direction_group", "topic", "pathway", "pval", "padj",
    "condition_topic_genes", "condition_topic_gene_symbols"
  ) %in% names(out)))
  expect_false("comparison_id" %in% names(out))
  expect_false("comparison_topic_genes" %in% names(out))
  row <- out[condition_id == "CondA" & direction_group == "All" & topic == 1L & pathway == "Reactome: Condition Path A"]
  expect_equal(row$query_size, 2L)
  expect_equal(row$overlap_hits, 2L)
  expect_equal(row$overlap_genes, "GeneA;GeneB")
  expect_equal(
    row$pval,
    stats::phyper(2 - 1, 10, 20000 - 10, 2, lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("TF topic assignment keeps direction-specific theta memberships", {
  theta <- matrix(
    c(
      0.62, 0.31, 0.07,
      0.12, 0.76, 0.12,
      0.34, 0.33, 0.33,
      0.20, 0.18, 0.62
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c(
        "CompA::KLF5::Target-Up",
        "CompA::KLF5::Target-Down",
        "CompA::SOX9::Target-Up",
        "CompB::Tbx21::Target-Up"
      ),
      c("Topic1", "Topic2", "Topic3")
    )
  )

  res <- .build_tf_topic_assignment_tables(
    theta = theta,
    doc_design = "comparison",
    membership_cutoff = 0.3,
    primary_margin_cutoff = 0.1
  )

  expect_equal(nrow(res$membership), 12L)
  expect_true(all(c("membership", "primary", "pass", "direction_summary") %in% names(res)))

  klf_up <- res$primary[doc_id == "CompA::KLF5::Target-Up"]
  klf_down <- res$primary[doc_id == "CompA::KLF5::Target-Down"]
  expect_equal(klf_up$primary_topic, "Topic1")
  expect_equal(klf_down$primary_topic, "Topic2")
  expect_false(klf_up$ambiguous)
  expect_false(klf_down$ambiguous)

  sox9 <- res$primary[doc_id == "CompA::SOX9::Target-Up"]
  expect_true(sox9$ambiguous)
  expect_equal(res$pass[doc_id == "CompA::SOX9::Target-Up", .N], 3L)

  summary <- res$direction_summary[comparison_id == "CompA" & tf == "KLF5"]
  expect_equal(summary$direction_topic_status, "different_topic")
  expect_equal(summary$up_primary_topic, "Topic1")
  expect_equal(summary$down_primary_topic, "Topic2")

  expect_equal(res$primary[tf == "Tbx21", tf_display], "Tbet")
})

test_that("TF topic assignment supports condition TF documents", {
  theta <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.15, 0.75, 0.10,
      0.20, 0.25, 0.55
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("CondA::KLF5", "CondB::KLF5", "CondB::Tbx21"),
      c("Topic1", "Topic2", "Topic3")
    )
  )

  res <- .build_tf_topic_assignment_tables(
    theta = theta,
    doc_design = "condition",
    membership_cutoff = 0.3,
    primary_margin_cutoff = 0.1
  )

  expect_equal(nrow(res$membership), 9L)
  expect_equal(nrow(res$direction_summary), 0L)
  expect_equal(res$primary[doc_id == "CondA::KLF5", comparison_id], "CondA")
  expect_equal(res$primary[doc_id == "CondA::KLF5", primary_topic], "Topic1")
  expect_equal(res$primary[doc_id == "CondB::KLF5", primary_topic], "Topic2")
  expect_equal(res$primary[doc_id == "CondB::Tbx21", tf_display], "Tbet")
  expect_true(all(is.na(res$primary$direction)))
})

test_that("TF coverage is counted from document theta assignment", {
  theta <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.20, 0.25, 0.55,
      0.25, 0.25, 0.25
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("CondA::KLF5", "CondB::KLF5", "CondA::SOX9"),
      c("Topic1", "Topic2", "Topic3")
    )
  )
  assign <- .build_tf_topic_assignment_tables(
    theta = theta,
    doc_design = "condition",
    membership_cutoff = 0.3,
    primary_margin_cutoff = 0.1
  )

  coverage <- .topic_item_coverage_from_tf_assignment(assign)

  expect_equal(coverage[unit == "TFs" & status == "Pass", count], 1)
  expect_equal(coverage[unit == "TFs" & status == "Fail", count], 1)
  expect_equal(coverage[unit == "TFs" & status == "Pass", total], 2)
  expect_equal(coverage[unit == "TFs" & status == "Pass", percent], 50)
  expect_equal(
    coverage[unit == "TFs" & status == "Pass", count_basis],
    "TFs assigned to at least one topic from raw theta documents"
  )
})

test_that("raw theta document heatmap writes comparison and condition PDF files", {
  theta_cmp <- matrix(
    c(
      0.60, 0.30, 0.10,
      0.10, 0.82, 0.08,
      0.24, 0.16, 0.60,
      0.55, 0.20, 0.25
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c(
        "CompA::KLF5::Target-Up",
        "CompA::KLF5::Target-Down",
        "CompB::SOX9::Target-Up",
        "CompA::SOX9::Target-Down"
      ),
      c("Topic1", "Topic2", "Topic3")
    )
  )
  theta_cond <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.20, 0.65, 0.15,
      0.25, 0.15, 0.60
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("CondA::KLF5", "CondB::KLF5", "CondB::SOX9"),
      c("Topic1", "Topic2", "Topic3")
    )
  )
  out_cmp <- tempfile(fileext = ".pdf")
  out_cond <- tempfile(fileext = ".pdf")

  cmp_path <- .plot_raw_theta_document_heatmap(
    theta = theta_cmp,
    out_file = out_cmp,
    doc_design = "comparison",
    title_prefix = "comparison aggr MultiVI"
  )
  cond_path <- .plot_raw_theta_document_heatmap(
    theta = theta_cond,
    out_file = out_cond,
    doc_design = "condition",
    title_prefix = "condition aggr LDA"
  )

  expect_identical(cmp_path, out_cmp)
  expect_identical(cond_path, out_cond)
  expect_gt(file.info(out_cmp)$size, 1000)
  expect_gt(file.info(out_cond)$size, 1000)
})

test_that("raw theta document heatmap hides dense row labels by default", {
  expect_false(formals(.plot_raw_theta_document_heatmap)$show_rownames)
})

test_that("document theta UMAP writes condition and selected TF metadata", {
  skip_if_not_installed("uwot")
  configured_colors <- .module3_topic_condition_colors(list(
    report = list(condition_colors = list(CondA = "#112233", CondB = "#aabbcc"))
  ))
  expect_equal(configured_colors, c(CondA = "#112233", CondB = "#AABBCC"))
  theta <- matrix(
    c(
      0.85, 0.10, 0.05,
      0.75, 0.20, 0.05,
      0.10, 0.80, 0.10,
      0.15, 0.75, 0.10,
      0.05, 0.15, 0.80,
      0.10, 0.20, 0.70,
      0.60, 0.10, 0.30,
      0.20, 0.60, 0.20
    ),
    nrow = 8L,
    byrow = TRUE,
    dimnames = list(
      c(
        "CondA::TF1", "CondB::TF1", "CondA::TF2", "CondB::TF2",
        "CondA::TF3", "CondB::TF3", "CondA::TF4", "CondB::TF4"
      ),
      paste0("Topic", 1:3)
    )
  )
  out <- file.path(tempdir(), "document_theta_umap_K3.pdf")
  expect_equal(
    .plot_document_theta_umap(
      theta,
      out_file = out,
      doc_design = "condition",
      selected_tfs = c("TF1", "TF3"),
      seed = 17L,
      n_neighbors = 3L,
      condition_colors = configured_colors
    ),
    out
  )
  expect_true(file.exists(out))
  coords <- data.table::fread(sub("[.]pdf$", ".csv", out))
  selected <- data.table::fread(sub("[.]pdf$", "_selected_tfs.csv", out))
  expect_equal(sort(unique(coords$group_label)), c("CondA", "CondB"))
  expect_equal(sort(coords[selected_tf == TRUE, unique(tf_display)]), c("TF1", "TF3"))
  expect_equal(sort(selected[selected == TRUE, tf_display]), c("TF1", "TF3"))
  expect_true(all(c("UMAP1", "UMAP2", "primary_topic", "primary_theta", "condition_color", "topic_color") %in% names(coords)))
  expect_equal(unique(coords[group_label == "CondA", condition_color]), "#112233")
  expect_equal(unique(coords[group_label == "CondB", condition_color]), "#AABBCC")
  expect_equal(unique(coords[primary_topic == "Topic1", topic_color]), "#E15759")
  qpdf <- Sys.which("qpdf")
  if (nzchar(qpdf)) {
    expect_equal(trimws(system2(qpdf, c("--show-npages", out), stdout = TRUE)), "1")
  }
})

test_that("topic-term assignment writer keeps dense term-topic assignment outputs", {
  score_mat <- matrix(
    c(
      1.0, 0.1, 0.4,
      0.2, 1.0, 0.5
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:G1", "GENE:G2", "PEAK:P1"))
  )
  topic_terms <- data.table::data.table(
    topic = c(1L, 2L),
    term_id = c("GENE:G1", "GENE:G2"),
    score = c(1, 1),
    in_topic = TRUE
  )
  out_csv <- tempfile(fileext = ".csv")

  res <- .write_topic_term_primary_assignment(
    score_mat = score_mat,
    topic_terms = topic_terms,
    assignment_file = out_csv
  )

  expect_identical(res, out_csv)
  expect_true(file.exists(out_csv))
  assignment <- data.table::fread(out_csv)
  expect_true(all(c("term_id", "term_group", "primary_topic", "in_any_topic", "max_score") %in% names(assignment)))
  expect_equal(assignment[term_id == "GENE:G1", primary_topic], "Topic1")
  expect_equal(assignment[term_id == "GENE:G2", primary_topic], "Topic2")
})

test_that("topic-term score heatmap writes derived score diagnostic", {
  phi <- matrix(
    c(
      0.60, 0.20, 0.20,
      0.20, 0.60, 0.20
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:G1", "GENE:G2", "GENE:COMMON"))
  )
  score_mat <- score_terms_normtop(phi, method = "normtop_specificity")
  topic_terms <- data.table::data.table(
    topic = c(1L, 2L),
    term_id = c("GENE:G1", "GENE:G2"),
    score = c(1, 1),
    in_topic = TRUE
  )
  out_file <- tempfile(fileext = ".pdf")

  res <- .plot_topic_term_phi_score_comparison_heatmap(
    phi = phi,
    score_mat = score_mat,
    topic_terms = topic_terms,
    out_file = out_file,
    topic_score_method = "normtop_specificity",
    title_prefix = "topic term test"
  )

  expect_identical(res, out_file)
  expect_gt(file.info(out_file)$size, 1000)
  expect_lt(score_mat["Topic1", "GENE:COMMON"], score_mat["Topic1", "GENE:G1"])
  expect_lt(score_mat["Topic2", "GENE:COMMON"], score_mat["Topic2", "GENE:G2"])
})

test_that("topic extraction standard output skips raw theta documents", {
  theta <- matrix(
    c(
      0.70, 0.20, 0.10,
      0.20, 0.65, 0.15,
      0.25, 0.15, 0.60
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("CondA::KLF5", "CondB::KLF5", "CondB::SOX9"),
      c("Topic1", "Topic2", "Topic3")
    )
  )
  phi <- matrix(
    c(
      0.6, 0.3, 0.1,
      0.2, 0.7, 0.1,
      0.1, 0.2, 0.7
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2", "Topic3"), c("GENE:G1", "GENE:G2", "PEAK:P1"))
  )
  dtm <- Matrix::Matrix(
    matrix(c(4, 1, 0, 1, 4, 0, 0, 1, 4), nrow = 3, byrow = TRUE),
    sparse = TRUE
  )
  rownames(dtm) <- rownames(theta)
  colnames(dtm) <- colnames(phi)
  edges_docs <- data.table::data.table(
    doc_id = rownames(theta),
    tf = c("KLF5", "KLF5", "SOX9"),
    gene_key = c("G1", "G2", "G1"),
    peak_id = c("P1", "P1", "P1")
  )
  topic_base <- list(theta = theta, phi = phi)

  out_enabled <- withr::local_tempdir()
  run_tfdocs_report_from_topic_base(
    topic_base = topic_base,
    dtm = dtm,
    edges_docs = edges_docs,
    out_dir = out_enabled,
    doc_design = "condition",
    in_topic_min_terms = 1L,
    top_n_terms = 3L,
    extraction_steps = c("tf_topic_assignment", "topic_term_heatmap"),
    title_prefix = "condition aggr LDA"
  )
  expect_false(file.exists(file.path(out_enabled, "raw_theta_documents_K3.pdf")))
  expect_false(file.exists(file.path(out_enabled, "topic_term_score_heatmap_K3.pdf")))
  expect_true(file.exists(file.path(out_enabled, "topic_term_phi_score_heatmap_K3.pdf")))
  expect_true(file.exists(file.path(out_enabled, "topic_term_primary_assignment.csv")))
  expect_false(dir.exists(file.path(out_enabled, "tf_topic_assignment")))
  expect_false(dir.exists(file.path(out_enabled, "topic_term_assignment")))
  expect_false(dir.exists(file.path(out_enabled, "doc_topic_heatmaps")))
  expect_false(dir.exists(file.path(out_enabled, "ldavis")))
  expect_false(file.exists(file.path(out_enabled, "tf_topic_membership.csv")))
  expect_true(file.exists(file.path(out_enabled, "tf_topic_membership_pass.csv")))
  expect_true(file.exists(file.path(out_enabled, "tf_topic_primary.csv")))
  expect_true(file.exists(file.path(out_enabled, "tf_direction_topic_summary.csv")))
  expect_false(file.exists(file.path(out_enabled, "tf_topic_assignment_heatmap.pdf")))
  expect_false(file.exists(file.path(out_enabled, "tf_topic_assignment_heatmaps.pdf")))
  expect_false(file.exists(file.path(out_enabled, "tf_primary_topic_dotplot.pdf")))
  expect_false(file.exists(file.path(out_enabled, "tf_direction_topic_map.pdf")))
  expect_false(file.exists(file.path(out_enabled, "tf_topic_assignment_browser.html")))

  out_disabled <- withr::local_tempdir()
  run_tfdocs_report_from_topic_base(
    topic_base = topic_base,
    dtm = dtm,
    edges_docs = edges_docs,
    out_dir = out_disabled,
    doc_design = "condition",
    in_topic_min_terms = 1L,
    top_n_terms = 3L,
    extraction_steps = "tf_topic_assignment",
    title_prefix = "condition aggr LDA"
  )
  expect_false(file.exists(file.path(out_disabled, "raw_theta_documents_K3.pdf")))
})

test_that("topic extraction defaults keep per-comparison pathway outputs flat", {
  expect_true(isTRUE(formals(run_tfdocs_report_from_topic_base)$pathway_per_comparison_flat))
  expect_null(formals(run_tfdocs_report_from_topic_base)$pathway_per_condition)
})

test_that("topic extraction skips removed marker feature outputs", {
  expect_false("topic_marker_heatmap" %in% .topic_extraction_step_names())

  theta <- matrix(
    c(0.80, 0.20, 0.25, 0.75),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("CondA::KLF5", "CondB::SOX9"), c("Topic1", "Topic2"))
  )
  phi <- matrix(
    c(0.70, 0.20, 0.10, 0.10, 0.25, 0.65),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:G1", "GENE:G2", "PEAK:P1"))
  )
  dtm <- Matrix::Matrix(matrix(c(3, 1, 0, 1, 2, 3), nrow = 2, byrow = TRUE), sparse = TRUE)
  rownames(dtm) <- rownames(theta)
  colnames(dtm) <- colnames(phi)
  edges_docs <- data.table::data.table(
    doc_id = rownames(theta),
    tf = c("KLF5", "SOX9"),
    gene_key = c("G1", "G2"),
    peak_id = c("P1", "P1")
  )
  out_dir <- withr::local_tempdir()

  run_tfdocs_report_from_topic_base(
    topic_base = list(theta = theta, phi = phi),
    dtm = dtm,
    edges_docs = edges_docs,
    out_dir = out_dir,
    doc_design = "condition",
    in_topic_min_terms = 1L,
    top_n_terms = 3L,
    extraction_steps = c("topic_terms"),
    title_prefix = "condition aggr LDA"
  )

  expect_false(file.exists(file.path(out_dir, "topic_marker_features.csv")))
  expect_false(file.exists(file.path(out_dir, "topic_marker_term_heatmap.pdf")))
})

test_that("TF topic assignment writer creates reusable CSV outputs without browser artifacts", {
  theta <- matrix(
    c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("CompA::Batf::Target-Up", "CompA::Batf::Target-Down"),
      c("Topic1", "Topic2", "Topic3")
    )
  )
  out_dir <- withr::local_tempdir()

  .write_tf_topic_assignment_outputs(
    theta = theta,
    out_dir = out_dir,
    doc_design = "comparison",
    membership_cutoff = 0.3,
    primary_margin_cutoff = 0.1
  )

  expect_false(dir.exists(file.path(out_dir, "tf_topic_assignment")))
  expect_false(file.exists(file.path(out_dir, "tf_topic_membership.csv")))
  expect_true(file.exists(file.path(out_dir, "tf_topic_membership_pass.csv")))
  expect_true(file.exists(file.path(out_dir, "tf_topic_primary.csv")))
  expect_true(file.exists(file.path(out_dir, "tf_direction_topic_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "tf_topic_tf_list.csv")))
  expect_true(file.exists(file.path(out_dir, "tf_topic_assignment_summary.csv")))
  expect_false(file.exists(file.path(out_dir, "tf_topic_assignment_browser.html")))

  primary <- data.table::fread(file.path(out_dir, "tf_topic_primary.csv"))
  expect_equal(primary[doc_id == "CompA::Batf::Target-Up", primary_topic], "Topic1")
  expect_equal(primary[doc_id == "CompA::Batf::Target-Down", primary_topic], "Topic2")
  tf_list <- data.table::fread(file.path(out_dir, "tf_topic_tf_list.csv"))
  expect_equal(tf_list$tf, "Batf")
  expect_equal(tf_list$tf_display, "Batf")
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
