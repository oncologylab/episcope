test_that("TFBS condition comparison uses active sites and summed FP scores", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  result <- predict_tfbs(compact, r_cutoff = 0.8, verbose = FALSE)

  plot <- plot_tfbs_condition_comparison(
    result,
    cond1 = "cond_a",
    cond2 = "cond_d",
    pseudocount = 1,
    verbose = FALSE
  )
  summary <- attr(plot, "tfbs_summary")
  tf_a <- summary[summary$tf == "TF_A", , drop = FALSE]

  expect_s3_class(plot, "ggplot")
  expect_equal(tf_a$n_tfbs_cond1, 1L)
  expect_equal(tf_a$n_tfbs_cond2, 1L)
  expect_equal(tf_a$delta_n_predicted_tfbs, 0L)
  expect_equal(tf_a$log2fc_fp_score, log2((1 + 1) / (9 + 1)))
  expect_equal(tf_a$max_tf_expression, 9)
  expect_equal(tf_a$log2fc_tf_expression, log2((1 + 1) / (9 + 1)))
  expect_equal(plot$labels$title, "cond_a v cond_d")
})

test_that("TFBS condition comparison aggregates FP scores over the active-site union", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  fp <- rownames(compact$matrices$fp_score)[1:2]
  compact$matrices$fp_bound[fp, ] <- FALSE
  compact$matrices$fp_bound[fp[[1L]], "cond_a"] <- TRUE
  compact$matrices$fp_bound[fp[[2L]], "cond_d"] <- TRUE
  compact$matrices$fp_score[fp, "cond_a"] <- c(2, 100)
  compact$matrices$fp_score[fp, "cond_d"] <- c(50, 4)
  tf_i <- match("TF_A", rownames(compact$matrices$gene_on))
  compact$matrices$gene_on[tf_i, c("cond_a", "cond_d")] <- TRUE

  plot <- plot_tfbs_condition_comparison(
    tibble::tibble(tf = "TF_A", fp_id = fp),
    multiomic_data = compact,
    cond1 = "cond_a",
    cond2 = "cond_d",
    pseudocount = 1,
    verbose = FALSE
  )
  summary <- attr(plot, "tfbs_summary")

  expect_equal(summary$fp_score_cond1, 102)
  expect_equal(summary$fp_score_cond2, 54)
  expect_equal(summary$log2fc_fp_score, log2(103 / 55))
})

test_that("TF-TF co-binding heatmap supports counts and JSD", {
  predicted <- tibble::tibble(
    tf = c("TF1", "TF1", "TF2", "TF2", "TF3"),
    fp_id = c("fp1", "fp2", "fp2", "fp3", "fp4")
  )

  count_plot <- plot_tf_tf_cobinding_heatmap(predicted, metric = "absolute", verbose = FALSE)
  count_matrix <- attr(count_plot, "cobinding_matrix")
  jsd_plot <- plot_tf_tf_cobinding_heatmap(predicted, metric = "jsd", verbose = FALSE)
  jsd <- attr(jsd_plot, "cobinding_matrix")

  expect_s3_class(count_plot, "ggplot")
  expect_equal(count_matrix["TF1", "TF2"], 1)
  expect_equal(diag(count_matrix), c(TF1 = 2, TF2 = 2, TF3 = 1))
  expect_equal(jsd, t(jsd))
  expect_equal(unname(diag(jsd)), rep(0, 3))
  expect_true(all(jsd >= 0 & jsd <= 1))

  incidence <- Matrix::Matrix(rbind(c(1, 1, 0, 0), c(0, 1, 1, 0), c(0, 0, 0, 1)), sparse = TRUE)
  rownames(incidence) <- c("A", "B", "C")
  expect_equal(craftgrn:::.module1_jsd_matrix(incidence)["A", "B"], sqrt(0.5), tolerance = 1e-12)
})

test_that("Module 1 Run Summary does not report zero TFs when scanning is skipped", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  cards <- craftgrn:::.module1_qc_run_cards(
    tibble::tibble(metric = character(), value = numeric()),
    compact,
    list(tf_summary_all = tibble::tibble()),
    predicted_rows = 10
  )

  expect_equal(cards$value[cards$label == "TFs with predicted binding"], "Not recorded")
  expect_equal(cards$value[cards$label == "Expressed TFs"], "5")
  expect_equal(cards$value[cards$label == "Raw ATAC peaks"], "Not provided")
  expect_equal(cards$label, c(
    "Raw ATAC peaks", "Raw footprints", "Aligned footprints",
    "Filtered footprints", "Predicted unique TFBS",
    "Expressed TFs", "TFs with predicted binding"
  ))
})

test_that("TFBS UMAP report provides interactive cluster choices", {
  skip_if_not_installed("uwot")
  predicted <- tibble::tibble(
    tf = rep(paste0("TF", 1:6), each = 3),
    fp_id = c(
      "p1", "p2", "p3", "p1", "p2", "p4", "p2", "p3", "p5",
      "p3", "p4", "p6", "p4", "p5", "p7", "p5", "p6", "p8"
    )
  )
  out <- tempfile(fileext = ".html")

  path <- build_tfbs_umap_report(
    predicted,
    output_file = out,
    top_variable_tfbs = 8L,
    cluster_range = 2:4,
    default_clusters = 3L,
    seed = 11L,
    verbose = FALSE
  )
  page <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(file.exists(path))
  expect_match(page, "Number of clusters", fixed = TRUE)
  expect_match(page, "option value=\"3\" selected", fixed = TRUE)
  expect_match(page, "TF UMAP from highly variable predicted TFBS", fixed = TRUE)
  expect_match(page, "predicted TFBS", fixed = TRUE)
})

test_that("TFBS Explorer writes compact exact co-binding controls", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  fp <- rownames(compact$matrices$fp_score)[1:4]
  predicted <- tibble::tibble(
    tf = c("TF_A", "TF_A", "TF_B", "TF_B", "TF_C", "TF_C"),
    fp_id = c(fp[[1L]], fp[[2L]], fp[[2L]], fp[[3L]], fp[[2L]], fp[[4L]])
  )
  out <- tempfile(fileext = ".html")

  path <- build_module1_tfbs_explorer(
    predicted,
    multiomic_data = compact,
    output_file = out,
    default_condition = "cond_a",
    verbose = FALSE
  )
  page <- paste(readLines(path, warn = FALSE), collapse = "\n")
  cache <- readRDS(file.path(dirname(out), "module1_qc_analysis.rds"))
  bits <- cache$tf_bits
  shared <- bitwAnd(as.integer(bits$TF_A), as.integer(bits$TF_B))
  shared_with_c <- bitwAnd(shared, as.integer(bits$TF_C))
  popcount <- vapply(0:255, function(value) sum(as.integer(intToBits(value))[1:8]), integer(1L))

  expect_true(file.exists(path))
  expect_match(page, "Binding", fixed = TRUE)
  expect_match(page, "Co-binding", fixed = TRUE)
  expect_match(page, "Motif", fixed = TRUE)
  expect_match(page, "multiple focal TFs", fixed = TRUE)
  expect_match(page, "Export SVG", fixed = TRUE)
  expect_match(page, "focal-search", fixed = TRUE)
  expect_match(page, "motif-search", fixed = TRUE)
  expect_match(page, "DecompressionStream", fixed = TRUE)
  expect_match(page, "new Worker", fixed = TRUE)
  expect_match(page, "Exact co-binding data ready.';applyCobinding()", fixed = TRUE)
  expect_false(grepl("D.tfBits=D.tfBits.map", page, fixed = TRUE))
  expect_equal(sum(popcount[shared_with_c + 1L]), 1L)
  expect_equal(cache$fingerprint$schema, "module1_qc_analysis_v8")
  expect_equal(cache$tf_counts["TF_A", "Overall"], 2)
  expect_equal(sum(cache$tf_counts[, "Overall"]), nrow(unique(predicted[c("tf", "fp_id")])))
  expect_true(sum(cache$tf_counts[, cache$conditions, drop = FALSE]) != sum(cache$tf_counts[, "Overall"]))
})

test_that("Module 1 PCA maps assay IDs and preserves transformed RNA", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  compact$samples <- tibble::tibble(
    id = paste0("raw_", 1:4),
    condition_id = paste0("cond_", letters[1:4]),
    condition_base = c("group_1", "group_1", "group_2", "group_2")
  )
  annotation <- craftgrn:::.module1_qc_sample_annotations(compact, c("raw_1", "cond_b"))
  expect_equal(annotation$display_label, c("cond_a", "cond_b"))
  expect_equal(annotation$biological_group, c("group_1", "group_1"))
  expect_true(all(annotation$matched))

  normalized <- matrix(c(-1, 0, 1, 2, 2, 3, 4, 5), nrow = 2)
  transformed <- craftgrn:::.module1_qc_transform_matrix(normalized, assay = "rna")
  expect_equal(transformed$transformation, "identity")
  expect_equal(transformed$matrix, normalized)

  counts <- matrix(c(0, 100, 1000, 10000), nrow = 2)
  transformed_counts <- craftgrn:::.module1_qc_transform_matrix(counts, assay = "rna")
  expect_equal(transformed_counts$transformation, "log2(x + 1)")
  expect_equal(transformed_counts$matrix, log2(counts + 1))
})

test_that("Motif explorer excludes motifs without retained canonical TF binding", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  fp <- rownames(compact$matrices$fp_score)[[1L]]
  compact$features$fp_motif <- tibble::tibble(
    fp_id = c(fp, fp),
    motif = c("REFERENCE_MOTIF", "ABSENT_MOTIF"),
    tf = c("ZZZ", "ABSENT")
  )
  predicted <- tibble::tibble(tf = c(sprintf("TF%02d", 1:24), "ZZZ"), fp_id = fp)
  out <- tempfile(fileext = ".html")

  build_module1_tfbs_explorer(predicted, multiomic_data = compact, output_file = out, verbose = FALSE)
  page <- paste(readLines(out, warn = FALSE), collapse = "\n")
  cache <- readRDS(file.path(dirname(out), "module1_qc_analysis.rds"))
  reference <- cache$motif_counts[cache$motif_counts$motif == "REFERENCE_MOTIF", ]
  absent <- cache$motif_counts[cache$motif_counts$motif == "ABSENT_MOTIF", ]

  expect_equal(sum(reference$top20), 20L)
  expect_equal(reference$rank[reference$tf == "ZZZ"], 25L)
  expect_equal(reference$canonical_status[reference$tf == "ZZZ"], "predicted_outside_top20")
  expect_equal(nrow(absent), 0L)
  expect_equal(cache$motif_summary, list(total = 2L, eligible = 1L, excluded = 1L))
  expect_false(grepl("not predicted", page, fixed = TRUE))
  expect_match(page, "excluded without retained canonical TF binding", fixed = TRUE)
})

test_that("Explorer cutoff prefers saved Module 1 run metadata", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  compact$project$config <- list(threshold_fp_tf_corr_r = 0.5)
  module1_dir <- tempfile("module1-cutoff-")
  dir.create(module1_dir)
  readr::write_csv(
    tibble::tibble(
      parameter = "r_cutoff",
      yaml_value = "0.3",
      command_value = "",
      effective_value = "0.3",
      source = "YAML"
    ),
    file.path(module1_dir, "module1_run_parameters.csv")
  )

  expect_equal(craftgrn:::.module1_explorer_r_cutoff(module1_dir, omics = compact), 0.3)
})

test_that("Module 1 config provenance keeps relevant values", {
  cfg <- craftgrn:::.module1_relevant_config(list(
    db = "JASPAR2024",
    ref_genome = "hg38",
    threshold_fp_tf_corr_r = 0.4,
    module3_backend = "lda",
    api_key = "do-not-keep"
  ))

  expect_equal(cfg$db, "JASPAR2024")
  expect_equal(cfg$threshold_fp_tf_corr_r, 0.4)
  expect_false("module3_backend" %in% names(cfg))
  expect_false("api_key" %in% names(cfg))
})
