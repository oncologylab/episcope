test_that("load_config resolves relative paths against the config directory", {
  project_dir <- tempfile("craftgrn-project-")
  dir.create(file.path(project_dir, "data"), recursive = TRUE)
  config_path <- file.path(project_dir, "project.yaml")
  writeLines(
    c(
      "base_dir: \".\"",
      "ref_genome: hg38",
      "db: JASPAR2024",
      "threshold_expr: 10",
      "threshold_fp_score: 2",
      "threshold_fp_tf_corr_r: 0.3",
      "link_window_bp: 100000",
      "threshold_rna_gene_corr_r: 0.3",
      "threshold_rna_gene_corr_p: null",
      "threshold_rna_gene_corr_fdr: null",
      "threshold_fp_gene_corr_r: 0.3",
      "threshold_fp_gene_corr_p: null",
      "threshold_fp_gene_corr_fdr: null",
      "sample_metadata: data/sample_metadata_strict.csv",
      "atac_master: data/atac_master.csv",
      "rna_mapped: data/rna_strict.csv"
    ),
    config_path
  )

  env <- new.env(parent = emptyenv())
  load_config(config_path, env = env)

  project_dir_n <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  expect_equal(get("base_dir", envir = env), project_dir_n)
  expect_equal(
    get("sample_metadata", envir = env),
    file.path(project_dir_n, "data/sample_metadata_strict.csv")
  )
  expect_equal(
    get("atac_master", envir = env),
    file.path(project_dir_n, "data/atac_master.csv")
  )
  expect_equal(
    get("rna_mapped", envir = env),
    file.path(project_dir_n, "data/rna_strict.csv")
  )
})

test_that("project config validation rejects unknown and malformed entries", {
  base <- list(
    db = "JASPAR2024",
    ref_genome = "hg38",
    threshold_expr = 10,
    threshold_fp_score = 2,
    threshold_fp_tf_corr_r = 0.3,
    link_window_bp = 100000,
    threshold_rna_gene_corr_r = 0.3,
    threshold_fp_gene_corr_r = 0.3
  )
  validate <- function(x) validate_config(
    required = character(),
    numeric_required = character(),
    project_config = x
  )
  expect_true(validate(base))

  expression_filter <- base
  expression_filter$topic_condition_gene_expression_min <- 10
  expect_true(validate(expression_filter))
  expression_filter$topic_condition_gene_expression_min <- -1
  expect_error(validate(expression_filter), "non-negative", fixed = TRUE)

  misspelled <- base
  misspelled$topic_link_prob_cutof <- 0.3
  expect_error(
    validate(misspelled),
    "topic_link_prob_cutoff",
    fixed = TRUE
  )

  malformed <- base
  malformed$topic_link_prob_cutoff <- 1.2
  expect_error(
    validate(malformed),
    "between 0 and 1",
    fixed = TRUE
  )
})

test_that("project config validation checks nested report state", {
  cfg <- list(
    report = list(
      defaults = list(
        condition_1_suffix = "_0_BCAA_Ctrl",
        condition_2_suffix = "_10_FBS_Ctrl",
        topic = 1L,
        trained_k = 30L,
        topic_space = "combined",
        short_condition_names = TRUE
      ),
      condition_order = c("HPAFII_0_BCAA_Ctrl", "HPAFII_10_FBS_Ctrl"),
      condition_colors = list(
        HPAFII_0_BCAA_Ctrl = "#3C5488",
        HPAFII_10_FBS_Ctrl = "#FF9DA7"
      )
    )
  )
  expect_true(validate_config(
    required = character(),
    numeric_required = character(),
    project_config = cfg
  ))

  bad_default <- cfg
  bad_default$report$defaults$condition_one <- "A"
  expect_error(
    validate_config(
      required = character(),
      numeric_required = character(),
      project_config = bad_default
    ),
    "Unknown `report.defaults` entries",
    fixed = TRUE
  )

  bad_color <- cfg
  bad_color$report$condition_colors[[1L]] <- "red"
  expect_error(
    validate_config(
      required = character(),
      numeric_required = character(),
      project_config = bad_color
    ),
    "six-digit hex",
    fixed = TRUE
  )

  bad_topic_space <- cfg
  bad_topic_space$report$defaults$topic_space <- "optimized"
  expect_error(
    validate_config(
      required = character(),
      numeric_required = character(),
      project_config = bad_topic_space
    ),
    "raw or combined",
    fixed = TRUE
  )

  expect_true(validate_config(
    required = character(),
    numeric_required = character(),
    project_config = list(
      report = list(
        defaults = list(trained_k = 30L, topic_space = "combined")
      )
    )
  ))
})

test_that("project resource and batch-correction configuration is validated", {
  cfg <- list(
    resources = list(max_memory_fraction = 0.8, reserve_memory_gb = 32, max_workers = 12L),
    batch_correction = list(
      enabled = TRUE,
      method = "ruvr",
      batch_column = "study_id",
      preserve_column = "condition_id",
      k_candidates = c(0L, 1L, 2L, 3L, 5L, 10L),
      effect_spearman_min = 0.95
    )
  )
  expect_true(validate_config(required = character(), numeric_required = character(), project_config = cfg))
  cfg$resources$max_memory_fraction <- 0.95
  expect_error(
    validate_config(required = character(), numeric_required = character(), project_config = cfg),
    "max_memory_fraction"
  )
})

test_that("nested Module 2 and Module 3 configuration is validated", {
  cfg <- list(
    module2 = list(
      max_distance_bp = 100000,
      threshold_fp_gene_corr_r = 0.3,
      threshold_tf_target_corr_r = 0.4,
      threshold_fp_target_corr_fdr = 0.05
    ),
    module3 = list(
      topic_benchmark_methods = c("condition_aggr_lda", "condition_aggr_multivi"),
      topic_k_grid = c(10L, 20L, 30L),
      topic_term_assignment_method = "gammafit_maxprob",
      topic_gammafit_thrP = list(lda = 0.8, multivi = 0.85),
      topic_condition_gene_weighting = "specificity",
      topic_condition_peak_weighting = "tf_expression",
      topic_condition_specificity_temperature = 0.5,
      topic_condition_specificity_floor = 0.1,
      topic_condition_specificity_expression_min = 10,
      topic_vae_device = "cpu",
      pathway_species = "human"
    )
  )
  validate <- function(x) validate_config(
    required = character(),
    numeric_required = character(),
    project_config = x
  )

  expect_true(validate(cfg))

  bad_module2 <- cfg
  bad_module2$module2$max_distnce_bp <- 1000
  expect_error(validate(bad_module2), "Unknown `module2` entries", fixed = TRUE)

  bad_module3 <- cfg
  bad_module3$module3$topic_vae_device <- "gpu"
  expect_error(validate(bad_module3), "module3_topic_vae_device", fixed = TRUE)

  bad_assignment <- cfg
  bad_assignment$module3$topic_term_assignment_method <- "maximum"
  expect_error(validate(bad_assignment), "module3_topic_term_assignment_method", fixed = TRUE)

  bad_gammafit <- cfg
  bad_gammafit$module3$topic_gammafit_thrP$multivi <- 1.2
  expect_error(validate(bad_gammafit), "GammaFit probability", fixed = TRUE)

  bad_specificity <- cfg
  bad_specificity$module3$topic_condition_specificity_floor <- 1
  expect_error(validate(bad_specificity), "specificity floors", fixed = TRUE)

  bad_peak_weighting <- cfg
  bad_peak_weighting$module3$topic_condition_peak_weighting <- "raw_expression"
  expect_error(
    validate(bad_peak_weighting),
    "module3_topic_condition_peak_weighting",
    fixed = TRUE
  )

  unknown_module3 <- cfg
  unknown_module3$module3$pathway_speces <- "human"
  expect_error(validate(unknown_module3), "Unknown `module3` entries", fixed = TRUE)
})

test_that("legacy run labels do not partially match nested config blocks", {
  cfg <- list(
    module2_run_label = "m1R0p3_100kb",
    threshold_tf_target_corr_r = 0.4,
    threshold_fp_target_corr_fdr = 0.05,
    module3_label_col = "condition_id",
    report_state = list(defaults = list(short_condition_names = TRUE))
  )
  expect_true(validate_config(
    required = character(),
    numeric_required = character(),
    project_config = cfg
  ))
  cfg$threshold_tf_target_corr_r <- 1.2
  expect_error(
    validate_config(
      required = character(),
      numeric_required = character(),
      project_config = cfg
    ),
    "Correlation cutoffs"
  )
})
