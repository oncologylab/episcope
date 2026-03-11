combined_theta_review_script <- normalizePath(file.path(
  Sys.getenv(
    "CRAFTGRN_REPO_DIR",
    unset = normalizePath(file.path(getwd(), "../.."), mustWork = FALSE)
  ),
  "dev/benchmark/03_build_topic_combined_theta_review_html.R"
), mustWork = FALSE)

test_that("combined theta review HTML has one shared K control and two plot regions", {
  testthat::skip_if_not(file.exists(combined_theta_review_script), message = "dev benchmark script is not included in built package")
  old_wd <- setwd(dirname(dirname(dirname(combined_theta_review_script))))
  on.exit(setwd(old_wd), add = TRUE)
  source(combined_theta_review_script, local = TRUE)

  phi_summary <- data.table::data.table(
    k = c(10L, 20L),
    selected_k = c(10L, 20L),
    setup_label = c("cond fp uniq", "diff fp uniq"),
    doc_design = c("condition", "comparison"),
    model_label = c("MultiVI", "LDA"),
    pearson_r = c(0.7, 0.6),
    spearman_r = c(0.4, 0.3),
    panel_label = c("cond fp uniq\nMultiVI", "diff fp uniq\nLDA")
  )
  mds_points <- data.table::data.table(
    group_label = c("HPAFII_10_FBS", "HPAFII_0.05_Glc"),
    display_label = c("10 FBS", "0.05 Glc"),
    color = c("#717171", "#5A8EBC"),
    MDS1 = c(0, 1),
    MDS2 = c(0, 1),
    k = c(10L, 10L),
    panel_label = c("cond fp uniq\nMultiVI", "cond fp uniq\nMultiVI")
  )
  phi_pairs <- data.table::data.table(
    k = rep(10L, 2L),
    selected_k = rep(10L, 2L),
    setup_label = rep("cond fp uniq", 2L),
    model_label = rep("MultiVI", 2L),
    theta_jsd = c(0.1, 0.8),
    phi_jsd = c(0.2, 0.7),
    panel_label = rep("cond fp uniq\nMultiVI", 2L)
  )

  html <- paste(
    .combined_theta_review_html_document(
      phi_summary = phi_summary,
      mds_points = mds_points,
      panels = "cond fp uniq\nMultiVI",
      phi_pairs = phi_pairs
    ),
    collapse = "\n"
  )

  expect_match(html, "id=\"kSlider\"", fixed = TRUE)
  expect_equal(length(gregexpr("id=\\\"kSlider\\\"", html)[[1L]]), 1L)
  expect_match(html, "correlationImg", fixed = TRUE)
  expect_match(html, "groupMdsImg", fixed = TRUE)
  expect_match(html, "IMAGE_MANIFEST", fixed = TRUE)
  expect_match(html, "theta_phi_topic_distance_correlation_k10.png", fixed = TRUE)
  expect_match(html, "theta_group_mds_k10.png", fixed = TRUE)
  expect_match(html, "grid-template-columns:minmax(0,1.3fr) minmax(0,3.7fr)", fixed = TRUE)
  expect_false(grepl("COR_SUMMARY", html, fixed = TRUE))
  expect_false(grepl("COR_PAIRS", html, fixed = TRUE))
  expect_false(grepl("renderCorrelation", html, fixed = TRUE))
  expect_false(grepl("renderMds", html, fixed = TRUE))
  expect_false(grepl("repelLabels", html, fixed = TRUE))
  expect_false(grepl("mdsPanelSelect", html, fixed = TRUE))
  expect_false(grepl("thetaPhiImage", html, fixed = TRUE))
})

test_that("theta group MDS labels and shapes use requested condition/comparison encoding", {
  testthat::skip_if_not(file.exists(combined_theta_review_script), message = "dev benchmark script is not included in built package")
  old_wd <- setwd(dirname(dirname(dirname(combined_theta_review_script))))
  on.exit(setwd(old_wd), add = TRUE)
  source(combined_theta_review_script, local = TRUE)

  expect_equal(.short_group_label("10 FBS"), "Full")
  expect_equal(.short_group_label("0 BCAA"), "0BCAA")
  expect_equal(.short_group_label("25 BCAA"), "25BCAA")
  expect_equal(
    .theta_mds_color_value("#DB7270", "0 Lys", "HPAFII_0_Lys"),
    "#D81B60"
  )
  expect_equal(
    .theta_mds_color_value("#B45D5C", "0 Lys Down", "HPAFII_0_Lys_vs_HPAFII_10_FBS::Target-Down"),
    "#AD1457"
  )
  expect_equal(
    .split_theta_mds_panel(
      panel_label = "diff fp uniq\nLDA",
      doc_design = "comparison",
      display_label = "0 BCAA Up",
      group_label = "HPAFII_0_BCAA_vs_HPAFII_10_FBS::Target-Up"
    ),
    "diff fp uniq\nLDA\nUp"
  )
  expect_equal(
    .split_theta_mds_panel(
      panel_label = "diff fp uniq\nLDA",
      doc_design = "comparison",
      display_label = "0 BCAA Down",
      group_label = "HPAFII_0_BCAA_vs_HPAFII_10_FBS::Target-Down"
    ),
    "diff fp uniq\nLDA\nDown"
  )

  dt <- data.table::data.table(
    group_label = c(
      "HPAFII_0_BCAA_vs_HPAFII_10_FBS::Target-Up",
      "HPAFII_0_BCAA_vs_HPAFII_10_FBS::Target-Down",
      "HPAFII_10_FBS"
    ),
    display_label = c("0 BCAA Up", "0 BCAA Down", "10 FBS"),
    color = c("#111111", "#222222", "#333333"),
    MDS1 = c(0, 1, 2),
    MDS2 = c(0, 1, 2),
    k = c(10L, 10L, 10L),
    doc_design = c("comparison", "comparison", "condition"),
    panel_label = c("diff fp uniq\nLDA", "diff fp uniq\nLDA", "cond fp uniq\nMultiVI")
  )
  tmp <- tempfile(fileext = ".png")
  expect_true(.plot_theta_group_mds_png(
    mds_points = dt,
    panels = unique(dt$panel_label),
    out_file = tmp,
    k = 10L
  ))
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})
