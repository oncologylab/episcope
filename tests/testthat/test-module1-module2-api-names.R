test_that("Module 1 and Module 2 expose publication-facing function names", {
  exports <- getNamespaceExports("craftgrn")

  expected <- c(
    "module1_correlate_TF_to_canonical_tfbs",
    "module1_filter_canonical_bound_tfbs",
    "module1_predict_full_tfbs",
    "output_predicted_tfbs",
    "predict_tf_targets",
    "module2_identify_candidate_links",
    "module2_link_fp_targets",
    "module2_output_predicted_links",
    "load_predicted_links",
    "query_predicted_links",
    "check_predicted_links",
    "report_top_tf_targets",
    "report_direct_tf_tf_regulations",
    "report_tf_tf_coregulations"
  )
  retired <- c(
    "module1_correlate_motif_supported_tfbs",
    "module1_select_prediction_footprints",
    "module1_predict_tfbs_from_correlations",
    "build_predicted_tfbs",
    "link_tf_targets",
    "module2_prepare_link_inputs",
    "module2_build_fp_target_candidates",
    "module2_assemble_links",
    "output_predicted_links",
    "load_module2_links",
    "query_module2_links",
    "validate_module2_links",
    "check_module2_links",
    "build_module2_reports",
    "export_top_tf_targets",
    "export_direct_tf_tf_browser",
    "export_tf_tf_connectivity_browser"
  )

  expect_true(all(expected %in% exports))
  expect_false(any(retired %in% exports))
})
