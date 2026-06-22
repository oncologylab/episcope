# craftgrn 0.1.7

* Updated DESCRIPTION references to include author-year text for the next CRAN
  update.
* Fixed a CRAN macOS test failure by comparing config-resolved relative paths
  against the canonical project directory path.

# craftgrn 0.1.6

* Revised DESCRIPTION references to avoid CRAN incoming spell-check notes and
  replaced the arXiv URL reference with the requested arXiv DOI format.

# craftgrn 0.1.5

* Addressed CRAN resubmission feedback by adding method references to
  DESCRIPTION, completing value documentation, replacing `\dontrun{}` examples,
  avoiding default writes to the working directory, and removing a fixed seed
  from package function code.

# craftgrn 0.1.4

* Added production-oriented Module 3 APIs for regulatory topic modeling:
  `run_topic_modeling()`, `module3_construct_docs()`,
  `module3_train_topic_models()`, `module3_extract_topics()`, and
  `build_module3_qc_report()`.
* Added reusable Module 3 topic-input caches and compact pass-only topic-link
  output as the default production behavior.
* Added `visualize_topic_modeling_results()` and
  `visualize_differential_grns()` for standard Module 3 result browsers.

# craftgrn 0.1.3

* Fixed Windows incoming pretest failures by normalizing path comparisons in
  report tests and hardening temporary topic-model fixture directory creation.
* Revised CRAN-facing package metadata wording to reduce incoming
  spell-check notes for abbreviated scientific terms.

# craftgrn 0.1.2

* Initial CRAN submission candidate.
* Added Module 1 TFBS prediction, Module 2 TF-target prediction, and Module 3
  regulatory topic modeling workflows.
* Added package-native QC reports and pkgdown documentation.
