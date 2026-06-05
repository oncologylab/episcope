# craftgrn 0.1.4

* Added production-oriented Module 3 APIs for regulatory topic modeling:
  `run_regulatory_topics()`, `module3_prepare_topic_inputs()`, and
  `build_module3_qc_report()`.
* Added reusable Module 3 topic-input caches and compact pass-only topic-link
  output as the default production behavior.
* Updated Module 3 benchmark/report discovery to read compact topic-link
  outputs.

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
