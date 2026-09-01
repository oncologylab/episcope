# craftgrn 0.1.7

* Added Module 1 condition-comparison, TF-TF co-binding, and interactive TFBS
  UMAP visualization utilities, plus a provenance-aware QC Run Summary.
* Redesigned the Module 1 QC report to match the report plan, with exact run
  metrics, full-width raincloud distributions, metadata-aware labeled PCA,
  threshold-colored correlations, and responsive Binding, Co-binding, and
  Motif explorers embedded in one HTML file. Motif rankings now retain the
  true top 20 and report canonical TF references with their actual rank or
  not-predicted status.
* Refined the Module 1 QC visual hierarchy with focused distribution, PCA, and
  count review modes, cleaner violin summaries, full-width PCA panels, larger
  correlation plots, compact guidance, and an automatically populated
  co-binding view.
* Made the Module 1 ATAC master input optional, retained raw-footprint
  provenance for regenerated QC reports, and added responsive 1080p and 4K
  report layouts with ATAC-aware diagnostics.
* Corrected Module 1 correlation QC to persist all-evaluated Pearson, Spearman,
  and best-R histograms with the actual run cutoffs, including negative
  correlations and an explicit legacy recomputation path.
* Redesigned the Module 2 QC report around the Module 1 review hierarchy with
  sticky navigation, accessible tabbed diagnostics, compact condition and
  correlation views, clearer candidate and final-link evidence, and a
  collapsed technical appendix that preserves all detailed checks.
* Rebuilt the Module 2 top-TF report as a condition-aware regulon browser with
  full or selected target sets, optional supporting TFs, correlation controls,
  condition activity and expression views, and standalone SVG export.
* Updated Module 3 topic optimization to merge small, similar topics while
  preserving Gene/Peak assignments and improving TF-term versus document-theta
  correspondence. Merged-topic theta QC uses the maximum member-topic value.
* Added condition-specific Gene weighting, optional condition term-IDF and
  final Peak/Gene token balancing, configurable WarpLDA alpha mass, and new VAE
  initialization, validation, regularization, and paired-term controls.
* Standardized Module 3 input-design QC as the default topic-assignment report,
  with Gene and all-Peak UMAPs when available, condition-composition review,
  optional Module 2 TF-target evidence, and full-universe retention funnels.
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
