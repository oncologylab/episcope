## Test environments

* Local Ubuntu 24.04.4 LTS, R 4.5.3

## R CMD check results

0 errors | 0 warnings | 3 notes

* checking CRAN incoming feasibility ... NOTE
  New submission

  This is the first CRAN submission for craftgrn.

* checking for future file timestamps ... NOTE
  unable to verify current time

  This appears to be a local check environment clock-verification issue.

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-mno-omit-leaf-frame-pointer'

  This flag is injected by the local compiler/toolchain environment.

## Resubmission

This is a resubmission with the version bumped to 0.1.7, addressing CRAN
feedback from Benjamin Altmann.

In this resubmission:

* Added method references to the DESCRIPTION field using CRAN auto-linking
  formats for DOI and HTTPS references. The complete CraftGRN manuscript is in
  preparation.
* Revised DESCRIPTION references after incoming pretest feedback to use the
  arXiv DOI form, include author-year reference text, and avoid incoming
  spell-check notes.
* Added missing value documentation for exported functions, including
  `run_app()`. The current public Module 3 step pages
  `module3_train_topic_models.Rd` and `module3_extract_topics.Rd` include
  value sections; old internal aliases `train_topic_models` and
  `extract_regulatory_topics` are no longer exported or documented.
* Replaced `\dontrun{}` examples with `\donttest{}` and made the examples
  self-contained or guarded when they require local project files.
* Removed default writes to the working directory from Module 1 public writing
  functions. `predict_tfbs()` and `module1_predict_full_tfbs()` now default to
  in-memory operation and require an explicit output directory when
  `write_outputs = TRUE`. README and vignette examples now write under
  `tempdir()`.
* Removed the fixed RNG seed from the Module 1 QC sampling helper by replacing
  random subsampling with deterministic evenly spaced row selection.
* Added regression coverage for the no-default-write Module 1 behavior.

## Package data

A large local GeneHancer CSV file used in earlier development is excluded from
the package source with `.Rbuildignore`. Current Module 2 workflows accept
regulatory priors as user-supplied inputs rather than bundling the full
GeneHancer database.
