## Test environments

* Local Ubuntu 24.04.4 LTS, R 4.5.3

## R CMD check results

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-mno-omit-leaf-frame-pointer'

  This flag is injected by the local compiler/toolchain environment.

## Resubmission

This is a resubmission with the version bumped to 0.1.4.

The previous incoming pretests failed on Windows test fixtures only. Version
0.1.2 had path-separator-sensitive tests and temporary topic model fixture
directories that were not robustly created before CSV writes. Version 0.1.3
fixed those issues, but one Module 3 test fixture still used the historical
deep benchmark folder layout, which exceeded the practical Windows temporary
path length during CRAN incoming checks. This resubmission uses the current
short Module 3 output layout for synthetic test fixtures while preserving
legacy benchmark path discovery in package code.

This version also completes the package-native Module 3 regulatory topic
workflow, including reusable topic input caches, pass-only compact topic-link
outputs, and a Module 3 QC report. The package metadata wording was also revised
to reduce incoming spell-check notes for abbreviated scientific terms.

## Package data

A large local GeneHancer CSV file used in earlier development is excluded from
the package source with `.Rbuildignore`. Current Module 2 workflows accept
regulatory priors as user-supplied inputs rather than bundling the full
GeneHancer database.
