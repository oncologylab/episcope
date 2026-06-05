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

The previous incoming pretest for 0.1.2 failed on Windows because two tests
compared temporary paths with platform-specific separators, and temporary topic
model fixture directories were not robustly created before CSV writes. This
resubmission normalizes path comparisons and explicitly creates the fixture
directories before writing test files.

This version also completes the package-native Module 3 regulatory topic
workflow, including reusable topic input caches, pass-only compact topic-link
outputs, and a Module 3 QC report. The package metadata wording was also revised
to reduce incoming spell-check notes for abbreviated scientific terms.

## Package data

A large local GeneHancer CSV file used in earlier development is excluded from
the package source with `.Rbuildignore`. Current Module 2 workflows accept
regulatory priors as user-supplied inputs rather than bundling the full
GeneHancer database.
