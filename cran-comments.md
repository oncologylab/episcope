## Test environments

* Local Ubuntu 24.04.4 LTS, R 4.5.3

## R CMD check results

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-mno-omit-leaf-frame-pointer'

  This flag is injected by the local compiler/toolchain environment.

## First submission

This is the first CRAN submission of craftgrn.

## Package data

A large local GeneHancer CSV file used in earlier development is excluded from
the package source with `.Rbuildignore`. Current Module 2 workflows accept
regulatory priors as user-supplied inputs rather than bundling the full
GeneHancer database.
