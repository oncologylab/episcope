## /!\ Note: if needs to change the name of app during development,
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
##
golem::fill_desc(
  pkg_name = "episcope", # The name of the golem package (lowercase, no underscores or periods)
  pkg_title = "Episcope: An Integrative Epigenomic Analysis Framework", # One-line title describing the package
  pkg_description = "Robust and reproducible framework for epigenetic data analysis. It integrates functionalities for ChromHMM analysis, DNA methylation assessment, differentiation state profiling, paired ATAC-seq and RNA-seq analysis, and functional mapping to facilitate comprehensive epigenomic studies.",
  authors = c(
    person(
      given = "Yaoxiang",
      family = "Li",
      email = "yaoxiang@example.com",
      role = c("aut", "cre")
    ),
    person(
      given = "Chunling",
      family = "Yi",
      email = "chunling@example.com",
      role = c("aut")
    )
  ),
  repo_url = NULL,
  pkg_version = "0.1.1",
  set_options = TRUE
)

## Install the required dev dependencies ----
golem::install_dev_deps()

# Add dependencies ----
usethis::use_package("basilisk")
usethis::use_package("reticulate")



# Core runtime Imports (all are used directly in episcope code)
imports <- c(
  "readr", "dplyr", "tidyr", "tibble",
  "stringr", "purrr", "rlang",
  "igraph", "visNetwork",
  "htmltools", "htmlwidgets", "jsonlite",
  "scales"
)

# Add them to DESCRIPTION: Imports
for (pkg in imports) usethis::use_package(pkg, type = "Imports")

# (Optional)
suggests <- c("testthat", "knitr", "rmarkdown")
for (pkg in suggests) usethis::use_package(pkg, type = "Suggests")


render_tf_hub_delta_network(
  comp_csv  = "Z:/episcope/inst/extdata/lighting/AsPC1_Glc_vs_AsPC1_Ctrl_delta_links_filtered_lda_K20.csv",
  input_tf  = "HNF4A",
  out_html  = "AsPC1_Glc_vs_AsPC1_Ctrl_delta_links_HNF4A.html",
  edge_filter_min = 2,
  edge_filter_on  = "either",
  gene_fc_thresh  = 1.5,
  de_reference    = "str_over_ctrl",
  motif_db = "jaspar2024",
  ring_tf_only_direct = TRUE
)



## Create Common Files ----
## See ?usethis for more information
usethis::use_gpl3_license() # Set GPL-3 License
golem::use_readme_rmd(open = FALSE)
devtools::build_readme()
# Note that `contact` is required since usethis version 2.1.5
# If your {usethis} version is older, you can remove that param
usethis::use_code_of_conduct(contact = "Yaoxiang Li")
usethis::use_lifecycle_badge("Experimental")
usethis::use_news_md(open = FALSE)


## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Favicon ----
# If you want to change the favicon
golem::use_favicon() # path = "path/to/ico". Can be an online file.
# golem::remove_favicon() # Uncomment to remove the default favicon

## Add helper functions ----
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)


# go to dev/02_dev.R
rstudioapi::navigateToFile("dev/02_dev.R")
