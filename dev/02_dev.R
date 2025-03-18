# 1. Dependencies
attachment::att_amend_desc()

# 2. Modules - skipping for now
# golem::add_module(name = "chromhmm_analysis", with_test = TRUE)
# golem::add_module(name = "methylation_analysis", with_test = TRUE)

# 3. Helper functions
golem::add_fct("helpers", with_test = TRUE)
golem::add_utils("helpers", with_test = TRUE)

# 4. (Optional) External resources - skipping for now
## Creates .js and .css files at inst/app/www
# golem::add_js_file("script")
# golem::add_js_handler("handlers")
# golem::add_css_file("custom")
# golem::add_sass_file("custom")
# golem::add_any_file("file.json")


# 5. (Optional) Add internal datasets - skipping for now
# usethis::use_data_raw(name = "my_dataset", open = FALSE)

# 6. Testing - skipping for now
# usethis::use_test("app")

# 7. Documentation & Vignettes
usethis::use_vignette("episcope")
devtools::build_vignettes()

# 8. Code coverage - optional
# usethis::use_coverage()
# covrpage::covrpage()

# 9. CI - skipping if not on Git
# usethis::use_github()
# usethis::use_github_action_check_standard()

rstudioapi::navigateToFile("dev/03_deploy.R")

