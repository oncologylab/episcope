testthat::skip_if_not_installed("shiny")
testthat::skip_if_not_installed("golem")

test_that("app_ui returns a Shiny tag list and keeps request argument", {
  ui <- app_ui()
  golem::expect_shinytaglist(ui)
  # Check that formals have not been removed
  fmls <- formals(app_ui)
  for (i in c("request")) {
    expect_true(i %in% names(fmls))
  }
})

test_that("app_server keeps the standard Shiny server signature", {
  server <- app_server
  expect_type(server, "closure")
  # Check that formals have not been removed
  fmls <- formals(app_server)
  for (i in c("input", "output", "session")) {
    expect_true(i %in% names(fmls))
  }
})

test_that(
  "app_sys resolves packaged app configuration files",
  {
    expect_true(
      app_sys("golem-config.yml") != ""
    )
  }
)

test_that(
  "golem config reads production and development flags",
  {
    config_file <- app_sys("golem-config.yml")
    skip_if(config_file == "")

    expect_true(
      get_golem_config(
        "app_prod",
        config = "production",
        file = config_file
      )
    )
    expect_false(
      get_golem_config(
        "app_prod",
        config = "dev",
        file = config_file
      )
    )
  }
)

test_that(
  "Shiny application starts without runtime errors",
  {
    golem::expect_running(sleep = 5)
  }
)
