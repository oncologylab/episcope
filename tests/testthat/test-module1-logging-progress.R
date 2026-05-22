test_that("progress helpers are quiet when verbose is FALSE", {
  out <- capture.output(
    msg <- capture.output(
      .with_progress(
        total = 3L,
        message = "quiet progress",
        verbose = FALSE,
        code = {
          for (i in seq_len(3L)) .progress_update()
        }
      ),
      type = "message"
    ),
    type = "output"
  )

  expect_equal(out, character())
  expect_equal(msg, character())
})

test_that("progress helpers run user code and return its value", {
  value <- .with_progress(
    total = 2L,
    message = "return progress",
    verbose = FALSE,
    code = {
      .progress_update()
      .progress_update()
      "done"
    }
  )

  expect_equal(value, "done")
})

test_that("progress helpers tolerate extra updates and explicit done", {
  value <- .with_progress(
    total = 1L,
    message = "extra progress",
    verbose = FALSE,
    code = {
      .progress_update()
      .progress_update()
      .progress_done()
      TRUE
    }
  )

  expect_true(value)
})
