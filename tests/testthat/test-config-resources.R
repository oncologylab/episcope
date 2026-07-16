test_that("safe worker planning respects cores and current memory", {
  plan <- .safe_worker_plan(
    requested = 12L,
    estimated_bytes_per_worker = 10 * 1024^3,
    max_memory_fraction = 0.8,
    reserve_memory_gb = 20,
    available_bytes = 100 * 1024^3,
    physical_cores = 8L
  )
  expect_equal(plan$workers, 8L)
  expect_true(plan$allowed)

  constrained <- .safe_worker_plan(
    requested = 12L,
    estimated_bytes_per_worker = 30 * 1024^3,
    max_memory_fraction = 0.8,
    reserve_memory_gb = 20,
    available_bytes = 100 * 1024^3,
    physical_cores = 8L
  )
  expect_equal(constrained$workers, 2L)
})

test_that("safe worker planning refuses an unaffordable worker", {
  plan <- .safe_worker_plan(
    requested = 4L,
    estimated_bytes_per_worker = 50 * 1024^3,
    available_bytes = 40 * 1024^3,
    reserve_memory_gb = 16,
    physical_cores = 4L
  )
  expect_equal(plan$workers, 0L)
  expect_false(plan$allowed)
})

test_that("resource preflight enforces the live safe budget", {
  expect_true(.resource_preflight(
    estimated_bytes = 20,
    stage = "test",
    available_bytes = 100,
    max_memory_fraction = 0.8,
    reserve_memory_gb = 0,
    verbose = FALSE
  )$allowed)
  expect_error(
    .resource_preflight(
      estimated_bytes = 90,
      stage = "test",
      available_bytes = 100,
      max_memory_fraction = 0.8,
      reserve_memory_gb = 0,
      verbose = FALSE
    ),
    "refused"
  )
})
