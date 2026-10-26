# Runtime dispatch only: mocked fork execution never launches unsafe workers.

ml3_runtime_environment <- function() {
  environment <- new.env(parent = asNamespace("latents"))
  sys.source(file.path("features", "ml-three-step", "helpers.R"), envir = environment)
  environment
}

test_that("ML3 uses serial execution when requested or on unsafe platforms", {
  environment <- ml3_runtime_environment()
  withr::local_options(latents.ml_three_step.serial = TRUE)
  testthat::local_mocked_bindings(.latents_fork_safe = function(...) stop("not needed"),
                                  .package = "latents")
  testthat::local_mocked_bindings(mclapply = function(...) stop("unsafe fork"),
                                  .package = "parallel")
  expect_identical(environment$ml3_map_cases(list(a = 1L, b = 2L), identity),
                   list(a = 1L, b = 2L))
  options(latents.ml_three_step.serial = FALSE)
  environment$.Platform <- list(OS.type = "windows")
  expect_identical(environment$ml3_map_cases(list(1L), identity), list(1L))
})

test_that("ML3 automatically avoids Accelerate forks and retains safe-fork controls", {
  environment <- ml3_runtime_environment()
  environment$.Platform <- list(OS.type = "unix")
  withr::local_options(latents.ml_three_step.serial = FALSE)
  safe <- FALSE
  seen <- NULL
  testthat::local_mocked_bindings(.latents_fork_safe = function(...) safe,
                                  .package = "latents")
  testthat::local_mocked_bindings(mclapply = function(x, fun, mc.cores, mc.preschedule) {
    seen <<- list(cores = mc.cores, preschedule = mc.preschedule)
    lapply(x, fun)
  }, .package = "parallel")
  expect_identical(environment$ml3_map_cases(list(1L, 2L), identity), list(1L, 2L))
  expect_null(seen)
  safe <- TRUE
  expect_identical(environment$ml3_map_cases(list(a = 1L, b = 2L), identity),
                   list(a = 1L, b = 2L))
  expect_identical(seen, list(cores = 2L, preschedule = FALSE))
})

test_that("actual fixed-seed ML3 case fits are identical across runtime routes", {
  environment <- ml3_runtime_environment()
  environment$.Platform <- list(OS.type = "unix")
  cases <- readRDS(file.path("features", "ml-three-step", "kit", "cases.rds"))[1:2]
  testthat::local_mocked_bindings(.latents_fork_safe = function(...) TRUE,
                                  .package = "latents")
  testthat::local_mocked_bindings(mclapply = function(x, fun, mc.cores, mc.preschedule) {
    lapply(x, fun)
  }, .package = "parallel")
  withr::local_seed(526)
  before <- .Random.seed
  withr::local_options(latents.ml_three_step.serial = TRUE)
  serial <- environment$ml3_map_cases(cases, environment$ml3_fit)
  expect_identical(.Random.seed, before)
  options(latents.ml_three_step.serial = FALSE)
  mapped <- environment$ml3_map_cases(cases, environment$ml3_fit)
  expect_identical(.Random.seed, before)
  expect_identical(mapped, serial)
  expect_true(all(vapply(serial, function(fit) isTRUE(fit$converged), logical(1))))
})
