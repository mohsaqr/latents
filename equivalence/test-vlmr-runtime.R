# Runtime dispatch only: mocked fork execution never launches unsafe workers.

vlmr_runtime_environment <- function() {
  environment <- new.env(parent = asNamespace("latents"))
  sys.source(file.path("features", "vlmr", "helpers.R"), envir = environment)
  sys.source(file.path("features", "vlmr", "reference.R"), envir = environment)
  environment
}

test_that("VLMR uses serial execution when requested or on unsafe platforms", {
  environment <- vlmr_runtime_environment()
  withr::local_options(latents.vlmr.serial = TRUE)
  testthat::local_mocked_bindings(.latents_fork_safe = function(...) stop("not needed"),
                                  .package = "latents")
  testthat::local_mocked_bindings(mclapply = function(...) stop("unsafe fork"),
                                  .package = "parallel")
  expect_identical(environment$vlmr_map_cases(list(a = 1L, b = 2L), identity),
                   list(a = 1L, b = 2L))
  options(latents.vlmr.serial = FALSE)
  environment$.Platform <- list(OS.type = "windows")
  expect_identical(environment$vlmr_map_cases(list(1L), identity), list(1L))
})

test_that("VLMR automatically avoids Accelerate forks and preserves safe worker controls", {
  environment <- vlmr_runtime_environment()
  environment$.Platform <- list(OS.type = "unix")
  withr::local_options(latents.vlmr.serial = FALSE)
  safe <- FALSE
  seen <- NULL
  testthat::local_mocked_bindings(.latents_fork_safe = function(...) safe,
                                  .package = "latents")
  testthat::local_mocked_bindings(mclapply = function(x, fun, mc.cores, mc.preschedule) {
    seen <<- list(cores = mc.cores, preschedule = mc.preschedule)
    lapply(x, fun)
  }, .package = "parallel")
  expect_identical(environment$vlmr_map_cases(list(1L, 2L), identity), list(1L, 2L))
  expect_null(seen)
  safe <- TRUE
  expect_identical(environment$vlmr_map_cases(list(a = 1L, b = 2L), identity, cores = 3L),
                   list(a = 1L, b = 2L))
  expect_identical(seen, list(cores = 3L, preschedule = FALSE))
})

test_that("a real seeded VLMR fit and reference preserve results across runtime routes", {
  environment <- vlmr_runtime_environment()
  environment$.Platform <- list(OS.type = "unix")
  withr::local_seed(326)
  data <- data.frame(y1 = rep(c(-2, 2), each = 60L) + rnorm(120L),
                     y2 = rep(c(-2, 2), each = 60L) + rnorm(120L))
  case <- list(data = data, external = FALSE, config = "runtime",
               spec = list(model = "EEI", range = 1:2, true_k = 2L, vars = names(data),
                           weights = NULL, missing = FALSE, n = nrow(data)))
  environment$vlmr_cases <- function() list(runtime = case)
  withr::local_envvar(VLMR_FIT_CACHE = NA_character_)
  testthat::local_mocked_bindings(.latents_fork_safe = function(...) TRUE,
                                  .package = "latents")
  testthat::local_mocked_bindings(mclapply = function(x, fun, mc.cores, mc.preschedule) {
    lapply(x, fun)
  }, .package = "parallel")
  before <- .Random.seed
  withr::local_options(latents.vlmr.serial = TRUE)
  serial <- environment$vlmr_compute_all("runtime", draws = 200L)
  expect_identical(.Random.seed, before)
  options(latents.vlmr.serial = FALSE)
  mapped <- environment$vlmr_compute_all("runtime", draws = 200L)
  expect_identical(.Random.seed, before)
  expect_identical(mapped, serial)
  expect_length(serial$errors, 0L)
  expect_equal(nrow(serial$table), 1L)
  expect_true(is.finite(serial$table$lr))
  expect_identical(serial$table$ref_draws, 200L)
  expect_lt(abs(serial$table$lr - serial$table$ref_lr), 1e-9)
})
