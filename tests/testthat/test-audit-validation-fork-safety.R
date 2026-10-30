skip_on_cran()

validation_fork_guard <- function() {
  path <- testthat::test_path("..", "..", "validation", "fork-runtime.R")
  skip_if_not(file.exists(path), "validation sources are outside the installed package")
  definition <- Filter(function(x) {
    is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(x[[2L]], as.name(".validation_fork_cores"))
  }, as.list(parse(path)))[[1L]]
  eval(definition[[3L]], baseenv())
}

test_that("validation preflight refuses unsafe forks without changing core counts", {
  guard <- validation_fork_guard()
  expect_identical(guard(4L, "unix", "Linux", "reference BLAS"), 4L)
  expect_identical(guard(4, "unix", "Darwin", "libRblas.dylib"), 4)
  expect_identical(guard(1L, "windows", "Windows", "BLAS"), 1L)
  expect_identical(guard(1, "unix", "Darwin", "Accelerate/vecLib"), 1)
  expect_error(guard(2L, "windows", "Windows", "BLAS"),
                 class = "latents_validation_unsafe_fork")
  expect_error(guard(2L, "unix", "Darwin", "Accelerate/vecLib"),
                 class = "latents_validation_unsafe_fork")
  expect_error(guard(2L, "unix", "Darwin", "libRblas.vecLib.dylib"),
                 "reference BLAS")
  expect_error(guard(2L, "unix", "Darwin", "Accelerate"), "one-core")
  invisible(lapply(list(0, -1, 1.5, NA_real_, Inf, c(1, 2)), function(cores) {
    expect_error(guard(cores, "unix", "Linux", "BLAS"), "positive whole")
  }))
})

test_that("validation preflight preserves explicit arguments and the RNG state", {
  guard <- validation_fork_guard()
  withr::local_seed(708)
  before <- .Random.seed
  supplied <- list(mc.cores = 3L, mc.set.seed = TRUE, mc.preschedule = FALSE)
  recorded <- (function(mc.cores, mc.set.seed, mc.preschedule) {
    list(mc.cores = mc.cores, mc.set.seed = mc.set.seed, mc.preschedule = mc.preschedule)
  })(guard(supplied$mc.cores, "unix", "Linux", "reference BLAS"),
     supplied$mc.set.seed, supplied$mc.preschedule)
  expect_identical(recorded, supplied)
  expect_error(guard(3L, "unix", "Darwin", "vecLib"),
                 class = "latents_validation_unsafe_fork")
  expect_identical(.Random.seed, before)
})

test_that("every retained study guards each direct mclapply core request", {
  files <- c("additive-diagnostic.R", "additive-simulation.R", "bootstrap-lrt-size.R",
    "covariate-missing-coverage.R", "covariate-slopes-coverage.R",
    "families-simulation.R", "growth-simulation.R", "lta-simulation.R",
    "mixture-regression-recovery.R", "montecarlo-recovery.R",
    "pool-imputations-coverage.R", "staged-bootstrap-coverage.R",
    "transition-inference-coverage.R", paste0("review-2026-10-02-followup/",
      c("simulation-repro-additive.R", "simulation-repro-families.R",
        "simulation-repro-lta.R", "transition-bootstrap-validation.R")))
  paths <- file.path(testthat::test_path("..", "..", "validation"), files)
  skip_if_not(all(file.exists(paths)), "validation sources are outside the installed package")
  calls <- list()
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    head <- x[[1L]]
    if (is.call(head) && identical(head[[1L]], as.name("::")) &&
        identical(head[[2L]], as.name("parallel")) &&
        identical(head[[3L]], as.name("mclapply"))) {
      calls[[length(calls) + 1L]] <<- x
    }
    invisible(lapply(as.list(x)[-1L], walk))
  }
  invisible(lapply(paths, function(path) {
    expressions <- parse(path)
    expect_true(any(vapply(expressions, function(x) {
      is.call(x) && identical(x[[1L]], as.name("source")) &&
        grepl("fork-runtime.R", paste(deparse(x), collapse = ""), fixed = TRUE)
    }, logical(1))))
    invisible(lapply(expressions, walk))
  }))
  expect_length(calls, 18L)
  invisible(lapply(calls, function(call) {
    requested <- call$mc.cores
    expect_true(is.call(requested))
    expect_identical(requested[[1L]], as.name(".validation_fork_cores"))
    expect_length(requested, 2L)
  }))
})
