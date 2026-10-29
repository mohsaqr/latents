# Cached controls and mutation tests for the actual kit wrappers.
# Run from the repository root: Rscript equivalence/validation/check-package-kit-gates.R
# No model fitting, reference-package installation, or report mutation is needed.
library(testthat)
source(file.path("equivalence", "package-kit-runner.R"))

kits <- c("lca-packages", "lpa-packages", "lta-packages", "regression-packages",
           "pool-packages", "enumeration-packages", "growth-published", "stata-published",
           "mplus-growth", "latentgold-coverage")
compare_names <- c("compare_all_lca", "compare_all_lpa", "compare_all_lta_packages",
                    "compare_all_regression", "compare_all_pool", "compare_all_enumeration",
                    "compare_all_growth_published", "compare_all_stata", "compare_all_mplus_growth",
                    "compare_all")
directory_names <- c("packages_dir", "lpa_packages_dir", "lta_packages_dir",
                      "regression_packages_dir", "pool_packages_dir", "enumeration_packages_dir",
                      "growth_published_dir", "stata_dir", "mplus_growth_dir", "coverage_dir")

#' Original comparison test body with only the expensive fitter replaced
#' @param i Kit index.
#' @param mutation Function applied to a copy of its saved comparison.
#' @return The wrapper's test body and its evaluation environment.
wrapper_case <- function(i, mutation = identity) {
  directory <- file.path("equivalence", kits[[i]])
  env <- new.env(parent = globalenv())
  env$comparison <- mutation(utils::read.csv(file.path(directory, "comparison.csv")))
  env[[compare_names[[i]]]] <- function(...) env$comparison
  env[[directory_names[[i]]]] <- directory
  expression <- tail(parse(file.path("equivalence", paste0("test-", kits[[i]], ".R"))), 1L)[[1L]]
  stopifnot(identical(expression[[1L]], as.name("test_that")))
  list(body = expression[[3L]], env = env)
}

invisible(lapply(seq_along(kits), function(i) {
  test_that(paste(kits[[i]], "accepts complete retained comparisons"), {
    control <- wrapper_case(i)
    eval(control$body, control$env)
  })
  test_that(paste(kits[[i]], "rejects missing and nonfinite comparisons"), {
    missing <- wrapper_case(i, function(d) { d$agrees[] <- NA; d })
    expect_error(eval(missing$body, missing$env), "agreement flags must be logical and complete")
    invisible(lapply(c(NA_real_, Inf, -Inf), function(value) {
      broken <- wrapper_case(i, function(d) {
        at <- which(d$n_values > 0L)[[1L]]
        d$difference[[at]] <- value
        d
      })
      expect_error(eval(broken$body, broken$env), "finite nonnegative differences")
    }))
    empty <- wrapper_case(i, function(d) d[FALSE, ])
    expect_error(eval(empty$body, empty$env), "comparison must contain results")
    incomplete <- wrapper_case(i, function(d) d[-which(d$agrees)[[1L]], ])
    expect_s3_class(capture_expectation(withCallingHandlers(
      eval(incomplete$body, incomplete$env),
      expectation_success = function(cnd) invokeRestart("muffle_expectation"))),
      "expectation_failure")
  })
}))

test_that("parameter-count defects cannot be registered as local-maximum differences", {
  comparison <- utils::read.csv("equivalence/mplus-growth/comparison.csv")
  known <- utils::read.csv("equivalence/mplus-growth/known-differences.csv")
  at <- which(comparison$case == "runner_gmm_k3" & comparison$quantity == "number of parameters")
  expect_length(at, 1L)
  comparison$latents[at] <- 9999
  comparison$difference[at] <- abs(comparison$latents[at] - comparison$reference_value[at])
  comparison$agrees[at] <- FALSE
  expect_identical(apply_register(comparison, known)$comparison$status[at], "disagrees")
  # Even a future broad register entry must not override this invariant.
  known$quantity[] <- "."
  expect_identical(apply_register(comparison, known)$comparison$status[at], "disagrees")
})

test_that("zero-value reference failures remain explicit failures", {
  comparison <- data.frame(case = "no_fit", reference = "reference", quantity = "a solution",
                            difference = NA_real_, tolerance = NA_real_, n_values = 0L,
                            agrees = FALSE)
  known <- data.frame(case = character(), reference = character(), quantity = character(),
                      reason = character())
  expect_identical(apply_register(comparison, known)$comparison$status, "disagrees")
  comparison$agrees <- TRUE
  expect_error(apply_register(comparison, known), "missing reference solution cannot agree")
})

test_that("the runner rejects incomplete evidence before writing reports", {
  directory <- tempfile("latents-kit-gate-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  known <- data.frame(case = character(), reference = character(), quantity = character(),
                      reason = character())
  utils::write.csv(known, file.path(directory, "known-differences.csv"), row.names = FALSE)
  broken <- data.frame(case = "synthetic", reference = "reference", quantity = "estimate",
                        difference = NA_real_, tolerance = 1e-6, n_values = 1L, agrees = NA)
  expect_error(run_package_kit(directory, function() broken, character()),
                "agreement flags must be logical and complete")
  expect_false(file.exists(file.path(directory, "comparison.csv")))
  utils::write.csv(broken, file.path(directory, "comparison.csv"), row.names = FALSE)
  before <- tools::md5sum(file.path(directory, "comparison.csv"))
  expect_error(reapply_package_register(directory), "agreement flags must be logical and complete")
  expect_identical(tools::md5sum(file.path(directory, "comparison.csv")), before)
})
