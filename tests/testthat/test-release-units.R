# Fits that do not depend on the units of the data, in every family:
# SQUAREM in lta() and the additive family, the quasi-Newton finish of the
# general transition engine, and the replication tolerance of every family.
# A change of units shifts every start's log likelihood by one constant, so
# the gaps between starts, each start's iterations and the tolerance must
# be the same in any units.

units_gaps <- function(fit) {
  starts <- get_results(fit, "starts")
  starts$log_likelihood - fit$log_likelihood
}

# Scaling by a power of two is exact, so any difference left is the fit's own
# dependence on the units, not rounding of the rescaled data.
expect_same_path <- function(raw, scaled, gap = 1e-8) {
  expect_identical(get_results(scaled, "starts")$iterations,
                   get_results(raw, "starts")$iterations)
  expect_lt(max(abs(units_gaps(scaled) - units_gaps(raw))), gap)
  expect_equal(scaled$replication_tolerance, raw$replication_tolerance)
  expect_identical(scaled$n_best_replicated, raw$n_best_replicated)
}

with_local_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  code
}

# Two persistent, well-separated states: a likelihood with a clear maximum.
units_panel <- with_local_seed(8, {
  path <- unlist(lapply(seq_len(80L), function(i) {
    Reduce(function(previous, step) {
      if (stats::runif(1) < 0.85) previous else 3L - previous
    }, seq_len(4L), sample(1:2, 1L), accumulate = TRUE)
  }))
  data.frame(id = rep(seq_len(80L), each = 5L), time = rep(1:5, 80L),
             y1 = stats::rnorm(400L, c(0, 2.5)[path]),
             y2 = stats::rnorm(400L, c(0, 1.5)[path], 0.8))
})

# Group classes with distinct means and a group-level spread.
units_groups <- with_local_seed(9, {
  class <- rep(sample(1:2, 50L, replace = TRUE), each = 6L)
  intercept <- rep(stats::rnorm(50L, 0, 0.7), each = 6L)
  data.frame(id = rep(seq_len(50L), each = 6L), time = rep(0:5, 50L),
             y1 = c(0, 2.5)[class] + intercept + stats::rnorm(300L),
             y2 = c(0, 1.5)[class] + 0.5 * intercept + stats::rnorm(300L, 0, 0.8))
})

rescale <- function(data, vars, factor) {
  data[vars] <- data[vars] * factor
  data
}

test_that("lta() takes the same path in any units, under SQUAREM and plain EM", {
  skip_on_cran()
  vars <- c("y1", "y2")
  fit <- function(data, ...) {
    quietly(lta(data, vars, "id", n_profiles = 2, time = "time", n_starts = 3,
                seed = 2, ...))
  }
  # The homogeneous engine, then the general one (mover-stayer), whose
  # quasi-Newton finish used to stop short at small units.
  # The general engine's L-BFGS-B finish stops within about 1e-5 of the
  # maximum (in any units: starts at one maximum differ by that much at the
  # same scale), far inside the replication tolerance; EM is exact.
  engines <- list(homogeneous = list(engine = list(), gap = 1e-8),
                  general = list(engine = list(mover_stayer = TRUE), gap = 1e-5))
  invisible(lapply(engines, function(engine) {
    invisible(lapply(c("squarem", "none"), function(acceleration) {
      raw <- do.call(fit, c(list(units_panel, acceleration = acceleration), engine$engine))
      small <- do.call(fit, c(list(rescale(units_panel, vars, 1 / 1024),
                                   acceleration = acceleration), engine$engine))
      expect_same_path(raw, small, engine$gap)
      expect_identical(raw$replication_tolerance, 1e-6 * 2 * nrow(units_panel))
    }))
  }))
})

test_that("the additive family takes the same SQUAREM path in any units", {
  skip_on_cran()
  vars <- c("y1", "y2")
  fit <- function(data) {
    quietly(multilpa(data, vars, "id", n_group_classes = 2, family = "additive",
                     n_starts = 3, seed = 4))
  }
  raw <- fit(units_groups)
  large <- fit(rescale(units_groups, vars, 1024))
  expect_same_path(raw, large)
  expect_identical(raw$replication_tolerance, 1e-6 * 2 * nrow(units_groups))
})

test_that("every family counts replication within 1e-6 per observed value", {
  skip_on_cran()
  vars <- c("y1", "y2")
  full <- quietly(multilpa(units_groups, vars, "id", n_profiles = 2, n_group_classes = 2,
                           family = "full_cross_level", n_starts = 2, seed = 1))
  # The rows and the group means.
  expect_identical(full$replication_tolerance,
                   1e-6 * (2 * nrow(units_groups) + 2 * 50))
  covariate_data <- transform(units_panel, z = rep(seq(-1, 1, length.out = 80L), each = 5L))
  covariates <- quietly(multilpa(covariate_data, vars, "id", 2, 1,
                                 profile_covariates = "z", n_starts = 2, seed = 1))
  expect_identical(covariates$replication_tolerance, 1e-6 * 2 * nrow(units_panel))
  regression_data <- transform(units_panel, x = sin(seq_len(nrow(units_panel))))
  regression <- quietly(mixture_regression(y1 ~ x, regression_data, n_classes = 2,
                                           n_starts = 2, seed = 1, vcov_type = "none"))
  expect_identical(regression$replication_tolerance, 1e-6 * nrow(units_panel))
  scaled <- quietly(mixture_regression(y1 ~ x, rescale(regression_data, "y1", 1e3),
                                       n_classes = 2, n_starts = 2, seed = 1,
                                       vcov_type = "none"))
  expect_identical(scaled$n_best_replicated, regression$n_best_replicated)
  growth <- quietly(mixture_regression(y1 ~ time, units_groups, n_classes = 2, id = "id",
                                       class_level = "group", random = ~ 1,
                                       n_starts = 2, seed = 1, vcov_type = "none"))
  expect_identical(growth$replication_tolerance, 1e-6 * nrow(units_groups))
})

test_that("sensitivity() reads the fit's own replication tolerance", {
  skip_on_cran()
  growth <- quietly(mixture_regression(y1 ~ time, units_groups, n_classes = 2, id = "id",
                                       class_level = "group", random = ~ 1,
                                       n_starts = 2, seed = 1, vcov_type = "none"))
  table <- quietly(sensitivity(growth, seeds = 1:2))
  expect_identical(attr(table, "tolerance"), growth$replication_tolerance)
})
