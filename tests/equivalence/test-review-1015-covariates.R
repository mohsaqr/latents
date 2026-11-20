# Comparisons with independent references moved from tests/testthat/test-review-1015-covariates.R,
# so that they run in CI's equivalence job rather than under R CMD check.

review_r3step_data <- function() {
  set.seed(21)
  g <- rep(seq_len(50), each = 12)
  x <- stats::rnorm(600)
  truth <- 1L + as.integer(stats::runif(600) < stats::plogis(-0.3 + 1.2 * x))
  data.frame(
    g = g, x = x,
    a = stats::rnorm(600, ifelse(truth == 2L, 1.2, -1.2)),
    b = stats::rnorm(600, ifelse(truth == 2L, 1.2, -1.2))
  )
}

test_that("r3step() reaches the maximum of its likelihood", {
  data <- review_r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)
  result <- r3step(fit, data, "x")
  # The step-three likelihood written from its definition (Vermunt 2010):
  # modal assignments, D[t, s] the posterior mass of true class t among units
  # assigned to s over the mass of t, and a two-class logit prior on x.
  posteriors <- fit$subject_posteriors
  modal <- max.col(posteriors, ties.method = "first")
  expect_identical(as.integer(fit$subject_profiles), modal)
  errors <- vapply(1:2, \(s) colSums(posteriors * (modal == s)), numeric(2)) /
    colSums(posteriors)
  by_assigned <- t(errors)[modal, , drop = FALSE]
  design <- cbind(1, data$x)
  objective <- function(theta) {
    eta <- cbind(design %*% theta, 0)
    prior <- exp(eta - apply(eta, 1L, max))
    prior <- prior / rowSums(prior)
    -sum(log(rowSums(prior * by_assigned)))
  }
  tight <- stats::optim(c(0, 0), objective, method = "BFGS",
                        control = list(maxit = 5000L, reltol = 1e-15))
  tight <- stats::optim(tight$par, objective, method = "Nelder-Mead",
                        control = list(maxit = 20000L, reltol = 1e-15))
  # Before the repair the default BFGS tolerance left it 1e-4 away. 1e-6
  # relative, the original tolerance: two optimisers, each to its own criterion.
  expect_reference(result$estimate, tight$par,
                   reference = "independent: stats::optim() maximisation of the step-three likelihood",
                   covers = "r3step()", quantity = "R3STEP logit coefficients",
                   tolerance = 1e-6, measure = "relative", key = c("(Intercept)", "x"))
})
