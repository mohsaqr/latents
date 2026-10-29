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
  pieces <- .multilpa_level_assignments(fit, "individuals")
  by_assigned <- t(.multilpa_error_matrix(pieces))[pieces$modal, , drop = FALSE]
  design <- cbind(1, data$x)
  objective <- function(theta) {
    prior <- .multilpa_logit_prior(matrix(theta, 2L, 1L), design)
    -sum(log(rowSums(prior * by_assigned)))
  }
  tight <- stats::optim(c(0, 0), objective, method = "BFGS",
                        control = list(maxit = 5000L, reltol = 1e-15))
  tight <- stats::optim(tight$par, objective, method = "Nelder-Mead",
                        control = list(maxit = 20000L, reltol = 1e-15))
  # Before the repair the default BFGS tolerance left it 1e-4 away.
  expect_equal(result$estimate, tight$par, tolerance = 1e-6)
})
