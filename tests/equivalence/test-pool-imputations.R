# Comparisons with independent references moved from tests/testthat/test-pool-imputations.R,
# so that they run in CI's equivalence job rather than under R CMD check.

.pool_fixture <- function(n_groups = 40L, size = 8L, seed = 5L) {
  set.seed(seed)
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size), z = stats::rnorm(n))
  frame$w <- stats::ave(stats::rnorm(n), frame$g)
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z))
  frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  frame$y2 <- stats::rnorm(n, c(0, 2)[state])
  frame
}

# Stand-in imputations: the same data with a few covariate values redrawn,
# which is all the pooling step can see of an imputation model.
.pool_imputed <- function(frame, m = 3L, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(m), function(index) {
    copy <- frame
    copy$z[seq_len(20L)] <- stats::rnorm(20L)
    copy
  })
}

.pool_quietly <- function(expr) {
  quietly(expr, classes = c(.multilpa_expected_warnings, "latents_extreme_coefficients"))
}

test_that("Rubin's rules agree with mice::pool.scalar for every parameter", {
  skip_if_not_installed("mice")
  pooled <- .pool_quietly(pool_imputations(
    .pool_imputed(.pool_fixture()), c("y1", "y2"), "g", 2L,
    n_group_classes = 1L, profile_covariates = "z", n_starts = 1, seed = 1))
  estimates <- as.data.frame(pooled)
  per_imputation <- get_results(pooled, "imputations")
  key <- function(table) paste(table$level, table$outcome, table$term, table$parameter)
  invisible(lapply(seq_len(nrow(estimates)), function(row) {
    rows <- per_imputation[key(per_imputation) == key(estimates)[row], ]
    reference <- mice::pool.scalar(rows$estimate, rows$standard_error^2, n = Inf)
    expect_equal(estimates$estimate[row], reference$qbar, tolerance = 1e-12)
    expect_equal(estimates$standard_error[row]^2, reference$t, tolerance = 1e-12)
    expect_equal(estimates$df[row], reference$df, tolerance = 1e-10)
    expect_equal(estimates$fmi[row], reference$fmi, tolerance = 1e-10)
  }))
})
