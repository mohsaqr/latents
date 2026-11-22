# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

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
  skip_without_reference("mice")
  pooled <- .pool_quietly(pool_imputations(
    .pool_imputed(.pool_fixture()), c("y1", "y2"), "g", 2L,
    n_group_classes = 1L, profile_covariates = "z", n_starts = 1, seed = 1))
  keys <- c("level", "outcome", "term", "parameter")
  # The package's pooled table, with the total variance beside the standard error.
  estimates <- transform(as.data.frame(pooled), variance = standard_error^2)
  # mice::pool.scalar() on each parameter's per-imputation estimates and variances.
  per_imputation <- get_results(pooled, "imputations")
  reference <- do.call(rbind, lapply(split(per_imputation, per_imputation[keys], drop = TRUE), \(rows) {
    scalar <- mice::pool.scalar(rows$estimate, rows$standard_error^2, n = Inf)
    data.frame(rows[1L, keys], estimate = scalar$qbar, variance = scalar$t,
               df = scalar$df, fmi = scalar$fmi, row.names = NULL)
  }))
  covers <- c("pool_imputations()", "as.data.frame.latents_pooled()")
  expect_reference(estimates, reference[c(keys, "estimate", "variance")], reference = "mice",
                   covers = covers, quantity = "pooled estimate and total variance",
                   tolerance = 1e-12, measure = "relative", key = keys)
  expect_reference(estimates, reference[c(keys, "df", "fmi")], reference = "mice",
                   covers = covers, quantity = "degrees of freedom and fraction of missing information",
                   tolerance = 1e-10, measure = "relative", key = keys)
})
