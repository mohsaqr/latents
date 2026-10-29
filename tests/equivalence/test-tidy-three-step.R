# Comparisons with independent references moved from tests/testthat/test-tidy-three-step.R,
# so that they run in CI's equivalence job rather than under R CMD check.

.tidy_step_data <- function(seed = 11L, n_groups = 50L, per = 10L,
                            separation = 1.1, difference = 10) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per)
  n <- length(group)
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  data.frame(g = group,
             a = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             b = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             x = stats::rnorm(n),
             y = stats::rnorm(n, ifelse(truth == 2L, difference, 0)))
}

# The fixture separates its classes clearly, so one start is enough, and most
# tests fit the same data: each fit is computed once and reused.
.tidy_step_fit <- local({
  cache <- list()
  function(data, n_profiles = 2L) {
    hit <- Filter(function(entry) entry$n_profiles == n_profiles &&
                    identical(entry$data, data), cache)
    if (length(hit)) return(hit[[1L]]$fit)
    fit <- multilpa(data, c("a", "b"), "g", n_profiles = n_profiles,
                    n_group_classes = 1, n_starts = 1, seed = 1)
    cache[[length(cache) + 1L]] <<- list(data = data, n_profiles = n_profiles,
                                         fit = fit)
    fit
  }
})

test_that("the pairwise standard error agrees with a cluster bootstrap", {
  skip_on_cran()
  # Weak separation on purpose: the BCH weights for the two classes are then
  # strongly anticorrelated, so a standard error that ignored the covariance
  # between the two class means would miss by more than the resampling noise.
  data <- .tidy_step_data(separation = 0.8)
  fit <- .tidy_step_fit(data)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  weights <- .multilpa_step_weights(pieces, "bch")
  analytic <- three_step(fit, data, "y", contrast = "pairs")

  # Resample whole groups, holding the step-one weights fixed, which is exactly
  # the assumption the cluster-robust variance makes. The influence-function
  # standard error should reproduce the spread this generates.
  set.seed(99)
  rows_by_group <- split(seq_len(nrow(data)), pieces$group_index)
  replicates <- replicate(600L, {
    picked <- sample(seq_along(rows_by_group), replace = TRUE)
    rows <- unlist(rows_by_group[picked], use.names = FALSE)
    taken <- weights[rows, , drop = FALSE]
    means <- colSums(taken * data$y[rows]) / colSums(taken)
    means[[2L]] - means[[1L]]
  })

  ratio <- stats::sd(replicates) / analytic$standard_error
  expect_gt(ratio, 0.90)
  expect_lt(ratio, 1.16)
  # Treating the two class means as independent would inflate the ratio past
  # that bound, which is what makes this a calibration and not a tautology.
  means <- three_step(fit, data, "y")
  independent <- sqrt(sum(means$standard_error^2))
  expect_gt(stats::sd(replicates) / independent, 1.16)
})
