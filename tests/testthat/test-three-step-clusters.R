## How many independent groups a cluster-robust three-step standard error needs.
##
## The meat of a cluster-robust sandwich is the sum of one outer product per
## group, formed from that group's summed influence contributions. Those
## contributions sum to zero at the estimate -- that sum *is* the stationarity
## condition the estimate solves -- so the meat has rank at most `G - 1`. With
## one group it is exactly zero, and version 0.10.0 read a standard error of
## 5.9e-17 off it and reported identical confidence endpoints. With as many
## groups as reported quantities it is singular, and some contrast among them is
## handed a variance of exactly zero. Both are refused here.

.clusters_one_group <- function() {
  set.seed(28)
  data <- data.frame(g = 1, a = stats::rnorm(200), y = stats::rnorm(200))
  list(data = data,
       fit = multilpa(data, "a", "g", n_profiles = 1, n_group_classes = 1,
                      n_starts = 1, seed = 1))
}

.clusters_many_groups <- function(n_groups = 40L, per = 8L, seed = 404L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per)
  n <- length(group)
  data <- data.frame(g = group, a = stats::rnorm(n),
                     y = stats::rnorm(n, 3, 2))
  list(data = data,
       fit = multilpa(data, "a", "g", n_profiles = 1, n_group_classes = 1,
                      n_starts = 2, seed = 1))
}

## The influence contributions `three_step()` clusters, recomputed here so the
## rank claim the guard rests on is checked against the real quantities rather
## than against a hand-made matrix.
.clusters_influence <- function(fit, data, outcome, method = "bch") {
  pieces <- .multilpa_level_assignments(fit, "individuals")
  weights <- .multilpa_step_weights(pieces, method)
  totals <- colSums(weights)
  values <- data[[outcome]]
  estimates <- as.vector(crossprod(weights, values)) / totals
  sweep(weights * outer(values, estimates, "-"), 2L, totals, "/")
}

test_that("a single independent group is refused, not given a zero standard error", {
  parts <- .clusters_one_group()
  expect_equal(parts$fit$n_groups, 1)

  expect_error(three_step(parts$fit, parts$data, "y"),
               class = "multilpa_too_few_groups")
  expect_error(three_step(parts$fit, parts$data, "y", method = "modal"),
               class = "multilpa_too_few_groups")
  expect_error(three_step(parts$fit, parts$data, "y", method = "proportional"),
               class = "multilpa_too_few_groups")
  expect_error(three_step(parts$fit, parts$data, "y", contrast = "pairs"),
               class = "multilpa_too_few_groups")
})

test_that("the sole group's influence contributions really do cancel to nothing", {
  parts <- .clusters_one_group()
  influence <- .clusters_influence(parts$fit, parts$data, "y")
  clustered <- rowsum(influence, parts$fit$group_index, reorder = FALSE)

  expect_equal(nrow(clustered), 1L)
  # Not "small": zero to within the rounding of summing 200 numbers, and many
  # orders of magnitude below the contributions it was built from.
  expect_lt(abs(drop(clustered)), 1e-14)
  expect_lt(abs(drop(clustered)) / max(abs(influence)), 1e-10)
})

test_that("the unclustered variance is offered explicitly and is the right size", {
  parts <- .clusters_one_group()
  result <- three_step(parts$fit, parts$data, "y", vcov_type = "independent")

  expect_identical(attr(result, "vcov_type"), "independent")
  expect_named(result, c("level", "method", "class", "estimate",
                         "standard_error", "conf_low", "conf_high",
                         "effective_n"))
  # A one-profile fit gives every unit the weight 1, so the influence
  # contribution of observation i is (y_i - ybar) / n and the unclustered
  # standard error is exactly sd(y) * sqrt((n - 1) / n) / sqrt(n). That is an
  # algebraic identity, not an approximation, so it is asserted tightly.
  n <- nrow(parts$data)
  identity <- stats::sd(parts$data$y) * sqrt((n - 1) / n) / sqrt(n)
  expect_equal(result$standard_error, identity, tolerance = 1e-12)
  # And it sits within a percent of the order of magnitude the review names as
  # the sane one, sd(y) / sqrt(n), instead of 5.9e-17.
  reference <- stats::sd(parts$data$y) / sqrt(n)
  expect_lt(abs(result$standard_error / reference - 1), 0.01)
  expect_gt(result$conf_high - result$conf_low, 0.1)
})

test_that("enough groups give a standard error of a plausible size", {
  skip_on_cran()
  parts <- .clusters_many_groups()
  result <- three_step(parts$fit, parts$data, "y")

  expect_identical(attr(result, "vcov_type"), "cluster")
  expect_equal(result, three_step(parts$fit, parts$data, "y",
                                  vcov_type = "cluster"))
  # The outcome is drawn independently of the grouping, so there is no true
  # intraclass correlation and the cluster-robust standard error should track
  # the independent reference sd(y) / sqrt(n). The CR0 variance estimate from
  # G groups has a relative standard deviation of about sqrt(2 / (G - 1)),
  # here 0.23, so about 11% on the standard error; the band below is roughly
  # three of those either way and is a magnitude check, not a calibration.
  reference <- stats::sd(parts$data$y) / sqrt(nrow(parts$data))
  ratio <- result$standard_error / reference
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
  expect_equal(result$conf_high - result$conf_low,
               2 * stats::qnorm(0.975) * result$standard_error)
})

test_that("the refusal threshold is the rank the cluster sums can reach", {
  skip_on_cran()
  parts <- .clusters_many_groups()
  influence <- .clusters_influence(parts$fit, parts$data, "y")
  # Every column sums to zero over the units, which is what caps the rank of
  # any clustering of them at one less than the number of clusters.
  expect_lt(max(abs(colSums(influence))) / max(abs(influence)), 1e-10)
  ranks <- vapply(2:5, function(n_clusters) {
    folded <- rowsum(influence, rep_len(seq_len(n_clusters), nrow(influence)),
                     reorder = FALSE)
    qr(folded)$rank
  }, integer(1))
  expect_true(all(ranks <= (2:5) - 1L))

  # The guard passes exactly when the units outnumber the quantities, which is
  # the same as asking the meat to have full rank for them.
  expect_true(.multilpa_require_clusters(5L, 4L, "A variance"))
  expect_error(.multilpa_require_clusters(4L, 4L, "A variance"),
               class = "multilpa_too_few_groups")
  expect_error(.multilpa_require_clusters(1L, 1L, "A variance"),
               class = "multilpa_too_few_groups")
  expect_error(.multilpa_require_clusters(0L, 1L, "A variance"),
               class = "multilpa_too_few_groups")
  # A one-row score matrix cannot make a sandwich, whoever asks for it.
  expect_error(.multilpa_cross_product(matrix(c(1, 2, 3), 1L, 3L)),
               class = "multilpa_too_few_groups")
  expect_equal(.multilpa_cross_product(matrix(c(1, 2), 2L, 1L)),
               matrix(5, 1L, 1L))
})

test_that("two groups cannot carry a variance for two classes", {
  skip_on_cran()
  set.seed(77)
  group <- rep(seq_len(2), each = 60)
  n <- length(group)
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  data <- data.frame(
    g = group,
    a = stats::rnorm(n, ifelse(truth == 2L, 1.5, -1.5)),
    b = stats::rnorm(n, ifelse(truth == 2L, 1.5, -1.5)),
    y = stats::rnorm(n, ifelse(truth == 2L, 5, 0)))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)
  expect_equal(fit$n_groups, 2)

  # Two groups leave one independent contribution for a two-by-two covariance,
  # so one contrast between the class means would have variance exactly zero.
  expect_error(three_step(fit, data, "y"), class = "multilpa_too_few_groups")
  expect_error(three_step(fit, data, "y", contrast = "pairs"),
               class = "multilpa_too_few_groups")
  # Treating the 120 observations as independent is a different assumption and
  # is available under its own name.
  unclustered <- three_step(fit, data, "y", vcov_type = "independent")
  expect_identical(attr(unclustered, "vcov_type"), "independent")
  expect_true(all(unclustered$standard_error > 1e-3))
})

test_that("r3step refuses a robust covariance its groups cannot support", {
  skip_on_cran()
  set.seed(29)
  one_group <- data.frame(
    g = 1, x = stats::rnorm(200),
    a = c(stats::rnorm(100, -2), stats::rnorm(100, 2)))
  fit <- multilpa(one_group, "a", "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)
  expect_error(r3step(fit, one_group, "x", vcov_type = "robust"),
               class = "multilpa_too_few_groups")
  # The model-based alternative is still available and is labelled as such.
  observed <- r3step(fit, one_group, "x", vcov_type = "observed")
  expect_identical(attr(observed, "vcov_type"), "observed")
  expect_true(all(observed$standard_error > 1e-3))

  # Three groups against three coefficients: the sandwich is singular, not
  # merely imprecise, so it is refused for the same reason.
  set.seed(30)
  group <- rep(seq_len(3), each = 60)
  n <- length(group)
  truth <- 1L + as.integer(stats::runif(n) < stats::plogis(-0.3 + 1.2 * stats::rnorm(n)))
  three <- data.frame(
    g = group, x = stats::rnorm(n), w = stats::rnorm(n),
    a = stats::rnorm(n, ifelse(truth == 2L, 1.5, -1.5)),
    b = stats::rnorm(n, ifelse(truth == 2L, 1.5, -1.5)))
  small <- multilpa(three, c("a", "b"), "g", n_profiles = 2,
                    n_group_classes = 1, n_starts = 4, seed = 1)
  expect_equal(small$n_groups, 3)
  expect_error(r3step(small, three, c("x", "w"), vcov_type = "robust"),
               class = "multilpa_too_few_groups")
  # Two coefficients against three groups is within the rule and goes through.
  expect_true(all(r3step(small, three, "x", vcov_type = "robust")$standard_error > 0))
})

test_that("a robust r3step with many groups is unchanged and comparable", {
  skip_on_cran()
  set.seed(21)
  group <- rep(seq_len(50), each = 12L)
  n <- length(group)
  covariate <- stats::rnorm(n)
  truth <- 1L + as.integer(stats::runif(n) < stats::plogis(-0.3 + 1.2 * covariate))
  data <- data.frame(
    g = group, x = covariate,
    a = stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)),
    b = stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)

  robust <- r3step(fit, data, "x", vcov_type = "robust")
  observed <- r3step(fit, data, "x", vcov_type = "observed")
  expect_identical(attr(robust, "vcov_type"), "robust")
  expect_equal(robust$estimate, observed$estimate)
  expect_true(all(robust$standard_error > 0))
  # With 50 groups and no true clustering the two variances should agree to
  # within a factor of two; the point is that the robust one is an estimate,
  # not that it equals the model-based one.
  ratio <- robust$standard_error / observed$standard_error
  expect_true(all(ratio > 0.5 & ratio < 2))
})
