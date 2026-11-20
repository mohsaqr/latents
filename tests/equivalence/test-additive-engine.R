# Additive model engine: sufficient statistics, exact density, group
# posteriors. References come from helper-additive-dense.R and never use
# the sufficient-statistic formulas under test. The package side is the
# internal engine behind multilpa(family = "additive"); it has no exported
# accessor at a parameter point of the test's choosing.

additive_parameters <- list(
  means = rbind(c(-1, 0.5), c(1.4, -0.8)),
  between = rbind(c(0.3, 0.7), c(0.8, 0.2)),
  within = c(0.6, 1.1), weights = c(0.4, 0.6))
additive_sizes <- c(1L, 2L, 3L, 5L, 7L, 4L, 6L, 2L, 8L, 3L, 5L, 4L)
additive_data <- additive_draw(additive_parameters, additive_sizes, seed = 20260930)
additive_vars <- c("y1", "y2")

# Dense log densities, one row per group and one column per class.
additive_engine_dense <- function(blocks) {
  t(vapply(blocks, \(x) vapply(1:2, \(h) additive_dense_density(
    x, additive_parameters$means[h, ], additive_parameters$between[h, ],
    additive_parameters$within), numeric(1)), numeric(2)))
}

test_that("compact density equals the full-covariance Cholesky density", {
  stats <- latents:::.additive_prepare(additive_data, additive_vars, "group")
  compact <- latents:::.additive_log_density(stats, additive_parameters)
  dense <- additive_engine_dense(additive_blocks(additive_data, additive_vars, "group"))
  # Two factorisations of one Gaussian density; rounding only.
  expect_reference(unname(compact), dense,
                   reference = "independent: dense Cholesky density (helper-additive-dense.R)",
                   covers = "multilpa()", quantity = "group log density, equal within variances",
                   tolerance = 1e-10, measure = "relative",
                   note = "Package side: the internal sufficient-statistic density (.additive_log_density) of multilpa(family = 'additive'); rows are groups, columns classes.")
  expect_identical(rownames(compact), stats$ids)
})

test_that("group posteriors are normalized and match Bayes' rule on dense densities", {
  stats <- latents:::.additive_prepare(additive_data, additive_vars, "group")
  expectation <- latents:::.additive_expectation(stats, additive_parameters)
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  joint <- additive_engine_dense(blocks) +
    rep(log(additive_parameters$weights), each = length(blocks))
  reference <- exp(joint) / rowSums(exp(joint))
  expect_reference(unname(expectation$posterior), reference,
                   reference = "independent: Bayes' rule on dense Cholesky densities (helper-additive-dense.R)",
                   covers = "multilpa()", quantity = "group class posterior",
                   tolerance = 1e-10, measure = "relative",
                   note = "Package side: the internal E-step (.additive_expectation) of multilpa(family = 'additive').")
  expect_equal(rowSums(expectation$posterior), rep(1, length(blocks)),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_reference(expectation$log_likelihood, sum(log(rowSums(exp(joint)))),
                   reference = "independent: Bayes' rule on dense Cholesky densities (helper-additive-dense.R)",
                   covers = "multilpa()", quantity = "log likelihood at the test's parameters",
                   tolerance = 1e-10, measure = "relative",
                   note = "Package side: the internal E-step (.additive_expectation) of multilpa(family = 'additive').")
  expect_equal(sum(expectation$group_log_likelihood), expectation$log_likelihood)
})

test_that("intercept moments match joint-Gaussian conditioning", {
  stats <- latents:::.additive_prepare(additive_data, additive_vars, "group")
  expectation <- latents:::.additive_expectation(stats, additive_parameters)
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  # One row per group, class and indicator, on both sides.
  cells <- expand.grid(indicator = seq_along(additive_vars), class = 1:2,
                       group = seq_along(blocks))
  reference <- lapply(seq_along(blocks), \(j) lapply(1:2, \(h) additive_dense_moments(
    blocks[[j]], additive_parameters$means[h, ],
    additive_parameters$between[h, ], additive_parameters$within)))
  expected <- data.frame(cells,
    mean = vapply(seq_len(nrow(cells)), \(i) {
      reference[[cells$group[i]]][[cells$class[i]]]["mean", cells$indicator[i]]
    }, numeric(1)),
    variance = vapply(seq_len(nrow(cells)), \(i) {
      reference[[cells$group[i]]][[cells$class[i]]]["variance", cells$indicator[i]]
    }, numeric(1)))
  actual <- data.frame(cells,
    mean = expectation$conditional_mean[cbind(cells$group, cells$indicator, cells$class)],
    variance = expectation$conditional_variance[cbind(cells$group, cells$indicator, cells$class)])
  # The original bound was absolute: max |difference| < 1e-10.
  expect_reference(actual, expected,
                   reference = "independent: joint-Gaussian conditioning (helper-additive-dense.R)",
                   covers = "multilpa()", quantity = "posterior moments of the group intercept",
                   tolerance = 1e-10, measure = "absolute", key = c("group", "class", "indicator"),
                   note = "Package side: the internal E-step (.additive_expectation) of multilpa(family = 'additive').")
})

# Moved from tests/testthat/test-additive-engine.R: comparisons with independent references.

test_that("one-indicator density equals numerical integration over the intercept", {
  x <- c(0.4, -0.2, 1.3, 0.9)
  data <- data.frame(g = c(1, 1, 1, 1, 2, 2), y = c(x, 0, 1))
  stats <- latents:::.additive_prepare(data, "y", "g")
  parameters <- list(means = matrix(0.3), between = matrix(0.5),
                     within = 0.8, weights = 1)
  integral <- stats::integrate(\(b) vapply(b, \(v) {
    prod(stats::dnorm(x, v, sqrt(0.8))) * stats::dnorm(v, 0.3, sqrt(0.5))
  }, numeric(1)), -Inf, Inf, rel.tol = 1e-11, abs.tol = 0)
  # 1e-9: integrate() is asked for 1e-11 relative on the integral, and its
  # error estimate is not a bound.
  expect_reference(unname(latents:::.additive_log_density(stats, parameters)[1L, 1L]),
                   log(integral$value), reference = "stats",
                   covers = "multilpa()", quantity = "one-indicator group log density",
                   tolerance = 1e-9, measure = "relative",
                   note = "Package side: the internal density (.additive_log_density) of multilpa(family = 'additive'); reference: stats::integrate() over the group intercept.")
})
