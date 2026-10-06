# Additive model engine: sufficient statistics, exact density, group
# posteriors. References come from helper-additive-dense.R and never use
# the sufficient-statistic formulas under test.

additive_parameters <- list(
  means = rbind(c(-1, 0.5), c(1.4, -0.8)),
  between = rbind(c(0.3, 0.7), c(0.8, 0.2)),
  within = c(0.6, 1.1), weights = c(0.4, 0.6))
additive_sizes <- c(1L, 2L, 3L, 5L, 7L, 4L, 6L, 2L, 8L, 3L, 5L, 4L)
additive_data <- additive_draw(additive_parameters, additive_sizes, seed = 20260930)
additive_vars <- c("y1", "y2")

test_that("compact density equals the full-covariance Cholesky density", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  compact <- .additive_log_density(stats, additive_parameters)
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  dense <- t(vapply(blocks, \(x) vapply(1:2, \(h) additive_dense_density(
    x, additive_parameters$means[h, ], additive_parameters$between[h, ],
    additive_parameters$within), numeric(1)), numeric(2)))
  expect_equal(unname(compact), dense, tolerance = 1e-10)
  expect_identical(rownames(compact), stats$ids)
})

test_that("group posteriors are normalized and match Bayes' rule on dense densities", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  expectation <- .additive_expectation(stats, additive_parameters)
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  joint <- t(vapply(blocks, \(x) vapply(1:2, \(h) additive_dense_density(
    x, additive_parameters$means[h, ], additive_parameters$between[h, ],
    additive_parameters$within), numeric(1)), numeric(2))) +
    rep(log(additive_parameters$weights), each = length(blocks))
  reference <- exp(joint) / rowSums(exp(joint))
  expect_equal(unname(expectation$posterior), reference, tolerance = 1e-10)
  expect_equal(rowSums(expectation$posterior), rep(1, length(blocks)),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(expectation$log_likelihood, sum(log(rowSums(exp(joint)))),
               tolerance = 1e-10)
  expect_equal(sum(expectation$group_log_likelihood), expectation$log_likelihood)
})

test_that("intercept moments match joint-Gaussian conditioning", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  expectation <- .additive_expectation(stats, additive_parameters)
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  differences <- unlist(lapply(seq_along(blocks), \(j) lapply(1:2, \(h) {
    reference <- additive_dense_moments(
      blocks[[j]], additive_parameters$means[h, ],
      additive_parameters$between[h, ], additive_parameters$within)
    c(reference["mean", ] - expectation$conditional_mean[j, , h],
      reference["variance", ] - expectation$conditional_variance[j, , h])
  })))
  expect_lt(max(abs(differences)), 1e-10)
})
