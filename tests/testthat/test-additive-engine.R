# Additive model engine: sufficient statistics, exact density, group
# posteriors. The comparisons with dense reference densities are in
# tests/equivalence/.

additive_parameters <- list(
  means = rbind(c(-1, 0.5), c(1.4, -0.8)),
  between = rbind(c(0.3, 0.7), c(0.8, 0.2)),
  within = c(0.6, 1.1), weights = c(0.4, 0.6))
additive_sizes <- c(1L, 2L, 3L, 5L, 7L, 4L, 6L, 2L, 8L, 3L, 5L, 4L)
additive_data <- additive_draw(additive_parameters, additive_sizes, seed = 20260930)
additive_vars <- c("y1", "y2")

test_that("sufficient statistics match per-group means and scatter", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  expect_identical(stats$ids, unique(additive_data$group))
  expect_identical(stats$sizes, additive_sizes)
  expect_identical(stats$n_obs, sum(additive_sizes))
  expect_equal(unname(stats$averages),
               do.call(rbind, lapply(blocks, colMeans)), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(unname(stats$scatter),
               do.call(rbind, lapply(blocks, \(x) apply(x, 2L, \(v) sum((v - mean(v))^2)))),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(stats$pooled_scatter, colSums(stats$scatter))
})

test_that("one-indicator density equals numerical integration over the intercept", {
  x <- c(0.4, -0.2, 1.3, 0.9)
  data <- data.frame(g = c(1, 1, 1, 1, 2, 2), y = c(x, 0, 1))
  stats <- .additive_prepare(data, "y", "g")
  parameters <- list(means = matrix(0.3), between = matrix(0.5),
                     within = 0.8, weights = 1)
  integral <- stats::integrate(\(b) vapply(b, \(v) {
    prod(stats::dnorm(x, v, sqrt(0.8))) * stats::dnorm(v, 0.3, sqrt(0.5))
  }, numeric(1)), -Inf, Inf, rel.tol = 1e-11, abs.tol = 0)
  expect_equal(.additive_log_density(stats, parameters)[1L, 1L],
               log(integral$value), tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("marginal intercept moments follow the law of total variance", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  e <- .additive_expectation(stats, additive_parameters)
  # Independent route: E[B^2] mixed over classes, minus the squared mean.
  second <- Reduce(`+`, lapply(1:2, \(h) e$posterior[, h] *
    (e$conditional_variance[, , h] + e$conditional_mean[, , h]^2)))
  first <- Reduce(`+`, lapply(1:2, \(h) e$posterior[, h] * e$conditional_mean[, , h]))
  expect_equal(e$intercept_mean, first, tolerance = 1e-12)
  expect_equal(e$intercept_variance, second - first^2, tolerance = 1e-10)
})

test_that("likelihood is invariant to row order and to class labels", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  base <- .additive_expectation(stats, additive_parameters)
  shuffled <- additive_data[rev(seq_len(nrow(additive_data))), ]
  reordered <- .additive_expectation(
    .additive_prepare(shuffled, additive_vars, "group"), additive_parameters)
  expect_equal(reordered$log_likelihood, base$log_likelihood, tolerance = 1e-12)
  expect_equal(reordered$posterior[stats$ids, ], base$posterior, tolerance = 1e-12)
  # A non-self-inverse 3-cycle on three classes.
  three <- list(means = rbind(c(-1, 0.5), c(1.4, -0.8), c(0, 0)),
                between = rbind(c(0.3, 0.7), c(0.8, 0.2), c(0.1, 0.4)),
                within = c(0.6, 1.1), weights = c(0.2, 0.5, 0.3))
  cycle <- c(2L, 3L, 1L)
  cycled <- list(means = three$means[cycle, ], between = three$between[cycle, ],
                 within = three$within, weights = three$weights[cycle])
  a <- .additive_expectation(stats, three)
  b <- .additive_expectation(stats, cycled)
  expect_equal(b$log_likelihood, a$log_likelihood, tolerance = 1e-12)
  expect_equal(b$posterior, a$posterior[, cycle], tolerance = 1e-12)
  expect_equal(b$conditional_mean, a$conditional_mean[, , cycle], tolerance = 1e-12)
})

test_that("zero between variance reduces to independent Gaussian ratings", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  zero <- additive_parameters
  zero$between[, ] <- 0
  e <- .additive_expectation(stats, zero)
  blocks <- additive_blocks(additive_data, additive_vars, "group")
  reference <- t(vapply(blocks, \(x) vapply(1:2, \(h) sum(stats::dnorm(
    x, rep(zero$means[h, ], each = nrow(x)),
    rep(sqrt(zero$within), each = nrow(x)), log = TRUE)), numeric(1)), numeric(2)))
  expect_equal(unname(e$log_density), reference, tolerance = 1e-10)
  expect_equal(e$conditional_variance, array(0, dim(e$conditional_variance)),
               ignore_attr = TRUE)
})

test_that("singleton groups leave within and between variances on a ridge", {
  singles <- additive_data[!duplicated(additive_data$group), ]
  mixed <- .additive_prepare(additive_data, additive_vars, "group")
  shifted <- additive_parameters
  shifted$within <- shifted$within + 0.05
  shifted$between <- shifted$between - 0.05
  # g01 has one rating: moving variance from between to within leaves its
  # density unchanged, while replicated groups detect the shift.
  expect_identical(mixed$sizes[1L], 1L)
  expect_equal(.additive_log_density(mixed, shifted)["g01", ],
               .additive_log_density(mixed, additive_parameters)["g01", ],
               tolerance = 1e-12)
  expect_gt(abs(.additive_expectation(mixed, shifted)$log_likelihood -
                  .additive_expectation(mixed, additive_parameters)$log_likelihood),
            1e-4)
  expect_error(.additive_prepare(singles, additive_vars, "group"),
               class = "latents_unidentified")
})

test_that("raw-rating contract refusals are classed", {
  bad_missing <- additive_data
  bad_missing$y1[3L] <- NA
  expect_error(.additive_prepare(bad_missing, additive_vars, "group"),
               class = "latents_bad_data")
  bad_factor <- additive_data
  bad_factor$y2 <- factor(bad_factor$y2 > 0)
  expect_error(.additive_prepare(bad_factor, additive_vars, "group"),
               class = "latents_bad_data")
  bad_constant <- transform(additive_data, y2 = 1)
  expect_error(.additive_prepare(bad_constant, additive_vars, "group"),
               class = "latents_bad_data")
  bad_id <- additive_data
  bad_id$group[1L] <- NA
  expect_error(.additive_prepare(bad_id, additive_vars, "group"),
               class = "latents_bad_data")
  expect_error(.additive_prepare(additive_data, c("y1", "zz"), "group"),
               class = "latents_bad_data")
  expect_error(.additive_prepare(additive_data, additive_vars, "y1"),
               class = "latents_bad_data")
})

test_that("malformed parameters are refused by class", {
  stats <- .additive_prepare(additive_data, additive_vars, "group")
  broken <- list(
    negative_between = modifyList(additive_parameters,
                                  list(between = -additive_parameters$between)),
    zero_within = modifyList(additive_parameters, list(within = c(0, 1))),
    weights_off = modifyList(additive_parameters, list(weights = c(0.5, 0.6))),
    wrong_width = modifyList(additive_parameters,
                             list(means = cbind(additive_parameters$means, 0))))
  invisible(lapply(broken, \(p) expect_error(.additive_log_density(stats, p),
                                             class = "latents_bad_start")))
})
