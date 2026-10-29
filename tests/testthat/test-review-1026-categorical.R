# Review 2026-10-26: categorical, ordinal, count, local-dependence, direct
# effect, bivariate-residual and sampling-weight mathematics. Calibration
# against hand computation, numDeriv; invariants. The poLCA and Rmpfr
# comparisons are in equivalence/test-review-1026-oracles.R.
# The comparisons with independent references are in tests/equivalence/test-review-1026-categorical.R.

# More than 1.5 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

review_two_class_items <- function(n, seed) {
  withr::local_seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(.4, .6))
  probabilities <- list(rbind(c(.8, .2), c(.2, .8)),
                        rbind(c(.7, .2, .1), c(.1, .3, .6)),
                        rbind(c(.6, .4), c(.3, .7)),
                        rbind(c(.5, .3, .2), c(.2, .2, .6)))
  items <- lapply(probabilities, function(p) vapply(class, function(k)
    sample(seq_len(ncol(p)), 1L, prob = p[k, ]), integer(1)))
  stats::setNames(as.data.frame(items), paste0("i", 1:4))
}

test_that("Poisson log-ratio offset recombines to the hand-computed log pmf", {
  withr::local_seed(1)
  n <- 50
  y <- stats::rpois(n, 4)
  y[c(3, 7)] <- NA
  o <- sample(1:4, n, TRUE)
  o[5] <- NA
  extra <- .latents_set_count_model(
    .latents_prepare_extra(data.frame(o = o, y = y), "o", "y", missing = "fiml"))
  parameters <- list(ordinal_intercepts = list(c(.2, -.3, .5)),
                     ordinal_locations = matrix(c(.4, -.6, 0), 3),
                     count_means = matrix(c(2, 4.5, 9), 3))
  plain <- .latents_extra_log_density(extra, parameters, n, 3)
  offset <- .latents_extra_log_density(extra, parameters, n, 3, offset = TRUE)
  ordinal <- .latents_ordinal_log_probabilities(parameters$ordinal_intercepts[[1]],
                                                parameters$ordinal_locations[, 1])
  hand <- t(vapply(seq_len(n), function(i) {
    (if (is.na(o[i])) 0 else ordinal[, o[i]]) +
      (if (is.na(y[i])) 0 else y[i] * log(parameters$count_means[, 1]) -
         parameters$count_means[, 1] - lgamma(y[i] + 1))
  }, numeric(3)))
  expect_equal(plain, hand, tolerance = 1e-12)
  expect_equal(offset$log_density + offset$offset, plain, tolerance = 1e-12)
  posterior <- function(log_density) {
    joint <- log_density + matrix(log(c(.2, .3, .5)), n, 3, byrow = TRUE)
    weight <- exp(joint - apply(joint, 1L, max))
    weight / rowSums(weight)
  }
  expect_equal(posterior(offset$log_density), posterior(plain), tolerance = 1e-14)
})

test_that("categorical log density is equivariant under a relabelling of profiles", {
  withr::local_seed(11)
  codes <- matrix(sample(1:3, 60, TRUE), 20, 3)
  probabilities <- lapply(1:3, function(j) {
    p <- matrix(stats::runif(9), 3, 3)
    p / rowSums(p)
  })
  permutation <- c(3L, 1L, 2L)
  original <- .multilpa_categorical_log_density(codes, probabilities)
  permuted <- .multilpa_categorical_log_density(
    codes, lapply(probabilities, function(p) p[permutation, ]))
  expect_equal(permuted, original[, permutation], tolerance = 1e-15)
})

test_that("bounded multinomial shares satisfy the Karush-Kuhn-Tucker conditions", {
  counts <- c(0, 1e-3, 5, 2, 0)
  share <- .multilpa_bound_probabilities(counts, 0.05)
  expect_equal(sum(share), 1, tolerance = 1e-15)
  expect_true(all(share >= 0.05))
  free <- share > 0.05 + 1e-12
  multiplier <- counts[free] / share[free]
  expect_equal(diff(range(multiplier)), 0, tolerance = 1e-12)
  expect_true(all(counts[!free] / 0.05 <= multiplier[1] + 1e-12))
  expect_equal(share[free], c(5, 2) / 7 * (1 - 3 * 0.05), tolerance = 1e-15)
})

test_that("bivariate residual chi-squares match the hand-computed Pearson tables", {
  data <- review_two_class_items(300, 4)
  fit <- lca(data, vars = names(data), n_classes = 2, n_starts = 2, seed = 1)
  posterior <- fit$subject_posteriors
  responses <- fit$response_probabilities
  overall <- get_results(fit, "residuals", by = "overall")
  by_profile <- get_results(fit, "residuals", by = "profile")
  hand <- function(first, second, profile = NULL) {
    weight <- if (is.null(profile)) rep(1, nrow(data)) else posterior[, profile]
    observed <- stats::xtabs(weight ~ factor(data[[first]]) + factor(data[[second]]))
    size <- if (is.null(profile)) colSums(posterior) else
      replace(numeric(2), profile, sum(weight))
    expected <- Reduce(`+`, lapply(1:2, function(k)
      size[k] * outer(responses[[first]][k, ], responses[[second]][k, ])))
    sum((observed - expected)^2 / expected)
  }
  pick <- function(table, first, second, label) {
    subset(table, indicator_1 == first & indicator_2 == second & profile == label)
  }
  expect_equal(pick(overall, "i1", "i2", "overall")$statistic, hand("i1", "i2"),
               tolerance = 1e-12)
  expect_equal(pick(overall, "i2", "i4", "overall")$statistic, hand("i2", "i4"),
               tolerance = 1e-12)
  expect_equal(pick(by_profile, "i1", "i3", "profile_2")$statistic,
               hand("i1", "i3", profile = 2), tolerance = 1e-12)
  expect_identical(pick(overall, "i2", "i4", "overall")$df, 4L)
})

test_that("sampling weights: relative within-unit check and scale invariance", {
  unit <- c(1, 1, 2, 2, 3)
  expect_error(.latents_sampling_weights(
    data.frame(w = c(1e-3, 1.0001e-3, 2e-3, 2e-3, 1e-3)), "w", unit, 3),
    class = "latents_bad_weights")
  # 0.3 - 0.2 and 0.1 differ by rounding only.
  accepted <- .latents_sampling_weights(
    data.frame(w = c(0.1, 0.3 - 0.2, 0.2, 0.2, 0.5)), "w", unit, 3)
  expect_equal(accepted, c(0.375, 0.75, 1.875), tolerance = 1e-14)
  expect_equal(.latents_sampling_weights(data.frame(w = c(1e300, 1e300, 1e308, 1e308, 1)),
                                         "w", unit, 3),
               c(3e-8, 3 - 3e-8, 3e-308), tolerance = 1e-12)
  withr::local_seed(5)
  n <- 300
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = stats::rbinom(n, 1, c(.2, .8)[class]),
                     b = stats::rbinom(n, 1, c(.3, .9)[class]),
                     c = stats::rbinom(n, 1, c(.1, .6)[class]),
                     w = sample(1:3, n, TRUE))
  weighted <- lca(data, vars = c("a", "b", "c"), n_classes = 2, weights = "w",
                  n_starts = 2, seed = 1, tol = 1e-12)
  data$w <- data$w * 1e-6
  scaled <- lca(data, vars = c("a", "b", "c"), n_classes = 2, weights = "w",
                n_starts = 2, seed = 1, tol = 1e-12)
  expect_equal(scaled$log_likelihood, weighted$log_likelihood, tolerance = 1e-10)
  expect_equal(scaled$response_probabilities, weighted$response_probabilities,
               tolerance = 1e-10)
  # Normalized weights sum to n, so n / sum(w) times the replicated-data log
  # likelihood is the weighted one.
  data$w <- data$w * 1e6
  replicated <- lca(data[rep(seq_len(n), round(data$w)), ], vars = c("a", "b", "c"),
                    n_classes = 2, n_starts = 2, seed = 1, tol = 1e-12)
  expect_equal(weighted$log_likelihood * sum(round(data$w)) / n,
               replicated$log_likelihood, tolerance = 1e-8)
})

test_that("direct-effect M-step solves the weighted normal equations", {
  withr::local_seed(9)
  n <- 100
  x <- cbind(stats::rnorm(n), stats::rnorm(n))
  x[c(2, 9), 1] <- NA
  z <- cbind(z1 = stats::rnorm(n), z2 = stats::runif(n))
  extra <- list(direct = list(covariates = z, effects = data.frame(
    indicator = c("a", "a", "b"), covariate = c("z1", "z2", "z1"),
    column = c(1L, 1L, 2L), covariate_column = c(1L, 2L, 1L))))
  posterior <- matrix(stats::runif(2 * n), n)
  posterior <- posterior / rowSums(posterior)
  means <- rbind(c(0, 1), c(2, -1))
  variances <- rbind(c(1, 2), c(.5, 1.5))
  effects <- .latents_direct_maximize(x, extra, posterior, means, variances)
  normal <- function(column, covariates) {
    observed <- !is.na(x[, column])
    precision <- sweep(posterior[observed, ], 2L, variances[, column], "/")
    design <- z[observed, covariates, drop = FALSE]
    right <- rowSums(precision * outer(x[observed, column], means[, column], "-"))
    solve(crossprod(design, design * rowSums(precision)), crossprod(design, right))
  }
  expect_equal(effects, c(normal(1, 1:2), normal(2, 1)), tolerance = 1e-12)
  scores <- .latents_direct_scores(x, extra, posterior,
                                   list(means = means, variances = variances,
                                        direct_effects = effects))
  expect_lt(max(abs(scores)), 1e-10)
})

test_that("a direct effect whose covariate is constant where its target is observed is refused", {
  withr::local_seed(4)
  n <- 200
  class <- sample(1:2, n, TRUE)
  data <- data.frame(y1 = stats::rnorm(n, c(0, 3)[class]),
                     y2 = stats::rnorm(n, c(0, 2)[class]), z = 0.1)
  data$z[1:30] <- stats::rnorm(30)
  data$y1[1:30] <- NA
  # Before the repair this fitted y1:z = -3.6e13 and class means of 3.6e12.
  expect_error(lpa(data, vars = c("y1", "y2"), n_profiles = 2, direct_effects = "y1 ~ z",
                   missing = "fiml", n_starts = 2, seed = 1),
               class = "latents_bad_data")
  # The same covariate identifies an effect on the fully observed target.
  fit <- lpa(data, vars = c("y1", "y2"), n_profiles = 2, direct_effects = "y2 ~ z",
             missing = "fiml", n_starts = 2, seed = 1)
  expect_lt(abs(unname(fit$direct_effects)), 10)
})

test_that("local-dependence table labels a joint category nobody gave", {
  withr::local_seed(7)
  n <- 300
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = stats::rbinom(n, 1, c(.2, .8)[class]),
                     b = stats::rbinom(n, 1, c(.3, .9)[class]),
                     c = sample(1:3, n, TRUE),
                     e = stats::rbinom(n, 1, c(.2, .7)[class]))
  data$c[data$a == 1 & data$c == 3] <- 2
  fit <- lca(data, vars = c("a", "b", "c", "e"), n_classes = 2, local_dependence = "a:c",
             n_starts = 2, seed = 1)
  table <- get_results(fit, "local_dependence")
  expect_false(anyNA(table$category_1))
  expect_false(anyNA(table$independent))
  unseen <- subset(table, category_1 == "1" & category_2 == "3")
  expect_identical(nrow(unseen), 2L)
  # Each profile's joint distribution sums to one, and so does its product of
  # marginals.
  expect_equal(as.vector(tapply(table$probability, table$profile, sum)), c(1, 1),
               tolerance = 1e-12)
  expect_equal(as.vector(tapply(table$independent, table$profile, sum)), c(1, 1),
               tolerance = 1e-12)
})

test_that("huge native ordinal values survive fitting, prediction and redrawing", {
  values <- 1e17 + c(0, 16, 32)
  withr::local_seed(6)
  n <- 200
  class <- sample(1:2, n, TRUE)
  code <- vapply(class, function(k) sample(1:3, 1L, prob = if (k == 1) c(.7, .2, .1) else
    c(.1, .3, .6)), integer(1))
  y <- stats::rnorm(n, c(0, 2)[class])
  native <- lpa(data.frame(o = values[code], y = y), vars = c("o", "y"), ordinal = "o",
                n_profiles = 2, n_starts = 2, seed = 1)
  coded <- lpa(data.frame(o = code, y = y), vars = c("o", "y"), ordinal = "o",
               n_profiles = 2, n_starts = 2, seed = 1)
  expect_equal(native$log_likelihood, coded$log_likelihood, tolerance = 1e-12)
  expect_identical(lengths(native$extra_data$ordinal_levels), c(o = 3L))
  expect_equal(predict(native, newdata = data.frame(o = values[code], y = y)[1:5, ]),
               predict(coded, newdata = data.frame(o = code, y = y)[1:5, ]),
               tolerance = 1e-12)
  drawn <- .latents_draw_extra_columns(native$extra_data)
  expect_identical(.latents_as_column_type(drawn$o, values[code]), values[code])
  expect_identical(.latents_as_column_type(factor(c("2", "1")), c(1, 2, 2)), c(2, 1))
})

test_that("a probability on its upper bound has no error, whatever the category order", {
  withr::local_seed(26)
  n <- 400L
  class <- sample(1:2, n, TRUE)
  draw <- function(p) ifelse(stats::runif(n) < p, "yes", "no")
  # Item `u` is "no" in every member of class 1: one probability of that row
  # sits on the floor and its complement on one minus the floor.
  items <- c("u", "v", "w", "z", "s")
  # The other items separate the classes sharply, so EM drives the floored
  # probability onto its bound rather than stopping in its flat approach.
  data <- data.frame(unit = seq_len(n),
                     u = ifelse(class == 1L, "no", draw(0.7)),
                     v = draw(c(0.05, 0.95)[class]), w = draw(c(0.1, 0.9)[class]),
                     z = draw(c(0.05, 0.9)[class]), s = draw(c(0.9, 0.1)[class]))
  bound_errors <- function(levels) {
    ordered <- data
    ordered[items] <- lapply(ordered[items], factor, levels = levels)
    fit <- quietly(multilpa(ordered, vars = items, id = "unit",
                            n_profiles = 2, n_group_classes = 1, categorical = items,
                            n_starts = 10, seed = 1, tol = 1e-12, max_iter = 20000))
    responses <- get_results(fit, "responses")
    held <- subset(responses, indicator == "u" &
                     (probability < 1e-6 | probability > 1 - 1e-6))
    list(held = held, others = subset(responses, indicator != "u"))
  }
  no_last <- bound_errors(c("yes", "no"))
  yes_last <- bound_errors(c("no", "yes"))
  expect_identical(nrow(no_last$held), 2L)
  expect_identical(nrow(yes_last$held), 2L)
  # Before the fix one of these orders reported 0, with the interval [1, 1].
  expect_true(all(is.na(no_last$held$probability_standard_error)))
  expect_true(all(is.na(yes_last$held$probability_standard_error)))
  # Every probability off the bound keeps a positive error.
  expect_true(all(no_last$others$probability_standard_error > 0))
  expect_true(all(yes_last$others$probability_standard_error > 0))
})

test_that("standard errors do not depend on where an unused category sits", {
  withr::local_seed(27)
  n <- 600L
  class <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  pick <- function(p) vapply(seq_len(n), function(i) {
    sample(c("low", "mid", "high"), 1L, prob = p[[class[[i]]]])
  }, character(1))
  # Item `u` never takes "high" in class 1: a row whose third category sits on
  # the floor while "low" and "mid" are interior.
  data <- data.frame(unit = seq_len(n),
                     u = pick(list(c(0.7, 0.3, 0), c(0.2, 0.3, 0.5))),
                     v = pick(list(c(0.1, 0.2, 0.7), c(0.7, 0.2, 0.1))),
                     w = pick(list(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))),
                     z = pick(list(c(0.2, 0.6, 0.2), c(0.5, 0.1, 0.4))))
  items <- c("u", "v", "w", "z")
  errors <- function(levels) {
    ordered <- data
    ordered[items] <- lapply(ordered[items], factor, levels = levels)
    fit <- quietly(multilpa(ordered, vars = items, id = "unit", n_profiles = 2,
                            n_group_classes = 1, categorical = items, n_starts = 10,
                            seed = 1, tol = 1e-12, max_iter = 20000))
    responses <- get_results(fit, "responses")
    # Classes in a fixed order, by their probability of "low" on `v`.
    order_by <- subset(responses, indicator == "v" & category == "low")
    rank <- rank(order_by$probability)[match(responses$profile, order_by$profile)]
    responses$class <- rank
    shares <- subset(parameter_inference(fit, boundary = "fix"),
                     level == "profile")
    shares <- shares$standard_error[order(rank[match(
      as.integer(sub("profile_", "", shares$outcome)), responses$profile)])]
    list(responses = responses[order(responses$indicator, responses$category,
                                     responses$class), ], shares = shares)
  }
  last <- errors(c("low", "mid", "high"))
  first <- errors(c("high", "low", "mid"))
  held <- subset(last$responses, indicator == "u" & probability < 1e-6)
  expect_identical(nrow(held), 1L)
  interior <- subset(last$responses, indicator == "u" & class == held$class &
                       probability > 1e-6)
  # Before the fix, with the unused category last, these were NA.
  expect_true(all(is.finite(interior$probability_standard_error)))
  expect_equal(last$responses$probability_standard_error,
               first$responses$probability_standard_error, tolerance = 1e-5)
  expect_equal(last$shares, first$shares, tolerance = 1e-5)
})
