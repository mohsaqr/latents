# Sampling weights: pseudo maximum likelihood over the independent units.

weighted_rows <- function(seed = 1, n = 100) {
  set.seed(seed)
  data.frame(a = c(rnorm(0.6 * n), rnorm(0.4 * n, 3)),
             b = c(rnorm(0.6 * n), rnorm(0.4 * n, 2)),
             w = sample(1:3, n, TRUE))
}

# Repeat every unit as many times as its integer weight, as separate units.
duplicate_units <- function(data, id, weight) {
  rows <- split(seq_len(nrow(data)), data[[id]])
  copies <- lapply(names(rows), function(unit) {
    taken <- rows[[unit]]
    times <- data[[weight]][taken[1L]]
    block <- data[rep(taken, times), , drop = FALSE]
    block[[id]] <- paste(unit, rep(seq_len(times), each = length(taken)), sep = "_")
    block
  })
  do.call(rbind, copies)
}

# The pseudo log likelihood is on the scale of the units, so it is the
# duplicated data's log likelihood times units / sum of weights.
on_duplicate_scale <- function(fit, weights) {
  fit$log_likelihood * sum(weights) / length(weights)
}

numeric_gradient <- function(objective, theta, step = 1e-6) {
  vapply(seq_along(theta), function(index) {
    up <- theta
    down <- theta
    up[index] <- up[index] + step
    down[index] <- down[index] - step
    (objective(up) - objective(down)) / (2 * step)
  }, numeric(1))
}

test_that("integer weights reproduce the duplicated data under every structure", {
  data <- weighted_rows()
  duplicated <- data[rep(seq_len(nrow(data)), data$w), ]
  vapply(c("VVI", "EEI", "VVV", "VEI"), function(structure) {
    weighted <- quietly(lpa(data, c("a", "b"), 2, weights = "w", model = structure,
                            n_starts = 1, seed = 1, tol = 1e-12),
                        "latents_single_level")
    repeated <- quietly(lpa(duplicated, c("a", "b"), 2, model = structure,
                            n_starts = 1, seed = 1, tol = 1e-12),
                        "latents_single_level")
    expect_equal(on_duplicate_scale(weighted, data$w), repeated$log_likelihood,
                 tolerance = 1e-10)
    expect_equal(weighted$means, repeated$means, tolerance = 1e-6)
    expect_equal(weighted$variances, repeated$variances, tolerance = 1e-6)
    TRUE
  }, logical(1))
})

test_that("integer weights reproduce duplication for categorical and FIML fits", {
  data <- weighted_rows(2)
  data$c <- factor(sample(c("x", "y", "z"), nrow(data), TRUE))
  data$a[c(3, 9, 40)] <- NA
  duplicated <- data[rep(seq_len(nrow(data)), data$w), ]
  weighted <- quietly(lpa(data, c("a", "b", "c"), 2, categorical = "c",
                          missing = "fiml", weights = "w", n_starts = 1, seed = 1,
                          tol = 1e-12), "latents_single_level")
  repeated <- quietly(lpa(duplicated, c("a", "b", "c"), 2, categorical = "c",
                          missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12),
                      "latents_single_level")
  expect_equal(on_duplicate_scale(weighted, data$w), repeated$log_likelihood,
               tolerance = 1e-10)
  expect_equal(weighted$response_probabilities, repeated$response_probabilities,
               tolerance = 1e-6)
})

test_that("two-level weights weight groups, and reproduce duplicated groups", {
  set.seed(3)
  data <- data.frame(id = rep(seq_len(30), each = 5), a = rnorm(150), b = rnorm(150))
  data$a <- data$a + rep(sample(c(0, 3), 30, TRUE), each = 5) * rbinom(150, 1, 0.7)
  data$w <- rep(sample(1:3, 30, TRUE), each = 5)
  unit_weights <- data$w[!duplicated(data$id)]
  weighted <- multilpa(data, c("a", "b"), "id", 2, 2, weights = "w", n_starts = 1,
                       seed = 1, tol = 1e-12)
  repeated <- multilpa(duplicate_units(data, "id", "w"), c("a", "b"), "id", 2, 2,
                       n_starts = 1, seed = 1, tol = 1e-12)
  expect_equal(on_duplicate_scale(weighted, unit_weights), repeated$log_likelihood,
               tolerance = 1e-10)
  expect_equal(sort(weighted$means), sort(repeated$means), tolerance = 1e-5)
  expect_equal(sum(weighted$sampling_weights), 30)
  # Classification stays per unit: each row's posterior is a probability.
  expect_equal(unname(rowSums(weighted$subject_posteriors)), rep(1, 150))
  # Weighted counts: the scaled weights sum to the 30 groups, each of 5 rows.
  counts <- get_results(weighted, "counts")
  expect_equal(sum(subset(counts, level == "groups")$effective_count), 30)
  expect_equal(sum(subset(counts, level == "individuals")$effective_count), 150)
})

test_that("membership covariates take weights and reproduce duplication", {
  set.seed(4)
  n <- 200
  z <- rnorm(n)
  membership <- rbinom(n, 1, plogis(1.2 * z))
  data <- data.frame(a = rnorm(n, 3 * membership), b = rnorm(n, 2 * membership),
                     z = z, w = sample(1:3, n, TRUE))
  duplicated <- data[rep(seq_len(n), data$w), ]
  weighted <- quietly(lpa(data, c("a", "b"), 2, profile_covariates = "z",
                          weights = "w", n_starts = 1, seed = 1, tol = 1e-12),
                      "latents_single_level")
  repeated <- quietly(lpa(duplicated, c("a", "b"), 2, profile_covariates = "z",
                          n_starts = 1, seed = 1, tol = 1e-12),
                      "latents_single_level")
  expect_equal(on_duplicate_scale(weighted, data$w), repeated$log_likelihood,
               tolerance = 1e-10)
  expect_equal(weighted$profile_coefficients, repeated$profile_coefficients,
               tolerance = 1e-5)
  expect_identical(get_results(weighted, "model")$weights, "w")
})

test_that("unit weights change nothing and weights are scale invariant", {
  data <- weighted_rows(5)
  data$one <- 1
  data$scaled <- data$w * 7.5
  plain <- quietly(lpa(data, c("a", "b"), 2, n_starts = 1, seed = 1, tol = 1e-12),
                   "latents_single_level")
  unit <- quietly(lpa(data, c("a", "b"), 2, weights = "one", n_starts = 1, seed = 1,
                      tol = 1e-12), "latents_single_level")
  expect_equal(unit$log_likelihood, plain$log_likelihood, tolerance = 1e-12)
  expect_equal(parameter_inference(unit)$standard_error,
               parameter_inference(plain, vcov_type = "robust")$standard_error,
               tolerance = 1e-8)
  once <- quietly(lpa(data, c("a", "b"), 2, weights = "w", n_starts = 1, seed = 1,
                      tol = 1e-12), "latents_single_level")
  scaled <- quietly(lpa(data, c("a", "b"), 2, weights = "scaled", n_starts = 1,
                        seed = 1, tol = 1e-12), "latents_single_level")
  expect_equal(scaled$log_likelihood, once$log_likelihood, tolerance = 1e-12)
  expect_equal(scaled$bic, once$bic, tolerance = 1e-12)
})

test_that("weighted scores are the gradient of the weighted pseudo likelihood", {
  data <- weighted_rows(6)
  data$w <- runif(nrow(data), 0.3, 3)
  fits <- list(
    quietly(lpa(data, c("a", "b"), 2, weights = "w", model = "VVV", n_starts = 1,
                seed = 1), "latents_single_level"),
    quietly(lpa(data, c("a", "b"), 2, weights = "w", model = "VEI", n_starts = 1,
                seed = 1), "latents_single_level"))
  vapply(fits, function(fit) {
    prepared <- .multilpa_inference_matrix(fit, data)
    set.seed(7)
    theta <- unname(.multilpa_coefficients(fit, "unconstrained"))
    theta <- theta + rnorm(length(theta), 0, 0.05)
    objective <- function(point) {
      -.multilpa_expectation(prepared$x, fit$group_index, .multilpa_decode(point, fit),
                             prepared$codes, unname(fit$sampling_weights))$log_likelihood
    }
    analytic <- .multilpa_score(theta, prepared$x, fit, prepared$codes)
    expect_equal(analytic, numeric_gradient(objective, theta), tolerance = 1e-5)
    by_unit <- .multilpa_group_scores(theta, prepared$x, fit, prepared$codes)
    expect_equal(-colSums(by_unit), analytic, tolerance = 1e-10)
    TRUE
  }, logical(1))
  set.seed(8)
  groups <- data.frame(id = rep(seq_len(60), each = 5))
  driver <- rnorm(60)
  classes <- rbinom(60, 1, plogis(1.5 * driver))
  groups$u <- driver[groups$id]
  groups$a <- rnorm(300, 3 * rbinom(300, 1, ifelse(classes[groups$id] == 1, 0.8, 0.2)))
  groups$b <- rnorm(300)
  groups$w <- rep(runif(60, 0.3, 3), each = 5)
  fit <- multilpa(groups, c("a", "b"), "id", 2, 2, group_covariates = "u",
                  weights = "w", n_starts = 1, seed = 1)
  x <- sweep(as.matrix(groups[c("a", "b")]), 2L, fit$center, "-")
  theta <- .multilpa_cov_encode(fit) + rnorm(length(.multilpa_cov_encode(fit)), 0, 0.05)
  objective <- function(point) {
    pieces <- .multilpa_cov_decode(point, fit)
    -.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                               fit$profile_design, fit$group_design, pieces$beta,
                               pieces$gamma, NULL,
                               unname(fit$sampling_weights))$log_likelihood
  }
  expect_equal(unname(-colSums(.multilpa_cov_group_scores(theta, x, fit))),
               numeric_gradient(objective, theta), tolerance = 1e-5)
})

test_that("weighted fits default to robust errors and refuse the others", {
  data <- weighted_rows(9)
  fit <- quietly(lpa(data, c("a", "b"), 2, weights = "w", n_starts = 1, seed = 1),
                 "latents_single_level")
  inference <- parameter_inference(fit)
  expect_identical(attr(inference, "vcov_type"), "robust")
  expect_true(all(is.finite(inference$standard_error)))
  expect_error(parameter_inference(fit, vcov_type = "observed"),
               class = "latents_unsupported_weights")
  expect_error(vcov(fit, vcov_type = "opg"), class = "latents_unsupported_weights")
  boot <- parameter_inference(fit, method = "bootstrap", iter = 3, n_starts = 1,
                              seed = 1)
  expect_identical(nrow(boot), nrow(inference))
  expect_true("w" %in% names(get_results(fit, "data")))
})

test_that("weights are validated and unsupported pairings refused", {
  data <- weighted_rows(10)
  data$negative <- -data$w
  data$text <- as.character(data$w)
  data$missing_weight <- replace(data$w, 4, NA)
  invalid <- c("not_a_column", "negative", "text", "missing_weight")
  vapply(invalid, function(column) {
    expect_error(quietly(lpa(data, c("a", "b"), 2, weights = column, n_starts = 1),
                         "latents_single_level"),
                 class = "latents_bad_weights")
    TRUE
  }, logical(1))
  grouped <- data.frame(id = rep(seq_len(20), each = 5), a = rnorm(100),
                        b = rnorm(100), w = runif(100))
  expect_error(multilpa(grouped, c("a", "b"), "id", 2, 1, weights = "w", n_starts = 1),
               class = "latents_bad_weights")
  expect_error(quietly(lpa(data, c("a", "b"), 2, weights = "w", noise = TRUE, n_starts = 1),
                       "latents_single_level"),
               class = "latents_unsupported_weights")
  expect_error(quietly(lpa(data, c("a", "b"), 2, weights = "w",
                           prior = prior_control(), n_starts = 1), "latents_single_level"),
               class = "latents_unsupported_weights")
  fit <- quietly(lpa(data, c("a", "b"), 2, weights = "w", n_starts = 1, seed = 1),
                 "latents_single_level")
  smaller <- quietly(lpa(data, c("a", "b"), 1, weights = "w", n_starts = 1), "latents_single_level")
  expect_error(bootstrap_lrt(smaller, fit, iter = 2),
               class = "latents_unsupported_weights")
  expect_error(three_step(fit, data, "a"), class = "latents_unsupported_weights")
  expect_error(r3step(fit, data, "b"), class = "latents_unsupported_weights")
})

test_that("a weighted fit says so when printed", {
  data <- weighted_rows(11)
  fit <- quietly(lpa(data, c("a", "b"), 2, weights = "w", n_starts = 1, seed = 1),
                 "latents_single_level")
  expect_output(print(fit), "Sampling weights `w`.*Kish effective n")
  expect_identical(get_results(fit, "model")$weights, "w")
  plain <- quietly(lpa(data, c("a", "b"), 2, n_starts = 1, seed = 1),
                   "latents_single_level")
  expect_identical(get_results(plain, "model")$weights, NA_character_)
})

# Two-level data with one integer weight per group, and the same data with
# every group repeated that many times.
weighted_groups <- function(seed, n_groups = 50, size = 6) {
  set.seed(seed)
  membership <- rbinom(n_groups, 1, 0.4)
  data <- data.frame(id = rep(seq_len(n_groups), each = size))
  data$a <- 2 * membership[data$id] + rnorm(n_groups, 0, 0.5)[data$id] +
    2 * rbinom(nrow(data), 1, 0.5) + rnorm(nrow(data), 0, 0.6)
  data$b <- membership[data$id] + rnorm(nrow(data))
  data$w <- sample(1:3, n_groups, TRUE)[data$id]
  data
}

test_that("weighted transition fits reproduce duplicated persons", {
  skip_on_cran()
  data <- course_engagement
  set.seed(12)
  persons <- unique(data$student)
  person_weight <- sample(1:3, length(persons), TRUE)
  data$w <- person_weight[match(data$student, persons)]
  repeated <- duplicate_units(data, "student", "w")
  activity <- c("browse", "lectures")
  weighted <- lta(data, activity, "student", n_profiles = 2, time = "sequence",
                  weights = "w", n_starts = 1, seed = 1, tol = 1e-10)
  expect_s3_class(weighted, "multilpa_lta")
  plain <- lta(repeated, activity, "student", n_profiles = 2, time = "sequence",
               n_starts = 1, seed = 1, tol = 1e-10)
  expect_equal(on_duplicate_scale(weighted, person_weight), plain$log_likelihood,
               tolerance = 1e-8)
  covariate <- lta(data, activity, "student", n_profiles = 2, time = "sequence",
                   transition_covariates = "previous_grade", weights = "w",
                   n_starts = 1, seed = 1, tol = 1e-10)
  covariate_plain <- lta(repeated, activity, "student", n_profiles = 2,
                         time = "sequence", transition_covariates = "previous_grade",
                         n_starts = 1, seed = 1, tol = 1e-10)
  expect_equal(on_duplicate_scale(covariate, person_weight),
               covariate_plain$log_likelihood, tolerance = 1e-8)
  expect_equal(covariate$transition_coefficients,
               covariate_plain$transition_coefficients, tolerance = 1e-4)
  expect_output(print(weighted), "Sampling weights `w`")
  expect_identical(attr(parameter_inference(weighted), "vcov_type"), "robust")
  expect_error(parameter_inference(weighted, vcov_type = "observed"),
               class = "latents_unsupported_weights")
})

test_that("weighted transition scores are the gradient of the pseudo likelihood", {
  skip_on_cran()
  data <- course_engagement
  set.seed(13)
  persons <- unique(data$student)
  data$w <- runif(length(persons), 0.3, 3)[match(data$student, persons)]
  fit <- lta(data, c("browse", "lectures"), "student", n_profiles = 2,
             time = "sequence", weights = "w", n_starts = 1, seed = 1)
  theta <- coef(fit) + rnorm(length(coef(fit)), 0, 0.03)
  objective <- function(point) {
    -.lta_expectation(fit$x, fit$codes, fit$layout, fit$designs,
                      .lta_unpack(point, fit), fit$occasion_of_row,
                      unname(fit$sampling_weights))$log_likelihood
  }
  expect_equal(unname(-colSums(.lta_group_scores(theta, fit))),
               numeric_gradient(objective, unname(theta)), tolerance = 1e-5)
})

test_that("weighted group-class families reproduce duplicated groups", {
  skip_on_cran()
  data <- weighted_groups(14)
  unit_weights <- data$w[!duplicated(data$id)]
  repeated <- duplicate_units(data, "id", "w")
  vapply(c("additive", "additive_dispersion", "restricted_cross_level",
           "full_cross_level"), function(family) {
    arguments <- list(vars = c("a", "b"), id = "id", n_group_classes = 2,
                      family = family, n_starts = 1, seed = 1, tol = 1e-12)
    if (endsWith(family, "cross_level")) arguments$n_profiles <- 2
    weighted <- quietly(do.call(multilpa, c(list(data = data, weights = "w"),
                                            arguments)))
    plain <- quietly(do.call(multilpa, c(list(data = repeated), arguments)))
    expect_equal(on_duplicate_scale(weighted, unit_weights), plain$log_likelihood,
                 tolerance = 1e-9)
    expect_equal(sort(weighted$group_probabilities), sort(plain$group_probabilities),
                 tolerance = 1e-6)
    TRUE
  }, logical(1))
})

test_that("weighted additive scores are the gradient and inference is robust", {
  skip_on_cran()
  data <- weighted_groups(15, n_groups = 60)
  data$w <- runif(60, 0.3, 3)[data$id]
  fit <- quietly(multilpa(data, c("a", "b"), "id", n_group_classes = 2,
                          family = "additive", weights = "w", n_starts = 1,
                          seed = 1))
  stats <- fit$sufficient_statistics
  class_names <- names(fit$group_probabilities)
  parameters <- list(means = unname(fit$means), between = unname(fit$between_variances),
                     within = unname(fit$within_variances),
                     weights = unname(fit$group_probabilities))
  set.seed(16)
  theta <- .additive_pack(parameters, fit$structure, class_names, fit$vars)
  theta <- theta + rnorm(length(theta), 0, 0.03)
  unpack <- function(point) .additive_unpack(point, 2L, 2L, fit$structure)
  objective <- function(point) .additive_expectation(stats, unpack(point))$log_likelihood
  scores <- .additive_group_scores(stats, unpack(theta), fit$structure, class_names,
                                   fit$vars)
  expect_equal(unname(colSums(scores)), numeric_gradient(objective, unname(theta)),
               tolerance = 1e-5)
  expect_error(vcov(fit, type = "observed"), class = "latents_unsupported_weights")
})

test_that("weighted mixture regressions reproduce duplication at every nesting", {
  set.seed(17)
  n <- 300
  x <- rnorm(n)
  class <- rbinom(n, 1, 0.4)
  data <- data.frame(x = x, y = ifelse(class == 1, 1 + 2 * x, -1 - x) +
                       rnorm(n, 0, 0.5),
                     w = sample(1:3, n, TRUE), id = rep(seq_len(60), each = 5))
  data$group_weight <- sample(1:3, 60, TRUE)[data$id]
  rows <- data[rep(seq_len(n), data$w), ]
  weighted <- mixture_regression(y ~ x, data, 2, weights = "w", n_starts = 1,
                                 seed = 1, tol = 1e-12)
  plain <- mixture_regression(y ~ x, rows, 2, n_starts = 1, seed = 1, tol = 1e-12)
  expect_equal(on_duplicate_scale(weighted, data$w), plain$log_likelihood,
               tolerance = 1e-10)
  expect_equal(weighted$params$beta, plain$params$beta, tolerance = 1e-6)
  expect_identical(weighted$inference$vcov_type, "robust")
  expect_output(print(weighted), "Sampling weights `w`")
  groups <- duplicate_units(data, "id", "group_weight")
  group_weights <- data$group_weight[!duplicated(data$id)]
  weighted_group <- mixture_regression(y ~ x, data, 2, id = "id",
                                       class_level = "group",
                                       weights = "group_weight", n_starts = 1,
                                       seed = 1, tol = 1e-12, vcov_type = "none")
  plain_group <- mixture_regression(y ~ x, groups, 2, id = "id",
                                    class_level = "group", n_starts = 1, seed = 1,
                                    tol = 1e-12, vcov_type = "none")
  expect_equal(on_duplicate_scale(weighted_group, group_weights),
               plain_group$log_likelihood, tolerance = 1e-9)
  expect_error(mixture_regression(y ~ x, data, 2, weights = "w",
                                  vcov_type = "opg", n_starts = 1),
               class = "latents_unsupported_weights")
  expect_error(enumerate_regressions(y ~ x, data, n_classes = 1:2, bootstrap = 2,
                                     weights = "w", n_starts = 1),
               class = "latents_unsupported_weights")
})
