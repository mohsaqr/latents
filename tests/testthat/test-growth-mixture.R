# Growth mixture models: mixture_regression(random = ).

growth_fixture_data <- function(n = 160L, seed = 21L) {
  set.seed(seed)
  times <- 0:5
  class <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  covariance <- matrix(c(1, 0.2, 0.2, 0.1), 2L)
  effects <- t(chol(covariance)) %*% matrix(stats::rnorm(2L * n), 2L)
  out <- data.frame(id = rep(seq_len(n), each = length(times)),
                    time = rep(times, n), x = stats::rnorm(n * length(times)),
                    age = rep(stats::rnorm(n), each = length(times)))
  person_class <- class[out$id]
  out$y <- c(2, 6)[person_class] + c(0.5, -0.4)[person_class] * out$time +
    0.3 * out$x + effects[1L, out$id] + effects[2L, out$id] * out$time +
    stats::rnorm(nrow(out), 0, c(0.8, 0.6)[person_class])
  out[stats::runif(nrow(out)) > 0.15, ]
}

growth_fit <- function(data, ...) {
  mixture_regression(y ~ time + x, data, n_classes = 2, id = "id",
                     class_level = "group", common = ~ x, membership = ~ age,
                     random = ~ 1 + time, n_starts = 2, seed = 1, ...)
}

test_that("estimation coordinates round-trip and EM never decreases", {
  skip_on_cran()
  data <- growth_fixture_data()
  fit <- growth_fit(data, vcov_type = "none")
  theta <- latents:::.growth_pack(fit$spec, fit$params)
  expect_identical(names(theta), latents:::.growth_names(fit$spec))
  expect_identical(length(theta), fit$n_parameters)
  expect_equal(latents:::.growth_pack(fit$spec,
                                      latents:::.growth_unpack(fit$spec, theta, fit$params)),
               theta, tolerance = 1e-9)
  expect_true(all(diff(fit$history) >= -1e-7 * (1 + abs(fit$log_likelihood))))
  expect_true(fit$converged)
  expect_lt(fit$inference$score_norm %||% max(abs(colSums(
    latents:::.growth_scores(fit$spec, fit$stats, fit$params)))), 1e-2)
})

test_that("the fit does not depend on the order of persons in the data", {
  skip_on_cran()
  data <- growth_fixture_data()
  fit <- growth_fit(data, vcov_type = "none")
  set.seed(9)
  persons <- sample(unique(data$id))
  shuffled <- data[order(match(data$id, persons), data$time), ]
  refit <- growth_fit(shuffled, vcov_type = "none")
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-7)
  ours <- get_results(fit, "classes", vcov_type = "observed")
  theirs <- get_results(refit, "classes", vcov_type = "observed")
  expect_equal(sort(theirs$share), sort(ours$share), tolerance = 1e-5)
})

test_that("integer sampling weights equal duplicated persons", {
  skip_on_cran()
  data <- growth_fixture_data(n = 100L)
  set.seed(2)
  person_weight <- sample(1:3, 100L, replace = TRUE)
  data$w <- person_weight[data$id]
  copies <- rep(seq_len(100L), person_weight)
  duplicated_data <- do.call(rbind, lapply(seq_along(copies), function(j) {
    rows <- data[data$id == copies[j], ]
    rows$id <- j
    rows
  }))
  weighted <- growth_fit(data, weights = "w", tol = 1e-12)
  duplicated_fit <- growth_fit(duplicated_data, tol = 1e-12)
  # Weights are normalized to the number of persons.
  expect_equal(weighted$log_likelihood,
               duplicated_fit$log_likelihood * 100 / sum(person_weight), tolerance = 1e-6)
  expect_equal(sort(unlist(weighted$params$beta)), sort(unlist(duplicated_fit$params$beta)),
               tolerance = 1e-4)
  expect_identical(weighted$inference$vcov_type, "robust")
})

test_that("growth tables are tidy and internally consistent", {
  skip_on_cran()
  data <- growth_fixture_data()
  fit <- growth_fit(data)
  classes <- get_results(fit, "classes")
  expect_equal(sum(classes$share), 1, tolerance = 1e-10)
  expect_equal(sum(classes$expected_persons), length(unique(data$id)), tolerance = 1e-8)
  expect_identical(sum(classes$assigned_persons), length(unique(data$id)))
  assignments <- get_results(fit, "assignments")
  posterior_columns <- grep("^posterior_", names(assignments))
  expect_equal(rowSums(assignments[posterior_columns]), rep(1, nrow(assignments)),
               tolerance = 1e-10)
  coefficients <- get_results(fit, "coefficients")
  expect_setequal(coefficients$class, c("class_1", "class_2", "common"))
  random <- get_results(fit, "random")
  expect_identical(nrow(random), 2L * 5L)
  expect_true(all(random$conf_low[random$parameter == "correlation"] > -1 &
                    random$conf_high[random$parameter == "correlation"] < 1))
  expect_true(all(random$conf_low[random$parameter == "variance"] > 0))
  trajectories <- get_results(fit, "trajectories")
  expect_identical(nrow(trajectories), 2L * 60L)
  expect_true(all(trajectories$conf_low <= trajectories$estimate &
                    trajectories$estimate <= trajectories$conf_high))
  individual <- get_results(fit, "individual")
  expect_identical(nrow(individual), nrow(data))
  expect_named(individual, c("id", "time", "class", "observed", "class_mean", "predicted"))
  fit_row <- get_results(fit, "fit")
  expect_identical(nrow(fit_row), 1L)
  expect_identical(fit_row$n_parameters, fit$n_parameters)
  expect_named(get_results(fit, "all"), latents:::.growth_tables())
  expect_identical(as.data.frame(fit), coefficients)
  expect_s3_class(summary(fit), "summary_latents_growth_mixture")
  expect_output(print(fit), "Growth mixture model")
  expect_identical(nobs(fit), length(unique(data$id)))
})

test_that("a proportional covariance is the shared matrix times class scales", {
  skip_on_cran()
  data <- growth_fixture_data()
  fit <- growth_fit(data, random_covariance = "proportional", vcov_type = "none")
  covariances <- latents:::.growth_covariances(fit$spec, fit$params)
  expect_identical(fit$params$random_scale[2L], 1)
  expect_equal(covariances[[2L]], fit$params$random_covariance[[1L]], tolerance = 1e-12)
  expect_equal(covariances[[1L]] / covariances[[2L]],
               matrix(fit$params$random_scale[1L]^2, 2L, 2L), tolerance = 1e-10)
})

test_that("growth mixture requests outside the model are refused by class", {
  data <- growth_fixture_data(n = 60L)
  expect_error(mixture_regression(y ~ time, data, 2, random = ~ 1),
               class = "latents_bad_argument")
  expect_error(mixture_regression(y ~ time, data, 2, id = "id", random = ~ 1),
               class = "latents_bad_argument")
  data$k <- stats::rpois(nrow(data), 3)
  expect_error(mixture_regression(k ~ time, data, 2, family = "poisson", id = "id",
                                  class_level = "group", random = ~ 1),
               class = "latents_bad_argument")
  expect_error(mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                  random = ~ 1, random_diagonal = NA),
               class = "latents_bad_argument")
  expect_error(get_results(mixture_regression(y ~ time, data, 2, id = "id",
                                              class_level = "group", random = ~ 1,
                                              n_starts = 1, seed = 1),
                           "trajectories", time = "missing_column"),
               class = "latents_bad_argument")
})

test_that("every growth plot draws", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data <- growth_fixture_data(n = 80L)
  fit <- growth_fit(data)
  invisible(lapply(c("trajectories", "individuals", "coefficients", "random",
                      "classification", "posteriors"), function(w) {
    expect_s3_class(plot(fit, what = w), "ggplot")
  }))
  expect_s3_class(plot(fit, facet = TRUE, spread = 0, max_persons = 0), "ggplot")
  expect_error(plot(fit, spread = 1))
  intercept_only <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                       random = ~ 1, n_starts = 1, seed = 1)
  expect_s3_class(plot(intercept_only, what = "random"), "ggplot")
})

test_that("growth estimates recover the generating trajectories", {
  skip_on_cran()
  data <- growth_fixture_data(n = 400L, seed = 7L)
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", random = ~ 1 + time,
                            n_starts = 3, seed = 1)
  coefficients <- get_results(fit, "coefficients")
  intercepts <- sort(coefficients$estimate[coefficients$term == "(Intercept)"])
  expect_equal(intercepts, c(2, 6), tolerance = 0.1)
  shares <- sort(get_results(fit, "classes")$share)
  expect_equal(shares, c(0.4, 0.6), tolerance = 0.12)
})

test_that("variable names stand in for one-sided formulas", {
  skip_on_cran()
  data <- growth_fixture_data(n = 100L)
  by_formula <- mixture_regression(y ~ time + x, data, n_classes = 2, id = "id",
                                   class_level = "group", common = ~ x,
                                   membership = ~ age, random = ~ 1 + time,
                                   n_starts = 1, seed = 1, vcov_type = "none")
  by_name <- mixture_regression(y ~ time + x, data, n_classes = 2, id = "id",
                                class_level = "group", common = "x",
                                membership = "age", random = "time",
                                n_starts = 1, seed = 1, vcov_type = "none")
  expect_identical(by_name$log_likelihood, by_formula$log_likelihood)
  expect_identical(coef(by_name), coef(by_formula))
  intercept_only <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                       random = "intercept", n_starts = 1, seed = 1,
                                       vcov_type = "none")
  expect_identical(colnames(intercept_only$spec$random_design), "(Intercept)")
  expect_error(mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                  random = 3), class = "latents_bad_argument")
  expect_error(mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                  random = ""), class = "latents_bad_argument")
})

test_that("compare_models() compares fits to the same data in one tidy table", {
  skip_on_cran()
  data <- growth_fixture_data(n = 100L)
  trajectories <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                     n_starts = 1, seed = 1)
  growth <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                               random = "time", n_starts = 1, seed = 1)
  comparison <- compare_models(trajectories = trajectories, growth = growth)
  expect_s3_class(comparison, "latents_comparison")
  expect_s3_class(comparison, "data.frame")
  expect_identical(comparison$model, c("trajectories", "growth"))
  expect_identical(comparison$random, c("none", "intercept + time"))
  expect_equal(sum(comparison$bic_weight), 1, tolerance = 1e-12)
  expect_identical(min(comparison$delta_bic), 0)
  expect_equal(comparison$log_likelihood,
               c(trajectories$log_likelihood, growth$log_likelihood))
  # The growth mixture model nests the trajectory model: never a lower likelihood.
  expect_gte(comparison$log_likelihood[2L], comparison$log_likelihood[1L] - 1e-6)
  expect_output(print(comparison), "Weight")
  other <- data
  other$y <- other$y + 1
  shifted <- mixture_regression(y ~ time, other, 2, id = "id", class_level = "group",
                                n_starts = 1, seed = 1)
  expect_error(compare_models(trajectories, shifted), class = "latents_incomparable_models")
  expect_error(compare_models(trajectories), class = "latents_bad_argument")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(plot(comparison), "ggplot")
  }
})

test_that("result tables print readably and stay numeric", {
  skip_on_cran()
  data <- growth_fixture_data(n = 100L)
  fit <- growth_fit(data)
  coefficients <- get_results(fit, "coefficients")
  expect_s3_class(coefficients, "latents_table")
  expect_true(is.numeric(coefficients$estimate) && is.numeric(coefficients$p_value))
  printed <- utils::capture.output(print(coefficients))
  expect_true(any(grepl("95% CI", printed, fixed = TRUE)))
  expect_true(any(grepl("<.001", printed, fixed = TRUE)))
  expect_false(any(grepl("e-0", printed, fixed = TRUE)))
  plain <- as.data.frame(coefficients)
  expect_identical(class(plain), "data.frame")
  expect_identical(plain$estimate, coefficients$estimate)
})

test_that("a degenerate random-effect covariance warns and withholds its errors", {
  skip_on_cran()
  # Random intercepts with no variance in the generating model: the fitted
  # variance runs to its floor.
  set.seed(4)
  n <- 120L
  data <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(0:4, n))
  data$y <- rep(c(0, 3), n / 2)[data$id] + stats::rnorm(nrow(data), 0, 0.5)
  expect_warning(
    fit <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                              random = "intercept", n_starts = 2, seed = 1),
    class = "latents_random_boundary")
  expect_true(get_results(fit, "fit")$random_boundary)
  random <- get_results(fit, "random")
  expect_true(any(random$boundary))
  expect_true(all(is.na(random$std_error[random$boundary])))
  expect_true(all(is.na(random$conf_low[random$boundary])))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(plot(fit, what = "random"), "ggplot")
  }
})

test_that("trajectory models without random effects share the trajectory views", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data <- growth_fixture_data(n = 80L)
  trajectories <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                                     n_starts = 1, seed = 1)
  invisible(lapply(c("trajectories", "individuals", "coefficients", "classification"),
                   function(w) expect_s3_class(plot(trajectories, what = w), "ggplot")))
  expect_s3_class(plot(trajectories), "ggplot")
  # Without random effects a person's own curve is their class trajectory.
  individual <- latents:::.growth_individual_table(trajectories)
  expect_identical(individual$predicted, individual$class_mean)
  # Binomial trajectories are drawn on the probability scale.
  set.seed(2)
  data$k <- stats::rbinom(nrow(data), 5, stats::plogis(rep(c(-1, 1), 40)[data$id] +
                                                         0.2 * data$time))
  binomial <- mixture_regression(cbind(k, 5 - k) ~ time, data, 2, family = "binomial",
                                 id = "id", class_level = "group", n_starts = 1, seed = 1)
  band <- latents:::.growth_trajectory_table(binomial, 0.95,
                                             latents:::.trajectory_inference(binomial))
  expect_true(all(band$conf_low >= 0 & band$conf_high <= 1 &
                    band$conf_low <= band$estimate & band$estimate <= band$conf_high))
  expect_s3_class(plot(binomial), "ggplot")
  row_level <- mixture_regression(y ~ time, data, 2, n_starts = 1, seed = 1)
  expect_error(plot(row_level, what = "trajectories"), class = "latents_bad_argument")
})

test_that("a degenerate random intercept is advised away, not re-suggested", {
  skip_on_cran()
  set.seed(4)
  n <- 120L
  data <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(0:4, n))
  data$y <- rep(c(0, 3), n / 2)[data$id] + stats::rnorm(nrow(data), 0, 0.5)
  warning_text <- tryCatch(
    mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
                       random = "intercept", n_starts = 2, seed = 1),
    latents_random_boundary = function(w) conditionMessage(w))
  expect_match(warning_text, "drop `random`")
  expect_no_match(warning_text, 'random = "intercept"', fixed = TRUE)
})

test_that("growth_scores has the documented structure", {
  expect_identical(names(growth_scores),
                   c("student", "wave", "score", "motivation", "trajectory"))
  expect_identical(nrow(growth_scores), 1844L)
  expect_identical(length(unique(growth_scores$student)), 300L)
  expect_identical(levels(growth_scores$trajectory), c("improving", "stable", "declining"))
  expect_false(anyNA(growth_scores))
  constant <- tapply(growth_scores$motivation, growth_scores$student,
                     function(v) length(unique(v)) == 1L)
  expect_true(all(constant))
})

test_that("recovery compares the classes with a known classification", {
  skip_on_cran()
  fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3, id = "student",
                            class_level = "group", random = "wave",
                            random_covariance = "equal", n_starts = 3, seed = 1)
  recovery <- get_results(fit, "recovery", data = growth_scores, truth = "trajectory")
  expect_identical(sum(recovery$n), 300L)
  expect_equal(as.vector(tapply(recovery$share, recovery$assigned, sum)), rep(1, 3),
               tolerance = 1e-12)
  # Each class is dominated by one generating kind.
  expect_true(all(tapply(recovery$share, recovery$assigned, max) > 0.9))
  random <- get_results(fit, "random")
  expect_identical(unique(random$class), "all")
  expect_error(get_results(fit, "recovery"), class = "latents_bad_argument")
})

test_that("simulation draws from the fitted growth model", {
  skip_on_cran()
  data <- growth_fixture_data(n = 100L)
  fit <- growth_fit(data, vcov_type = "none")
  set.seed(42)
  before <- .Random.seed
  draws <- simulate(fit, nsim = 300, seed = 7)
  expect_identical(.Random.seed, before)
  expect_identical(dim(draws), c(nrow(data), 300L))
  expect_identical(simulate(fit, nsim = 2, seed = 7), simulate(fit, nsim = 2, seed = 7))
  # The average draw approaches the model-implied marginal mean of each row.
  spec <- fit$spec
  prior <- exp(latents:::.mixture_log_softmax(spec$w, fit$params$gamma))
  design <- cbind(spec$x, spec$z)
  implied <- vapply(seq_len(spec$n), function(r) {
    sum(prior[spec$group_index[r], ] * vapply(seq_len(spec$n_classes), function(k) {
      sum(design[r, ] * latents:::.growth_coefficients(fit$params, k))
    }, numeric(1)))
  }, numeric(1))
  expect_lt(mean(abs(rowMeans(draws) - implied)), 0.25)
})

test_that("growth fits enumerate over classes and support the bootstrap test", {
  skip_on_cran()
  data <- growth_fixture_data(n = 60L)
  enumeration <- enumerate_regressions(y ~ time, data, n_classes = 1:2, id = "id",
                                       class_level = "group", random = "intercept",
                                       random_covariance = "equal", n_starts = 1,
                                       bootstrap = 2, bootstrap_starts = 1, seed = 1)
  table <- as.data.frame(enumeration)
  expect_identical(table$n_classes, 1:2)
  expect_true(all(is.finite(table$bic)))
  expect_identical(sum(table$best_bic), 1L)
  expect_true(is.finite(table$blrt_statistic[2L]))
  expect_true(table$blrt_p_value[2L] > 0 && table$blrt_p_value[2L] <= 1)
})

test_that("distal outcomes of trajectory models use the three-step correction", {
  skip_on_cran()
  fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3, id = "student",
                            class_level = "group", random = "wave",
                            random_covariance = "equal", n_starts = 3, seed = 1)
  means <- three_step(fit, data = growth_scores, outcome = "motivation")
  expect_s3_class(means, "data.frame")
  expect_identical(unique(means$level), "persons")
  expect_identical(nrow(means), 3L)
  # Motivation was generated to raise the improving class: its class mean is
  # the largest, the declining class's the smallest.
  recovery <- as.data.frame(get_results(fit, "recovery", data = growth_scores,
                                        truth = "trajectory"))
  improving <- recovery$assigned[recovery$trajectory == "improving" &
                                   recovery$share > 0.5]
  declining <- recovery$assigned[recovery$trajectory == "declining" &
                                   recovery$share > 0.5]
  class_number <- function(label) as.integer(sub("class_", "", label))
  expect_identical(means$class[which.max(means$estimate)], class_number(improving))
  expect_identical(means$class[which.min(means$estimate)], class_number(declining))
  pairs <- three_step(fit, data = growth_scores, outcome = "motivation", contrast = "pairs")
  expect_identical(nrow(pairs), 3L)
  expect_error(three_step(fit, data = growth_scores, outcome = "score"),
               class = "latents_bad_outcome")
  expect_error(three_step(fit, data = growth_scores, outcome = "wave"),
               class = "latents_bad_outcome")
  with_membership <- mixture_regression(score ~ wave, growth_scores, n_classes = 2,
                                        id = "student", class_level = "group",
                                        random = "intercept", membership = "motivation",
                                        n_starts = 1, seed = 1)
  expect_error(three_step(with_membership, data = growth_scores, outcome = "motivation"),
               class = "latents_unsupported_three_step")
  expect_error(r3step(fit, data = growth_scores, covariates = "motivation"),
               class = "latents_unsupported_three_step")
  trajectories <- mixture_regression(score ~ wave, growth_scores, n_classes = 3,
                                     id = "student", class_level = "group",
                                     n_starts = 2, seed = 1)
  expect_identical(nrow(three_step(trajectories, data = growth_scores,
                                   outcome = "motivation")), 3L)
})

test_that("prediction refuses clearly", {
  data <- growth_fixture_data(n = 60L)
  fit <- growth_fit(data, vcov_type = "none")
  expect_error(predict(fit), class = "latents_unsupported_prediction")
})
