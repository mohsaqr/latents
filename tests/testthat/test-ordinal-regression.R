# Ordinal regression mixtures, mixture_regression(family = "ordinal") (phase
# E7 of validation/ENGINE_DESIGN.md): a proportional-odds (cumulative logit)
# regression within each class under the regression structures. Here: a
# likelihood written with plogis() differences alone, recovery of simulated
# truths, and two categories as the logistic mixture. The comparisons with
# MASS::polr() and numerical derivatives are in tests/equivalence/.

# Rows with one ordinal outcome each, from two classes whose slopes differ in
# sign; categories 1..4 as an ordered factor.
ordinal_rows <- function(seed = 11, n = 2000) {
  latents:::.mixture_with_seed(seed, {
    x <- stats::rnorm(n)
    w <- stats::runif(n)
    class <- 1L + (stats::runif(n) < 0.4)
    cuts <- rbind(c(-1, 0, 1), c(-0.5, 1, 2))
    latent <- c(2, -1.5)[class] * x - 0.5 * w + stats::rlogis(n)
    y <- rowSums(latent > cuts[class, , drop = FALSE]) + 1
    data.frame(y = factor(y, levels = 1:4, labels = c("never", "rarely", "often", "always"),
                          ordered = TRUE),
               x = x, w = w, class = class)
  })
}

# Persons with five ordinal outcomes each (waves 0..4), whose class sets a
# rising or a falling trajectory; categories are the numbers 1..4.
ordinal_persons <- function(seed = 5, persons = 300, waves = 5) {
  latents:::.mixture_with_seed(seed, {
    class <- 1L + (stats::runif(persons) < 0.4)
    data <- data.frame(person = rep(seq_len(persons), each = waves),
                       wave = rep(seq_len(waves) - 1, persons))
    data$class <- class[data$person]
    cuts <- rbind(c(-1, 0.5, 2), c(-2, -0.5, 1))
    latent <- c(0.8, -0.6)[data$class] * data$wave + stats::rlogis(nrow(data))
    data$y <- rowSums(latent > cuts[data$class, , drop = FALSE]) + 1
    data
  })
}

test_that("a simulated two-class ordinal mixture is recovered", {
  skip_on_cran()
  data <- ordinal_rows()
  fit <- mixture_regression(y ~ x + w, data, n_classes = 2, family = "ordinal",
                            common = ~ w, n_starts = 3, seed = 1)
  coefficients <- get_results(fit, "coefficients")
  # One outcome per row identifies the classes only through the slopes, so
  # the estimates are judged against their own standard errors.
  slopes <- subset(coefficients, term %in% c("x", "w"))
  expect_true(all(abs(slopes$estimate - c(2, -1.5, -0.5)) < 3 * slopes$std_error))
  expect_equal(get_results(fit, "classes")$share, c(0.6, 0.4), tolerance = 0.15)
  expect_true(all(diff(fit$params$thresholds) > 0))
  expect_true(all(is.finite(coefficients$std_error)))
  expect_true(all(is.na(subset(coefficients, startsWith(term, "threshold:"))$p_adjusted)))
  # The maximum is at least as likely as the truth.
  truth <- fit$params
  truth$beta[] <- c(2, -1.5)
  truth$common[] <- -0.5
  truth$thresholds[] <- c(-1, 0, 1, -0.5, 1, 2)
  truth$gamma[1L, 2L] <- log(0.4 / 0.6)
  expect_gt(fit$log_likelihood,
            latents:::.mixture_expectation(fit$spec, truth)$log_likelihood)
})

test_that("persons with repeated ordinal outcomes are classified (class_level = group)", {
  data <- ordinal_persons()
  fit <- mixture_regression(y ~ wave, data, n_classes = 2, family = "ordinal",
                            id = "person", class_level = "group", n_starts = 3, seed = 1)
  coefficients <- get_results(fit, "coefficients")
  expect_equal(subset(coefficients, term == "wave")$estimate, c(0.8, -0.6),
               tolerance = 0.2)
  thresholds <- subset(coefficients, startsWith(term, "threshold:"))
  expect_equal(thresholds$estimate, c(-1, 0.5, 2, -2, -0.5, 1), tolerance = 0.25)
  recovery <- get_results(fit, "recovery", data = data, truth = "class")
  expect_true(all(subset(recovery, n > 50)$share > 0.9))
  # The class mean of an ordinal outcome is its expected category score.
  fitted <- get_results(fit, "fitted")
  expect_true(all(fitted$fitted > 1 & fitted$fitted < 4))
  trajectory <- latents:::.growth_trajectory_table(
    fit, 0.95, latents:::.trajectory_inference(fit))
  expect_true(all(trajectory$conf_low <= trajectory$estimate &
                    trajectory$estimate <= trajectory$conf_high))
  expect_true(all(trajectory$estimate >= 1 & trajectory$estimate <= 4))
})

test_that("two categories are the logistic mixture with thresholds as negated intercepts", {
  hours <- study_hours
  ordinal <- mixture_regression(factor(passed) ~ hours, hours, n_classes = 2,
                                family = "ordinal", id = "student",
                                class_level = "group", n_starts = 2, seed = 1)
  binomial <- mixture_regression(passed ~ hours, hours, n_classes = 2,
                                 family = "binomial", id = "student",
                                 class_level = "group", n_starts = 2, seed = 1)
  expect_equal(ordinal$log_likelihood, binomial$log_likelihood, tolerance = 1e-10)
  expect_equal(unname(ordinal$params$thresholds[1L, ]),
               unname(-binomial$params$beta["(Intercept)", ]), tolerance = 1e-6)
  expect_equal(unname(ordinal$params$beta["hours", ]),
               unname(binomial$params$beta["hours", ]), tolerance = 1e-6)
})

test_that("category probabilities sum to one and relabelling leaves the likelihood", {
  skip_on_cran()
  data <- ordinal_rows(n = 600)
  fit <- mixture_regression(y ~ x, data, n_classes = 2, family = "ordinal",
                            n_starts = 2, seed = 1)
  probabilities <- predict(fit, type = "probabilities")
  totals <- stats::aggregate(probability ~ row + class, data = probabilities, FUN = sum)
  expect_identical(nrow(totals), 2L * nrow(data))
  expect_equal(totals$probability, rep(1, nrow(totals)), tolerance = 1e-12)
  expect_setequal(unique(probabilities$category), levels(data$y))
  # The expected score is the probability-weighted category number.
  response <- predict(fit, type = "class_response")
  scored <- stats::aggregate(probability * match(category, levels(data$y)) ~ row + class,
                             data = probabilities, FUN = sum)
  expect_equal(scored[[3L]], response$fitted, tolerance = 1e-12)
  swapped <- latents:::.mixture_permute(fit$spec, fit$params, c(2L, 1L))
  expect_equal(latents:::.mixture_expectation(fit$spec, swapped)$log_likelihood,
               fit$log_likelihood, tolerance = 1e-10)
  # New rows: posteriors given the outcome, in the fitted categories.
  new_rows <- data.frame(x = c(-1, 1), y = c("always", "never"))
  posterior <- predict(fit, new_rows, type = "posterior")
  expect_equal(posterior$probability_class_1 + posterior$probability_class_2,
               c(1, 1), tolerance = 1e-12)
})

test_that("unusable ordinal requests are refused with classed conditions", {
  data <- ordinal_rows(n = 300)
  fit_ordinal <- function(...) {
    mixture_regression(data = data, family = "ordinal", n_starts = 1, ...)
  }
  expect_error(fit_ordinal(formula = w ~ x, n_classes = 2), class = "latents_bad_data")
  expect_error(fit_ordinal(formula = as.character(y) ~ x, n_classes = 2),
               class = "latents_bad_data")
  expect_error(fit_ordinal(formula = rep(2, 300) ~ x, n_classes = 1),
               class = "latents_bad_data")
  expect_error(fit_ordinal(formula = y ~ 0 + x, n_classes = 2),
               class = "latents_bad_argument")
  expect_error(fit_ordinal(formula = y ~ x, n_classes = 2, variance = "equal"),
               class = "latents_bad_argument")
  expect_error(fit_ordinal(formula = I(y > "rarely") ~ x, n_classes = 2),
               class = "latents_bad_data")
  expect_error(fit_ordinal(formula = factor(y > "rarely") ~ x, n_classes = 2),
               class = "latents_not_identified")
  expect_error(fit_ordinal(formula = y ~ 1, n_classes = 2),
               class = "latents_not_identified")
  persons <- ordinal_persons(persons = 60)
  expect_error(mixture_regression(y ~ wave, persons, n_classes = 2, family = "ordinal",
                                  id = "person", class_level = "group",
                                  random = "wave", n_starts = 1),
               class = "latents_bad_argument")
  fit <- fit_ordinal(formula = y ~ x, n_classes = 1)
  expect_error(predict(fit, data.frame(x = 0, y = "sometimes"), type = "posterior"),
               class = "latents_bad_data")
  gaussian <- mixture_regression(w ~ x, data, n_classes = 1, n_starts = 1)
  expect_error(predict(gaussian, type = "probabilities"), class = "latents_bad_argument")
})

test_that("simulate() is seeded and draws the fitted categories", {
  skip_on_cran()
  data <- ordinal_rows(n = 600)
  fit <- mixture_regression(y ~ x, data, n_classes = 2, family = "ordinal",
                            n_starts = 2, seed = 1)
  draws <- simulate(fit, nsim = 2, seed = 4)
  expect_identical(simulate(fit, nsim = 2, seed = 4), draws)
  expect_false(identical(simulate(fit, nsim = 2, seed = 5), draws))
  expect_true(is.ordered(draws$sim_1))
  expect_identical(levels(draws$sim_2), levels(data$y))
  persons <- ordinal_persons(persons = 80)
  numbered <- mixture_regression(y ~ wave, persons, n_classes = 2, family = "ordinal",
                                 id = "person", class_level = "group",
                                 n_starts = 1, seed = 1)
  values <- unlist(simulate(numbered, nsim = 3, seed = 2))
  expect_true(is.numeric(values) && all(values %in% 1:4))
  expect_output(print(numbered), "Threshold 1\\|2")
})
