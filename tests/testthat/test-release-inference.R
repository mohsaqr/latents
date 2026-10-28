# Regression tests for the release review of 2026-10-27, cluster "inference":
# INF-1, INF-2, INF-3 (covariate fits and the multilpa dots), INF-4, INF-6
# and INF-7. Each test fails on the code before the fix.

.release_engagement_vars <- c("browse", "lectures", "forum_read", "forum_post",
                              "attendance")

# A third profile that lives in group 1 only (6 rows of it), the design of the
# review's probe p16.
.release_thin_data <- function() {
  set.seed(8)
  n_groups <- 30L
  per <- 8L
  g <- rep(seq_len(n_groups), each = per)
  z <- 1L + as.integer(stats::runif(n_groups * per) > 0.5)
  z[g == 1L][1:6] <- 3L
  mu <- rbind(c(-2, -2), c(2, 2), c(8, -8))
  data.frame(g = g,
             a = stats::rnorm(n_groups * per, mu[z, 1L], 0.7),
             b = stats::rnorm(n_groups * per, mu[z, 2L], 0.7),
             y = stats::rnorm(n_groups * per, c(0, 1, 5)[z]))
}

test_that("INF-1: a class carried by one group has no standard error in any method", {
  data <- .release_thin_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 3, n_group_classes = 1,
                  n_starts = 3, seed = 1)
  # The profile whose rows all sit in group 1.
  in_group_one <- tapply(data$g == 1L, fit$subject_profiles, all)
  thin <- as.integer(names(in_group_one)[in_group_one])
  expect_length(thin, 1L)

  lapply(c("modal", "proportional", "bch"), function(method) {
    expect_warning(means <- three_step(fit, data, "y", method = method),
                   class = "latents_small_classes")
    table <- as.data.frame(means)
    thin_row <- table$class == thin
    expect_true(is.na(table$standard_error[thin_row]))
    expect_true(is.na(table$conf_low[thin_row]) && is.na(table$conf_high[thin_row]))
    # The mean itself is still reported, and the other classes keep theirs.
    expect_true(is.finite(table$estimate[thin_row]))
    expect_true(all(table$standard_error[!thin_row] > 0.01))

    expect_warning(pairs <- three_step(fit, data, "y", method = method,
                                       contrast = "pairs"),
                   class = "latents_small_classes")
    pairs <- as.data.frame(pairs)
    involved <- pairs$class == thin | pairs$reference_class == thin
    expect_true(all(is.na(pairs$standard_error[involved])))
    expect_true(all(is.na(pairs$p_value[involved])))
    expect_true(all(is.na(pairs$p_adjusted[involved])))
    expect_true(all(is.finite(pairs$p_value[!involved])))
  })

  # Under `vcov_type = "independent"` the class has six carrying observations,
  # so its standard error is estimated and nothing is withheld.
  expect_no_warning(
    independent <- three_step(fit, data, "y", method = "bch", vcov_type = "independent"))
  expect_true(all(as.data.frame(independent)$standard_error > 0))
})

test_that("INF-1: ordinary data keep every standard error and raise no warning", {
  set.seed(11)
  g <- rep(seq_len(20), each = 8)
  truth <- 1L + as.integer(stats::runif(160) > 0.5)
  data <- data.frame(g = g, a = stats::rnorm(160, ifelse(truth == 2L, 1.5, -1.5)),
                     b = stats::rnorm(160, ifelse(truth == 2L, 1.5, -1.5)),
                     y = stats::rnorm(160, ifelse(truth == 2L, 2, 0)))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 2, seed = 1)
  expect_no_warning(means <- three_step(fit, data, "y"))
  expect_true(all(as.data.frame(means)$standard_error > 0))
})

test_that("INF-1: the carrier rule counts units by absolute weight", {
  pieces <- list(group_index = c(1L, 1L, 2L, 2L, 3L))
  # Class 1 has weight in group 1 only (a negative BCH weight elsewhere is
  # counted by its size); class 2 is spread over every group.
  carrying <- cbind(c(1, 1, 0, 0, 0), c(0.5, 0.5, 1, 1, 1))
  clustered <- matrix(1, nrow = 3L, ncol = 2L)
  expect_warning(kept <- .multilpa_step_thin_classes(clustered, carrying, pieces,
                                                     "cluster"),
                 class = "latents_small_classes")
  expect_true(all(is.na(kept[, 1L])))
  expect_identical(kept[, 2L], c(1, 1, 1))
  # A negative weight in another group is still a carrier.
  carrying[3L, 1L] <- -1
  expect_no_warning(.multilpa_step_thin_classes(clustered, carrying, pieces, "cluster"))
  # Per observation, the first class now has three carrying observations.
  expect_no_warning(.multilpa_step_thin_classes(matrix(1, 5L, 2L), carrying, pieces,
                                                "independent"))
})

test_that("INF-2: three_step() and r3step() name the families they refuse", {
  general <- lta(course_engagement, vars = .release_engagement_vars, id = "student",
                 n_profiles = 2, time = "sequence", transitions = "occasion",
                 n_starts = 2, seed = 1)
  expect_error(three_step(general, course_engagement, "previous_grade"),
               class = "latents_unsupported_three_step")
  expect_error(r3step(general, course_engagement, "previous_grade"),
               class = "latents_unsupported_three_step")
  refused <- tryCatch(three_step(general, course_engagement, "previous_grade"),
                      error = identity)
  expect_false(inherits(refused, "latents_not_a_fit"))

  covariate <- lta(course_engagement, vars = .release_engagement_vars, id = "student",
                   n_profiles = 2, time = "sequence",
                   transition_covariates = "previous_grade", n_starts = 2, seed = 1)
  expect_error(three_step(covariate, course_engagement, "previous_grade"),
               class = "latents_unsupported_three_step")

  set.seed(4)
  n_groups <- 60L
  g <- rep(seq_len(n_groups), each = 8L)
  h <- rep(sample(1:2, n_groups, TRUE), each = 8L)
  d <- data.frame(g = g,
                  a = stats::rnorm(n_groups)[g] + stats::rnorm(480, ifelse(h == 1L, 1, -1)),
                  b = stats::rnorm(n_groups)[g] + stats::rnorm(480, ifelse(h == 1L, 0.5, -0.5)),
                  y = stats::rnorm(n_groups)[g] + h, x = stats::rnorm(n_groups)[g])
  additive <- multilpa(d, c("a", "b"), "g", n_group_classes = 2, family = "additive",
                       n_starts = 2, seed = 1)
  expect_error(three_step(additive, d, "y", level = "groups"),
               class = "latents_unsupported_three_step")
  expect_error(r3step(additive, d, "x", level = "groups"),
               class = "latents_unsupported_three_step")

  # Not a fit at all is still `latents_not_a_fit`.
  expect_error(three_step(d, d, "y"), class = "latents_not_a_fit")
})

test_that("INF-2: parameter_inference() does not tell a mixture regression it is not a fit", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(120)
  regression <- mixture_regression(y ~ x, data = d, n_classes = 1, n_starts = 0)
  refused <- expect_error(parameter_inference(regression), class = "latents_bad_argument")
  expect_false(inherits(refused, "latents_not_a_fit"))
  expect_match(conditionMessage(refused), "coefficients", fixed = TRUE)
})

.release_covariate_data <- function() {
  set.seed(21)
  n <- 240L
  g <- rep(seq_len(40), each = 6L)
  x <- stats::rnorm(n)
  z <- 1L + as.integer(stats::runif(n) < stats::plogis(-0.3 + 0.8 * x))
  data.frame(g = g, x = x, a = stats::rnorm(n, ifelse(z == 1L, 1.5, -1.5)),
             b = stats::rnorm(n, ifelse(z == 1L, 1.2, -1), 1.2))
}

test_that("INF-3: vcov() of a covariate fit refuses or names what it would ignore", {
  d <- .release_covariate_data()
  fit <- multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  profile_covariates = "x", n_starts = 2, seed = 1)
  observed <- vcov(fit)
  expect_warning(misnamed <- vcov(fit, type = "robust"),
                 class = "latents_unused_argument")
  expect_identical(misnamed, observed)
  expect_error(vcov(fit, method = "bootstrap"), class = "latents_unsupported_inference")
  expect_no_warning(vcov(fit, method = "wald"))
  expect_false(isTRUE(all.equal(vcov(fit, vcov_type = "robust"), observed)))
  expect_error(confint(fit, type = "robust"), class = "latents_bad_argument")
})

test_that("INF-3: vcov() and confint() of a multilpa fit class a misnamed argument", {
  d <- .release_covariate_data()
  fit <- multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 2, seed = 1)
  expect_error(vcov(fit, type = "robust"), class = "latents_bad_argument")
  expect_error(confint(fit, type = "robust"), class = "latents_bad_argument")
  # The arguments parameter_inference() takes still pass, partial names included.
  expect_identical(vcov(fit, vcov_type = "robust"), vcov(fit, vcov = "robust"))
  expect_identical(nrow(confint(fit, vcov_type = "robust")), nrow(confint(fit)))
})

test_that("INF-4: sensitivity() of an lta() fit reproduces it under its own seed", {
  fit <- lta(course_engagement, vars = .release_engagement_vars, id = "student",
             n_profiles = 3, time = "sequence", n_starts = 2, seed = 1, tol = 1e-3)
  seeds <- as.data.frame(sensitivity(fit, data = course_engagement, seeds = 1:2))
  own <- subset(seeds, seed == 1)
  expect_equal(own$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  expect_equal(own$iterations, fit$iterations)
  expect_equal(own$agreement, 1)
})

test_that("INF-6: the inference verbs class invalid scalar arguments", {
  d <- .release_covariate_data()
  fit <- multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 2, seed = 1)
  bad <- "latents_bad_argument"
  expect_error(parameter_inference(fit, level = 1.5), class = bad)
  expect_error(parameter_inference(fit, step = -1), class = bad)
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 1), class = bad)
  expect_error(parameter_inference(fit, method = "bootstrap", n_starts = 0), class = bad)
  expect_error(parameter_inference(fit, data = as.matrix(d)), class = "latents_bad_data")
  expect_error(confint(fit, level = 2), class = bad)

  d$w <- stats::rnorm(nrow(d))
  expect_error(r3step(fit, d, "nope"), class = "latents_bad_data")
  expect_error(r3step(fit, d, "w", ci_level = 2), class = bad)
  expect_error(r3step(fit, as.matrix(d), "w"), class = "latents_bad_data")
  expect_error(r3step(fit, d, character()), class = bad)
  d$text <- "a"
  expect_error(r3step(fit, d, "text"), class = bad)

  covariate <- multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                        profile_covariates = "x", n_starts = 2, seed = 1)
  expect_error(parameter_inference(covariate, level = 2), class = bad)
  expect_error(pool_imputations(list(d, d), c("a", "b"), "g", 2, level = 2), class = bad)
})

test_that("INF-6: r3step() with transition covariates classes its argument checks", {
  moves <- lta(course_engagement, vars = .release_engagement_vars, id = "student",
               n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
  expect_error(r3step(moves, course_engagement, transition_covariates = "nope"),
               class = "latents_bad_data")
  expect_error(r3step(moves, course_engagement, transition_covariates = "sequence",
                      ci_level = 2), class = "latents_bad_argument")
  expect_error(r3step(moves, as.matrix(course_engagement),
                      transition_covariates = "sequence"), class = "latents_bad_data")
})

test_that("INF-7: coef() and confint() of pooled imputations agree with its table", {
  set.seed(1)
  completed <- lapply(seq_len(2), function(i) {
    within(subset(course_engagement, student <= 40),
           previous_grade <- previous_grade + stats::rnorm(length(previous_grade), sd = 0.1))
  })
  pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"), "student",
                             n_profiles = 2, n_group_classes = 1,
                             profile_covariates = "previous_grade", n_starts = 2, seed = 1)
  table <- as.data.frame(pooled)
  estimates <- coef(pooled)
  expect_type(estimates, "double")
  expect_identical(names(estimates), rownames(vcov(pooled)))
  expect_equal(unname(estimates), table$estimate)

  intervals <- confint(pooled)
  expect_identical(dim(intervals), c(nrow(table), 2L))
  expect_identical(colnames(intervals), c("2.5%", "97.5%"))
  expect_equal(unname(intervals[, 1L]), table$conf_low)
  expect_equal(unname(intervals[, 2L]), table$conf_high)
  narrow <- confint(pooled, level = 0.8)
  expect_true(all(narrow[, 2L] - narrow[, 1L] <= intervals[, 2L] - intervals[, 1L]))

  first <- names(estimates)[[1L]]
  expect_identical(confint(pooled, parm = 1), confint(pooled, parm = first))
  expect_error(confint(pooled, parm = "nope"), class = "latents_bad_argument")
  expect_error(confint(pooled, parm = 0), class = "latents_bad_argument")
  expect_error(confint(pooled, level = 2), class = "latents_bad_argument")
})

test_that("INF-6: an invalid covariance type or adjustment is a classed refusal", {
  d <- .release_covariate_data()
  fit <- multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 2, seed = 1)
  covariate <- multilpa(d, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                        profile_covariates = "x", n_starts = 2, seed = 1)
  expect_error(parameter_inference(fit, vcov_type = "sandwich"),
               class = "latents_bad_argument")
  expect_error(parameter_inference(fit, adjust = "nope"), class = "latents_bad_argument")
  expect_error(parameter_inference(fit, boundary = "nope"), class = "latents_bad_argument")
  expect_error(vcov(fit, scale = "nope"), class = "latents_bad_argument")
  expect_error(parameter_inference(covariate, vcov_type = "sandwich"),
               class = "latents_bad_argument")
  # A misspelt method is a bad argument, not an unimplemented bootstrap.
  expect_error(parameter_inference(covariate, method = "nope"),
               class = "latents_bad_argument")
  expect_error(parameter_inference(covariate, method = "boot"),
               class = "latents_unsupported_inference")
  # Valid abbreviations still match, as match.arg() matched them.
  expect_identical(parameter_inference(fit, vcov_type = "rob"),
                   parameter_inference(fit, vcov_type = "robust"))
  regression <- mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                   n_starts = 1, seed = 1)
  expect_error(get_results(regression, "coefficients", vcov_type = "sandwich"),
               class = "latents_bad_argument")
  expect_error(get_results(regression, "recovery", by = "nope"),
               class = "latents_bad_argument")
})

test_that("INF-2: a verb names the family it does not take, not 'not a fit'", {
  regression <- mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                   n_starts = 1, seed = 1)
  condition <- expect_error(get_tna(regression), class = "latents_not_a_fit")
  expect_false(grepl("must be a fitted model of this package", conditionMessage(condition),
                     fixed = TRUE))
  # Anything that is not a fit at all is still told so, for every family.
  expect_error(parameter_inference(data.frame(a = 1)), class = "latents_not_a_fit")
})
