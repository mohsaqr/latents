# Pre-release review 2026-10-27, cluster "regression": mixture regressions,
# growth mixtures and additive fits. Each block fails on the code before the
# fix (REG-1..8, API-1/2/4/5/6, INF-3).

release_three_class_data <- function(seed = 1L, n = 120L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  class <- sample(1:3, n, TRUE)
  data.frame(y = c(0, 1, 2)[class] + c(1, -1, 0.3)[class] * x + stats::rnorm(n, sd = 0.9),
             x = x, w = 1)
}

release_quiet <- function(expr) {
  withCallingHandlers(expr, latents_degenerate_start = function(w) invokeRestart("muffleWarning"),
                      latents_unconverged = function(w) invokeRestart("muffleWarning"))
}

release_same_fit <- function(plain, padded) {
  expect_equal(padded$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
  expect_equal(padded$starts, plain$starts, tolerance = 1e-10)
  expect_equal(coef(padded), coef(plain), tolerance = 1e-8)
  expect_identical(padded$n_parameters, plain$n_parameters)
}

test_that("REG-1: a zero-weight row changes no start, likelihood or estimate (observation)", {
  data <- release_three_class_data()
  fit <- function(d) release_quiet(mixture_regression(
    y ~ x, d, n_classes = 3, weights = "w", n_starts = 3, seed = 11, vcov_type = "none"))
  plain <- fit(data)
  padded <- fit(rbind(data, data.frame(y = c(1e6, -1e6), x = c(30, 0), w = 0)))
  release_same_fit(plain, padded)
  # The absent rows keep their row and posterior.
  expect_identical(nrow(as.data.frame(get_results(padded, "assignments"))), nrow(data) + 2L)
  expect_identical(nobs(padded), nobs(plain))
})

test_that("REG-1: a zero-weight group changes nothing under group and two-level nesting", {
  set.seed(5)
  groups <- 30L
  size <- 6L
  id <- rep(seq_len(groups), each = size)
  x <- stats::rnorm(groups * size)
  class <- rep(stats::rbinom(groups, 1, 0.5), each = size)
  data <- data.frame(id = id, x = x, w = 1,
                     y = ifelse(class == 1, 2 + 1.5 * x, -1 - x) +
                       stats::rnorm(groups * size, sd = 0.5))
  pad <- data.frame(id = 99L, x = c(25, -25, 0), w = 0, y = c(1e6, -1e6, 5e5))
  group_fit <- function(d) release_quiet(mixture_regression(
    y ~ x, d, n_classes = 2, id = "id", class_level = "group", weights = "w",
    n_starts = 2, seed = 3, vcov_type = "none"))
  release_same_fit(group_fit(data), group_fit(rbind(data, pad)))
  two_level <- function(d) release_quiet(mixture_regression(
    y ~ x, d, n_classes = 2, id = "id", n_group_classes = 2, weights = "w",
    n_starts = 2, seed = 3, vcov_type = "none"))
  release_same_fit(two_level(data), two_level(rbind(data, pad)))
})

test_that("REG-1: fits without zero weights are unchanged by the carrying restriction", {
  data <- release_three_class_data(seed = 2L)
  fit <- release_quiet(mixture_regression(y ~ x, data, n_classes = 2, weights = "w",
                                          n_starts = 2, seed = 4, vcov_type = "none"))
  spec <- fit$spec
  expect_identical(latents:::.mixture_carrying_spec(spec), spec)
})

test_that("REG-2: a factor level or ordinal category on zero-weight rows only adds no parameter", {
  set.seed(7)
  n <- 300L
  x <- stats::rnorm(n)
  f <- factor(sample(c("a", "b"), n, TRUE), levels = c("a", "b", "c"))
  eta <- 0.8 * x + 0.5 * (f == "b")
  data <- data.frame(o = cut(eta + stats::rlogis(n), c(-Inf, -1, 0.5, 1.5, Inf), labels = FALSE),
                     x = x, f = f, w = stats::runif(n, 0.5, 1.5),
                     y = 1 + x + 0.5 * (f == "b") + stats::rnorm(n))
  level_pad <- data.frame(o = 2L, x = c(0.1, -0.2), f = factor("c", levels = c("a", "b", "c")),
                          w = 0, y = 3)
  plain <- mixture_regression(y ~ x + f, data, n_classes = 1, weights = "w", n_starts = 0)
  # The default robust covariance used to refuse the padded fit as rank deficient.
  expect_warning(padded <- mixture_regression(y ~ x + f, rbind(data, level_pad), n_classes = 1,
                                              weights = "w", n_starts = 0),
                 class = "latents_rows_dropped")
  expect_identical(padded$n_parameters, plain$n_parameters)
  expect_equal(get_results(padded, "fit")$bic, get_results(plain, "fit")$bic)
  expect_equal(as.data.frame(padded)$std_error, as.data.frame(plain)$std_error,
               tolerance = 1e-8)

  category_pad <- data.frame(o = 5L, x = c(0.1, -0.2), f = factor("a", levels = c("a", "b", "c")),
                             w = 0, y = 3)
  ordinal <- mixture_regression(o ~ x + f, data, n_classes = 1, family = "ordinal",
                                weights = "w", n_starts = 0, vcov_type = "none")
  expect_warning(ordinal_padded <- mixture_regression(
    o ~ x + f, rbind(data, category_pad), n_classes = 1, family = "ordinal",
    weights = "w", n_starts = 0, vcov_type = "none"), class = "latents_rows_dropped")
  expect_identical(ordinal_padded$n_parameters, ordinal$n_parameters)
  expect_equal(ordinal_padded$log_likelihood, ordinal$log_likelihood, tolerance = 1e-10)
})

test_that("REG-2: zero weights that leave a factor one level, or no unit, are refused with a class", {
  set.seed(8)
  n <- 80L
  data <- data.frame(y = stats::rnorm(n), x = stats::rnorm(n),
                     f = factor(rep(c("a", "b"), each = n / 2L)))
  # No unit of positive weight: the weights refusal, not base R's contrasts error.
  expect_error(mixture_regression(y ~ x + f, transform(data, w = 0), n_classes = 1,
                                  weights = "w", n_starts = 0),
               class = "latents_bad_weights")
  # A subgroup selected by zero weights whose indicator is in the formula:
  # refused with latents_bad_data, and no rows-dropped warning beforehand.
  subgroup <- transform(data, w = as.numeric(f == "a"))
  dropped <- 0L
  condition <- withCallingHandlers(
    tryCatch(mixture_regression(y ~ x + f, subgroup, n_classes = 1, weights = "w",
                                n_starts = 0, vcov_type = "none"), error = identity),
    latents_rows_dropped = function(w) {
      dropped <<- dropped + 1L
      invokeRestart("muffleWarning")
    })
  expect_s3_class(condition, "latents_bad_data")
  expect_null(conditionCall(condition))
  expect_identical(dropped, 0L)
  # The unweighted fit with a one-level factor now gets the same class.
  expect_error(mixture_regression(y ~ x + f, subset(data, f == "a"), n_classes = 1,
                                  n_starts = 0),
               class = "latents_bad_data")
  expect_error(mixture_regression(y ~ x, subset(data, f == "a"), n_classes = 2,
                                  membership = ~ f, n_starts = 1, seed = 1),
               class = "latents_bad_data")
})

test_that("REG-3: a residual SD held on its floor has no standard error, not zero", {
  set.seed(11)
  persons <- 120L
  waves <- 6L
  id <- rep(seq_len(persons), each = waves)
  time <- rep(0:(waves - 1L), persons)
  class <- rep(stats::rbinom(persons, 1, 0.5), each = waves)
  intercept <- rep(stats::rnorm(persons), each = waves)
  # An autoregressive residual process and no white noise: the residual SD of
  # one class sits on its floor.
  process <- unlist(lapply(seq_len(persons), function(i) {
    Reduce(function(previous, shock) exp(-0.5) * previous + shock,
           stats::rnorm(waves - 1L, sd = 1.5 * sqrt(1 - exp(-1))),
           accumulate = TRUE, stats::rnorm(1, sd = 1.5))
  }))
  data <- data.frame(id = id, t = time,
                     y = ifelse(class == 1, 10 + time, 14 - 0.8 * time) + intercept + process)
  fit <- suppressWarnings(mixture_regression(
    y ~ t, data, n_classes = 2, id = "id", class_level = "group", random = "intercept",
    residual_correlation = "ar", residual_time = "t", n_starts = 1, seed = 1))
  held <- fit$residual_boundary
  skip_if(length(held) == 0L || !any(startsWith(held, "log_sigma")),
          "the simulated fit did not hold a residual SD on its floor")
  classes <- get_results(fit, "classes")
  floored <- paste0("class_", sub("^log_sigma\\.class_", "", held[startsWith(held, "log_sigma")]))
  expect_true(all(is.na(classes$residual_sd_std_error[classes$class %in% floored])))
  expect_true(all(is.finite(classes$residual_sd_std_error[!classes$class %in% floored])))
  expect_true(all(is.na(diag(vcov(fit))[held])))
  intervals <- confint(fit)
  expect_true(all(is.na(intervals$conf_low[intervals$parameter %in% held])))
  # Every other coordinate keeps its interval.
  free <- !intervals$parameter %in% held
  expect_true(all(is.finite(intervals$conf_low[free])))
  # The refusal names the residual SD, not a random-effect covariance.
  refusal <- tryCatch(parameter_inference(fit), latents_boundary_fit = identity)
  expect_s3_class(refusal, "latents_boundary_fit")
  expect_match(conditionMessage(refusal), "residual standard deviation")
})

test_that("REG-5: zero-weight groups are not counted among the assigned groups", {
  set.seed(5)
  groups <- 40L
  size <- 8L
  g <- rep(seq_len(groups), each = size)
  group_class <- rep(stats::rbinom(groups, 1, 0.5), each = size)
  class <- stats::rbinom(groups * size, 1, ifelse(group_class == 1, 0.85, 0.2))
  x <- stats::rnorm(groups * size)
  data <- data.frame(y = ifelse(class == 1, 2 + 1.5 * x, -1 - x) +
                       stats::rnorm(groups * size, sd = 0.5),
                     x = x, g = g, w = ifelse(g %in% 1:6, 0, 1))
  fit <- release_quiet(mixture_regression(y ~ x, data, n_classes = 2, id = "g",
                                          n_group_classes = 2, weights = "w",
                                          n_starts = 2, seed = 2, vcov_type = "none"))
  table <- get_results(fit, "group_classes")
  first <- table[table$class == "class_1", "n_groups_assigned"]
  expect_identical(sum(first), 34L)
  expect_identical(sum(first), nobs(fit))
})

test_that("REG-6: confint() refuses a level outside (0, 1)", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1, seed = 1)
  expect_error(confint(fit, level = -0.5), class = "latents_bad_argument")
  expect_error(confint(fit, level = 2), class = "latents_bad_argument")
  expect_error(confint(fit, level = c(0.9, 0.95)), class = "latents_bad_argument")
  growth <- mixture_regression(score ~ wave, growth_scores, n_classes = 2, id = "student",
                               class_level = "group", random = "intercept",
                               n_starts = 1, seed = 1)
  expect_error(confint(growth, level = -0.5), class = "latents_bad_argument")
  expect_error(confint(growth, level = 2), class = "latents_bad_argument")
})

test_that("REG-7: the regression surface classes its refusals", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1, seed = 1)
  growth <- mixture_regression(score ~ wave, growth_scores, n_classes = 2, id = "student",
                               class_level = "group", random = "intercept",
                               n_starts = 1, seed = 1)
  expect_error(simulate(fit, nsim = 0), class = "latents_bad_argument")
  expect_error(simulate(growth, nsim = 0), class = "latents_bad_argument")
  expect_error(get_results(fit, level = 2), class = "latents_bad_argument")
  expect_error(get_results(growth, level = 2), class = "latents_bad_argument")
  expect_error(enumerate_regressions(score ~ hours, study_hours, n_classes = 0),
               class = "latents_bad_argument")
  expect_error(enumerate_regressions(score ~ hours, study_hours, n_classes = 1, bootstrap = -1),
               class = "latents_bad_argument")
  expect_error(get_results(fit, "recovery", data = study_hours, truth = "nope"),
               class = "latents_bad_data")
})

test_that("REG-8: a `common` with no predictor term is refused", {
  expect_error(mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                  common = "intercept", n_starts = 0),
               class = "latents_bad_argument")
  expect_error(mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                  common = ~ 1, n_starts = 0),
               class = "latents_bad_argument")
})

test_that("API-1: enumerate_regressions() names an unknown argument, with no call", {
  error <- tryCatch(enumerate_regressions(score ~ hours, study_hours, n_classes = 1:2, foo = 1),
                    error = identity)
  expect_s3_class(error, "latents_bad_argument")
  expect_null(conditionCall(error))
  # The backstop: an error raised inside a cell never carries the data, even
  # base R's own "unused argument" error of the do.call().
  bounded <- tryCatch(latents:::.mixture_try_fit(score ~ hours, study_hours, 2L, 1L,
                                                 list(foo = 1)), error = identity)
  expect_s3_class(bounded, "error")
  expect_lt(nchar(paste(deparse(conditionCall(bounded)), collapse = "")), 1000L)
})

test_that("API-2: a data error shared by every cell is raised; partial failures are printed", {
  expect_error(enumerate_regressions(score ~ hours, study_hours, n_classes = 1:2,
                                     n_starts = 1, common = ~ nope),
               class = "latents_bad_data")
  small <- utils::head(study_hours, 4L)
  grid <- release_quiet(enumerate_regressions(score ~ hours, small, n_classes = c(1, 5),
                                              n_starts = 0, seed = 1))
  expect_identical(sum(!is.na(as.data.frame(grid)$note)), 1L)
  expect_output(print(grid), "1 model\\(s\\) could not be fitted")
})

test_that("API-4: regression, growth and additive verbs refuse arguments they do not use", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1, seed = 1)
  growth <- mixture_regression(score ~ wave, growth_scores, n_classes = 2, id = "student",
                               class_level = "group", random = "intercept",
                               n_starts = 1, seed = 1)
  expect_error(get_results(fit, "coefficients", levl = 0.5), class = "latents_bad_argument")
  expect_error(as.data.frame(fit, levl = 0.5), class = "latents_bad_argument")
  expect_error(confint(fit, levl = 0.5), class = "latents_bad_argument")
  expect_error(summary(fit, levl = 0.5), class = "latents_bad_argument")
  expect_error(get_results(growth, "coefficients", levl = 0.5), class = "latents_bad_argument")
  expect_error(as.data.frame(growth, levl = 0.5), class = "latents_bad_argument")
  expect_error(confint(growth, levl = 0.5), class = "latents_bad_argument")
  expect_error(summary(growth, levl = 0.5), class = "latents_bad_argument")
  # The arguments the tables do take still work.
  expect_s3_class(get_results(fit, "coefficients", level = 0.9), "data.frame")
  expect_s3_class(get_results(growth, "classes", level = 0.9), "data.frame")
})

test_that("INF-3: vcov() warns about a covariance argument under the wrong name", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1, seed = 1)
  growth <- mixture_regression(score ~ wave, growth_scores, n_classes = 2, id = "student",
                               class_level = "group", random = "intercept",
                               n_starts = 1, seed = 1)
  expect_warning(vcov(fit, vcov_type = "robust"), class = "latents_unused_argument")
  expect_warning(vcov(growth, vcov_type = "robust"), class = "latents_unused_argument")
  expect_no_warning(vcov(fit, type = "robust"))
})

test_that("API-5: a `what` the object does not offer is a classed refusal", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1, seed = 1)
  growth <- mixture_regression(score ~ wave, growth_scores, n_classes = 2, id = "student",
                               class_level = "group", random = "intercept",
                               n_starts = 1, seed = 1)
  expect_error(get_results(fit, "nope"), class = "latents_bad_argument")
  expect_error(get_results(growth, "classification"), class = "latents_bad_argument")
  grid <- enumerate_regressions(score ~ hours, study_hours, n_classes = 1:2,
                                n_starts = 1, seed = 1)
  expect_error(get_results(grid, "nope"), class = "latents_bad_argument")
  everything <- get_results(grid, "all")
  expect_identical(names(everything), c("fit", "tree"))
  expect_identical(everything$fit, as.data.frame(grid))
  expect_error(get_results(grid, "fit", foo = 1), class = "latents_bad_argument")
  expect_error(get_results(grid, "model", n_classes = "2"), class = "latents_bad_argument")
})

test_that("API-6: the summaries of regression and growth fits give their tables", {
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 2, n_starts = 1, seed = 1)
  growth <- mixture_regression(score ~ wave, growth_scores, n_classes = 2, id = "student",
                               class_level = "group", random = "intercept",
                               n_starts = 1, seed = 1)
  fit_summary <- summary(fit, level = 0.9)
  expect_identical(get_results(fit_summary, "classes"),
                   get_results(fit, "classes", level = 0.9))
  expect_identical(as.data.frame(fit_summary), get_results(fit, "coefficients", level = 0.9))
  expect_identical(names(get_results(fit_summary, "all")), names(fit_summary))
  growth_summary <- summary(growth)
  expect_identical(get_results(growth_summary, "random"), get_results(growth, "random"))
  expect_identical(as.data.frame(growth_summary), get_results(growth, "coefficients"))
  expect_error(get_results(growth_summary, "nope"), class = "latents_bad_argument")
  expect_error(get_results(fit_summary, "classes", level = 0.5), class = "latents_bad_argument")
})

test_that("API-4/5/6, INF-3: the additive fit refuses unused arguments and serves its summary", {
  set.seed(1)
  ratings <- data.frame(
    team = rep(seq_len(40L), each = 6L),
    climate = rep(stats::rnorm(40L, rep(c(-1, 1), each = 20L), 0.5), each = 6L) +
      stats::rnorm(240L),
    support = rep(stats::rnorm(40L, rep(c(1, -1), each = 20L), 0.7), each = 6L) +
      stats::rnorm(240L))
  fit <- multilpa(ratings, c("climate", "support"), "team", n_group_classes = 2,
                  family = "additive", n_starts = 1, seed = 1)
  # Forty teams make weak classes; that warning is the fixture's, not the test's.
  weak <- c(.multilpa_expected_warnings, "latents_weak_class")
  quietly(expect_warning(vcov(fit, vcov_type = "robust"), class = "latents_unused_argument"),
          classes = weak)
  expect_error(confint(fit, levl = 0.5), class = "latents_bad_argument")
  expect_error(confint(fit, level = 2), class = "latents_bad_argument")
  expect_error(get_results(fit, "parameters", levl = 0.5), class = "latents_bad_argument")
  expect_error(get_results(fit, "nope"), class = "latents_bad_argument")
  expect_error(get_results(fit, level = 2), class = "latents_bad_argument")
  expect_error(summary(fit, levl = 0.5), class = "latents_bad_argument")
  fit_summary <- quietly(summary(fit), classes = weak)
  expect_identical(get_results(fit_summary, "group_classes"), fit_summary$group_classes)
  expect_identical(as.data.frame(fit_summary), get_results(fit_summary))
})

test_that("INF-5: confint() of an additive fit refuses a parm that is not an estimate", {
  set.seed(1)
  ratings <- data.frame(
    team = rep(seq_len(40L), each = 6L),
    climate = rep(stats::rnorm(40L, rep(c(-1, 1), each = 20L), 0.5), each = 6L) +
      stats::rnorm(240L),
    support = rep(stats::rnorm(40L, rep(c(1, -1), each = 20L), 0.7), each = 6L) +
      stats::rnorm(240L))
  fit <- multilpa(ratings, c("climate", "support"), "team", n_group_classes = 2,
                  family = "additive", n_starts = 1, seed = 1)
  weak <- c(.multilpa_expected_warnings, "latents_weak_class")
  bounds <- quietly(confint(fit), classes = weak)
  expect_identical(quietly(confint(fit, parm = 2), classes = weak),
                   bounds[2L, , drop = FALSE])
  expect_identical(quietly(confint(fit, parm = rownames(bounds)[3L]), classes = weak),
                   bounds[3L, , drop = FALSE])
  expect_error(quietly(confint(fit, parm = 999), classes = weak),
               class = "latents_bad_argument")
  expect_error(quietly(confint(fit, parm = "nope"), classes = weak),
               class = "latents_bad_argument")
  expect_error(quietly(confint(fit, type = "sandwich"), classes = weak),
               class = "latents_bad_argument")
  expect_error(vcov(fit, scale = "nope"), class = "latents_bad_argument")
})
