# Counterexamples from the November review: designs, identification and output.
skip_on_cran()

fix1105_regression_data <- function(n = 60L, seed = 110511L) {
  withr::local_seed(seed)
  id <- rep(seq_len(n), each = 4L)
  time <- rep(0:3, n)
  state <- rep(sample(0:1, n, TRUE), each = 4L)
  data.frame(id = id, time = time,
    y = 3 * state + (.7 - 1.1 * state) * time +
      rep(rnorm(n, sd = .7), each = 4L) + rnorm(length(id), sd = .5),
    member = factor(rep(ifelse(seq_len(n) %% 2L == 0L, "a", "b"), each = 4L),
                      levels = c("a", "b", "unused")),
    truth = state)
}

test_that("unused membership and group-membership levels do not become parameters", {
  data <- fix1105_regression_data()
  fit <- suppressWarnings(mixture_regression(y ~ time, data, 2, id = "id",
    class_level = "group", membership = ~ member, seed = 91, n_starts = 2))
  clean <- droplevels(data)
  reference <- mixture_regression(y ~ time, clean, 2, id = "id",
    class_level = "group", membership = ~ member, seed = 91, n_starts = 2)
  expect_identical(colnames(fit$spec$w), c("(Intercept)", "memberb"))
  expect_identical(fit$n_parameters, reference$n_parameters)
  expect_equal(fit$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_equal(get_results(fit, "fit")$bic, get_results(reference, "fit")$bic)
  expect_true(all(is.finite(vcov(fit))))
  expect_equal(predict(fit, data, type = "response")$fitted,
                predict(reference, clean, type = "response")$fitted)
  grouped <- function(d) suppressWarnings(mixture_regression(y ~ time, d, 2,
    id = "id", n_group_classes = 2, group_membership = ~ member,
    seed = 91, n_starts = 1, vcov_type = "none"))
  with_unused <- grouped(data)
  without_unused <- grouped(clean)
  expect_identical(colnames(with_unused$spec$v), c("(Intercept)", "memberb"))
  expect_identical(with_unused$n_parameters, without_unused$n_parameters)
  expect_equal(with_unused$log_likelihood, without_unused$log_likelihood, tolerance = 1e-10)
})

test_that("unused random factor levels are dropped with their prediction coding", {
  data <- fix1105_regression_data(n = 80)
  fit <- suppressWarnings(mixture_regression(y ~ time, data, 1, id = "id",
    class_level = "group", random = ~ 0 + member, random_diagonal = TRUE,
    seed = 71, n_starts = 1))
  expect_identical(colnames(fit$spec$random_design), c("membera", "memberb"))
  expect_identical(fit$spec$random_xlevels$member, c("a", "b"))
  expect_identical(fit$n_parameters, 5L)
  expect_true(all(is.finite(predict(fit, data, type = "response")$fitted)))
  expect_error(predict(fit, transform(data[1:4, ], member = "unused")),
                class = "latents_bad_data")
})

test_that("a membership level carried only by zero-weight units has no parameter", {
  data <- droplevels(fix1105_regression_data())
  zero <- data[data$id == 1L, ]
  zero$id <- max(data$id) + 1L
  weighted <- rbind(data, zero)
  weighted$member <- factor(c(as.character(data$member), rep("rare", nrow(zero))))
  weighted$w <- c(rep(1, nrow(data)), rep(0, nrow(zero)))
  expect_warning(fit <- mixture_regression(y ~ time, weighted, 2, id = "id",
    class_level = "group", membership = ~ member, weights = "w",
    seed = 91, n_starts = 2, vcov_type = "none"), class = "latents_rows_dropped")
  reference <- mixture_regression(y ~ time, data, 2, id = "id",
    class_level = "group", membership = ~ member, seed = 91,
    n_starts = 2, vcov_type = "none")
  expect_identical(fit$n_parameters, reference$n_parameters)
  expect_identical(colnames(fit$spec$w), colnames(reference$spec$w))
  expect_equal(fit$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_equal(get_results(fit, "fit")$bic, get_results(reference, "fit")$bic)
})

test_that("growth covariance identification respects the actual designs and restrictions", {
  data <- fix1105_regression_data()
  two <- subset(data, time < 2)
  fit <- function(d, ...) suppressWarnings(mixture_regression(y ~ time, d, 1,
    id = "id", class_level = "group", random = "time", seed = 62,
    n_starts = 1, vcov_type = "none", ...))
  expect_error(fit(two), class = "latents_not_identified")
  expect_error(fit(transform(two, time = time * 1e-6 + 2000)),
                class = "latents_not_identified")
  diagonal <- fit(two, random_diagonal = TRUE)
  expect_s3_class(diagonal, "latents_growth_mixture")
  expect_identical(diagonal$n_parameters, 5L)
  varying <- two
  varying$time <- varying$time + .4 * (varying$id %% 2L)
  expect_s3_class(fit(varying), "latents_growth_mixture")
  expect_s3_class(fit(subset(data, time < 3)), "latents_growth_mixture")
  proportional <- suppressWarnings(mixture_regression(y ~ time, two, 2,
    id = "id", class_level = "group", random = "time",
    random_covariance = "proportional", variance = "equal", seed = 62,
    n_starts = 1, vcov_type = "none"))
  expect_s3_class(proportional, "latents_growth_mixture")
  zero <- subset(data, time < 2 | id == 1L)
  zero$w <- as.numeric(zero$id != 1L)
  expect_error(fit(zero, weights = "w"), class = "latents_not_identified")
})

test_that("binary and ordinal repeated mixtures need enough response-count patterns", {
  data <- fix1105_regression_data()
  withr::local_seed(110512)
  data$b <- rbinom(nrow(data), 1, ifelse(data$truth == 0, .2, .8))
  data$o <- ifelse(data$b == 0, 0L, 2L)
  fit <- function(family, formula, d) suppressWarnings(mixture_regression(formula,
    d, 2, family = family, id = "id", class_level = "group", seed = 5,
    n_starts = 1, vcov_type = "none"))
  expect_error(fit("binomial", b ~ 1, subset(data, time < 2)),
                class = "latents_not_identified")
  expect_error(fit("ordinal", o ~ 1, subset(data, time < 2)),
                class = "latents_not_identified")
  expect_s3_class(fit("binomial", b ~ 1, subset(data, time < 3)),
                   "latents_mixture_regression")
})

test_that("class-named IDs survive regression and growth output and predictions", {
  data <- droplevels(fix1105_regression_data())
  data$class <- data$id + 1000L
  regression <- mixture_regression(y ~ time, data, 2, id = "class",
    class_level = "group", seed = 32, n_starts = 2, vcov_type = "none")
  rows <- get_results(regression, "assignments")
  groups <- get_results(regression, "groups")
  expect_identical(rows$id_class, data$class)
  expect_identical(groups$id_class, sort(unique(data$class)))
  expect_true(all(rows$class %in% 1:2))
  expect_false(anyDuplicated(names(rows)) > 0L)
  growth <- suppressWarnings(mixture_regression(y ~ time, data, 1, id = "class",
    class_level = "group", random = "intercept", seed = 32, n_starts = 1))
  assignments <- get_results(growth, "assignments")
  individuals <- get_results(growth, "individual")
  effects <- get_results(growth, "random_effects")
  expect_identical(assignments$id_class, sort(unique(data$class)))
  expect_identical(individuals$id_class, data$class)
  expect_true(all(individuals$class == "class_1"))
  expect_identical(effects$id_class, sort(unique(data$class)))
  expect_identical(predict(growth, type = "posterior")$id_class,
                    assignments$id_class)
  expect_identical(predict(growth, type = "random_effects")$id_class,
                    assignments$id_class)
  interval <- get_results(growth, "classes")
  expect_identical(interval$share_conf_low, 1)
  expect_identical(interval$share_conf_high, 1)
})

test_that("cluster identifiers and occupied identifier prefixes remain distinct", {
  data <- droplevels(fix1105_regression_data())
  data$class <- 500L + ceiling(data$id / 5L)
  fit <- suppressWarnings(mixture_regression(y ~ time, data, 1, id = "id",
    cluster = "class", class_level = "group", random = "intercept",
    seed = 41, n_starts = 1, vcov_type = "none"))
  tab <- get_results(fit, "assignments")
  expect_identical(tab$id, sort(unique(data$id)))
  expect_identical(tab$cluster_class, data$class[match(tab$id, data$id)])
  expect_true(all(tab$class == "class_1"))
  expect_false(anyDuplicated(names(tab)) > 0L)
  expect_identical(get_results(fit, "clusters")$cluster_class, sort(unique(data$class)))
  data$id_class <- data$time
  collision <- mixture_regression(y ~ time + id_class, data, 1, id = "class",
    class_level = "group", vcov_type = "none", n_starts = 1, seed = 4)
  assignments <- get_results(collision, "assignments")
  expect_identical(assignments$id_class.1, data$class)
})

test_that("individual trajectory plots select the original collision-safe IDs", {
  skip_if_not_installed("ggplot2")
  data <- droplevels(fix1105_regression_data(n = 35))
  data$class <- data$id + 1000L
  fit <- suppressWarnings(mixture_regression(y ~ time, data, 1, id = "class",
    class_level = "group", random = "intercept", seed = 32, n_starts = 1,
    vcov_type = "none"))
  chosen <- c(1001L, 1005L, 1012L)
  individual <- plot(fit, what = "individuals", persons = chosen)
  expect_setequal(individual$data$person, chosen)
  expect_identical(individual$data$person, individual$data$id_class)
  expect_false(anyNA(individual$data$person_label))
  expect_s3_class(ggplot2::ggplot_build(individual), "ggplot_built")
})

test_that("recovery refuses another row order and two-level labels count rows", {
  data <- droplevels(fix1105_regression_data())
  fit <- mixture_regression(y ~ time, data, 2, seed = 12, n_starts = 2,
    vcov_type = "none")
  correct <- get_results(fit, "recovery", data = data, truth = "truth")
  expect_identical(sum(correct$n), nrow(data))
  expect_error(get_results(fit, "recovery", data = data[rev(seq_len(nrow(data))), ],
    truth = "truth"), class = "latents_bad_inference_data")
  grouped <- suppressWarnings(mixture_regression(y ~ time, data, 2, id = "id",
    n_group_classes = 2, seed = 12, n_starts = 1, vcov_type = "none"))
  expect_output(print(get_results(grouped, "classes")), "Rows assigned")
  expect_output(print(get_results(grouped, "recovery", data = data, truth = "truth")),
                "Rows")
})

test_that("exponentiated intercepts are not labeled as slope ratios", {
  data <- droplevels(fix1105_regression_data())
  withr::local_seed(110513)
  data$count <- rpois(nrow(data), exp(.4 + .2 * data$time))
  fit <- mixture_regression(count ~ time, data, 1, family = "poisson", seed = 7)
  coefficients <- get_results(fit, "coefficients")
  expect_true(all(is.na(coefficients$exp_estimate[coefficients$term == "(Intercept)"])))
  expect_equal(coefficients$exp_estimate[coefficients$term == "time"],
                exp(coefficients$estimate[coefficients$term == "time"]))
  membership <- mixture_regression(y ~ time, data, 2, id = "id", class_level = "group",
    membership = ~ member, seed = 9, n_starts = 2)
  table <- get_results(membership, "membership")
  expect_true(all(is.na(table$odds_ratio[table$term == "(Intercept)"])))
  expect_equal(table$odds_ratio[table$term != "(Intercept)"],
                exp(table$estimate[table$term != "(Intercept)"]))
})

test_that("correlated growth starts describe the fitted model rather than its initializer", {
  data <- droplevels(fix1105_regression_data(n = 35L, seed = 110514L))
  fit <- suppressWarnings(mixture_regression(y ~ time, data, 1, id = "id",
    class_level = "group", random = "intercept", residual_correlation = "ar",
    seed = 52, n_starts = 1, vcov_type = "none"))
  starts <- get_results(fit, "starts")
  expect_true(is.na(fit$n_best_replicated))
  expect_equal(starts$log_likelihood[starts$selected], fit$log_likelihood, tolerance = 1e-10)
  expect_identical(starts$kind, "correlated")
  expect_true(is.data.frame(fit$initialization_starts))
  initial <- mixture_regression(y ~ time, data, 1, id = "id", class_level = "group",
    random = "intercept", seed = 52, n_starts = 1, vcov_type = "none")
  expect_equal(fit$initialization_starts$log_likelihood[fit$initialization_starts$selected],
                initial$log_likelihood, tolerance = 1e-9)
})
