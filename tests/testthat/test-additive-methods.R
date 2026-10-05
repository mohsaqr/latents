# multilpa(family = "additive"): scores, information, standard errors,
# tables, S3 methods and plots.

additive_methods_truth <- list(
  means = rbind(c(-1, 0.5), c(1.2, -0.7)),
  between = rbind(c(0.3, 0.4), c(0.5, 0.2)),
  within = c(0.8, 1.1), weights = c(0.45, 0.55))
# 200 groups, as in the reference simulation cells, so no class falls under
# the weak-class threshold for want of groups.
additive_methods_data <- additive_draw(additive_methods_truth,
                                       rep(c(3L, 5L, 8L, 10L), 50L), seed = 4)
# A well-separated fit with about 100 groups per class, as in the reference
# simulation cells, where the weak-class warning must stay silent.
additive_methods_fit_large <- function() {
  truth <- list(means = rbind(c(-1, -1), c(1, 1)),
                between = rbind(c(0.25, 0.5), c(0.5, 0.25)),
                within = c(1, 1), weights = c(0.5, 0.5))
  multilpa(additive_draw(truth, rep(10L, 200), seed = 1), c("y1", "y2"),
           "group", n_group_classes = 2, family = "additive", seed = 3)
}
additive_methods_fit <- multilpa(additive_methods_data, c("y1", "y2"), "group",
                                 n_group_classes = 2, family = "additive",
                                 tol = 1e-10, seed = 2)

test_that("pack and unpack are inverse on both between-variance restrictions", {
  p <- additive_methods_truth
  equal <- modifyList(p, list(between = rbind(c(0.3, 0.4), c(0.3, 0.4))))
  classes <- c("group_class_1", "group_class_2")
  lapply(list(list(p, .additive_structure("additive", "varying")),
              list(equal, .additive_structure("additive", "equal"))), \(case) {
    theta <- .additive_pack(case[[1]], case[[2]], classes, c("y1", "y2"))
    back <- .additive_unpack(theta, 2L, 2L, case[[2]])
    expect_equal(back$means, case[[1]]$means)
    expect_equal(back$between, case[[1]]$between)
    expect_equal(back$within, matrix(case[[1]]$within, 2L, 2L, byrow = TRUE))
    expect_equal(back$weights, case[[1]]$weights)
  })
})

test_that("one class, balanced groups: the SE of the mean is exact", {
  data <- data.frame(group = rep(seq_len(30), each = 5L),
                     y = rep(seq(-2, 2, length.out = 30), each = 5L) +
                       rep(c(-1, -0.4, 0, 0.5, 0.9), 30) *
                       rep(c(1, 0.8, 1.2), length.out = 150))
  fit <- multilpa(data, "y", "group", n_group_classes = 1, family = "additive",
                  tol = 1e-13, max_iter = 5000, seed = 1)
  expect_false(any(fit$between_zero))
  table <- get_results(fit)
  mean_row <- subset(table, parameter == "mean")
  # Var(mu) = (W + n T) / (J n) for a balanced one-way random-effects model.
  exact <- sqrt((fit$within_variances[1, 1] + 5 * fit$between_variances[1, 1]) / (30 * 5))
  expect_equal(mean_row$standard_error, unname(exact), tolerance = 1e-6)
})

test_that("natural covariance is the delta-method map; names agree", {
  fit <- additive_methods_fit
  expect_identical(names(coef(fit)), rownames(vcov(fit)))
  expect_identical(names(coef(fit, scale = "unconstrained")),
                   rownames(vcov(fit, scale = "unconstrained")))
  table <- get_results(fit)
  expect_identical(nrow(table), length(coef(fit)))
  expect_equal(table$estimate, unname(coef(fit)))
  expect_equal(table$standard_error, unname(sqrt(diag(vcov(fit)))),
               tolerance = 1e-10)
  expect_equal(unname(confint(fit)), cbind(table$conf_low, table$conf_high))
  # Variance intervals are on the log scale, so they stay positive.
  expect_true(all(subset(table, parameter == "variance")$conf_low > 0))
  weights <- subset(table, parameter == "weight")
  expect_equal(sum(weights$estimate), 1)
  expect_equal(weights$standard_error[1], weights$standard_error[2],
               tolerance = 1e-8)
})

test_that("robust and OPG covariances follow their definitions", {
  fit <- additive_methods_fit
  robust <- parameter_inference(fit, vcov_type = "robust")
  scores <- attr(robust, "group_scores")
  information <- -attr(robust, "hessian")
  inverse <- solve(information)
  expect_equal(attr(robust, "covariance_unconstrained"),
               inverse %*% crossprod(scores) %*% inverse, tolerance = 1e-8,
               ignore_attr = TRUE)
  opg <- parameter_inference(fit, vcov_type = "opg")
  expect_equal(attr(opg, "covariance_unconstrained"), solve(crossprod(scores)),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_identical(nrow(scores), fit$n_groups)
  # Scores sum to (nearly) zero at a converged maximum: below the 1e-2
  # Newton decrement at which parameter_inference() warns.
  expect_lt(attr(robust, "scaled_score"), 1e-2)
})

test_that("group-class labels do not depend on the start that won", {
  a <- multilpa(additive_methods_data, c("y1", "y2"), "group",
                n_group_classes = 2, family = "additive", tol = 1e-10, seed = 2)
  b <- multilpa(additive_methods_data, c("y1", "y2"), "group",
                n_group_classes = 2, family = "additive", tol = 1e-10, seed = 77,
                n_starts = 5)
  expect_equal(coef(a), coef(b), tolerance = 1e-5)
  expect_gte(a$group_probabilities[1], a$group_probabilities[2])
  # Three classes: the reported order is by weight whatever the start.
  three <- list(means = rbind(c(-2, -2), c(0, 0), c(2, 2)),
                between = matrix(0.2, 3, 2), within = c(1, 1),
                weights = c(0.2, 0.5, 0.3))
  data <- additive_draw(three, rep(8L, 150), seed = 9)
  fits <- lapply(c(1, 5), \(seed) multilpa(data, c("y1", "y2"), "group",
                                           n_group_classes = 3, family = "additive",
                                           between_variance = "equal", seed = seed))
  expect_equal(coef(fits[[1]]), coef(fits[[2]]), tolerance = 1e-5)
  expect_true(all(diff(fits[[1]]$group_probabilities) <= 0))
})

test_that("tables have one row per unit and agree with each other", {
  fit <- additive_methods_fit
  tables <- get_results(fit, "all")
  expect_named(tables, .additive_tables())
  expect_true(all(vapply(tables, is.data.frame, logical(1))))
  groups <- tables$groups
  expect_identical(nrow(groups), fit$n_groups)
  expect_identical(groups$group, unique(additive_methods_data$group))
  expect_equal(unname(rowSums(groups[, c("probability_group_class_1",
                                         "probability_group_class_2")])),
               rep(1, fit$n_groups), tolerance = 1e-12)
  expect_identical(nrow(tables$intercepts), fit$n_groups * 2L)
  expect_identical(nrow(tables$conditional_intercepts), fit$n_groups * 2L * 2L)
  assignments <- tables$assignments
  expect_identical(nrow(assignments), nrow(additive_methods_data))
  expect_identical(assignments$group, additive_methods_data$group)
  expect_identical(assignments$group_class,
                   groups$group_class[match(assignments$group, groups$group)])
  expect_true(all(assignments$assignment_level == "group"))
  classification <- tables$classification
  expect_equal(as.vector(tapply(classification$mean_posterior,
                               classification$assigned, sum)), c(1, 1),
               tolerance = 1e-12)
  expect_equal(sum(tables$group_classes$weight), 1)
  expect_identical(sum(tables$group_classes$n_assigned), fit$n_groups)
  expect_identical(tables$fit$n_groups, nobs(fit))
  expect_equal(tables$fit$bic, BIC(fit))
  expect_equal(tables$fit$aic, AIC(fit))
  expect_identical(attr(logLik(fit), "df"), fit$n_parameters)
  # Marginal intercept mean is the posterior mix of the conditional ones.
  conditional <- tables$conditional_intercepts
  mixed <- tapply(conditional$posterior * conditional$intercept_mean,
                  paste(conditional$group, conditional$indicator), sum)
  marginal <- stats::setNames(tables$intercepts$intercept_mean,
                              paste(tables$intercepts$group,
                                    tables$intercepts$indicator))
  expect_equal(as.vector(mixed[names(marginal)]), unname(marginal),
               tolerance = 1e-12)
})

test_that("summary, print and as.data.frame return the documented objects", {
  fit <- additive_methods_fit
  summary_object <- summary(fit)
  expect_s3_class(summary_object, "summary_multilpa_additive")
  expect_named(summary_object, c("fit", "parameters", "group_classes"))
  expect_identical(as.data.frame(summary_object), get_results(fit))
  expect_identical(as.data.frame(fit), get_results(fit))
  expect_identical(as.data.frame(fit, what = "groups"), get_results(fit, "groups"))
  expect_output(print(summary_object), "Group classes")
  expect_output(print(fit), "get_results")
})

test_that("inference is refused, by class, where Wald inference does not apply", {
  balanced <- data.frame(group = rep(seq_len(6), each = 4L),
                         y = rep(seq(-0.05, 0.05, length.out = 6), each = 4L) +
                           c(-1, -0.3, 0.3, 1))
  boundary <- suppressWarnings(multilpa(balanced, "y", "group",
                                        n_group_classes = 1, family = "additive",
                                        seed = 1))
  expect_true(boundary$between_zero[1, 1])
  expect_error(vcov(boundary), class = "latents_boundary_fit")
  expect_error(parameter_inference(boundary), class = "latents_boundary_fit")
  expect_message(table <- get_results(boundary),
                 class = "latents_no_standard_errors")
  expect_true(all(is.na(table$standard_error)))
  expect_identical(table$estimate, unname(coef(boundary)))
  unconverged <- suppressWarnings(multilpa(
    additive_methods_data, c("y1", "y2"), "group", n_group_classes = 2,
    family = "additive", max_iter = 2, seed = 1))
  expect_error(vcov(unconverged), class = "latents_no_converge")
  expect_error(parameter_inference(additive_methods_fit, method = "bootstrap"),
               class = "latents_unsupported_inference")
  altered <- additive_methods_data
  altered$y1[1] <- altered$y1[1] + 1
  expect_error(parameter_inference(additive_methods_fit, data = altered),
               class = "latents_bad_inference_data")
  expect_s3_class(parameter_inference(additive_methods_fit,
                                      data = additive_methods_data), "data.frame")
})

test_that("robust covariance needs more groups than parameters", {
  truth <- list(means = matrix(c(0, 0), 1), between = matrix(c(4, 4), 1),
                within = c(1, 1), weights = 1)
  few <- additive_draw(truth, rep(6L, 5L), seed = 3)
  fit <- multilpa(few, c("y1", "y2"), "group", n_group_classes = 1,
                  family = "additive", seed = 1)
  expect_false(any(fit$between_zero))
  expect_identical(fit$n_parameters, 6L)
  expect_error(vcov(fit, type = "robust"), class = "latents_too_few_groups")
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

test_that("plots are ggplot objects with intervals only when inference applies", {
  skip_if_not_installed("ggplot2")
  fit <- additive_methods_fit
  views <- c("means", "variances", "intercepts")
  plots <- lapply(views, \(view) plot(fit, what = view))
  expect_true(all(vapply(plots, ggplot2::is_ggplot, logical(1))))
  layers <- vapply(plots[1:2], \(p) any(vapply(p$layers, \(layer)
    inherits(layer$geom, "GeomErrorbar"), logical(1))), logical(1))
  expect_true(all(layers))
  balanced <- data.frame(group = rep(seq_len(6), each = 4L),
                         y = rep(seq(-0.05, 0.05, length.out = 6), each = 4L) +
                           c(-1, -0.3, 0.3, 1))
  boundary <- suppressWarnings(multilpa(balanced, "y", "group",
                                        n_group_classes = 1, family = "additive",
                                        seed = 1))
  bare <- plot(boundary)
  expect_false(any(vapply(bare$layers, \(layer)
    inherits(layer$geom, "GeomErrorbar"), logical(1))))
  expect_match(bare$labels$caption, "No intervals")
})

test_that("Newton polishing reaches the maximum EM stops short of", {
  # A rare, weakly separated class: EM's likelihood-change rule stops with a
  # scaled score above 1e-2 here; the Newton finish must bring it to ~0.
  truth <- list(means = rbind(c(-0.4, -0.4), c(0.4, 0.4)),
                between = matrix(0.4, 2, 2), within = c(1, 1),
                weights = c(0.9, 0.1))
  data <- additive_draw(truth, rep(10L, 200), seed = 5)
  fit <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                  family = "additive", between_variance = "equal", seed = 1)
  expect_gt(fit$newton_steps, 0L)
  expect_warning(inference <- parameter_inference(fit),
                 class = "latents_weak_class")
  expect_lt(attr(inference, "scaled_score"), 1e-6)
  stats <- fit$sufficient_statistics
  structure <- .additive_structure("additive", "equal")
  em_only <- .additive_em(stats, .additive_initialize(stats, 2L, structure, 1e-6, 1L),
                          structure, 1e-6, 1000L, 1e-8, matrix(FALSE, 2, 2))
  expect_gte(fit$log_likelihood, em_only$expectation$log_likelihood - 1e-8)
})

test_that("effective groups equal class counts under perfect separation", {
  truth <- list(means = rbind(c(-6, -6), c(6, 6)), between = matrix(0.2, 2, 2),
                within = c(1, 1), weights = c(0.7, 0.3))
  data <- additive_draw(truth, rep(10L, 300), seed = 1)
  fit <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                  family = "additive", seed = 1)
  classes <- get_results(fit, "group_classes")
  expect_equal(classes$effective_groups, classes$count, tolerance = 1e-6)
  expect_equal(get_results(fit, "fit")$min_effective_groups,
               min(classes$effective_groups))
  one <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 1,
                  family = "additive", seed = 1)
  expect_true(is.na(get_results(one, "group_classes")$effective_groups))
})

test_that("a rare, weakly separated class raises latents_weak_class", {
  truth <- list(means = rbind(c(-0.4, -0.4), c(0.4, 0.4)),
                between = matrix(0.4, 2, 2), within = c(1, 1),
                weights = c(0.9, 0.1))
  data <- additive_draw(truth, rep(10L, 200), seed = 5)
  fit <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                  family = "additive", between_variance = "equal", seed = 1)
  expect_warning(get_results(fit), class = "latents_weak_class")
  expect_warning(parameter_inference(fit), class = "latents_weak_class")
  expect_lt(get_results(fit, "fit")$min_effective_groups, 50)
  # Well separated, 100 groups per class: no warning.
  expect_no_warning(get_results(additive_methods_fit_large()))
})
