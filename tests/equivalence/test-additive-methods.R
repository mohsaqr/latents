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

test_that("analytic scores equal numerical derivatives away from the optimum", {
  skip_if_not_installed("numDeriv")
  stats <- .additive_prepare(additive_methods_data, c("y1", "y2"), "group")
  classes <- c("group_class_1", "group_class_2")
  point <- list(means = rbind(c(-0.6, 0.1), c(0.9, -0.2)),
                between = rbind(c(0.2, 0.6), c(0.2, 0.6)),
                within = c(1.3, 0.7), weights = c(0.3, 0.7))
  lapply(c("varying", "equal"), \(between_variance) {
    structure <- .additive_structure("additive", between_variance)
    theta <- .additive_pack(point, structure, classes, c("y1", "y2"))
    analytic <- colSums(.additive_group_scores(stats, point, structure,
                                               classes, c("y1", "y2")))
    numeric_gradient <- numDeriv::grad(\(v) .additive_expectation(
      stats, .additive_unpack(v, 2L, 2L, structure))$log_likelihood, theta)
    expect_equal(unname(analytic), numeric_gradient, tolerance = 1e-7)
  })
})

test_that("observed information matches an independent numerical Hessian", {
  skip_if_not_installed("numDeriv")
  fit <- additive_methods_fit
  stats <- fit$sufficient_statistics
  classes <- names(fit$group_probabilities)
  theta <- coef(fit, scale = "unconstrained")
  structure <- fit$structure
  information <- .additive_information(stats, theta, structure, classes,
                                       fit$vars, 1e-4)
  reference <- -numDeriv::jacobian(\(v) colSums(.additive_group_scores(
    stats, .additive_unpack(v, 2L, 2L, structure), structure, classes,
    fit$vars)), theta)
  expect_equal(unname(information), reference, tolerance = 1e-6)
  hessian <- numDeriv::hessian(\(v) .additive_expectation(
    stats, .additive_unpack(v, 2L, 2L, structure))$log_likelihood, theta)
  expect_equal(unname(information), -hessian, tolerance = 1e-3)
})
