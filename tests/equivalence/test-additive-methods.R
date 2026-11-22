# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

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
  skip_without_reference("numDeriv")
  stats <- latents:::.additive_prepare(additive_methods_data, c("y1", "y2"), "group")
  classes <- c("group_class_1", "group_class_2")
  point <- list(means = rbind(c(-0.6, 0.1), c(0.9, -0.2)),
                between = rbind(c(0.2, 0.6), c(0.2, 0.6)),
                within = c(1.3, 0.7), weights = c(0.3, 0.7))
  sides <- lapply(c("varying", "equal"), \(between_variance) {
    structure <- latents:::.additive_structure("additive", between_variance)
    theta <- latents:::.additive_pack(point, structure, classes, c("y1", "y2"))
    analytic <- colSums(latents:::.additive_group_scores(stats, point, structure,
                                                         classes, c("y1", "y2")))
    numeric_gradient <- numDeriv::grad(\(v) latents:::.additive_expectation(
      stats, latents:::.additive_unpack(v, 2L, 2L, structure))$log_likelihood, theta)
    list(actual = data.frame(between_variance = between_variance, parameter = names(theta),
                             score = unname(analytic)),
         expected = data.frame(between_variance = between_variance, parameter = names(theta),
                               score = numeric_gradient))
  })
  # numDeriv's error is absolute near zero; a relative measure would judge
  # the rounding noise of a near-zero score. The original tolerance, 1e-7.
  expect_reference(do.call(rbind, lapply(sides, `[[`, "actual")),
                   do.call(rbind, lapply(sides, `[[`, "expected")), reference = "numDeriv",
                   covers = "multilpa()", quantity = "score off the optimum, additive family",
                   tolerance = 1e-7, measure = "scaled", key = c("between_variance", "parameter"),
                   note = "Package side: the internal analytic group scores (.additive_group_scores) behind the standard errors of multilpa(family = 'additive'); reference: numDeriv::grad() of the package's likelihood.")
})

test_that("observed information matches an independent numerical Hessian", {
  skip_without_reference("numDeriv")
  fit <- additive_methods_fit
  stats <- fit$sufficient_statistics
  classes <- names(fit$group_probabilities)
  theta <- coef(fit, scale = "unconstrained")
  structure <- fit$structure
  information <- latents:::.additive_information(stats, theta, structure, classes,
                                                 fit$vars, 1e-4)
  reference <- -numDeriv::jacobian(\(v) colSums(latents:::.additive_group_scores(
    stats, latents:::.additive_unpack(v, 2L, 2L, structure), structure, classes,
    fit$vars)), theta)
  dimnames(reference) <- list(names(theta), names(theta))
  note <- "Package side: the internal observed information (.additive_information, central differences of the analytic scores, step 1e-4) behind vcov() of multilpa(family = 'additive')."
  # Information entries run down to zero, where numDeriv's error is absolute:
  # scaled. The original tolerance, 1e-6.
  expect_reference(information, reference, reference = "numDeriv",
                   covers = c("multilpa()", "coef.multilpa_additive()", "vcov.multilpa_additive()"),
                   quantity = "observed information, Jacobian of the analytic scores",
                   tolerance = 1e-6, measure = "scaled",
                   note = paste(note, "Reference: numDeriv::jacobian() of the analytic scores."))
  hessian <- numDeriv::hessian(\(v) latents:::.additive_expectation(
    stats, latents:::.additive_unpack(v, 2L, 2L, structure))$log_likelihood, theta)
  dimnames(hessian) <- list(names(theta), names(theta))
  # 1e-3, the original tolerance: the package differentiates the scores with a
  # step of 1e-4, numDeriv differentiates the likelihood twice; both carry
  # truncation error. Scaled, because entries run down to zero.
  expect_reference(information, -hessian, reference = "numDeriv",
                   covers = c("multilpa()", "coef.multilpa_additive()", "vcov.multilpa_additive()"),
                   quantity = "observed information, Hessian of the likelihood",
                   tolerance = 1e-3, measure = "scaled",
                   note = paste(note, "Reference: numDeriv::hessian() of the package's likelihood."))
})
