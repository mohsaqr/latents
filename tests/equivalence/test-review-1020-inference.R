# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1020-inference.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("direct-effect inference transforms the full covariance across predictor origins", {
  set.seed(802)
  data <- data.frame(g = rep(1:20, each = 20), z = stats::rnorm(400),
                     w = stats::rnorm(400))
  data$y <- 0.8 + 0.7 * data$z - 0.4 * data$w + stats::rnorm(400)
  fit <- function(d, ...) quietly(multilpa(d, "y", "g", 1, 1,
    direct_effects = "y ~ z + w", n_starts = 1, seed = 3, tol = 1e-12, ...))
  original <- fit(data)
  shifted_data <- transform(data, z = z + 2026, w = w - 1000)
  shifted <- fit(shifted_data)
  # Coordinates are intercept, log variance, z slope, w slope.
  transform <- diag(4)
  transform[1, 3:4] <- c(-2026, 1000)
  inverse_transform <- diag(4)
  inverse_transform[1, 3:4] <- -transform[1, 3:4]
  invisible(lapply(c("observed", "robust", "opg"), function(type) {
    a <- parameter_inference(original, vcov_type = type)
    b <- parameter_inference(shifted, vcov_type = type)
    expected <- transform %*% attr(a, "covariance_unconstrained") %*% t(transform)
    expect_equal(unname(attr(b, "covariance_unconstrained")), unname(expected),
                 tolerance = 2e-7)
    expect_equal(b$standard_error[b$parameter == "direct_effect"],
                 a$standard_error[a$parameter == "direct_effect"], tolerance = 2e-7)
    expect_equal(b$p_value[b$parameter == "direct_effect"],
                 a$p_value[a$parameter == "direct_effect"], tolerance = 2e-7)
    if (type != "observed") {
      expect_equal(unname(attr(b, "group_scores")),
                   unname(attr(a, "group_scores") %*% inverse_transform), tolerance = 2e-6)
    }
    # The cached observed Hessian must work for a second covariance type too.
    expect_equal(unname(attr(b, "hessian")),
                 unname(t(inverse_transform) %*% attr(a, "hessian") %*%
                          inverse_transform), tolerance = 2e-7)
  }))
  reference <- stats::lm(y ~ z + w, shifted_data)
  mle_variance <- sum(stats::residuals(reference)^2) / nrow(data)
  reference_covariance <- summary(reference)$cov.unscaled * mle_variance
  inferred <- parameter_inference(shifted)
  # The one-class fit is this regression. 2e-7 relative, the original
  # tolerance: the fit stops at tol = 1e-12 on the likelihood.
  expect_reference(unname(attr(inferred, "covariance")[c(1, 3, 4), c(1, 3, 4)]),
                   unname(reference_covariance), reference = "stats",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "covariance of intercept and direct effects, shifted predictors",
                   tolerance = 2e-7, measure = "relative")
  expect_reference(inferred$estimate[c(1, 3, 4)], unname(stats::coef(reference)), reference = "stats",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "intercept and direct effects, shifted predictors",
                   tolerance = 2e-7, measure = "relative", key = c("(Intercept)", "z", "w"))
  # Independently form the normal-model sandwich in the original predictor
  # coordinates, then apply the known intercept transformation.
  reference <- stats::lm(y ~ z + w, data)
  residual <- stats::residuals(reference)
  variance <- mean(residual^2)
  design <- stats::model.matrix(reference)
  unit_scores <- rowsum(cbind(design * (residual / variance),
                              (residual^2 / variance - 1) / 2), data$g)
  inverse <- matrix(0, 4, 4)
  inverse[1:3, 1:3] <- summary(reference)$cov.unscaled * variance
  inverse[4, 4] <- 2 / nrow(data)
  reference_covariance <- inverse %*% crossprod(unit_scores) %*% inverse * 20 / 19
  order <- c(1, 4, 2, 3)
  reference_covariance <- transform %*% reference_covariance[order, order] %*% t(transform)
  expect_reference(unname(attr(parameter_inference(shifted, vcov_type = "robust"),
                               "covariance_unconstrained")),
                   unname(reference_covariance),
                   reference = "independent: CR1 sandwich of the normal regression, base R",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "robust covariance, unconstrained scale, shifted predictors",
                   tolerance = 2e-7, measure = "relative")
})

test_that("direct-effect charts use target-specific observations and preserve public covariates", {
  set.seed(804)
  data <- data.frame(g = rep(1:40, each = 10), z = stats::rnorm(400),
                     w = stats::rnorm(400))
  data$z[201:400] <- data$z[201:400] + 1e5
  data$a <- 0.8 + 0.7 * data$z + stats::rnorm(400)
  data$a[201:400] <- NA_real_
  data$b <- -0.3 - 0.4 * data$z + 0.2 * data$w + stats::rnorm(400)
  fit <- quietly(multilpa(data, c("a", "b"), "g", 1, 1,
    direct_effects = c("a ~ z", "b ~ z + w"), missing = "fiml", n_starts = 1,
    seed = 1, tol = 1e-12))
  expect_true(fit$converged)
  inference <- parameter_inference(fit)
  references <- list(a = stats::lm(a ~ z, data), b = stats::lm(b ~ z + w, data))
  # One row per indicator and coefficient: the intercept (the profile mean)
  # and the direct effects of that indicator, in the order lm() reports them.
  rows_of <- function(indicator) {
    which((inference$parameter == "mean" & inference$term == indicator) |
            (inference$parameter == "direct_effect" &
               startsWith(inference$term, paste0(indicator, ":"))))
  }
  coefficient_table <- function(values_of) {
    do.call(rbind, lapply(names(references), function(indicator) {
      data.frame(indicator = indicator,
                 coefficient = names(stats::coef(references[[indicator]])),
                 estimate = values_of(indicator))
    }))
  }
  covariance_table <- function(matrix_of) {
    do.call(rbind, lapply(names(references), function(indicator) {
      values <- matrix_of(indicator)
      terms <- names(stats::coef(references[[indicator]]))
      data.frame(indicator = indicator,
                 row = rep(terms, times = length(terms)),
                 column = rep(terms, each = length(terms)),
                 covariance = as.vector(values))
    }))
  }
  # Each Gaussian indicator of a one-profile fit is its own regression on the
  # rows where it is observed. 2e-6 relative, the original tolerance.
  expect_reference(coefficient_table(\(indicator) inference$estimate[rows_of(indicator)]),
                   coefficient_table(\(indicator) unname(stats::coef(references[[indicator]]))),
                   reference = "stats", covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "intercepts and direct effects by indicator, FIML",
                   tolerance = 2e-6, measure = "relative", key = c("indicator", "coefficient"))
  expect_reference(covariance_table(\(indicator) {
                     rows <- rows_of(indicator)
                     unname(attr(inference, "covariance")[rows, rows])
                   }),
                   covariance_table(\(indicator) {
                     reference <- references[[indicator]]
                     summary(reference)$cov.unscaled * mean(stats::residuals(reference)^2)
                   }),
                   reference = "stats", covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "covariance of intercepts and direct effects by indicator, FIML",
                   tolerance = 2e-6, measure = "relative", key = c("indicator", "row", "column"))
  expect_identical(colnames(fit$extra_data$direct$covariates), c("z", "w"))
  expect_identical(fit$extra_data$direct$effects$covariate_column, c(1L, 1L, 2L))
  expect_identical(latents:::.multilpa_refit_arguments(fit)$direct_effects, c("a ~ z", "b ~ z + w"))
  expect_equal(predict(fit, newdata = data, type = "posterior"),
               predict(fit, type = "posterior"), tolerance = 1e-12)
})
