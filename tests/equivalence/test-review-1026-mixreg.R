# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1026-mixreg.R,
# so that they run in CI's equivalence job rather than under R CMD check.

.review_mixreg_data <- function(n = 240L, seed = 2026L) {
  latents:::.latents_local_seed(seed)
  frame <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  class_two <- stats::rbinom(n, 1L, stats::plogis(0.3 + 0.8 * frame$z)) == 1L
  frame$y <- ifelse(class_two, -1 - frame$x, 1 + 2 * frame$x) +
    stats::rnorm(n, sd = ifelse(class_two, 1, 0.5))
  frame$count <- stats::rpois(n, exp(ifelse(class_two, 1.5 - 0.3 * frame$x,
                                             0.2 + 0.8 * frame$x)))
  frame$exposure <- stats::runif(n, 1, 3)
  frame$successes <- stats::rbinom(n, 5L, stats::plogis(-0.2 + 0.7 * frame$x))
  frame$failures <- 5L - frame$successes
  frame$ordinal <- cut(ifelse(class_two, -1.5 * frame$x, 2 * frame$x) +
                         stats::rlogis(n), c(-Inf, -1, 0.5, 1.5, Inf),
                       ordered_result = TRUE)
  frame
}

.review_mixreg_se <- function(fit) sqrt(diag(vcov(fit)))

# One class's regression entries of a vector named like coef(), renamed to the
# bare term labels glm() uses.
.review_class_terms <- function(values, k = 1L) {
  prefix <- sprintf("coefficient.class_%d.", k)
  kept <- values[startsWith(names(values), prefix)]
  stats::setNames(unname(kept), substring(names(kept), nchar(prefix) + 1L))
}

test_that("one-class gaussian, binomial and poisson fits are glm()", {
  frame <- .review_mixreg_data()
  n <- nrow(frame)
  # A one-class fit against glm(): coefficients, log likelihood, standard errors.
  against_glm <- function(fit, reference, family, coefficient_tolerance, error_tolerance,
                          error_factor = 1) {
    expect_reference(.review_class_terms(coef(fit)), stats::coef(reference),
                     reference = "stats",
                     covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                     quantity = sprintf("coefficients, one %s class", family),
                     tolerance = coefficient_tolerance, measure = "relative")
    expect_reference(as.numeric(logLik(fit)), as.numeric(stats::logLik(reference)),
                     reference = "stats",
                     covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                     quantity = sprintf("log likelihood, one %s class", family),
                     tolerance = 1e-10, measure = "relative")
    expect_reference(.review_class_terms(.review_mixreg_se(fit)),
                     sqrt(diag(stats::vcov(reference)) * error_factor),
                     reference = "stats",
                     covers = c("mixture_regression()", "vcov.latents_mixture_regression()"),
                     quantity = sprintf("standard errors, one %s class", family),
                     tolerance = error_tolerance, measure = "relative")
  }
  # The maximum-likelihood variance divides by n, glm()'s dispersion by n - p.
  against_glm(mixture_regression(y ~ x + z, frame, 1, n_starts = 0),
              stats::glm(y ~ x + z, data = frame), "gaussian",
              coefficient_tolerance = 1e-10, error_tolerance = 1e-8,
              error_factor = (n - 3) / n)
  # Binomial and Poisson: IRLS stops at glm()'s deviance criterion, and the
  # observed information is differenced from analytic scores; the original
  # tolerances, 1e-8 and 1e-5.
  against_glm(mixture_regression(cbind(successes, failures) ~ x, frame, 1,
                                 family = "binomial", n_starts = 0),
              stats::glm(cbind(successes, failures) ~ x, stats::binomial, frame), "binomial",
              coefficient_tolerance = 1e-8, error_tolerance = 1e-5)
  against_glm(mixture_regression(count ~ x + offset(log(exposure)), frame, 1,
                                 family = "poisson", n_starts = 0),
              stats::glm(count ~ x + offset(log(exposure)), stats::poisson, frame),
              "poisson with an offset",
              coefficient_tolerance = 1e-8, error_tolerance = 1e-5)
})

test_that("one-class negative-binomial and ordinal fits are glm.nb() and polr()", {
  skip_without_reference("MASS")
  frame <- .review_mixreg_data()
  latents:::.latents_local_seed(5L)
  frame$overdispersed <- MASS::rnegbin(nrow(frame), exp(0.5 + 0.4 * frame$x), theta = 2)
  reference <- MASS::glm.nb(overdispersed ~ x, frame,
                            control = stats::glm.control(epsilon = 1e-12, maxit = 100))
  fit <- mixture_regression(overdispersed ~ x, frame, 1, family = "negative_binomial",
                            n_starts = 0)
  # 1e-6, the original tolerance: glm.nb alternates IRLS and a theta step.
  expect_reference(.review_class_terms(coef(fit)), stats::coef(reference), reference = "MASS",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "coefficients, one negative-binomial class", tolerance = 1e-6,
                   measure = "relative")
  # latents reports log(dispersion); glm.nb reports theta = 1 / dispersion.
  expect_reference(exp(coef(fit)[["log_dispersion.class_1"]]), 1 / reference$theta,
                   reference = "MASS",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "dispersion 1 / theta, one class", tolerance = 1e-6,
                   measure = "relative")
  expect_reference(as.numeric(logLik(fit)), as.numeric(stats::logLik(reference)),
                   reference = "MASS",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one negative-binomial class", tolerance = 1e-9,
                   measure = "relative")

  ordinal_reference <- MASS::polr(ordinal ~ x + z, frame, Hess = TRUE,
                                  control = list(reltol = 1e-14))
  ordinal_fit <- mixture_regression(ordinal ~ x + z, frame, 1, family = "ordinal",
                                    n_starts = 0)
  ordinal_terms <- as.data.frame(get_results(ordinal_fit, "coefficients"))
  # Scaled: the outcome is generated without z, so its coefficient (and x's,
  # which mixes two classes of opposite sign) is near zero, where a relative
  # measure would judge the optimisers' stopping noise. 1e-6, the original.
  expect_reference(subset(ordinal_terms, !startsWith(term, "threshold:")),
                   data.frame(term = names(stats::coef(ordinal_reference)),
                              estimate = unname(stats::coef(ordinal_reference))),
                   reference = "MASS",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()"),
                   quantity = "coefficients, one ordinal class", tolerance = 1e-6,
                   measure = "scaled", key = "term")
  expect_reference(subset(ordinal_terms, startsWith(term, "threshold:")),
                   data.frame(term = paste0("threshold:", names(ordinal_reference$zeta)),
                              estimate = unname(ordinal_reference$zeta)),
                   reference = "MASS",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()"),
                   quantity = "thresholds, one ordinal class", tolerance = 1e-6,
                   measure = "relative", key = "term")
  expect_reference(as.numeric(logLik(ordinal_fit)),
                   as.numeric(stats::logLik(ordinal_reference)), reference = "MASS",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one ordinal class", tolerance = 1e-9,
                   measure = "relative")
})

test_that("observed-information errors equal a numerical Hessian of the likelihood", {
  skip_without_reference("numDeriv")
  frame <- .review_mixreg_data()
  check <- function(fit, setting) {
    spec <- fit$spec
    params <- fit$params
    theta <- latents:::.mixture_pack(spec, params)
    log_likelihood <- function(point) {
      latents:::.mixture_expectation(spec, latents:::.mixture_unpack(spec, point, params))$log_likelihood
    }
    hessian <- numDeriv::hessian(log_likelihood, theta)
    expect_length(theta, fit$n_parameters)
    # 1e-5, the original tolerance: Richardson extrapolation of second differences.
    expect_reference(.review_mixreg_se(fit), sqrt(diag(solve(-hessian))),
                     reference = "independent: numDeriv derivative of the package's own likelihood (self-consistency)",
                     covers = c("mixture_regression()", "vcov.latents_mixture_regression()"),
                     quantity = sprintf("standard errors, %s", setting), tolerance = 1e-5,
                     measure = "relative",
                     note = "Reference: the inverse of numDeriv::hessian() of the package's likelihood at the estimate.")
  }
  check(mixture_regression(y ~ x, frame, 2, membership = ~z, n_starts = 2, seed = 1),
        "two gaussian classes with a membership covariate")
  check(mixture_regression(count ~ x, frame, 2, family = "poisson", n_starts = 2,
                           seed = 1),
        "two poisson classes")
  check(mixture_regression(ordinal ~ x, frame, 2, family = "ordinal", n_starts = 2,
                           seed = 1),
        "two ordinal classes")
})

test_that("the robust covariance of a one-class Poisson fit is the GLM sandwich", {
  frame <- .review_mixreg_data()
  reference <- stats::glm(count ~ x, stats::poisson, frame)
  fit <- mixture_regression(count ~ x, frame, 1, family = "poisson", n_starts = 0,
                            vcov_type = "robust")
  design <- stats::model.matrix(reference)
  scores <- design * (frame$count - stats::fitted(reference))
  n <- nrow(frame)
  bread <- stats::vcov(reference)
  sandwich <- bread %*% crossprod(scores) %*% bread * n / (n - 1)
  robust <- vcov(fit)
  dimnames(robust) <- rep(list(sub("coefficient.class_1.", "", rownames(robust), fixed = TRUE)), 2L)
  # 1e-5, the original tolerance: glm()'s IRLS stops at its deviance criterion.
  expect_reference(robust, sandwich,
                   reference = "independent: GLM sandwich from glm()'s bread and score outer products",
                   covers = c("mixture_regression()", "vcov.latents_mixture_regression()"),
                   quantity = "robust covariance, one poisson class", tolerance = 1e-5,
                   measure = "relative")
})

test_that("the large-origin least-squares step keeps an empty or tiny class solvable", {
  frame <- .review_mixreg_data(n = 60L)
  frame$shifted <- frame$x + 1e7
  spec <- mixture_regression(y ~ shifted, frame, 2, n_starts = 0, vcov_type = "none")$spec
  start <- latents:::.mixture_null_params(spec)
  empty <- cbind(rep(1, spec$n), rep(0, spec$n))
  updated <- latents:::.mixture_update_regression(spec, start, empty)
  expect_true(all(is.finite(updated$beta)))
  expect_equal(unname(updated$beta[, 2L]), c(0, 0))
  reference <- stats::lm(y ~ x, frame)
  expect_reference(unname(updated$beta[2L, 1L]), unname(stats::coef(reference)[2L]),
                   reference = "stats", covers = "mixture_regression()",
                   quantity = "slope of the least-squares step at origin 1e7", tolerance = 1e-8,
                   measure = "relative",
                   note = "Package side: the internal regression M step (.mixture_update_regression) of mixture_regression() with every row in class 1; reference: stats::lm() on the unshifted predictor.")
  # One row cannot determine an intercept and a slope.
  single <- cbind(rep(1, spec$n), c(1, rep(0, spec$n - 1L)))
  expect_true(all(is.finite(latents:::.mixture_update_regression(spec, start, single)$beta)))
})

test_that("predictions use prior class probabilities and posteriors reproduce the fit", {
  frame <- .review_mixreg_data()
  fit <- mixture_regression(y ~ x, frame, 2, membership = ~z, n_starts = 2, seed = 1,
                            vcov_type = "none")
  new_rows <- data.frame(x = c(-1, 0, 2), z = c(0.5, -1, 0), y = c(0, 1, 3))
  gamma <- fit$params$gamma
  logits <- cbind(1, new_rows$z) %*% gamma
  priors <- exp(logits) / rowSums(exp(logits))
  means <- cbind(1, new_rows$x) %*% fit$params$beta
  reference <- "independent: class priors, means and posteriors from the estimates, base R"
  covers <- c("mixture_regression()", "predict.latents_mixture_regression()")
  # testthat's default tolerance, the original test's.
  expect_reference(predict(fit, new_rows),
                   data.frame(row = 1:3, fitted = rowSums(priors * means)),
                   reference = reference, covers = covers,
                   quantity = "prior-weighted prediction", tolerance = sqrt(.Machine$double.eps),
                   measure = "relative", key = "row")
  expect_reference(predict(fit, new_rows, type = "class_response"),
                   data.frame(row = rep(1:3, 2), class = rep(c("class_1", "class_2"), each = 3),
                              prior = as.vector(priors)),
                   reference = reference, covers = covers,
                   quantity = "prior class probabilities", tolerance = sqrt(.Machine$double.eps),
                   measure = "relative", key = c("row", "class"))
  densities <- vapply(1:2, function(k) {
    stats::dnorm(new_rows$y, means[, k], sqrt(fit$params$sigma2[k]))
  }, numeric(3))
  posterior <- priors * densities / rowSums(priors * densities)
  expect_reference(predict(fit, new_rows, type = "posterior"),
                   data.frame(row = 1:3, probability_class_2 = posterior[, 2]),
                   reference = reference, covers = covers,
                   quantity = "posterior of class 2", tolerance = 1e-10,
                   measure = "relative", key = "row")
  # Two routes to one result inside latents: a unit check.
  training <- predict(fit, frame, type = "posterior")
  expect_equal(training$probability_class_1, unname(fit$expectation$tau[, 1]),
               tolerance = 1e-10)
})

test_that("the negative-binomial dispersion step never falls back to a lower Poisson limit", {
  # A weighted profile in alpha with a local maximum at zero (non-positive
  # score there) and a higher interior maximum near alpha = 26.8.
  y <- c(0, 0, 16, 0, 0, 6, 0)
  mu <- c(0.261227, 0.908894, 15.731, 1.18125, 2.08141, 3.39144, 0.352174)
  weight <- c(0.0263183, 0.570606, 0.205777, 0.938473, 0.189234, 0.0251619, 0.88033)
  expect_lte(sum(weight * ((y - mu)^2 - y)), 0)
  spec <- list(family = "negative_binomial", variance = "varying", n = length(y),
               n_classes = 1L, y = y, trials = rep(1, length(y)),
               x = matrix(1, length(y), 1L, dimnames = list(NULL, "(Intercept)")),
               z = matrix(0, length(y), 0L), offset = log(mu))
  params <- list(beta = matrix(0, 1L, 1L), common = numeric(), sigma2 = 20)
  tau <- matrix(weight, ncol = 1L)
  objective <- function(alpha) {
    point <- params
    point$sigma2 <- alpha
    latents:::.mixture_regression_objective(spec, point, tau)
  }
  updated <- latents:::.mixture_update_dispersion(spec, params, tau)
  expect_gte(objective(updated$sigma2), objective(params$sigma2))
  # stats::optimize() as an independent maximiser of the same profile;
  # 1e-6, the original tolerance (optimize's own default tolerance).
  interior <- stats::optimize(objective, c(1, 1000), maximum = TRUE)
  expect_reference(objective(updated$sigma2), interior$objective,
                   reference = "independent: stats optimiser on the package's own objective (self-consistency)", covers = "mixture_regression()",
                   quantity = "weighted profile at the dispersion step's alpha", tolerance = 1e-6,
                   measure = "relative",
                   note = "Package side: the internal dispersion M step (.mixture_update_dispersion) of mixture_regression(family = 'negative_binomial'); reference: stats::optimize() of the package's weighted objective over alpha in [1, 1000].")
})

test_that("standard errors do not depend on a predictor's origin (MR-F5, MR-F6)", {
  fit_at <- function(shift) {
    mixture_regression(score ~ hours, transform(study_hours, hours = hours + shift),
                       n_classes = 2, n_starts = 10, seed = 1, tol = 1e-12,
                       max_iter = 5000, vcov_type = "observed")
  }
  coefficients <- function(fit) {
    table <- as.data.frame(get_results(fit, "coefficients"))
    slopes <- table[table$term == "hours", ]
    table[order(match(table$class, slopes$class[order(slopes$estimate)]), table$term), ]
  }
  original <- coefficients(fit_at(0))
  # A calendar year: every standard error used to be NA, even with one class.
  shifted <- coefficients(fit_at(2000))
  expect_true(all(is.finite(shifted$std_error)))
  # The same model at two origins: a unit check of invariance.
  slope <- original$term == "hours"
  expect_equal(shifted$std_error[slope], original$std_error[slope], tolerance = 1e-5)
  # One class reproduces the closed-form ML errors at the shifted origin.
  one <- mixture_regression(score ~ hours, transform(study_hours, hours = hours + 2000),
                            n_classes = 1, n_starts = 0, seed = 1)
  reference <- stats::lm(score ~ I(hours + 2000), study_hours)
  design <- stats::model.matrix(reference)
  closed_form <- sqrt(diag(mean(stats::residuals(reference)^2) * solve(crossprod(design))))
  expect_reference(as.data.frame(get_results(one, "coefficients")),
                   data.frame(term = c("(Intercept)", "hours"), std_error = unname(closed_form)),
                   reference = "independent: closed-form ML standard errors, sigma^2 (X'X)^-1 from lm() residuals",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()"),
                   quantity = "standard errors, one class at origin 2000", tolerance = 1e-6,
                   measure = "relative", key = "term")
})
