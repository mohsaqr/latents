# ---- Review 2026-10-26: mixture-regression mathematics ----
#
# Calibration: the one-class limit is the ordinary GLM (glm(), MASS::glm.nb(),
# MASS::polr()), and the observed-information standard errors agree with a
# numerical Hessian of the log likelihood. Invariants: class relabelling,
# translating or rescaling a predictor. Regressions: the large-origin Gaussian
# least-squares shortcut with an empty or rank-deficient class, and an aliased
# constant column whose rounded mean is not exact.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

.review_mixreg_data <- function(n = 240L, seed = 2026L) {
  .latents_local_seed(seed)
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

.review_mixreg_se <- function(fit) sqrt(diag(fit$inference$vcov))

test_that("one-class gaussian, binomial and poisson fits are glm()", {
  frame <- .review_mixreg_data()
  n <- nrow(frame)
  gaussian_glm <- stats::glm(y ~ x + z, data = frame)
  gaussian_fit <- mixture_regression(y ~ x + z, frame, 1, n_starts = 0)
  expect_equal(drop(gaussian_fit$params$beta), stats::coef(gaussian_glm),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(gaussian_fit$log_likelihood, as.numeric(stats::logLik(gaussian_glm)),
               tolerance = 1e-10)
  # The maximum-likelihood variance divides by n, glm()'s dispersion by n - p.
  expect_equal(unname(.review_mixreg_se(gaussian_fit)[1:3]),
               unname(sqrt(diag(stats::vcov(gaussian_glm)) * (n - 3) / n)),
               tolerance = 1e-8)

  binomial_glm <- stats::glm(cbind(successes, failures) ~ x, stats::binomial, frame)
  binomial_fit <- mixture_regression(cbind(successes, failures) ~ x, frame, 1,
                                     family = "binomial", n_starts = 0)
  expect_equal(drop(binomial_fit$params$beta), stats::coef(binomial_glm),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(binomial_fit$log_likelihood, as.numeric(stats::logLik(binomial_glm)),
               tolerance = 1e-10)
  expect_equal(unname(.review_mixreg_se(binomial_fit)),
               unname(sqrt(diag(stats::vcov(binomial_glm)))), tolerance = 1e-5)

  poisson_glm <- stats::glm(count ~ x + offset(log(exposure)), stats::poisson, frame)
  poisson_fit <- mixture_regression(count ~ x + offset(log(exposure)), frame, 1,
                                    family = "poisson", n_starts = 0)
  expect_equal(drop(poisson_fit$params$beta), stats::coef(poisson_glm),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(poisson_fit$log_likelihood, as.numeric(stats::logLik(poisson_glm)),
               tolerance = 1e-10)
  expect_equal(unname(.review_mixreg_se(poisson_fit)),
               unname(sqrt(diag(stats::vcov(poisson_glm)))), tolerance = 1e-5)
})

test_that("one-class negative-binomial and ordinal fits are glm.nb() and polr()", {
  skip_if_not_installed("MASS")
  frame <- .review_mixreg_data()
  .latents_local_seed(5L)
  frame$overdispersed <- MASS::rnegbin(nrow(frame), exp(0.5 + 0.4 * frame$x), theta = 2)
  reference <- MASS::glm.nb(overdispersed ~ x, frame,
                            control = stats::glm.control(epsilon = 1e-12, maxit = 100))
  fit <- mixture_regression(overdispersed ~ x, frame, 1, family = "negative_binomial",
                            n_starts = 0)
  expect_equal(drop(fit$params$beta), stats::coef(reference), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(fit$params$sigma2, 1 / reference$theta, tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)), tolerance = 1e-9)

  ordinal_reference <- MASS::polr(ordinal ~ x + z, frame, Hess = TRUE,
                                  control = list(reltol = 1e-14))
  ordinal_fit <- mixture_regression(ordinal ~ x + z, frame, 1, family = "ordinal",
                                    n_starts = 0)
  expect_equal(drop(ordinal_fit$params$beta), stats::coef(ordinal_reference),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(drop(ordinal_fit$params$thresholds), ordinal_reference$zeta,
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(ordinal_fit$log_likelihood, as.numeric(stats::logLik(ordinal_reference)),
               tolerance = 1e-9)
})

test_that("observed-information errors equal a numerical Hessian of the likelihood", {
  skip_if_not_installed("numDeriv")
  frame <- .review_mixreg_data()
  check <- function(fit) {
    spec <- fit$spec
    params <- fit$params
    theta <- .mixture_pack(spec, params)
    log_likelihood <- function(point) {
      .mixture_expectation(spec, .mixture_unpack(spec, point, params))$log_likelihood
    }
    hessian <- numDeriv::hessian(log_likelihood, theta)
    expect_length(theta, fit$n_parameters)
    expect_equal(unname(.review_mixreg_se(fit)), sqrt(diag(solve(-hessian))),
                 tolerance = 1e-5)
  }
  check(mixture_regression(y ~ x, frame, 2, membership = ~z, n_starts = 2, seed = 1))
  check(mixture_regression(count ~ x, frame, 2, family = "poisson", n_starts = 2,
                           seed = 1))
  check(mixture_regression(ordinal ~ x, frame, 2, family = "ordinal", n_starts = 2,
                           seed = 1))
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
  expect_equal(unname(fit$inference$vcov), unname(sandwich), tolerance = 1e-5)
})

test_that("relabelling the classes leaves the likelihood and the fitted means unchanged", {
  frame <- .review_mixreg_data()
  fit <- mixture_regression(y ~ x, frame, 2, membership = ~z, n_starts = 2, seed = 1,
                            vcov_type = "none")
  swapped <- .mixture_permute(fit$spec, fit$params, c(2L, 1L))
  original <- .mixture_expectation(fit$spec, fit$params)
  relabelled <- .mixture_expectation(fit$spec, swapped)
  expect_equal(relabelled$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  expect_equal(relabelled$tau, original$tau[, 2:1], tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(.mixture_class_means(fit$spec, swapped),
               .mixture_class_means(fit$spec, fit$params)[, 2:1], ignore_attr = TRUE)
})

test_that("translating a predictor by 1e7 changes only the intercepts", {
  frame <- .review_mixreg_data()
  frame$shifted <- frame$x + 1e7
  compare <- function(formula, shifted_formula, family) {
    near <- mixture_regression(formula, frame, 2, family = family, n_starts = 2,
                               seed = 1, vcov_type = "none")
    far <- mixture_regression(shifted_formula, frame, 2, family = family,
                              n_starts = 2, seed = 1, vcov_type = "none")
    expect_equal(far$log_likelihood, near$log_likelihood, tolerance = 1e-9)
    expect_equal(unname(far$params$beta[2L, ]), unname(near$params$beta[2L, ]),
                 tolerance = 1e-8)
    expect_equal(.mixture_linear_predictors(far$spec, far$params),
                 .mixture_linear_predictors(near$spec, near$params), tolerance = 1e-8)
  }
  compare(y ~ x, y ~ shifted, "gaussian")
  compare(count ~ x, count ~ shifted, "poisson")
})

test_that("rescaling a predictor rescales its slopes and nothing else", {
  frame <- .review_mixreg_data()
  frame$scaled <- frame$x * 1e6
  compare <- function(response, family) {
    unit <- mixture_regression(stats::reformulate("x", response), frame, 2,
                               family = family, n_starts = 2, seed = 1,
                               vcov_type = "none")
    large <- mixture_regression(stats::reformulate("scaled", response), frame, 2,
                                family = family, n_starts = 2, seed = 1,
                                vcov_type = "none")
    expect_equal(large$log_likelihood, unit$log_likelihood, tolerance = 1e-9)
    slope <- function(fit, term) unname(fit$params$beta[term, ])
    expect_equal(slope(large, "scaled") * 1e6, slope(unit, "x"), tolerance = 1e-6)
  }
  compare("y", "gaussian")
  compare("count", "poisson")
  # The ordinal M-step solved a raw-ridge system whose ridge, scaled by the
  # slope's curvature, swamped the thresholds' and stalled short of the MLE.
  compare("ordinal", "ordinal")
})

test_that("the large-origin least-squares step keeps an empty or tiny class solvable", {
  frame <- .review_mixreg_data(n = 60L)
  frame$shifted <- frame$x + 1e7
  spec <- mixture_regression(y ~ shifted, frame, 2, n_starts = 0, vcov_type = "none")$spec
  start <- .mixture_null_params(spec)
  empty <- cbind(rep(1, spec$n), rep(0, spec$n))
  updated <- .mixture_update_regression(spec, start, empty)
  expect_true(all(is.finite(updated$beta)))
  expect_equal(unname(updated$beta[, 2L]), c(0, 0))
  reference <- stats::lm(y ~ x, frame)
  expect_equal(unname(updated$beta[2L, 1L]), unname(stats::coef(reference)[2L]),
               tolerance = 1e-8)
  # One row cannot determine an intercept and a slope.
  single <- cbind(rep(1, spec$n), c(1, rep(0, spec$n - 1L)))
  expect_true(all(is.finite(.mixture_update_regression(spec, start, single)$beta)))
})

test_that("an overfitted large-origin gaussian mixture fails only with a classed condition", {
  .latents_local_seed(17L)
  frame <- data.frame(x = round(stats::rnorm(20), 1) + 1e7)
  frame$y <- round(1 + 0.5 * (frame$x - 1e7) + stats::rnorm(20))
  outcome <- tryCatch(
    suppressWarnings(mixture_regression(y ~ x, frame, 4, n_starts = 10, seed = 17,
                                        vcov_type = "none")),
    latents_no_valid_start = function(condition) "no valid start")
  expect_true(identical(outcome, "no valid start") ||
                inherits(outcome, "latents_mixture_regression"))
})

test_that("an aliased constant column with an inexact mean adds nothing", {
  frame <- .review_mixreg_data(n = 99L)
  frame$constant <- 0.1
  # The premise: the column mean rounds, so a spread-is-zero test misses it.
  column <- as.matrix(frame["constant"])
  expect_gt(max(abs(column - colMeans(column))), 0)
  for_family <- function(response, family) {
    plain <- mixture_regression(stats::reformulate("x", response), frame, 1,
                                family = family, n_starts = 0, vcov_type = "none")
    aliased <- mixture_regression(stats::reformulate(c("x", "constant"), response),
                                  frame, 1, family = family, n_starts = 0,
                                  vcov_type = "none")
    expect_equal(aliased$log_likelihood, plain$log_likelihood, tolerance = 1e-9)
    expect_equal(unname(aliased$params$beta["x", ]), unname(plain$params$beta["x", ]),
                 tolerance = 1e-8)
  }
  for_family("y", "gaussian")
  for_family("count", "poisson")
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
  expect_equal(predict(fit, new_rows)$fitted, rowSums(priors * means),
               ignore_attr = TRUE)
  by_class <- predict(fit, new_rows, type = "class_response")
  expect_equal(by_class$prior, as.vector(priors), ignore_attr = TRUE)
  densities <- vapply(1:2, function(k) {
    stats::dnorm(new_rows$y, means[, k], sqrt(fit$params$sigma2[k]))
  }, numeric(3))
  posterior <- priors * densities / rowSums(priors * densities)
  scored <- predict(fit, new_rows, type = "posterior")
  expect_equal(scored$probability_class_2, posterior[, 2], tolerance = 1e-10)
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
    .mixture_regression_objective(spec, point, tau)
  }
  updated <- .mixture_update_dispersion(spec, params, tau)
  expect_gte(objective(updated$sigma2), objective(params$sigma2))
  interior <- stats::optimize(objective, c(1, 1000), maximum = TRUE)
  expect_equal(objective(updated$sigma2), interior$objective, tolerance = 1e-6)
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
  slope <- original$term == "hours"
  expect_equal(shifted$std_error[slope], original$std_error[slope], tolerance = 1e-5)
  # One class reproduces the closed-form ML errors at the shifted origin.
  one <- mixture_regression(score ~ hours, transform(study_hours, hours = hours + 2000),
                            n_classes = 1, n_starts = 0, seed = 1)
  reference <- stats::lm(score ~ I(hours + 2000), study_hours)
  design <- stats::model.matrix(reference)
  closed_form <- sqrt(diag(mean(stats::residuals(reference)^2) * solve(crossprod(design))))
  expect_equal(as.data.frame(get_results(one, "coefficients"))$std_error,
               unname(closed_form), tolerance = 1e-6)
})
