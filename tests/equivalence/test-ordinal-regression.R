# Ordinal regression mixtures, mixture_regression(family = "ordinal") (phase
# E7 of validation/ENGINE_DESIGN.md): a proportional-odds (cumulative logit)
# regression within each class under the regression structures. References:
# one class is MASS::polr(method = "logistic"); a likelihood written with
# plogis() differences alone; numerical derivatives of the likelihood;
# recovery of simulated truths; two categories are the logistic mixture.

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

# The mixture likelihood from its definition: class shares times
# plogis(t_y - eta) - plogis(t_y-1 - eta), sharing no code with the engine.
ordinal_dense <- function(fit, data, design) {
  params <- fit$params
  shares <- exp(params$gamma[1L, ]) / sum(exp(params$gamma[1L, ]))
  y <- as.integer(data$y)
  sum(log(rowSums(vapply(seq_len(fit$spec$n_classes), function(k) {
    eta <- drop(design %*% c(params$beta[, k], params$common))
    bounds <- c(-Inf, params$thresholds[, k], Inf)
    shares[k] * (stats::plogis(bounds[y + 1L] - eta) - stats::plogis(bounds[y] - eta))
  }, numeric(nrow(data))))))
}

test_that("one class is MASS::polr(method = \"logistic\")", {
  skip_if_not_installed("MASS")
  data <- ordinal_rows()
  data$o <- 0.3 * data$w
  fit <- mixture_regression(y ~ x + w + offset(o), data, n_classes = 1,
                            family = "ordinal", n_starts = 1, seed = 1)
  reference <- MASS::polr(y ~ x + w + offset(o), data = data, method = "logistic",
                          Hess = TRUE, control = list(reltol = 1e-14, maxit = 1000))
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)),
               tolerance = 1e-8)
  expect_equal(unname(fit$params$thresholds[, 1L]), unname(reference$zeta),
               tolerance = 1e-5)
  expect_equal(unname(fit$params$beta[, 1L]), unname(stats::coef(reference)),
               tolerance = 1e-5)
  # Both are the observed information at the same estimate, polr's by
  # optimHess() on its analytic gradient and ours by differencing the
  # analytic scores, so the errors agree to differencing error.
  coefficients <- get_results(fit, "coefficients")
  expect_identical(coefficients$term,
                   c("threshold:never|rarely", "threshold:rarely|often",
                     "threshold:often|always", "x", "w"))
  expect_equal(coefficients$std_error,
               unname(sqrt(diag(stats::vcov(reference)))[
                 c(names(reference$zeta), names(stats::coef(reference)))]),
               tolerance = 1e-3)
})

test_that("the analytic scores are the gradient of the log likelihood", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  check <- function(fit) {
    theta <- latents:::.mixture_pack(fit$spec, fit$params) + 0.03
    value <- function(v) {
      latents:::.mixture_expectation(
        fit$spec, latents:::.mixture_unpack(fit$spec, v, fit$params))$log_likelihood
    }
    params <- latents:::.mixture_unpack(fit$spec, theta, fit$params)
    scores <- latents:::.mixture_unit_scores(
      fit$spec, params, latents:::.mixture_expectation(fit$spec, params))
    expect_equal(unname(colSums(scores)), numDeriv::grad(value, theta), tolerance = 1e-6)
    # The packed coordinates round-trip to ordered thresholds.
    expect_equal(latents:::.mixture_pack(fit$spec, params), theta, tolerance = 1e-12)
  }
  check(mixture_regression(y ~ x + w, ordinal_rows(), n_classes = 2, family = "ordinal",
                           common = ~ w, n_starts = 2, seed = 1))
  check(mixture_regression(y ~ wave, ordinal_persons(), n_classes = 2,
                           family = "ordinal", id = "person", class_level = "group",
                           n_starts = 2, seed = 1))
})
