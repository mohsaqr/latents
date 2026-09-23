# Inference for a fit that holds part of its measurement model fixed.
#
# The references here are deliberately independent of the package: the
# conditional log likelihood is rebuilt by exhaustive enumeration of the latent
# assignments with direct dnorm arithmetic (helper-independent-likelihood.R),
# the held blocks are pasted back in by hand, and the information matrix comes
# from second central differences of that likelihood. Nothing in the reference
# path calls .multilpa_score(), .multilpa_decode() or optimHess().

#' Two-level data with two well-separated profiles and two group classes
#' @return A list with the data frame, the indicator matrix and the group vector.
fixed_measurement_data <- function() {
  set.seed(4021)
  n_groups <- 24L
  size <- 6L
  dat <- data.frame(g = rep(seq_len(n_groups), each = size))
  class_of_group <- rep(c(1, 2), length.out = n_groups)
  probability <- c(0.25, 0.75)[class_of_group][dat$g]
  profile <- rbinom(nrow(dat), 1, probability)
  dat$a <- rnorm(nrow(dat), 5 * profile, 1)
  dat$b <- rnorm(nrow(dat), 3 * profile, 1)
  list(data = dat, x = as.matrix(dat[c("a", "b")]), group = dat$g)
}

#' The conditional log likelihood of the two-profile, two-class design
#'
#' Built from held means and free coordinates without touching the package's
#' encoder or decoder, so it can calibrate both.
#'
#' @param means Profiles-by-indicators matrix of held means.
#' @param log_variances Four log variances, profile major.
#' @param profile_logit One logit per group class.
#' @param group_logit The single group-class logit.
#' @param x Indicator matrix.
#' @param group Group identifiers in row order.
#' @return The log likelihood as a single number.
enumerated_conditional_loglik <- function(means, log_variances, profile_logit,
                                          group_logit, x, group) {
  stopifnot("`means` must be two profiles by two indicators" =
              is.matrix(means) && identical(dim(means), c(2L, 2L)),
            "`log_variances` must give one value per profile and indicator" =
              is.numeric(log_variances) && length(log_variances) == 4L,
            "`profile_logit` must give one logit per group class" =
              is.numeric(profile_logit) && length(profile_logit) == 2L,
            "`group_logit` must be a single logit" =
              is.numeric(group_logit) && length(group_logit) == 1L)
  variances <- matrix(exp(log_variances), 2L, 2L, byrow = TRUE)
  profile_logits <- cbind(matrix(profile_logit, 2L, 1L), 0)
  group_logits <- c(group_logit, 0)
  enumerate_latent_assignments(x, group, list(
    means = means, variances = variances,
    profile_probabilities = exp(profile_logits) / rowSums(exp(profile_logits)),
    group_probabilities = exp(group_logits) / sum(exp(group_logits))
  ))$log_likelihood
}

#' Second central differences of a scalar function
#' @param objective Function of a numeric vector returning one number.
#' @param point Where to differentiate.
#' @param step Positive finite-difference step.
#' @return A symmetric numerical Hessian of `objective` at `point`.
central_difference_hessian <- function(objective, point, step = 5e-4) {
  stopifnot("`point` must be numeric" = is.numeric(point),
            "`step` must be a single positive number" =
              is.numeric(step) && length(step) == 1L && step > 0)
  base <- objective(point)
  pairs <- expand.grid(row = seq_along(point), column = seq_along(point))
  values <- vapply(seq_len(nrow(pairs)), function(index) {
    row <- pairs$row[index]
    column <- pairs$column[index]
    if (column > row) return(NA_real_)
    shift <- function(signs) {
      moved <- point
      moved[c(row, column)] <- moved[c(row, column)] + step * signs
      objective(moved)
    }
    if (row == column) {
      moved_up <- point
      moved_up[row] <- moved_up[row] + step
      moved_down <- point
      moved_down[row] <- moved_down[row] - step
      (objective(moved_up) - 2 * base + objective(moved_down)) / step^2
    } else {
      (shift(c(1, 1)) - shift(c(1, -1)) - shift(c(-1, 1)) + shift(c(-1, -1))) /
        (4 * step^2)
    }
  }, numeric(1))
  hessian <- matrix(values, length(point), length(point))
  lower <- lower.tri(hessian)
  hessian[upper.tri(hessian)] <- t(hessian)[upper.tri(hessian)]
  hessian
}

#' Central-difference gradient of a scalar function
#' @param objective Function of a numeric vector returning one number.
#' @param point Where to differentiate.
#' @param step Positive finite-difference step.
#' @return A numeric gradient the same length as `point`.
central_difference_gradient <- function(objective, point, step = 1e-5) {
  vapply(seq_along(point), function(index) {
    moved_up <- point
    moved_up[index] <- moved_up[index] + step
    moved_down <- point
    moved_down[index] <- moved_down[index] - step
    (objective(moved_up) - objective(moved_down)) / (2 * step)
  }, numeric(1))
}

test_that("a held-measurement fit reports exactly the parameters it estimated", {
  prepared <- fixed_measurement_data()
  stage <- multilpa(prepared$data, c("a", "b"), "g", 2, 1, n_starts = 4,
                    seed = 11, tol = 1e-12)
  held <- multilpa(prepared$data, c("a", "b"), "g", 2, 2, n_starts = 3,
                   start = starting_values(stage, what = "measurement"),
                   fixed = "measurement", seed = 11, tol = 1e-12)
  expect_true(held$converged)
  expect_equal(held$fixed, c("means", "variances"))
  # The held blocks are the stage's own values, unmoved by the second stage.
  expect_equal(held$means, stage$means)
  expect_equal(held$variances, stage$variances)

  information <- parameter_inference(held, prepared$data)
  # Nothing measurement-level is reported, because nothing measurement-level was
  # estimated here.
  expect_false(any(information$level == "measurement"))
  expect_equal(nrow(information), held$n_group_classes * held$n_profiles +
                 held$n_group_classes)
  expect_equal(attr(information, "fixed"), c("means", "variances"))

  # The estimation-scale covariance is exactly n_parameters square: a held block
  # contributes no row and no column to the information matrix.
  estimation <- vcov(held, data = prepared$data, scale = "unconstrained")
  expect_equal(dim(estimation), rep(held$n_parameters, 2L))
  expect_equal(length(.multilpa_free_index(held, "unconstrained")),
               held$n_parameters)
  expect_equal(nrow(vcov(held, data = prepared$data)), nrow(information))
  expect_equal(nrow(confint(held, data = prepared$data)), nrow(information))
  expect_true(all(is.finite(information$standard_error)))
  expect_true(all(information$standard_error > 0))
})

test_that("held-measurement standard errors match an independent numerical information", {
  skip_on_cran()
  prepared <- fixed_measurement_data()
  stage <- multilpa(prepared$data, c("a", "b"), "g", 2, 1, n_starts = 4,
                    seed = 11, tol = 1e-12)
  held <- multilpa(prepared$data, c("a", "b"), "g", 2, 2, n_starts = 3,
                   start = starting_values(stage, what = "measurement"),
                   fixed = "measurement", seed = 11, tol = 1e-12)
  information <- parameter_inference(held, prepared$data)
  free <- .multilpa_free_index(held, "unconstrained")
  free_natural <- .multilpa_free_index(held, "natural")
  fitted_free <- unname(coef(held, scale = "unconstrained"))[free]
  held_log_variances <- as.vector(t(log(held$variances)))

  loglik <- function(values) {
    enumerated_conditional_loglik(held$means, held_log_variances, values[1:2],
                                  values[[3L]], prepared$x, prepared$group)
  }
  # The held blocks really are in the likelihood, at the values they were held
  # at: an independent enumeration of the same model reproduces the fit.
  expect_equal(loglik(fitted_free), held$log_likelihood, tolerance = 1e-10)

  # Analytic conditional score against numerical differentiation of the
  # independent conditional log likelihood, away from the optimum.
  centers <- colMeans(prepared$x)
  centred <- held
  centred$means <- sweep(held$means, 2L, centers, "-")
  restore <- .multilpa_restore_held(coef(centred, scale = "unconstrained"), free)
  centred_x <- sweep(prepared$x, 2L, centers, "-")
  point <- fitted_free + c(0.13, -0.09, 0.21)
  analytic_score <- -.multilpa_score(restore(point), centred_x, held)[free]
  numerical_score <- central_difference_gradient(loglik, point)
  expect_equal(analytic_score, numerical_score, tolerance = 1e-7)
  expect_lt(max(abs(analytic_score - numerical_score)), 1e-7)

  # Analytic observed information against second central differences of the
  # independent negative conditional log likelihood. Achieved agreement in the
  # session that wrote this test: 1.49e-07 relative.
  numerical_information <- central_difference_hessian(
    function(values) -loglik(values), fitted_free)
  analytic_information <- unname(attr(information, "hessian"))
  expect_equal(numerical_information, analytic_information, tolerance = 1e-6)

  # The delta-method map restricted to the free coordinates, against a
  # numerical Jacobian of the natural coefficients.
  natural_from_free <- function(values) {
    model <- held
    decoded <- .multilpa_decode(restore(values), held)
    model[names(decoded)] <- decoded
    unname(coef(model))[free_natural]
  }
  numerical_jacobian <- vapply(seq_along(fitted_free), function(index) {
    moved_up <- fitted_free
    moved_up[index] <- moved_up[index] + 1e-6
    moved_down <- fitted_free
    moved_down[index] <- moved_down[index] - 1e-6
    (natural_from_free(moved_up) - natural_from_free(moved_down)) / 2e-6
  }, numeric(length(free_natural)))
  expect_equal(numerical_jacobian,
               unname(.multilpa_inference_jacobian(held)[free_natural, free, drop = FALSE]),
               tolerance = 1e-8)

  # Standard errors assembled entirely from the independent pieces. Achieved
  # agreement in the session that wrote this test: 1.28e-07 relative.
  independent_errors <- sqrt(pmax(diag(
    numerical_jacobian %*% solve(numerical_information) %*% t(numerical_jacobian)), 0))
  expect_equal(information$standard_error, independent_errors, tolerance = 1e-6)
  expect_lt(max(abs(information$standard_error - independent_errors) /
                  information$standard_error), 1e-6)
})

test_that("holding only the means leaves the variances estimated and correct", {
  skip_on_cran()
  prepared <- fixed_measurement_data()
  stage <- multilpa(prepared$data, c("a", "b"), "g", 2, 1, n_starts = 2,
                    seed = 11, tol = 1e-12)
  held <- multilpa(prepared$data, c("a", "b"), "g", 2, 2, n_starts = 2,
                   start = starting_values(stage, what = "measurement"),
                   fixed = "means", seed = 11, tol = 1e-12)
  expect_true(held$converged)
  expect_equal(held$fixed, "means")
  expect_equal(held$means, stage$means)

  information <- parameter_inference(held, prepared$data)
  # The variances are estimated here and are reported; the means are not.
  expect_setequal(subset(information, level == "measurement")$parameter, "variance")
  expect_equal(length(.multilpa_free_index(held, "unconstrained")),
               held$n_parameters)

  free <- .multilpa_free_index(held, "unconstrained")
  fitted_free <- unname(coef(held, scale = "unconstrained"))[free]
  loglik <- function(values) {
    enumerated_conditional_loglik(held$means, values[1:4], values[5:6],
                                  values[[7L]], prepared$x, prepared$group)
  }
  expect_equal(loglik(fitted_free), held$log_likelihood, tolerance = 1e-10)

  numerical_information <- central_difference_hessian(
    function(values) -loglik(values), fitted_free)
  expect_equal(numerical_information, unname(attr(information, "hessian")),
               tolerance = 1e-5)
  expect_true(all(is.finite(information$standard_error)))
})

test_that("an unheld fit takes the unconstrained path unchanged", {
  set.seed(4100)
  dat <- data.frame(g = rep(seq_len(30), each = 8), y = rnorm(240, 3, 2))
  plain <- multilpa(dat, "y", "g", 1, 1, n_starts = 1, seed = 4)
  expect_equal(plain$fixed, character(0))
  # Every coordinate is free, on both scales, and the free map is the identity.
  expect_equal(.multilpa_free_index(plain, "unconstrained"),
               seq_along(coef(plain, scale = "unconstrained")))
  expect_equal(.multilpa_free_index(plain, "natural"),
               seq_along(coef(plain)))

  information <- parameter_inference(plain, dat)
  expect_equal(nrow(information), length(coef(plain)))
  # Closed-form maximum-likelihood errors for a single Gaussian: the fixed-block
  # machinery must not have moved the ordinary answer by so much as an ulp.
  # The single profile and single group class are certain, hence the two trailing
  # probabilities of one with no error at all.
  variance <- mean((dat$y - mean(dat$y))^2)
  expect_equal(information$estimate, c(mean(dat$y), variance, 1, 1))
  expect_equal(information$standard_error,
               c(sqrt(variance / nrow(dat)), variance * sqrt(2 / nrow(dat)), 0, 0),
               tolerance = 1e-7)
  expect_equal(unname(confint(plain, data = dat)),
               unname(cbind(information$conf_low, information$conf_high)))
  expect_equal(rownames(vcov(plain, data = dat)), names(coef(plain)))

  # A fit made with an explicit empty `fixed` is the same fit, not a different
  # code path: identical estimates, identical errors, identical intervals.
  explicit <- multilpa(dat, "y", "g", 1, 1, n_starts = 1, seed = 4,
                       fixed = character())
  expect_equal(parameter_inference(explicit, dat), information)
})

test_that("inference refuses what a held fit cannot answer", {
  prepared <- fixed_measurement_data()
  stage <- multilpa(prepared$data, c("a", "b"), "g", 2, 1, n_starts = 4,
                    seed = 11, tol = 1e-12)
  held <- multilpa(prepared$data, c("a", "b"), "g", 2, 2, n_starts = 3,
                   start = starting_values(stage, what = "measurement"),
                   fixed = "measurement", seed = 11, tol = 1e-12)

  # A held coefficient is in coef(), because it is part of the model, but it was
  # not estimated here and has no interval.
  expect_error(confint(held, parm = "measurement.mean.profile_1.a",
                       data = prepared$data),
               class = "latents_held_parameter")

  # Holding every parameter the model has leaves nothing to report.
  set.seed(4055)
  small <- data.frame(g = rep(seq_len(12), each = 5), a = rnorm(60), b = rnorm(60))
  single <- multilpa(small, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 2)
  everything <- multilpa(small, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 2,
                         start = starting_values(single, what = "measurement"),
                         fixed = "measurement")
  expect_equal(everything$n_parameters, 0L)
  expect_error(parameter_inference(everything, small),
               class = "latents_no_free_parameters")
  expect_error(vcov(everything, data = small),
               class = "latents_no_free_parameters")

  # A `fixed` field naming something that is not a measurement block is refused
  # rather than silently treated as free.
  corrupted <- held
  corrupted$fixed <- c("means", "profile_probabilities")
  expect_error(.multilpa_free_index(corrupted, "unconstrained"),
               class = "latents_bad_fixed")
  expect_error(parameter_inference(corrupted, prepared$data),
               class = "latents_bad_fixed")
})

test_that("fit_staged makes good on its promise of conditional standard errors", {
  skip_on_cran()
  set.seed(4200)
  dat <- data.frame(g = rep(seq_len(40), each = 8))
  class_of_group <- rep(c(1, 2), length.out = 40)
  profile <- rbinom(nrow(dat), 1, c(0.2, 0.8)[class_of_group][dat$g])
  dat$a <- rnorm(nrow(dat), 4 * profile, 1)
  dat$b <- rnorm(nrow(dat), 3 * profile, 1)
  staged <- fit_staged(dat, c("a", "b"), "g", n_profiles = 2,
                       n_group_classes = 2, n_starts = 3, seed = 6)
  information <- parameter_inference(staged, dat)
  expect_false(any(information$level == "measurement"))
  expect_true(all(is.finite(information$standard_error)))
  expect_equal(dim(vcov(staged, data = dat, scale = "unconstrained")),
               rep(staged$n_parameters, 2L))
  expect_equal(nrow(confint(staged, data = dat)), nrow(information))
  expect_true(staged$n_parameters_with_measurement > staged$n_parameters)
})
