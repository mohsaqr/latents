# Step-one uncertainty is checked against a separate Gaussian likelihood and
# step-three estimators; none of the reference derivatives use package internals.
skip_on_cran()

.step_one_case <- local({
  made <- NULL
  function() {
    if (!is.null(made)) return(made)
    set.seed(532)
    n <- 700L
    truth <- sample.int(2L, n, replace = TRUE, prob = c(0.57, 0.43))
    data <- data.frame(
      a = stats::rnorm(n, c(-1.1, 1.3)[truth], c(0.8, 1.1)[truth]),
      b = stats::rnorm(n, c(-0.8, 0.9)[truth], c(1.1, 0.85)[truth]),
      y = stats::rnorm(n, c(-0.5, 2)[truth], c(0.9, 1.2)[truth]),
      category = factor(ifelse(stats::runif(n) < c(0.28, 0.70)[truth], "yes", "no")))
    fit <- lpa(data, c("a", "b"), n_profiles = 2, model = "VVI", n_starts = 3,
               seed = 1, tol = 1e-12, max_iter = 10000)
    made <<- list(data = data, fit = fit)
    made
  }
})

# Independent coordinates: four means, four log standard deviations, one
# profile log odds per group class, then the group-class log odds if needed.
# This likelihood also integrates both latent levels for the group check.
.step_one_reference <- function(case, level = "individuals") {
  observed <- as.matrix(subset(case$data, select = c(a, b)))
  fit <- case$fit
  n_types <- fit$n_group_classes
  point <- c(as.vector(t(fit$means)), as.vector(t(log(sqrt(fit$variances)))),
             stats::qlogis(as.vector(subset(fit$profile_probabilities, select = 1L))),
             if (n_types > 1L) stats::qlogis(fit$group_probabilities[1L]))
  informative <- rowSums(!is.na(observed)) > 0L
  group_kept <- as.vector(rowsum(as.numeric(informative), fit$group_index,
                                 reorder = FALSE)) > 0L
  rows <- if (identical(level, "groups")) {
    match(seq_len(fit$n_groups), fit$group_index)[group_kept]
  } else which(informative)
  state <- function(theta) {
    means <- matrix(theta[seq_len(4L)], 2L, byrow = TRUE)
    sds <- matrix(exp(theta[4L + seq_len(4L)]), 2L, byrow = TRUE)
    prior <- stats::plogis(theta[8L + seq_len(n_types)])
    eta <- if (n_types == 1L) 1 else
      c(stats::plogis(theta[11L]), stats::plogis(-theta[11L]))
    log_density <- vapply(seq_len(2L), function(class) {
      density_a <- stats::dnorm(as.vector(subset(observed, select = a)), means[class, 1L],
                                sds[class, 1L], log = TRUE)
      density_b <- stats::dnorm(as.vector(subset(observed, select = b)), means[class, 2L],
                                sds[class, 2L], log = TRUE)
      replace(density_a, is.na(density_a), 0) + replace(density_b, is.na(density_b), 0)
    }, numeric(nrow(observed)))
    largest <- apply(log_density, 1L, max)
    density <- exp(log_density - largest)
    conditional <- lapply(seq_len(n_types), function(type) {
      joint <- sweep(density, 2L, c(prior[type], 1 - prior[type]), "*")
      marginal <- rowSums(joint)
      list(posterior = joint / marginal,
           evidence = as.vector(rowsum(largest + log(marginal), fit$group_index,
                                        reorder = FALSE)))
    })
    log_group <- sweep(vapply(conditional, `[[`, numeric(fit$n_groups), "evidence"),
                        2L, log(eta), "+")
    group_max <- apply(log_group, 1L, max)
    group_joint <- exp(log_group - group_max)
    group_posterior <- group_joint / rowSums(group_joint)
    subject_posterior <- Reduce(`+`, lapply(seq_len(n_types), function(type) {
      conditional[[type]]$posterior * as.vector(subset(group_posterior, select = type))[
        fit$group_index]
    }))
    list(log_likelihood = sum(group_max + log(rowSums(group_joint))),
         posteriors = if (identical(level, "groups")) {
           group_posterior[group_kept, , drop = FALSE]
         } else subject_posterior[informative, , drop = FALSE])
  }
  log_likelihood <- function(theta) state(theta)$log_likelihood
  list(point = point, covariance = solve(-numDeriv::hessian(log_likelihood, point)),
       posteriors = function(theta) state(theta)$posteriors, rows = rows)
}

# A separate EM implementation with a parameter-based stopping rule. At
# perturbed theta1 the assigned classes remain those of the original fit.
.step_one_estimator <- function(case, reference, method, outcome, fixed, variance = "varying") {
  modal <- max.col(reference$posteriors(reference$point), ties.method = "first")
  categorical <- identical(outcome, "category")
  values <- if (categorical) as.numeric(case$data$category[reference$rows] == "no") else
    case$data$y[reference$rows]
  center <- if (categorical) 0 else mean(values)
  scale <- if (categorical) 1 else stats::sd(values)
  values <- (values - center) / scale
  start <- if (categorical) subset(fixed, category == "no")$estimate else
    (fixed$estimate - center) / scale
  start_variance <- if (categorical) NULL else fixed$variance / scale^2
  function(theta) {
    posterior <- reference$posteriors(theta)
    errors <- vapply(seq_len(2L), function(assigned) {
      colSums(posterior * (modal == assigned)) / colSums(posterior)
    }, numeric(2L))
    if (!identical(method, "ml")) {
      weights <- switch(method,
                        bch = solve(errors)[modal, , drop = FALSE],
                        proportional = posterior,
                        modal = outer(modal, seq_len(2L), "==") * 1)
      estimate <- as.vector(crossprod(weights, values)) / colSums(weights)
    } else {
      prior <- colMeans(posterior)
      fixed_weights <- t(subset(errors, select = modal)) *
        matrix(prior, length(values), 2L, byrow = TRUE)
      estimate <- start
      variances <- start_variance
      iteration <- 0L
      change <- Inf
      while (iteration < 10000L && change > 1e-13) {
        iteration <- iteration + 1L
        likelihood <- if (categorical) {
          outer(values, estimate, function(y, p) ifelse(y == 1, p, 1 - p))
        } else vapply(seq_len(2L), function(class) {
          stats::dnorm(values, estimate[class], sqrt(variances[class]))
        }, numeric(length(values)))
        joint <- fixed_weights * likelihood
        responsibility <- joint / rowSums(joint)
        updated <- as.vector(crossprod(responsibility, values)) / colSums(responsibility)
        updated_variance <- if (categorical) NULL else {
          squares <- responsibility * outer(values, updated, "-")^2
          if (identical(variance, "equal")) rep(sum(squares) / length(values), 2L) else
            colSums(squares) / colSums(responsibility)
        }
        change <- max(abs(c(updated - estimate, updated_variance - variances)))
        estimate <- updated
        variances <- updated_variance
      }
      stopifnot("the independent step-three EM must converge" = change <= 1e-13)
    }
    estimate <- center + scale * estimate
    if (categorical) as.vector(rbind(estimate, 1 - estimate)) else estimate
  }
}

.step_one_gradient <- function(estimator, point) {
  size <- length(estimator(point))
  vapply(seq_along(point), function(j) {
    h <- 1e-4 * max(1, abs(point[j]))
    difference <- function(step) {
      up <- point
      down <- point
      up[j] <- up[j] + step
      down[j] <- down[j] - step
      (estimator(up) - estimator(down)) / (2 * step)
    }
    (4 * difference(h / 2) - difference(h)) / 3
  }, numeric(size))
}

# With two classes, their marginal SEs and their difference SE identify
# their covariance. For a binary outcome p(no) determines p(yes) exactly.
.step_one_fixed_covariance <- function(fixed, pairs, categorical) {
  errors <- if (categorical) subset(fixed, category == "no")$standard_error else
    fixed$standard_error
  difference <- if (categorical) subset(pairs, category == "no")$standard_error else
    pairs$standard_error
  covariance <- diag(errors^2)
  covariance[1L, 2L] <- covariance[2L, 1L] <-
    (sum(errors^2) - difference^2) / 2
  if (!categorical) return(covariance)
  transform <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  transform %*% covariance %*% t(transform)
}

test_that("both levels use the kept posteriors and include their step-one uncertainty", {
  skip_if_not_installed("numDeriv")
  set.seed(613)
  n_groups <- 90L
  per_group <- 7L
  group_class <- rep(seq_len(2L), length.out = n_groups)
  truth <- 1L + as.integer(stats::runif(n_groups * per_group) <
                             c(0.15, 0.85)[rep(group_class, each = per_group)])
  data <- data.frame(
    g = rep(seq_len(n_groups), each = per_group),
    a = stats::rnorm(n_groups * per_group, c(-1.6, 1.6)[truth]),
    b = stats::rnorm(n_groups * per_group, c(-1.3, 1.3)[truth]),
    y = rep(stats::rnorm(n_groups, c(-0.6, 1.1)[group_class], 1.4), each = per_group),
    category = factor(rep(ifelse(stats::runif(n_groups) < c(0.3, 0.7)[group_class],
                                  "yes", "no"), each = per_group)))
  # One whole group and three other rows are unclassified and must be
  # dropped from D and step three, while remaining in the fitted likelihood.
  data$a[seq_len(10L)] <- data$b[seq_len(10L)] <- NA_real_
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                  n_starts = 3, seed = 1, tol = 1e-11, max_iter = 10000, missing = "fiml")
  case <- list(data = data, fit = fit)
  invisible(lapply(c("individuals", "groups"), function(level) {
    reference <- .step_one_reference(case, level)
    grid <- expand.grid(method = c("bch", "proportional", "ml"),
                        outcome = c("y", "category"), stringsAsFactors = FALSE)
    invisible(lapply(seq_len(nrow(grid)), function(row) {
      method <- grid$method[row]
      outcome <- grid$outcome[row]
      run <- function(step_one, contrast = "none") {
        three_step(fit, data, outcome, level = level, method = method,
                   step_one = step_one, contrast = contrast)
      }
      fixed <- run("fixed")
      corrected <- run("estimated")
      estimator <- .step_one_estimator(case, reference, method, outcome, fixed)
      expect_equal(estimator(reference$point), fixed$estimate, tolerance = 1e-8)
      gradient <- .step_one_gradient(estimator, reference$point)
      covariance <- .step_one_fixed_covariance(fixed, run("fixed", "pairs"),
                                               identical(outcome, "category"))
      expected <- sqrt(diag(covariance + gradient %*% reference$covariance %*% t(gradient)))
      expect_equal(corrected$standard_error / expected, rep(1, length(expected)),
                   tolerance = 1e-5)
      expect_true(all(corrected$standard_error >= fixed$standard_error))
    }))
  }))
})

test_that("fixed step one is the default for every result path", {
  case <- .step_one_case()
  grid <- expand.grid(method = c("bch", "proportional", "modal", "ml"),
                      outcome = c("y", "category"),
                      contrast = c("none", "pairs", "overall"),
                      stringsAsFactors = FALSE)
  invisible(lapply(seq_len(nrow(grid)), function(row) {
    args <- list(x = case$fit, data = case$data, outcome = grid$outcome[row],
                 method = grid$method[row], contrast = grid$contrast[row])
    default <- do.call(three_step, args)
    fixed <- do.call(three_step, c(args, list(step_one = "fixed")))
    expect_identical(fixed, default)
    expect_identical(attr(fixed, "step_one"), "fixed")
  }))
})

test_that("estimated step one agrees with an independent likelihood and re-solve", {
  skip_if_not_installed("numDeriv")
  case <- .step_one_case()
  reference <- .step_one_reference(case)
  expect_equal(unname(reference$posteriors(reference$point)),
               unname(case$fit$subject_posteriors), tolerance = 1e-10)
  grid <- expand.grid(method = c("bch", "proportional", "ml"),
                      outcome = c("y", "category"), stringsAsFactors = FALSE)
  invisible(lapply(seq_len(nrow(grid)), function(row) {
    method <- grid$method[row]
    outcome <- grid$outcome[row]
    categorical <- identical(outcome, "category")
    run <- function(step_one, contrast = "none") {
      three_step(case$fit, case$data, outcome, method = method,
                 step_one = step_one, contrast = contrast)
    }
    fixed <- run("fixed")
    pairs <- run("fixed", "pairs")
    corrected <- run("estimated")
    estimator <- .step_one_estimator(case, reference, method, outcome, fixed)
    expect_equal(estimator(reference$point), fixed$estimate, tolerance = 1e-8)
    gradient <- .step_one_gradient(estimator, reference$point)
    fixed_covariance <- .step_one_fixed_covariance(fixed, pairs, categorical)
    covariance <- fixed_covariance + gradient %*% reference$covariance %*% t(gradient)
    expected <- sqrt(diag(covariance))
    expect_equal(corrected$standard_error / expected, rep(1, length(expected)),
                 tolerance = 1e-5)
    expect_identical(attr(corrected, "step_one"), "estimated")
    expect_identical(corrected$estimate, fixed$estimate)
    expect_true(all(corrected$standard_error >= fixed$standard_error))
    pair_contrast <- if (categorical) rbind(c(-1, 0, 1, 0), c(0, -1, 0, 1)) else
      matrix(c(-1, 1), 1L)
    pair_covariance <- pair_contrast %*% covariance %*% t(pair_contrast)
    expected_pairs <- sqrt(diag(pair_covariance))
    expect_equal(run("estimated", "pairs")$standard_error / expected_pairs,
                 rep(1, length(expected_pairs)), tolerance = 1e-5)
    estimates <- fixed$estimate
    overall_contrast <- if (categorical) {
      matrix(c(-1 / estimates[1L], 1 / estimates[2L],
               1 / estimates[3L], -1 / estimates[4L]), 1L)
    } else pair_contrast
    difference <- if (categorical) {
      log(estimates[3L] / estimates[4L]) - log(estimates[1L] / estimates[2L])
    } else diff(estimates)
    expected_overall <- difference^2 /
      as.numeric(overall_contrast %*% covariance %*% t(overall_contrast))
    expect_equal(run("estimated", "overall")$statistic / expected_overall, 1,
                 tolerance = 1e-5)
  }))
})

test_that("modal estimates and uncertainty are unchanged when step one is estimated", {
  case <- .step_one_case()
  invisible(lapply(c("y", "category"), function(outcome) {
    invisible(lapply(c("none", "pairs", "overall"), function(contrast) {
      fixed <- three_step(case$fit, case$data, outcome, method = "modal", contrast = contrast)
      estimated <- three_step(case$fit, case$data, outcome, method = "modal", contrast = contrast,
                              step_one = "estimated")
      expect_identical(attr(estimated, "step_one"), "estimated")
      attr(estimated, "step_one") <- "fixed"
      expect_identical(estimated, fixed)
    }))
  }))
})

test_that("estimated step one needs only the outcome in data", {
  case <- .step_one_case()
  full <- three_step(case$fit, case$data, "y", step_one = "estimated")
  # The step-one inference reads the fit's own indicators, as vcov() does.
  expect_warning(
    outcome_only <- three_step(case$fit, subset(case$data, select = "y"), "y",
                               step_one = "estimated"),
    class = "latents_unverified_alignment")
  expect_identical(outcome_only, full)
})

test_that("a step-one covariance of deficient rank is used, not refused", {
  # Person centring makes the profile-weighted means sum to zero within every
  # group, so the step-one covariance has zero directions.
  set.seed(31)
  n <- 240L
  truth <- sample(1:2, n, TRUE, c(0.6, 0.4))
  shift <- rep(stats::rnorm(40L), each = 6L)
  data <- data.frame(g = rep(seq_len(40L), each = 6L),
                     a = 5 + shift + stats::rnorm(n, c(-1.2, 1.8)[truth], c(1, 0.8)[truth]),
                     b = -3 + shift + stats::rnorm(n, c(1, -1.5)[truth], c(0.9, 1.1)[truth]),
                     y = stats::rnorm(n))
  fit <- multilpa(data, c("a", "b"), "g", 2, 1, centering = "person", n_starts = 2,
                  seed = 1, tol = 1e-10)
  values <- eigen(vcov(fit, scale = "unconstrained"), symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(abs(values)), 1e-12)
  fixed <- three_step(fit, data, "y", method = "proportional")
  estimated <- three_step(fit, data, "y", method = "proportional", step_one = "estimated")
  expect_true(all(is.finite(estimated$standard_error)))
  expect_true(all(estimated$standard_error >= fixed$standard_error))
})

test_that("the correction does not depend on the units of the data", {
  set.seed(5)
  n <- 500L
  truth <- sample(1:2, n, TRUE)
  data <- data.frame(a = stats::rnorm(n, c(0, 2)[truth]), b = stats::rnorm(n, c(0, 1.5)[truth]),
                     y = stats::rnorm(n, c(0, 1)[truth]))
  data$b <- data$b + 0.5 * data$a
  rescaled <- transform(data, b = b * 1e4)
  full <- function(frame) {
    lpa(frame, c("a", "b"), 2, model = "VVV", n_starts = 3, seed = 1, tol = 1e-12,
        max_iter = 10000)
  }
  original_fit <- full(data)
  rescaled_fit <- full(rescaled)
  invisible(lapply(c("bch", "proportional", "ml"), function(method) {
    expect_equal(three_step(rescaled_fit, rescaled, "y", method = method,
                            step_one = "estimated")$standard_error,
                 three_step(original_fit, data, "y", method = method,
                            step_one = "estimated")$standard_error,
                 tolerance = 1e-6, label = paste("full covariance,", method))
  }))
  # A direct effect is in indicator units per covariate unit.
  data$z <- stats::rnorm(n, 0.5, 1.5)
  data$a <- data$a + 0.8 * data$z
  rescaled <- transform(data, z = z * 1e6)
  direct <- function(frame) {
    lpa(frame, c("a", "b"), 2, n_starts = 3, seed = 1, tol = 1e-12, max_iter = 10000,
        direct_effects = "a ~ z")
  }
  expect_equal(three_step(direct(rescaled), rescaled, "y", method = "proportional",
                          step_one = "estimated")$standard_error,
               three_step(direct(data), data, "y", method = "proportional",
                          step_one = "estimated")$standard_error,
               tolerance = 1e-5, label = "direct effect")
})

test_that("estimated step one refuses unsupported fits and preserves inference errors", {
  case <- .step_one_case()
  expect_error(three_step(case$fit, case$data, "y", step_one = "unknown"),
               class = "latents_bad_argument")
  boundary <- case$fit
  boundary$boundary <- TRUE
  expect_error(three_step(boundary, case$data, "y", step_one = "estimated"),
               'step_one = "fixed"', class = "latents_boundary_fit", fixed = TRUE)
  expect_no_error(three_step(boundary, case$data, "y", step_one = "fixed"))
  repeated <- transform(subset(case$data, select = c(a, b, y)),
                        id = rep(seq_len(140L), each = 5L),
                        time = rep(seq_len(5L), 140L))
  transition <- lta(repeated, c("a", "b"), "id", time = "time", n_profiles = 2,
                    n_starts = 1, seed = 1)
  expect_error(three_step(transition, repeated, "y", step_one = "estimated"),
               'step_one = "fixed"', class = "latents_unsupported_three_step", fixed = TRUE)
  trajectory <- mixture_regression(y ~ time, repeated, n_classes = 2, id = "id",
                                    class_level = "group", n_starts = 1, seed = 1,
                                    vcov_type = "none")
  expect_error(three_step(trajectory, repeated, "y", step_one = "estimated"),
               'step_one = "fixed"', class = "latents_unsupported_three_step", fixed = TRUE)
})

test_that("equal-variance ML includes step-one uncertainty", {
  skip_if_not_installed("numDeriv")
  case <- .step_one_case()
  reference <- .step_one_reference(case)
  fixed <- three_step(case$fit, case$data, "y", method = "ml", variance = "equal")
  pairs <- three_step(case$fit, case$data, "y", method = "ml", variance = "equal",
                      contrast = "pairs")
  corrected <- three_step(case$fit, case$data, "y", method = "ml", variance = "equal",
                          step_one = "estimated")
  estimator <- .step_one_estimator(case, reference, "ml", "y", fixed, "equal")
  gradient <- .step_one_gradient(estimator, reference$point)
  covariance <- .step_one_fixed_covariance(fixed, pairs, FALSE) +
    gradient %*% reference$covariance %*% t(gradient)
  expect_equal(corrected$standard_error / sqrt(diag(covariance)), c(1, 1), tolerance = 1e-5)
})

test_that("interior LCA fits support estimated step one for both outcome types", {
  set.seed(917)
  n <- 600L
  truth <- sample.int(2L, n, replace = TRUE)
  data <- as.data.frame(stats::setNames(lapply(seq_len(5L), function(item) {
    factor(stats::rbinom(n, 1L, c(0.18, 0.78)[truth]), levels = c(0, 1))
  }), paste0("u", seq_len(5L))))
  data$y <- stats::rnorm(n, c(-0.4, 1.3)[truth], 1.2)
  data$category <- factor(ifelse(stats::runif(n) < c(0.35, 0.65)[truth], "yes", "no"))
  fit <- lca(data, paste0("u", seq_len(5L)), n_classes = 2, n_starts = 3, seed = 1)
  invisible(lapply(c("y", "category"), function(outcome) {
    invisible(lapply(c("bch", "proportional", "ml"), function(method) {
      fixed <- three_step(fit, data, outcome, method = method)
      estimated <- three_step(fit, data, outcome, method = method, step_one = "estimated")
      expect_identical(estimated$estimate, fixed$estimate)
      expect_true(all(is.finite(estimated$standard_error)))
      expect_true(all(estimated$standard_error >= fixed$standard_error))
    }))
  }))
})
