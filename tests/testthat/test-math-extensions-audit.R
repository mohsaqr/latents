test_that("univariate random-intercept mixtures retain profile-by-indicator shapes", {
  set.seed(4)
  data <- data.frame(g = rep(seq_len(12L), each = 6L),
                     y = rnorm(72L, rep(c(-2, 2), 36L), 0.5))
  fits <- lapply(c("varying", "equal"), function(variance_model)
    fit_random_intercept(data, "y", "g", 2L, n_starts = 1L,
                         variance_model = variance_model, seed = 1L))
  invisible(lapply(fits, function(fit) {
    expect_identical(dim(fit$means), c(2L, 1L))
    expect_identical(dim(fit$variances), c(2L, 1L))
    expect_equal(rowSums(fit$subject_posteriors), rep(1, 72L), tolerance = 1e-12)
    expect_true(is.finite(fit$log_likelihood))
    expect_true(fit$converged)
    expect_equal(unname(sort(drop(fit$means))), c(-2, 2), tolerance = 0.2)
  }))
})

test_that("mixture random-intercept posteriors agree with adaptive integration", {
  x <- matrix(c(-0.5, 1, 0.2, 0.7, 0.9, -0.3), 3L, 2L)
  parameters <- list(means = matrix(c(-0.6, 0.8, -0.2, 0.6), 2L),
    variances = matrix(c(0.5, 0.8, 0.6, 0.9), 2L), random_sd = 0.35,
    profile_probabilities = c(0.3, 0.7))
  actual <- .ri_evaluate(x, c(1L, 1L, 2L), parameters, .ri_quadrature(61L), TRUE)
  reference <- lapply(list(1:2, 3L), function(rows) {
    density <- function(intercept) {
      vapply(seq_len(2L), function(k) parameters$profile_probabilities[k] *
        dnorm(x[rows, 1L], parameters$means[k, 1L] + intercept,
              sqrt(parameters$variances[k, 1L])) *
        dnorm(x[rows, 2L], parameters$means[k, 2L] + intercept,
              sqrt(parameters$variances[k, 2L])), numeric(length(rows))) |>
        matrix(nrow = length(rows), ncol = 2L)
    }
    integrate_weighted <- function(weight) {
      integrate(function(intercepts) vapply(intercepts, function(intercept) {
        components <- density(intercept)
        likelihood <- prod(rowSums(components)) * dnorm(intercept, sd = parameters$random_sd)
        if (likelihood == 0) 0 else likelihood * weight(intercept, components)
      }, numeric(1)), -Inf, Inf, rel.tol = 1e-11)$value
    }
    normalizer <- integrate_weighted(function(intercept, components) 1)
    mean <- integrate_weighted(function(intercept, components) intercept) / normalizer
    second <- integrate_weighted(function(intercept, components) intercept^2) / normalizer
    posteriors <- vapply(seq_along(rows), function(row)
      integrate_weighted(function(intercept, components)
        components[row, 1L] / sum(components[row, ])) / normalizer, numeric(1))
    list(mean = mean, sd = sqrt(second - mean^2), posterior = posteriors)
  })
  expect_equal(unname(actual$random_intercept_mean),
               vapply(reference, `[[`, numeric(1), "mean"), tolerance = 1e-10)
  expect_equal(unname(actual$random_intercept_sd),
               vapply(reference, `[[`, numeric(1), "sd"), tolerance = 1e-10)
  expect_equal(actual$subject_posteriors[, 1L],
               unlist(lapply(reference, `[[`, "posterior")), tolerance = 1e-10)
})

test_that("pooled Gaussian residuals exclude association explained by profile means", {
  residual <- as.matrix(expand.grid(a = c(-1, 1), b = c(-1, 1)))
  data <- as.data.frame(rbind(residual - 5, residual + 5))
  object <- structure(list(n_observations = 8L, n_profiles = 2L,
    vars = c("a", "b"), means = matrix(c(-5, 5, -5, 5), 2L),
    variances = matrix(1, 2L, 2L), covariance_model = "diagonal",
    subject_posteriors = cbind(rep(c(1, 0), each = 4L), rep(c(0, 1), each = 4L))),
    class = "multilpa")
  result <- bivariate_residuals(object, data, "overall")
  expect_gt(cor(data$a, data$b), 0.9)
  expect_equal(result$observed, 0)
  expect_equal(result$expected, 0)
  expect_equal(result$p_value, 1)
  # Pool fitted within-profile covariances as well as empirical moments.
  data$b <- data$a
  object$covariance_model <- "full"
  object$covariances <- array(1, c(2L, 2L, 2L))
  resolved <- bivariate_residuals(object, data, "overall")
  expect_equal(resolved$observed, 1)
  expect_equal(resolved$expected, 1)
  expect_equal(resolved$residual, 0)
})

test_that("categorical profile residuals compare to their own product probabilities", {
  patterns <- expand.grid(a = c("no", "yes"), b = c("no", "yes"),
                          stringsAsFactors = FALSE)
  probabilities <- rbind(c(0.8, 0.2), c(0.2, 0.8))
  component <- vapply(seq_len(2L), function(k)
    as.vector(outer(probabilities[k, ], probabilities[k, ])), numeric(4L))
  marginal <- rowMeans(component)
  indices <- rep(seq_len(4L), round(1000 * marginal))
  data <- patterns[indices, , drop = FALSE]
  object <- structure(list(n_observations = 1000L, n_profiles = 2L,
    vars = c("a", "b"), continuous = character(), categorical = c("a", "b"),
    categorical_levels = list(a = c("no", "yes"), b = c("no", "yes")),
    response_probabilities = list(a = probabilities, b = probabilities),
    subject_posteriors = (component / rowSums(component))[indices, , drop = FALSE]),
    class = "multilpa")
  profile <- bivariate_residuals(object, data)
  overall <- bivariate_residuals(object, data, "overall")
  expect_equal(profile$statistic, c(0, 0), tolerance = 1e-20)
  expect_equal(profile$effective_n, c(500, 500))
  expect_equal(overall$statistic, 0, tolerance = 1e-20)
  # Preserve all modelled categories even when a pair has no observation in one.
  data$b[data$b == "yes"] <- NA_character_
  incomplete <- bivariate_residuals(object, data, "overall")
  # Complete-pair posteriors need not keep equal class shares.
  complete <- !is.na(data$b)
  counts <- colSums(object$subject_posteriors[complete, , drop = FALSE])
  expected <- counts[1L] * outer(probabilities[1L, ], probabilities[1L, ]) +
    counts[2L] * outer(probabilities[2L, ], probabilities[2L, ])
  observed <- matrix(c(340, 160, 0, 0), 2L)
  expect_equal(incomplete$statistic, sum((observed - expected)^2 / expected))
  expect_equal(incomplete$effective_n, 500)
  expect_equal(incomplete$df, 1)
  data$a[1L] <- "unknown"
  expect_error(bivariate_residuals(object, data), "levels absent")
})

test_that("Gaussian residuals use the complete rows of each pair", {
  data <- data.frame(a = c(-2, -1, 1, 2, NA), b = c(1, -1, -1, 1, 3))
  object <- structure(list(n_observations = 5L, n_profiles = 1L,
    vars = c("a", "b"), means = matrix(0, 1L, 2L),
    variances = matrix(1, 1L, 2L), covariance_model = "diagonal",
    subject_posteriors = matrix(1, 5L, 1L)), class = "multilpa")
  result <- bivariate_residuals(object, data)
  expect_equal(result$effective_n, 4)
  expect_equal(result$observed, 0)
  expect_equal(result$p_value, 1)
  data$a[1L] <- NA_real_
  expect_true(is.na(bivariate_residuals(object, data)$p_value))
  data$a[] <- NA_real_
  empty <- bivariate_residuals(object, data)
  expect_equal(empty$effective_n, 0)
  expect_true(is.na(empty$observed))
})

test_that("LMR rejects equally sized fits on different data or group layouts", {
  object <- structure(list(n_observations = 4L, n_groups = 2L,
    vars = "y", group = "g", group_values = 1:2,
    group_index = c(1L, 1L, 2L, 2L), indicator_data = matrix(1:4, 4L),
    n_parameters = 2L, log_likelihood = -10), class = "multilpa")
  alternative <- object
  alternative$n_parameters <- 3L
  alternative$log_likelihood <- -9
  expect_equal(lmr_lrt(object, alternative)$statistic, 2)
  alternative$indicator_data[1L, 1L] <- 9L
  expect_error(lmr_lrt(object, alternative), class = "multilpa_incomparable_models")
  alternative$indicator_data <- object$indicator_data
  alternative$group_index <- c(1L, 2L, 1L, 2L)
  expect_error(lmr_lrt(object, alternative), class = "multilpa_incomparable_models")
  alternative$group_index <- object$group_index
  invisible(lapply(c("covariance_model", "variance_model", "min_variance", "min_probability"),
    function(field) {
      changed <- alternative
      changed[[field]] <- "different"
      expect_error(lmr_lrt(object, changed), class = "multilpa_incomparable_models")
    }))
  object$n_profiles <- 1L
  object$n_group_classes <- 2L
  alternative$n_profiles <- 3L
  alternative$n_group_classes <- 1L
  expect_error(lmr_lrt(object, alternative), class = "multilpa_bad_nesting")
  object$n_groups <- alternative$n_groups <- 1L
  expect_error(lmr_lrt(object, alternative, n = "groups"), "greater than one")
})

test_that("three-step methods reject nonfinite outcomes and empty modal classes", {
  object <- structure(list(n_observations = 4L, n_profiles = 2L,
    subject_posteriors = matrix(rep(c(0.8, 0.2), each = 4L), 4L),
    subject_profiles = rep(1L, 4L), group_index = c(1L, 1L, 2L, 2L)),
    class = "multilpa")
  data <- data.frame(y = c(1, 2, 3, Inf))
  expect_error(three_step(object, data, "y"), "must be finite")
  expect_error(r3step(object, data, "y"), "must be finite")
  data$y[4L] <- 4
  expect_error(three_step(object, data, "y", method = "modal"),
               class = "multilpa_inseparable_classes")
  class(object) <- "multilpa_random_intercept"
  expect_error(three_step(object, data, "y"), "discrete group-class")
  expect_error(bivariate_residuals(object, data), "discrete group-class")
})
