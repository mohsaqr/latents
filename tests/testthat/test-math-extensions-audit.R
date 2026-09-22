test_that("pooled Gaussian residuals exclude association explained by profile means", {
  residual <- as.matrix(expand.grid(a = c(-1, 1), b = c(-1, 1)))
  data <- as.data.frame(rbind(residual - 5, residual + 5))
  object <- structure(list(n_observations = 8L, n_profiles = 2L,
    vars = c("a", "b"), means = matrix(c(-5, 5, -5, 5), 2L),
    variances = matrix(1, 2L, 2L), covariance_model = "diagonal",
    subject_posteriors = cbind(rep(c(1, 0), each = 4L), rep(c(0, 1), each = 4L))),
    class = "multilpa")
  result <- get_data(object, "residuals", data = data, by = "overall")
  expect_gt(cor(data$a, data$b), 0.9)
  expect_equal(result$observed, 0)
  expect_equal(result$expected, 0)
  expect_equal(result$p_value, 1)
  # Pool fitted within-profile covariances as well as empirical moments.
  data$b <- data$a
  object$covariance_model <- "full"
  object$covariances <- array(1, c(2L, 2L, 2L))
  resolved <- get_data(object, "residuals", data = data, by = "overall")
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
  profile <- get_data(object, "residuals", data = data)
  overall <- get_data(object, "residuals", data = data, by = "overall")
  expect_equal(profile$statistic, c(0, 0), tolerance = 1e-20)
  expect_equal(profile$effective_n, c(500, 500))
  expect_equal(overall$statistic, 0, tolerance = 1e-20)
  # Preserve all modelled categories even when a pair has no observation in one.
  data$b[data$b == "yes"] <- NA_character_
  incomplete <- get_data(object, "residuals", data = data, by = "overall")
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
  expect_error(get_data(object, "residuals", data = data), "levels absent")
})

test_that("Gaussian residuals use the complete rows of each pair", {
  data <- data.frame(a = c(-2, -1, 1, 2, NA), b = c(1, -1, -1, 1, 3))
  object <- structure(list(n_observations = 5L, n_profiles = 1L,
    vars = c("a", "b"), means = matrix(0, 1L, 2L),
    variances = matrix(1, 1L, 2L), covariance_model = "diagonal",
    subject_posteriors = matrix(1, 5L, 1L)), class = "multilpa")
  result <- get_data(object, "residuals", data = data)
  expect_equal(result$effective_n, 4)
  expect_equal(result$observed, 0)
  expect_equal(result$p_value, 1)
  data$a[1L] <- NA_real_
  expect_true(is.na(get_data(object, "residuals", data = data)$p_value))
  data$a[] <- NA_real_
  empty <- get_data(object, "residuals", data = data)
  expect_equal(empty$effective_n, 0)
  expect_true(is.na(empty$observed))
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
})
