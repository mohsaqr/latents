# Ordinal (adjacent-category logit) and count (Poisson) indicators.

extra_rows <- function(seed = 1, n = 600) {
  set.seed(seed)
  profile <- rbinom(n, 1, 0.4) + 1
  probabilities <- exp(.latents_ordinal_log_probabilities(c(0.5, 0, -1), c(1.2, 0)))
  levels <- c("never", "rarely", "often", "always")
  data.frame(
    y = rnorm(n, c(0, 1.5)[profile]),
    o = factor(vapply(profile, function(h) sample(levels, 1L, prob = probabilities[h, ]),
                      character(1)), levels = levels),
    k = rpois(n, c(1, 5)[profile]))
}

extra_fit <- function(data, ...) {
  quietly(lpa(data, c("y", "o", "k"), 2, ordinal = "o", count = "k", n_starts = 3,
              seed = 1, ...), "latents_single_level")
}

numeric_gradient <- function(objective, theta, step = 1e-6) {
  vapply(seq_along(theta), function(index) {
    up <- theta
    down <- theta
    up[index] <- up[index] + step
    down[index] <- down[index] - step
    (objective(up) - objective(down)) / (2 * step)
  }, numeric(1))
}

test_that("the ordinal density is the adjacent-category logit", {
  log_probability <- .latents_ordinal_log_probabilities(c(0.2, -0.5), c(0.7, 0))
  expect_equal(rowSums(exp(log_probability)), c(1, 1))
  # log P(k) / P(k - 1) = a_k - a_(k-1) + location.
  expect_equal(diff(log_probability[1L, ]), diff(c(0, 0.2, -0.5)) + 0.7)
  expect_equal(diff(log_probability[2L, ]), diff(c(0, 0.2, -0.5)))
  data <- data.frame(o = c(1, 3, NA, 2), k = c(0, 4, 2, NA))
  extra <- .latents_prepare_extra(data, "o", "k", "fiml")
  parameters <- list(ordinal_intercepts = list(c(0.2, -0.5)),
                     ordinal_locations = matrix(c(0.7, 0), 2L),
                     count_means = matrix(c(1.5, 4), 2L))
  density <- .latents_extra_log_density(extra, parameters, 4L, 2L)
  by_hand <- t(vapply(seq_len(4L), function(i) {
    ordinal <- if (is.na(data$o[i])) c(0, 0) else log_probability[, data$o[i]]
    count <- if (is.na(data$k[i])) c(0, 0) else dpois(data$k[i], c(1.5, 4), log = TRUE)
    ordinal + count
  }, numeric(2)))
  expect_equal(density, by_hand, tolerance = 1e-14)
})

test_that("the ordinal M-step maximizes its concave objective", {
  set.seed(2)
  counts <- matrix(runif(15, 1, 40), 3L, 5L)
  solved <- .latents_ordinal_maximize(counts, numeric(4), numeric(3))
  objective <- function(theta) {
    -sum(counts * .latents_ordinal_log_probabilities(theta[1:4], c(theta[5:6], 0)))
  }
  reference <- optim(numeric(6), objective, method = "BFGS",
                     control = list(reltol = 1e-15, maxit = 10000))
  ours <- objective(c(solved$intercepts, solved$locations[1:2]))
  expect_lte(ours, reference$value + 1e-10)
  expect_equal(c(solved$intercepts, solved$locations[1:2]), reference$par,
               tolerance = 1e-5)
  expect_identical(solved$locations[3L], 0)
})

test_that("fits agree with Latent GOLD 6.1", {
  skip_on_cran()
  fixture <- readRDS(test_path("fixtures", "latentgold-ordinal.rds"))
  vapply(names(fixture), function(name) {
    case <- fixture[[name]]
    fit <- quietly(multilpa(case$data, case$vars, if (case$two_level) "g" else NULL,
                            case$n_profiles,
                            n_group_classes = if (case$two_level) case$n_group_classes else 1L,
                            ordinal = case$vars[case$types == "ordinal"],
                            count = case$vars[case$types == "count"],
                            n_starts = 10, seed = 1, tol = 1e-12),
                   "latents_single_level")
    # Latent GOLD prints four decimals.
    expect_lt(abs(fit$log_likelihood - case$latent_gold_ll), 1e-4)
    expect_equal(fit$n_parameters, case$latent_gold_npar)
    TRUE
  }, logical(1))
})

test_that("a fit recovers its parameters and reports them in tables", {
  data <- extra_rows()
  fit <- extra_fit(data)
  expect_equal(fit$n_parameters, 11)
  ordinal <- get_results(fit, "ordinal")
  expect_identical(nrow(ordinal), 8L)
  expect_equal(as.vector(tapply(ordinal$probability, ordinal$profile, sum)), c(1, 1))
  counts <- get_results(fit, "count_means")
  expect_equal(sort(counts$mean), c(1, 5), tolerance = 0.1)
  expect_true(all(is.finite(counts$standard_error)))
  expect_output(print(fit), "ordinal o .*count k")
  inference <- parameter_inference(fit)
  expect_setequal(unique(inference$parameter),
                  c("mean", "variance", "ordinal_intercept", "ordinal_location",
                    "count_mean", "probability"))
  location <- subset(inference, parameter == "ordinal_location")
  expect_true(location$conf_low < 1.2 && 1.2 < location$conf_high)
})

test_that("analytic scores are the gradient, with missing values and weights", {
  data <- extra_rows(3)
  data$o[c(4, 9)] <- NA
  data$k[11] <- NA
  data$w <- runif(nrow(data), 0.5, 2)
  fit <- extra_fit(data, missing = "fiml", weights = "w")
  prepared <- .multilpa_inference_matrix(fit, data)
  theta <- unname(.multilpa_coefficients(fit, "unconstrained"))
  theta[1:2] <- theta[1:2] - prepared$centers
  set.seed(4)
  theta <- theta + rnorm(length(theta), 0, 0.03)
  objective <- function(point) {
    -.multilpa_expectation(prepared$x, fit$group_index, .multilpa_decode(point, fit),
                           prepared$codes, unname(fit$sampling_weights),
                           fit$extra_data)$log_likelihood
  }
  analytic <- .multilpa_score(theta, prepared$x, fit, prepared$codes)
  expect_equal(analytic, numeric_gradient(objective, theta), tolerance = 1e-5)
  expect_equal(-colSums(.multilpa_group_scores(theta, prepared$x, fit, prepared$codes)),
               analytic, tolerance = 1e-10)
})

test_that("integer weights reproduce duplicated rows", {
  data <- extra_rows(5, n = 300)
  data$w <- sample(1:3, nrow(data), TRUE)
  weighted <- extra_fit(data, weights = "w", tol = 1e-12)
  repeated <- extra_fit(data[rep(seq_len(nrow(data)), data$w), ], tol = 1e-12)
  expect_equal(weighted$log_likelihood * sum(data$w) / nrow(data),
               repeated$log_likelihood, tolerance = 1e-10)
  # The two fits may label the profiles in either order.
  expect_equal(sort(weighted$count_means), sort(repeated$count_means), tolerance = 1e-6)
})

test_that("prediction, bootstrap and simulation carry the new indicators", {
  data <- extra_rows(6, n = 400)
  fit <- extra_fit(data)
  predicted <- predict(fit, data)
  expect_equal(predicted$profile, unname(fit$subject_profiles))
  expect_error(predict(fit, transform(data, o = "sometimes")),
               class = "latents_bad_data")
  simulated <- .multilpa_draw_indicators(fit, rep(1:2, 200))
  expect_identical(levels(simulated$o), levels(data$o))
  expect_true(all(simulated$k >= 0))
  smaller <- quietly(lpa(data, c("y", "o", "k"), 1, ordinal = "o", count = "k"),
                     "latents_single_level")
  test <- quietly(bootstrap_lrt(smaller, fit, iter = 3, n_starts = 1, seed = 1))
  expect_true(is.finite(as.data.frame(test)$p_value))
  boot <- parameter_inference(fit, method = "bootstrap", iter = 3, n_starts = 1,
                              seed = 1)
  expect_identical(nrow(boot), nrow(parameter_inference(fit)))
})

test_that("invalid ordinal and count indicators are refused", {
  data <- extra_rows(7, n = 100)
  data$negative <- -data$k
  data$fraction <- data$k + 0.5
  data$text <- as.character(data$o)
  expect_error(quietly(lpa(data, c("y", "negative"), 2, count = "negative"),
                       "latents_single_level"), class = "latents_bad_data")
  expect_error(quietly(lpa(data, c("y", "fraction"), 2, count = "fraction"),
                       "latents_single_level"), class = "latents_bad_data")
  expect_error(quietly(lpa(data, c("y", "text"), 2, ordinal = "text"),
                       "latents_single_level"), class = "latents_bad_data")
  expect_error(quietly(lpa(data, c("y", "o"), 2, ordinal = "o", categorical = "o"),
                       "latents_single_level"), class = "latents_bad_argument")
  expect_error(quietly(lpa(data, "y", 2, count = "k"), "latents_single_level"),
               class = "latents_bad_argument")
  expect_error(quietly(lpa(data, c("y", "k"), 2, count = "k",
                           profile_covariates = "y"), "latents_single_level"),
               class = "latents_unsupported_indicator")
  expect_error(quietly(lpa(data, c("y", "k"), 2, count = "k", noise = TRUE),
                       "latents_single_level"),
               class = "latents_unsupported_indicator")
})
