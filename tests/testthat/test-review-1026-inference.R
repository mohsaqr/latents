# Review 2026-10-26, inference partition: Wald information, robust and OPG
# covariances, Rubin's rules, label alignment and seed scoping, checked against
# independent computations rather than against the package's own arithmetic.

review_two_profile_data <- function(seed = 11L, n_groups = 30L, per_group = 6L) {
  set.seed(seed)
  data <- data.frame(g = rep(seq_len(n_groups), each = per_group))
  membership <- stats::rbinom(nrow(data), 1L, 0.4)
  data$a <- stats::rnorm(nrow(data), ifelse(membership == 1L, 2, 0),
                         ifelse(membership == 1L, 0.7, 1))
  data$b <- stats::rnorm(nrow(data), ifelse(membership == 1L, -1, 0.5), 1)
  data
}

# Independent observed-data log likelihood of a two-profile, one-group-class
# diagonal (VVI) mixture in the package's unconstrained order: means by
# profile, log variances by profile, the profile-1 logit against profile 2.
review_vvi_loglik <- function(theta, x, rows = seq_len(nrow(x))) {
  means <- matrix(theta[1:4], 2L, byrow = TRUE)
  log_variances <- matrix(theta[5:8], 2L, byrow = TRUE)
  weights <- exp(c(theta[9], 0) - max(theta[9], 0))
  weights <- weights / sum(weights)
  density <- vapply(1:2, function(k) {
    weights[k] *
      stats::dnorm(x[rows, 1L], means[k, 1L], sqrt(exp(log_variances[k, 1L]))) *
      stats::dnorm(x[rows, 2L], means[k, 2L], sqrt(exp(log_variances[k, 2L])))
  }, numeric(length(rows)))
  sum(log(rowSums(matrix(density, ncol = 2L))))
}

test_that("observed, robust (CR1) and OPG covariances match independent numDeriv computations", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data <- review_two_profile_data()
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          variance_model = "varying", n_starts = 3, seed = 1, tol = 1e-10))
  x <- as.matrix(data[, c("a", "b")])
  theta <- unname(coef(fit, scale = "unconstrained"))
  expect_equal(review_vvi_loglik(theta, x), fit$log_likelihood, tolerance = 1e-10)

  observed <- parameter_inference(fit)
  hessian <- numDeriv::hessian(function(t) -review_vvi_loglik(t, x), theta)
  expect_lt(max(abs(hessian - unname(attr(observed, "hessian")))) / max(abs(hessian)), 1e-6)
  covariance <- solve(hessian)
  natural <- function(t) {
    weights <- exp(c(t[9], 0)) / sum(exp(c(t[9], 0)))
    c(t[1:4], exp(t[5:8]), weights, 1)
  }
  jacobian <- numDeriv::jacobian(natural, theta)
  delta_se <- sqrt(diag(jacobian %*% covariance %*% t(jacobian)))
  expect_equal(observed$standard_error, delta_se, tolerance = 1e-6)

  # Per-group scores of the independent likelihood, then the CR1 sandwich
  # G / (G - 1) * A^-1 B A^-1 (Cameron and Miller 2015) and the OPG B^-1.
  n_groups <- max(data$g)
  scores <- t(vapply(seq_len(n_groups), function(group) {
    numDeriv::grad(function(t) review_vvi_loglik(t, x, which(data$g == group)), theta)
  }, numeric(length(theta))))
  meat <- crossprod(scores)
  robust <- parameter_inference(fit, vcov_type = "robust")
  sandwich_covariance <- covariance %*% meat %*% covariance * n_groups / (n_groups - 1)
  expect_lt(max(abs(sandwich_covariance - unname(attr(robust, "covariance_unconstrained")))) /
              max(abs(sandwich_covariance)), 1e-6)
  expect_equal(attr(robust, "scaling_correction"),
               sum(diag(covariance %*% meat)) / length(theta), tolerance = 1e-6)
  opg <- parameter_inference(fit, vcov_type = "opg")
  expect_lt(max(abs(solve(meat) - unname(attr(opg, "covariance_unconstrained")))) /
              max(abs(solve(meat))), 1e-6)

  # Intervals: symmetric for means, log scale for variances, logit scale for
  # interior probabilities.
  critical <- stats::qnorm(0.975)
  mean_rows <- observed$parameter == "mean"
  expect_equal(observed$conf_low[mean_rows],
               observed$estimate[mean_rows] - critical * observed$standard_error[mean_rows])
  variance_rows <- observed$parameter == "variance"
  expect_equal(observed$conf_low[variance_rows],
               observed$estimate[variance_rows] *
                 exp(-critical * observed$standard_error[variance_rows] /
                       observed$estimate[variance_rows]))
  interior <- observed$parameter == "probability" & observed$estimate < 1
  p <- observed$estimate[interior]
  expect_equal(observed$conf_high[interior],
               stats::plogis(stats::qlogis(p) +
                               critical * observed$standard_error[interior] / (p * (1 - p))))
})

test_that("boundary = 'fix' inverts the Hessian of the free coordinates only", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set.seed(31)
  n <- 200L
  state <- sample.int(2L, n, TRUE)
  data <- data.frame(g = rep(1:40, each = 5), y = c(-2, 2)[state] + stats::rnorm(n))
  data$c <- factor(ifelse(state == 1L, sample(c("b", "c"), n, TRUE),
                          sample(c("a", "b", "c"), n, TRUE)))
  fit <- quietly(multilpa(data, c("y", "c"), "g", n_profiles = 2, n_group_classes = 1,
                          categorical = "c", n_starts = 3, seed = 1, tol = 1e-12,
                          max_iter = 5000))
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
  held <- parameter_inference(fit, boundary = "fix")
  expect_identical(attr(held, "fixed_at_bound"), "measurement.response.profile_1.c:a")

  theta <- unname(coef(fit, scale = "unconstrained"))
  free <- which(!grepl("profile_1\\.c:a$", names(coef(fit, scale = "unconstrained"))))
  code <- match(as.character(data$c), colnames(fit$response_probabilities$c))
  loglik <- function(t) {
    response <- t(apply(cbind(matrix(t[5:8], 2L, byrow = TRUE), 0), 1L,
                        function(z) exp(z - max(z)) / sum(exp(z - max(z)))))
    weights <- exp(c(t[9], 0)) / sum(exp(c(t[9], 0)))
    sum(log(rowSums(vapply(1:2, function(k) {
      weights[k] * stats::dnorm(data$y, t[k], sqrt(exp(t[2 + k]))) * response[k, code]
    }, numeric(n)))))
  }
  expect_equal(loglik(theta), fit$log_likelihood, tolerance = 1e-10)
  hessian <- numDeriv::hessian(function(t) {
    full <- theta
    full[free] <- t
    -loglik(full)
  }, theta[free])
  covariance <- solve(hessian)
  expect_lt(max(abs(covariance - unname(attr(held, "covariance_unconstrained")))) /
              max(abs(covariance)), 1e-5)
})

test_that("score-difference information recovers an analytic Hessian under both schemes", {
  set.seed(3)
  x <- stats::rnorm(40, 1, 2)
  # Normal log likelihood in (mu, eta = log sigma); score and -Hessian by hand.
  score <- function(theta) {
    sigma2 <- exp(2 * theta[2])
    c(sum(x - theta[1]) / sigma2, -length(x) + sum((x - theta[1])^2) / sigma2)
  }
  theta <- c(mu = 0.7, eta = 0.4)
  sigma2 <- exp(2 * theta[[2]])
  analytic <- matrix(c(length(x) / sigma2, 2 * sum(x - theta[[1]]) / sigma2,
                       2 * sum(x - theta[[1]]) / sigma2,
                       2 * sum((x - theta[[1]])^2) / sigma2), 2L)
  invisible(lapply(c("central", "five_point"), function(scheme) {
    invisible(lapply(c("plus1", "max1", "absolute"), function(rule) {
      information <- .inference_score_information(score, theta, 1e-4, scheme, rule)
      expect_equal(unname(information), analytic, tolerance = 1e-7)
      expect_identical(dimnames(information), list(names(theta), names(theta)))
    }))
  }))
  expect_error(.inference_invert(diag(c(1, -1))), class = "latents_singular_information")
  expect_error(.inference_invert(diag(c(1, 1e-12))), class = "latents_singular_information")
  expect_equal(.inference_invert(matrix(c(2, 1, 1, 2), 2L)), solve(matrix(c(2, 1, 1, 2), 2L)))
})

test_that("Rubin's rules match hand computation and mice::pool.scalar", {
  set.seed(4)
  m <- 5L
  estimates <- matrix(stats::rnorm(3 * m, c(1, 2, 0.5), 0.2), 3L)
  variances <- matrix(stats::runif(3 * m, 0.01, 0.05), 3L)
  estimates[3, ] <- 1
  variances[3, ] <- 0
  tables <- lapply(seq_len(m), function(i) {
    data.frame(level = "measurement", outcome = "profile_1", term = c("a", "b", "c"),
               parameter = "mean", estimate = estimates[, i],
               standard_error = sqrt(variances[, i]), statistic = c(1, NA, NA))
  })
  pooled <- .multilpa_rubin(tables, 0.95)
  q_bar <- mean(estimates[1, ])
  u_bar <- mean(variances[1, ])
  b <- stats::var(estimates[1, ])
  total <- u_bar + (1 + 1 / m) * b
  riv <- (1 + 1 / m) * b / u_bar
  df <- (m - 1) * (1 + 1 / riv)^2
  expect_equal(pooled$estimate[1], q_bar)
  expect_equal(pooled$standard_error[1], sqrt(total))
  expect_equal(pooled$df[1], df)
  expect_equal(pooled$fmi[1], (riv + 2 / (df + 3)) / (riv + 1))
  expect_equal(pooled$conf_low[1], q_bar - stats::qt(0.975, df) * sqrt(total))
  expect_equal(pooled$p_value[1], 2 * stats::pt(-abs(q_bar / sqrt(total)), df))
  expect_true(is.na(pooled$p_value[2]))
  # A structural constant: no variance, no missing information, infinite df.
  expect_identical(c(pooled$riv[3], pooled$fmi[3], pooled$df[3]), c(0, 0, Inf))
  skip_if_not_installed("mice")
  reference <- lapply(1:2, function(j) {
    mice::pool.scalar(estimates[j, ], variances[j, ], n = Inf, k = 1)
  })
  expect_equal(pooled$standard_error[1:2]^2, vapply(reference, `[[`, numeric(1), "t"))
  expect_equal(pooled$df[1:2], vapply(reference, `[[`, numeric(1), "df"))
  expect_equal(pooled$riv[1:2], vapply(reference, `[[`, numeric(1), "r"))
  expect_equal(pooled$fmi[1:2], vapply(reference, `[[`, numeric(1), "fmi"))
})

test_that("alignment recovers a known permutation and relabelling leaves the SE table invariant", {
  skip_on_cran()
  set.seed(8)
  n <- 240L
  state <- sample.int(3L, n, TRUE, prob = c(0.3, 0.3, 0.4))
  data <- data.frame(g = rep(1:40, each = 6),
                     a = c(-2, 0, 2)[state] + stats::rnorm(n, sd = c(1, 0.5, 1.5)[state]),
                     b = c(1, -1, 0)[state] + stats::rnorm(n))
  permutation <- c(3L, 1L, 2L)
  invisible(lapply(c("diagonal", "full"), function(covariance_model) {
    fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 3, n_group_classes = 1,
                            variance_model = "varying", covariance_model = covariance_model,
                            n_starts = 3, seed = 1))
    shuffled <- .multilpa_permute_profiles(fit, permutation)
    aligned <- .multilpa_pool_align(shuffled, fit)
    expect_identical(aligned$profile_order, order(permutation))
    expect_identical(coef(aligned$fit), coef(fit))
    expect_identical(aligned$fit$subject_profiles, fit$subject_profiles)
    original <- parameter_inference(fit)
    relabelled <- parameter_inference(shuffled)
    old_label <- function(outcome) {
      profile <- grepl("^profile_", outcome)
      outcome[profile] <- paste0("profile_",
                                 permutation[as.integer(sub("^profile_", "", outcome[profile]))])
      outcome
    }
    position <- match(paste(original$parameter, original$outcome, original$term),
                      paste(relabelled$parameter, old_label(relabelled$outcome),
                            relabelled$term))
    expect_false(anyNA(position))
    expect_equal(relabelled$standard_error[position], original$standard_error,
                 tolerance = 1e-6)
  }))
  # Equal means: only the variance features can tell the profiles apart.
  set.seed(9)
  state <- sample.int(2L, 300L, TRUE)
  spread <- data.frame(g = rep(1:50, each = 6),
                       a = stats::rnorm(300, 0, c(0.5, 3)[state]),
                       b = stats::rnorm(300, 0, c(0.5, 3)[state]))
  fit <- quietly(multilpa(spread, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          variance_model = "varying", n_starts = 3, seed = 1))
  signature <- .multilpa_profile_signature(fit)
  expect_true(all(is.finite(signature)))
  expect_identical(.multilpa_pool_align(.multilpa_permute_profiles(fit, 2:1), fit)$profile_order,
                   2:1)
})

test_that("standard errors are translation invariant and scale equivariant", {
  skip_on_cran()
  set.seed(17)
  n <- 240L
  state <- sample.int(2L, n, TRUE)
  data <- data.frame(g = rep(1:40, each = 6))
  data$a <- c(-1.5, 1.5)[state] + stats::rnorm(n)
  data$b <- c(1, -1)[state] + stats::rnorm(n)
  moved <- data
  moved$a <- moved$a * 1000 + 1e6
  fit_of <- function(frame) quietly(multilpa(frame, c("a", "b"), "g", n_profiles = 2,
                                             n_group_classes = 1, n_starts = 3, seed = 1,
                                             tol = 1e-11, variance_model = "varying",
                                             covariance_model = "full"))
  base <- parameter_inference(fit_of(data))
  scaled <- parameter_inference(fit_of(moved))
  power <- ifelse(base$parameter == "mean" & base$term %in% "a", 1,
                  ifelse(base$term %in% "a:a", 2, ifelse(base$term %in% "b:a", 1, 0)))
  expect_equal(scaled$standard_error / 1000^power, base$standard_error, tolerance = 1e-4)
})

test_that("a seeded bootstrap restores the caller's stream and is reproducible", {
  skip_on_cran()
  set.seed(1)
  data <- data.frame(g = rep(1:12, each = 6),
                     a = c(stats::rnorm(36, -1.5), stats::rnorm(36, 1.5)),
                     b = stats::rnorm(72))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          n_starts = 2, seed = 1))
  set.seed(99)
  expected <- stats::runif(3)
  set.seed(99)
  first <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3,
                                       n_starts = 1, seed = 7))
  expect_identical(stats::runif(3), expected)
  second <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3,
                                        n_starts = 1, seed = 7))
  expect_identical(attr(second, "replicates"), attr(first, "replicates"))
  set.seed(5)
  after_error <- stats::runif(1)
  set.seed(5)
  failing <- function() {
    .latents_local_seed(3)
    stop("boom")
  }
  expect_error(failing(), "boom")
  expect_identical(stats::runif(1), after_error)
})

test_that("starting variances keep profiles in rows with one indicator", {
  start <- list(covariances = array(c(2, 5), c(1, 1, 2)))
  expect_identical(.multilpa_start_variances(start), matrix(c(2, 5), 2L, 1L))
  start <- list(covariances = array(c(1, 0, 0, 3, 4, 0, 0, 6), c(2, 2, 2)))
  expect_identical(.multilpa_start_variances(start), matrix(c(1, 4, 3, 6), 2L, 2L))
})

# Regressions from the 2026-10-26 inference review (INF-1, INF-2).

test_that("covariate fits keep response standard errors for reserved category names", {
  set.seed(5)
  n <- 200L
  data <- data.frame(g = rep(1:40, each = 5), w = stats::rnorm(n))
  state <- stats::rbinom(n, 1, 0.5)
  data$a <- stats::rnorm(n, 2 * state)
  data$c <- factor(ifelse(stats::runif(n) < ifelse(state == 1, 0.7, 0.3), "x, y",
                          ifelse(stats::runif(n) < 0.5, "p:q", "z")))
  fit <- quietly(multilpa(data, c("a", "c"), "g", n_profiles = 2, n_group_classes = 1,
                          categorical = "c", profile_covariates = "w", n_starts = 3,
                          seed = 1))
  expect_false(anyNA(get_results(fit, "responses")$probability_standard_error))
})

test_that("a CSV round trip of a person-centred fit's data still aligns", {
  set.seed(1)
  data <- data.frame(g = rep(1:20, each = 5), a = stats::rnorm(100), b = stats::rnorm(100))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          centering = "person", n_starts = 1, seed = 1))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data, path, row.names = FALSE)
  expect_no_error(.multilpa_check_alignment(fit, utils::read.csv(path)))
})

test_that("a CSV round trip still refuses a changed person-centred value", {
  set.seed(1)
  data <- data.frame(g = rep(1:20, each = 5), a = stats::rnorm(100), b = stats::rnorm(100))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          centering = "person", n_starts = 1, seed = 1))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data, path, row.names = FALSE)
  changed <- utils::read.csv(path)
  changed$a[7L] <- changed$a[7L] + 1e-9
  expect_error(.multilpa_check_alignment(fit, changed),
               class = "latents_bad_inference_data")
  swapped <- utils::read.csv(path)
  swapped[c(1L, 2L), ] <- swapped[c(2L, 1L), ]
  expect_error(.multilpa_check_alignment(fit, swapped),
               class = "latents_bad_inference_data")
})
