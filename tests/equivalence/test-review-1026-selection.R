# Comparisons with independent references moved from tests/testthat/test-review-1026-selection.R,
# so that they run in CI's equivalence job rather than under R CMD check.

selection_fit <- function() {
  withr::local_seed(42)
  n <- 240
  truth <- sample(1:3, n, TRUE, prob = c(0.5, 0.3, 0.2))
  data <- data.frame(a = stats::rnorm(n, c(-2, 0, 2.5)[truth]),
                     b = stats::rnorm(n, c(0, 2, -1.5)[truth]))
  list(data = data, fit = quietly(lpa(data, c("a", "b"), n_profiles = 3,
                                      n_starts = 5, seed = 1)))
}

distal_case <- function() {
  withr::local_seed(11)
  n_groups <- 40
  group <- rep(seq_len(n_groups), each = 10)
  truth <- 1L + as.integer(stats::runif(length(group)) > 0.45)
  data <- data.frame(g = group,
                     a = stats::rnorm(length(group), ifelse(truth == 2L, 1, -1)),
                     b = stats::rnorm(length(group), ifelse(truth == 2L, 1, -1)),
                     y = stats::rnorm(length(group), ifelse(truth == 2L, 10, 0),
                                      ifelse(truth == 2L, 2, 1)))
  fit <- quietly(multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                          n_starts = 4, seed = 1))
  list(data = data, fit = fit)
}

test_that("every information criterion equals its definition (Nylund et al. 2007; Masyn 2013)", {
  case <- selection_fit()
  fit <- case$fit
  posterior <- .multilpa_classified_posteriors(fit)
  n <- nrow(posterior)
  log_lik <- fit$log_likelihood
  q <- fit$n_parameters
  entropy <- -sum(posterior[posterior > 0] * log(posterior[posterior > 0]))
  hand <- c(deviance = -2 * log_lik, aic = -2 * log_lik + 2 * q,
            kic = -2 * log_lik + 3 * (q + 1), bic = -2 * log_lik + q * log(n),
            sabic = -2 * log_lik + q * log((n + 2) / 24),
            caic = -2 * log_lik + q * (log(n) + 1),
            awe = -2 * (log_lik - entropy) + 2 * q * (1.5 + log(n)),
            icl = -2 * log_lik + q * log(n) + 2 * entropy,
            clc = -2 * log_lik + 2 * entropy)
  long <- .multilpa_information_criteria(fit, format = "long")
  individual <- subset(long, is.na(convention) | convention == "individuals")
  expect_equal(individual$value[match(names(hand), individual$criterion)],
               unname(hand), tolerance = 1e-12)
  wide <- .multilpa_information_criteria(fit)
  expect_equal(wide$bic_individual, hand[["bic"]], tolerance = 1e-12)
  expect_equal(wide$icl_individual, hand[["icl"]], tolerance = 1e-12)

  # mclust's convention is BIC = 2L - q log n, with its own parameter count.
  skip_if_not_installed("mclust")
  expect_identical(as.integer(q),
                   as.integer(mclust::nMclustParams("VVI", d = 2, G = 3)))
  expect_equal(-wide$bic_individual,
               mclust::bic("VVI", loglik = log_lik, n = n, d = 2, G = 3),
               tolerance = 1e-12)
})

test_that("ML step three matches an independent likelihood and sandwich (Bakk et al. 2013)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  case <- distal_case()
  posterior <- case$fit$subject_posteriors
  modal <- case$fit$subject_profiles
  y <- case$data$y
  groups <- case$data$g
  k <- ncol(posterior)
  n <- length(y)
  errors <- crossprod(posterior, outer(modal, seq_len(k), "==") * 1) / colSums(posterior)
  log_weights <- log(t(errors[, modal]) * matrix(colMeans(posterior), n, k, byrow = TRUE))
  unit_loglik <- function(theta, equal) {
    sds <- if (equal) rep(exp(theta[k + 1]), k) else exp(theta[k + seq_len(k)])
    joint <- log_weights + vapply(seq_len(k), function(class) {
      stats::dnorm(y, theta[class], sds[class], log = TRUE)
    }, numeric(n))
    peak <- apply(joint, 1L, max)
    peak + log(rowSums(exp(joint - peak)))
  }
  invisible(lapply(c(FALSE, TRUE), function(equal) {
    start <- c(tapply(y, modal, mean),
               if (equal) log(stats::sd(y) / 2) else log(tapply(y, modal, stats::sd)))
    objective <- function(theta) -sum(unit_loglik(theta, equal))
    found <- stats::optim(start, objective, method = "BFGS",
                          control = list(reltol = 1e-15, maxit = 5000))
    found <- stats::optim(found$par, objective, method = "Nelder-Mead",
                          control = list(reltol = 1e-15, maxit = 20000))
    found <- stats::optim(found$par, objective, method = "BFGS",
                          control = list(reltol = 1e-15, maxit = 5000))
    information <- -numDeriv::hessian(function(theta) sum(unit_loglik(theta, equal)),
                                      found$par)
    scores <- numDeriv::jacobian(function(theta) unit_loglik(theta, equal), found$par)
    influence <- (scores %*% solve(information))[, seq_len(k)]
    n_groups <- length(unique(groups))
    clustered <- rowsum(influence, groups) * sqrt(n_groups / (n_groups - 1))
    result <- three_step(case$fit, case$data, "y", method = "ml",
                         variance = if (equal) "equal" else "varying")
    variances <- if (equal) rep(exp(2 * found$par[k + 1]), k) else
      exp(2 * found$par[k + seq_len(k)])
    expect_lt(max(abs(result$estimate - found$par[seq_len(k)])), 1e-6)
    expect_lt(max(abs(result$variance / variances - 1)), 1e-6)
    expect_lt(max(abs(result$standard_error / sqrt(colSums(clustered^2)) - 1)), 1e-5)
  }))
})

test_that("ICC(1) equals the one-way ANOVA estimator with n0 for unbalanced groups", {
  withr::local_seed(5)
  sizes <- c(3, 8, 15, 4, 20, 6, 2, 11)
  groups <- rep(seq_along(sizes), sizes)
  n <- length(groups)
  values <- stats::rnorm(length(sizes), 0, 0.8)[groups] + stats::rnorm(n)
  anova_table <- stats::anova(stats::lm(values ~ factor(groups)))
  between <- anova_table[["Mean Sq"]][1]
  within <- anova_table[["Mean Sq"]][2]
  n0 <- (n - sum(sizes^2) / n) / (length(sizes) - 1)
  expect_equal(.multilpa_icc(values, groups),
               (between - within) / (between + (n0 - 1) * within), tolerance = 1e-12)
  # Invariant to the order of the rows and to the group labels.
  order <- sample(n)
  expect_equal(.multilpa_icc(values[order], groups[order]), .multilpa_icc(values, groups),
               tolerance = 1e-12)
  expect_equal(.multilpa_icc(values, letters[groups]), .multilpa_icc(values, groups),
               tolerance = 1e-12)
  # Distinct native identifiers that print alike are distinct groups.
  native <- c(0.3, 0.1 + 0.2)[1 + groups %% 2]
  expect_equal(.multilpa_icc(values, native), .multilpa_icc(values, 1 + groups %% 2),
               tolerance = 1e-12)
  expect_identical(.multilpa_icc(rep(0.1, n), groups), NA_real_)
  expect_equal(.multilpa_icc(stats::rnorm(length(sizes))[groups], groups), 1)
})

test_that("the weighted chi-square tail agrees with independent integrals", {
  # Positive weights: P(X1 - 0.3 X2 > 1) by conditioning on X2.
  conditioned <- stats::integrate(function(x) {
    stats::pchisq(1 + 0.3 * x, 1, lower.tail = FALSE) * stats::dchisq(x, 1)
  }, 0, Inf, rel.tol = 1e-13)$value
  expect_lt(abs(.latents_imhof_upper(1, c(1, -0.3)) - conditioned), 1e-12)
  # A negative statistic, split at the kink of the conditional probability.
  kink <- 0.5 / 0.3
  negative <- stats::pchisq(kink, 1) + stats::integrate(function(x) {
    stats::pchisq(-0.5 + 0.3 * x, 1, lower.tail = FALSE) * stats::dchisq(x, 1)
  }, kink, Inf, rel.tol = 1e-13)$value
  expect_lt(abs(.latents_imhof_upper(-0.5, c(1, -0.3)) - negative), 1e-12)
  expect_equal(.latents_imhof_upper(3, c(1, 1, 1)), stats::pchisq(3, 3, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("VLMR weights equal the eigenvalues of Vuong's W from independent derivatives", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  withr::local_seed(7)
  n <- 300
  class <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
  data <- data.frame(y = stats::rnorm(n, c(-1.5, 0.5, 2.5)[class], c(0.7, 1, 0.8)[class]))
  grid <- quietly(enumerate_lpa(data, "y", n_profiles = 2:3, model = "VVI",
                                n_starts = 10, seed = 1, tol = 1e-12))
  test <- get_results(grid, "lmr")
  coordinates <- function(fit) {
    k <- fit$n_profiles
    shares <- as.vector(fit$profile_probabilities)
    c(as.vector(fit$means), log(sqrt(as.vector(fit$variances))), log(shares[-k] / shares[k]))
  }
  unit_loglik <- function(theta, k) {
    eta <- c(theta[2 * k + seq_len(k - 1)], 0)
    log_share <- eta - max(eta)
    log_share <- log_share - log(sum(exp(log_share)))
    joint <- vapply(seq_len(k), function(j) {
      log_share[j] + stats::dnorm(data$y, theta[j], exp(theta[k + j]), log = TRUE)
    }, numeric(n))
    peak <- apply(joint, 1L, max)
    peak + log(rowSums(exp(joint - peak)))
  }
  null <- get_candidate(grid, 2)
  alternative <- get_candidate(grid, 3)
  null_point <- coordinates(null)
  alternative_point <- coordinates(alternative)
  expect_equal(sum(unit_loglik(null_point, 2)), null$log_likelihood, tolerance = 1e-10)
  scores_1 <- numDeriv::jacobian(function(theta) unit_loglik(theta, 3), alternative_point)
  scores_0 <- numDeriv::jacobian(function(theta) unit_loglik(theta, 2), null_point)
  hessian_1 <- numDeriv::hessian(function(theta) sum(unit_loglik(theta, 3)), alternative_point)
  hessian_0 <- numDeriv::hessian(function(theta) sum(unit_loglik(theta, 2)), null_point)
  cross <- crossprod(scores_1, scores_0)
  w <- rbind(cbind(-crossprod(scores_1) %*% solve(hessian_1), -cross %*% solve(hessian_0)),
             cbind(t(cross) %*% solve(hessian_1), crossprod(scores_0) %*% solve(hessian_0)))
  weights <- Re(eigen(w, only.values = TRUE)$values)
  expect_equal(test$reference_mean, sum(weights), tolerance = 1e-4)
  expect_equal(test$reference_sd, sqrt(2 * sum(weights^2)), tolerance = 1e-4)
  statistic <- 2 * (alternative$log_likelihood - null$log_likelihood)
  expect_equal(test$lr, statistic, tolerance = 1e-12)
  expect_equal(test$vlmr_p, .latents_imhof_upper(statistic, weights[abs(weights) > 1e-9]),
               tolerance = 1e-4)
  added <- alternative$n_parameters - null$n_parameters
  expect_equal(test$adjusted_lr, statistic / (1 + 1 / (added * log(n))), tolerance = 1e-12)
})
