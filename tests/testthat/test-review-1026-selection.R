# Review 2026-10-26 (selection partition): calibration of the model-selection
# statistics and the three-step distal-outcome estimators against hand
# computations from the fitted posteriors and likelihoods, plus invariants.

selection_fit <- function() {
  withr::local_seed(42)
  n <- 240
  truth <- sample(1:3, n, TRUE, prob = c(0.5, 0.3, 0.2))
  data <- data.frame(a = stats::rnorm(n, c(-2, 0, 2.5)[truth]),
                     b = stats::rnorm(n, c(0, 2, -1.5)[truth]))
  list(data = data, fit = quietly(lpa(data, c("a", "b"), n_profiles = 3,
                                      n_starts = 5, seed = 1)))
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

test_that("entropy, AvePP and odds of correct classification equal their definitions", {
  fit <- selection_fit()$fit
  posterior <- .multilpa_classified_posteriors(fit)
  n <- nrow(posterior)
  k <- ncol(posterior)
  entropy <- -sum(posterior[posterior > 0] * log(posterior[posterior > 0]))
  table <- .multilpa_entropy_table(fit)
  individual <- subset(table, level == "individuals")
  expect_equal(individual$entropy_sum, entropy, tolerance = 1e-12)
  expect_equal(individual$relative_entropy, 1 - entropy / (n * log(k)), tolerance = 1e-12)

  modal <- max.col(posterior, ties.method = "first")
  avepp <- vapply(seq_len(k), function(class) mean(posterior[modal == class, class]),
                  numeric(1))
  share <- colMeans(posterior)
  occ <- (avepp / (1 - avepp)) / (share / (1 - share))
  quality <- .multilpa_classification_table(fit)
  expect_equal(quality$n_modal, tabulate(modal, nbins = k))
  expect_equal(quality$average_posterior, avepp, tolerance = 1e-12)
  expect_equal(quality$estimated_proportion, unname(share), tolerance = 1e-12)
  expect_equal(quality$odds_correct_classification, unname(occ), tolerance = 1e-10)
})

test_that("classification statistics are invariant to relabelling the classes", {
  fit <- selection_fit()$fit
  posterior <- .multilpa_classified_posteriors(fit)
  order <- c(3L, 1L, 2L)
  relabelled <- posterior[, order]
  expect_equal(.multilpa_entropy_sum(relabelled), .multilpa_entropy_sum(posterior),
               tolerance = 1e-14)
  expect_equal(.multilpa_relative_entropy(relabelled),
               .multilpa_relative_entropy(posterior), tolerance = 1e-14)
  expect_equal(.multilpa_average_posterior_matrix(relabelled),
               .multilpa_average_posterior_matrix(posterior)[order, order],
               tolerance = 1e-14)
  original <- .multilpa_classification_summary(posterior, "individuals")
  permuted <- .multilpa_classification_summary(relabelled, "individuals")
  expect_equal(permuted$odds_correct_classification,
               original$odds_correct_classification[order], tolerance = 1e-12)
  # The rows' order does not matter either.
  rows <- rev(seq_len(nrow(posterior)))
  expect_equal(.multilpa_entropy_sum(posterior[rows, ]), .multilpa_entropy_sum(posterior),
               tolerance = 1e-12)
})

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

test_that("BCH means and clustered standard errors equal a hand computation (Vermunt 2010)", {
  case <- distal_case()
  posterior <- case$fit$subject_posteriors
  modal <- case$fit$subject_profiles
  y <- case$data$y
  groups <- case$data$g
  k <- ncol(posterior)
  # D[t, s] = P(assigned s | true t), from the posteriors.
  errors <- t(vapply(seq_len(k), function(true) {
    vapply(seq_len(k), function(assigned) {
      sum(posterior[, true] * (modal == assigned)) / sum(posterior[, true])
    }, numeric(1))
  }, numeric(k)))
  weights <- solve(errors)[modal, ]
  means <- colSums(weights * y) / colSums(weights)
  influence <- sweep(weights * outer(y, means, "-"), 2L, colSums(weights), "/")
  n_groups <- length(unique(groups))
  clustered <- rowsum(influence, groups) * sqrt(n_groups / (n_groups - 1))
  result <- three_step(case$fit, case$data, "y")
  expect_equal(result$estimate, unname(means), tolerance = 1e-12)
  expect_equal(result$standard_error, unname(sqrt(colSums(clustered^2))), tolerance = 1e-12)
  expect_equal(result$effective_n, unname(colSums(weights)^2 / colSums(weights^2)),
               tolerance = 1e-12)
  pairs <- three_step(case$fit, case$data, "y", contrast = "pairs")
  expect_equal(pairs$estimate, unname(means[2] - means[1]), tolerance = 1e-12)
  expect_equal(pairs$standard_error,
               sqrt(sum((clustered[, 2] - clustered[, 1])^2)), tolerance = 1e-12)
  stored <- subset(get_results(case$fit, "bch_weights"), level == "individuals")
  expect_equal(stored$weight, as.vector(t(weights)), tolerance = 1e-12)
  # Proportional and modal estimates are the posterior- and modal-weighted means.
  expect_equal(three_step(case$fit, case$data, "y", method = "proportional")$estimate,
               unname(colSums(posterior * y) / colSums(posterior)), tolerance = 1e-12)
  expect_equal(three_step(case$fit, case$data, "y", method = "modal")$estimate,
               as.vector(tapply(y, modal, mean)), tolerance = 1e-12)
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

test_that("an unexpected error in the equal-variance BFGS search is not swallowed", {
  # Regression: the BFGS candidates were wrapped in error = function(e) NULL,
  # which hid every failure. The first log-likelihood evaluation of an
  # equal-variance fit is a BFGS start; an error there must reach the caller.
  pieces <- list(posteriors = rbind(cbind(rep(0.9, 30), 0.1), cbind(rep(0.1, 30), 0.9)),
                 modal = rep(1:2, each = 30), n_classes = 2L)
  values <- c(rep(c(-2, -1, 0), 10), rep(c(3, 4, 5), 10))
  calls <- 0L
  original <- .multilpa_step_ml_log_likelihood
  local_mocked_bindings(.multilpa_step_ml_log_likelihood = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) {
      stop(errorCondition("injected", class = "latents_test_injected"))
    }
    original(...)
  })
  expect_error(.multilpa_step_ml_fit(pieces, values, "equal"),
               class = "latents_test_injected")
})

test_that("a constant outcome is refused by every step-three method", {
  case <- distal_case()
  data <- case$data
  # Exactly constant values whose naive mean could round (the constancy
  # pattern repaired elsewhere): var() is exactly zero for them.
  data$y <- 0.1
  expect_error(three_step(case$fit, data, "y", method = "ml"),
               class = "latents_bad_outcome")
  constant <- three_step(case$fit, data, "y")
  expect_identical(constant$estimate, rep(0.1, 2))
  expect_identical(constant$standard_error, c(0, 0))
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

test_that("the bootstrap p-value counts ties and the observed statistic", {
  replicates <- data.frame(valid = TRUE, statistic = as.numeric(1:9))
  expect_equal(.latents_blrt_p_value(replicates, 6.5, 9L)$p_value, (1 + 3) / 10)
  expect_equal(.latents_blrt_p_value(replicates, 6, 9L)$p_value, (1 + 4) / 10)
  expect_equal(.latents_blrt_p_value(replicates, 100, 9L)$p_value, 1 / 10)
  # An invalid replicate is left out (round 2: p over the valid replicates,
  # NA only below half valid); replicate 2 (statistic 2) no longer counts.
  replicates$valid[2] <- FALSE
  expect_warning(kept <- .latents_blrt_p_value(replicates, 6, 9L),
                 class = "latents_failed_replicates")
  expect_equal(kept$p_value, (1 + 4) / 9)
  replicates$valid[1:5] <- FALSE
  expect_warning(withheld <- .latents_blrt_p_value(replicates, 6, 9L),
                 class = "latents_failed_replicates")
  expect_identical(withheld$p_value, NA_real_)
})

test_that("the knee is the point furthest below the chord in rescaled units", {
  x <- 1:6
  y <- c(100, 60, 40, 35, 33, 32)
  knee <- .latents_kneedle(x, y)
  xs <- (x - 1) / 5
  ys <- (y - 32) / 68
  slope <- ys[6] - ys[1]
  expect_equal(knee$distance, ((ys[1] + slope * xs) - ys) / sqrt(1 + slope^2),
               tolerance = 1e-14)
  expect_identical(knee$at, 3L)
  # Rescaling the criterion does not move the knee.
  expect_identical(.latents_kneedle(x, 1e6 + 7 * y)$at, knee$at)
})

test_that("comparison identity ignores row names that only one engine keeps", {
  # Regression: a homogeneous transition fit keeps the frame's row names on
  # its indicators and the general engine does not, so compare_models()
  # refused two fits of the same data as "fitted to different data".
  data <- subset(course_engagement, student <= 12)
  homogeneous <- quietly(lta(data, c("browse", "lectures"), "student", n_profiles = 2,
                             time = "sequence", n_starts = 1, seed = 1, max_iter = 200))
  occasion <- quietly(lta(data, c("browse", "lectures"), "student", n_profiles = 2,
                          time = "sequence", n_starts = 1, seed = 1, max_iter = 200,
                          transitions = "occasion"))
  expect_s3_class(occasion, "multilpa_lta")
  compared <- compare_models(homogeneous, occasion)
  expect_identical(nrow(as.data.frame(compared)), 2L)
  # Different values are still refused.
  changed <- data
  changed$browse[1] <- changed$browse[1] + 1e-9
  other <- quietly(lta(changed, c("browse", "lectures"), "student", n_profiles = 2,
                       time = "sequence", n_starts = 1, seed = 1, max_iter = 200,
                       transitions = "occasion"))
  expect_error(compare_models(homogeneous, other), class = "latents_incomparable_models")
})

test_that("r3step reports NA, not 0, for a coefficient without positive variance", {
  frame <- data.frame(level = "profile", outcome = "profile_1", term = "x")
  covariance <- diag(c(0.04, -1e-9, 0.09))
  expect_warning(
    table <- .multilpa_r3step_frame(c(1, 2, 3), covariance, c("(Intercept)", "x", "z"),
                                    1L, 2L, "profile", 0.95, "observed", "none"),
    class = "latents_singular_information")
  expect_equal(table$standard_error, c(0.2, NA, 0.3))
  expect_true(is.na(table$p_value[[2L]]))
  expect_true(is.na(table$conf_low[[2L]]))
})
