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

imhof_note <- "Package side: the internal weighted chi-square tail (.latents_imhof_upper) behind the VLMR and adjusted LMR p-values of get_results(grid, \"lmr\")."

test_that("every information criterion equals its definition (Nylund et al. 2007; Masyn 2013)", {
  case <- selection_fit()
  fit <- case$fit
  posterior <- latents:::.multilpa_classified_posteriors(fit)
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
  definitions <- "independent: information criteria from their definitions (Nylund et al. 2007; Masyn 2013), base R"
  covers <- c("lpa()", "get_results.multilpa()")
  long <- get_results(fit, "information_criteria", format = "long")
  individual <- subset(long, is.na(convention) | convention == "individuals",
                       select = c(criterion, value))
  expect_reference(individual, data.frame(criterion = names(hand), value = unname(hand)),
                   reference = definitions, covers = covers,
                   quantity = "information criteria, individual convention, long table",
                   tolerance = 1e-12, measure = "relative", key = "criterion")
  wide <- get_results(fit, "information_criteria")
  expect_reference(c(bic = wide$bic_individual, icl = wide$icl_individual), hand[c("bic", "icl")],
                   reference = definitions, covers = covers,
                   quantity = "BIC and ICL, individual convention, wide table",
                   tolerance = 1e-12, measure = "relative")

  # mclust's convention is BIC = 2L - q log n, with its own parameter count.
  skip_without_reference("mclust")
  expect_reference(as.integer(q), as.integer(mclust::nMclustParams("VVI", d = 2, G = 3)),
                   reference = "mclust", covers = c("lpa()", "logLik.multilpa()"),
                   quantity = "number of parameters, VVI, three profiles, two indicators",
                   tolerance = 0)
  expect_reference(-wide$bic_individual,
                   mclust::bic("VVI", loglik = log_lik, n = n, d = 2, G = 3),
                   reference = "mclust", covers = covers,
                   quantity = "BIC in mclust's sign convention, VVI, three profiles",
                   tolerance = 1e-12, measure = "relative")
})

test_that("ML step three matches an independent likelihood and sandwich (Bakk et al. 2013)", {
  skip_on_cran()
  skip_without_reference("numDeriv")
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
  # The independent ML step three, for one variance model: means and
  # variances from optim() on the likelihood above, cluster-robust standard
  # errors from numDeriv's information and scores; and the package's table.
  step_three <- function(equal) {
    label <- if (equal) "equal" else "varying"
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
    result <- as.data.frame(three_step(case$fit, case$data, "y", method = "ml", variance = label))
    list(package = data.frame(variance_model = label, class = result$class,
                              estimate = result$estimate, variance = result$variance,
                              standard_error = result$standard_error),
         reference = data.frame(variance_model = label, class = seq_len(k),
                                estimate = unname(found$par[seq_len(k)]),
                                variance = if (equal) rep(exp(2 * found$par[k + 1]), k) else
                                  unname(exp(2 * found$par[k + seq_len(k)])),
                                standard_error = sqrt(colSums(clustered^2))))
  }
  varying <- step_three(FALSE)
  equal <- step_three(TRUE)
  package <- rbind(varying$package, equal$package)
  reference <- rbind(varying$reference, equal$reference)
  keys <- c("variance_model", "class")
  covers <- c("three_step()", "as.data.frame.latents_table()")
  likelihood <- "independent: ML step-three likelihood (Bakk et al. 2013), base R, maximised by stats::optim()"
  # optim() stops at its own criterion: the original bounds, 1e-6 absolute on
  # the means and 1e-6 relative on the variances.
  expect_reference(package, reference[c(keys, "estimate")], reference = likelihood, covers = covers,
                   quantity = "class means of the distal outcome", tolerance = 1e-6,
                   measure = "absolute", key = keys)
  expect_reference(package, reference[c(keys, "variance")], reference = likelihood, covers = covers,
                   quantity = "class variances of the distal outcome", tolerance = 1e-6,
                   measure = "relative", key = keys)
  # numDeriv's Hessian and Jacobian at optim's point: 1e-5 relative, the
  # original bound on |ratio - 1|.
  expect_reference(package, reference[c(keys, "standard_error")], reference = "numDeriv", covers = covers,
                   quantity = "cluster-robust standard errors of the class means", tolerance = 1e-5,
                   measure = "relative", key = keys,
                   note = "Reference: the sandwich of numDeriv::hessian() and numDeriv::jacobian() of the independent step-three likelihood, clustered by group.")
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
  expect_reference(latents:::.multilpa_icc(values, groups),
                   (between - within) / (between + (n0 - 1) * within),
                   reference = "stats", covers = c("descriptives()", "descriptives.data.frame()"),
                   quantity = "ICC(1) from the one-way ANOVA mean squares, unbalanced groups",
                   tolerance = 1e-12, measure = "relative",
                   note = "Package side: the internal ICC (.multilpa_icc) behind the icc column of descriptives(); reference: stats::anova(stats::lm()) mean squares with n0.")
  # Invariant to the order of the rows and to the group labels.
  order <- sample(n)
  expect_equal(latents:::.multilpa_icc(values[order], groups[order]), latents:::.multilpa_icc(values, groups),
               tolerance = 1e-12)
  expect_equal(latents:::.multilpa_icc(values, letters[groups]), latents:::.multilpa_icc(values, groups),
               tolerance = 1e-12)
  # Distinct native identifiers that print alike are distinct groups.
  native <- c(0.3, 0.1 + 0.2)[1 + groups %% 2]
  expect_equal(latents:::.multilpa_icc(values, native), latents:::.multilpa_icc(values, 1 + groups %% 2),
               tolerance = 1e-12)
  expect_identical(latents:::.multilpa_icc(rep(0.1, n), groups), NA_real_)
  expect_equal(latents:::.multilpa_icc(stats::rnorm(length(sizes))[groups], groups), 1)
})

test_that("the weighted chi-square tail agrees with independent integrals", {
  # Positive weights: P(X1 - 0.3 X2 > 1) by conditioning on X2.
  conditioned <- stats::integrate(function(x) {
    stats::pchisq(1 + 0.3 * x, 1, lower.tail = FALSE) * stats::dchisq(x, 1)
  }, 0, Inf, rel.tol = 1e-13)$value
  # The Imhof integral's accuracy is absolute: the original bound, 1e-12.
  expect_reference(latents:::.latents_imhof_upper(1, c(1, -0.3)), conditioned,
                   reference = "independent: conditioning on one chi-square, stats::integrate()",
                   covers = "get_results.multilpa_enumeration()",
                   quantity = "P(X1 - 0.3 X2 > 1)", tolerance = 1e-12, measure = "absolute",
                   note = imhof_note)
  # A negative statistic, split at the kink of the conditional probability.
  kink <- 0.5 / 0.3
  negative <- stats::pchisq(kink, 1) + stats::integrate(function(x) {
    stats::pchisq(-0.5 + 0.3 * x, 1, lower.tail = FALSE) * stats::dchisq(x, 1)
  }, kink, Inf, rel.tol = 1e-13)$value
  expect_reference(latents:::.latents_imhof_upper(-0.5, c(1, -0.3)), negative,
                   reference = "independent: conditioning on one chi-square, stats::integrate()",
                   covers = "get_results.multilpa_enumeration()",
                   quantity = "P(X1 - 0.3 X2 > -0.5)", tolerance = 1e-12, measure = "absolute",
                   note = imhof_note)
  expect_reference(latents:::.latents_imhof_upper(3, c(1, 1, 1)), stats::pchisq(3, 3, lower.tail = FALSE),
                   reference = "stats", covers = "get_results.multilpa_enumeration()",
                   quantity = "P(chi-square 3 > 3), three unit weights", tolerance = 1e-12,
                   measure = "relative", note = imhof_note)
})

test_that("VLMR weights equal the eigenvalues of Vuong's W from independent derivatives", {
  skip_on_cran()
  skip_without_reference("numDeriv")
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
  expect_reference(as.numeric(logLik(null)), sum(unit_loglik(null_point, 2)),
                   reference = "independent: univariate normal mixture likelihood, base R",
                   covers = c("enumerate_lpa()", "get_candidate()", "logLik.multilpa()"),
                   quantity = "log likelihood of the two-profile candidate, VVI",
                   tolerance = 1e-10, measure = "relative")
  scores_1 <- numDeriv::jacobian(function(theta) unit_loglik(theta, 3), alternative_point)
  scores_0 <- numDeriv::jacobian(function(theta) unit_loglik(theta, 2), null_point)
  hessian_1 <- numDeriv::hessian(function(theta) sum(unit_loglik(theta, 3)), alternative_point)
  hessian_0 <- numDeriv::hessian(function(theta) sum(unit_loglik(theta, 2)), null_point)
  cross <- crossprod(scores_1, scores_0)
  w <- rbind(cbind(-crossprod(scores_1) %*% solve(hessian_1), -cross %*% solve(hessian_0)),
             cbind(t(cross) %*% solve(hessian_1), crossprod(scores_0) %*% solve(hessian_0)))
  weights <- Re(eigen(w, only.values = TRUE)$values)
  statistic <- 2 * (alternative$log_likelihood - null$log_likelihood)
  # The weights from numDeriv's scores and Hessians of the likelihood above;
  # 1e-4 relative, the original test's tolerance (numerical second
  # derivatives against the package's analytic ones).
  expect_reference(c(reference_mean = test$reference_mean, reference_sd = test$reference_sd),
                   c(reference_mean = sum(weights), reference_sd = sqrt(2 * sum(weights^2))),
                   reference = "numDeriv",
                   covers = c("enumerate_lpa()", "get_results.multilpa_enumeration()"),
                   quantity = "mean and standard deviation of the VLMR reference distribution",
                   tolerance = 1e-4, measure = "relative",
                   note = "Reference: eigenvalues of Vuong's W from numDeriv::jacobian()/hessian() of an independent normal mixture likelihood.")
  expect_equal(test$lr, statistic, tolerance = 1e-12)
  expect_reference(test$vlmr_p, latents:::.latents_imhof_upper(statistic, weights[abs(weights) > 1e-9]),
                   reference = "numDeriv",
                   covers = c("enumerate_lpa()", "get_results.multilpa_enumeration()"),
                   quantity = "VLMR p-value, two against three profiles", tolerance = 1e-4,
                   measure = "relative",
                   note = "Reference: the package's Imhof tail (checked against stats and integrate() above) at the eigenvalues of Vuong's W from numDeriv derivatives of an independent likelihood.")
  added <- alternative$n_parameters - null$n_parameters
  expect_equal(test$adjusted_lr, statistic / (1 + 1 / (added * log(n))), tolerance = 1e-12)
})
