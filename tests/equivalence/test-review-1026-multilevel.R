# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-review-1026-multilevel.R,
# so that they run in CI's equivalence job rather than under R CMD check.

review_ml_data <- function(seed = 11, n_groups = 12, size = 4) {
  withr::local_seed(seed)
  g <- rep(seq_len(n_groups), each = size)
  h <- rep(sample(1:2, n_groups, TRUE), each = size)
  z <- ifelse(stats::runif(n_groups * size) < ifelse(h == 1, 0.8, 0.2), 1, 2)
  data.frame(
    g = g,
    y1 = stats::rnorm(n_groups * size, c(0, 2.5)[z]),
    y2 = stats::rnorm(n_groups * size, c(0, 1.5)[z], c(1, 0.6)[z]),
    wt = rep(sample(1:3, n_groups, TRUE), each = size)
  )
}

review_ml_covariate_data <- function(seed = 3, n_groups = 20, size = 6) {
  withr::local_seed(seed)
  g <- rep(seq_len(n_groups), each = size)
  gz <- rep(stats::rnorm(n_groups), each = size)
  hc <- rep(ifelse(stats::runif(n_groups) < stats::plogis(1.2 * gz[!duplicated(g)]),
                   2, 1), each = size)
  x <- stats::rnorm(n_groups * size)
  z <- ifelse(stats::runif(n_groups * size) < stats::plogis(-1 + 2 * (hc == 2) + 0.8 * x),
              2, 1)
  data.frame(
    g = g, x = x, gz = gz,
    wt = rep(stats::runif(n_groups, 0.5, 2), each = size),
    y1 = stats::rnorm(n_groups * size, c(-1.5, 1.5)[z], 0.8),
    y2 = stats::rnorm(n_groups * size, c(-1, 1)[z], c(0.7, 1)[z])
  )
}

# Row densities of a diagonal Gaussian mixture, computed with dnorm() alone.
review_ml_densities <- function(y, means, variances) {
  vapply(seq_len(nrow(means)), function(k) {
    Reduce(`*`, lapply(seq_len(ncol(y)), function(j) {
      stats::dnorm(y[, j], means[k, j], sqrt(variances[k, j]))
    }))
  }, numeric(nrow(y)))
}

review_ml_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, c("latents_boundary", "latents_small_classes",
                      "latents_thin_profiles"))) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("the nested likelihood and both posteriors equal brute-force enumeration", {
  data <- review_ml_data()
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", n_profiles = 2,
                                  n_group_classes = 2, n_starts = 3, seed = 1))
  y <- as.matrix(data[c("y1", "y2")])
  density <- review_ml_densities(y, unname(as.matrix(fit$means)),
                                 unname(as.matrix(fit$variances)))
  profile_p <- unname(as.matrix(fit$profile_probabilities))
  omega <- unname(fit$group_probabilities)
  configurations <- as.matrix(expand.grid(rep(list(1:2), 4)))
  # Every (group class, profile of each member) assignment of every group.
  per_group <- lapply(seq_len(12), function(j) {
    rows <- which(data$g == j)
    joint <- lapply(1:2, function(h) {
      apply(configurations, 1L, function(cfg) {
        omega[h] * prod(profile_p[h, cfg] * density[cbind(rows, cfg)])
      })
    })
    total <- sum(unlist(joint))
    list(
      log_likelihood = log(total),
      group = vapply(joint, sum, numeric(1)) / total,
      subject = t(vapply(seq_along(rows), function(i) {
        vapply(1:2, function(k) {
          sum(vapply(joint, function(v) sum(v[configurations[, i] == k]), numeric(1)))
        }, numeric(1)) / total
      }, numeric(2))))
  })
  enumeration <- "independent: brute-force enumeration of class and profile assignments, stats::dnorm()"
  # The same probabilities summed in another order: rounding only.
  expect_reference(as.numeric(logLik(fit)),
                   sum(vapply(per_group, `[[`, numeric(1), "log_likelihood")),
                   reference = enumeration, covers = c("multilpa()", "logLik.multilpa()"),
                   quantity = "nested log likelihood", tolerance = 1e-12, measure = "relative")
  expect_reference(unname(fit$group_posteriors),
                   do.call(rbind, lapply(per_group, `[[`, "group")),
                   reference = enumeration, covers = "multilpa()",
                   quantity = "group class posteriors", tolerance = 1e-12, measure = "relative")
  expect_reference(unname(fit$subject_posteriors),
                   do.call(rbind, lapply(per_group, `[[`, "subject")),
                   reference = enumeration, covers = "multilpa()",
                   quantity = "profile posteriors of the members", tolerance = 1e-12, measure = "relative")
})

test_that("one group class is the single-level mixture, at least as high as mclust", {
  skip_without_reference("mclust")
  data <- review_ml_data(seed = 21, n_groups = 15, size = 5)
  two_level <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", n_profiles = 2,
                                        n_group_classes = 1, n_starts = 3, seed = 1,
                                        tol = 1e-12))
  single <- review_ml_quiet(lpa(data, c("y1", "y2"), 2, n_starts = 3, seed = 1,
                                tol = 1e-12))
  expect_equal(two_level$log_likelihood, single$log_likelihood, tolerance = 1e-10)
  density <- review_ml_densities(as.matrix(data[c("y1", "y2")]),
                                 unname(as.matrix(single$means)),
                                 unname(as.matrix(single$variances)))
  proportions <- unname(as.vector(single$profile_probabilities))
  expect_reference(as.numeric(logLik(single)), sum(log(density %*% proportions)),
                   reference = "independent: diagonal Gaussian mixture likelihood, stats::dnorm()",
                   covers = c("lpa()", "logLik.multilpa()"),
                   quantity = "single-level log likelihood at the reported estimates",
                   tolerance = 1e-12, measure = "relative")
  # mclust's own EM for VVI, started from a split on y1 and run to a tight
  # tolerance, reaches the same maximum.
  y <- as.matrix(data[c("y1", "y2")])
  reference <- mclust::meVVI(y, mclust::unmap(ifelse(y[, "y1"] > 1.25, 2, 1)),
                             control = mclust::emControl(
                               tol = c(1e-12, sqrt(.Machine$double.eps))))
  # 1e-7, the original tolerance: two EMs stopped by their own criteria.
  expect_reference(as.numeric(logLik(single)), reference$loglik, reference = "mclust",
                   covers = c("lpa()", "logLik.multilpa()"),
                   quantity = "VVI log likelihood, two profiles",
                   tolerance = 1e-7, measure = "relative",
                   note = "Reference: mclust::meVVI() started from a split on y1.")
  expect_reference(as.numeric(single$n_parameters), as.numeric(mclust::nMclustParams("VVI", 2, 2)),
                   reference = "mclust", covers = "lpa()",
                   quantity = "number of parameters, VVI, two profiles",
                   tolerance = 0, measure = "absolute")
})

# Soft class counts with two covariates, one far from zero (x2), for the
# membership Newton step.
review_newton_counts <- function() {
  withr::local_seed(4)
  n <- 200
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n, 5, 0.3)
  design <- cbind(1, x1, x2)
  eta <- cbind(0.5 + x1, -0.3 + 0.5 * (x2 - 5), 0)
  probabilities <- exp(eta) / rowSums(exp(eta))
  list(x1 = x1, x2 = x2, design = design, counts = probabilities * stats::runif(n, 0.5, 2))
}

test_that("the membership Newton step reaches an independent optimizer's minimum", {
  newton <- review_newton_counts()
  solved <- latents:::.multilpa_weighted_logits(newton$design, newton$counts, matrix(0, 3, 2))
  expect_true(solved$converged)
  objective <- function(b) {
    -sum(newton$counts * latents:::.multilpa_log_softmax(newton$design, matrix(b, 3, 2)))
  }
  reference <- stats::optim(rep(0, 6), objective, method = "BFGS",
                            control = list(reltol = 1e-16, maxit = 5000))
  # One-sided: the Newton step must be at least as low as BFGS's minimum.
  expect_lte(objective(as.vector(solved$coefficients)), reference$value + 1e-10)
  # 1e-2, the original tolerance: x2 sits at 5 with SD 0.3, so BFGS stops far
  # from the minimum along the ill-conditioned intercept/x2 direction.
  # Scaled: some coefficients are near zero, where BFGS's error is absolute.
  expect_reference(as.vector(solved$coefficients), reference$par, reference = "stats",
                   covers = c("multilpa()", "lta()"),
                   quantity = "membership logit coefficients, three categories",
                   tolerance = 1e-2, measure = "scaled",
                   key = paste(rep(c("(Intercept)", "x1", "x2"), 2L), rep(c("class 1", "class 2"), each = 3L)),
                   note = "Package side: the internal weighted multinomial-logit Newton step (.multilpa_weighted_logits) of the covariate M-steps of multilpa() and lta(); reference: stats::optim(method = 'BFGS') on the weighted multinomial log likelihood.")
})

test_that("the membership Newton step with two categories is glm's weighted logit", {
  newton <- review_newton_counts()
  binary <- cbind(newton$counts[, 1], newton$counts[, 2] + newton$counts[, 3])
  two <- latents:::.multilpa_weighted_logits(newton$design, binary, matrix(0, 3, 1))
  x1 <- newton$x1
  x2 <- newton$x2
  glm_fit <- stats::glm(binary[, 1] / rowSums(binary) ~ x1 + x2,
                        family = stats::quasibinomial(), weights = rowSums(binary),
                        control = stats::glm.control(epsilon = 1e-14, maxit = 100))
  expect_reference(setNames(as.vector(two$coefficients), names(stats::coef(glm_fit))),
                   stats::coef(glm_fit), reference = "stats",
                   covers = c("multilpa()", "lta()"),
                   quantity = "membership logit coefficients, two categories",
                   tolerance = 1e-8, measure = "relative",
                   note = "Package side: the internal weighted logit Newton step (.multilpa_weighted_logits) of the covariate M-steps of multilpa() and lta(); reference: stats::glm(family = quasibinomial()) with weights.")
})

test_that("the covariate likelihood equals an independent two-level computation", {
  data <- review_ml_covariate_data()
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                  profile_covariates = "x", group_covariates = "gz",
                                  n_starts = 2, seed = 1, tol = 1e-10))
  softmax <- function(eta) {
    e <- exp(cbind(eta, 0) - apply(cbind(eta, 0), 1L, max))
    e / rowSums(e)
  }
  beta <- unname(fit$profile_coefficients)
  gamma <- unname(fit$group_coefficients)
  density <- review_ml_densities(as.matrix(data[c("y1", "y2")]),
                                 unname(fit$means), unname(fit$variances))
  first <- !duplicated(data$g)
  omega <- softmax(cbind(1, data$gz[first]) %*% gamma)
  per_group <- vapply(seq_len(20), function(j) {
    rows <- which(data$g == j)
    sum(vapply(1:2, function(h) {
      prior <- softmax(cbind(matrix(as.numeric(1:2 == h), length(rows), 2L, byrow = TRUE),
                             data$x[rows]) %*% beta)
      omega[j, h] * prod(rowSums(prior * density[rows, ]))
    }, numeric(1)))
  }, numeric(1))
  expect_reference(as.numeric(logLik(fit)), sum(log(per_group)),
                   reference = "independent: two-level covariate mixture likelihood, stats::dnorm()",
                   covers = c("multilpa()", "logLik.multilpa_covariates()"),
                   quantity = "log likelihood with profile and group covariates",
                   tolerance = 1e-12, measure = "relative")
  # (G - 1)(1 + q) + (G + p)(K - 1) membership and 2 K J measurement parameters.
  expect_reference(as.numeric(fit$n_parameters), 2 + 3 + 8,
                   reference = "independent: parameter count from the model definition",
                   covers = "multilpa()", quantity = "number of parameters with covariates",
                   tolerance = 0, measure = "absolute")
  expect_identical(nobs(fit), 20L)
  expect_equal(BIC(fit), fit$bic, tolerance = 1e-12)
})

test_that("covariate scores and observed information match numerical derivatives", {
  skip_without_reference("numDeriv")
  data <- review_ml_covariate_data()
  check <- function(fit, type, label) {
    x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
    theta <- latents:::.multilpa_cov_encode(fit)
    log_likelihood <- function(values) {
      pieces <- latents:::.multilpa_cov_decode(values, fit)
      latents:::.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                                          fit$profile_design, fit$group_design, pieces$beta,
                                          pieces$gamma, NULL, unname(fit$sampling_weights),
                                          fit$extra_data)$log_likelihood
    }
    withr::local_seed(9)
    moved <- theta + stats::rnorm(length(theta), sd = 0.1)
    inverse <- solve(-numDeriv::hessian(log_likelihood, theta))
    reference <- if (identical(type, "observed")) inverse else {
      scores <- latents:::.multilpa_cov_group_scores(theta, x, fit, NULL)
      inverse %*% (crossprod(scores) * nrow(scores) / (nrow(scores) - 1)) %*% inverse
    }
    cells <- expand.grid(row = seq_along(theta), column = seq_along(theta))
    list(score = list(
           actual = data.frame(model = label, parameter = seq_along(theta),
                               score = unname(colSums(latents:::.multilpa_cov_group_scores(moved, x, fit, NULL)))),
           expected = data.frame(model = label, parameter = seq_along(theta),
                                 score = numDeriv::grad(log_likelihood, moved))),
         covariance = list(
           actual = data.frame(model = label, cells,
                               covariance = as.vector(unname(latents:::.multilpa_cov_covariance(fit, data, 1e-4, type)))),
           expected = data.frame(model = label, cells, covariance = as.vector(reference))))
  }
  sides <- list(
    check(review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                   profile_covariates = "x", group_covariates = "gz",
                                   n_starts = 2, seed = 1)), "observed",
          "profile and group covariates, observed"),
    check(review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                   profile_covariates = "x", covariance_model = "full",
                                   n_starts = 2, seed = 1)), "observed",
          "profile covariate, full covariance, observed"),
    check(review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                   profile_covariates = "x", group_covariates = "gz",
                                   weights = "wt", n_starts = 2, seed = 1)), "robust",
          "profile and group covariates, weighted, robust"))
  stack <- function(part, side) do.call(rbind, lapply(sides, \(s) s[[part]][[side]]))
  note <- "Package side: the internal analytic group scores and covariance (.multilpa_cov_group_scores, .multilpa_cov_covariance) behind vcov() and parameter_inference() of multilpa() with covariates."
  # numDeriv's error is absolute near zero; a relative measure would judge
  # the rounding noise of a near-zero score. The original tolerance, 1e-6.
  expect_reference(stack("score", "actual"), stack("score", "expected"), reference = "numDeriv",
                   covers = "multilpa()", quantity = "covariate-model score off the optimum",
                   tolerance = 1e-6, measure = "scaled", key = c("model", "parameter"),
                   note = paste(note, "Reference: numDeriv::grad() of the package's likelihood."))
  # 1e-5, the original tolerance: the package's information is a central
  # difference of the analytic scores with step 1e-4, the reference the
  # inverse of numDeriv's Hessian (with the sandwich for robust). Relative:
  # the covariances are of order 1e-5 to 1, so a scaled measure would judge
  # most of them on an absolute 1e-5, looser than the original test.
  expect_reference(stack("covariance", "actual"), stack("covariance", "expected"),
                   reference = "numDeriv",
                   covers = c("multilpa()", "vcov.multilpa_covariates()", "parameter_inference.multilpa_covariates()"),
                   quantity = "covariate-model covariance of the estimates",
                   tolerance = 1e-5, measure = "relative", key = c("model", "row", "column"),
                   note = paste(note, "Reference: the inverse of minus numDeriv::hessian() of the package's likelihood; for robust, the sandwich with the analytic scores' cross-product, n / (n - 1)."))
})
