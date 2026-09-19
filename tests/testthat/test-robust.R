robust_fixture_data <- function(seed = 11L, n_groups = 30L, per_group = 8L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per_group)
  group_class <- rep(c(1L, 2L), length.out = n_groups)
  prevalence <- rbind(c(0.85, 0.15), c(0.2, 0.8))
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = prevalence[group_class[j], ])
  }, integer(1))
  means <- rbind(c(-2, -2), c(2, 2))
  y <- t(vapply(profile, function(k) stats::rnorm(2L, means[k, ], 1), numeric(2)))
  data.frame(g = group, a = y[, 1L], b = y[, 2L])
}

centered_matrix <- function(fit, dat) {
  x <- as.matrix(dat[, fit$indicators, drop = FALSE])
  sweep(x, 2L, colMeans(x, na.rm = TRUE), "-")
}

centered_theta <- function(fit, dat) {
  x <- as.matrix(dat[, fit$indicators, drop = FALSE])
  centered <- fit
  centered$means <- sweep(fit$means, 2L, colMeans(x, na.rm = TRUE), "-")
  .multilpa_coefficients(centered, "unconstrained")
}

test_that("per-group scores sum to the aggregate gradient in every model family", {
  dat <- robust_fixture_data()
  incomplete <- dat
  set.seed(2)
  incomplete$a[sample(nrow(dat), 25L)] <- NA
  incomplete$b[sample(nrow(dat), 18L)] <- NA
  specifications <- list(
    list(label = "diagonal varying", data = dat, arguments = list()),
    list(label = "diagonal equal", data = dat,
         arguments = list(variance_model = "equal")),
    list(label = "full varying", data = dat,
         arguments = list(covariance_model = "full")),
    list(label = "full equal", data = dat,
         arguments = list(covariance_model = "full", variance_model = "equal")),
    list(label = "diagonal fiml", data = incomplete,
         arguments = list(missing = "fiml")),
    list(label = "full fiml varying", data = incomplete,
         arguments = list(covariance_model = "full", missing = "fiml")),
    list(label = "full fiml equal", data = incomplete,
         arguments = list(covariance_model = "full", variance_model = "equal",
                          missing = "fiml")))
  invisible(lapply(specifications, function(specification) {
    fit <- do.call(multilpa, c(list(data = specification$data,
      indicators = c("a", "b"), group = "g", n_profiles = 2,
      n_group_classes = 2, n_starts = 5, seed = 5), specification$arguments))
    theta <- centered_theta(fit, specification$data)
    x <- centered_matrix(fit, specification$data)
    aggregate_gradient <- -.multilpa_score(theta, x, fit)
    scores <- .multilpa_group_scores(theta, x, fit)
    expect_identical(dim(scores), c(fit$n_groups, length(theta)),
                     info = specification$label)
    expect_equal(unname(colSums(scores)), unname(aggregate_gradient),
                 tolerance = 1e-8, info = specification$label)
  }))
})

test_that("per-group scores are localized to their own group", {
  dat <- robust_fixture_data()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  theta <- centered_theta(fit, dat)
  # Perturb one group's data only. Centering would spread that change across
  # every group, so the comparison uses the uncentered indicators directly.
  baseline_x <- as.matrix(dat[, c("a", "b")])
  shifted_x <- baseline_x
  shifted_x[dat$g == 1L, "a"] <- shifted_x[dat$g == 1L, "a"] + 0.5
  baseline <- .multilpa_group_scores(theta, baseline_x, fit)
  moved <- .multilpa_group_scores(theta, shifted_x, fit)
  expect_false(isTRUE(all.equal(baseline[1L, ], moved[1L, ])))
  expect_equal(baseline[-1L, ], moved[-1L, ])
})

test_that("robust standard errors reproduce genuine Mplus MLR results", {
  reference <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-robust.rds"))
  dat <- reference$data
  fit <- multilpa(dat, c("y1", "y2"), "clus", 2, 2, variance_model = "varying",
                    start = starting_values(reference), n_starts = 1,
                    max_iter = 5000, tol = 1e-13)
  information <- parameter_inference(fit, dat, vcov_type = "robust")
  expect_identical(attr(information, "vcov_type"), "robust")
  q <- fit$n_parameters
  n_measurement <- q - 3L
  # Mplus reports a group logit, a within-class intercept at cb = 2, and the
  # cb = 1 minus cb = 2 contrast, with variances rather than log variances.
  ordering <- c(1, 2, 5, 6, 3, 4, 7, 8)
  derivatives <- c(rep(1, 4), as.vector(t(fit$variances)))
  jacobian <- matrix(0, q, q)
  jacobian[cbind(seq_len(n_measurement), ordering)] <- derivatives[ordering]
  jacobian[q - 2L, q] <- 1
  jacobian[q - 1L, q - 1L] <- 1
  jacobian[q, c(q - 2L, q - 1L)] <- c(1, -1)
  standard_errors <- sqrt(diag(jacobian %*% attr(information, "covariance_unconstrained") %*%
                                 t(jacobian)))
  expect_lt(max(abs(standard_errors - reference$mplus_standard_errors)), 1e-6)
  expect_lt(abs(attr(information, "scaling_correction") - reference$mplus_scaling), 1e-6)
  indices <- information_criteria(fit)
  # A criterion that has no sample-size convention -- AIC counts parameters, not
  # units -- carries NA there rather than a sentinel level, so it is selected
  # by NA and not by a magic string.
  select <- function(criterion, convention) {
    wanted <- if (is.na(convention)) is.na(indices$convention) else
      !is.na(indices$convention) & indices$convention == convention
    indices$value[indices$criterion == criterion & wanted]
  }
  expect_lt(abs(select("aic", NA_character_) - reference$mplus_aic), 1e-3)
  expect_lt(abs(select("bic", "individuals") - reference$mplus_bic), 1e-3)
  expect_lt(abs(select("sabic", "individuals") - reference$mplus_sabic), 1e-3)
})

test_that("robust and observed covariances differ only through the score product", {
  dat <- robust_fixture_data()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  observed <- parameter_inference(fit, dat)
  robust <- parameter_inference(fit, dat, vcov_type = "robust")
  expect_identical(attr(observed, "vcov_type"), "observed")
  # These are diagnostics of the fit and travel as attributes, so they are
  # checked as attributes; reading them off the frame with `$` would silently
  # compare NULL with NULL and pass whatever the code did.
  expect_null(attr(observed, "group_scores"))
  expect_true(is.matrix(attr(robust, "group_scores")))
  expect_true(is.na(attr(observed, "scaling_correction")))
  expect_equal(attr(observed, "hessian"), attr(robust, "hessian"))
  expect_equal(observed$estimate, robust$estimate)
  # The sandwich is symmetric and positive semidefinite in free coordinates.
  free <- attr(robust, "covariance_unconstrained")
  expect_equal(free, t(free))
  expect_gt(min(eigen(free, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_true(all(robust$standard_error >= 0))
  expect_gt(attr(robust, "scaling_correction"), 0)
  # The data are generated from the fitted family, so the sandwich should be
  # close to the observed information; a large ratio would signal a bug.
  ratio <- robust$standard_error / observed$standard_error
  expect_true(all(ratio > 0.5 & ratio < 2))
  expect_equal(unname(vcov(fit, data = dat, vcov_type = "robust")),
               unname(attr(robust, "covariance")))
  expect_false(isTRUE(all.equal(unname(attr(observed, "covariance")),
                                unname(attr(robust, "covariance")))))
})

test_that("robust inference refuses fits with too few groups", {
  set.seed(4)
  # Four indicators and three groups: identified enough to fit, but far too few
  # independent units for the cross-product matrix to have full rank.
  small <- data.frame(g = rep(seq_len(3L), each = 20L),
                      a = stats::rnorm(60L), b = stats::rnorm(60L))
  fit <- multilpa(small, c("a", "b"), "g", 2, 1, n_starts = 4, seed = 1)
  expect_gt(fit$n_parameters, fit$n_groups)
  expect_error(parameter_inference(fit, small, vcov_type = "robust"),
               class = "multilpa_too_few_groups")
})

test_that("the cross-product matrix rejects non-finite scores", {
  expect_error(.multilpa_cross_product(matrix(c(1, NA, 2, 3), 2L, 2L)),
               class = "multilpa_bad_scores")
  expect_equal(.multilpa_cross_product(matrix(c(1, 2), 2L, 1L)),
               matrix(5, 1L, 1L))
})

test_that("robust standard errors are invariant to indicator location", {
  dat <- robust_fixture_data()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  robust <- parameter_inference(fit, dat, vcov_type = "robust")
  shifted <- dat
  shifted$a <- shifted$a + 100
  shifted_fit <- multilpa(shifted, c("a", "b"), "g", 2, 2,
                            start = local({
                              start <- starting_values(fit)
                              start$means[, 1L] <- start$means[, 1L] + 100
                              start
                            }), n_starts = 1, max_iter = 5000, tol = 1e-13)
  shifted_robust <- parameter_inference(shifted_fit, shifted, vcov_type = "robust")
  expect_equal(unname(shifted_robust$standard_error),
               unname(robust$standard_error), tolerance = 1e-6)
  expect_equal(attr(shifted_robust, "scaling_correction"), attr(robust, "scaling_correction"),
               tolerance = 1e-6)
})
