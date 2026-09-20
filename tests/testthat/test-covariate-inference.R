.covariate_fit <- function(data, ...) {
  fit_covariates(data, c("y1", "y2"), "g", 2, 2, "z", "w",
                 n_starts = 10, seed = 812, tol = 1e-12, ...)
}

test_that("the analytic scores agree with a numerical gradient", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  x <- sweep(as.matrix(fixture$data[fit$vars]), 2L, fit$center, "-")
  theta <- multilpa:::.multilpa_cov_encode(fit)
  likelihood <- function(parameters) {
    pieces <- multilpa:::.multilpa_cov_decode(parameters, fit)
    multilpa:::.multilpa_cov_expectation(
      x, fit$group_index, pieces$parameters, fit$profile_design,
      fit$group_design, pieces$beta, pieces$gamma)$log_likelihood
  }
  analytic <- unname(colSums(multilpa:::.multilpa_cov_group_scores(theta, x, fit)))
  numerical <- vapply(seq_along(theta), function(j) {
    step <- 1e-5 * max(1, abs(theta[j]))
    up <- theta; up[j] <- up[j] + step
    down <- theta; down[j] <- down[j] - step
    (likelihood(up) - likelihood(down)) / (2 * step)
  }, numeric(1))

  # The formula, not the plumbing: every block of the score must match.
  expect_equal(analytic, numerical, tolerance = 1e-4)
  # and at the maximum the gradient is zero
  expect_true(max(abs(analytic)) < 1e-3)
})

test_that("standard errors reproduce a genuine Mplus covariate run", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  inference <- parameter_inference(fit, fixture$data)
  membership <- subset(inference, parameter == "coefficient")
  value <- function(term) membership$estimate[membership$term == term]
  error <- function(term) membership$standard_error[membership$term == term]

  # Mplus 9 covariates.out, MODEL RESULTS. Its group-class reference is the
  # opposite of this package's, so the group-level signs are mirrored.
  expect_equal(value("z"), -1.008, tolerance = 5e-4)
  expect_equal(error("z"), 0.117, tolerance = 5e-3)
  expect_equal(value("group_class_1"), 1.662, tolerance = 5e-4)
  expect_equal(error("group_class_1"), 0.152, tolerance = 5e-3)
  expect_equal(abs(value("w")), 0.740, tolerance = 5e-4)
  expect_equal(error("w"), 0.311, tolerance = 5e-3)
  expect_equal(abs(value("(Intercept)")), 0.214, tolerance = 5e-3)
  expect_equal(error("(Intercept)"), 0.278, tolerance = 5e-3)
})

test_that("the group-class contrast matches Mplus through vcov()", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  covariance <- vcov(fit, fixture$data)
  rows <- grep("^profile[.]coefficient[.]profile_1[.]group_class", rownames(covariance))
  contrast <- c(-1, 1)
  inference <- parameter_inference(fit, fixture$data)
  estimates <- inference$estimate[grepl("^group_class", inference$term) &
                                    inference$level == "profile"]

  # Mplus reports CW#1 ON CB#1 = -3.233 (S.E. 0.226); this package carries one
  # intercept per group class, whose difference is the same quantity.
  expect_equal(sum(contrast * estimates), -3.233, tolerance = 1e-3)
  expect_equal(sqrt(drop(contrast %*% covariance[rows, rows] %*% contrast)),
               0.226, tolerance = 5e-3)
})

test_that("the table is tidy and covers every level", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  inference <- parameter_inference(fit, fixture$data)

  expect_s3_class(inference, "data.frame")
  expect_named(inference, c("level", "outcome", "term", "parameter", "estimate",
                            "standard_error", "statistic", "p_value",
                            "p_adjusted", "conf_low", "conf_high"))
  expect_setequal(unique(inference$level), c("measurement", "profile", "group"))
  expect_setequal(unique(inference$parameter), c("mean", "variance", "coefficient"))
  expect_equal(nrow(inference), fit$n_parameters)
  expect_true(all(inference$standard_error > 0))
  expect_true(all(is.finite(inference$estimate)))
  # means are reported in input units, not the centred scale used internally
  expect_equal(inference$estimate[inference$parameter == "mean"],
               as.vector(t(fit$means)))
  # variances are positive and their intervals do not straddle zero downward
  variances <- subset(inference, parameter == "variance")
  expect_true(all(variances$estimate > 0))
  expect_true(all(is.na(variances$p_value)))
  expect_true(all(variances$conf_low > 0))
  # the interval is the estimate plus or minus the usual multiple
  membership <- subset(inference, parameter == "coefficient")
  expect_equal(membership$conf_high - membership$conf_low,
               2 * stats::qnorm(0.975) * membership$standard_error)
})

test_that("vcov is a named, symmetric matrix over the same parameters", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  covariance <- vcov(fit, fixture$data)

  expect_equal(dim(covariance), c(fit$n_parameters, fit$n_parameters))
  expect_equal(covariance, t(covariance))
  expect_true(all(diag(covariance) > 0))
  expect_identical(rownames(covariance), colnames(covariance))
  expect_identical(rownames(covariance), names(coef(fit)))
})

test_that("vcov agrees with the tidy table on whichever scale is asked for", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  inference <- parameter_inference(fit, fixture$data)
  natural <- vcov(fit, fixture$data)
  estimation <- vcov(fit, fixture$data, scale = "unconstrained")

  # The default is the scale the tidy table reports on, so the two agree
  # exactly; that is what makes confint() and vcov() mean the same thing.
  expect_equal(inference$standard_error, unname(sqrt(diag(natural))))
  # On the estimation scale variances are logs, so they differ from the tidy
  # table by exactly the delta-method factor and nothing else does.
  is_variance <- inference$parameter == "variance"
  expect_equal(inference$standard_error[is_variance],
               unname(sqrt(diag(estimation))[is_variance]) *
                 inference$estimate[is_variance])
  expect_equal(diag(natural)[!is_variance], diag(estimation)[!is_variance])
})

test_that("a named multiplicity correction is applied and nothing is corrected by stealth", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  plain <- parameter_inference(fit, fixture$data)
  corrected <- parameter_inference(fit, fixture$data, adjust = "BH")

  expect_identical(attr(plain, "adjust"), "none")
  expect_equal(plain$p_adjusted, plain$p_value)
  expect_identical(attr(corrected, "adjust"), "BH")
  expect_equal(corrected$p_value, plain$p_value)
  # The family is the tests the table reports, not its rows: a variance has no
  # test and must not inflate the correction.
  tested <- !is.na(plain$p_value)
  expect_equal(corrected$p_adjusted[tested],
               stats::p.adjust(plain$p_value[tested], method = "BH"))
  expect_true(all(is.na(corrected$p_adjusted[!tested])))
  expect_true(all(corrected$p_adjusted[tested] >= plain$p_value[tested]))
})

test_that("the robust sandwich runs and differs from the observed information", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  observed <- parameter_inference(fit, fixture$data)
  robust <- parameter_inference(fit, fixture$data, vcov_type = "robust")

  expect_equal(robust$estimate, observed$estimate)
  expect_true(all(robust$standard_error > 0))
  expect_false(isTRUE(all.equal(robust$standard_error, observed$standard_error)))
  expect_identical(attr(robust, "vcov_type"), "robust")
  expect_identical(attr(observed, "vcov_type"), "observed")
})

test_that("a broken contract is refused rather than answered", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)

  expect_error(parameter_inference(fit, fixture$data, level = 1),
               "`level` must be a single number")
  wrong <- fixture$data
  wrong$y1 <- wrong$y1 + 1
  expect_error(parameter_inference(fit, wrong),
               class = "multilpa_bad_inference_data")
  expect_error(parameter_inference(fit, fixture$data[c("y1", "z", "w", "g")]),
               class = "multilpa_bad_inference_data")
  stale <- fit
  stale$profile_design <- NULL
  expect_error(parameter_inference(stale, fixture$data),
               class = "multilpa_unsupported_inference")
})

test_that("confidence level widens the interval as asked", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  narrow <- parameter_inference(fit, fixture$data, level = 0.90)
  wide <- parameter_inference(fit, fixture$data, level = 0.99)

  expect_true(all(wide$conf_high - wide$conf_low >
                    narrow$conf_high - narrow$conf_low))
  expect_equal(attr(narrow, "level"), 0.90)
})
