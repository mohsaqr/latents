# Multiple-imputation pooling (`pool_imputations()`): Rubin's rules against
# mice's own implementation, label alignment across imputations, and the
# relabelling of a covariate fit, whose membership logits change reference
# category when the labels move.

skip_on_cran()

.pool_fixture <- function(n_groups = 40L, size = 8L, seed = 5L) {
  set.seed(seed)
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size), z = stats::rnorm(n))
  frame$w <- stats::ave(stats::rnorm(n), frame$g)
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z))
  frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  frame$y2 <- stats::rnorm(n, c(0, 2)[state])
  frame
}

# Stand-in imputations: the same data with a few covariate values redrawn,
# which is all the pooling step can see of an imputation model.
.pool_imputed <- function(frame, m = 3L, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(m), function(index) {
    copy <- frame
    copy$z[seq_len(20L)] <- stats::rnorm(20L)
    copy
  })
}

.pool_quietly <- function(expr) {
  quietly(expr, classes = c(.multilpa_expected_warnings, "latents_extreme_coefficients"))
}

test_that("Rubin's rules agree with mice::pool.scalar for every parameter", {
  skip_if_not_installed("mice")
  pooled <- .pool_quietly(pool_imputations(
    .pool_imputed(.pool_fixture()), c("y1", "y2"), "g", 2L,
    n_group_classes = 1L, profile_covariates = "z", n_starts = 1, seed = 1))
  estimates <- as.data.frame(pooled)
  per_imputation <- get_results(pooled, "imputations")
  key <- function(table) paste(table$level, table$outcome, table$term, table$parameter)
  invisible(lapply(seq_len(nrow(estimates)), function(row) {
    rows <- per_imputation[key(per_imputation) == key(estimates)[row], ]
    reference <- mice::pool.scalar(rows$estimate, rows$standard_error^2, n = Inf)
    expect_equal(estimates$estimate[row], reference$qbar, tolerance = 1e-12)
    expect_equal(estimates$standard_error[row]^2, reference$t, tolerance = 1e-12)
    expect_equal(estimates$df[row], reference$df, tolerance = 1e-10)
    expect_equal(estimates$fmi[row], reference$fmi, tolerance = 1e-10)
  }))
})

test_that("identical imputations pool to the single fit, with no missing information", {
  frame <- .pool_fixture()
  pooled <- .pool_quietly(pool_imputations(
    list(frame, frame, frame), c("y1", "y2"), "g", 2L, n_group_classes = 1L,
    profile_covariates = "z", n_starts = 1, seed = 1))
  single <- parameter_inference(.pool_quietly(multilpa(
    frame, c("y1", "y2"), "g", 2L, 1L, profile_covariates = "z",
    n_starts = 1, seed = 1)))
  estimates <- as.data.frame(pooled)
  expect_equal(estimates$estimate, single$estimate)
  expect_equal(estimates$standard_error, single$standard_error)
  expect_true(all(estimates$between == 0 & is.infinite(estimates$df) &
                    estimates$fmi == 0))
})

test_that("relabelling a covariate fit leaves its likelihood and inference intact", {
  frame <- .pool_fixture(n_groups = 60L)
  # Three profiles, two group classes, slopes by group class and a group
  # covariate: every block the relabelling has to move is present.
  fit <- .pool_quietly(multilpa(frame, c("y1", "y2"), "g", 3L, 2L,
                                profile_covariates = "z", group_covariates = "w",
                                profile_slopes = "group_class", n_starts = 1,
                                seed = 1, max_iter = 300))
  moved <- .multilpa_cov_permute(fit, c(3L, 1L, 2L), c(2L, 1L))
  x <- sweep(as.matrix(frame[c("y1", "y2")]), 2L, fit$center, "-")
  likelihood <- function(object) {
    pieces <- .multilpa_cov_decode(.multilpa_cov_encode(object), object)
    .multilpa_cov_expectation(x, object$group_index, pieces$parameters,
                              object$profile_design, object$group_design,
                              pieces$beta, pieces$gamma, NULL)$log_likelihood
  }
  expect_equal(likelihood(moved), fit$log_likelihood, tolerance = 1e-10)
  # The relabelled posteriors are the original's columns, reordered.
  expect_equal(unname(moved$subject_posteriors),
               unname(fit$subject_posteriors[, c(3L, 1L, 2L)]))
  # Moving back restores every coefficient exactly.
  back <- .multilpa_cov_permute(moved, c(2L, 3L, 1L), c(2L, 1L))
  expect_equal(back$profile_coefficients, fit$profile_coefficients, tolerance = 1e-12)
  expect_equal(back$group_coefficients, fit$group_coefficients, tolerance = 1e-12)
  expect_equal(back$means, fit$means)
})

test_that("alignment undoes a scrambled imputation before pooling", {
  frame <- .pool_fixture()
  fit <- .pool_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 1L,
                                profile_covariates = "z", n_starts = 1, seed = 1))
  scrambled <- .multilpa_cov_permute(fit, c(2L, 1L), 1L)
  expect_false(isTRUE(all.equal(scrambled$profile_coefficients,
                                fit$profile_coefficients)))
  alignment <- .multilpa_pool_align(scrambled, fit)
  expect_identical(alignment$profile_order, c(2L, 1L))
  expect_equal(alignment$fit$profile_coefficients, fit$profile_coefficients,
               tolerance = 1e-12)
  expect_equal(parameter_inference(alignment$fit)$standard_error,
               parameter_inference(fit)$standard_error, tolerance = 1e-8)
})

test_that("imputations that come back relabelled are aligned inside the verb", {
  skip_if_not_installed("mice")
  set.seed(5)
  n <- 480L
  frame <- data.frame(g = rep(seq_len(60L), each = 8L), z = stats::rnorm(n))
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z))
  frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  frame$y2 <- stats::rnorm(n, c(0, 2)[state])
  frame$z[stats::runif(n) < stats::plogis(-2 + 0.8 * frame$y1)] <- NA
  imputed <- mice::mice(frame, m = 5, printFlag = FALSE, seed = 1)
  pooled <- .pool_quietly(pool_imputations(imputed, c("y1", "y2"), "g", 2L,
                                           n_group_classes = 1L,
                                           profile_covariates = "z",
                                           n_starts = 2, seed = 1))
  # On these data one imputation's EM labels its profiles the other way round
  # (measured: imputation 5). Pooled without alignment, its slope enters with
  # the opposite sign and the between-imputation variance swamps the rest.
  expect_gte(summary(pooled)$relabelled, 1L)
  slope <- subset(as.data.frame(pooled), term == "z")
  expect_gt(abs(slope$estimate), 1)
  expect_lt(slope$fmi, 0.5)
})

test_that("covariate-free fits are pooled too", {
  frame <- .pool_fixture()
  imputed <- lapply(1:2, function(index) {
    copy <- frame
    copy$y1[index * 3L] <- copy$y1[index * 3L] + 0.5
    copy
  })
  pooled <- .pool_quietly(pool_imputations(imputed, c("y1", "y2"), "g", 2L,
                                           n_group_classes = 1L, n_starts = 1,
                                           seed = 1))
  expect_identical(pooled$model, "multilpa")
  expect_true(all(is.finite(as.data.frame(pooled)$standard_error)))
  expect_identical(nrow(get_results(pooled, "fits")), 2L)
})

test_that("the methods return what they document", {
  pooled <- .pool_quietly(pool_imputations(
    .pool_imputed(.pool_fixture()), c("y1", "y2"), "g", 2L,
    n_group_classes = 1L, profile_covariates = "z", n_starts = 1, seed = 1))
  expect_s3_class(summary(pooled), "data.frame")
  expect_identical(nrow(summary(pooled)), 1L)
  expect_identical(dim(vcov(pooled)),
                   rep(nrow(as.data.frame(pooled)), 2L))
  expect_true(all(c("df", "riv", "fmi", "within", "between", "p_adjusted") %in%
                    names(as.data.frame(pooled))))
  expect_output(print(pooled), "Pooled over 3 imputations")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(pooled), pooled)
  expect_error(plot(pooled, level = "group"), class = "latents_bad_argument")
})

test_that("bad input and a failing imputation are refused by class", {
  frame <- .pool_fixture()
  call_with <- function(imputed) {
    pool_imputations(imputed, c("y1", "y2"), "g", 2L, n_group_classes = 1L,
                     profile_covariates = "z", n_starts = 1, seed = 1)
  }
  expect_error(call_with(frame), class = "latents_bad_data")
  expect_error(call_with(list(frame)), class = "latents_bad_data")
  expect_error(call_with(list(frame, frame[-1L, ])), class = "latents_bad_data")
  broken <- frame
  broken$z[1L] <- NA
  # Rubin's rules need every imputation, so one failure stops the pooling and
  # is named rather than dropped.
  expect_error(.pool_quietly(call_with(list(frame, broken))),
               class = "latents_pooling_failed")
  expect_error(.pool_quietly(call_with(list(frame, broken))), "Imputation 2")
})
