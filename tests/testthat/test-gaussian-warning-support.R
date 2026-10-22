#' A small thin Gaussian profile with covariate and categorical columns
#' @return An 80-row data frame with three clearly separated Gaussian rows.
gaussian_warning_support_data <- function() {
  set.seed(73)
  data <- as.data.frame(matrix(rnorm(80L * 5L), 80L, 5L))
  names(data) <- letters[seq_len(5L)]
  data[seq_len(3L), ] <- data[seq_len(3L), ] + 10
  data$z <- rnorm(80L)
  data$categ <- c(rep("rare", 3L), rep("common", 77L))
  data
}

#' Capture Gaussian-support warnings from a public fit
#' @param expr A public fitting expression.
#' @return The fit and its emitted thin-profile warning conditions.
gaussian_warning_support_capture <- function(expr) {
  warnings <- list()
  fit <- withCallingHandlers(expr, warning = function(condition) {
    if (inherits(condition, "latents_thin_profiles")) {
      warnings[[length(warnings) + 1L]] <<- condition
    }
    if (inherits(condition, .multilpa_expected_warnings)) invokeRestart("muffleWarning")
  })
  list(fit = fit, thin = warnings)
}

test_that("retained empty rows cannot hide ordinary or covariate Gaussian warnings", {
  data <- gaussian_warning_support_data()
  empty <- data[rep(1L, 300L), ]
  empty[, letters[seq_len(5L)]] <- NA_real_
  # Covariates remain observed even when no measurement is present.
  empty$z <- rep(c(-1, 1), 150L)
  invisible(lapply(c(FALSE, TRUE), function(covariate) {
    fit <- function(frame) gaussian_warning_support_capture(lpa(
      frame, letters[seq_len(5L)], 2L,
      profile_covariates = if (covariate) "z" else character(),
      n_starts = 1L, seed = 1L, max_iter = 1000L,
      missing = "fiml", drop_empty = FALSE))
    original <- fit(data)
    appended <- fit(rbind(data, empty))
    expect_true(original$fit$converged)
    expect_true(appended$fit$converged)
    expect_length(original$thin, 1L)
    expect_length(appended$thin, 1L)
    expect_equal(min(original$fit$gaussian_diagnostics$gaussian_mass), 3, tolerance = 1e-7)
    expect_equal(min(appended$fit$gaussian_diagnostics$gaussian_mass), 3, tolerance = 1e-7)
    expect_gt(min(appended$fit$effective_profile_counts), 10)
    expect_equal(appended$fit$gaussian_diagnostics$gaussian_mass,
                 unname(colSums(appended$fit$subject_posteriors[seq_len(80L), ])))
  }))
})

test_that("categorical-only rows cannot hide ordinary or covariate Gaussian warnings", {
  data <- gaussian_warning_support_data()
  categorical_only <- data[rep(1L, 300L), ]
  categorical_only[, letters[seq_len(5L)]] <- NA_real_
  categorical_only$z <- rep(c(-1, 1), 150L)
  invisible(lapply(c(FALSE, TRUE), function(covariate) {
    captured <- gaussian_warning_support_capture(lpa(
      rbind(data, categorical_only), c(letters[seq_len(5L)], "categ"), 2L,
      categorical = "categ", profile_covariates = if (covariate) "z" else character(),
      n_starts = 1L, seed = 1L, max_iter = 1000L, missing = "fiml"))
    fit <- captured$fit
    expect_true(fit$converged)
    expect_length(captured$thin, 1L)
    expect_gt(min(fit$effective_profile_counts), 10)
    expect_equal(min(fit$gaussian_diagnostics$gaussian_mass), 3, tolerance = 1e-7)
    expect_equal(fit$gaussian_diagnostics$gaussian_mass,
                 unname(colSums(fit$subject_posteriors[seq_len(80L), ])))
    expect_true(all(fit$gaussian_diagnostics$positive_rows <= 80L))
  }))
})
