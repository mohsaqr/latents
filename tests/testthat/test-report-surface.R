# The reporting surface gathers other verbs. These tests pin the two ways that
# can go wrong without a word being said: an argument that is documented as
# forwarded but is dropped on the way, and a view that is promised for every fit
# but only exists for one family.

# One fit, built once, for the tests that only read it. Seeded, so the tables
# and the panels are the same on every run.
.surface_data <- function() {
  set.seed(41)
  n <- 144L
  high <- rep(c(FALSE, TRUE), each = n / 2L)
  data.frame(school = rep(seq_len(12L), each = 12L),
             a = ifelse(high, 4, -4) + stats::rnorm(n, sd = 0.3),
             b = ifelse(high, 4, -4) + stats::rnorm(n, sd = 0.3),
             z = stats::rnorm(n))
}

.surface_fit <- function() {
  multilpa(.surface_data(), c("a", "b"), "school", n_profiles = 2L,
           n_group_classes = 1L, n_starts = 2L, seed = 1)
}

test_that("diagnostics() forwards `by` instead of dropping it", {
  skip_on_cran()
  fit <- .surface_fit()
  pooled <- get_data(diagnostics(fit, by = "overall"), "residuals")
  # The point of the test: the table has to be the pooled one. Before the fix
  # `by` was discarded and this returned the per-profile residuals, which read
  # as pooled ones and are not.
  expect_identical(pooled, get_data(fit, "residuals", by = "overall"))
  expect_identical(unique(pooled$profile), "overall")
  expect_false(any(grepl("profile_", pooled$profile, fixed = TRUE)))
  # The default is unchanged, and is still the per-profile table.
  by_profile <- get_data(diagnostics(fit), "residuals")
  expect_identical(by_profile, get_data(fit, "residuals"))
  expect_false(identical(pooled, by_profile))
})

test_that("report() forwards `by` through to the residual it prints", {
  skip_on_cran()
  fit <- .surface_fit()
  printed <- utils::capture.output(report(fit, plots = FALSE, by = "overall"))
  residual_line <- grep("Largest residual", printed, value = TRUE)
  expect_length(residual_line, 1L)
  expect_true(grepl("overall", residual_line, fixed = TRUE))
  expect_false(grepl("profile_", residual_line, fixed = TRUE))
})

test_that("an argument the reporting surface cannot forward is refused", {
  skip_on_cran()
  fit <- .surface_fit()
  expect_error(diagnostics(fit, nonsense = 1),
               class = "multilpa_bad_argument")
  expect_error(report(fit, plots = FALSE, nonsense = 1),
               class = "multilpa_bad_argument")
  # The refusal names the argument, so the caller knows which one was refused.
  expect_error(diagnostics(fit, nonsense = 1), regexp = "nonsense")
  # A value outside the contract of an argument that *is* forwarded still
  # refuses rather than falling back to the default.
  expect_error(diagnostics(fit, by = "nonsense"))
})

test_that("plot(diagnostics()) draws for a covariate fit", {
  skip_on_cran()
  fit <- multilpa(.surface_data(), c("a", "b"), "school",
                        n_profiles = 2L, n_group_classes = 1L,
                        profile_covariates = "z", n_starts = 1L, seed = 1)
  draw({
    expect_invisible(plot(diagnostics(fit)))
    expect_invisible(plot(fit, what = "entropy"))
    expect_invisible(plot(fit, what = "posteriors"))
  })
  expect_true(all(c("entropy", "posteriors") %in% .multilpa_supported_views(fit)))
})

test_that("plot(diagnostics()) draws for a random-intercept fit", {
  skip_on_cran()
  fit <- fit_random_intercept(.surface_data(), c("a", "b"), "school",
                              n_profiles = 2L, n_starts = 2L, seed = 1)
  draw(expect_invisible(plot(diagnostics(fit))))
})

test_that("a single-profile fit refuses the case diagnostics by class", {
  skip_on_cran()
  # One profile leaves the two panels nothing to separate. The refusal is the
  # contract; before the fix the drawing code failed on a zero-width axis with
  # "`breaks` must be increasing", which reads as an internal fault.
  fits <- list(
    multilpa(.surface_data(), c("a", "b"), "school", n_profiles = 1L,
             n_group_classes = 1L, n_starts = 2L, seed = 1),
    fit_random_intercept(.surface_data(), c("a", "b"), "school",
                         n_profiles = 1L, n_starts = 1L, seed = 1))
  draw(invisible(lapply(fits, function(fit) {
    expect_error(plot(diagnostics(fit)), class = "multilpa_nothing_to_plot")
    # The refusal names what can be drawn instead, so it is a signpost and not
    # just a wall.
    expect_error(plot(diagnostics(fit)), regexp = "profiles")
    expect_false(any(c("entropy", "posteriors") %in%
                       .multilpa_supported_views(fit)))
  })))
})

test_that("report() draws every view it offers, for every family", {
  skip_on_cran()
  data <- .surface_data()
  fits <- list(
    covariates = multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                                n_group_classes = 1L, profile_covariates = "z",
                                n_starts = 1L, seed = 1),
    random_intercept = fit_random_intercept(data, c("a", "b"), "school",
                                            n_profiles = 2L, n_starts = 2L,
                                            seed = 1),
    one_profile = multilpa(data, c("a", "b"), "school", n_profiles = 1L,
                           n_group_classes = 1L, n_starts = 2L, seed = 1))
  draw(invisible(lapply(fits, function(fit) {
    expect_identical(utils::capture.output(report(fit, plots = TRUE))[1L],
                     utils::capture.output(print(summary(fit)))[1L])
  })))
})
