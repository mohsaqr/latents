# A posterior can round a hair past 1 (seen on a two-level fit of VaSStra's
# engagement data: 2 of 1136 rows). Reported uncertainty and plotted certainty
# stay inside [0, 1]; the stored estimates are left as computed.

skip_on_cran()

rounded_fit <- function() {
  vars <- c("browse", "lectures")
  fit <- lpa(course_engagement, vars, n_profiles = 2, n_starts = 1, seed = 1)
  fit$subject_posteriors[1L, ] <- c(1 + .Machine$double.eps, 0)
  fit$subject_profiles[1L] <- 1L
  fit
}

test_that("uncertainty is never negative when a posterior rounds past 1", {
  fit <- rounded_fit()
  uncertainty <- get_results(fit, "assignments")$uncertainty
  expect_gte(min(uncertainty), 0)
  expect_identical(uncertainty[1L], 0)
  # The stored posterior is not altered.
  expect_gt(fit$subject_posteriors[1L, 1L], 1)
})

test_that("per-case certainty and entropy stay on their [0, 1] axes", {
  cases <- .gg_case_values(rounded_fit())
  expect_true(all(cases$certainty >= 0 & cases$certainty <= 1))
  expect_true(all(cases$entropy >= 0 & cases$entropy <= 1))
})

test_that("the classification plots keep every case", {
  skip_if_not_installed("ggplot2")
  fit <- rounded_fit()
  for_view <- \(what) ggplot2::ggplot_build(plot(fit, what = what))
  expect_no_warning(for_view("posteriors"))
  expect_no_warning(for_view("entropy"))
})
