# Random k-means starts on tightly clustered data can leave a cluster empty.
# Base R warns "empty cluster" for each such start; the fit already reports
# failed starts itself, so that warning is not passed on.

# More than 1.5 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

tight_clusters <- function() {
  set.seed(2026L)
  level <- rep(1:3, each = 72L)
  data.frame(
    student = rep(seq_len(36L), times = 6L),
    views = level * 6 + stats::rnorm(216L, sd = 0.08),
    sessions = level * 3 + stats::rnorm(216L, sd = 0.08),
    duration = level * 12 + stats::rnorm(216L, sd = 0.08)
  )
}

warning_messages <- function(expr) {
  messages <- character()
  withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  messages
}

test_that("base k-means' empty-cluster warning does not reach the caller", {
  messages <- warning_messages(
    lpa(tight_clusters(), c("views", "sessions", "duration"), n_profiles = 4,
        model = "EEI", n_starts = 25, seed = 123))
  expect_false(any(startsWith(messages, "empty cluster")))
})
