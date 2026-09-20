# Two ways a descriptive summary used to be quietly wrong: an NA stratifier
# poisoned every stratum's row selector, and the class a model assigned was
# written into the frame being described, over the top of whatever column
# already had that name.

.separated_data <- function(seed = 41L, n_groups = 12L, per = 12L) {
  set.seed(seed)
  n <- n_groups * per
  half <- rep(c(-4, 4), each = n / 2L)
  data.frame(g = rep(seq_len(n_groups), each = per),
             profile = half + stats::rnorm(n, sd = 0.3),
             group_class = half + stats::rnorm(n, sd = 0.3),
             b = half + stats::rnorm(n, sd = 0.3))
}

test_that("a missing stratifier neither invents missingness nor loses a value", {
  # The review's example verbatim: `y` is complete, and every one of its four
  # values is observed in exactly one stratum.
  described <- descriptives(data.frame(by = c("A", "A", "B", NA), y = 1:4),
                            vars = "y", by = "by")

  expect_identical(nrow(described), 3L)
  expect_identical(described$by, c("A", "B", NA))
  # `y` has no missing value anywhere, so no stratum may report one.
  expect_identical(described$n_missing, c(0L, 0L, 0L))
  expect_identical(described$n, c(2L, 1L, 1L))
  # The fourth observation is reported, not lost: it is the missing stratum's
  # only value.
  expect_equal(described$mean, c(1.5, 3, 4))
  expect_equal(described$min, c(1, 3, 4))
  expect_equal(described$max, c(2, 3, 4))
})

test_that("strata account for every input row, whatever the stratifier's type", {
  set.seed(7)
  data <- data.frame(
    numeric_by = c(rep(c(0.1, 0.2), each = 4L), rep(NA_real_, 2L)),
    factor_by = factor(c(rep(c("low", "high"), each = 4L), NA, NA),
                       levels = c("low", "high", "never_seen")),
    y = stats::rnorm(10),
    z = c(stats::rnorm(8), NA, NA))

  lapply(c("numeric_by", "factor_by"), function(stratifier) {
    described <- descriptives(data, vars = c("y", "z"), by = stratifier)
    accounted <- vapply(split(described$n + described$n_missing,
                              described$variable), sum, integer(1))
    # One row of the input is counted once for each variable, in exactly one
    # stratum, and the two missing values of `z` are still the only missing
    # values reported anywhere.
    expect_identical(unname(accounted), c(10L, 10L))
    expect_identical(sum(subset(described, variable == "z")$n_missing), 2L)
    expect_identical(sum(subset(described, variable == "y")$n_missing), 0L)
    expect_true(anyNA(described[[stratifier]]))
  })
  # The stratifier keeps its own type rather than being coerced to its label.
  expect_type(descriptives(data, vars = "y", by = "numeric_by")$numeric_by,
              "double")
  expect_s3_class(descriptives(data, vars = "y", by = "factor_by")$factor_by,
                  "factor")
})

test_that("a complete stratifier is summarized exactly as before", {
  described <- descriptives(data.frame(by = c("A", "A", "B", "B"), y = 1:4),
                            vars = "y", by = "by")
  expect_identical(described$by, c("A", "B"))
  expect_identical(described$n, c(2L, 2L))
  expect_equal(described$mean, c(1.5, 3.5))
})

test_that("an indicator named `profile` is summarized as itself", {
  data <- .separated_data()
  fit <- multilpa(data, c("profile", "b"), "g", n_profiles = 2,
                  n_group_classes = 1, n_starts = 2, seed = 1)
  described <- descriptives(fit, by = "profile")
  indicator <- subset(described, variable == "profile")
  truth <- vapply(split(data$profile, fit$subject_profiles), mean, numeric(1))

  # The class means of the caller's own `profile` column, not the class labels.
  expect_equal(indicator$mean, unname(truth))
  expect_true(all(indicator$sd > 0))
  # The stratum column still carries the assignment it splits on.
  expect_identical(sort(unique(described$profile)),
                   sort(unique(fit$subject_profiles)))
  expect_identical(sum(indicator$n + indicator$n_missing), nrow(data))
})

test_that("an indicator named `group_class` is summarized as itself", {
  data <- .separated_data()
  fit <- multilpa(data, c("group_class", "b"), "g", n_profiles = 2,
                  n_group_classes = 2, n_starts = 4, seed = 1)
  described <- descriptives(fit, by = "group_class")
  indicator <- subset(described, variable == "group_class")
  assigned <- fit$group_classes[fit$group_index]
  truth <- vapply(split(data$group_class, assigned), mean, numeric(1))

  expect_equal(indicator$mean, unname(truth))
  expect_true(all(indicator$sd > 0))
  expect_identical(sum(indicator$n + indicator$n_missing), nrow(data))
})

test_that("a stratifier that cannot be told apart from the summary is refused", {
  # `mean` is a column of the summary itself, so a stratum column of that name
  # would be unaddressable in the result.
  expect_error(descriptives(data.frame(mean = c("a", "b"), y = 1:2),
                            vars = "y", by = "mean"),
               class = "multilpa_bad_data")
})
