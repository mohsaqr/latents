# The numbers a view draws, read from the plot object: a plot must show the
# fit's values, in the right place.

# The assigned profile of each group (rows) at each position (columns), as the
# sequence view lays them out.
sequence_codes <- function(plot) {
  cells <- plot$data
  codes <- matrix(NA_integer_, max(cells$row), max(cells$column))
  codes[cbind(cells$row, cells$column)] <- cells$profile
  codes
}

test_that("single-profile response lines carry every indicator's probability", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(group = rep(1:2, each = 4),
                     a = c(1, 1, 2, 2, 1, 1, 1, 2),
                     b = c(1, 2, 2, 2, 1, 2, 2, 2))
  fit <- multilpa(data, c("a", "b"), "group", 1, 1,
                  categorical = c("a", "b"), n_starts = 1)
  drawn <- expect_plot(plot(fit, what = "responses"))$data
  expect_equal(drawn$position, 1:2)
  expect_equal(drawn$mean, c(3 / 8, 6 / 8))
})

test_that("sequence images align profile rows and group classes with factor identifiers", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(group = factor(rep(c("z", "a", "m"), each = 4),
                                    levels = c("m", "z", "a", "unused")),
                     time = rep(1:4, 3),
                     score = c(-2, -2, 2, 2, rep(-2, 4), rep(2, 4)))
  start <- list(means = matrix(c(-2, 2), 2), variances = matrix(.25, 2, 1),
                 profile_probabilities = rbind(c(.9, .1), c(.1, .9)),
                 group_probabilities = c(.5, .5))
  fit <- multilpa(data, "score", "group", 2, 2, n_starts = 1,
                  start = start, max_iter = 0, time = "time")
  drawn <- expect_plot(plot(fit, what = "sequences"))
  # Groups sort by their DECLARED factor level order (m, z, a), not
  # alphabetically: class 1 holds z then a, class 2 holds m. The unused factor
  # level is not an observed group and must not acquire a row.
  expect_equal(sequence_codes(drawn),
               rbind(c(1L, 1L, 2L, 2L), rep(1L, 4), rep(2L, 4)))
})

test_that("one-group sequence images preserve dimensions and text time labels", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(group = 1, time = c("early", "middle", "late"),
                     score = c(-1, .2, 2))
  fit <- multilpa(data, "score", "group", 1, 1, n_starts = 1, time = "time")
  drawn <- expect_plot(plot(fit, what = "sequences"))
  expect_equal(sequence_codes(drawn), matrix(1L, 1L, 3L))
  labels <- ggplot2::ggplot_build(drawn)$layout$panel_params[[1L]]$x$get_labels()
  expect_setequal(labels, c("early", "middle", "late"))
})

test_that("sequence plots preserve groups with identical printed numeric identifiers", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(group = rep(c(1, 1 + 1e-15), each = 3),
                     time = rep(1:3, 2), score = c(-2, -1, -3, 2, 1, 3))
  start <- list(means = matrix(c(-2, 2), 2), variances = matrix(1, 2, 1),
                profile_probabilities = matrix(c(.5, .5), 1), group_probabilities = 1)
  fit <- multilpa(data, "score", "group", 2, 1, n_starts = 1,
                  start = start, max_iter = 0, time = "time")
  # Covariate fits store native group values but do not have a group_ids field.
  fit$group_ids <- NULL
  drawn <- expect_plot(plot(fit, what = "sequences"))
  expect_equal(sequence_codes(drawn), rbind(rep(1L, 3), rep(2L, 3)))
})
