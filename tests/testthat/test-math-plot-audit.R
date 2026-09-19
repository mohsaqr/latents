test_that("single-profile response lines carry every indicator's probability", {
  data <- data.frame(group = rep(1:2, each = 4),
                     a = c(1, 1, 2, 2, 1, 1, 1, 2),
                     b = c(1, 2, 2, 2, 1, 2, 2, 2))
  fit <- multilpa(data, c("a", "b"), "group", 1, 1,
                  categorical = c("a", "b"), n_starts = 1)
  captured <- list()
  local_mocked_bindings(lines = function(x, y, ...) {
    captured[[length(captured) + 1L]] <<- list(x = x, y = y)
    invisible(NULL)
  }, .package = "graphics")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, what = "responses"))
  expect_length(captured, 1L)
  expect_equal(captured[[1L]]$x, 1:2)
  expect_equal(unname(captured[[1L]]$y), c(3 / 8, 6 / 8))
})

test_that("sequence images align profile rows and group classes with factor identifiers", {
  data <- data.frame(group = factor(rep(c("z", "a", "m"), each = 4),
                                    levels = c("m", "z", "a", "unused")),
                     time = rep(1:4, 3),
                     score = c(-2, -2, 2, 2, rep(-2, 4), rep(2, 4)))
  start <- list(means = matrix(c(-2, 2), 2), variances = matrix(.25, 2, 1),
                 profile_probabilities = rbind(c(.9, .1), c(.1, .9)),
                 group_probabilities = c(.5, .5))
  fit <- multilpa(data, "score", "group", 2, 2, n_starts = 1,
                  start = start, max_iter = 0, time = "time")
  captured <- NULL
  local_mocked_bindings(image = function(x, y, z, ...) {
    captured <<- list(x = x, y = y, z = z, extras = list(...))
    invisible(NULL)
  }, .package = "graphics")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, what = "sequences", palette = "black"))
  expect_equal(captured$x, 1:4)
  expect_equal(captured$y, 1:3)
  # Class 1 contains a then z; class 2 contains m. The unused factor level is
  # not an observed group and must not acquire an image row.
  expect_equal(captured$z, cbind(rep(1L, 4), c(1L, 1L, 2L, 2L), rep(2L, 4)))
  expect_identical(captured$extras$col, c("black", "black"))
})

test_that("one-group sequence images preserve dimensions and text time labels", {
  data <- data.frame(group = 1, time = c("early", "middle", "late"),
                     score = c(-1, .2, 2))
  fit <- multilpa(data, "score", "group", 1, 1, n_starts = 1, time = "time")
  captured <- NULL
  local_mocked_bindings(image = function(x, y, z, ...) {
    captured <<- list(x = x, y = y, z = z)
    invisible(NULL)
  }, .package = "graphics")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, what = "sequences"))
  expect_equal(captured$x, 1:3)
  expect_equal(captured$y, 1L)
  expect_equal(captured$z, matrix(1L, 3L, 1L))
})

test_that("sequence plots preserve groups with identical printed numeric identifiers", {
  data <- data.frame(group = rep(c(1, 1 + 1e-15), each = 3),
                     time = rep(1:3, 2), score = c(-2, -1, -3, 2, 1, 3))
  start <- list(means = matrix(c(-2, 2), 2), variances = matrix(1, 2, 1),
                profile_probabilities = matrix(c(.5, .5), 1), group_probabilities = 1)
  fit <- multilpa(data, "score", "group", 2, 1, n_starts = 1,
                  start = start, max_iter = 0, time = "time")
  # Covariate fits store native group values but do not have a group_ids field.
  fit$group_ids <- NULL
  captured <- NULL
  local_mocked_bindings(image = function(x, y, z, ...) {
    captured <<- z
    invisible(NULL)
  }, .package = "graphics")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, what = "sequences"))
  expect_equal(captured, cbind(rep(1L, 3), rep(2L, 3)))
})
