.sequence_fixture <- function(n_groups = 12L, positions = 8L, seed = 4L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = positions)
  time <- rep(seq_len(positions), times = n_groups)
  class_of_group <- rep(1:2, length.out = n_groups)
  profile <- 1L + as.integer(
    stats::runif(length(group)) > c(0.85, 0.15)[class_of_group[group]])
  data.frame(school = group, wave = time,
             a = stats::rnorm(length(group), c(-2.5, 2.5)[profile], 0.6),
             b = stats::rnorm(length(group), c(-2, 2)[profile], 0.7))
}

.sequence_fit <- function(data, ...) {
  multilpa(data, c("a", "b"), "school", n_profiles = 2, n_group_classes = 2,
           n_starts = 4, seed = 1, ...)
}

test_that("a fit carries the ordering only when asked, and refuses it otherwise", {
  data <- .sequence_fixture()
  bare <- .sequence_fit(data)
  timed <- .sequence_fit(data, time = "wave")

  expect_null(bare$time)
  expect_null(bare$time_values)
  expect_identical(timed$time, "wave")
  expect_identical(timed$time_values, data$wave)
  expect_error(sequences(bare), class = "multilpa_no_time")
  expect_error(sequence_summary(bare), class = "multilpa_no_time")
  # The ordering is metadata: carrying it must not move the estimates.
  expect_equal(bare$log_likelihood, timed$log_likelihood)
  expect_equal(bare$means, timed$means)
})

test_that("the long form is one tidy row per observation, in sequence order", {
  data <- .sequence_fixture()
  fit <- .sequence_fit(data, time = "wave")
  long <- sequences(fit)

  expect_s3_class(long, "data.frame")
  expect_named(long, c("group", "group_class", "time", "profile"))
  expect_equal(nrow(long), nrow(data))
  expect_false(anyNA(long))
  expect_true(all(long$profile %in% seq_len(fit$n_profiles)))
  expect_true(all(long$group_class %in% seq_len(fit$n_group_classes)))
  # ordered by group, then by time within group
  expect_false(is.unsorted(long$group))
  expect_true(all(vapply(split(long$time, long$group), is.unsorted, logical(1)) == FALSE))
  # one class per group, never two
  expect_true(all(vapply(split(long$group_class, long$group),
                         function(v) length(unique(v)), integer(1)) == 1L))
})

test_that("the wide form is one row per group and one column per position", {
  data <- .sequence_fixture()
  fit <- .sequence_fit(data, time = "wave")
  wide <- sequences(fit, format = "wide")
  long <- sequences(fit)

  expect_equal(dim(wide), c(fit$n_groups, length(unique(data$wave))))
  expect_named(wide, as.character(sort(unique(data$wave))))
  expect_setequal(row.names(wide), as.character(unique(long$group)))
  expect_true(all(vapply(wide, is.factor, logical(1))))
  expect_identical(levels(wide[[1L]]), as.character(seq_len(fit$n_profiles)))
  # the two shapes must agree cell by cell
  first <- as.character(long$group[1L])
  expect_equal(as.integer(unlist(wide[first, ])),
               long$profile[long$group == long$group[1L]])
})

test_that("an unbalanced group leaves NA at the positions it never reached", {
  data <- .sequence_fixture()
  trimmed <- data[!(data$school == 1L & data$wave > 5L), , drop = FALSE]
  fit <- .sequence_fit(trimmed, time = "wave")
  wide <- sequences(fit, format = "wide")

  expect_equal(nrow(sequences(fit)), nrow(trimmed))
  expect_equal(sum(is.na(wide)), 3L)
  expect_true(all(is.na(unlist(wide["1", c("6", "7", "8")]))))
})

test_that("the summary counts each class's groups and their lengths", {
  data <- .sequence_fixture()
  trimmed <- data[!(data$school == 1L & data$wave > 5L), , drop = FALSE]
  fit <- .sequence_fit(trimmed, time = "wave")
  summary_table <- sequence_summary(fit)
  long <- sequences(fit)

  expect_named(summary_table, c("group_class", "groups", "observations",
                                "mean_length", "median_length", "shortest",
                                "complete"))
  expect_equal(nrow(summary_table), fit$n_group_classes)
  expect_equal(sum(summary_table$groups), fit$n_groups)
  expect_equal(sum(summary_table$observations), nrow(long))
  expect_equal(min(summary_table$shortest), 5L)
  expect_equal(sum(summary_table$complete), fit$n_groups - 1L)
})

test_that("a broken ordering is refused rather than silently reshaped", {
  data <- .sequence_fixture()
  expect_error(.sequence_fit(data, time = "absent"),
               "`time` must name a column of `data`")
  missing_time <- data
  missing_time$wave[3L] <- NA
  expect_error(.sequence_fit(missing_time, time = "wave"),
               class = "multilpa_bad_time")
  repeated <- data
  repeated$wave[2L] <- repeated$wave[1L]
  expect_error(.sequence_fit(repeated, time = "wave"),
               class = "multilpa_bad_time")
})

test_that("the sequence panel draws and respects the no-ordering contract", {
  data <- .sequence_fixture()
  fit <- .sequence_fit(data, time = "wave")
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)

  grDevices::png(file, width = 900, height = 650)
  expect_silent(plot(fit, what = "sequences"))
  grDevices::dev.off()
  expect_gt(file.size(file), 1000)
  expect_error(plot(.sequence_fit(data), what = "sequences"),
               class = "multilpa_no_time")
})

test_that("the ordering is imposed, not inherited from the input row order", {
  data <- .sequence_fixture()
  set.seed(99)
  shuffled <- data[sample(nrow(data)), , drop = FALSE]
  row.names(shuffled) <- NULL
  fit <- .sequence_fit(shuffled, time = "wave")
  long <- sequences(fit)

  # The input is deliberately out of order, so a stable sort on group alone
  # would leave time scrambled within each group.
  expect_true(is.unsorted(shuffled$wave))
  expect_false(is.unsorted(long$group))
  expect_true(all(vapply(split(long$time, long$group),
                         function(v) identical(v, sort(v)), logical(1))))
  # and the wide form must put each observation in its own column
  wide <- sequences(fit, format = "wide")
  first <- as.character(long$group[1L])
  expect_equal(as.integer(unlist(wide[first, ])),
               long$profile[long$group == long$group[1L]])
})
