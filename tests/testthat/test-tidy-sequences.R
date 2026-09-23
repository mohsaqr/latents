# The tidy contract of the sequence and transition verbs.
#
# These tests are about shape, not about the estimates: a result must carry
# every identifier it describes, name its columns so they can be used unquoted,
# say honestly what it summarizes, and let a caller restrict it with an
# argument rather than with brackets.

.tidy_panel <- function(n_groups = 8L, positions = 6L, seed = 11L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = positions)
  class_of_group <- rep(1:2, length.out = n_groups)
  profile <- 1L + as.integer(
    stats::runif(length(group)) > c(0.85, 0.15)[class_of_group[group]])
  data.frame(school = group, wave = rep(seq_len(positions), times = n_groups),
             a = stats::rnorm(length(group), c(-2.5, 2.5)[profile], 0.6),
             b = stats::rnorm(length(group), c(-2, 2)[profile], 0.7))
}

.tidy_fit <- function(data, ...) {
  multilpa(data, c("a", "b"), "school", n_profiles = 2, n_group_classes = 2,
           n_starts = 4, seed = 1, time = "wave", ...)
}

test_that("the wide form attributes every row without relying on row order", {
  data <- .tidy_panel()
  fit <- .tidy_fit(data)
  wide <- get_results(fit, "sequences", format = "wide")

  # The defect this guards: a wide form whose only link to a person is the
  # position of the row, with integer column names no one can write unquoted.
  expect_true("group" %in% names(wide))
  expect_false(anyNA(wide$group))
  expect_equal(anyDuplicated(wide$group), 0L)
  expect_identical(names(wide), make.names(names(wide)))
  expect_identical(names(wide), make.unique(names(wide)))

  # A permuted fit describes the same groups, whatever order they arrive in.
  set.seed(3)
  shuffled <- data[sample(nrow(data)), , drop = FALSE]
  row.names(shuffled) <- NULL
  permuted <- get_results(.tidy_fit(shuffled), "sequences", format = "wide")
  expect_identical(permuted$group, wide$group)
  expect_identical(names(permuted), names(wide))
})

test_that("two groups that print alike stay two groups in every shape", {
  data <- data.frame(school = rep(c(1, 1 + 1e-15), each = 3),
                     wave = rep(1:3, 2), a = c(-2, -1, -3, 2, 1, 3),
                     b = c(-2, -1, -3, 2, 1, 3))
  fit <- multilpa(data, c("a", "b"), "school", n_profiles = 1,
                  n_group_classes = 1, n_starts = 1, time = "wave")
  wide <- get_results(fit, "sequences", format = "wide")

  expect_equal(nrow(wide), 2L)
  # Row names coerce to character and would merge these two; a column does not.
  expect_identical(wide$group, unique(data$school))
  expect_false(isTRUE(all.equal(wide$group[1L], wide$group[2L],
                                tolerance = 0)))
  expect_equal(get_results(fit, "sequence_summary")$groups, 2L)
})

test_that("a gapped panel is summarized as gapped, not merely as short", {
  data <- .tidy_panel()
  # group 1 stops early; group 2 keeps going but skips a wave in the middle
  gapped <- data[!(data$school == 1L & data$wave > 4L) &
                   !(data$school == 2L & data$wave == 3L), , drop = FALSE]
  fit <- .tidy_fit(gapped)
  summary_table <- get_results(fit, "sequence_summary")
  classes <- get_results(fit, "sequences")

  expect_equal(sum(summary_table$observations), nrow(gapped))
  # exactly one group has a hole inside its own span, and it is not the
  # truncated one, which the old summary could not tell apart
  expect_equal(sum(summary_table$gaps), 1L)
  gapped_class <- unique(classes$group_class[classes$group == 2L])
  expect_equal(summary_table$gaps[summary_table$group_class == gapped_class], 1L)
  expect_equal(sum(summary_table$complete), fit$n_groups - 2L)
  expect_false(anyNA(summary_table))
})

test_that("a group class with no groups is still a row, with NA lengths", {
  data <- .tidy_panel()
  fit <- .tidy_fit(data)
  # Reassign every group to the first class, leaving the second empty.
  fit$group_classes <- rep(1L, fit$n_groups)
  summary_table <- get_results(fit, "sequence_summary")

  expect_equal(nrow(summary_table), fit$n_group_classes)
  expect_equal(summary_table$groups, c(8L, 0L))
  expect_true(is.na(summary_table$mean_length[2L]))
  expect_true(is.na(summary_table$shortest[2L]))
  expect_true(is.na(summary_table$longest[2L]))
  expect_equal(summary_table$complete[2L], 0L)
  expect_equal(summary_table$gaps[2L], 0L)
})

test_that("sequence_summary refuses a fit without an ordering", {
  data <- .tidy_panel()
  bare <- multilpa(data, c("a", "b"), "school", n_profiles = 2,
                   n_group_classes = 2, n_starts = 2, seed = 1)
  expect_error(get_results(bare, "sequence_summary"), class = "latents_no_time")
  expect_error(get_results(bare, "sequences", format = "wide"), class = "latents_no_time")
})

test_that("the occasion columns are named after the ordering column", {
  data <- .tidy_panel(n_groups = 4L, positions = 3L)
  names(data)[names(data) == "wave"] <- "occasion number"
  fit <- multilpa(data, c("a", "b"), "school", n_profiles = 2,
                  n_group_classes = 2, n_starts = 2, seed = 1,
                  time = "occasion number")
  wide <- get_results(fit, "sequences", format = "wide")

  # A non-syntactic time column still yields syntactic occasion names.
  expect_identical(names(wide),
                   c("group", "group_class", "occasion.number_1",
                     "occasion.number_2", "occasion.number_3"))
})
