skip_on_cran()

# The order and blocks of the sequence view: ties in prevalence, noise.

test_that("profiles whose shares tie up to rounding are ranked by number", {
  # 3 + 4e-16 and 3 effective members are the same share; profile 1 comes
  # first whichever sum came out a bit larger.
  key <- .gg_profile_key(list(n_profiles = 2L,
                              effective_profile_counts = c(3, 3 + 4.4e-16)))
  expect_identical(key$rank, c(1L, 2L))
  key <- .gg_profile_key(list(n_profiles = 2L,
                              effective_profile_counts = c(3 + 4.4e-16, 3)))
  expect_identical(key$rank, c(1L, 2L))
  # A real difference still orders by share.
  key <- .gg_profile_key(list(n_profiles = 2L,
                              effective_profile_counts = c(2, 3)))
  expect_identical(key$rank, c(2L, 1L))
})

test_that("noise counts as a state when a group's most frequent one is found", {
  key <- data.frame(profile = c(1L, 2L, 0L), rank = c(1L, 2L, 3L))
  codes <- rbind(c(0L, 0L, 0L, 2L), c(1L, 1L, 0L, NA), c(0L, 0L, 0L, 0L))
  dominant <- .gg_sequence_dominant(codes, key)
  expect_identical(dominant$rank, c(3, 1, 3))
  expect_equal(dominant$share, c(3 / 4, 2 / 3, 1))
})

test_that("a one-class noise fit blocks its rows without an NA strip", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  n_groups <- 40L
  per <- 6L
  frame <- data.frame(g = rep(seq_len(n_groups), each = per),
                      t = rep(seq_len(per), n_groups))
  state <- unlist(lapply(seq_len(n_groups), \(i) {
    cumsum(stats::runif(per) > 0.7) %% 2 + 1
  }))
  frame$x <- stats::rnorm(nrow(frame), ifelse(state == 2, 2, -2))
  frame$y <- stats::rnorm(nrow(frame), ifelse(state == 2, 1, -1))
  # Two groups of pure outliers, classified as noise at every occasion.
  set.seed(9)
  outliers <- frame$g %in% c(5, 12)
  frame$x[outliers] <- stats::runif(sum(outliers), -40, 40)
  frame$y[outliers] <- stats::runif(sum(outliers), -40, 40)
  fit <- multilpa(frame, c("x", "y"), "g", n_profiles = 2, n_group_classes = 1,
                  time = "t", noise = TRUE, n_starts = 3, seed = 1)
  plot <- plot_sequences(fit)
  strips <- plot$data$class_label
  expect_false(anyNA(strips))
  expect_true("Mostly noise" %in% levels(strips))
  # The all-noise groups sit in the noise block.
  noise_rows <- unique(subset(plot$data, class_label == "Mostly noise")$row)
  expect_gte(length(noise_rows), 2L)
})

test_that("rows within a block hold less and less of their main profile", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  frame <- data.frame(g = rep(1:40, each = 6), t = rep(1:6, 40))
  state <- unlist(lapply(1:40, \(i) cumsum(stats::runif(6) > 0.7) %% 2 + 1))
  frame$x <- stats::rnorm(nrow(frame), ifelse(state == 2, 2, -2))
  frame$y <- stats::rnorm(nrow(frame), ifelse(state == 2, 1, -1))
  fit <- multilpa(frame, c("x", "y"), "g", n_profiles = 2, n_group_classes = 1,
                  time = "t", n_starts = 3, seed = 1)
  cells <- subset(plot_sequences(fit)$data, !empty)
  rows <- split(cells, cells$row)
  share <- vapply(rows, \(r) max(table(r$profile)) / nrow(r), numeric(1))
  block <- vapply(rows, \(r) as.character(r$class_label[1L]), character(1))
  ordered <- vapply(split(share, block), \(s) !is.unsorted(rev(s)), logical(1))
  expect_true(all(ordered))
})

test_that("an occasion no group answered keeps its own empty column", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  frame <- data.frame(g = rep(1:30, each = 4), t = rep(1:4, 30))
  state <- unlist(lapply(1:30, \(i) cumsum(stats::runif(4) > 0.6) %% 2 + 1))
  frame$x <- stats::rnorm(nrow(frame), ifelse(state == 2, 2, -2))
  frame$y <- stats::rnorm(nrow(frame), ifelse(state == 2, 1, -1))
  frame$x[frame$t == 2] <- NA
  frame$y[frame$t == 2] <- NA
  fit <- multilpa(frame, c("x", "y"), "g", n_profiles = 2, n_group_classes = 1,
                  time = "t", missing = "fiml", n_starts = 3, seed = 1)
  cells <- plot_sequences(fit)$data
  expect_setequal(unique(cells$column), 1:4)
  expect_identical(unique(subset(cells, empty)$column), 2L)
  # Every row is one unit's assignments, occasion by occasion.
  assigned <- get_results(fit, "assignments")
  assigned <- assigned[order(assigned$g, assigned$t), ]
  cells <- cells[order(cells$row, cells$column), ]
  drawn <- vapply(split(paste0(cells$column, ":",
                               ifelse(cells$empty, "-", cells$profile)),
                        cells$row), paste, character(1), collapse = " ")
  units <- vapply(split(paste0(assigned$t, ":",
                               ifelse(is.na(assigned$profile), "-", assigned$profile)),
                        assigned$g), paste, character(1), collapse = " ")
  expect_identical(sort(unname(drawn)), sort(unname(units)))
})

