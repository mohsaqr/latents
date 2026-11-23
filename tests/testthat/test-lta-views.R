skip_on_cran()

# A general transition fit (here: ordinal indicators) draws the profile and
# sequence views of the profile models from its one measurement block.

.lta_view_data <- function(seed = 3L, n_groups = 40L, per = 6L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(n_groups), each = per),
                      t = rep(seq_len(per), n_groups))
  state <- unlist(lapply(seq_len(n_groups), \(i) {
    cumsum(stats::runif(per) > 0.7) %% 2 + 1
  }))
  rating <- function(shift) {
    as.integer(cut(stats::rnorm(length(state), ifelse(state == 2, shift, -shift)),
                   c(-Inf, -1, 0, 1, Inf)))
  }
  frame$a <- rating(1.2)
  frame$b <- rating(1.2)
  frame$c <- rating(0.8)
  frame
}

.lta_view_fit <- function(...) {
  lta(.lta_view_data(), c("a", "b", "c"), id = "g", time = "t", n_profiles = 2,
      ordinal = c("a", "b", "c"), n_starts = 3, seed = 1, ...)
}

test_that("plot_profiles() draws an ordinal lta's expected categories", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  expect_s3_class(fit, "multilpa_lta")
  plot <- plot_profiles(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "Profile means across 3 ordinal indicators")

  categories <- get_results(fit, "ordinal")
  categories$weighted <- as.numeric(categories$category) * categories$probability
  expected <- aggregate(weighted ~ profile + indicator, data = categories, FUN = sum)
  expected$profile <- as.integer(sub("profile_", "", expected$profile))
  drawn <- merge(plot$data, expected, by = c("profile", "indicator"))
  expect_identical(nrow(drawn), 6L)
  expect_equal(drawn$mean, drawn$weighted, tolerance = 1e-12)
})

test_that("plot_sequences() draws each unit's modal profile at each occasion", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  plot <- plot_sequences(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(nrow(plot$data), 240L)
  assigned <- get_results(fit, "assignments")
  expect_identical(as.vector(table(plot$data$profile)),
                   as.vector(table(assigned$profile)))
  # Each drawn row is one unit's sequence, cell by cell in occasion order:
  # the rows, read left to right, are the units' assigned sequences.
  cells <- plot$data[order(plot$data$row, plot$data$column), ]
  drawn <- vapply(split(cells$profile, cells$row), paste, character(1),
                  collapse = "")
  assigned <- assigned[order(assigned$g, assigned$t), ]
  units <- vapply(split(assigned$profile, assigned$g), paste, character(1),
                  collapse = "")
  expect_identical(sort(unname(drawn)), sort(unname(units)))
})

test_that("plot() of an lta fit checks its arguments and names unused ones", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  expect_error(plot(fit, what = "sequences", labels = "yes"),
               class = "latents_bad_argument")
  expect_error(plot(fit, what = "sequences", cell_labels = NA),
               class = "latents_bad_argument")
  expect_warning(plot(fit, what = "transitions", statistic = "median"),
                 regexp = "`statistic`", class = "latents_unused_argument")
  expect_warning(plot(fit, what = "profiles", cell_labels = FALSE),
                 class = "latents_unused_argument")
  expect_no_warning(plot(fit, what = "sequences", labels = FALSE,
                         cell_labels = FALSE))
})

test_that("a mover-stayer fit's sequences name the stayer class", {
  skip_if_not_installed("ggplot2")
  set.seed(5)
  n_groups <- 80L
  per <- 6L
  frame <- data.frame(g = rep(seq_len(n_groups), each = per),
                      t = rep(seq_len(per), n_groups))
  state <- unlist(lapply(seq_len(n_groups), \(i) {
    if (i <= n_groups / 2) rep(i %% 2 + 1, per) else
      cumsum(stats::runif(per) > 0.4) %% 2 + 1
  }))
  frame$a <- stats::rnorm(nrow(frame), ifelse(state == 2, 2, -2))
  frame$b <- stats::rnorm(nrow(frame), ifelse(state == 2, 1, -1))
  fit <- lta(frame, c("a", "b"), id = "g", time = "t", n_profiles = 2,
             mover_stayer = TRUE, n_starts = 3, seed = 1)
  expect_s3_class(fit, "multilpa_lta")
  plot <- plot_sequences(fit)
  expect_setequal(levels(plot$data$class_label), c("Class 1", "Stayers"))
})

test_that("plot(fit) still draws the transitions by default", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  expect_identical(plot(fit)$labels$title, "Transition probabilities")
  expect_identical(plot_transitions(fit)$labels$title, "Transition probabilities")
})

test_that("a profile argument the lta views do not take is refused", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  expect_error(plot_profiles(fit, scale = "standardized"),
               class = "latents_bad_argument")
})

test_that("occasion-specific measurement has no single profile view", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit(measurement = "occasion")
  expect_error(plot_profiles(fit), class = "latents_nothing_to_plot")
  expect_error(plot_sequences(fit), class = "latents_nothing_to_plot")
})

test_that("sequence rows are blocked and sorted by each group's most frequent profile", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  plot <- plot_sequences(fit)
  key <- .gg_profile_key(.lta_plot_view(fit))
  cells <- plot$data
  observed <- subset(cells, !empty)
  # Each row's most frequent profile, by display rank, never decreases down the rows.
  # Ties go to the more prevalent profile, the lower display rank.
  modal <- aggregate(profile ~ row, data = observed, FUN = \(p) {
    counts <- table(key$rank[match(p, key$profile)])
    min(as.integer(names(counts)[counts == max(counts)]))
  })
  expect_false(is.unsorted(modal$profile))
  expect_match(plot$labels$subtitle, "sorted by their most frequent profile")
  expect_true(all(grepl("^Mostly ", levels(cells$class_label))))
})

test_that(".gg_sequence_dominant() counts the modal profile and its share", {
  key <- data.frame(profile = 1:3, rank = c(2L, 1L, 3L))
  codes <- rbind(c(1L, 1L, 2L, NA), c(3L, 3L, 3L, 3L), c(NA, NA, NA, NA))
  dominant <- .gg_sequence_dominant(codes, key)
  expect_identical(dominant$rank, c(2, 3, NA))
  expect_equal(dominant$share, c(2 / 3, 1, NA))
})

test_that("plot() of an lta fit keeps the title in its old position", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  expect_identical(plot(fit, "My title")$labels$title, "My title")
  expect_identical(plot(fit, "My title", "Sub")$labels$subtitle, "Sub")
})

test_that("plot_all() draws what a general transition fit has", {
  skip_if_not_installed("ggplot2")
  fit <- .lta_view_fit()
  drawn <- plot_all(fit)
  expect_true(length(drawn) >= 3L)
  by_occasion <- .lta_view_fit(measurement = "occasion")
  drawn <- expect_message(plot_all(by_occasion), "profiles, sequences")
  expect_true(length(drawn) >= 1L)
})

