skip_on_cran()

# A fit whose indicators are all ordinal (or all categorical) has no profile
# means: its primary table is its response table, and the profile views draw
# its category probabilities.

.ordinal_view_data <- function(seed = 3L, n_groups = 60L, per = 6L) {
  set.seed(seed)
  n <- n_groups * per
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  rating <- function(shift) {
    latent <- stats::rnorm(n, ifelse(truth == 2L, shift, -shift))
    as.integer(cut(latent, c(-Inf, -1, 0, 1, Inf)))
  }
  data.frame(g = rep(seq_len(n_groups), each = per),
             a = rating(1.2), b = rating(1.2), c = rating(0.8))
}

.ordinal_view_fit <- function() {
  multilpa(.ordinal_view_data(), c("a", "b", "c"), "g", n_profiles = 2,
           n_group_classes = 1, ordinal = c("a", "b", "c"), n_starts = 4,
           seed = 1)
}

test_that("an all-ordinal fit's primary table is its ordinal table", {
  fit <- .ordinal_view_fit()
  ordinal <- get_results(fit, "ordinal")
  expect_identical(get_results(fit), ordinal)
  expect_identical(as.data.frame(fit), ordinal)
  expect_gt(nrow(ordinal), 0L)
})

test_that("an all-categorical fit's primary table is its response table", {
  fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports",
                            "walking"), n_classes = 2, n_starts = 2, seed = 1)
  expect_identical(get_results(fit), get_results(fit, "responses"))
  expect_gt(nrow(get_results(fit)), 0L)
})

test_that("a fit with continuous indicators keeps the profile table", {
  fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                  id = "student", n_profiles = 2, n_group_classes = 1,
                  n_starts = 2, seed = 1)
  expect_identical(get_results(fit), get_results(fit, "profiles"))
})

test_that("plot_profiles() draws each ordinal profile's expected category", {
  skip_if_not_installed("ggplot2")
  fit <- .ordinal_view_fit()
  plot <- plot_profiles(fit)
  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "Profile means across 3 ordinal indicators")

  categories <- get_results(fit, "ordinal")
  categories$weighted <- as.numeric(categories$category) * categories$probability
  expected <- aggregate(weighted ~ profile + indicator, data = categories,
                        FUN = sum)
  drawn <- merge(plot$data, expected,
                 by.x = c("profile", "indicator"), by.y = c("profile", "indicator"))
  expect_identical(nrow(drawn), 6L)
  expect_equal(drawn$mean, drawn$weighted, tolerance = 1e-12)
  # Expected categories lie within the rating scale.
  expect_true(all(drawn$mean >= 1 & drawn$mean <= 4))
})

test_that("plot_responses() draws ordinal category probabilities", {
  skip_if_not_installed("ggplot2")
  fit <- .ordinal_view_fit()
  plot <- plot_responses(fit, category = "first")
  expect_s3_class(plot, "ggplot")
  expect_identical(nrow(plot$data), 6L)
  first <- subset(get_results(fit, "ordinal"), category == "1")
  drawn <- data.frame(profile = plot$data$profile,
                      indicator = c("a", "b", "c")[plot$data$position],
                      drawn = plot$data$mean)
  matched <- merge(drawn, first, by = c("profile", "indicator"))
  expect_identical(nrow(matched), 6L)
  # Each profile's probability of the lowest answer on each indicator, in
  # its own place: a swap of profiles or indicators fails here.
  expect_equal(matched$drawn, matched$probability, tolerance = 1e-12)
  expect_identical(plot$labels$subtitle, "2 profiles, 3 ordinal indicators")
  heatmap <- plot_heatmap(fit)
  expect_s3_class(heatmap, "ggplot")
  # Every tile carries its profile's probability of that category.
  tiles <- data.frame(profile = heatmap$data$profile,
                      indicator = as.character(heatmap$data$indicator),
                      category = as.character(heatmap$data$category_label),
                      drawn = heatmap$data$probability)
  table <- get_results(fit, "ordinal")
  matched <- merge(tiles, table, by = c("profile", "indicator", "category"))
  expect_identical(nrow(matched), nrow(table))
  expect_equal(matched$drawn, matched$probability, tolerance = 1e-12)
})

test_that("a fit without categorical or ordinal indicators still refuses responses", {
  skip_if_not_installed("ggplot2")
  fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                  id = "student", n_profiles = 2, n_group_classes = 1,
                  n_starts = 2, seed = 1)
  expect_error(plot_responses(fit), class = "latents_no_categorical")
})
