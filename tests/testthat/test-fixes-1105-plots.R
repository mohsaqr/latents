# Repairs found by reviewing the 2026-11-05 plot previews: blank legend keys,
# an empty group-class axis, collision-named person panels, weighted heatmap
# wording, plot_individuals(persons =) and the categorical boundary flag.

skip_on_cran()
skip_if_not_installed("ggplot2")

categorical_grid <- function(seed, exact) {
  set.seed(seed)
  n <- 250L
  class <- sample(1:2, n, TRUE)
  draw <- function(p) ifelse(runif(n) < p[class], "yes", "no")
  data <- data.frame(
    a = if (exact) ifelse(class == 1L, "yes", "no") else draw(c(.8, .25)),
    b = if (exact) ifelse(class == 1L, "up", "down") else draw(c(.75, .2)),
    c = draw(c(.85, .2)), d = draw(c(.9, .3)), e = draw(c(.75, .2)))
  quietly(enumerate_lca(data, names(data), n_classes = 2:3, n_starts = 2, seed = 17))
}

trust_statuses <- c("best replicated", "best found once", "replication unavailable",
                    "on a boundary", "not converged", "not fitted")

test_that("the trust view draws a legend key for every status", {
  trust <- expect_plot(plot_trust(categorical_grid(1106, exact = TRUE)))
  expect_true(isTRUE(trust$layers[[1L]]$show.legend))
  # Only "on a boundary" occurs here; the other five keep their key.
  expect_identical(as.character(unique(trust$data$status)), "on a boundary")
  expect_identical(ggplot2::get_guide_data(trust, "fill")$.label, trust_statuses)
})

test_that("a single-level grid has no group-class axis; a multilevel grid keeps it", {
  single <- categorical_grid(1106, exact = TRUE)
  views <- list(trust = plot_trust(single),
                heatmap = plot_selection(single, what = "heatmap"))
  lapply(views, \(view) {
    expect_null(expect_plot(view)$labels$y)
    expect_null(ggplot2::ggplot_build(view)$layout$panel_params[[1L]]$y$breaks)
  })
  multilevel <- quietly(enumerate_classes(
    subset(course_engagement, student <= 30), c("browse", "lectures"), "student",
    n_profiles = 1:2, n_group_classes = 1:2, model = "VVI", n_starts = 1, seed = 1))
  expect_identical(plot_trust(multilevel)$labels$y, "Number of group classes")
  expect_identical(plot_selection(multilevel, what = "heatmap")$labels$y,
                   "Number of group classes")
})

test_that("the group-class axis is dropped only when every candidate has one", {
  labels_of <- function(counts) {
    ggplot2::ggplot(data.frame(x = 1, y = "1")) + .latents_group_class_axis(counts)
  }
  expect_null(labels_of(c(1L, 1L, 1L))$labels$y)
  expect_identical(labels_of(c(1L, 2L))$labels$y, "Number of group classes")
  # One level that is not one group class still says how many there are.
  expect_identical(labels_of(2L)$labels$y, "Number of group classes")
})

test_that("the boundary flag separates interior from bounded categorical candidates", {
  grid <- categorical_grid(2211, exact = FALSE)
  candidates <- get_candidates(grid)
  # Invariant: the table's flag is the reporting predicate of each fit.
  expect_identical(candidates$boundary,
                   unname(vapply(grid$fits, .latents_boundary_active, logical(1))))
  expect_false(subset(candidates, n_profiles == 2L)$boundary)
  status <- as.character(plot_trust(grid)$data$status)
  expect_false(status[candidates$n_profiles == 2L] == "on a boundary")
  bounded <- get_candidates(categorical_grid(1106, exact = TRUE))
  expect_true(all(bounded$boundary))
})

test_that("the sensitivity plot keeps a key for statuses with no point", {
  set.seed(110514)
  class <- rep(1:2, each = 60L)
  data <- data.frame(a = rnorm(120, c(-1, 1)[class]), b = rnorm(120, c(1, -1)[class]))
  fit <- quietly(lpa(data, c("a", "b"), 2, n_starts = 1, seed = 9))
  drawn <- plot(quietly(sensitivity(fit, seeds = 1:2, n_starts = 1)))
  expect_true(isTRUE(drawn$layers[[1L]]$show.legend))
})

test_that("the weighted heatmap says weighted and uses weighted population moments", {
  set.seed(110517)
  class <- rep(1:2, each = 60L)
  data <- data.frame(a = rnorm(120, c(-2, 2)[class], .7),
                     b = rnorm(120, c(1, -1)[class], .8), w = c(1, 3)[class])
  data <- rbind(data, data.frame(a = 100, b = -100, w = 0))
  weighted <- quietly(lpa(data, c("a", "b"), 2, weights = "w", n_starts = 2, seed = 17))
  heatmap <- expect_plot(plot_heatmap(weighted))
  expect_match(heatmap$labels$subtitle, "weighted mean, in weighted standard deviations")
  # Against an independent implementation: stats::cov.wt(method = "ML") gives
  # the weighted mean and the weighted population covariance (denominator the
  # weight total). The zero-weight outlier enters with weight 0.
  moments <- stats::cov.wt(data[c("a", "b")], wt = data$w, method = "ML")
  centre <- moments$center
  spread <- sqrt(diag(moments$cov))
  raw <- as.data.frame(get_results(weighted, "profiles"))
  expected <- (raw$mean - centre[raw$indicator]) / spread[raw$indicator]
  standardized <- as.data.frame(get_results(weighted, "profiles", scale = "standardized"))
  expect_equal(standardized$mean, unname(expected), tolerance = 1e-10)
  tiles <- heatmap$data
  expect_equal(sort(tiles$mean), sort(unname(expected)), tolerance = 1e-10)
  unweighted <- quietly(lpa(subset(data, w > 0), c("a", "b"), 2, n_starts = 2, seed = 17))
  expect_match(plot_heatmap(unweighted)$labels$subtitle,
               "observed mean, in observed standard deviations")
})

test_that("person panels name a colliding identifier as the assignment table does", {
  set.seed(110515)
  id <- rep(1001:1040, each = 4L)
  state <- rep(sample(0:1, 40, TRUE), each = 4L)
  time <- rep(0:3, 40L)
  data <- data.frame(class = id, time = time,
    y = 3 * state + (.7 - 1.1 * state) * time +
      rep(rnorm(40, sd = .6), each = 4L) + rnorm(length(id), sd = .5))
  fit <- quietly(mixture_regression(y ~ time, data, 2, id = "class",
    class_level = "group", random = "intercept", n_starts = 2, seed = 32))
  chosen <- c(1001L, 1005L, 1012L)
  drawn <- expect_plot(plot_individuals(fit, persons = chosen))
  expect_setequal(unique(drawn$data$person), chosen)
  expect_true(all(startsWith(levels(drawn$data$person_label), "id_class 10")))
  expect_true("id_class" %in% names(as.data.frame(get_results(fit, "assignments"))))
  # The verb draws exactly what the fit's own view draws.
  own <- plot(fit, what = "individuals", persons = chosen)
  expect_identical(levels(drawn$data$person_label), levels(own$data$person_label))
  expect_length(unique(plot_individuals(fit)$data$person), 12L)
  expect_identical(unique(plot_individuals(fit, persons = "1001")$data$person), 1001L)
  expect_error(plot_individuals(fit, persons = 0), class = "latents_bad_argument")
  expect_error(plot_individuals(fit, persons = c(1L, 2L)), class = "latents_bad_argument")
  # An identifier with no collision keeps its own name.
  names(data)[names(data) == "class"] <- "student"
  plain <- quietly(mixture_regression(y ~ time, data, 2, id = "student",
    class_level = "group", random = "intercept", n_starts = 2, seed = 32))
  expect_true(all(startsWith(
    levels(plot_individuals(plain, persons = chosen)$data$person_label), "student 10")))
})
