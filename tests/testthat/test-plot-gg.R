skip_if_not_installed("ggplot2")

gg_vars <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
gg_grid <- local({
  set.seed(1)
  enumerate_lpa(as.data.frame(scale(iris[gg_vars])), gg_vars,
                n_profiles = 2:4, model = c("EEE", "VVI"), n_starts = 3)
})
gg_fit <- candidate_fit(gg_grid, n_profiles = 4, model = "EEE")

built_layer <- function(plot, layer) ggplot2::layer_data(plot, layer)

test_that("every view builds a ggplot without a warning", {
  fit_views <- c("profiles", "bars", "heatmap", "raincloud", "parallel",
                 "pairs", "sizes", "entropy", "posteriors", "avepp")
  lapply(fit_views, \(view) {
    plot <- .gg_plot(gg_fit, view)
    expect_s3_class(plot, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(plot))
  })
  lapply(c("enumeration", "tree"), \(view) {
    plot <- .gg_plot(gg_grid, view)
    expect_s3_class(plot, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(plot))
  })
})

test_that("the profile key orders by posterior share and matches counts", {
  key <- .gg_profile_key(gg_fit)
  counts <- get_results(gg_fit, what = "counts")
  counts <- counts[counts$level == "individuals", ]
  expect_equal(sum(key$share), 1)
  expect_false(is.unsorted(rev(key$share)))
  expect_setequal(key$profile, counts$class)
  expect_equal(key$share[match(counts$class, key$profile)],
               counts$effective_proportion)
})

test_that("bars start at zero and show the fitted means", {
  plot <- .gg_plot(gg_fit, "bars")
  bars <- built_layer(plot, 2L)
  expect_true(all(abs(pmin(bars$ymin, bars$ymax)) < 1e-12 |
                    abs(pmax(bars$ymin, bars$ymax)) < 1e-12))
  means <- get_results(gg_fit, what = "profiles")$mean
  expect_equal(sort(bars$y), sort(means))
})

test_that("the heatmap is the mixture's own standardized distance", {
  means <- get_results(gg_fit, what = "profiles")
  key <- .gg_profile_key(gg_fit)
  pi <- key$share[match(means$profile, key$profile)]
  rows <- means$indicator == "Petal.Length"
  overall <- sum(pi[rows] * means$mean[rows])
  spread <- sqrt(sum(pi[rows] * (means$variance[rows] + means$mean[rows]^2)) -
                   overall^2)
  expected <- (means$mean[rows] - overall) / spread
  plot <- .gg_plot(gg_fit, "heatmap")
  drawn <- plot$data$fill_value[plot$data$indicator == "Petal.Length"]
  expect_equal(sort(drawn), sort(expected))
})

test_that("case entropy is consistent with the reported relative entropy", {
  cases <- .gg_case_values(gg_fit)
  entropy <- get_results(gg_fit, what = "entropy")
  reported <- entropy$relative_entropy[entropy$level == "individuals"]
  expect_equal(1 - mean(cases$entropy), reported)
  expect_true(all(cases$certainty >= 1 / 4 - 1e-12 & cases$certainty <= 1))
})

test_that("an ellipse has the covariance it was drawn from", {
  sigma <- matrix(c(2, 0.6, 0.6, 1), 2L)
  ellipse <- .gg_ellipse(c(1, -1), sigma, n = 2000L)
  expect_equal(colMeans(ellipse), c(x = 1, y = -1), tolerance = 1e-2)
  # Points on the ellipse sit at the chi-square radius in Mahalanobis units.
  centred <- sweep(as.matrix(ellipse), 2L, c(1, -1))
  distance <- rowSums((centred %*% solve(sigma)) * centred)
  expect_equal(distance, rep(stats::qchisq(0.95, 2L), 2000L),
               tolerance = 1e-8)
})

test_that("label spreading keeps order and the minimum gap", {
  y <- c(1, 1.01, 1.02, 3)
  spread <- .gg_spread(y, gap = 0.5)
  expect_identical(order(spread), order(y))
  expect_true(all(diff(sort(spread)) >= 0.5 - 1e-12))
})

test_that("tree bars are posterior shares and every row sums to one", {
  plot <- .gg_plot(gg_grid, "tree")
  bars <- built_layer(plot, 2L)
  widths <- bars$xmax - bars$xmin
  bars$row <- round((bars$ymin + bars$ymax) / 2)
  row_total <- tapply(widths, list(bars$PANEL, bars$row), sum)
  expect_equal(as.vector(row_total), rep(1, length(row_total)),
               tolerance = 1e-8)
  fit_three <- candidate_fit(gg_grid, n_profiles = 3, model = "EEE")
  counts <- get_results(fit_three, what = "counts")
  shares <- counts$effective_proportion[counts$level == "individuals"]
  eee_three <- sort(widths[bars$PANEL == 1L & abs(bars$row) == 3])
  expect_equal(eee_three, sort(shares), tolerance = 1e-8)
})

test_that("views refuse the wrong input with a classed error", {
  expect_error(.gg_plot(gg_fit, "tree"))
  expect_error(.gg_plot(gg_grid, "profiles"))
  single <- enumerate_lpa(as.data.frame(scale(iris[gg_vars])), gg_vars,
                          n_profiles = 3, model = "EEE", n_starts = 2,
                          seed = 1)
  expect_error(.gg_plot(single, "tree"), class = "latents_bad_argument")
})
