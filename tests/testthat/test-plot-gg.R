skip_if_not_installed("ggplot2")

gg_vars <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
gg_grid <- local({
  set.seed(1)
  enumerate_lpa(as.data.frame(scale(iris[gg_vars])), gg_vars,
                n_profiles = 2:4, model = c("EEE", "VVI"), n_starts = 3)
})
gg_fit <- candidate_fit(gg_grid, n_profiles = 4, model = "EEE")

test_that("every view of a fit and a grid builds without a warning", {
  fit_views <- c("profiles", "bars", "heatmap", "raincloud", "parallel",
                 "pairs", "probabilities", "sizes", "entropy", "posteriors",
                 "avepp")
  invisible(lapply(fit_views, \(view) expect_plot(plot(gg_fit, what = view))))
  expect_plot(plot(gg_grid))
  expect_plot(plot(gg_grid, what = "tree"))
})

test_that("the profile key orders by posterior share and matches counts", {
  key <- .gg_profile_key(gg_fit)
  counts <- get_results(gg_fit, what = "counts")
  counts <- counts[counts$level == "individuals", ]
  expect_equal(sum(key$share), 1)
  expect_identical(key$profile, seq_len(gg_fit$n_profiles))
  expect_equal(key$share[counts$class], counts$effective_proportion)
  # Largest first: rank 1 is the largest share.
  expect_false(is.unsorted(rev(key$share[order(key$rank)])))
  expect_identical(levels(key$label), as.character(key$label[order(key$rank)]))
})

test_that("bars start at zero and show the fitted means", {
  bars <- plot_layer(plot(gg_fit, what = "bars"), 2L)
  expect_true(all(abs(pmin(bars$ymin, bars$ymax)) < 1e-12 |
                    abs(pmax(bars$ymin, bars$ymax)) < 1e-12))
  expect_equal(sort(bars$y), sort(as.vector(gg_fit$means)))
})

test_that("the heatmap is the package's one standardization", {
  drawn <- plot(gg_fit, what = "heatmap")$data
  table <- get_results(gg_fit, what = "profiles", scale = "standardized")
  expect_equal(drawn$mean[order(drawn$profile, drawn$position)],
               table$mean[order(table$profile,
                                match(table$indicator, gg_vars))])
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

test_that("a diagonal structure's ellipses come from its variances", {
  vvi <- candidate_fit(gg_grid, n_profiles = 3, model = "VVI")
  expect_null(vvi$covariances)
  expect_equal(.gg_pair_covariance(vvi, 2L, 1L, 3L),
               diag(vvi$variances[2L, c(1L, 3L)]))
  expect_equal(.gg_pair_covariance(gg_fit, 2L, 1L, 3L),
               gg_fit$covariances[c(1L, 3L), c(1L, 3L), 2L])
})

test_that("label spreading keeps order and the minimum gap", {
  y <- c(1, 1.01, 1.02, 3)
  spread <- .gg_spread(y, gap = 0.5)
  expect_identical(order(spread), order(y))
  expect_true(all(diff(sort(spread)) >= 0.5 - 1e-12))
  expect_equal(mean(spread), mean(y))
})

test_that("tree bars are posterior shares and every row sums to one", {
  tree <- plot(gg_grid, what = "tree")
  bars <- plot_layer(tree, 2L)
  widths <- bars$xmax - bars$xmin
  bars$row <- round((bars$ymin + bars$ymax) / 2)
  row_total <- tapply(widths, list(bars$PANEL, bars$row), sum)
  expect_equal(as.vector(row_total), rep(1, length(row_total)),
               tolerance = 1e-8)
  fit_three <- candidate_fit(gg_grid, n_profiles = 3, model = "EEE")
  counts <- get_results(fit_three, what = "counts")
  shares <- counts$effective_proportion[counts$level == "individuals"]
  eee <- which(levels(tree$layers[[2L]]$data$panel) == "EEE")
  eee_three <- sort(widths[bars$PANEL == eee & abs(bars$row) == 3])
  expect_equal(eee_three, sort(shares), tolerance = 1e-8)
})

test_that("a tree needs two numbers of profiles for some model", {
  single <- enumerate_lpa(as.data.frame(scale(iris[gg_vars])), gg_vars,
                          n_profiles = 3, model = "EEE", n_starts = 2,
                          seed = 1)
  expect_error(plot(single, what = "tree"), class = "latents_nothing_to_plot")
})

test_that("the plots of several views print each one", {
  plots <- .gg_plots(list(sizes = plot(gg_fit, what = "sizes"),
                          avepp = plot(gg_fit, what = "avepp")))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE, after = FALSE)
  expect_invisible(print(plots))
  expect_identical(print(plots), plots)
})
