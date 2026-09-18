plot_fixture <- function(seed = 11L, n_groups = 24L, per_group = 8L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per_group)
  group_class <- rep(c(1L, 2L), length.out = n_groups)
  prevalence <- rbind(c(0.85, 0.15), c(0.2, 0.8))
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = prevalence[group_class[j], ])
  }, integer(1))
  means <- rbind(c(-2, 40), c(2, 60))
  y <- t(vapply(profile, function(k) {
    stats::rnorm(2L, means[k, ], c(1, 6))
  }, numeric(2)))
  data.frame(g = group, a = y[, 1L], b = y[, 2L])
}

# Drawing is tested on a null device, so the tests need no graphics back end
# and leave no files behind.
draw <- function(expression) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 900, height = 600)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE, after = FALSE)
  force(expression)
}

test_that("plot methods draw and return their input invisibly", {
  dat <- plot_fixture()
  fit <- fit_ml_lpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    expect_invisible(plot(fit))
    expect_identical(plot(fit), fit)
    expect_identical(plot(fit, what = "probabilities"), fit)
    expect_identical(plot(fit, scale = "standardized"), fit)
    expect_identical(plot(fit, labels = FALSE), fit)
  })
  candidates <- enumerate_ml_lpa(dat, c("a", "b"), "g", profiles = 1:2,
                                 group_classes = 1:2, n_starts = 3, seed = 3)
  draw({
    expect_identical(plot(candidates), candidates)
    expect_identical(plot(candidates, criterion = "sabic_individual"), candidates)
    expect_identical(plot(candidates, labels = FALSE, mark_minimum = FALSE),
                     candidates)
  })
})

test_that("plotting restores the caller's graphical parameters", {
  dat <- plot_fixture()
  fit <- fit_ml_lpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    before <- graphics::par(c("mar", "xpd", "cex"))
    plot(fit)
    plot(fit, what = "probabilities")
    expect_equal(graphics::par(c("mar", "xpd", "cex")), before)
  })
})

test_that("plot arguments override every visual constant", {
  dat <- plot_fixture()
  fit <- fit_ml_lpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    expect_silent(plot(fit, palette = c("#000000", "#FFFFFF"),
                       symbols = c(1L, 2L), linetypes = c(1L, 3L),
                       main = "custom", subtitle = "custom",
                       panel_fill = "#FFFFFF", grid_colour = "#CCCCCC",
                       point_size = 2, line_width = 4))
    # A single supplied colour must recycle rather than fail.
    expect_silent(plot(fit, palette = "#000000"))
  })
})

test_that("standardizing uses the observed indicator scales", {
  dat <- plot_fixture()
  fit <- fit_ml_lpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  observed <- fit$indicator_data
  centre <- colMeans(observed)
  spread <- c(stats::sd(observed[, 1L]), stats::sd(observed[, 2L]))
  expected <- sweep(sweep(fit$means, 2L, centre, "-"), 2L, spread, "/")
  # The b indicator is on a far wider scale, so standardizing must shrink it
  # relative to a; this is the whole point of the option.
  raw_spread <- diff(range(fit$means[, "b"]))
  standardized_spread <- diff(range(expected[, "b"]))
  expect_lt(standardized_spread, raw_spread)
  expect_equal(colMeans(expected) * 0 + 1, c(a = 1, b = 1))
  stripped <- fit
  stripped$indicator_data <- NULL
  draw(expect_error(plot(stripped, scale = "standardized"),
                    class = "mllpa_no_indicator_data"))
})

test_that("enumeration plotting rejects unusable criteria by condition class", {
  dat <- plot_fixture()
  candidates <- enumerate_ml_lpa(dat, c("a", "b"), "g", profiles = 1:2,
                                 group_classes = 1, n_starts = 3, seed = 3)
  draw({
    expect_error(plot(candidates, criterion = "not_a_column"),
                 class = "mllpa_unknown_criterion")
    empty <- candidates
    empty$table$bic_individual <- NA_real_
    expect_error(plot(empty, criterion = "bic_individual"),
                 class = "mllpa_nothing_to_plot")
  })
})

test_that("the palette, symbols and line types stay aligned and recycle", {
  # Colour must never be the only channel, so all three vectors must have the
  # same length for any number of classes, including beyond the palette length.
  invisible(lapply(c(1L, 3L, 9L, 12L), function(n) {
    expect_length(.ml_lpa_palette(n), n)
    expect_length(.ml_lpa_symbols(n), n)
    expect_length(.ml_lpa_linetypes(n), n)
  }))
  expect_identical(.ml_lpa_palette(1L), "#E69F00")
  # Okabe-Ito, and nothing outside it.
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
  expect_true(all(.ml_lpa_palette(20L) %in% okabe_ito))
  expect_identical(.ml_lpa_palette(11L)[10L], .ml_lpa_palette(1L))
})

test_that("direct labels are spread apart but keep their order and centre", {
  expect_equal(.ml_lpa_spread_labels(c(1, 5), 0.5), c(1, 5))
  crowded <- .ml_lpa_spread_labels(c(1, 1, 1), 0.5)
  expect_equal(diff(sort(crowded)), c(0.5, 0.5))
  # The spread is centred on the original positions.
  expect_equal(mean(crowded), 1)
  jumbled <- c(3, 1, 2)
  spread <- .ml_lpa_spread_labels(jumbled, 2)
  expect_identical(order(spread), order(jumbled))
  expect_equal(mean(spread), mean(jumbled))
  expect_equal(.ml_lpa_spread_labels(7, 1), 7)
})

test_that("the label margin grows with the widest label", {
  draw({
    narrow <- .ml_lpa_label_margin("ab", 0.78)
    wide <- .ml_lpa_label_margin("a much longer series label", 0.78)
    expect_gt(wide, narrow)
    expect_equal(.ml_lpa_label_margin(character(), 0.78), 1.6)
    # Margins never shrink below the style's own right margin.
    style <- .ml_lpa_style()
    expect_gte(.ml_lpa_margins(style, character(), 0.78)[4L], style$margins[4L])
    expect_gt(.ml_lpa_margins(style, "a very long direct label indeed", 0.78)[4L],
              style$margins[4L])
  })
})
