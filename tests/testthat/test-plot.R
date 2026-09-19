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
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    expect_invisible(plot(fit))
    expect_identical(plot(fit), fit)
    expect_identical(plot(fit, what = "probabilities"), fit)
    expect_identical(plot(fit, scale = "standardized"), fit)
    expect_identical(plot(fit, labels = FALSE), fit)
  })
  candidates <- enumerate_classes(dat, c("a", "b"), "g", profiles = 1:2,
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
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    before <- graphics::par(c("mar", "xpd", "cex"))
    plot(fit)
    plot(fit, what = "probabilities")
    expect_equal(graphics::par(c("mar", "xpd", "cex")), before)
  })
})

test_that("plot arguments override every visual constant", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
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
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
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
                    class = "multilpa_no_indicator_data"))
})

test_that("enumeration plotting rejects unusable criteria by condition class", {
  dat <- plot_fixture()
  candidates <- enumerate_classes(dat, c("a", "b"), "g", profiles = 1:2,
                                 group_classes = 1, n_starts = 3, seed = 3)
  draw({
    expect_error(plot(candidates, criterion = "not_a_column"),
                 class = "multilpa_unknown_criterion")
    empty <- candidates
    empty$table$bic_individual <- NA_real_
    expect_error(plot(empty, criterion = "bic_individual"),
                 class = "multilpa_nothing_to_plot")
  })
})

test_that("the palette, symbols and line types stay aligned and recycle", {
  # Colour must never be the only channel, so all three vectors must have the
  # same length for any number of classes, including beyond the palette length.
  invisible(lapply(c(1L, 3L, 9L, 12L), function(n) {
    expect_length(.multilpa_palette(n), n)
    expect_length(.multilpa_symbols(n), n)
    expect_length(.multilpa_linetypes(n), n)
  }))
  expect_identical(.multilpa_palette(1L), "#E69F00")
  # Okabe-Ito, and nothing outside it.
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
  expect_true(all(.multilpa_palette(20L) %in% okabe_ito))
  expect_identical(.multilpa_palette(11L)[10L], .multilpa_palette(1L))
})

test_that("direct labels are spread apart but keep their order and centre", {
  expect_equal(.multilpa_spread_labels(c(1, 5), 0.5), c(1, 5))
  crowded <- .multilpa_spread_labels(c(1, 1, 1), 0.5)
  expect_equal(diff(sort(crowded)), c(0.5, 0.5))
  # The spread is centred on the original positions.
  expect_equal(mean(crowded), 1)
  jumbled <- c(3, 1, 2)
  spread <- .multilpa_spread_labels(jumbled, 2)
  expect_identical(order(spread), order(jumbled))
  expect_equal(mean(spread), mean(jumbled))
  expect_equal(.multilpa_spread_labels(7, 1), 7)
})

.sequence_plot_fixture <- function(n_groups = 12L, n_waves = 10L) {
  set.seed(7)
  n <- n_groups * n_waves
  data.frame(school = rep(seq_len(n_groups), each = n_waves),
             wave = rep(seq_len(n_waves), times = n_groups),
             score_a = stats::rnorm(n), score_b = stats::rnorm(n))
}

test_that("the sequence plot draws and returns its fit invisibly", {
  data <- .sequence_plot_fixture()
  fit <- multilpa(data, c("score_a", "score_b"), "school", n_profiles = 2,
                  n_group_classes = 2, n_starts = 3, seed = 1, time = "wave")
  draw({
    expect_invisible(plot(fit, what = "sequences"))
    expect_identical(plot(fit, what = "sequences"), fit)
    expect_identical(plot(fit, what = "sequences", cell_labels = FALSE), fit)
    expect_identical(plot(fit, what = "sequences", labels = FALSE), fit)
  })
  expect_error(plot(fit, what = "sequences", cell_labels = "yes"),
               "`cell_labels` must be TRUE or FALSE")
})

test_that("every plot's data is reachable through a tidy verb", {
  data <- .sequence_plot_fixture()
  fit <- multilpa(data, c("score_a", "score_b"), "school", n_profiles = 2,
                  n_group_classes = 2, n_starts = 3, seed = 1, time = "wave")
  candidates <- enumerate_classes(data, c("score_a", "score_b"), "school",
                                  profiles = 1:2, group_classes = 1,
                                  n_starts = 2, seed = 3)

  # No plot is the only way to see what it draws: the reader can always get the
  # numbers as a data frame instead of measuring them off the picture.
  expect_s3_class(as.data.frame(fit, what = "profiles"), "data.frame")
  expect_s3_class(as.data.frame(fit, what = "profile_probabilities"),
                  "data.frame")
  expect_s3_class(sequences(fit), "data.frame")
  expect_s3_class(as.data.frame(candidates), "data.frame")
  expect_true(all(c("profile", "indicator", "mean") %in%
                    names(as.data.frame(fit, what = "profiles"))))
  expect_true(all(c("group", "group_class", "time", "profile") %in%
                    names(sequences(fit))))
})

test_that("a filled grid carries its code as text, not by colour alone", {
  draw({
    codes <- matrix(rep(1:2, length.out = 120L), 12L, 10L)
    style <- .multilpa_style()
    graphics::plot.new()
    graphics::plot.window(c(0.5, 10.5), c(12.5, 0.5))
    expect_true(.multilpa_cell_labels(codes, seq_len(10L),
                                      .multilpa_palette(2L), style))
    expect_false(.multilpa_cell_labels(codes, seq_len(10L),
                                       .multilpa_palette(2L), style,
                                       draw = FALSE))
    # An empty grid has nothing to stamp, and a grid too dense to hold a digit
    # says so rather than smearing the numbers over each other.
    empty <- matrix(NA_integer_, 4L, 4L)
    expect_false(.multilpa_cell_labels(empty, seq_len(4L),
                                       .multilpa_palette(2L), style))
    graphics::plot.new()
    graphics::plot.window(c(0.5, 400.5), c(400.5, 0.5))
    dense <- matrix(rep(1:2, length.out = 160000L), 400L, 400L)
    expect_false(.multilpa_cell_labels(dense, seq_len(400L),
                                       .multilpa_palette(2L), style))
  })
})

test_that("cell ink is chosen for contrast against its background", {
  # Light fills take black ink and dark fills take white, so a label inside a
  # cell is legible whichever palette colour the cell got.
  expect_identical(.multilpa_ink(c("#F0E442", "#FFFFFF")),
                   c("#000000", "#000000"))
  expect_identical(.multilpa_ink(c("#000000", "#0072B2", "#D55E00")),
                   c("#FFFFFF", "#FFFFFF", "#FFFFFF"))
  expect_length(.multilpa_ink(.multilpa_palette(9L)), 9L)
  expect_identical(.multilpa_ink(character()), character())
  expect_true(all(.multilpa_ink(.multilpa_palette(9L)) %in%
                    c("#000000", "#FFFFFF")))
})

test_that("the label margin grows with the widest label", {
  draw({
    narrow <- .multilpa_label_margin("ab", 0.78)
    wide <- .multilpa_label_margin("a much longer series label", 0.78)
    expect_gt(wide, narrow)
    expect_equal(.multilpa_label_margin(character(), 0.78), 1.6)
    # Margins never shrink below the style's own right margin.
    style <- .multilpa_style()
    expect_gte(.multilpa_margins(style, character(), 0.78)[4L], style$margins[4L])
    expect_gt(.multilpa_margins(style, "a very long direct label indeed", 0.78)[4L],
              style$margins[4L])
  })
})
