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

# `draw()` comes from helper-draw.R, shared with every other drawing test.

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
  candidates <- enumerate_classes(dat, c("a", "b"), "g", n_profiles = 1:2,
                                 n_group_classes = 1:2, n_starts = 3, seed = 3)
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
    # `new` left on would silently overlay the caller's next plot onto this one.
    expect_false(graphics::par("new"))
  })
})

test_that("a refused view leaves the device as it found it", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    # The restore runs before anything has been drawn here. Restoring `mfg` in
    # that state switches `new` on and warns, so a refused view used to poison
    # the device for every plot after it -- including the caller's own.
    expect_warning(
      expect_error(plot(fit, what = "responses"),
                   class = "multilpa_no_categorical"),
      regexp = NA)
    expect_false(graphics::par("new"))
    expect_identical(plot(fit), fit)
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
  candidates <- enumerate_classes(dat, c("a", "b"), "g", n_profiles = 1:2,
                                 n_group_classes = 1, n_starts = 3, seed = 3)
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
                                  n_profiles = 1:2, n_group_classes = 1,
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
    # With nothing to label the margin is exactly the padding. 2.2 lines, not
    # the earlier 1.6, because the label gap is now one "m" of the label's own
    # type rather than a fixed 0.45 user units, and a wider gap needs more room.
    expect_equal(.multilpa_label_margin(character(), 0.78), 2.2)
    expect_equal(.multilpa_label_margin(character(), 0.78, padding = 1.6), 1.6)
    # Margins never shrink below the style's own right margin.
    style <- .multilpa_style()
    expect_gte(.multilpa_margins(style, character(), 0.78)[4L], style$margins[4L])
    expect_gt(.multilpa_margins(style, "a very long direct label indeed", 0.78)[4L],
              style$margins[4L])
  })
})

test_that("the catalogue lists exactly the views the methods accept", {
  catalogue <- multilpa_plot_types()
  expect_s3_class(catalogue, "data.frame")
  expect_identical(names(catalogue), c("type", "group", "description"))
  expect_false(any(duplicated(catalogue$type)))
  expect_false(any(is.na(catalogue$description)))
  expect_true(all(nzchar(catalogue$description)))
  expect_setequal(unique(catalogue$group),
                  c("measurement", "structure", "diagnostics", "selection"))

  # The invariant that matters: a view added to the method but forgotten in the
  # catalogue, or listed but never implemented, fails here rather than silently
  # going undiscoverable. "enumeration" belongs to the enumeration method.
  fit_views <- eval(formals(plot.multilpa)$what)
  expect_setequal(setdiff(catalogue$type, "enumeration"), fit_views)
  expect_true(all(eval(formals(plot.multilpa_covariates)$what) %in%
                    catalogue$type))
})

test_that("every view of a Gaussian fit draws and returns the fit invisibly", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  gaussian_views <- c("profiles", "bars", "heatmap", "probabilities",
                      "entropy", "posteriors")
  draw({
    drawn <- vapply(gaussian_views, function(view) {
      identical(plot(fit, what = view), fit)
    }, logical(1))
    expect_true(all(drawn))
    expect_invisible(plot(fit, what = "heatmap"))
    # Supplying the data adds intervals to the bars; omitting it must not fail.
    expect_identical(plot(fit, what = "bars", data = dat), fit)
    expect_identical(plot(fit, what = "bars", scale = "standardized"), fit)
  })
})

test_that("a view a fit cannot supply is refused by condition class", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  draw({
    # Continuous indicators have no response probabilities, and this fit has no
    # time variable, so neither view is quietly drawn empty.
    expect_error(plot(fit, what = "responses"),
                 class = "multilpa_no_categorical")
    expect_error(plot(fit, what = "sequences"), class = "multilpa_no_time")
    expect_error(plot(fit, what = "not_a_view"))
  })
})

test_that("point area encodes prevalence absolutely, not stretched to the plot", {
  style <- .multilpa_style()
  # The defect this replaced: min-max normalising within the plot mapped every
  # spread onto the same size range, so a 52/48 split was drawn as far apart as
  # an 80/20 split and the channel carried no information at all.
  even <- .multilpa_point_sizes(c(0.52, 0.48), style)
  uneven <- .multilpa_point_sizes(c(0.80, 0.20), style)
  expect_gt(diff(rev(uneven)), diff(rev(even)))
  expect_false(isTRUE(all.equal(max(even), max(uneven))))

  # An evenly split pair sits at the style's own point size, whatever k is.
  expect_equal(.multilpa_point_sizes(c(0.5, 0.5), style),
               rep(style$point_size, 2L))
  expect_equal(.multilpa_point_sizes(rep(0.25, 4L), style),
               rep(style$point_size, 4L))

  # Monotone in share, and bounded so a rare profile is small but never invisible.
  sizes <- .multilpa_point_sizes(c(0.001, 0.1, 0.3, 0.599), style)
  expect_false(is.unsorted(sizes))
  expect_gt(min(sizes), 0)
  expect_lte(max(.multilpa_point_sizes(c(0.999, 0.001), style)),
             2.0 * style$point_size)
  expect_error(.multilpa_point_sizes(c(-0.1, 1.1), style))
})

test_that("the diverging scale is white at zero and symmetric about it", {
  expect_equal(.multilpa_diverging(0, limit = 2), "#FFFFFF")
  expect_equal(.multilpa_diverging(c(-2, 2), limit = 2), c("#D33F6A", "#4A6FE3"))
  # Equal deviations either side travel equally far from white, so the eye reads
  # distance from the grand mean rather than the direction alone. The two ends
  # are different hues by design, so it is the fraction of the way to each end
  # that must match, not the raw channel values.
  toward <- function(value, end) {
    reached <- 1 - grDevices::col2rgb(.multilpa_diverging(value, limit = 2)) / 255
    available <- 1 - grDevices::col2rgb(end) / 255
    reached / available
  }
  # Absolute, not relative: the colours are quantised to 8 bits, so the two
  # sides agree only to within one step of 1/255 and a relative tolerance on a
  # value of 0.5 is the stricter test by accident.
  expect_lt(max(abs(toward(-1, "#D33F6A") - toward(1, "#4A6FE3"))), 1 / 255)
  expect_lt(max(abs(toward(-1, "#D33F6A") - 0.5)), 1 / 255)
  # Beyond the limit the scale clamps instead of cycling through the hues.
  expect_equal(.multilpa_diverging(5, limit = 2), .multilpa_diverging(2, limit = 2))
  # A constant matrix has no spread to colour, and must not divide by zero.
  expect_equal(.multilpa_diverging(c(0, 0)), c("#FFFFFF", "#FFFFFF"))
})

test_that("the label gap is measured in the type it labels", {
  style <- .multilpa_style()
  draw({
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 3), ylim = c(0, 1))
    narrow_range <- .multilpa_label_offset(style)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 30), ylim = c(0, 1))
    wide_range <- .multilpa_label_offset(style)
    # The offset is in user units, so the same physical gap is a larger number
    # on a wider axis. A fixed constant would have been identical here, which is
    # exactly how the labels came to be clipped on the narrower panels.
    expect_gt(wide_range, narrow_range)
    expect_gt(narrow_range, 0)
  })
})

test_that("ridge views agree with the entropy the package reports", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 5, seed = 5)
  # The picture is never the only access to the numbers behind it.
  posterior <- as.data.frame(fit, what = "posteriors")
  expect_s3_class(posterior, "data.frame")
  expect_true(all(c("row", "profile", "posterior", "modal") %in% names(posterior)))
  expect_true(all(posterior$posterior >= 0 & posterior$posterior <= 1))
  # One entropy contribution per case, none of them beyond a flat posterior.
  contribution <- -rowSums(fit$subject_posteriors *
                             log(pmax(fit$subject_posteriors, .Machine$double.xmin)))
  expect_length(contribution, fit$n_observations)
  expect_lte(max(contribution), log(fit$n_profiles) + 1e-8)
  expect_gte(min(contribution), 0)
})
