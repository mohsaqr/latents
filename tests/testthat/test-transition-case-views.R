# The case views (raincloud, parallel, pairs) of a transition fit from lta()
# with homogeneous transitions: each answered occasion is a case, counted for
# its most likely state; occasions with nothing observed are not drawn.

.case_view_data <- function() {
  set.seed(11)
  n_units <- 40L
  n_occasions <- 5L
  frame <- data.frame(unit = rep(seq_len(n_units), each = n_occasions),
                      occasion = rep(seq_len(n_occasions), n_units))
  state <- unlist(lapply(seq_len(n_units), \(i) {
    cumsum(stats::runif(n_occasions) > 0.75) %% 2 + 1
  }))
  frame$a <- stats::rnorm(nrow(frame), ifelse(state == 2, 2, -2))
  frame$b <- stats::rnorm(nrow(frame), ifelse(state == 2, 1, -1))
  # Whole occasions unanswered, and one occasion with only `a` answered.
  frame[c(4L, 17L, 33L, 90L), c("a", "b")] <- NA
  frame$b[c(8L, 120L)] <- NA
  frame
}

.case_view_fit <- function(frame = .case_view_data()) {
  lta(frame, c("a", "b"), id = "unit", time = "occasion", n_profiles = 2,
      missing = "fiml", n_starts = 3, seed = 1)
}

.layer_with <- function(plot, column) {
  frames <- Filter(is.data.frame, lapply(plot$layers, \(layer) layer$data))
  Filter(\(frame) column %in% names(frame), frames)[[1L]]
}

test_that("a transition fit draws its cases without warnings", {
  skip_if_not_installed("ggplot2")
  fit <- .case_view_fit()
  expect_s3_class(fit, "multilpa_transitions")
  views <- list(plot_raincloud(fit), plot_parallel(fit), plot_pairs(fit),
                plot_profiles(fit, what = "raincloud"),
                plot(fit, what = "pairs"))
  invisible(lapply(views, \(view) {
    expect_s3_class(view, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(view))
  }))
})

test_that("raincloud medians and cases follow the assigned states", {
  skip_if_not_installed("ggplot2")
  frame <- .case_view_data()
  fit <- .case_view_fit(frame)
  drawn <- plot_raincloud(fit)
  key <- .gg_profile_key(fit)
  state_of <- \(labels) key$profile[match(as.character(labels),
                                          as.character(key$label))]
  box <- .layer_with(drawn, "median")
  box$profile <- state_of(box$profile_label)
  points <- Filter(\(frame) all(c("x", "y") %in% names(frame)) &&
                     !any(c("height", "median") %in% names(frame)),
                   Filter(is.data.frame,
                          lapply(drawn$layers, \(layer) layer$data)))[[1L]]
  points$profile <- state_of(points$profile_label)

  assigned <- get_results(fit, "assignments", data = frame)
  expected <- do.call(rbind, lapply(c("a", "b"), \(name) {
    values <- assigned[[name]]
    keep <- is.finite(values) & !is.na(assigned$profile)
    by_state <- split(values[keep], assigned$profile[keep])
    data.frame(indicator = name, profile = as.integer(names(by_state)),
               median = vapply(by_state, stats::median, numeric(1)),
               n = lengths(by_state))
  }))
  matched <- merge(expected, subset(box, select = c(indicator, profile,
                                                     median)),
                   by = c("indicator", "profile"), suffixes = c("", "_drawn"))
  expect_identical(nrow(matched), nrow(expected))
  expect_equal(matched$median_drawn, matched$median, tolerance = 1e-12)

  counted <- as.data.frame(table(indicator = points$indicator,
                                 profile = points$profile),
                           stringsAsFactors = FALSE)
  counted$profile <- as.integer(counted$profile)
  counted <- merge(expected, counted, by = c("indicator", "profile"))
  expect_identical(as.integer(counted$Freq), as.integer(counted$n))
  # Every answered value is drawn once, and nothing else.
  expect_identical(nrow(points), sum(is.finite(frame$a)) + sum(is.finite(frame$b)))
})

test_that("plot_all() of a transition fit includes the case views", {
  skip_if_not_installed("ggplot2")
  fit <- .case_view_fit()
  expect_true(all(c("raincloud", "parallel", "pairs") %in%
                    .multilpa_supported_views(fit)))
  drawn <- plot_all(fit)
  expect_true(all(c("raincloud", "parallel", "pairs") %in% names(drawn)))
})

# The case layers of a parallel plot: those coloured by profile, which carry
# the certainty of each case.
.parallel_case_layers <- function(plot) {
  Filter(\(layer) is.data.frame(layer$data) &&
           "certainty" %in% names(layer$data), plot$layers)
}

test_that("parallel draws every answered value, breaking lines at gaps", {
  skip_if_not_installed("ggplot2")
  frame <- .case_view_data()
  drawn <- plot_parallel(.case_view_fit(frame))
  layers <- .parallel_case_layers(drawn)
  # an occasion that answered only `a` is a point, not an invisible line
  geoms <- vapply(layers, \(layer) class(layer$geom)[1L], character(1))
  expect_setequal(geoms, c("GeomLine", "GeomPoint"))
  values <- do.call(rbind, lapply(layers, \(layer) layer$data))
  expect_identical(nrow(values),
                   sum(is.finite(frame$a)) + sum(is.finite(frame$b)))
  expect_no_warning(ggplot2::ggplot_build(drawn))
})

test_that("a missing middle indicator breaks the line instead of bridging it", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  frame <- data.frame(a = c(stats::rnorm(40, -2), stats::rnorm(40, 2)),
                      b = c(stats::rnorm(40, -1), stats::rnorm(40, 1)),
                      c = c(stats::rnorm(40, -1), stats::rnorm(40, 1)))
  frame$b[c(3L, 50L)] <- NA
  fit <- lpa(frame, c("a", "b", "c"), n_profiles = 2, missing = "fiml",
             n_starts = 2, seed = 1)
  lines <- Filter(\(layer) inherits(layer$geom, "GeomLine"),
                  .parallel_case_layers(plot_parallel(fit)))[[1L]]$data
  steps <- tapply(lines$position, lines$run, \(p) all(diff(sort(p)) == 1L))
  expect_true(all(steps))
  # the two cases with `b` missing appear as two points each
  points <- Filter(\(layer) inherits(layer$geom, "GeomPoint"),
                   .parallel_case_layers(plot_parallel(fit)))[[1L]]$data
  expect_identical(nrow(points), 4L)
})

test_that("pairs points sit in the panel of the case's own state", {
  skip_if_not_installed("ggplot2")
  frame <- .case_view_data()
  fit <- .case_view_fit(frame)
  drawn <- plot_pairs(fit)
  key <- .gg_profile_key(fit)
  scatter <- Filter(\(frame) all(c("x", "y", "uncertainty") %in%
                                   names(frame)),
                    lapply(drawn$layers, \(layer) layer$data))[[1L]]
  drawn_counts <- table(key$profile[match(as.character(scatter$profile_label),
                                          as.character(key$label))])
  assigned <- get_results(fit, "assignments", data = frame)
  both <- is.finite(assigned$a) & is.finite(assigned$b) &
    !is.na(assigned$profile)
  expect_identical(as.vector(drawn_counts),
                   as.vector(table(assigned$profile[both])))
  # the diagonal densities carry the state markers
  peaks <- Filter(\(layer) inherits(layer$geom, "GeomPoint") &&
                    is.data.frame(layer$data) &&
                    "height" %in% names(layer$data), drawn$layers)
  expect_length(peaks, 1L)
  expect_true("shape" %in% names(peaks[[1L]]$mapping))
})

test_that("a general transition fit still refuses the case views", {
  skip_if_not_installed("ggplot2")
  frame <- .case_view_data()
  general <- quietly(lta(frame, c("a", "b"), id = "unit", time = "occasion",
                         n_profiles = 2, missing = "fiml", mover_stayer = TRUE,
                         n_starts = 2, seed = 1))
  expect_s3_class(general, "multilpa_lta")
  expect_error(plot_raincloud(general), class = "latents_nothing_to_plot")
  expect_error(plot_pairs(general), class = "latents_nothing_to_plot")
})
