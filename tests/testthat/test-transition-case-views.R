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
