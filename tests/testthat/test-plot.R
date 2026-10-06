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

test_that("plot methods return ggplot objects for every view", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  expect_plot(plot(fit))
  expect_plot(plot(fit, what = "probabilities"))
  expect_plot(plot(fit, scale = "standardized"))
  expect_plot(plot(fit, labels = FALSE))
  candidates <- enumerate_classes(dat, c("a", "b"), "g", n_profiles = 1:2,
                                 n_group_classes = 1:2, n_starts = 1, seed = 3)
  expect_plot(plot(candidates))
  expect_plot(plot(candidates, criterion = "sabic_individual"))
  expect_plot(plot(candidates, labels = FALSE, mark_minimum = FALSE))
  expect_plot(plot(candidates, what = "tree"))
})

test_that("plot refuses the styling arguments it no longer takes", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  # Styling is ggplot2's job now; a dropped argument is named, not ignored.
  expect_error(plot(fit, palette = c("#000000", "#FFFFFF")),
               class = "latents_bad_argument")
  expect_error(plot(fit, point_size = 2), class = "latents_bad_argument")
  expect_error(plot(fit, labels = "yes"), "`labels` must be TRUE or FALSE")
  # Titles are still the caller's to set.
  custom <- plot(fit, main = "custom title", subtitle = "custom subtitle")
  expect_identical(custom$labels$title, "custom title")
  expect_identical(custom$labels$subtitle, "custom subtitle")
})

test_that("standardizing uses the observed indicator scales", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  observed <- fit$indicator_data
  centre <- colMeans(observed)
  spread <- c(stats::sd(observed[, 1L]), stats::sd(observed[, 2L]))
  expected <- sweep(sweep(fit$means, 2L, centre, "-"), 2L, spread, "/")
  # The b indicator is on a far wider scale, so standardizing must shrink it
  # relative to a; this is the whole point of the option.
  raw_spread <- diff(range(fit$means[, "b"]))
  standardized_spread <- diff(range(expected[, "b"]))
  expect_lt(standardized_spread, raw_spread)
  skip_if_not_installed("ggplot2")
  drawn <- expect_plot(plot(fit, scale = "standardized"))$data
  expect_equal(drawn$mean, expected[cbind(drawn$profile, drawn$position)])
  stripped <- fit
  stripped$indicator_data <- NULL
  expect_error(plot(stripped, scale = "standardized"),
               class = "latents_no_indicator_data")
})

test_that("enumeration plotting rejects unusable criteria by condition class", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  dat <- plot_fixture()
  candidates <- enumerate_classes(dat, c("a", "b"), "g", n_profiles = 1:2,
                                 n_group_classes = 1, n_starts = 1, seed = 3)
  expect_error(plot(candidates, criterion = "not_a_column"),
               class = "latents_unknown_criterion")
  expect_error(plot(candidates, criterion = "structure"),
               class = "latents_unknown_criterion")
  empty <- candidates
  empty$table$bic_individual <- NA_real_
  expect_error(plot(empty, criterion = "bic_individual"),
               class = "latents_nothing_to_plot")
})

test_that("the palette, symbols and line types stay aligned and recycle", {
  skip_on_cran()
  # The base-graphics plots that remain (bootstrap, pooled, mixture
  # regression) still draw with these; colour must never be the only channel.
  invisible(lapply(c(1L, 3L, 9L, 12L), function(n) {
    expect_length(.multilpa_palette(n), n)
    expect_length(.multilpa_symbols(n), n)
    expect_length(.multilpa_linetypes(n), n)
  }))
  expect_identical(.multilpa_palette(1L), "#E69F00")
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
  expect_true(all(.multilpa_palette(20L) %in% okabe_ito))
  expect_identical(.multilpa_palette(11L)[10L], .multilpa_palette(1L))
  expect_false("#F0E442" %in% .multilpa_palette(6L))
  # The ggplot2 views draw from the same palette, yellow last.
  expect_true(all(.gg_okabe_ito %in% okabe_ito))
  expect_identical(.gg_okabe_ito[length(.gg_okabe_ito)], "#F0E442")
})

# Twelve schools do not identify the group-class split in this fixture: the
# likelihood is flat along it. Plain EM stops partway along that ridge; SQUAREM
# follows it to the boundary and rightly warns of a near-empty group class. The
# plot tests are about drawing, so their fits use plain EM.
.sequence_plot_fixture <- function(n_groups = 12L, n_waves = 10L) {
  set.seed(7)
  n <- n_groups * n_waves
  data.frame(school = rep(seq_len(n_groups), each = n_waves),
             wave = rep(seq_len(n_waves), times = n_groups),
             score_a = stats::rnorm(n), score_b = stats::rnorm(n))
}

test_that("the sequence plot draws every group at every position", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data <- .sequence_plot_fixture()
  fit <- multilpa(data, c("score_a", "score_b"), "school", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1, time = "wave",
                  acceleration = "none")
  drawn <- expect_plot(plot(fit, what = "sequences"))$data
  expect_identical(nrow(drawn), 12L * 10L)
  expect_setequal(drawn$profile, fit$subject_profiles)
  expect_plot(plot(fit, what = "sequences", cell_labels = FALSE))
  expect_plot(plot(fit, what = "sequences", labels = FALSE))
  expect_error(plot(fit, what = "sequences", cell_labels = "yes"),
               "`cell_labels` must be TRUE or FALSE")
})

test_that("every plot's data is reachable through a tidy verb", {
  skip_on_cran()
  data <- .sequence_plot_fixture()
  fit <- multilpa(data, c("score_a", "score_b"), "school", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1, time = "wave",
                  acceleration = "none")
  candidates <- enumerate_classes(data, c("score_a", "score_b"), "school",
                                  n_profiles = 1:2, n_group_classes = 1,
                                  n_starts = 1, seed = 3)
  # No plot is the only way to see what it draws: the reader can always get the
  # numbers as a data frame instead of measuring them off the picture.
  expect_s3_class(get_results(fit, "profiles"), "data.frame")
  expect_s3_class(get_results(fit, "profile_probabilities"), "data.frame")
  expect_s3_class(get_results(fit, "sequences"), "data.frame")
  expect_s3_class(as.data.frame(candidates), "data.frame")
  expect_true(all(c("profile", "indicator", "mean") %in%
                    names(get_results(fit, "profiles"))))
  expect_true(all(c("group", "group_class", "time", "profile") %in%
                    names(get_results(fit, "sequences"))))
})

test_that("the catalogue lists exactly the views the methods accept", {
  skip_on_cran()
  catalogue <- plot_views()
  expect_s3_class(catalogue, "data.frame")
  expect_identical(names(catalogue), c("type", "group", "description"))
  expect_false(any(duplicated(catalogue$type)))
  expect_false(any(is.na(catalogue$description)))
  expect_true(all(nzchar(catalogue$description)))
  expect_setequal(unique(catalogue$group),
                  c("measurement", "structure", "diagnostics", "selection",
                    "every"))
  expect_identical(catalogue$group[catalogue$type == "all"], "every")
  # A view added to a method but forgotten in the catalogue, or listed but
  # never implemented, fails here. The catalogue is the union over every
  # method; "enumeration" and "tree" belong to the enumeration method.
  fit_methods <- list(plot.multilpa, plot.multilpa_covariates,
                      plot.multilpa_transitions)
  fit_views <- unique(unlist(lapply(fit_methods, function(method)
    eval(formals(method)$what)), use.names = FALSE))
  grid_views <- eval(formals(plot.multilpa_enumeration)$what)
  expect_setequal(catalogue$type, c(fit_views, grid_views))
  expect_setequal(grid_views, c("enumeration", "tree"))
})

test_that("a view a fit cannot supply is refused by condition class", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  expect_error(plot(fit, what = "responses"), class = "latents_no_categorical")
  expect_error(plot(fit, what = "sequences"), class = "latents_no_time")
  expect_error(plot(fit, what = "not_a_view"))
})

test_that("case views agree with the entropy the package reports", {
  skip_on_cran()
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  posterior <- get_results(fit, "posteriors")
  expect_true(all(c("row", "profile", "posterior", "modal") %in% names(posterior)))
  skip_if_not_installed("ggplot2")
  drawn <- expect_plot(plot(fit, what = "entropy"))$data
  expect_length(drawn$value, fit$n_observations)
  expect_true(all(drawn$value >= 0 & drawn$value <= 1 + 1e-12))
  entropy <- get_results(fit, "entropy")
  expect_equal(1 - mean(drawn$value),
               entropy$relative_entropy[entropy$level == "individuals"])
  certainty <- expect_plot(plot(fit, what = "posteriors"))$data
  expect_equal(certainty$value, unname(apply(fit$subject_posteriors, 1L, max)))
})

test_that("bar intervals come from the fit's own data, not from the caller", {
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5,
                  tol = 1e-10)
  supplied <- .multilpa_mean_error_matrix(fit, dat)
  carried <- .multilpa_mean_error_matrix(fit, NULL)
  expect_identical(carried, supplied)
  expect_identical(dim(carried), c(2L, 2L))
  expect_true(all(carried > 0))
  skip_if_not_installed("ggplot2")
  bars <- expect_plot(plot(fit, what = "bars"))
  expect_equal(bars$data$upper - bars$data$mean,
               stats::qnorm(0.975) *
                 carried[cbind(bars$data$profile, bars$data$position)])
  expect_equal(plot(fit, what = "bars", data = dat)$data, bars$data)
})

test_that("a transition fit gets bars without whiskers, not an error", {
  skip_if_not_installed("ggplot2")
  data <- .sequence_plot_fixture()
  # The transition bars carry point estimates only, even where inference
  # would refuse: on these pure-noise data a probability sits at its bound and
  # parameter_inference() raises latents_boundary_fit.
  moves <- lta(data, c("score_a", "score_b"), "school",
               n_profiles = 2, time = "wave", n_starts = 1, seed = 1)
  expect_error(parameter_inference(moves), class = "latents_boundary_fit")
  bars <- expect_plot(plot(moves, what = "bars"))
  expect_true(all(is.na(bars$data$upper)))
})

test_that("the avepp view draws the average posterior matrix", {
  skip_on_cran()
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  averages <- .multilpa_average_posterior_matrix(fit$subject_posteriors)
  expect_equal(unname(rowSums(averages)), rep(1, ncol(averages)))
  table <- get_results(fit, "average_posteriors")
  individuals <- table[table$level == "individuals", , drop = FALSE]
  expect_equal(individuals$average_posterior, as.vector(t(averages)))
  skip_if_not_installed("ggplot2")
  drawn <- expect_plot(plot(fit, what = "avepp"))$data
  expect_equal(drawn$average_posterior,
               averages[cbind(drawn$row, drawn$column)])
})

test_that("a one-profile fit still draws its sizes", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 5)
  expect_plot(plot(fit, what = "sizes"))
  expect_error(plot(fit, what = "posteriors"), class = "latents_nothing_to_plot")
})

test_that("what = \"all\" returns every drawable view, named", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  dat <- plot_fixture()
  fit <- multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5)
  # The views this fit lacks the ingredients for are never attempted, so
  # nothing is refused and nothing needs saying.
  everything <- expect_plots(expect_no_message(plot(fit, what = "all")))
  expect_true(all(c("profiles", "bars", "pairs", "avepp") %in%
                    names(everything)))
  expect_false(any(c("responses", "sequences") %in% names(everything)))
  # A view that refuses for a reason only it knows is named, not dropped.
  one_profile <- multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 5)
  fewer <- expect_plots(plot(one_profile, what = "all"))
  expect_false(any(c("entropy", "posteriors") %in% names(fewer)))
})
