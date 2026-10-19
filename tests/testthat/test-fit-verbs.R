# The named fit plot functions on every fitted-model class: each draws the
# view of the fit's own plot() that .latents_named_plot_registry() maps, or
# refuses with a classed condition naming the views that plot() offers.

fit_verb_cache <- new.env(parent = emptyenv())

# One small fit of every fitted-model class (and of each kind within a class),
# built once for the file.
fit_verb_fits <- function() {
  if (!is.null(fit_verb_cache$fits)) return(fit_verb_cache$fits)
  set.seed(21)
  persons <- 60L
  person_class <- sample.int(2L, persons, replace = TRUE)
  growth <- data.frame(id = rep(seq_len(persons), each = 5L),
                       time = rep(0:4, persons))
  growth$y <- c(2, 6)[person_class[growth$id]] +
    c(0.5, -0.4)[person_class[growth$id]] * growth$time +
    stats::rnorm(persons)[growth$id] + stats::rnorm(nrow(growth), 0, 0.7)

  set.seed(1)
  ratings <- data.frame(
    team = rep(seq_len(40L), each = 6L),
    climate = rep(stats::rnorm(40L, rep(c(-1, 1), each = 20L), 0.5), each = 6L) +
      stats::rnorm(240L),
    support = stats::rnorm(240L))

  set.seed(5)
  team_type <- rep(1:2, each = 40L)
  profile_means <- rbind(c(-1.5, -1), c(1.5, 1))
  cross <- do.call(rbind, lapply(seq_along(team_type), function(j) {
    shares <- if (team_type[j] == 1L) c(0.8, 0.2) else c(0.25, 0.75)
    k <- sample.int(2L, 6L, replace = TRUE, prob = shares)
    data.frame(team = sprintf("t%02d", j),
               y1 = stats::rnorm(6L, profile_means[k, 1L], 0.6),
               y2 = stats::rnorm(6L, profile_means[k, 2L], 0.6))
  }))

  set.seed(5)
  imputed_base <- data.frame(g = rep(seq_len(30L), each = 6L),
                             z = stats::rnorm(180L))
  state <- 2L - stats::rbinom(180L, 1L, stats::plogis(1.2 * imputed_base$z))
  imputed_base$y1 <- stats::rnorm(180L, c(0, 2.5)[state])
  imputed_base$y2 <- stats::rnorm(180L, c(0, 2)[state])
  imputations <- lapply(seq_len(3L), function(index) {
    copy <- imputed_base
    copy$z[seq_len(10L)] <- stats::rnorm(10L)
    copy
  })

  set.seed(8)
  panel <- do.call(rbind, lapply(seq_len(40L), function(j) {
    state <- cumsum(stats::runif(4L) < 0.3) %% 2L + 1L
    data.frame(id = j, time = seq_len(4L),
               y = stats::rnorm(4L, c(-1.5, 1.5)[state], 0.6))
  }))

  engagement <- subset(course_engagement, student <= 40)
  fits <- quietly(list(
    growth_mixture = mixture_regression(
      y ~ time, growth, n_classes = 2, id = "id", class_level = "group",
      random = ~ 1, n_starts = 1, seed = 1),
    trajectory_regression = mixture_regression(
      y ~ time, growth, n_classes = 2, id = "id", class_level = "group",
      n_starts = 1, seed = 1),
    mixture_regression = mixture_regression(score ~ hours, study_hours, 2,
                                            n_starts = 1, seed = 1),
    two_level_regression = mixture_regression(
      score ~ hours, study_hours, 2, id = "student", n_group_classes = 2,
      n_starts = 1, seed = 1, vcov_type = "none"),
    additive = multilpa(ratings, c("climate", "support"), "team",
                        n_group_classes = 2, family = "additive",
                        n_starts = 1, seed = 1),
    dispersion = multilpa(ratings, c("climate", "support"), "team",
                          n_group_classes = 2, family = "dispersion",
                          n_starts = 1, seed = 1),
    cross_level = multilpa(cross, c("y1", "y2"), "team", n_profiles = 2,
                           n_group_classes = 2, family = "full_cross_level",
                           n_starts = 1, seed = 1),
    pooled = pool_imputations(imputations, c("y1", "y2"), "g", 2L,
                              n_group_classes = 1L, profile_covariates = "z",
                              n_starts = 1, seed = 1),
    general_transitions = lta(panel, "y", "id", n_profiles = 2, time = "time",
                              transitions = "occasion", n_starts = 1, seed = 1),
    transitions = lta(engagement, c("browse", "lectures"), "student",
                      n_profiles = 2, time = "sequence", n_starts = 1, seed = 1),
    lpa = multilpa(engagement, c("browse", "lectures"), "student", 2, 2,
                   n_starts = 1, seed = 1),
    lca = lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
              n_classes = 2, n_starts = 1, seed = 1, tol = 1e-10),
    covariates = multilpa(imputed_base, c("y1", "y2"), "g", 2, 1,
                          profile_covariates = "z", n_starts = 1, seed = 1,
                          tol = 1e-10)),
    classes = c(.multilpa_expected_warnings, "latents_extreme_coefficients",
                "latents_inexact_optimum"))
  fit_verb_cache$fits <- fits
  fits
}

# Every named fit function with each of its `what` values: the function, its
# name, the `what` (`""` for a function without one) and the call itself.
fit_verb_calls <- function() {
  grouped <- list(
    plot_profiles = c("lines", "bars", "heatmap", "raincloud", "parallel",
                      "pairs", "all"),
    plot_diagnostics = c("all", "sizes", "posteriors", "entropy", "avepp",
                         "group_sizes", "group_posteriors", "group_entropy",
                         "group_avepp"))
  single <- c("plot_responses", "plot_sequences", "plot_prevalence",
              "plot_transitions")
  rbind(
    data.frame(verb = rep(names(grouped), lengths(grouped)),
               what = unlist(grouped, use.names = FALSE)),
    data.frame(verb = single, what = ""))
}

# Call one named function on one fit, returning the plot or the condition.
call_fit_verb <- function(fit, verb, what) {
  arguments <- if (nzchar(what)) list(fit, what = what) else list(fit)
  tryCatch(suppressMessages(do.call(verb, arguments)),
           error = function(condition) condition)
}

# Two drawings are the same plot: the same data in every layer and the same
# titles.
same_fit_drawing <- function(named, generic) {
  expect_equal(named$data, generic$data)
  expect_equal(lapply(named$layers, \(layer) layer$data),
               lapply(generic$layers, \(layer) layer$data))
  expect_identical(named$labels$title, generic$labels$title)
  expect_identical(named$labels$subtitle, generic$labels$subtitle)
}

test_that("every fit class is a kind the named functions know", {
  skip_on_cran()
  fits <- fit_verb_fits()
  kinds <- vapply(fits, .latents_plot_kind, character(1))
  expect_identical(unname(kinds), c(
    "growth_mixture", "trajectory_regression", "mixture_regression",
    "mixture_regression", "group_class_family", "group_class_family",
    "cross_level", "pooled", "general_transitions",
    rep("profile_model", 4L)))
  expect_true(is.na(.latents_plot_kind(data.frame(a = 1))))
})

test_that("the registry maps each function once, onto views the fit offers", {
  skip_on_cran()
  registry <- .latents_named_plot_registry()
  expect_s3_class(registry, "data.frame")
  expect_named(registry, c("kind", "verb", "what", "view"))
  expect_false(anyDuplicated(registry[c("kind", "verb", "what")]) > 0L)
  # Every mapped verb and `what` is one the function has.
  calls <- fit_verb_calls()
  expect_true(all(paste(registry$verb, registry$what) %in%
                    paste(calls$verb, calls$what)))
  fits <- fit_verb_fits()
  kinds <- vapply(fits, .latents_plot_kind, character(1))
  # The views each kind's fits offer between them: the profile models are
  # three classes, and not every one offers every view (only a transition fit
  # draws transitions).
  offered <- lapply(split(fits, kinds), function(same_kind) {
    unique(unlist(lapply(same_kind, .latents_plot_offered_views)))
  })
  invisible(Map(function(kind_name, views) {
    rows <- subset(registry, kind == kind_name)
    expect_gt(nrow(rows), 0L)
    # A single-view plot() takes no `what`, and its row says so with NA.
    expect_identical(is.na(rows$view), rep(length(views) == 0L, nrow(rows)))
    expect_true(all(rows$view %in% c(views, NA_character_)))
  }, names(offered), offered))
})

test_that("each named function on each fit draws or refuses with a classed condition", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fits <- fit_verb_fits()
  calls <- fit_verb_calls()
  outcomes <- do.call(rbind, Map(function(fit_name, fit) {
    do.call(rbind, Map(function(verb, what) {
      result <- call_fit_verb(fit, verb, what)
      data.frame(fit = fit_name, kind = .latents_plot_kind(fit), verb = verb,
                 what = what,
                 drawn = inherits(result, c("ggplot", "latents_plots")),
                 classed = inherits(result, "error") &&
                   any(grepl("^latents_", class(result))),
                 refused = inherits(result, "latents_nothing_to_plot"))
    }, calls$verb, calls$what))
  }, names(fits), fits))
  # Never an unclassed error.
  expect_true(all(outcomes$drawn | outcomes$classed))
  # Away from profile, class and transition models, a function draws exactly
  # when the registry maps it (its "all" when it maps any of its views), and
  # refuses with latents_nothing_to_plot otherwise.
  registry <- .latents_named_plot_registry()
  other <- subset(outcomes, kind != "profile_model")
  mapped <- paste(other$kind, other$verb, other$what) %in%
    paste(registry$kind, registry$verb, registry$what) |
    (other$what == "all" & paste(other$kind, other$verb) %in%
       paste(registry$kind, registry$verb))
  expect_identical(other$drawn, mapped)
  expect_true(all(other$refused[!mapped]))
  # The profile models draw what they drew before this change.
  lpa <- subset(outcomes, fit == "lpa")
  expect_true(all(subset(lpa, verb == "plot_profiles")$drawn))
  expect_true(all(subset(lpa, verb == "plot_diagnostics")$drawn))
  expect_identical(subset(outcomes, verb == "plot_transitions")$drawn,
                   outcomes$fit[outcomes$verb == "plot_transitions"] %in%
                     c("general_transitions", "transitions"))
})

test_that("a mapped function draws exactly what plot(what = ) draws", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fits <- fit_verb_fits()
  registry <- subset(.latents_named_plot_registry(), kind != "profile_model")
  checked <- Map(function(fit_name, fit) {
    rows <- subset(registry, kind == .latents_plot_kind(fit))
    Map(function(verb, what, view) {
      named <- call_fit_verb(fit, verb, what)
      expect_s3_class(named, "ggplot")
      generic <- if (is.na(view)) plot(fit) else plot(fit, what = view)
      same_fit_drawing(named, generic)
      # Titles reach the mapped view.
      titled <- do.call(verb, c(list(fit), if (nzchar(what)) list(what = what),
                                list(main = "custom", subtitle = "sub")))
      expect_identical(titled$labels$title, "custom")
      expect_identical(titled$labels$subtitle, "sub")
      fit_name
    }, rows$verb, rows$what, rows$view)
  }, names(fits), fits)
  expect_gt(length(unlist(checked)), 10L)
  # "all" draws the mapped view and names the others.
  expect_message(trajectories <- plot_profiles(fits$growth_mixture, what = "all"),
                 "Not drawn for this model: bars, heatmap")
  expect_s3_class(trajectories, "latents_plots")
  expect_named(trajectories, "lines")
  same_fit_drawing(trajectories$lines,
                   plot(fits$growth_mixture, what = "trajectories"))
  expect_message(posteriors <- plot_diagnostics(fits$mixture_regression),
                 "sizes, entropy, avepp")
  expect_named(posteriors, "posteriors")
  same_fit_drawing(posteriors$posteriors,
                   plot(fits$mixture_regression, what = "posteriors"))
})

test_that("the mapping is the documented one, view by view", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fits <- fit_verb_fits()
  # Written out by hand, independently of the registry, so a changed or
  # dropped row fails here.
  same_fit_drawing(plot_profiles(fits$growth_mixture),
                   plot(fits$growth_mixture, what = "trajectories"))
  same_fit_drawing(plot_diagnostics(fits$growth_mixture, what = "posteriors"),
                   plot(fits$growth_mixture, what = "posteriors"))
  same_fit_drawing(plot_profiles(fits$trajectory_regression),
                   plot(fits$trajectory_regression, what = "trajectories"))
  same_fit_drawing(plot_diagnostics(fits$trajectory_regression,
                                    what = "posteriors"),
                   plot(fits$trajectory_regression, what = "posteriors"))
  same_fit_drawing(plot_profiles(fits$mixture_regression),
                   plot(fits$mixture_regression, what = "fitted"))
  same_fit_drawing(plot_diagnostics(fits$mixture_regression, what = "posteriors"),
                   plot(fits$mixture_regression, what = "posteriors"))
  same_fit_drawing(plot_profiles(fits$two_level_regression),
                   plot(fits$two_level_regression, what = "fitted"))
  same_fit_drawing(plot_profiles(fits$additive),
                   plot(fits$additive, what = "means"))
  same_fit_drawing(plot_profiles(fits$dispersion),
                   plot(fits$dispersion, what = "means"))
  same_fit_drawing(plot_profiles(fits$cross_level),
                   plot(fits$cross_level, what = "profiles"))
  same_fit_drawing(plot_prevalence(fits$cross_level),
                   plot(fits$cross_level, what = "composition"))
  same_fit_drawing(plot_profiles(fits$pooled),
                   plot(fits$pooled, what = "measurement"))
  same_fit_drawing(plot_prevalence(fits$pooled),
                   plot(fits$pooled, what = "profile"))
  same_fit_drawing(plot_transitions(fits$general_transitions),
                   plot(fits$general_transitions))
  # The profile models' mapping is the one their tests have always checked.
  same_fit_drawing(plot_prevalence(fits$lpa),
                   plot(fits$lpa, what = "probabilities"))
  same_fit_drawing(plot_transitions(fits$transitions),
                   plot(fits$transitions, what = "transitions"))
  same_fit_drawing(plot_profiles(fits$lpa, what = "heatmap"),
                   plot(fits$lpa, what = "heatmap"))
})

test_that("refusals are classed and name the views of the fit's plot()", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fits <- fit_verb_fits()
  # A function or `what` with no view of that fit.
  refused <- expect_error(plot_profiles(fits$growth_mixture, what = "bars"),
                          class = "latents_nothing_to_plot")
  expect_s3_class(refused, "latents_bad_argument")
  expect_match(conditionMessage(refused), "\"trajectories\", \"individuals\"",
               fixed = TRUE)
  expect_match(conditionMessage(refused), "plot_diagnostics(what = \"posteriors\")",
               fixed = TRUE)
  expect_error(plot_diagnostics(fits$growth_mixture, what = "sizes"),
               class = "latents_nothing_to_plot")
  expect_error(plot_diagnostics(fits$additive), class = "latents_nothing_to_plot")
  expect_error(plot_transitions(fits$cross_level),
               class = "latents_nothing_to_plot")
  expect_error(plot_sequences(fits$pooled), class = "latents_nothing_to_plot")
  expect_error(plot_responses(fits$mixture_regression),
               class = "latents_nothing_to_plot")
  # A general lta() fit points to the function that draws it.
  expect_error(plot_profiles(fits$general_transitions), "plot_transitions()",
               fixed = TRUE)
  expect_error(plot_profiles(fits$general_transitions),
               class = "latents_bad_argument")
  # A trajectory view is for fits that classify persons only: an observation
  # regression draws its fitted lines instead.
  expect_identical(plot_profiles(fits$mixture_regression)$labels$title,
                   plot(fits$mixture_regression, what = "fitted")$labels$title)
})

test_that("an argument the mapped view does not take is refused by name", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fits <- fit_verb_fits()
  check <- function(expr, argument) {
    refused <- expect_error(expr, class = "latents_bad_argument")
    expect_false(inherits(refused, "latents_nothing_to_plot"))
    expect_match(conditionMessage(refused), sprintf("`%s`", argument),
                 fixed = TRUE)
  }
  check(plot_profiles(fits$growth_mixture, intervals = FALSE), "intervals")
  check(plot_profiles(fits$additive, scale = "standardized"), "scale")
  check(plot_profiles(fits$cross_level, statistic = "median"), "statistic")
  check(plot_profiles(fits$pooled, labels = FALSE), "labels")
  check(plot_prevalence(fits$cross_level, labels = FALSE), "labels")
  # Also under "all": the argument is not silently dropped.
  check(plot_profiles(fits$growth_mixture, what = "all", scale = "standardized"),
        "scale")
  # Not a fitted model at all.
  not_fit <- expect_error(plot_profiles(data.frame(a = 1)),
                          class = "latents_bad_argument")
  expect_false(inherits(not_fit, "latents_nothing_to_plot"))
})
