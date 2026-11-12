# The plotting verbs: plot() with every `what` its method offers, and each
# named plot_*() verb, on a fit of every family, an enumeration grid and a
# pooled result. A picture has no reference package, so the checks are: every
# view either draws to a device or refuses with a `latents_` condition (a view
# the fit has no data for), the views a family is documented to have do draw,
# an unknown view is refused by class, and the numbers behind two views are
# the fit's own tables.

sweep_module_plots <- function() {
  two_level <- sweep_two_level()
  panel <- sweep_longitudinal()$gaussian
  items <- sweep_categorical()$nominal
  #' Draw to a null device; a list of plots is drawn element by element
  draw <- function(value) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    if (inherits(value, "ggplot") || !is.list(value)) print(value) else lapply(value, print)
    invisible(TRUE)
  }
  outcome_of <- function(expression) {
    tryCatch({
      draw(expression())
      "drawn"
    }, error = function(e) {
      if (any(grepl("^(latents|multilpa)_", class(e)))) "refused by class" else
        sprintf("error: %s", conditionMessage(e))
    })
  }
  #' A view draws, or the fit has no such view and says so by class
  draws_or_refuses <- function(label, expression) {
    outcome <- outcome_of(expression)
    sweep_true(sprintf("%s draws, or is refused by class", label),
               outcome %in% c("drawn", "refused by class"),
               note = if (identical(outcome, "drawn")) NA_character_ else outcome)
  }
  must_draw <- function(label, expression) {
    outcome <- outcome_of(expression)
    sweep_true(sprintf("%s draws", label), identical(outcome, "drawn"),
               note = if (identical(outcome, "drawn")) NA_character_ else outcome)
  }
  #' An unknown view is refused by class; a method with one view has no `what`
  #' and reports the argument as unused, by class too.
  unknown_view <- function(object) {
    signalled <- character()
    outcome <- withCallingHandlers(
      tryCatch({
        draw(plot(object, what = "no-such-view"))
        "drawn"
      }, error = function(e) {
        if (any(grepl("^(latents|multilpa)_", class(e)))) "refused by class" else
          sprintf("unclassed error: %s", conditionMessage(e))
      }),
      warning = function(w) {
        signalled <<- c(signalled, class(w))
        invokeRestart("muffleWarning")
      })
    ok <- identical(outcome, "refused by class") ||
      (identical(outcome, "drawn") && "latents_unused_argument" %in% signalled)
    sweep_true("an unknown view is refused, or reported as unused, by class", ok,
               note = if (ok) NA_character_ else outcome)
  }
  fits <- list(
    profiles = function() lpa(two_level, c("y1", "y2", "y3"), 2L, n_starts = 3L, seed = 1L),
    classes = function() lca(items, names(items), 2L, n_starts = 3L, seed = 1L),
    two_level = function() multilpa(two_level, c("y1", "y2"), "id", n_profiles = 2L,
                                    n_group_classes = 2L, n_starts = 3L, seed = 1L),
    transitions = function() lta(panel$data, panel$vars, id = panel$id, time = panel$time,
                                 n_profiles = 2L, n_starts = 3L, seed = 1L),
    general_transitions = function() lta(panel$data, panel$vars, id = panel$id,
                                         time = panel$time, n_profiles = 2L,
                                         transitions = "occasion", n_starts = 3L, seed = 1L),
    additive = function() multilpa(two_level, c("y1", "y2"), "id", n_group_classes = 2L,
                                   family = "additive", n_starts = 3L, seed = 1L),
    cross_level = function() multilpa(two_level, c("y1", "y2"), "id", n_profiles = 2L,
                                      n_group_classes = 2L, family = "full_cross_level",
                                      n_starts = 3L, seed = 1L),
    regression = function() mixture_regression(hours ~ week, latents::study_hours,
                                               n_classes = 2L, n_starts = 3L, seed = 1L),
    growth = function() mixture_regression(score ~ wave, latents::growth_scores, n_classes = 2L,
                                           id = "student", class_level = "group",
                                           random = ~ 1, n_starts = 3L, seed = 1L),
    grid = function() enumerate_lpa(two_level, c("y1", "y2", "y3"), n_profiles = 1:3,
                                    n_starts = 3L, seed = 1L),
    pooled = function() {
      sets <- sweep_seeded(3L, function() lapply(1:3, function(i) {
        transform(two_level, y1 = y1 + stats::rnorm(length(y1), sd = 0.2))
      }))
      pool_imputations(sets, c("y1", "y2"), NULL, n_profiles = 2L, n_starts = 3L, seed = 1L)
    })
  # The view every family is documented to draw by default.
  verbs <- c("plot_profiles", "plot_responses", "plot_sequences", "plot_prevalence",
             "plot_diagnostics", "plot_transitions", "plot_sizes", "plot_entropy",
             "plot_posteriors", "plot_avepp", "plot_heatmap", "plot_bars", "plot_raincloud",
             "plot_parallel", "plot_pairs", "plot_trajectories", "plot_individuals",
             "plot_classification", "plot_coefficients", "plot_criteria", "plot_deltas",
             "plot_trust", "plot_tree", "plot_selection", "plot_enumeration", "plot_all")
  method_views <- function(object) {
    method <- Find(Negate(is.null), lapply(class(object), function(class) {
      utils::getS3method("plot", class, optional = TRUE)
    }))
    views <- tryCatch(eval(formals(method)$what), error = function(e) NULL)
    if (is.character(views)) views else character()
  }
  family_cases <- lapply(names(fits), function(name) {
    sweep_define(sprintf("views[fit=%s]", name), function() {
      object <- fits[[name]]()
      views <- method_views(object)
      rbind(
        must_draw("plot() with its default view", function() plot(object)),
        do.call(rbind, lapply(views, function(view) {
          draws_or_refuses(sprintf("plot(what = \"%s\")", view),
                           function() plot(object, what = view))
        })),
        unknown_view(object),
        do.call(rbind, lapply(verbs, function(verb) {
          draws_or_refuses(sprintf("%s()", verb), function() get(verb)(object))
        })))
    })
  })

  argument_cases <- list(
    sweep_define("plot_profiles-arguments", function() {
      object <- fits$profiles()
      options <- sweep_grid(statistic = c("mean", "median", "mode"),
                            scale = c("raw", "standardized"),
                            what = c("lines", "bars", "heatmap", "raincloud", "parallel",
                                     "pairs", "all"))
      rbind(
        do.call(rbind, lapply(options, function(o) {
          draws_or_refuses(sprintf("plot_profiles(statistic = %s, scale = %s, what = %s)",
                                   o$statistic, o$scale, o$what),
                           function() plot_profiles(object, statistic = o$statistic,
                                                    scale = o$scale, what = o$what))
        })),
        must_draw("plot_profiles(intervals = FALSE, labels = FALSE)",
                  function() plot_profiles(object, intervals = FALSE, labels = FALSE)),
        must_draw("plot_profiles(main, subtitle)",
                  function() plot_profiles(object, main = "Profiles", subtitle = "two")),
        sweep_refuses("an unknown statistic", function() plot_profiles(object, statistic = "max")),
        sweep_refuses("an unknown scale", function() plot_profiles(object, scale = "log")),
        sweep_refuses("an unknown layout", function() plot_profiles(object, what = "radar")))
    }),
    sweep_define("plot_diagnostics-arguments", function() {
      object <- fits$two_level()
      views <- eval(formals(plot_diagnostics)$what)
      rbind(
        do.call(rbind, lapply(views, function(view) {
          must_draw(sprintf("plot_diagnostics(what = \"%s\")", view),
                    function() plot_diagnostics(object, what = view))
        })),
        sweep_refuses("an unknown diagnostic", function() plot_diagnostics(object, what = "roc")))
    }),
    sweep_define("plot_responses-arguments", function() {
      object <- fits$classes()
      categories <- levels(items[[1L]])
      rbind(
        do.call(rbind, lapply(c("last", categories), function(category) {
          draws_or_refuses(sprintf("plot_responses(category = \"%s\")", category),
                           function() plot_responses(object, category = category))
        })),
        must_draw("plot_responses()", function() plot_responses(object)),
        must_draw("plot_responses(intervals = FALSE, labels = FALSE)",
                  function() plot_responses(object, intervals = FALSE, labels = FALSE)))
    }),
    sweep_define("grid-verbs-arguments", function() {
      grid <- fits$grid()
      criteria <- eval(formals(plot_enumeration)$what)
      selections <- eval(formals(plot_selection)$what)
      rbind(
        do.call(rbind, lapply(criteria, function(criterion) {
          must_draw(sprintf("plot_enumeration(what = \"%s\")", criterion),
                    function() plot_enumeration(grid, what = criterion))
        })),
        must_draw("plot_enumeration(combine = FALSE, labels = FALSE, mark_minimum = FALSE)",
                  function() plot_enumeration(grid, combine = FALSE, labels = FALSE,
                                              mark_minimum = FALSE)),
        do.call(rbind, lapply(selections, function(view) {
          draws_or_refuses(sprintf("plot_selection(what = \"%s\")", view),
                           function() plot_selection(grid, what = view))
        })),
        must_draw("plot_tree()", function() plot_tree(grid)),
        sweep_true("plot_views() lists the grid's views",
                   is.data.frame(plot_views(grid)) && nrow(plot_views(grid)) > 0L),
        sweep_refuses("an unknown criterion", function() plot_enumeration(grid, what = "dic")),
        sweep_refuses("plot_views() of a fit", function() plot_views(fits$profiles())))
    }),
    sweep_define("numbers-behind-the-views", function() {
      object <- fits$profiles()
      drawn <- ggplot2::ggplot_build(plot(object, what = "profiles"))$data
      points <- Find(function(layer) all(c("x", "y") %in% names(layer)) &&
                       !"label" %in% names(layer) && !"ymin" %in% names(layer), drawn)
      means <- sweep_table(object, "profiles")$mean
      sizes <- ggplot2::ggplot_build(plot(object, what = "sizes"))$data[[1L]]
      counts <- subset(sweep_table(object, "counts"), level == "individuals")$effective_count
      # The bars may run along either axis; the other holds the class positions.
      lengths <- if (max(sizes$x) > max(sizes$y)) sizes$x else sizes$y
      rbind(
        sweep_check("the profile plot's points are the table's means", sort(points$y),
                    sort(means), 1e-10),
        sweep_check("the size plot's bars are the table's effective counts",
                    sort(lengths), sort(counts), 1e-10))
    }))
  c(family_cases, argument_cases)
}
