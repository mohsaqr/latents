#' Draw one view of a fit through its own plot method
#'
#' The named plot functions list their arguments so an editor completes them;
#' the drawing stays in the fit's `plot()` method, so a named function and
#' `plot(fit, what = )` cannot draw different things. Only the arguments the
#' caller supplied are passed on, so each method's own defaults apply. An
#' argument the fit's method does not take is refused by name rather than
#' dropped.
#'
#' @param x A fitted model.
#' @param view The `what =` value to draw, or `NA` for a method that draws a
#'   single view and takes no `what`.
#' @param arguments Named list of the arguments the caller supplied.
#' @param verb The calling function's name, for messages.
#' @param what The named function's own `what` (`""` for a function without
#'   one), for messages.
#' @return The method's ggplot object.
#' @noRd
.multilpa_plot_named_view <- function(x, view, arguments, verb, what = "") {
  method <- .multilpa_plot_method(x)
  if (!is.na(view) && !view %in% .latents_plot_offered_views(x)) {
    .latents_refuse_named_plot(x, verb, what)
  }
  unused <- setdiff(names(arguments), names(formals(method)))
  if (length(unused) > 0L) {
    stop(errorCondition(sprintf(
      "A `%s` fit's %s() does not use %s.", class(x)[1L], verb,
      paste(sprintf("`%s`", unused), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  do.call(method, c(list(x), if (!is.na(view)) list(what = view), arguments))
}

#' The kind of fitted model, for the named plot functions
#'
#' Profile, class and transition models share one set of views. The other
#' fitted models each have views of their own; a mixture regression that
#' classifies persons (`class_level = "group"`) draws trajectories, as a
#' growth mixture does, so it is a kind of its own.
#' @param x Any object.
#' @return One string: `"profile_model"`, `"growth_mixture"`,
#'   `"trajectory_regression"`, `"mixture_regression"`,
#'   `"group_class_family"`, `"cross_level"`, `"pooled"` or
#'   `"general_transitions"`; `NA` for an object that is not a fitted model.
#' @noRd
.latents_plot_kind <- function(x) {
  if (.multilpa_any_fit(x)) return("profile_model")
  if (inherits(x, "latents_growth_mixture")) return("growth_mixture")
  if (inherits(x, "latents_mixture_regression")) {
    return(if (identical(x$spec$nesting, "group")) "trajectory_regression" else
      "mixture_regression")
  }
  kinds <- c(multilpa_additive = "group_class_family",
             multilpa_cross_level = "cross_level", latents_pooled = "pooled",
             multilpa_lta = "general_transitions")
  matched <- kinds[intersect(class(x), names(kinds))]
  if (length(matched) == 0L) NA_character_ else unname(matched[1L])
}

#' Which view of a fit's own plot() each named plot function draws
#'
#' The one table behind the named fit plot functions. Each row says that the
#' named function `verb`, given `what`, draws `plot(x, what = view)` for a
#' fit of that `kind` (see `.latents_plot_kind()`). A function, kind and
#' `what` with no row are refused with `.latents_refuse_named_plot()`; a
#' mapping is listed only where the view shows the same quantity, so most
#' kinds map few functions. Profile, class and transition models draw
#' `plot_diagnostics()` from their posteriors at both levels in that function
#' itself, so they have no `plot_diagnostics` rows.
#' @return A `data.frame`, one row per mapped function, kind and `what`:
#'   `kind`, `verb`, `what` (`""` for a function without `what`) and `view`
#'   (the plot() view; `NA` when that plot() draws a single view and takes no
#'   `what`).
#' @noRd
.latents_named_plot_registry <- function() {
  rows <- function(kind, verb, what, view) {
    data.frame(kind = kind, verb = verb, what = what, view = view,
               stringsAsFactors = FALSE)
  }
  regressions <- c("growth_mixture", "trajectory_regression",
                   "mixture_regression")
  rbind(
    rows("profile_model", "plot_profiles",
         c("lines", "bars", "heatmap", "raincloud", "parallel", "pairs"),
         c("profiles", "bars", "heatmap", "raincloud", "parallel", "pairs")),
    rows("profile_model", c("plot_responses", "plot_sequences",
                            "plot_prevalence", "plot_transitions"), "",
         c("responses", "sequences", "probabilities", "transitions")),
    # The class's curve of the outcome: over time for trajectories, over the
    # first numeric predictor for a regression of observations.
    rows(regressions, "plot_profiles", "lines",
         c("trajectories", "trajectories", "fitted")),
    # Each unit's posterior probability of its most likely class.
    rows(regressions, "plot_diagnostics", "posteriors", "posteriors"),
    rows("group_class_family", "plot_profiles", "lines", "means"),
    rows("cross_level", c("plot_profiles", "plot_prevalence"), c("lines", ""),
         c("profiles", "composition")),
    # The pooled profile means and variances ("measurement"), and the pooled
    # probability of each profile within each group class ("profile"): the
    # prevalence that plot_prevalence() draws for a profile model.
    rows("pooled", c("plot_profiles", "plot_prevalence"), c("lines", ""),
         c("measurement", "profile")),
    rows("general_transitions", "plot_transitions", "", NA_character_))
}

#' The view a named plot function draws for a fit, or `NULL`
#' @param x A fitted model.
#' @param verb The named function.
#' @param what Its `what` (`""` for a function without one).
#' @return The plot() view (a string, `NA` for a single-view method), or
#'   `NULL` when the table maps none.
#' @noRd
.latents_named_plot_lookup <- function(x, verb, what = "") {
  registry <- .latents_named_plot_registry()
  row <- which(registry$kind %in% .latents_plot_kind(x) &
                 registry$verb == verb & registry$what == what)
  if (length(row) == 0L) NULL else registry$view[row[1L]]
}

#' The views a fit's own plot() offers
#' @param x A fitted model.
#' @return A character vector of `what` values; empty for a method that
#'   draws a single view.
#' @noRd
.latents_plot_offered_views <- function(x) {
  if (inherits(x, "latents_mixture_regression")) return(.mixture_plot_views(x))
  method <- .multilpa_plot_method(x)
  arguments <- if (is.null(method)) list() else formals(method)
  choices <- arguments$what
  # Methods that keep their 0.9.21 argument order hold the views in the
  # released argument (`level` for pooled fits, `criterion` for family and
  # transition grids) and default `what` to NULL.
  if ("what" %in% names(arguments) && is.null(choices)) {
    choices <- arguments$level %||% arguments$criterion
  }
  if (is.null(choices)) character() else setdiff(eval(choices), "all")
}

#' Refuse a named plot function a fit has no view for, naming its views
#'
#' The message lists the views of the fit's own plot() and, for the model
#' kinds other than profile, class and transition models, which named
#' functions draw them.
#' @param x A fitted model.
#' @param verb The named function.
#' @param what Its `what` (`""` for a function without one, `"all"` when no
#'   view of the function is mapped).
#' @return Does not return; raises `latents_nothing_to_plot` (also of class
#'   `latents_bad_argument`, which these refusals had before).
#' @noRd
.latents_refuse_named_plot <- function(x, verb, what = "") {
  asked <- function(verb, what) {
    ifelse(nzchar(what), sprintf("%s(what = \"%s\")", verb, what),
           sprintf("%s()", verb))
  }
  offered <- .latents_plot_offered_views(x)
  offered_text <- if (length(offered) == 0L) "its plot() draws a single view" else
    sprintf("its plot() offers %s",
            paste(sprintf("\"%s\"", offered), collapse = ", "))
  kind <- .latents_plot_kind(x)
  registry <- .latents_named_plot_registry()
  mapped <- if (identical(kind, "profile_model")) registry[0L, ] else
    registry[registry$kind %in% kind, , drop = FALSE]
  named_text <- if (nrow(mapped) == 0L) "" else
    sprintf("; of the named plot functions, %s", paste(
      sprintf("%s draws %s", asked(mapped$verb, mapped$what),
              ifelse(is.na(mapped$view), "it", sprintf("\"%s\"", mapped$view))),
      collapse = " and "))
  stop(errorCondition(sprintf(
    "%s has no view of a `%s` fit: %s%s.", asked(verb, what), class(x)[1L],
    offered_text, named_text),
    class = c("latents_nothing_to_plot", "latents_bad_argument"), call = NULL))
}

#' Draw the view a named plot function maps to, refusing an unmapped one
#'
#' Profile, class and transition models pass each view only the arguments
#' that view uses (`uses`), as the grouped functions always have. The other
#' kinds pass every supplied argument, so the fit's plot() refuses by name
#' any it does not take.
#' @param x A fitted model.
#' @param verb The named function.
#' @param what Its `what` (`""` for a function without one).
#' @param arguments Named list of the arguments the caller supplied.
#' @param uses The arguments this view uses, or `NULL` for all.
#' @param warn_unused `TRUE` warns (`latents_unused_argument`) about supplied
#'   arguments outside `uses`.
#' @return The method's ggplot object.
#' @noRd
.latents_draw_named <- function(x, verb, what, arguments, uses = NULL,
                                warn_unused = FALSE) {
  view <- .latents_named_plot_lookup(x, verb, what)
  if (is.null(view)) .latents_refuse_named_plot(x, verb, what)
  if (identical(.latents_plot_kind(x), "profile_model") && !is.null(uses)) {
    # A view asked for by name says which supplied arguments it leaves unused;
    # under what = "all" each view takes its share without a word.
    if (isTRUE(warn_unused)) {
      .latents_warn_unused_arguments(
        arguments[setdiff(names(arguments), uses)],
        sprintf("%s(what = \"%s\")", verb, what),
        "Those arguments belong to the function's other views.")
    }
    arguments <- arguments[intersect(names(arguments), uses)]
  }
  .multilpa_plot_named_view(x, view, arguments, verb, what)
}

#' The registered plot method selected by an object's S3 class chain
#' @param x An object to plot.
#' @return A function, or `NULL` when no class has a plot method.
#' @noRd
.multilpa_plot_method <- function(x) {
  methods <- lapply(class(x), function(class_name) {
    utils::getS3method("plot", class_name, optional = TRUE)
  })
  methods <- Filter(Negate(is.null), methods)
  if (length(methods) == 0L) NULL else methods[[1L]]
}

#' The arguments a named plot function's caller supplied
#'
#' @param call The function's `match.call()`.
#' @param env The function's environment, holding the matched values.
#' @return A named list of the supplied arguments other than `x`.
#' @noRd
.multilpa_supplied_arguments <- function(call, env) {
  supplied <- setdiff(names(call)[-1L], "x")
  mget(supplied, envir = env)
}

#' Refuse an object that is not a fitted model, with a classed error
#'
#' The named fit plot functions take every fitted model of the package (see
#' `.latents_plot_kind()`); which of them draws a view of which fit is in
#' `.latents_named_plot_registry()`. Grids, comparisons, bootstrap tests and
#' diagnostics have views of their own, so the refusal says where to go.
#' @param x The object passed to the named plot function.
#' @param verb The calling function's name, for the message.
#' @return `NULL`, invisibly, when `x` is a fitted model.
#' @noRd
.gg_require_fit <- function(x, verb, expression = NULL) {
  if (!is.na(.latents_plot_kind(x))) return(invisible(NULL))
  label <- .latents_argument_label(expression, "x")
  advice <- if (inherits(x, .latents_enumeration_classes)) {
    " Plot a grid with plot_enumeration(), plot_selection() or plot_tree()."
  } else if (inherits(x, "latents_sensitivity")) {
    sprintf(" Draw it with plot(%s).", label)
  } else if (is.data.frame(x)) {
    sprintf(paste(" The plot is drawn from the fitted model itself, not from its data",
                  "or a table taken from it: pass the fit, as in %s(fit)."), verb)
  } else if (any(grepl("^(latents|multilpa)_", class(x)))) {
    " Use plot() for this result's own views."
  } else ""
  stop(errorCondition(sprintf("%s() needs a fitted model, but `%s` is %s.%s",
                              verb, label, .latents_describe(x), advice),
                      class = c("latents_not_a_fit", "latents_bad_argument"), call = NULL))
}

#' Draw several views of a fit, naming the views that do not apply
#'
#' A grouped plot function's `what = "all"` draws each of its views that
#' `.latents_named_plot_registry()` maps for the fit, with the arguments that
#' view takes. A view the table does not map, that the fit has no ingredients
#' for (`.multilpa_supported_views()`), or that refuses with a classed
#' condition, is left out and named in one message, as
#' `plot(x, what = "all")` does. A fit for which the table maps none of the
#' function's views is refused, naming the views its plot() offers.
#' @param x A fitted model.
#' @param uses Named list, one element per `what` of the function: the
#'   arguments that view takes.
#' @param arguments Named list of the arguments the caller supplied.
#' @param verb The calling function's name, for messages.
#' @return A `latents_plots` list named by view.
#' @noRd
.multilpa_plot_named_views <- function(x, uses, arguments, verb) {
  views <- lapply(names(uses), \(what) .latents_named_plot_lookup(x, verb, what))
  names(views) <- names(uses)
  if (all(vapply(views, is.null, logical(1)))) {
    .latents_refuse_named_plot(x, verb, "all")
  }
  refused <- \(condition) NULL
  profile_model <- identical(.latents_plot_kind(x), "profile_model")
  supported <- if (profile_model) .multilpa_supported_views(x) else
    .latents_plot_offered_views(x)
  plots <- lapply(names(uses), \(what) {
    view <- views[[what]]
    if (is.null(view) || !view %in% supported) return(NULL)
    tryCatch(.latents_draw_named(x, verb, what, arguments, uses[[what]]),
      latents_no_plot = refused, latents_nothing_to_plot = refused,
      latents_no_time = refused, latents_no_categorical = refused,
      latents_no_continuous = refused, latents_no_indicator_data = refused)
  })
  names(plots) <- names(uses)
  drawn <- !vapply(plots, is.null, logical(1))
  if (!any(drawn)) {
    stop(errorCondition(sprintf("%s() has no view this model can draw.", verb),
                        class = "latents_nothing_to_plot", call = NULL))
  }
  if (!all(drawn)) {
    message(sprintf("Not drawn for this model: %s.",
                    paste(names(uses)[!drawn], collapse = ", ")))
  }
  .gg_plots(plots[drawn])
}

#' Draw one view of a grouped plot function, or all of them
#' @param x A fitted model.
#' @param what The caller's `what`, already matched.
#' @param uses As in `.multilpa_plot_named_views()`.
#' @param arguments Named list of the arguments the caller supplied.
#' @param verb The calling function's name, for messages.
#' @return A ggplot object, or a `latents_plots` list for `"all"`.
#' @noRd
.multilpa_plot_grouped <- function(x, what, uses, arguments, verb) {
  arguments$what <- NULL
  if (identical(what, "all")) {
    return(.multilpa_plot_named_views(x, uses, arguments, verb))
  }
  .latents_draw_named(x, verb, what, arguments, uses[[what]], warn_unused = TRUE)
}

#' @rdname latents-plots
#' @export
plot_profiles <- function(x, statistic = c("mean", "median", "mode"),
                          scale = c("raw", "standardized"), intervals = TRUE,
                          labels = TRUE, data = NULL, main = NULL,
                          subtitle = NULL,
                          what = c("lines", "bars", "heatmap", "raincloud",
                                   "parallel", "pairs", "all")) {
  .gg_require_fit(x, "plot_profiles", substitute(x))
  what <- match.arg(what)
  statistic <- match.arg(statistic)
  scale <- match.arg(scale)
  titles <- c("main", "subtitle")
  # The arguments each view takes; which plot() view each one draws is in
  # .latents_named_plot_registry().
  .multilpa_plot_grouped(x, what, list(
    lines = c("statistic", "scale", "intervals", "labels", "data", titles),
    bars = c("scale", "intervals", "data", titles),
    heatmap = titles,
    raincloud = c("scale", titles),
    parallel = titles,
    pairs = titles),
    .multilpa_supplied_arguments(match.call(), environment()), "plot_profiles")
}

#' @rdname latents-plots
#' @export
plot_responses <- function(x, category = "last", intervals = TRUE,
                           labels = TRUE, data = NULL, main = NULL,
                           subtitle = NULL) {
  .gg_require_fit(x, "plot_responses", substitute(x))
  .latents_draw_named(x, "plot_responses", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' @rdname latents-plots
#' @export
plot_sequences <- function(x, labels = TRUE, cell_labels = TRUE, main = NULL,
                           subtitle = NULL) {
  .gg_require_fit(x, "plot_sequences", substitute(x))
  .latents_draw_named(x, "plot_sequences", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' @rdname latents-plots
#' @export
plot_prevalence <- function(x, labels = TRUE, main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_prevalence", substitute(x))
  .latents_draw_named(x, "plot_prevalence", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' @rdname latents-plots
#' @export
plot_diagnostics <- function(x, what = c("all", "sizes", "posteriors",
                                         "entropy", "avepp", "group_sizes",
                                         "group_posteriors", "group_entropy",
                                         "group_avepp"),
                             main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_diagnostics", substitute(x))
  what <- match.arg(what)
  .gg_require()
  # The other fitted models draw their posteriors in their own plot(), as
  # .latents_named_plot_registry() maps them.
  if (!identical(.latents_plot_kind(x), "profile_model")) {
    individual <- c("sizes", "posteriors", "entropy", "avepp")
    return(.multilpa_plot_grouped(
      x, what, stats::setNames(rep(list(c("main", "subtitle")), 4L), individual),
      .multilpa_supplied_arguments(match.call(), environment()),
      "plot_diagnostics"))
  }
  # Profile, class and transition fits carry posteriors and effective counts;
  # a fit with two or more group classes carries them for its groups too, and
  # each view draws either level from the same code.
  views <- list(sizes = .gg_view_sizes, posteriors = .gg_view_posteriors,
                entropy = .gg_view_entropy, avepp = .gg_view_avepp)
  group_views <- lapply(views, \(view) \(x, main, subtitle)
                        view(x, main, subtitle, level = "groups"))
  names(group_views) <- paste0("group_", names(views))
  if (!identical(what, "all")) {
    if (what %in% names(group_views)) {
      .multilpa_require_case_posteriors(x, what, level = "groups")
      return(group_views[[what]](x, main, subtitle))
    }
    if (!identical(what, "sizes")) {
      .multilpa_refuse_noise(x, sprintf("plot_diagnostics(what = \"%s\")", what))
    }
    return(views[[what]](x, main, subtitle))
  }
  if (isTRUE(x$noise)) views <- views["sizes"] else
    if (ncol(x$subject_posteriors) < 2L) views <- views[c("sizes", "avepp")]
  if (.gg_has_group_classes(x) && !isTRUE(x$noise)) {
    views <- c(views, group_views)
  }
  .gg_plots(lapply(views, \(view) view(x, main, subtitle)))
}

#' @rdname latents-plots
#' @export
plot_transitions <- function(x, main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_transitions", substitute(x))
  .latents_draw_named(x, "plot_transitions", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' @rdname latents-plots
#' @export
plot_tree <- function(x, main = NULL, subtitle = NULL, n_group_classes = NULL,
                      n_profiles = NULL, model = NULL, candidates = NULL,
                      from = NULL, to = NULL,
                      level = c("individuals", "groups"), combine = FALSE) {
  .latents_require_grid_views(x, "tree", "plot_tree")
  # A family grid's candidates have group classes only, so its tree is theirs.
  if (identical(level, eval(formals(plot_tree)$level)) &&
      identical(.latents_grid_kind(x), "families")) level <- "groups"
  level <- .latents_check_tree_arguments(level, n_profiles, n_group_classes,
                                         model, combine)
  candidates <- .latents_tree_candidates(x, candidates, from, to, model)
  if (!is.null(from)) model <- NULL
  .gg_require()
  .gg_view_tree(x, level = level, n_profiles = n_profiles,
                n_group_classes = n_group_classes, model = model,
                candidates = candidates, combine = combine, main = main,
                subtitle = subtitle)
}

#' @rdname latents-plots
#' @export
plot_selection <- function(x, what = c("all", "deltas", "heatmap", "smallest",
                                       "trust", "diagnostics"),
                           using = "bic_groups", cutoff = 0.8, main = NULL,
                           subtitle = NULL) {
  views <- .latents_grid_views(x, "plot_selection")
  what <- match.arg(what)
  view <- if (identical(what, "all")) "selection" else what
  .latents_require_grid_views(x, view, "plot_selection")
  # Each kind of grid reads its own criterion unless `using` names one.
  plot_enumeration(x, what = view,
                   using = if (missing(using)) views$using else using,
                   cutoff = cutoff, main = main, subtitle = subtitle)
}
