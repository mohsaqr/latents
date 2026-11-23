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
    rows("general_transitions",
         c("plot_transitions", "plot_profiles", "plot_sequences"),
         c("", "lines", ""), c("transitions", "profiles", "sequences")))
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

#' Plot how the profiles differ on the indicators: each profile's mean, or
#' another summary, on every indicator
#'
#' `plot_profiles()` plots how the profiles (classes) of a fitted model differ
#' on its indicators: by default one line per profile through its mean of each
#' indicator, with 95% intervals. `what` chooses another view of them.
#'
#' @param x A fitted model, or pooled imputations. What each draws:
#'   * a profile or class fit ([lpa()], [lca()], [multilpa()], [multilca()]):
#'     every view it has the indicators for;
#'   * an [lta()] fit: `"lines"`, and `"bars"` and `"heatmap"` when its class
#'     is `multilpa_transitions`; nothing when it has only count indicators
#'     or `measurement = "occasion"`;
#'   * a growth mixture or mixture regression: `"lines"` only, each class's
#'     trajectory over time or fitted line over the first numeric predictor;
#'   * an additive or cross-level fit: `"lines"` only;
#'   * pooled imputations ([pool_imputations()]): `"lines"` only, the pooled
#'     means with their intervals beside each imputation's estimates.
#' @param statistic What `"lines"` joins. It matters only for a profile, class
#'   or `multilpa_lta` fit with only categorical or ordinal indicators:
#'   * `"mean"` (default): the probability-weighted category score;
#'   * `"median"`: the first category where cumulative probability reaches 0.5;
#'   * `"mode"`: the most probable category.
#' @param scale `"raw"` (default) draws each indicator in its own units;
#'   `"standardized"` in observed standard deviations from the observed mean
#'   (continuous indicators only). For `"lines"`, `"bars"` and `"raincloud"`
#'   of profile, class and `multilpa_transitions` fits.
#' @param intervals `TRUE` (default) draws the 95% intervals of
#'   [parameter_inference()] on `"lines"` and `"bars"` of continuous indicators
#'   (the subtitle says why a fit has none); `FALSE` draws none.
#' @param labels `TRUE` (default) writes each profile's name at the end of
#'   its line; `FALSE` uses a legend instead. For `what = "lines"`.
#' @param data The data frame the model was fitted to, used only to compute
#'   the intervals. `NULL` (default) uses the data stored in the fit.
#' @param main,subtitle Title and subtitle of the plot, each one string;
#'   `NULL` (default) keeps the view's own.
#' @param what Which view to draw. The linked pages say what a view needs:
#'   * `"lines"` (default): each profile's mean of every indicator, joined;
#'   * `"bars"`: the same means as bars ([plot_bars()]);
#'   * `"heatmap"`: the means in standard deviations ([plot_heatmap()]);
#'   * `"raincloud"`: each indicator's distribution ([plot_raincloud()]);
#'   * `"parallel"`: every case as a line, by profile ([plot_parallel()]);
#'   * `"pairs"`: the cases in every pair of indicators ([plot_pairs()]);
#'   * `"all"`: every view the fit supports, as a list.
#' @return A ggplot object; for `what = "all"` a `latents_plots` list of the
#'   views the fit can draw, named by view (a message names the others).
#'   Errors: `latents_nothing_to_plot` (the fit has no such view),
#'   `latents_no_continuous` (too few continuous indicators),
#'   `latents_bad_argument` (an argument no view of the fit takes, or
#'   `"lines"` for category labels that are not numbers), `latents_not_a_fit`.
#'   An argument of another view gives a `latents_unused_argument` warning.
#' @seealso [plot_responses()] for categorical indicators; [plot_all()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_profiles(fit)
#' plot_profiles(fit, scale = "standardized")
#' plot_profiles(fit, what = "bars")
#' @export
plot_profiles <- function(x, statistic = c("mean", "median", "mode"),
                          scale = c("raw", "standardized"), intervals = TRUE,
                          labels = TRUE, data = NULL, main = NULL,
                          subtitle = NULL,
                          what = c("lines", "bars", "heatmap", "raincloud",
                                   "parallel", "pairs", "all")) {
  .gg_require_fit(x, "plot_profiles", substitute(x))
  what <- .latents_match_arg(what, eval(formals(plot_profiles)$what), "what")
  statistic <- .latents_match_arg(statistic, eval(formals(plot_profiles)$statistic),
                                  "statistic")
  scale <- .latents_match_arg(scale, eval(formals(plot_profiles)$scale), "scale")
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

#' Plot each profile's probability of one response category on every categorical
#' or ordinal indicator
#'
#' `plot_responses()` plots how likely each profile (class) is to give one
#' response category. There is one line per profile, joining the probability of
#' the chosen category on every categorical or ordinal indicator, on a 0 to 1
#' axis. The categorical indicators carry 95% intervals, formed on the logit
#' scale, when the fit has standard errors.
#'
#' @inheritParams plot_profiles
#' @param x A profile, class or transition fit with categorical or ordinal
#'   indicators: from [lca()], [multilca()], [lpa()] or [multilpa()] (with
#'   or without membership covariates), or from [lta()] with categorical
#'   indicators. An [lta()] fit of class `multilpa_lta` (ordinal or count
#'   indicators, covariates, weights, a stayer class, second-order or
#'   occasion-specific transitions, occasion-specific measurement) is not
#'   drawn: read its probabilities with `get_results(x, "responses")` or
#'   `get_results(x, "ordinal")`.
#' @param category Which category's probability is drawn:
#'   * `"last"` (default): the last category of each indicator, which for a
#'     no/yes indicator is the probability of "yes";
#'   * `"first"`: the first category of each indicator;
#'   * one category label (such as `"yes"`) or one whole-number position
#'     (such as `2`), which every indicator must have.
#' @param intervals `TRUE` (default) draws the 95% intervals of
#'   [parameter_inference()]; a fit that has no standard errors is drawn
#'   without them, and the subtitle says why. A probability estimated at 0
#'   or 1 is held there and has no interval. `FALSE` draws no intervals and
#'   does not compute them.
#' @param labels `TRUE` (default) writes each profile's name at the end of
#'   its line; `FALSE` uses a legend instead.
#' @return A ggplot object, which draws the plot when printed. Conditions
#'   raised:
#'   * `latents_no_categorical`: a profile, class or transition fit without
#'     categorical or ordinal indicators;
#'   * `latents_unknown_category`: `category` is not a category of every
#'     indicator;
#'   * `latents_nothing_to_plot`: a model of another family, which has no
#'     response view; the message names the plot functions that apply;
#'   * `latents_not_a_fit`: `x` is not a fitted model.
#' @seealso [plot_profiles()] for continuous indicators; [plot_heatmap()]
#'   for every category at once.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
#'            n_classes = 2, n_starts = 2, seed = 1)
#' plot_responses(fit)
#' plot_responses(fit, category = "no")
#' @export
plot_responses <- function(x, category = "last", intervals = TRUE,
                           labels = TRUE, data = NULL, main = NULL,
                           subtitle = NULL) {
  .gg_require_fit(x, "plot_responses", substitute(x))
  .latents_draw_named(x, "plot_responses", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' Plot each unit's sequence of assigned profiles across the occasions as a row
#' of coloured cells
#'
#' `plot_sequences()` plots each group's (person's) sequence of profiles as a
#' row of coloured cells: one row per group, one cell per occasion, coloured by
#' the profile the observation at that occasion is assigned to (its most likely
#' one). Within a block the rows are sorted by each group's most frequent
#' profile, so groups with similar sequences sit together. An occasion at which
#' nothing was observed is shown grey, and one absent from the data is left
#' blank.
#'
#' @inheritParams plot_profiles
#' @param x A fit with ordered occasions:
#'   * a transition fit from [lta()], of either class
#'     (`multilpa_transitions` or `multilpa_lta`), except one with
#'     `measurement = "occasion"`, whose profiles differ by occasion;
#'   * a profile or class fit from [lpa()], [lca()], [multilpa()] or
#'     [multilca()] (with or without membership covariates) made with
#'     `time`.
#' @param labels `TRUE` (default) puts the rows in blocks under labelled
#'   strips: one block per group class when the fit has two or more group
#'   classes, otherwise one per most frequent profile. `FALSE` draws the
#'   rows in the same order as one block without strips.
#' @param cell_labels `TRUE` (default) prints the profile number in each
#'   cell when the grid has at most 60 rows and 40 occasions; a larger grid
#'   shows the profiles by colour, with a dash still marking an occasion that
#'   observed nothing when the grid also has noise or absent occasions.
#'   `FALSE` never prints the numbers.
#' @return A ggplot object, which draws the plot when printed. Conditions
#'   raised:
#'   * `latents_no_time`: a profile or class fit made without `time`;
#'   * `latents_nothing_to_plot`: a model with no sequence view (an
#'     [lta()] fit with `measurement = "occasion"`, a growth mixture, a
#'     mixture regression, an additive or cross-level fit, pooled
#'     imputations);
#'   * `latents_not_a_fit`: `x` is not a fitted model.
#' @seealso [plot_transitions()] for the transition matrix. The sequences as
#'   a table: `get_results(x, "sequences")` for a profile or class model and
#'   a transition fit of class `multilpa_transitions`;
#'   `get_results(x, "assignments")` (one row per observation, with its
#'   profile) for a general transition fit of class `multilpa_lta`.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
#'            id = "student", n_profiles = 2, time = "sequence",
#'            n_starts = 2, seed = 1)
#' plot_sequences(fit)
#' @export
plot_sequences <- function(x, labels = TRUE, cell_labels = TRUE, main = NULL,
                           subtitle = NULL) {
  .gg_require_fit(x, "plot_sequences", substitute(x))
  .latents_draw_named(x, "plot_sequences", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' Plot how common each profile is within each group class of a two-level model
#'
#' `plot_prevalence()` plots how common each profile is in each group class: the
#' quantity a two-level model adds to a one-level one. There is one line per
#' group class, labelled with the class's share of the groups, joining the
#' probability of each profile within that class.
#'
#' @inheritParams plot_profiles
#' @param x One of:
#'   * a profile or class fit from [lpa()], [lca()], [multilpa()] or
#'     [multilca()] without membership covariates: the profile
#'     probabilities within each group class. A fit with one group class
#'     (or none) is drawn as a single line, the overall profile shares;
#'   * a cross-level fit: the share of each profile among the members of
#'     each group class, as stacked bars;
#'   * pooled imputations ([pool_imputations()]): the pooled profile-level
#'     parameters with their intervals, beside each imputation's. They are
#'     the profile probabilities within each group class, or the membership
#'     coefficients of a model with profile covariates.
#' @param labels `TRUE` (default) writes each group class's name at the end
#'   of its line; `FALSE` uses a legend instead. For a profile or class
#'   fit; a cross-level fit or pooled imputations raise
#'   `latents_bad_argument` when it is given.
#' @return A ggplot object, which draws the plot when printed. Conditions
#'   raised:
#'   * `latents_nothing_to_plot`: a model with no such view. A fit with
#'     membership covariates has none because its prevalence varies with
#'     each unit's covariates; a transition fit from [lta()], a growth
#'     mixture, a mixture regression and an additive fit have none either.
#'     The message names the plot functions that apply;
#'   * `latents_not_a_fit`: `x` is not a fitted model.
#' @seealso [plot_sizes()] for the class sizes at each level;
#'   `get_results(x, "profile_probabilities")` for the table.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                 "student", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 2, seed = 1)
#' plot_prevalence(fit)
#' @export
plot_prevalence <- function(x, labels = TRUE, main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_prevalence", substitute(x))
  .latents_draw_named(x, "plot_prevalence", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' Plot how well a fitted model separates its classes: sizes, posteriors,
#' entropy and average posterior probabilities
#'
#' `plot_diagnostics()` plots how well a fitted model separates its classes. A
#' profile, class or transition fit has four views of how its observations are
#' classified into profiles: the class sizes, the posterior probability of each
#' observation's assigned profile, each observation's entropy and the average
#' posterior matrix. A fit with two or more group classes has the same four
#' views of how its groups are classified into group classes. Each view is also
#' a function of its own, which takes models of the other families too:
#' [plot_sizes()], [plot_posteriors()], [plot_entropy()] and [plot_avepp()].
#'
#' @inheritParams plot_profiles
#' @param x One of:
#'   * a profile or class fit from [lpa()], [lca()], [multilpa()] or
#'     [multilca()] (with or without membership covariates), or an [lta()]
#'     fit of class `multilpa_transitions`: every view;
#'   * a growth mixture or a mixture regression: `"posteriors"` only, each
#'     unit's posterior probability of its most likely class (`"all"` then
#'     holds that one view).
#'
#'   Any other model (an [lta()] fit of class `multilpa_lta`, an additive
#'   or cross-level fit, pooled imputations) raises
#'   `latents_nothing_to_plot`. [plot_sizes()], [plot_posteriors()],
#'   [plot_entropy()] and [plot_avepp()] draw those fits, pooled
#'   imputations excepted.
#' @param what Which view to draw:
#'   * `"all"` (default): every view the fit supports, as a list;
#'   * `"sizes"`: the effective number of observations in each profile (the
#'     sum of their posterior probabilities), as bars labelled with the
#'     count and its share;
#'   * `"posteriors"`: each observation's posterior probability of the
#'     profile it is assigned to, one strip per profile, with chance (one
#'     over the number of profiles) marked;
#'   * `"entropy"`: each observation's entropy relative to a flat
#'     posterior (0 is certain, 1 is flat), one strip per profile;
#'   * `"avepp"`: the average posterior probability matrix: for the
#'     observations assigned to each profile (rows), their mean posterior
#'     probability of every profile (columns). The diagonal is the
#'     classification accuracy;
#'   * `"group_sizes"`, `"group_posteriors"`, `"group_entropy"`,
#'     `"group_avepp"`: the same four views of the groups in their group
#'     classes.
#'
#'   `"all"` leaves out what the fit cannot draw: the group views of a fit
#'   with fewer than two group classes, `"posteriors"` and `"entropy"` of a
#'   fit with one profile, and every view but `"sizes"` of a fit with a
#'   noise component.
#' @return A ggplot object, which draws the plot when printed; for `"all"`
#'   a `latents_plots` list of ggplot objects named by view, which draws
#'   every plot when printed. Conditions raised:
#'   * `latents_no_group_classes`: a group view of a fit with fewer than
#'     two group classes;
#'   * `latents_unsupported_noise`: a view other than `"sizes"` of a fit
#'     with a noise component (such a fit always has one group class);
#'   * `latents_nothing_to_plot`: `"posteriors"` or `"entropy"` of a fit
#'     with one profile, or a model without the view asked for; the
#'     message names the views and plot functions that apply;
#'   * `latents_not_a_fit`: `x` is not a fitted model.
#' @seealso [diagnostics()] and [get_checks()] for the tables behind the
#'   views.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_diagnostics(fit, what = "avepp")
#' plot_diagnostics(fit, what = "sizes")
#' @export
plot_diagnostics <- function(x, what = c("all", "sizes", "posteriors",
                                         "entropy", "avepp", "group_sizes",
                                         "group_posteriors", "group_entropy",
                                         "group_avepp"),
                             main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_diagnostics", substitute(x))
  what <- .latents_match_arg(what, eval(formals(plot_diagnostics)$what), "what")
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

#' Plot the estimated probability of moving from each profile to each profile
#' between occasions
#'
#' `plot_transitions()` plots the estimated transition matrix of a latent
#' transition model as a grid of tiles: the profile at one occasion in the rows,
#' the profile at the next occasion in the columns, and in each tile the
#' probability of that move. Each row sums to one. There is one panel per group
#' class.
#'
#' @inheritParams plot_profiles
#' @param x A transition fit from [lta()], of either class:
#'   * `multilpa_transitions` (the basic transition model: continuous or
#'     categorical indicators and one first-order transition matrix for
#'     all occasions): the estimated matrix of each group class, the panel
#'     labelled with the class's share of the groups. A row with no
#'     expected moves out of its profile is not estimated from the data,
#'     and its label is put in parentheses;
#'   * `multilpa_lta` (any extension of the basic model: ordinal or count
#'     indicators, covariates, weights, a stayer class, second-order or
#'     occasion-specific transitions, occasion-specific measurement): the
#'     probabilities the model implies, averaged over the observed groups,
#'     with one panel per group class and, when the transitions differ by
#'     occasion, per occasion.
#' @return A ggplot object, which draws the plot when printed. Conditions
#'   raised:
#'   * `latents_nothing_to_plot`: a model without transitions, which is
#'     every model not fitted by [lta()], including a profile or class fit
#'     made with `time`; the message names the plot functions that apply;
#'   * `latents_not_a_fit`: `x` is not a fitted model.
#' @seealso [plot_sequences()] for the assigned profiles over time;
#'   `get_results(x, "transitions")` for the probabilities as a table;
#'   [get_tna()] hands the transitions to the tna package.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lta(course_engagement, c("browse", "lectures", "forum_read"),
#'            id = "student", n_profiles = 2, time = "sequence",
#'            n_starts = 2, seed = 1)
#' plot_transitions(fit)
#' @export
plot_transitions <- function(x, main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_transitions", substitute(x))
  .latents_draw_named(x, "plot_transitions", "",
                      .multilpa_supplied_arguments(match.call(), environment()))
}

#' Plot how the classes of each solution split or rearrange as one more class is
#' added
#'
#' `plot_tree()` plots a clustering tree of an enumeration grid, which shows
#' whether a solution with one more class splits an existing class or rearranges
#' them all. A tree has one row per solution and one bar per class, as wide as
#' the class's share of the units. A band between two rows shows where the units
#' of a class go in the next solution, as wide as the share of all units the two
#' classes have in common. Colour follows a class from row to row, and a new
#' colour marks a class that split off.
#'
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()], [enumerate_lca()] or [enumerate_regressions()]. A
#'   family grid (`enumerate_classes(family = )`) has group classes only,
#'   so its trees are of the groups.
#' @param main,subtitle Title and subtitle of the plot, each one string;
#'   `NULL` (default) keeps the view's own.
#' @param n_group_classes,n_profiles,model Which trees and rows to draw;
#'   `NULL` (default) keeps every converged candidate. At the individual
#'   level `n_group_classes` and `model` choose the trees and `n_profiles`
#'   the rows; at the group level `n_profiles` and `model` choose the trees
#'   and `n_group_classes` the rows. A single count for the rows compares
#'   that solution across the trees instead, as the rows of one tree:
#'   `plot_tree(x, n_profiles = 4, n_group_classes = c(2, 5))`.
#' @param candidates Two or more candidate names, drawn as the rows of one
#'   tree in the order given. Not with `n_group_classes`, `n_profiles` or
#'   `model`.
#' @param from,to Two solutions compared as a two-row tree, each named by its
#'   class counts `c(n_profiles, n_group_classes)` (one count means one group
#'   class) or by its candidate name:
#'   `plot_tree(x, from = c(3, 3), to = c(4, 3))`. Given together, not with
#'   `candidates`.
#' @param level Which units a tree divides:
#'   * `"individuals"` (default): the observations, into profiles (classes).
#'     A tree runs down the numbers of profiles at one number of group
#'     classes and one covariance structure;
#'   * `"groups"`: the groups, into group classes. A tree runs down the
#'     numbers of group classes at one number of profiles, so the grid needs
#'     two or more of them. The default for a family grid.
#' @param combine `FALSE` (default) returns one plot per tree; `TRUE` draws
#'   every tree as a panel of one plot.
#' @return A ggplot object when there is one tree or `combine = TRUE`;
#'   otherwise a `latents_plots` list of ggplot objects named by tree.
#'   Conditions raised:
#'   * `latents_no_tree`: a tree would have fewer than two rows;
#'   * `latents_unknown_candidate`: `candidates`, `from` or `to` names no
#'     converged candidate of the grid;
#'   * `latents_unsupported_view`: `level = "individuals"` on a family grid;
#'   * `latents_bad_argument`: `x` is not an enumeration grid, or the
#'     arguments are combined in a way the tree does not take.
#' @seealso [get_tree()] for the flows as a table; [plot_enumeration()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_tree(grid)
#' plot_tree(grid, from = 2, to = 3)
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

#' Plot the evidence for choosing the number of classes from an enumeration grid
#'
#' `plot_selection()` plots five views of an enumeration grid that bear on the
#' choice of a solution: how much a criterion improves with each added class,
#' the criterion over the whole grid, the smallest class, whether each estimate
#' can be trusted, and how sharply each candidate classifies. None of them
#' selects a model. Each view is also a function of its own, named below.
#'
#' @inheritParams plot_enumeration
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()], [enumerate_lca()] or [enumerate_regressions()]. A
#'   family grid (`enumerate_classes(family = )`) has no numbers of
#'   profiles to lay the views out by and draws `"diagnostics"` only.
#' @param what Which view to draw:
#'   * `"all"` (default): the five views below, as a list;
#'   * `"deltas"`: the decrease in the criterion `using` with each added
#'     class and each added group class ([plot_deltas()]);
#'   * `"heatmap"`: the criterion `using` over the numbers of classes and
#'     group classes, relative to its lowest value ([plot_heatmap()]);
#'   * `"smallest"`: the share of the smallest class of every candidate,
#'     with 5% drawn dashed ([plot_sizes()]);
#'   * `"trust"`: whether each candidate converged, whether it has an
#'     estimate at a bound, and how many starts reached its best log
#'     likelihood ([plot_trust()]);
#'   * `"diagnostics"`: the relative entropy and the lowest average
#'     posterior probability (AvePP) of every candidate, with `cutoff`
#'     drawn dashed ([plot_entropy()], [plot_avepp()]).
#' @param using The information criterion that `"deltas"` and `"heatmap"`
#'   read: one criterion column of the grid. When not given it is
#'   `"bic_groups"` on a profile or class grid and `"bic"` on a transition
#'   or regression grid.
#' @param cutoff The AvePP cutoff drawn dashed by `"diagnostics"`: one
#'   number above 0 and at most 1 (default 0.8).
#' @return A ggplot object, which draws the plot when printed; for `"all"`
#'   a `latents_plots` list of the five ggplot objects, named `deltas`,
#'   `heatmap`, `smallest`, `trust` and `diagnostics`, which draws every
#'   plot when printed. Conditions raised:
#'   * `latents_unsupported_view`: a view the kind of grid does not offer
#'     (every view but `"diagnostics"` on a family grid, `"all"`
#'     included); the message gives the reason and what the grid offers;
#'   * `latents_bad_argument`: `x` is not an enumeration grid, or `using`
#'     is not a criterion of the grid.
#' @seealso [get_deltas()], [get_knees()] and [get_best()] for the tables;
#'   [plot_enumeration()] for the information criteria; [plot_views()]
#'   lists the views a grid offers.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_selection(grid, what = "deltas")
#' plot_selection(grid, what = "trust")
#' @export
plot_selection <- function(x, what = c("all", "deltas", "heatmap", "smallest",
                                       "trust", "diagnostics"),
                           using = "bic_groups", cutoff = 0.8, main = NULL,
                           subtitle = NULL) {
  views <- .latents_grid_views(x, "plot_selection")
  what <- .latents_match_arg(what, eval(formals(plot_selection)$what), "what")
  view <- if (identical(what, "all")) "selection" else what
  .latents_require_grid_views(x, view, "plot_selection")
  # Each kind of grid reads its own criterion unless `using` names one.
  plot_enumeration(x, what = view,
                   using = if (missing(using)) views$using else using,
                   cutoff = cutoff, main = main, subtitle = subtitle)
}
