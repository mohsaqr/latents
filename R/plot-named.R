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
#' @param view The `what =` value to draw.
#' @param arguments Named list of the arguments the caller supplied.
#' @param verb The calling function's name, for messages.
#' @return The method's ggplot object.
#' @noRd
.multilpa_plot_named_view <- function(x, view, arguments, verb) {
  method <- .multilpa_plot_method(x)
  views <- if (is.null(method)) character() else
    setdiff(eval(formals(method)$what), "all")
  if (!view %in% views) {
    stop(errorCondition(sprintf(
      "%s() cannot draw a `%s` fit, whose plot() offers %s.", verb,
      class(x)[1L], paste(sprintf("\"%s\"", views), collapse = ", ")),
      class = "latents_nothing_to_plot", call = NULL))
  }
  unused <- setdiff(names(arguments), names(formals(method)))
  if (length(unused) > 0L) {
    stop(errorCondition(sprintf(
      "A `%s` fit's %s() does not use %s.", class(x)[1L], verb,
      paste(sprintf("`%s`", unused), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  do.call(method, c(list(x, what = view), arguments))
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

#' Refuse an object the named fit views cannot draw, with a classed error
#'
#' The named views draw profile, class and transition models. Every other
#' result of the package has views of its own under `plot()`, and a general
#' [lta()] fit draws its transition matrix, so the refusal says where to go.
#' @param x The object passed to the named plot function.
#' @param verb The calling function's name, for the message.
#' @return `NULL`, invisibly, when `x` is a fit the named views take.
#' @noRd
.gg_require_fit <- function(x, verb) {
  if (.multilpa_any_fit(x)) return(invisible(NULL))
  advice <- if (inherits(x, "multilpa_lta")) {
    " This lta() fit draws its transition matrix: use plot_transitions()."
  } else if (any(grepl("^(latents|multilpa)_", class(x)))) {
    " Use plot() for this result's own views."
  } else ""
  stop(errorCondition(sprintf(
    "%s() takes a profile, class or transition model of this package, not a `%s` object.%s",
    verb, class(x)[1L], advice), class = "latents_bad_argument", call = NULL))
}

#' Refuse a grid other than a profile or class enumeration
#' @param x The object passed to `plot_enumeration()` or `plot_tree()`.
#' @param verb The calling function's name, for the message.
#' @return `NULL`, invisibly, for a `multilpa_enumeration`.
#' @noRd
.gg_require_enumeration <- function(x, verb) {
  if (inherits(x, "multilpa_enumeration")) return(invisible(NULL))
  advice <- if (inherits(x, c("latents_transition_enumeration",
                              "latents_family_enumeration"))) {
    paste(" A transition or group-class family grid plots with",
          "`plot(x, what = \"bic\")`, `\"aic\"` or `\"bic_individual\"`.")
  } else ""
  stop(errorCondition(sprintf(paste(
    "%s() takes a `multilpa_enumeration` from enumerate_classes(),",
    "enumerate_lpa() or enumerate_lca(), not a `%s` object.%s"),
    verb, class(x)[1L], advice), class = "latents_bad_argument", call = NULL))
}

#' Draw several views of a fit, naming the views that do not apply
#'
#' A grouped plot function's `what = "all"` draws each of its views with the
#' arguments that view takes. A view the fit has no ingredients for
#' (`.multilpa_supported_views()`), or that refuses with a classed condition,
#' is left out and named in one message, as `plot(x, what = "all")` does.
#' @param x A fitted model.
#' @param views Named list: the plot() view for each name, and in `uses` the
#'   arguments that view takes.
#' @param arguments Named list of the arguments the caller supplied.
#' @param verb The calling function's name, for messages.
#' @return A `latents_plots` list named by view, or one ggplot.
#' @noRd
.multilpa_plot_named_views <- function(x, views, arguments, verb) {
  refused <- \(condition) NULL
  supported <- .multilpa_supported_views(x)
  plots <- lapply(views, \(view) {
    if (!view$what %in% supported) return(NULL)
    tryCatch(.multilpa_plot_named_view(
      x, view$what, arguments[intersect(names(arguments), view$uses)], verb),
      latents_no_plot = refused, latents_nothing_to_plot = refused,
      latents_no_time = refused, latents_no_categorical = refused,
      latents_no_continuous = refused, latents_no_indicator_data = refused)
  })
  drawn <- !vapply(plots, is.null, logical(1))
  if (!any(drawn)) {
    stop(errorCondition(sprintf("%s() has no view this model can draw.", verb),
                        class = "latents_nothing_to_plot", call = NULL))
  }
  if (!all(drawn)) {
    message(sprintf("Not drawn for this model: %s.",
                    paste(names(views)[!drawn], collapse = ", ")))
  }
  .gg_plots(plots[drawn])
}

#' Draw one view of a grouped plot function, or all of them
#' @param x A fitted model.
#' @param what The caller's `what`, already matched.
#' @param views As in `.multilpa_plot_named_views()`.
#' @param arguments Named list of the arguments the caller supplied.
#' @param verb The calling function's name, for messages.
#' @return A ggplot object, or a `latents_plots` list for `"all"`.
#' @noRd
.multilpa_plot_grouped <- function(x, what, views, arguments, verb) {
  arguments$what <- NULL
  if (identical(what, "all")) {
    return(.multilpa_plot_named_views(x, views, arguments, verb))
  }
  view <- views[[what]]
  .multilpa_plot_named_view(
    x, view$what, arguments[intersect(names(arguments), view$uses)], verb)
}

#' @rdname latents-plots
#' @export
plot_profiles <- function(x, statistic = c("mean", "median", "mode"),
                          scale = c("raw", "standardized"), intervals = TRUE,
                          labels = TRUE, data = NULL, main = NULL,
                          subtitle = NULL,
                          what = c("lines", "bars", "heatmap", "raincloud",
                                   "parallel", "pairs", "all")) {
  .gg_require_fit(x, "plot_profiles")
  what <- match.arg(what)
  statistic <- match.arg(statistic)
  scale <- match.arg(scale)
  titles <- c("main", "subtitle")
  .multilpa_plot_grouped(x, what, list(
    lines = list(what = "profiles", uses = c("statistic", "scale", "intervals",
                                             "labels", "data", titles)),
    bars = list(what = "bars", uses = c("scale", "intervals", "data", titles)),
    heatmap = list(what = "heatmap", uses = titles),
    raincloud = list(what = "raincloud", uses = c("scale", titles)),
    parallel = list(what = "parallel", uses = titles),
    pairs = list(what = "pairs", uses = titles)),
    .multilpa_supplied_arguments(match.call(), environment()), "plot_profiles")
}

#' @rdname latents-plots
#' @export
plot_responses <- function(x, category = "last", intervals = TRUE,
                           labels = TRUE, data = NULL, main = NULL,
                           subtitle = NULL) {
  .gg_require_fit(x, "plot_responses")
  .multilpa_plot_named_view(
    x, "responses", .multilpa_supplied_arguments(match.call(), environment()),
    "plot_responses")
}

#' @rdname latents-plots
#' @export
plot_sequences <- function(x, labels = TRUE, cell_labels = TRUE, main = NULL,
                           subtitle = NULL) {
  .gg_require_fit(x, "plot_sequences")
  .multilpa_plot_named_view(
    x, "sequences", .multilpa_supplied_arguments(match.call(), environment()),
    "plot_sequences")
}

#' @rdname latents-plots
#' @export
plot_prevalence <- function(x, labels = TRUE, main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_prevalence")
  .multilpa_plot_named_view(
    x, "probabilities",
    .multilpa_supplied_arguments(match.call(), environment()),
    "plot_prevalence")
}

#' @rdname latents-plots
#' @export
plot_diagnostics <- function(x, what = c("all", "sizes", "posteriors",
                                         "entropy", "avepp", "group_sizes",
                                         "group_posteriors", "group_entropy",
                                         "group_avepp"),
                             main = NULL, subtitle = NULL) {
  .gg_require_fit(x, "plot_diagnostics")
  what <- match.arg(what)
  .gg_require()
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
  if (inherits(x, "multilpa_lta")) {
    return(plot(x, main = main, subtitle = subtitle))
  }
  if (!inherits(x, "multilpa_transitions")) {
    stop(errorCondition(
      "plot_transitions() takes a transition fit from lta().",
      class = "latents_bad_argument", call = NULL))
  }
  plot(x, what = "transitions", main = main, subtitle = subtitle)
}

#' @rdname latents-plots
#' @export
plot_tree <- function(x, main = NULL, subtitle = NULL) {
  .gg_require_enumeration(x, "plot_tree")
  .gg_require()
  .gg_view_tree(x, main, subtitle)
}

#' @rdname latents-plots
#' @export
plot_selection <- function(x, what = c("all", "deltas", "heatmap", "smallest",
                                       "trust", "diagnostics"),
                           using = "bic_groups", cutoff = 0.8, main = NULL,
                           subtitle = NULL) {
  .gg_require_enumeration(x, "plot_selection")
  what <- match.arg(what)
  plot_enumeration(x, what = if (identical(what, "all")) "selection" else what,
                   using = using, cutoff = cutoff, main = main,
                   subtitle = subtitle)
}
