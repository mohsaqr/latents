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

#' @rdname latents-plots
#' @export
plot_profiles <- function(x, statistic = c("mean", "median", "mode"),
                          scale = c("raw", "standardized"), intervals = TRUE,
                          labels = TRUE, data = NULL, main = NULL,
                          subtitle = NULL) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  statistic <- match.arg(statistic)
  scale <- match.arg(scale)
  .multilpa_plot_named_view(
    x, "profiles", .multilpa_supplied_arguments(match.call(), environment()),
    "plot_profiles")
}

#' @rdname latents-plots
#' @export
plot_responses <- function(x, category = "last", intervals = TRUE,
                           labels = TRUE, data = NULL, main = NULL,
                           subtitle = NULL) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  .multilpa_plot_named_view(
    x, "responses", .multilpa_supplied_arguments(match.call(), environment()),
    "plot_responses")
}

#' @rdname latents-plots
#' @export
plot_prevalence <- function(x, labels = TRUE, main = NULL, subtitle = NULL) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  .multilpa_plot_named_view(
    x, "probabilities",
    .multilpa_supplied_arguments(match.call(), environment()),
    "plot_prevalence")
}

#' @rdname latents-plots
#' @export
plot_diagnostics <- function(x, what = c("all", "sizes", "posteriors",
                                         "entropy", "avepp"),
                             main = NULL, subtitle = NULL) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  what <- match.arg(what)
  .gg_require()
  # Profile, class and transition fits carry posteriors and effective counts.
  views <- list(sizes = .gg_view_sizes, posteriors = .gg_view_posteriors,
                entropy = .gg_view_entropy, avepp = .gg_view_avepp)
  if (!identical(what, "all")) {
    if (!identical(what, "sizes")) {
      .multilpa_refuse_noise(x, sprintf("plot_diagnostics(what = \"%s\")", what))
    }
    return(views[[what]](x, main, subtitle))
  }
  if (isTRUE(x$noise)) views <- views["sizes"] else
    if (ncol(x$subject_posteriors) < 2L) views <- views[c("sizes", "avepp")]
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
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"))
  .gg_require()
  .gg_view_tree(x, main, subtitle)
}
