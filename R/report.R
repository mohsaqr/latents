#' Every classification diagnostic, in one call
#'
#' `summary()` reports what a model estimated. This reports whether to believe
#' it: how sharply the posterior separates the classes, how big each class
#' actually is, how confidently each unit was assigned, and whether the
#' within-profile independence the model assumes survives contact with the data.
#'
#' It gathers [entropy_table()], [classification_table()],
#' [average_posteriors()] and [bivariate_residuals()], which return four
#' differently shaped tables, so the result is an object with a `print()` method
#' and an `as.data.frame()` accessor rather than four tables printed at once.
#'
#' @param x A fitted model of this package.
#' @param data Optional. The data the model was fitted to. A fit carries the
#'   columns it was built from, so this is only needed to override them.
#' @param plots `TRUE` also draws the classification plots, as a side effect,
#'   before returning. Equivalent to calling `plot()` on the result.
#' @param by Passed to [bivariate_residuals()]: `"profile"`, the default,
#'   assesses each profile separately, `"overall"` pools them. It is the one
#'   argument of a gathered verb this function forwards, because it is the one
#'   that changes what a gathered table means rather than which fit it is taken
#'   from. The `level` the classification tables use is not an argument here: it
#'   follows from whether the fit has discrete group classes.
#' @param ... For `diagnostics()`, nothing further is accepted. An argument this
#'   function cannot forward raises an error of class `multilpa_bad_argument`
#'   naming it, rather than being dropped on the way to a table that then means
#'   something other than what was asked for. For `plot()`, style overrides, as
#'   in [plot.multilpa()].
#' @return An object of class `multilpa_diagnostics`: a list with the elements
#'   `entropy`, `classification`, `posteriors` and `residuals`, each a tidy
#'   base `data.frame` as its own verb returns it, plus `fit` for the plot
#'   method. Reach the tables with `as.data.frame(result, what = )`, never with
#'   `$`. `residuals` is `NULL` for a model family that has no bivariate
#'   residuals, and asking for that table raises `multilpa_no_group_classes`.
#'
#'   `as.data.frame()` returns one of those tables: `"entropy"` exactly as
#'   [entropy_table()] returns it, `"classification"` as
#'   [classification_table()] does, `"posteriors"` as [average_posteriors()]
#'   does, and `"residuals"` as [bivariate_residuals()] does.
#'   `print()` returns the object invisibly, having printed one line per
#'   diagnostic: relative entropy, smallest class, lowest average posterior and
#'   largest residual, at each level the fit has. `plot()` returns the object
#'   invisibly, having drawn the case-level entropy and posterior panels.
#' @seealso [descriptives()] for the before-the-fit counterpart, [summary()]
#'   for what the model estimated rather than whether to trust it.
#' @examples
#' activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
#' fit <- multilpa(course_engagement, vars = activity, id = "student",
#'                 n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 4, seed = 1)
#' quality <- diagnostics(fit)
#' quality
#' as.data.frame(quality, what = "classification")
#' as.data.frame(diagnostics(fit, by = "overall"), what = "residuals")
#' @export
diagnostics <- function(x, data = NULL, plots = FALSE,
                        by = c("profile", "overall"), ...) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x),
            "`plots` must be TRUE or FALSE" = isTRUE(plots) || isFALSE(plots))
  .multilpa_reject_extra_arguments(
    list(...), "diagnostics()",
    paste("The only argument of a gathered verb it forwards is `by`;",
          "call that verb directly for anything else."))
  by <- match.arg(by)
  data <- .multilpa_resolve_data(x, data)
  has_groups <- !is.null(x$group_posteriors)
  level <- if (has_groups) "both" else "individuals"
  result <- list(
    entropy = entropy_table(x),
    classification = classification_table(x, level = level),
    posteriors = average_posteriors(x, level = level),
    # Residuals need a discrete group-class model; a random-intercept fit
    # refuses rather than returning a table of a different meaning.
    residuals = tryCatch(bivariate_residuals(x, data, by = by),
                         multilpa_no_group_classes = function(condition) NULL),
    fit = x)
  class(result) <- "multilpa_diagnostics"
  if (isTRUE(plots)) plot(result)
  result
}

#' @rdname diagnostics
#' @param row.names,optional Passed to the base generic; `row.names` is applied
#'   to the returned table.
#' @param what Which table: `"entropy"`, `"classification"`, `"posteriors"` or
#'   `"residuals"`.
#' @export
as.data.frame.multilpa_diagnostics <- function(x, row.names = NULL,
                                               optional = FALSE,
                                               what = c("entropy",
                                                        "classification",
                                                        "posteriors",
                                                        "residuals"), ...) {
  stopifnot(inherits(x, "multilpa_diagnostics"))
  what <- match.arg(what)
  result <- x[[what]]
  if (is.null(result)) {
    stop(errorCondition(
      "This fit has no bivariate residuals: they need a discrete group-class model.",
      class = "multilpa_no_group_classes", call = NULL))
  }
  row.names(result) <- row.names
  result
}

#' @rdname diagnostics
#' @export
print.multilpa_diagnostics <- function(x, ...) {
  stopifnot(inherits(x, "multilpa_diagnostics"))
  fit <- x$fit
  cat(sprintf("Classification quality: %d profiles", fit$n_profiles))
  if (!is.null(fit$n_group_classes)) cat(sprintf(", %d group classes", fit$n_group_classes))
  cat("\n\n")
  # `split()` orders alphabetically, which would print groups before
  # individuals on one line and after it on the next. The reporting order is
  # fixed once here so every line reads the same way.
  entropy <- x$entropy
  levels_in_order <- entropy$level
  line <- function(label, values) {
    cat(sprintf("  %-20s %s\n", label, paste(values, collapse = "    ")))
  }
  per_level <- function(frame, f) {
    parts <- lapply(levels_in_order, function(level) {
      f(frame[frame$level == level, , drop = FALSE], level)
    })
    unlist(parts, use.names = FALSE)
  }
  line("Relative entropy",
       sprintf("%s %.3f", entropy$level, entropy$relative_entropy))
  line("Smallest class", per_level(x$classification, function(rows, level) {
    row <- rows[which.min(rows$estimated_proportion), , drop = FALSE]
    sprintf("%s %.1f (%.1f%%)", level, row$estimated_n,
            100 * row$estimated_proportion)
  }))
  assigned <- x$posteriors[x$posteriors$assigned_class == x$posteriors$class, , drop = FALSE]
  line("Lowest avg posterior", per_level(assigned, function(rows, level) {
    sprintf("%s %.3f", level, min(rows$average_posterior))
  }))
  if (is.null(x$residuals)) {
    line("Largest residual", "not available for this model family")
  } else {
    worst <- x$residuals[which.max(abs(x$residuals$residual)), , drop = FALSE]
    line("Largest residual",
         sprintf("%.3f  (%s, %s; %s)", worst$residual, worst$indicator_1,
                 worst$indicator_2, worst$profile))
  }
  cat("\nTables: as.data.frame(x, what = \"entropy\" | \"classification\" |",
      "\n        \"posteriors\" | \"residuals\").  plot(x) draws them.\n")
  invisible(x)
}

#' Refuse arguments a forwarding verb cannot pass on
#'
#' A verb that gathers other verbs documents which of their arguments it
#' forwards. Anything else is named back to the caller instead of being dropped:
#' a discarded `by =` returns a table that answers a different question from the
#' one that was asked, and does so without a word.
#'
#' @param extra `list(...)` as the forwarding verb received it.
#' @param verb The verb's name, for the message.
#' @param advice One sentence saying what the verb does forward.
#' @return `NULL`, invisibly, when there is nothing to refuse.
#' @noRd
.multilpa_reject_extra_arguments <- function(extra, verb, advice) {
  if (length(extra) == 0L) return(invisible(NULL))
  supplied <- names(extra)
  if (is.null(supplied)) supplied <- rep("", length(extra))
  labels <- ifelse(nzchar(supplied), sprintf("`%s`", supplied), "an unnamed argument")
  stop(errorCondition(
    sprintf("%s does not use %s. %s", verb, paste(labels, collapse = ", "), advice),
    class = "multilpa_bad_argument", call = NULL))
}

#' @rdname diagnostics
#' @export
plot.multilpa_diagnostics <- function(x, ...) {
  stopifnot(inherits(x, "multilpa_diagnostics"))
  # Drawn from the fit's posteriors directly rather than through
  # `plot(fit, what = )`. The two panels need only the individual posteriors and
  # the effective profile counts, which every family carries, while not every
  # family's `plot()` method offers a `what` -- a covariate fit's takes
  # "profiles" and "sequences", a random-intercept fit's "profiles" and
  # "random_intercepts" -- so dispatching through the generic asked those fits
  # for a view their method had never heard of and failed on the match.
  .multilpa_plot_case_diagnostic(x$fit, what = "entropy", ...)
  .multilpa_plot_case_diagnostic(x$fit, what = "posteriors", ...)
  invisible(x)
}

#' Everything about a fit, in one call
#'
#' Prints the estimates, then the classification diagnostics, then draws every
#' plot the fit supports. A first look, not a substitute for the verbs: each
#' section is what [summary()], [diagnostics()] and [plot()] return, and the
#' tables are reachable from those rather than from here.
#'
#' @param x A fitted model of this package.
#' @param data Optional. The data the model was fitted to; a fit carries the
#'   columns it was built from.
#' @param plots `TRUE`, the default, draws the plots. `FALSE` prints only.
#' @param by Passed to [diagnostics()], and from there to
#'   [bivariate_residuals()]: `"profile"`, the default, assesses each profile
#'   separately, `"overall"` pools them.
#' @param ... Nothing further is accepted. An argument this function cannot
#'   forward raises an error of class `multilpa_bad_argument` naming it, before
#'   anything has been printed, rather than being dropped.
#' @return The fitted model, invisibly. Called for the printing and drawing.
#' @seealso [summary()], [diagnostics()], [descriptives()],
#'   [multilpa_plot_types()].
#' @examples
#' activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
#' fit <- multilpa(course_engagement, vars = activity, id = "student",
#'                 n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 4, seed = 1)
#' report(fit, plots = FALSE)
#' @export
report <- function(x, data = NULL, plots = TRUE,
                   by = c("profile", "overall"), ...) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x),
            "`plots` must be TRUE or FALSE" = isTRUE(plots) || isFALSE(plots))
  # Checked here as well as in diagnostics(), so a refused argument is refused
  # before three sections have already been printed.
  .multilpa_reject_extra_arguments(
    list(...), "report()",
    "It forwards `by` to diagnostics() and accepts nothing else.")
  by <- match.arg(by)
  data <- .multilpa_resolve_data(x, data)
  print(summary(x))
  cat("\n")
  print(descriptives(x))
  cat("\n")
  quality <- diagnostics(x, data = data, plots = FALSE, by = by)
  print(quality)
  if (isTRUE(plots)) {
    views <- .multilpa_supported_views(x)
    # A view can still refuse for a reason only the method knows -- a family
    # with no plot at all, a panel with nothing in it. `report()` is a first
    # look, so it moves on; it says which views it could not draw rather than
    # dropping them silently.
    refused <- unlist(lapply(views, function(view) {
      drawn <- tryCatch({
        if (is.na(view)) plot(x) else plot(x, what = view)
        NULL
      },
      multilpa_no_plot = function(condition) view,
      multilpa_nothing_to_plot = function(condition) view,
      multilpa_no_time = function(condition) view,
      multilpa_no_categorical = function(condition) view,
      multilpa_no_indicator_data = function(condition) view)
      drawn
    }))
    if (length(refused) > 0L) {
      cat(sprintf("\nNot drawn for this model: %s\n",
                  paste(ifelse(is.na(refused), "the default plot", refused),
                        collapse = ", ")))
    }
  }
  invisible(x)
}

#' The plot views a given fit can actually draw
#'
#' `multilpa_plot_types()` lists every view the methods accept; this filters it
#' to the ones this fit has the ingredients for, so `report()` does not stop on
#' the first refusal.
#'
#' @param x A fitted model of this package.
#' @return A character vector of `what` values.
#' @noRd
.multilpa_supported_views <- function(x) {
  # Ask the method that will be called what it accepts, rather than assuming
  # every class takes the full catalogue. A random-intercept fit's plot method
  # offers "profiles" and "random_intercepts" only, so report() must not hand it
  # a view it has never heard of.
  method <- utils::getS3method("plot", class(x)[1L], optional = TRUE)
  if (is.null(method)) return(character())
  # A method with no `what` draws one thing. `NA_character_` stands for "call
  # plot() with no view", which `report()` passes through as a bare plot(x).
  if (!"what" %in% names(formals(method))) return(NA_character_)
  candidates <- eval(formals(method)$what)
  Filter(function(view) {
    if (identical(view, "responses")) return(length(x$categorical) > 0L)
    if (identical(view, "sequences")) return(!is.null(x$time))
    if (view %in% c("profiles", "bars", "heatmap")) {
      return(!is.null(x$means) && ncol(x$means) > 0L)
    }
    if (identical(view, "probabilities")) return(!is.null(x$profile_probabilities))
    # The two case-level diagnostics separate cases by profile, so a
    # single-profile fit has nothing for them to separate: every case sits in
    # the one profile with probability one.
    if (view %in% c("entropy", "posteriors")) {
      return(!is.null(x$subject_posteriors) && ncol(x$subject_posteriors) > 1L)
    }
    TRUE
  }, candidates)
}

#' The supported views of a fit, phrased for a refusal message
#'
#' @param x A fitted model of this package.
#' @return A single string listing the views, or `"none"`.
#' @noRd
.multilpa_available_views_text <- function(x) {
  views <- .multilpa_supported_views(x)
  if (length(views) == 0L) return("none")
  paste(ifelse(is.na(views), "the default plot", views), collapse = ", ")
}
