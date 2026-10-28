# Checks of what a caller passed, shared by every exported verb.
#
# A refusal is a classed condition with no call (so it never carries data),
# and its message names the argument, says what was received, and when there
# is one, the fix: the misspelt column and its likely spelling, the indicator
# with missing values and `missing = "fiml"`, the fitted result passed where
# its data belonged (review 2026-10-27).

#' The classes of the package's fitted models and of its enumerations
#' @noRd
.latents_model_classes <- c("multilpa", "multilpa_covariates", "multilpa_transitions",
                            "multilpa_lta", "multilpa_additive", "multilpa_cross_level",
                            "latents_growth_mixture", "latents_mixture_regression")
#' @noRd
.latents_enumeration_classes <- c("multilpa_enumeration", "latents_transition_enumeration",
                                  "latents_family_enumeration",
                                  "latents_regression_enumeration")
#' The families `.multilpa_any_fit()` accepts, in words
#'
#' A verb that takes only these refused the others as not being "a fitted
#' model of this package" while naming them as one (review 2026-10-27).
#' @noRd
.latents_profile_fit_text <- paste(
  "a profile, class or transition fit from lpa(), lca(), multilpa(),",
  "multilca() or lta() with time-homogeneous transitions")

#' Classed `stopifnot()`
#'
#' The same contract as named `stopifnot()`: each argument is evaluated in
#' order and the first that is not all `TRUE` stops, with its name as the
#' message (the deparsed expression when unnamed). The condition carries
#' `class` and no call.
#' @param ... Named logical expressions.
#' @param .class The condition class.
#' @return `NULL`, invisibly.
#' @noRd
.latents_require <- function(..., .class = "latents_bad_argument") {
  checks <- as.list(substitute(list(...)))[-1L]
  labels <- names(checks) %||% rep("", length(checks))
  frame <- parent.frame()
  # Sequential on purpose: a later check may assume an earlier one held.
  invisible(lapply(seq_along(checks), function(index) {
    passed <- eval(checks[[index]], frame)
    if (!(is.logical(passed) && length(passed) > 0L && !anyNA(passed) && all(passed))) {
      message <- if (nzchar(labels[[index]])) labels[[index]] else
        sprintf("`%s` is not TRUE.", deparse(checks[[index]], width.cutoff = 60L)[[1L]])
      stop(errorCondition(message, class = .class, call = NULL))
    }
  }))
  invisible(NULL)
}

#' Classed `match.arg()` for a single choice
#'
#' The contract of `match.arg(arg, choices)` with `several.ok = FALSE`: `NULL`
#' or the whole default vector gives the first choice, an exact name wins, and
#' a unique abbreviation names its choice. Anything else -- an unknown or
#' ambiguous name, `NA`, several names, a non-character value -- raises
#' `latents_bad_argument` listing the choices, where `match.arg()` raised an
#' unclassed error that carried its own internal call.
#'
#' @param arg The caller's value.
#' @param choices The accepted values, the first being the default.
#' @param argument The argument's name, for the message.
#' @return A single element of `choices`.
#' @noRd
.latents_match_arg <- function(arg, choices, argument = "what") {
  if (is.null(arg) || identical(arg, choices)) return(choices[[1L]])
  single <- is.character(arg) && length(arg) == 1L && !is.na(arg)
  at <- if (single) pmatch(arg, choices) else NA_integer_
  if (!is.na(at)) return(choices[[at]])
  received <- if (single) sprintf("\"%s\"", arg) else .latents_describe(arg)
  stop(errorCondition(sprintf(
    "`%s` must be one of %s, not %s.%s", argument,
    paste(sprintf("\"%s\"", choices), collapse = ", "), received,
    if (single) .latents_did_you_mean(arg, choices) else ""),
    class = "latents_bad_argument", call = NULL))
}

#' What an object is, in a few words
#' @param x Any object.
#' @return One string, such as "a data frame with 200 rows and 4 columns".
#' @noRd
.latents_describe <- function(x) {
  results <- c(
    multilpa_enumeration = "an enumeration result (from enumerate_classes(), enumerate_lpa() or enumerate_lca())",
    latents_transition_enumeration = "an enumeration of transition models",
    latents_family_enumeration = "an enumeration across model families",
    latents_regression_enumeration = "an enumeration of mixture regressions (from enumerate_regressions())",
    multilpa_bootstrap_lrt = "a bootstrap likelihood ratio test (from bootstrap_lrt())",
    multilpa_diagnostics = "a diagnostics table (from diagnostics())",
    latents_comparison = "a model comparison (from compare_models())",
    latents_growth_mixture = "a fitted growth mixture",
    latents_mixture_regression = "a fitted mixture regression (from mixture_regression())",
    multilpa_transitions = "a fitted transition model (from lta())",
    multilpa_lta = "a fitted transition model (from lta())",
    multilpa_covariates = "a fitted model with membership covariates",
    multilpa_additive = "a fitted additive group-class model",
    multilpa_cross_level = "a fitted cross-level model",
    multilpa = "a fitted model (from lpa(), lca(), multilpa() or multilca())",
    latents_table = "a results table",
    latents_sensitivity = "a seed-sensitivity table (from sensitivity())")
  known <- names(results)[vapply(names(results), function(kind) inherits(x, kind), logical(1))]
  if (length(known) > 0L) return(results[[known[[1L]]]])
  if (is.null(x)) return("NULL")
  if (is.data.frame(x)) {
    return(sprintf("a data frame with %d %s and %d %s", nrow(x),
                   if (nrow(x) == 1L) "row" else "rows", ncol(x),
                   if (ncol(x) == 1L) "column" else "columns"))
  }
  if (is.matrix(x)) return(sprintf("a %s matrix [%d x %d]", typeof(x), nrow(x), ncol(x)))
  if (is.function(x)) return("a function")
  if (is.factor(x)) return(sprintf("a factor of length %d", length(x)))
  if (is.atomic(x)) {
    return(sprintf("a %s vector of length %d", class(x)[[1L]], length(x)))
  }
  sprintf("an object of class \"%s\"", class(x)[[1L]])
}

#' The name a caller gave an argument, when it is a plain name
#' @param expression The result of `substitute()` on the argument.
#' @param default The argument's own name.
#' @return The caller's name when `expression` is a symbol, else `default`.
#' @noRd
.latents_argument_label <- function(expression, default) {
  if (is.symbol(expression) && nzchar(as.character(expression))) {
    as.character(expression)
  } else default
}

#' Backquoted names joined for a message: `a`, `b` and `c`
#' @param names A character vector.
#' @return One string.
#' @noRd
.latents_name_list <- function(names) {
  quoted <- sprintf("`%s`", names)
  if (length(quoted) <= 1L) return(paste(quoted, collapse = ""))
  paste(paste(quoted[-length(quoted)], collapse = ", "), "and", quoted[[length(quoted)]])
}

#' A "did you mean" clause for misspelt names
#'
#' The closest available name by edit distance, offered only when it is at
#' most two edits (or a third of the name) away, or differs only in case.
#' @param absent The names that were not found.
#' @param available The names that exist.
#' @return A clause starting with a space, or "".
#' @noRd
.latents_did_you_mean <- function(absent, available) {
  if (length(absent) == 0L || length(available) == 0L) return("")
  suggestions <- vapply(absent, function(name) {
    same_case <- available[tolower(available) == tolower(name)]
    if (length(same_case) > 0L) return(same_case[[1L]])
    distance <- utils::adist(name, available)[1L, ]
    best <- which.min(distance)
    if (distance[[best]] <= max(2L, nchar(name) %/% 3L)) available[[best]] else NA_character_
  }, character(1))
  found <- !is.na(suggestions)
  if (!any(found)) return("")
  sprintf(" Did you mean %s?", .latents_name_list(unique(suggestions[found])))
}

#' Refuse anything but a data frame, saying what was passed instead
#'
#' A fitted result or an enumeration passed as `data` is the common case:
#' the message says so and asks for the data frame it was fitted on.
#' @param data The `data` argument.
#' @param expression `substitute(data)` in the verb, for the caller's name.
#' @return `NULL`, invisibly, or `latents_bad_data`.
#' @noRd
.latents_check_data_frame <- function(data, expression = NULL) {
  if (is.data.frame(data)) return(invisible(NULL))
  label <- .latents_argument_label(expression, "data")
  what <- .latents_describe(data)
  result <- inherits(data, c(.latents_model_classes, .latents_enumeration_classes,
                             "multilpa_bootstrap_lrt", "latents_comparison",
                             "multilpa_diagnostics"))
  message <- if (result) {
    sprintf(paste("`data` must be a data frame, but `%s` is %s.",
                  "Pass the data frame the model is to be fitted on."), label, what)
  } else {
    sprintf("`data` must be a data frame, but `%s` is %s.", label, what)
  }
  stop(errorCondition(message, class = "latents_bad_data", call = NULL))
}

#' Require a fitted result of this package, saying what was passed instead
#' @param x The object.
#' @param expression `substitute(x)` in the verb, for the caller's name.
#' @param accepted Classes accepted; `NULL` accepts any fitted model.
#' @param expected What is accepted, in words, for the message. When
#'   `accepted` is `NULL` (profile, class and homogeneous transition fits
#'   only), the default names those families, so a fit of another family is
#'   never told it must be "a fitted model of this package".
#' @param hint What to do instead, for a fitted model of a family the verb
#'   does not take; one sentence, or `NULL` for the general hints.
#' @return `NULL`, invisibly, or `latents_not_a_fit` (a `latents_bad_argument`),
#'   also for a fitted model of a family the verb does not take.
#' @noRd
.latents_require_fit <- function(x, expression = NULL, accepted = NULL,
                                 expected = NULL, hint = NULL) {
  ok <- if (is.null(accepted)) .multilpa_any_fit(x) else inherits(x, accepted)
  if (isTRUE(ok)) return(invisible(NULL))
  expected <- expected %||% if (is.null(accepted)) .latents_profile_fit_text else
    "a fitted model of this package"
  label <- .latents_argument_label(expression, "x")
  wants_enumeration <- length(accepted) > 0L && all(accepted %in% .latents_enumeration_classes)
  other_family <- !is.null(hint) && inherits(x, .latents_model_classes)
  hint <- if (other_family) {
    paste0(" ", hint)
  } else if (is.data.frame(x) && wants_enumeration) {
    " Enumerate models on the data first, for example with enumerate_lpa(), and pass the result."
  } else if (is.data.frame(x)) {
    " Fit a model to the data first, for example with lpa() or lca(), and pass the fit."
  } else if (wants_enumeration && inherits(x, .latents_model_classes)) {
    " It is already a single model, so there is nothing to pick from."
  } else if (inherits(x, .latents_enumeration_classes)) {
    " Pick one of its models with get_candidate() and pass that."
  } else ""
  # A fitted model of a family the verb does not take is a fit: it is told
  # which it is and what the verb needs, never that it must be "a fitted model
  # of this package" (INF-2).
  generic <- identical(expected, "a fitted model of this package")
  message <- if (inherits(x, .latents_model_classes) && !nzchar(hint)) {
    sprintf("`%s` is %s, which this verb does not take; it needs %s.", label,
            .latents_describe(x), if (generic) "another kind of model" else expected)
  } else {
    sprintf("`%s` must be %s, but it is %s.%s", label, expected, .latents_describe(x), hint)
  }
  stop(errorCondition(message, class = c("latents_not_a_fit", "latents_bad_argument"),
                      call = NULL))
}

#' Refuse indicators with missing or non-finite values, naming them
#'
#' Infinite and NaN values are refused under every `missing` setting, since
#' FIML treats only NA as missing; NA values are refused under
#' `missing = "error"`, with the counts per indicator and the remedy.
#' @param values A matrix or data frame of indicator values, columns named.
#' @param missing The fit's `missing` setting.
#' @param fiml_available Whether `missing = "fiml"` is a remedy for this model.
#' @return `NULL`, invisibly, or `latents_bad_data`.
#' @noRd
.latents_refuse_missing_values <- function(values, missing = "error",
                                           fiml_available = TRUE) {
  values <- as.data.frame(values, optional = TRUE)
  non_finite <- vapply(names(values), function(name) {
    column <- values[[name]]
    if (is.numeric(column)) sum(is.infinite(column) | is.nan(column)) else 0L
  }, numeric(1))
  if (any(non_finite > 0)) {
    hit <- non_finite[non_finite > 0]
    stop(errorCondition(sprintf(paste(
      "%s %s infinite or NaN values (%s). Recode them as NA or remove those rows."),
      if (length(hit) == 1L) "Indicator" else "Indicators",
      paste(.latents_name_list(names(hit)), if (length(hit) == 1L) "has" else "have"),
      paste(sprintf("%s: %d", names(hit), as.integer(hit)), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  if (!identical(missing, "error")) return(invisible(NULL))
  absent <- vapply(values, function(column) sum(is.na(column)), numeric(1))
  if (!any(absent > 0)) return(invisible(NULL))
  hit <- absent[absent > 0]
  remedy <- if (fiml_available) {
    "Use `missing = \"fiml\"` to fit every observed value, or remove the incomplete rows."
  } else "This model needs complete data: remove the incomplete rows."
  stop(errorCondition(sprintf("%s %s missing values (%s). %s",
    if (length(hit) == 1L) "Indicator" else "Indicators",
    paste(.latents_name_list(names(hit)), if (length(hit) == 1L) "has" else "have"),
    paste(sprintf("%s: %d", names(hit), as.integer(hit)), collapse = ", "), remedy),
    class = "latents_bad_data", call = NULL))
}

#' Compact structure of a result of this package
#'
#' A fitted model is a list of some ninety components, and `str()` walked all
#' of them: a host that shows `str()` before `print()` when an object is
#' clicked (Beatrina, and RStudio's Environment pane) showed internals first.
#' The result's own header line comes first instead, with the way to its
#' tables and to the raw list.
#' @param object A result of this package.
#' @param ... Ignored.
#' @return `NULL`, invisibly.
#' @noRd
.latents_str <- function(object, ...) {
  header <- utils::capture.output(print(object))
  header <- header[nzchar(trimws(header))]
  cat(sprintf("<%s> %s\n", class(object)[[1L]],
              if (length(header) > 0L) trimws(header[[1L]]) else .latents_describe(object)))
  cat(sprintf(" A list of %d components. Tables: get_results(x, what = ); raw list: str(unclass(x))\n",
              length(unclass(object))))
  invisible(NULL)
}

#' @exportS3Method utils::str
str.multilpa <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_covariates <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_transitions <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_lta <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_additive <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_cross_level <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.latents_growth_mixture <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.latents_mixture_regression <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_enumeration <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.latents_transition_enumeration <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.latents_family_enumeration <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.latents_regression_enumeration <- function(object, ...) .latents_str(object, ...)
#' @exportS3Method utils::str
str.multilpa_diagnostics <- function(object, ...) .latents_str(object, ...)
