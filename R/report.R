#' Diagnose a fitted model: convergence, replication, class sizes, separation
#' and local independence
#'
#' @description
#' `diagnostics()` examines the estimation and interpretation of a fitted
#' mixture model. It assesses convergence and replication across starting
#' values, describes the size and separation of the latent classes, and,
#' where supported, examines associations between indicators that are not
#' adequately represented by the model.
#'
#' Convergence and likelihood replication indicate whether different starting
#' values lead to the same numerical solution. Class sizes and posterior
#' membership probabilities describe how observations are distributed across
#' classes and how uncertain their classification remains. Bivariate residuals
#' compare observed indicator associations with those implied by the fitted
#' model, helping identify departures from the measurement specification.
#'
#' The function brings these results together in diagnostic tables and, where
#' available, classification plots. For multilevel models, classification
#' summaries distinguish observation-level profiles from group-level classes.
#' Conditions such as incomplete convergence, estimates at parameter bounds,
#' and classes with no assigned observations are identified for further
#' examination.
#'
#' @section Interpretation:
#' Numerical convergence, classification separation, and measurement fit
#' address different aspects of a model. A converged solution may contain
#' poorly separated classes, while sharply separated classes may arise under
#' an inadequate measurement specification. The diagnostics should therefore
#' be interpreted jointly, with attention to the study design and the
#' substantive meaning of the estimated classes.
#'
#' @details
#' The checks are grouped in five domains:
#'
#' * **Estimation**: convergence of the selected solution; replication of the
#'   best log likelihood across completed starts; failed or
#'   unfinished starts; and estimates at parameter bounds. Preliminary
#'   screening starts are reported separately from full optimisation runs.
#' * **Class size**, at each level the fit classifies: the smallest class,
#'   the largest gap between modal and model-estimated class proportions, and
#'   classes to which no unit is assigned.
#' * **Classification**: relative entropy, which
#'   describes classification precision and is not a criterion for the number
#'   of classes; the lowest average posterior
#'   probability of the assigned class and the lowest odds of correct
#'   classification.
#' * **Bivariate residuals**: discrepancies between observed and model-implied
#'   indicator associations, reported under the `"Local independence"` domain.
#'   Where tests are available, pairs are flagged at 0.05 after
#'   Benjamini-Hochberg adjustment. Residuals
#'   without an available test are reported descriptively.
#' * **Information**: log likelihood, number of parameters and BIC, for
#'   comparison with other fitted models using [compare_models()].
#'
#' The status labels are `"fail"` for a non-converged selected solution or a
#' class with effective membership below one; `"check"` for conditions requiring
#' review, including failed or unfinished starts, active parameter bounds,
#' classes with no modal assignments or residual indicator dependence; `"ok"`
#' when the checked condition is absent; `"info"` for a descriptive quantity;
#' and `"n/a"` for an unavailable or inapplicable quantity. An `"ok"` status
#' applies only to the specified check and does not establish overall model
#' adequacy. Classification summaries use the fitted posterior probabilities,
#' eligible observations and sampling weights. The available diagnostics depend
#' on the model family.
#'
#' Entropy, classification and average-posterior tables are supplied at the
#' classification levels supported by each family, together with a summary of
#' replication across starts (`"replication"`).
#' Group-classified regressions report the group level; hierarchical regression
#' and growth models report both levels. All-empty transition rows follow the
#' fitted `drop_empty` policy, and zero-weight units are never classified.
#' Bivariate residuals are gathered for families that support them.
#'
#' @param x A fitted model: a profile or class fit ([lpa()], [lca()],
#'   [multilpa()], [multilca()], with or without membership covariates), a
#'   transition fit ([lta()]), an additive or cross-level fit, a mixture
#'   regression or a growth mixture ([mixture_regression()]). Anything else,
#'   including pooled imputations and enumeration grids, raises
#'   `latents_not_a_fit`.
#' @param data The data frame the model was fitted to, used for the
#'   bivariate residuals. `NULL` (default) uses the columns stored in the
#'   fit.
#' @param plots Whether to draw the classification plots as well:
#'   * `TRUE` (default): for a profile or class fit and an [lta()] fit of
#'     class `multilpa_transitions`, draws the views of
#'     [plot_diagnostics()]: effective profile sizes, the posterior
#'     probability of each observation's assigned profile, each
#'     observation's entropy and the average posterior probability matrix,
#'     and the same four for the groups of a fit with two or more group
#'     classes. For any other model a message says that it has no
#'     classification plots here, and the tables are returned without
#'     them. Without the ggplot2 package a message says so, and the tables
#'     are returned;
#'   * `FALSE`: returns the tables and draws nothing.
#' @param by How the bivariate residuals are assessed:
#'   * `"profile"` (default): within each profile separately, one row per
#'     profile and pair of indicators;
#'   * `"overall"`: pooled over the profiles, one row per pair of
#'     indicators.
#' @param ... Not used. Any further argument to `diagnostics()`, or to
#'   `plot()` of its result, raises `latents_bad_argument`; style the
#'   returned plots with ggplot2.
#' @return An object of class `multilpa_diagnostics`. It prints the checks
#'   grouped by domain, with those that need review marked. Its tables are
#'   data frames read with [get_results()]:
#'   * `"entropy"` (the default, also `as.data.frame()`): one row per
#'     level of classification: `level`, `n_classes`, `n_units`,
#'     `entropy_sum`, `relative_entropy`;
#'   * `"classification"`: one row per level and class: `level`, `class`,
#'     `n_modal`, `proportion_modal`, `estimated_n`,
#'     `estimated_proportion`, `average_posterior`,
#'     `odds_correct_classification`;
#'   * `"average_posteriors"`: for the units assigned to each class, their
#'     mean posterior probability of every class: `level`,
#'     `assigned_class`, `class`, `n_assigned`, `average_posterior`;
#'   * `"residuals"`: the bivariate residuals, one row per pair of
#'     indicators (and per profile under `by = "profile"`): `profile`,
#'     `indicator_1`, `indicator_2`, `kind`, `observed`, `expected`,
#'     `residual`, `effective_n`, `statistic`, `df`, `p_value`,
#'     `p_adjusted` (Benjamini-Hochberg);
#'   * `"replication"`: one row on the starts: `n_starts`, `n_screened`,
#'     `n_converged`, `n_failed`, `n_best_replicated`, `share`,
#'     `tolerance`, `best_start`, `log_likelihood`;
#'   * `"checks"`: one row per check: `domain`, `check`, `level`, `value`
#'     (numeric), `display` (the value formatted), `criterion`, `status`
#'     (`"ok"`, `"check"`, `"fail"`, `"info"` or `"n/a"`), `detail` and
#'     `reference`;
#'   * `"all"`: a named list of the tables above.
#'
#'   `"residuals"` exists for profile and class fits and for [lta()] fits
#'   of class `multilpa_transitions`. Any other model (an [lta()] fit of
#'   class `multilpa_lta`, an additive or cross-level fit, a mixture
#'   regression, a growth mixture) has no such table: it is left out of
#'   `"all"`, and asking for it by name raises `latents_bad_argument`. A
#'   fit with a noise component has none either, and asking for it raises
#'   `latents_unsupported_noise`.
#'
#'   `plot()` of the result is [plot_diagnostics()] of the fit: a
#'   `latents_plots` list of ggplot objects, which draws every plot when
#'   printed. `summary()` returns the result itself.
#' @seealso [get_results()] for these and the other tables, [get_checks()]
#'   for the checks alone, [plot_diagnostics()] for the plots,
#'   [descriptives()] for the variables before fitting, and [summary()] for
#'   the estimates.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' quality <- diagnostics(fit, plots = FALSE)
#' quality
#' get_results(quality, "checks")
#' get_results(quality, what = "classification")
#' pooled <- diagnostics(fit, by = "overall", plots = FALSE)
#' get_results(pooled, what = "residuals")
#' @export
diagnostics <- function(x, data = NULL, plots = TRUE,
                        by = c("profile", "overall"), ...) {
  .latents_require_fit(x, substitute(x), .latents_model_classes,
                       "a fitted model of this package")
  .latents_require("`plots` must be TRUE or FALSE" = isTRUE(plots) || isFALSE(plots))
  .multilpa_reject_extra_arguments(
    list(...), "diagnostics()",
    paste("The only argument of a gathered table it forwards is `by`;",
          "ask get_results() for that table directly for anything else."))
  by <- .latents_match_arg(by, c("profile", "overall"), "by")
  # A profile, class or transition fit has every classification table; the
  # other families have their own classification and fit tables, from which
  # the checks are read, and no bivariate residuals.
  classified <- .multilpa_any_fit(x)
  residual_reason <- "not defined for this model family"
  result <- if (classified) list() else .latents_family_diagnostic_tables(x)
  if (classified) {
    data <- .multilpa_resolve_data(x, data)
    level <- if (.multilpa_has_group_classes(x)) "both" else "individuals"
    result$entropy <- .multilpa_entropy_table(x)
    result$classification <- .multilpa_classification_table(x, level = level)
    result$average_posteriors <- .multilpa_average_posteriors(x, level = level)
    # Residuals need a discrete group-class model and no noise component; a
    # fit without them reports the check as not available, with the reason.
    # The check reads BH-adjusted p-values, as it states; without `adjust`
    # the table came back unadjusted (up to 1.1.9), and correctly specified
    # profile models were flagged in up to half of simulated data sets
    # (validation/bivariate-residual-size.R).
    result$residuals <- tryCatch(
      .multilpa_bivariate_residuals(x, data, by = by, adjust = "BH"),
      latents_no_group_classes = function(condition) {
        residual_reason <<- conditionMessage(condition)
        NULL
      },
      latents_unsupported_noise = function(condition) {
        residual_reason <<- "not available for a fit with a noise component"
        NULL
      })
  }
  result$replication <- .multilpa_replication_table(x)
  result$checks <- .latents_diagnostic_checks(x, result$residuals, residual_reason)
  result$fit <- x
  class(result) <- "multilpa_diagnostics"
  if (isTRUE(plots) && .gg_available()) {
    if (classified) {
      print(plot(result))
    } else {
      message("This model family has no classification plots; plot() of the fit draws its own views.")
    }
  }
  result
}

#' Does a fit classify its groups into discrete group classes?
#'
#' Every fit of the profile family stores group posteriors, including a
#' single-level fit (one row per observation) and a fit with one group class
#' (one column). Neither classifies anything at the group level: the "groups"
#' row of a classification table was a single class holding every group, with
#' an average posterior of one by construction (review 2026-10-27). Only a
#' group posterior with two or more columns is a group classification.
#' @param x A fitted model of this package.
#' @return A single logical.
#' @noRd
.multilpa_has_group_classes <- function(x) {
  !isTRUE(x$single_level) && !is.null(x$group_posteriors) &&
    NCOL(x$group_posteriors) > 1L
}

#' @noRd
#' @export
as.data.frame.multilpa_diagnostics <-function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  stopifnot(inherits(x, "multilpa_diagnostics"))
  .multilpa_reject_extra_arguments(
    list(...), "as.data.frame()",
    "It coerces to the primary table; get_results(x, what = ) has the others.")
  result <- get_results(x)
  row.names(result) <- row.names
  result
}

#' @noRd
#' @export
summary.multilpa_diagnostics <- function(object, ...) {
  stopifnot("`object` must be a `multilpa_diagnostics` object" =
              inherits(object, "multilpa_diagnostics"))
  .multilpa_reject_extra_arguments(list(...), "summary()", .multilpa_summary_advice())
  # The diagnostics are already the tabular reading of a fit: printing them is
  # their summary, and get_results() reads their tables. summary.default
  # listed the internal components, the fit among them, as a Length/Class/Mode
  # table instead (review 2026-10-27).
  object
}

#' @noRd
#' @export
print.multilpa_diagnostics <-function(x, ...) {
  stopifnot(inherits(x, "multilpa_diagnostics"))
  fit <- x$fit
  # A diagnostics object saved before the checks existed prints as it did.
  if (!is.null(x$checks)) {
    header <- utils::capture.output(print(fit))
    header <- header[nzchar(trimws(header))]
    cat(sprintf("Diagnostics: %s\n", sub(":", ",", trimws(header[[1L]]), fixed = TRUE)))
    .latents_print_checks(x$checks)
    cat("\nEvery check with its reference: get_results(x, \"checks\"). Other tables: get_results(x, what = ).\n")
    return(invisible(x))
  }
  cat(sprintf("Classification quality: %d profile%s", fit$n_profiles,
              if (identical(as.integer(fit$n_profiles), 1L)) "" else "s"))
  one_group_class <- is.null(fit$n_group_classes) ||
    identical(as.integer(fit$n_group_classes), 1L)
  if (!one_group_class) cat(sprintf(", %d group classes", fit$n_group_classes))
  cat("\n\n")
  cat(sprintf("  %-20s %s\n", "Best likelihood",
              .multilpa_replication_phrase(fit, subject = "")))
  # `split()` orders alphabetically, which would print groups before
  # individuals on one line and after it on the next. The reporting order is
  # fixed once here so every line reads the same way.
  entropy <- x$entropy
  # With one group class the group level classifies nothing, so only the
  # individual level is worth a column.
  if (one_group_class) entropy <- entropy[entropy$level == "individuals", , drop = FALSE]
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
  averages <- x$average_posteriors
  assigned <- averages[averages$assigned_class == averages$class, , drop = FALSE]
  line("Lowest avg posterior", per_level(assigned, function(rows, level) {
    # A class to which no unit is assigned has no average posterior; the
    # minimum is over the classes that have members, and the empty ones are
    # named rather than turning the whole level into NA.
    empty <- rows$class[is.na(rows$average_posterior)]
    if (length(empty) == nrow(rows)) return(sprintf("%s none assigned", level))
    lowest <- sprintf("%s %.3f", level, min(rows$average_posterior, na.rm = TRUE))
    if (length(empty) == 0L) return(lowest)
    sprintf("%s (class %s: no %s assigned)", lowest, paste(empty, collapse = ", "),
            if (identical(level, "groups")) "group" else "unit")
  }))
  if (is.null(x$residuals)) {
    line("Largest residual",
         if (isTRUE(fit$noise)) "not available for a fit with a noise component"
         else "not available for this model family")
  } else {
    worst <- x$residuals[which.max(abs(x$residuals$residual)), , drop = FALSE]
    line("Largest residual",
         sprintf("%.3f  (%s, %s; %s)", worst$residual, worst$indicator_1,
                 worst$indicator_2, worst$profile))
  }
  cat("\nTables: get_results(x, what = \"entropy\" | \"classification\" |",
      "\n        \"average_posteriors\" | \"residuals\" | \"replication\" | \"all\").",
      "plot(x) draws them.\n")
  invisible(x)
}

#' How often the starts reached the best log likelihood
#'
#' A start counts when it ran to the end and its log likelihood is within
#' `replication_tolerance` of the best (the rule of the fit's own count). A
#' best solution reached by one start of several is a warning sign: the
#' likelihood may have higher maxima no start found.
#' @param x A fitted model.
#' @return A one-row data frame: `n_starts` (run to the end), `n_converged`,
#'   `n_failed` (no log likelihood), `n_best_replicated`, `share`, `tolerance`,
#'   `best_start`, `log_likelihood`.
#' @noRd
.multilpa_replication_table <- function(x) {
  starts <- x$starts
  # The restricted cross-level family searches its two factors separately and
  # records one row per factor (`part`), not per start: those rows are not
  # starts, so neither the count of starts nor the converged and failed
  # counts can be read from them (audit 2026-11-01: "2 of 2 starts" was
  # printed for a fit of five starts per factor).
  per_part <- is.data.frame(starts) && "part" %in% names(starts)
  n_starts <- x$n_starts %||%
    (if (is.null(starts) || per_part) NA_integer_ else nrow(starts))
  replicated <- x$n_best_replicated %||% NA_integer_
  counted <- is.data.frame(starts) && !per_part &&
    all(c("log_likelihood", "converged") %in% names(starts))
  n_screened <- 0L
  if (counted && "stage" %in% names(starts)) {
    # A family that screens starts (short runs that choose which to finish)
    # records them beside the finished ones; only a finished start can have
    # converged, failed or reached the best, so only those are counted.
    # Stages differ by family ("completed", "final", "full", "failed"); only
    # "screened" marks a short probe that was never meant to finish.
    finished <- !starts$stage %in% "screened"
    n_screened <- sum(!finished)
    starts <- starts[finished, , drop = FALSE]
    n_starts <- nrow(starts)
  }
  n_failed <- if (counted) {
    sum(!is.finite(starts$log_likelihood) |
          (if ("stage" %in% names(starts)) starts$stage %in% "failed" else FALSE))
  } else NA_integer_
  n_converged <- if (counted) {
    sum(is.finite(starts$log_likelihood) & starts$converged %in% TRUE &
          (if ("stage" %in% names(starts)) !starts$stage %in% "failed" else TRUE))
  } else NA_integer_
  # The finished starts are counted from the start record itself, so the
  # counts never exceed the total.
  if (counted) n_starts <- nrow(starts)
  data.frame(n_starts = as.integer(n_starts), n_screened = as.integer(n_screened),
             n_converged = as.integer(n_converged),
             n_failed = as.integer(n_failed), n_best_replicated = as.integer(replicated),
             share = if (isTRUE(n_starts > 0)) replicated / n_starts else NA_real_,
             tolerance = x$replication_tolerance %||% NA_real_,
             best_start = as.integer(x$best_start %||% NA_integer_),
             log_likelihood = x$log_likelihood)
}

#' The replication count as a phrase for printed output
#' @param x A fitted model.
#' @return One string.
#' @noRd
.multilpa_replication_phrase <- function(x, subject = "best log likelihood ") {
  table <- .multilpa_replication_table(x)
  if (is.na(table$n_best_replicated) || is.na(table$n_starts)) {
    return(sprintf("best start: %s", format(table$best_start)))
  }
  phrase <- sprintf("%sreached by %d of %d starts", subject, table$n_best_replicated,
                    table$n_starts)
  notes <- character()
  if (!is.na(table$n_converged) && table$n_starts > 1L) {
    unconverged <- table$n_starts - table$n_converged - table$n_failed
    notes <- c(notes, sprintf("%d converged", table$n_converged),
               if (unconverged > 0L) sprintf("%d stopped at max_iter", unconverged),
               if (table$n_failed > 0L) sprintf("%d failed", table$n_failed))
  }
  if (identical(x$search, "emEM") && !is.null(x$n_screen)) {
    notes <- c(notes, sprintf("emEM: %d screened", as.integer(x$n_screen)))
  } else if (isTRUE(table$n_screened > 0L)) {
    notes <- c(notes, sprintf("%d more screened", table$n_screened))
  }
  if (table$n_best_replicated == 1L && table$n_starts > 1L) {
    notes <- c(notes, "found once: refit with more starts to confirm it")
  }
  if (length(notes) > 0L) phrase <- sprintf("%s (%s)", phrase, paste(notes, collapse = "; "))
  phrase
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
    class = "latents_bad_argument", call = NULL))
}

#' Warn about arguments a released method used to ignore
#'
#' Up to latents 0.9.21 these `plot()` methods ignored any argument they did not
#' use. Refusing one now would break code written for that release, and
#' ignoring it silently would hide a call that does not do what it says, so the
#' method draws as before and says what it ignored.
#'
#' @param extra `list(...)` as the method received it.
#' @param verb The call, for the message.
#' @param advice What to do instead.
#' @return `NULL`, invisibly; warns with class `latents_unused_argument`.
#' @noRd
.latents_warn_unused_arguments <- function(extra, verb, advice) {
  if (length(extra) == 0L) return(invisible(NULL))
  supplied <- names(extra) %||% rep("", length(extra))
  labels <- ifelse(nzchar(supplied), sprintf("`%s`", supplied),
                   "an unnamed argument")
  warning(warningCondition(
    sprintf("%s ignores %s. %s", verb, paste(labels, collapse = ", "), advice),
    class = "latents_unused_argument", call = NULL))
  invisible(NULL)
}

#' Take a value given under an argument's 0.9.21 name
#'
#' Arguments renamed after latents 0.9.21 keep their old name as a second
#' formal, so code written for that release runs unchanged. Giving both names
#' is ambiguous and refused.
#'
#' @param current The value of the argument under its present name.
#' @param current_given Whether the caller supplied the present name.
#' @param former The value under the 0.9.21 name, `NULL` when not given.
#' @param names The present and the former name, for the message.
#' @return The value to use.
#' @noRd
.latents_former_argument <- function(current, current_given, former, names) {
  if (is.null(former)) return(current)
  if (current_given) {
    stop(errorCondition(sprintf(paste(
      "`%s` and `%s` are the same argument (`%s` is its name in latents",
      "0.9.21); give one of them."), names[1L], names[2L], names[2L]),
      class = "latents_bad_argument", call = NULL))
  }
  former
}

#' @noRd
#' @export
plot.multilpa_diagnostics <-function(x, ...) {
  stopifnot(inherits(x, "multilpa_diagnostics"))
  .multilpa_reject_extra_arguments(
    list(...), "plot()",
    "Style the returned plots with ggplot2, for example `+ ggplot2::theme()`.")
  plot_diagnostics(x$fit)
}

#' Print the estimates, descriptive statistics and diagnostic checks of a fitted
#' model, and draw its plots
#'
#' @description
#' `report()` prints a first look at a fitted profile, class or transition
#' model in one call: its estimates, the descriptive statistics of its
#' variables and its diagnostic checks, followed by its plots. It computes
#' nothing of its own; each part is the output of another function, named
#' under *Report contents*, which returns the same numbers as tables.
#'
#' The output goes to the R console and the graphics device. The function
#' returns the fitted model invisibly, so it can end a pipeline.
#'
#' @param x A profile, class or transition fit: from [lpa()], [lca()],
#'   [multilpa()] or [multilca()] (with or without membership covariates),
#'   or an [lta()] fit of class `multilpa_transitions` (the basic
#'   transition model: continuous or categorical indicators and one
#'   first-order transition matrix for all occasions, without covariates,
#'   weights or a stayer class). Any other model (an [lta()] fit of class
#'   `multilpa_lta`, a mixture regression, a growth mixture, an additive or
#'   cross-level fit) raises `latents_not_a_fit`; its own `summary()`,
#'   [diagnostics()] and `plot()` report it.
#' @param data The data frame the model was fitted to. `NULL` (default)
#'   uses the columns stored in the fit.
#' @param plots Whether to draw the plots:
#'   * `TRUE` (default): after the printed output, draws every plot view
#'     the fit has the ingredients for, one after another (see
#'     [plot_views()]), and then names the views that could not be drawn.
#'     Without the ggplot2 package a message says so and nothing is drawn;
#'   * `FALSE`: prints only.
#' @param rows How many rows of each table of `summary(x)` to print: one
#'   whole number of 0 or more (default 10). A longer table ends with the
#'   number of rows left out and the call that returns it whole.
#' @param by How [diagnostics()] assesses the bivariate residuals:
#'   * `"profile"` (default): within each profile separately;
#'   * `"overall"`: pooled over the profiles.
#' @param ... Not used. Any further argument raises `latents_bad_argument`
#'   before anything is printed.
#' @section Report contents:
#' The parts are printed in this order:
#'
#' 1. `summary(x)`: every table of the fit, each cut to `rows` rows;
#' 2. [descriptives()] of the fit: the number of observations, missing
#'    values, mean, standard deviation and range of each variable, and
#'    for grouped data its intraclass correlation;
#' 3. the checks of [diagnostics()]: estimation, class sizes,
#'    classification, bivariate residuals, and the log likelihood with its
#'    BIC, each with its status;
#' 4. with `plots = TRUE`, the plot views of the fit.
#'
#' The tables behind each part are returned by [get_results()],
#' [descriptives()] and [diagnostics()].
#' @return The fitted model `x`, invisibly. The function is called for its
#'   printed and plotted output. Conditions raised:
#'   * `latents_not_a_fit`: `x` is not a profile, class or
#'     `multilpa_transitions` fit;
#'   * `latents_bad_argument`: an argument `report()` does not take, or a
#'     `plots` or `by` with a value it does not take.
#' @seealso [summary()], [diagnostics()], [descriptives()], [plot_views()];
#'   [plot_all()] for the plots alone.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' report(fit, plots = FALSE)
#' @export
report <- function(x, data = NULL, plots = TRUE,
                   by = c("profile", "overall"), rows = 10L, ...) {
  .latents_require_fit(x, substitute(x), hint = paste(
    "report() reads the profile, class and transition tables; summary(),",
    "diagnostics() and plot() of the fit report this family."))
  .latents_require("`plots` must be TRUE or FALSE" = isTRUE(plots) || isFALSE(plots))
  # Checked here as well as in diagnostics(), so a refused argument is refused
  # before three sections have already been printed.
  .multilpa_reject_extra_arguments(
    list(...), "report()",
    "It forwards `by` to diagnostics() and `rows` to print(summary(x)).")
  by <- .latents_match_arg(by, c("profile", "overall"), "by")
  data <- .multilpa_resolve_data(x, data)
  print(summary(x), rows = rows)
  cat("\n")
  print(descriptives(x))
  cat("\n")
  quality <- diagnostics(x, data = data, plots = FALSE, by = by)
  print(quality)
  if (isTRUE(plots) && .gg_available()) {
    views <- .multilpa_supported_views(x)
    # A view can still refuse for a reason only the method knows -- a family
    # with no plot at all, a panel with nothing in it. `report()` is a first
    # look, so it moves on; it says which views it could not draw rather than
    # dropping them silently.
    refused <- unlist(lapply(views, function(view) {
      drawn <- tryCatch({
        print(if (is.na(view)) plot(x) else plot(x, what = view))
        NULL
      },
      latents_no_plot = function(condition) view,
      latents_nothing_to_plot = function(condition) view,
      latents_no_time = function(condition) view,
      latents_no_categorical = function(condition) view,
      latents_no_continuous = function(condition) view,
      latents_no_indicator_data = function(condition) view)
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
#' `plot_views()` lists every view the methods accept; this filters it
#' to the ones this fit has the ingredients for, so `report()` does not stop on
#' the first refusal.
#'
#' @param x A fitted model of this package.
#' @return A character vector of `what` values.
#' @noRd
.multilpa_supported_views <- function(x) {
  # A general transition fit keeps its measurement per occasion block. With
  # one block its profile and sequence views draw from that block, read here
  # as the plot method reads it; with one block per occasion only the
  # transition matrix has a single view.
  if (inherits(x, "multilpa_lta")) {
    if (length(x$measurement) != 1L) return("transitions")
    x <- .lta_plot_view(x)
  }
  # Ask the method that will be called what it accepts, rather than assuming
  # every class takes the full catalogue. A covariate fit's plot method offers
  # "profiles" and "sequences" only, so report() must not hand it a view it has
  # never heard of.
  method <- .multilpa_plot_method(x)
  if (is.null(method)) return(character())
  # A method with no `what` draws one thing. `NA_character_` stands for "call
  # plot() with no view", which `report()` passes through as a bare plot(x).
  if (!"what" %in% names(formals(method))) return(NA_character_)
  # `"all"` is a request to draw every view, not a view, so it is dropped
  # before the list is walked: leaving it in would make both `report()` and
  # `plot(what = "all")` call themselves.
  candidates <- setdiff(eval(formals(method)$what), "all")
  Filter(function(view) {
    if (identical(view, "responses")) return(length(.gg_response_blocks(x)) > 0L)
    if (identical(view, "sequences")) return(!is.null(x$time))
    if (identical(view, "profiles") &&
        length(.multilpa_continuous_names(x)) == 0L &&
        length(.gg_response_blocks(x)) > 0L) {
      return(.multilpa_scored_profiles(x))
    }
    # The heatmap draws the response probabilities when there are no means.
    if (identical(view, "heatmap") && length(.gg_response_blocks(x)) > 0L) {
      return(TRUE)
    }
    if (view %in% c("profiles", "bars", "heatmap")) {
      return(!is.null(x$means) && ncol(x$means) > 0L)
    }
    # The case-level views draw the stored indicators, which a centred fit
    # does not sit on, and the two multi-indicator views need two of them.
    if (view %in% c("raincloud", "parallel", "pairs")) {
      enough <- if (identical(view, "raincloud")) 1L else 2L
      return(!is.null(x$indicator_data) &&
               length(.multilpa_continuous_names(x)) >= enough &&
               !.multilpa_is_transformed(x))
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
