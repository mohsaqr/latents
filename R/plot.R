#' @rdname latents-plots
#' @export
plot.multilpa <- function(x, what = c("profiles", "bars", "heatmap", "raincloud",
                                      "parallel", "pairs", "responses",
                                      "probabilities", "sequences", "sizes",
                                      "entropy", "posteriors", "avepp", "all"),
                          data = NULL, scale = c("raw", "standardized"),
                          category = "last", labels = TRUE, intervals = TRUE,
                          cell_labels = TRUE, main = NULL, subtitle = NULL,
                          statistic = c("mean", "median", "mode"),
                          ...) {
  stopifnot("`x` must be a `multilpa` fit" = inherits(x, "multilpa"))
  .multilpa_check_plot_arguments(list(...), labels, intervals, cell_labels,
                                 category)
  what <- match.arg(what)
  scale <- match.arg(scale)
  statistic <- match.arg(statistic)
  .multilpa_check_plot_statistic(what, statistic)
  .gg_require()
  if (identical(what, "all")) {
    return(.gg_every_view(x, match.call(), parent.frame()))
  }
  .multilpa_draw_view(x, what, data = data, scale = scale,
                      category = category, labels = labels,
                      intervals = intervals, cell_labels = cell_labels,
                      main = main, subtitle = subtitle, statistic = statistic)
}

#' Validate where a profile summary statistic applies
#' @noRd
.multilpa_check_plot_statistic <- function(what, statistic) {
  stopifnot(is.character(what), is.character(statistic))
  if (!what %in% c("profiles", "all") && statistic != "mean") {
    stop(errorCondition("`statistic` applies only to `what = \"profiles\"`.",
                        class = "latents_bad_argument", call = NULL))
  }
  invisible(NULL)
}

#' Validate the arguments every fit's plot method shares
#' @noRd
.multilpa_check_plot_arguments <- function(extra, labels, intervals,
                                           cell_labels, category) {
  .multilpa_reject_extra_arguments(
    extra, "plot()",
    paste("Style the returned plot with ggplot2, for example",
          "`+ ggplot2::theme()` or `+ ggplot2::scale_colour_manual()`."))
  stopifnot(
    "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
    "`intervals` must be TRUE or FALSE" =
      isTRUE(intervals) || isFALSE(intervals),
    "`cell_labels` must be TRUE or FALSE" =
      isTRUE(cell_labels) || isFALSE(cell_labels),
    "`category` must be a single label or index" = length(category) == 1L)
  invisible(NULL)
}

#' Draw one view of a fit
#'
#' @noRd
.multilpa_draw_view <- function(x, what, data, scale, category, labels,
                                intervals, cell_labels, main, subtitle,
                                statistic = "mean") {
  if (what %in% c("entropy", "posteriors", "avepp")) {
    .multilpa_refuse_noise(x, sprintf("plot(what = \"%s\")", what),
                           class = "latents_nothing_to_plot")
  }
  with_errors <- isTRUE(intervals)
  switch(what,
    profiles = .gg_view_profiles(
      x, scale, labels,
      if (with_errors) .multilpa_mean_error_matrix(x, data), main, subtitle,
      statistic = statistic, intervals = with_errors),
    bars = .gg_view_bars(
      x, scale,
      if (with_errors) .multilpa_mean_error_matrix(x, data),
      main, subtitle),
    heatmap = .gg_view_heatmap(x, main, subtitle),
    raincloud = .gg_view_raincloud(x, scale, main, subtitle),
    parallel = .gg_view_parallel(x, main, subtitle),
    pairs = .gg_view_pairs(x, main, subtitle),
    responses = .gg_view_responses(
      x, category, labels,
      if (with_errors) .multilpa_response_error_table(x, data), main,
      subtitle),
    probabilities = .gg_view_probabilities(x, labels, main, subtitle),
    sequences = .gg_view_sequences(x, labels, cell_labels, main, subtitle),
    transitions = .gg_view_transitions(x, main, subtitle),
    sizes = .gg_view_sizes(x, main, subtitle),
    entropy = .gg_view_entropy(x, main, subtitle),
    posteriors = .gg_view_posteriors(x, main, subtitle),
    avepp = .gg_view_avepp(x, main, subtitle))
}

#' Every view a fit supports, as one list of plots
#'
#' `what = "all"` re-issues the caller's own `plot()` call once per supported
#' view, so every argument the caller gave is carried into each view. A view
#' can still refuse for a reason only it knows; a refusal moves on to the next
#' view and is named at the end instead of stopping halfway.
#'
#' @param x A fitted model of this package.
#' @param call The caller's `match.call()`, re-evaluated once per view.
#' @param env The caller's `parent.frame()`, where that call is evaluated.
#' @return A `latents_plots` list, named by view.
#' @noRd
.gg_every_view <- function(x, call, env) {
  # `match.call()` inside an S3 method names the method, and the methods are
  # registered rather than exported, so dispatch through the generic instead.
  call[[1L]] <- quote(plot)
  views <- .multilpa_supported_views(x)
  if (length(views) == 0L || anyNA(views)) {
    stop(errorCondition("This model family has no named plot views to draw.",
                        class = "latents_no_plot", call = NULL))
  }
  refused <- \(condition) NULL
  plots <- lapply(views, \(view) {
    call$what <- view
    if (view != "profiles") call$statistic <- NULL
    tryCatch(eval(call, env),
             latents_no_plot = refused,
             latents_nothing_to_plot = refused,
             latents_no_time = refused,
             latents_no_categorical = refused,
             latents_no_continuous = refused,
             latents_no_indicator_data = refused)
  })
  names(plots) <- views
  drawn <- !vapply(plots, is.null, logical(1))
  if (!all(drawn)) {
    message(sprintf("Not drawn for this model: %s.",
                    paste(views[!drawn], collapse = ", ")))
  }
  .gg_plots(plots[drawn])
}

#' Refuse a case-level diagnostic that has nothing to separate
#'
#' Both case-level views lay the cases of each profile out side by side. With
#' one profile every case is assigned to it with probability one, so the view
#' is refused rather than drawn as a single spike.
#'
#' @param x A fitted model of this package.
#' @param what The view that was asked for, for the message.
#' @return `NULL`, invisibly, when the view can be drawn.
#' @noRd
.multilpa_require_case_posteriors <- function(x, what, level = "individuals") {
  if (identical(level, "groups")) {
    if (!.gg_has_group_classes(x)) {
      stop(errorCondition(sprintf(
        "This fit has fewer than two group classes, so `what = \"%s\"` has nothing to separate.",
        what), class = "latents_no_group_classes", call = NULL))
    }
    return(invisible(NULL))
  }
  posteriors <- x$subject_posteriors
  if (is.null(posteriors) || !is.matrix(posteriors)) {
    stop(errorCondition(sprintf(
      "This fit carries no individual posteriors, so `what = \"%s\"` has nothing to draw. Views available for this fit: %s.",
      what, .multilpa_available_views_text(x)),
      class = "latents_nothing_to_plot", call = NULL))
  }
  if (ncol(posteriors) < 2L) {
    stop(errorCondition(sprintf(
      paste("A single-profile fit assigns every case to that profile with probability one,",
            "so `what = \"%s\"` has nothing to separate. Views available for this fit: %s."),
      what, .multilpa_available_views_text(x)),
      class = "latents_nothing_to_plot", call = NULL))
  }
  invisible(NULL)
}

#' The observations behind a case-level view, with each case's profile
#'
#' The indicators are taken as the fit stored them, which are the values it was
#' fitted to. A centred fit reports its means on the centred scale, where the
#' stored columns would not sit, so it is refused rather than drawn wrongly.
#' @param x A fitted model.
#' @param vars The continuous indicators.
#' @return A list with `values` (cases by indicators) and `profile`.
#' @noRd
.multilpa_plot_observations <- function(x, vars) {
  if (.multilpa_is_transformed(x)) {
    stop(errorCondition(paste(
      "This view draws the indicators as supplied, and this fit was",
      "centred or scaled, so its means are on a different scale."),
      class = c("latents_nothing_to_plot", "latents_bad_argument"),
      call = NULL))
  }
  values <- x$indicator_data
  if (is.null(values) || is.null(x$subject_profiles)) {
    stop(errorCondition("This fit did not retain its indicators, so they cannot be drawn.",
                        class = "latents_no_indicator_data", call = NULL))
  }
  values <- matrix(as.numeric(values), nrow(values), ncol(values),
                   dimnames = list(NULL, colnames(values)))
  # A noise component's cases belong to no profile, so they are not drawn.
  kept <- x$subject_profiles >= 1L & x$subject_profiles <= x$n_profiles
  list(values = values[kept, vars, drop = FALSE],
       profile = x$subject_profiles[kept])
}

#' @rdname latents-plots
#' @export
plot.multilpa_enumeration <- function(x, what = c("aic", "bic_groups",
                                               "bic_individual",
                                               "icl_individual"),
                                      criterion = NULL,
                                      combine = TRUE, labels = TRUE,
                                      mark_minimum = TRUE, main = NULL,
                                      subtitle = NULL, cutoff = 0.8,
                                      using = "bic_groups", ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"))
  # The released view names were matched with match.arg(), so their unique
  # abbreviations remain valid alongside the new criterion names.
  if (is.null(what)) what <- "enumeration"
  legacy <- if (is.character(what) && length(what) == 1L && !is.na(what)) {
    pmatch(what, c("enumeration", "tree"))
  } else NA_integer_
  if (!is.na(legacy)) what <- c("enumeration", "tree")[legacy]
  if (identical(what, "tree")) {
    # The tree takes its own arguments; any other is refused by name.
    tree_arguments <- list(...)
    .multilpa_reject_extra_arguments(
      tree_arguments[setdiff(names(tree_arguments), .latents_tree_arguments())],
      "plot(what = \"tree\")", paste(
        "The tree takes `n_group_classes`, `n_profiles`, `model`, `candidates`,",
        "`from`, `to`, `level`, `combine`, `main` and `subtitle`."))
    # In 0.9.21 `criterion` always had a value and the tree did not use it.
    .latents_warn_unused_arguments(
      if (is.null(criterion)) list() else list(criterion = criterion),
      "plot(what = \"tree\")", "The tree shows no criterion.")
    return(do.call(plot_tree, c(list(x), tree_arguments,
                                list(combine = if (missing(combine)) FALSE else combine,
                                     main = main, subtitle = subtitle))))
  }
  .multilpa_reject_extra_arguments(
    list(...), "plot()", .multilpa_enumeration_plot_advice())
  # 0.9.21 chose the view with `what = "enumeration"` and the criteria with
  # `criterion`; both spellings still draw the same plot.
  what <- .latents_former_argument(
    what, !missing(what) && !identical(what, "enumeration"), criterion,
    c("what", "criterion"))
  if ("tree" %in% what) {
    stop(errorCondition(
      "`what = \"tree\"` draws the tree on its own; ask for criteria separately.",
      class = "latents_bad_argument", call = NULL))
  }
  # The view's name before criteria were chosen through `what`.
  if (identical(what, "enumeration")) what <- eval(formals(plot_enumeration)$what)
  plot_enumeration(x, what = what, combine = combine, labels = labels,
                   mark_minimum = mark_minimum, main = main,
                   subtitle = subtitle, cutoff = cutoff, using = using)
}

#' The arguments a tree takes besides `combine`, `main` and `subtitle`
#' @return A character vector.
#' @noRd
.latents_tree_arguments <- function() {
  c("n_group_classes", "n_profiles", "model", "candidates", "from", "to", "level")
}

#' What an enumeration plot takes instead of a refused argument
#' @return One sentence for `.multilpa_reject_extra_arguments()`.
#' @noRd
.multilpa_enumeration_plot_advice <- function() {
  paste("Name the criteria with `what`, for example `what = \"bic_groups\"`",
        "or `what = \"all\"`, and style the plot with ggplot2.")
}

#' What a transition or group-class family grid's plot takes instead
#' @return One sentence for the unused-argument warning.
#' @noRd
.latents_grid_plot_advice <- function() {
  paste("Name the view with `what`, such as \"bic\", \"aic\", \"bic_individual\",",
        "\"all\" or a view listed by plot_views(x); style the plot with ggplot2.")
}

#' The view a transition or family grid's plot draws
#'
#' In 0.9.21 these methods matched `"bic"`, `"aic"` and `"bic_individual"`
#' with `match.arg()`: the full vector meant `"bic"`, and a unique
#' abbreviation such as `"bic_i"` names its criterion. Both still hold. A name
#' the grid offers in full (a criterion, `"all"` or a view) is kept, as are
#' several criteria; anything else is left for [plot_enumeration()] to refuse
#' by class.
#'
#' @param what The caller's `what`.
#' @param views `.latents_grid_views()` of the grid.
#' @return The `what` to pass to [plot_enumeration()].
#' @noRd
.latents_grid_criterion <- function(what, views) {
  if (!(is.character(what) && length(what) >= 1L && !anyNA(what))) {
    stop(errorCondition(sprintf(
      "`what` must name views of the grid, such as %s; not %s.",
      .latents_grid_offer_text(views), deparse1(what)),
      class = "latents_bad_argument", call = NULL))
  }
  choices <- c("bic", "aic", "bic_individual")
  if (identical(what, choices)) return(choices[1L])
  offered <- c(unname(views$criteria), "all", views$views)
  if (length(what) > 1L || what %in% offered) return(what)
  at <- pmatch(what, choices)
  if (is.na(at)) what else choices[at]
}

#' @rdname latents-plots
#' @export
plot_enumeration <- function(x, what = c("aic", "bic_groups",
                                         "bic_individual", "icl_individual"),
                             combine = TRUE, labels = TRUE,
                             mark_minimum = TRUE, main = NULL,
                             subtitle = NULL, cutoff = 0.8,
                             using = "bic_groups", ...) {
  grid_views <- .latents_grid_views(x, "plot_enumeration")
  # Each kind of grid has its own criteria and its own default criterion for
  # the selection views; an explicit `what` or `using` is taken as given.
  if (missing(what)) what <- grid_views$default
  if (missing(using)) using <- grid_views$using
  # The tree's own arguments reach it through `...`; no other view takes any.
  tree_arguments <- list(...)
  tree_only <- is.character(what) && identical(what, "tree")
  .multilpa_reject_extra_arguments(
    if (tree_only) tree_arguments[setdiff(names(tree_arguments),
                                          .latents_tree_arguments())] else
      tree_arguments,
    "plot_enumeration()", .multilpa_enumeration_plot_advice())
  stopifnot("`what` must name one or more criteria, or be \"all\"" =
              is.character(what) && length(what) >= 1L && !anyNA(what),
            "`combine` must be TRUE or FALSE" = isTRUE(combine) || isFALSE(combine),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`mark_minimum` must be TRUE or FALSE" =
              isTRUE(mark_minimum) || isFALSE(mark_minimum))
  .latents_require_grid_views(x, what, "plot_enumeration")
  .gg_require()
  if ("tree" %in% what) {
    if (!identical(what, "tree")) {
      stop(errorCondition(
        "`what = \"tree\"` draws the tree on its own; ask for other views separately.",
        class = "latents_bad_argument", call = NULL))
    }
    return(do.call(plot_tree, c(list(x), tree_arguments,
                                list(combine = if (missing(combine)) FALSE else combine,
                                     main = main, subtitle = subtitle))))
  }
  quality <- c("entropy", "avepp", "diagnostics")
  # Model selection: each view on its own plot; "selection" draws them all,
  # with classification quality, which is part of choosing a model.
  selection <- .latents_selection_views()
  if (any(what %in% c(selection, "selection"))) {
    if (!all(what %in% c(selection, "selection", quality))) {
      stop(errorCondition(paste(
        "Plot information criteria and the selection views separately: `what`",
        "takes criteria, or views such as \"deltas\", \"heatmap\" and \"selection\"."),
        class = "latents_bad_argument", call = NULL))
    }
    stopifnot("`cutoff` must be a single number between 0 and 1" =
                is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff) &&
                cutoff > 0 && cutoff <= 1)
    # "selection" is every selection view the grid offers, and diagnostics.
    views <- unique(unlist(lapply(what, function(view) {
      if (identical(view, "selection")) {
        intersect(c(selection, "diagnostics"), grid_views$views)
      } else view
    })))
    plots <- stats::setNames(lapply(views, function(view) {
      if (view %in% quality) {
        .gg_view_enumeration_quality(
          x, if (identical(view, "diagnostics")) c("entropy", "avepp") else view,
          labels, main, subtitle, cutoff)
      } else .gg_selection_view(x, view, using, main, subtitle)
    }), views)
    if (length(plots) == 1L) return(plots[[1L]])
    return(.gg_plots(plots))
  }
  # Classification quality: the same words as a single model's plot.
  if (any(what %in% quality)) {
    if (!all(what %in% quality)) {
      stop(errorCondition(paste(
        "Plot information criteria and classification quality separately:",
        "`what` takes criteria, or \"entropy\", \"avepp\" and \"diagnostics\"."),
        class = "latents_bad_argument", call = NULL))
    }
    stopifnot("`cutoff` must be a single number between 0 and 1" =
                is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff) &&
                cutoff > 0 && cutoff <= 1)
    measures <- if ("diagnostics" %in% what) c("entropy", "avepp") else what
    if (isFALSE(combine) && length(measures) > 1L) {
      return(.gg_plots(stats::setNames(lapply(measures, \(measure) {
        .gg_view_enumeration_quality(x, measure, labels, main, subtitle, cutoff)
      }), measures)))
    }
    return(.gg_view_enumeration_quality(x, measures, labels, main, subtitle,
                                        cutoff))
  }
  grid <- .latents_grid_table(x)
  offered <- unname(grid_views$criteria)
  criterion <- if (identical(what, "all")) offered else what
  unknown <- setdiff(criterion, offered)
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "%s is not an information criterion of this grid. Its `what` takes %s.",
      paste(sprintf("`%s`", unknown), collapse = ", "),
      .latents_grid_offer_text(grid_views)),
      class = c("latents_unknown_criterion", "latents_bad_argument"), call = NULL))
  }
  empty <- criterion[!vapply(criterion, \(name) {
    any(grid$converged %in% TRUE & is.finite(grid[[name]]))
  }, logical(1))]
  if (length(empty) > 0L) {
    stop(errorCondition(sprintf(
      "No converged candidate has a finite %s; nothing to plot.",
      paste(sprintf("`%s`", empty), collapse = ", ")),
      class = "latents_nothing_to_plot", call = NULL))
  }
  criteria <- .latents_grid_criteria_titles(x, grid, criterion)
  if (length(criteria) > 1L && isFALSE(combine)) {
    return(.gg_plots(stats::setNames(lapply(seq_along(criteria), \(index) {
      .gg_view_enumeration(x, criteria[index], labels, mark_minimum, main,
                           subtitle)
    }), criteria)))
  }
  .gg_view_enumeration(x, criteria, labels, mark_minimum, main, subtitle)
}

#' The criteria worth a panel, named for display
#'
#' A single-level grid has one unit per observation, so every criterion's
#' group-level and individual-level versions are the same numbers; drawing both
#' repeats a panel. Such a pair is kept once under the criterion's plain name.
#' A grid whose every candidate has one group class classifies nothing at the
#' group level, so the group version of a classification criterion (AWE, ICL,
#' CLC) is its penalty alone, and ICL over groups equals BIC over groups. When
#' both versions are asked for, only the individual one is kept, under the
#' criterion's plain name.
#' @param grid The enumeration table.
#' @param criterion Requested criterion columns.
#' @return The criteria to draw, named by their display titles.
#' @noRd
.multilpa_distinct_criteria <- function(grid, criterion) {
  classification <- c("awe", "icl", "clc")
  # An empty table (titles alone) has no candidates to judge.
  no_group_classes <- nrow(grid) > 0L && !is.null(grid$n_group_classes) &&
    all(grid$n_group_classes %in% 1L)
  if (no_group_classes) {
    stems <- sub("_groups$", "", criterion)
    redundant <- grepl("_groups$", criterion) & stems %in% classification &
      paste0(stems, "_individual") %in% criterion
    criterion <- criterion[!redundant]
  }
  base <- sub("_(groups|individual)$", "", criterion)
  level <- ifelse(grepl("_groups$", criterion), "groups",
                  ifelse(grepl("_individual$", criterion), "individuals", ""))
  same_at_both <- vapply(base, function(stem) {
    both <- paste0(stem, c("_groups", "_individual"))
    all(both %in% names(grid)) &&
      isTRUE(all.equal(grid[[both[1L]]], grid[[both[2L]]]))
  }, logical(1))
  plain <- same_at_both | !nzchar(level) |
    (no_group_classes & base %in% classification & level == "individuals")
  keep <- !duplicated(ifelse(same_at_both, base, criterion))
  titles <- ifelse(plain, toupper(base), sprintf("%s (%s)", toupper(base), level))
  stats::setNames(criterion[keep], titles[keep])
}

#' @rdname latents-plots
#' @export
plot.multilpa_covariates <- function(x, what = c("profiles", "bars", "heatmap",
                                                 "raincloud", "parallel",
                                                 "pairs", "responses",
                                                 "sequences", "sizes",
                                                 "entropy", "posteriors",
                                                 "avepp", "all"),
                                     data = NULL,
                                     scale = c("raw", "standardized"),
                                     category = "last", labels = TRUE,
                                     intervals = TRUE, cell_labels = TRUE,
                                     main = NULL, subtitle = NULL,
                                     statistic = c("mean", "median", "mode"), ...) {
  stopifnot("`x` must be a fitted `multilpa_covariates` model" =
              inherits(x, "multilpa_covariates"))
  .multilpa_check_plot_arguments(list(...), labels, intervals, cell_labels,
                                 category)
  # Asked for by name or by a unique abbreviation, as match.arg() would take
  # it, so "prob" gets the reason too rather than a bare match.arg() error.
  offered <- c(eval(formals(plot.multilpa_covariates)$what), "probabilities")
  asks_prevalence <- is.character(what) && length(what) == 1L && !is.na(what) &&
    identical(offered[pmatch(what, offered)], "probabilities")
  if (asks_prevalence) {
    stop(errorCondition(
      paste("A covariate model has no single profile prevalence: it varies with",
            "each unit's covariates. Use parameter_inference() for the logits."),
      class = "latents_nothing_to_plot", call = NULL))
  }
  what <- match.arg(what)
  scale <- match.arg(scale)
  statistic <- match.arg(statistic)
  .multilpa_check_plot_statistic(what, statistic)
  .gg_require()
  if (identical(what, "all")) {
    return(.gg_every_view(x, match.call(), parent.frame()))
  }
  .multilpa_draw_view(x, what, data = data, scale = scale,
                      category = category, labels = labels,
                      intervals = intervals, cell_labels = cell_labels,
                      main = main, subtitle = subtitle, statistic = statistic)
}

#' @rdname latents-plots
#' @export
plot_views <- function(x = NULL) {
  if (!is.null(x)) return(.latents_grid_view_table(x))
  data.frame(
    type = c("profiles", "bars", "heatmap", "raincloud", "parallel", "pairs",
             "responses", "probabilities", "sequences", "transitions",
             "sizes", "entropy", "posteriors", "avepp", "enumeration", "tree",
             "all"),
    group = c(rep("measurement", 7L), rep("structure", 3L),
              rep("diagnostics", 4L), rep("selection", 2L), "every"),
    description = c(
      "Profile means, medians or modes across indicators, one line per profile",
      "Profile means as grouped bars from zero, with 95% intervals",
      "Profile means in observed standard deviations from the observed mean",
      "Each indicator's distribution by assigned profile: density, box, observations",
      "Every case as a line across indicators, one panel per profile",
      "Scatter-plot matrix with each profile's 95% covariance ellipse",
      "Categorical response probabilities, one line per profile",
      "Profile prevalence within each group class, the two-level quantity",
      "Each group's profile at each occasion, one row per group",
      "Estimated transition matrix, one panel per group class (a transition fit)",
      "Effective number of cases in each profile, with its share",
      "Each case's entropy relative to a flat posterior, by profile",
      "Posterior probability of each case's assigned profile, by profile",
      "Average posterior probability: assigned profile by posterior profile",
      "Information criteria across a candidate grid (plot any enumeration; or name criteria in `what`; plot_views(grid) lists a grid's views)",
      "How classes split as more are added, one plot per tree (a profile, class or transition grid)",
      "Every view above that this fit has the ingredients for, in one call"
    ),
    stringsAsFactors = FALSE
  )
}

#' Measurement standard errors shaped like the means matrix
#'
#' @param x A fitted model.
#' @param data The data frame the model was fitted to, or `NULL` to use the
#'   columns the fit carries.
#' @return A profiles-by-indicators matrix of standard errors; or a single
#'   string saying why there are none, for the plot's subtitle.
#' @noRd
.multilpa_mean_error_matrix <- function(x, data) {
  # A fit carries the columns it was built from, so the intervals cost the
  # caller nothing to ask for and are drawn by default. Families whose standard
  # errors are not implemented refuse by condition class; the means are still
  # worth drawing without whiskers, so the refusal is caught and reported as
  # a reason rather than propagated out of a plot call. The reason travels to
  # the subtitle, so a plot without whiskers says why.
  because <- function(reason) function(condition) reason
  errors <- tryCatch(
    .multilpa_measurement_errors(x, .multilpa_resolve_data(x, data)),
    latents_unsupported_inference = because("no intervals for this model"),
    latents_unsupported_noise = because("no intervals with a noise component"),
    latents_singular_information = because("no intervals: singular information"),
    latents_incomplete_fit = because("no intervals: incomplete fit"),
    latents_bad_inference_data = because("no intervals: data do not reproduce the fit"),
    latents_boundary_fit = because("no intervals: an estimate is at its bound"),
    latents_no_converge = because("no intervals: the fit did not converge"))
  if (is.character(errors)) return(errors)
  vars <- .multilpa_continuous_names(x)
  cells <- expand.grid(indicator = vars, profile = seq_len(x$n_profiles),
                       stringsAsFactors = FALSE)
  values <- .multilpa_match_error(errors, "mean", cells$profile, cells$indicator)
  if (all(is.na(values))) return("no intervals for this model")
  matrix(values, x$n_profiles, length(vars), byrow = TRUE)
}

#' Response-probability standard errors for the responses plot
#' @param x A fitted model.
#' @param data The fitting data, or `NULL` for the columns the fit carries.
#' @return The measurement rows of the inference table; or a string saying
#'   why there are none.
#' @noRd
.multilpa_response_error_table <- function(x, data) {
  because <- function(reason) function(condition) reason
  # Probabilities on their bound are held there, so every other one keeps its
  # interval instead of one bound-active probability removing them all.
  tryCatch(
    {
      inference <- parameter_inference(x, .multilpa_resolve_data(x, data),
                                       boundary = "fix")
      inference <- inference[inference$level == "measurement", , drop = FALSE]
      attr(inference, "held") <- length(attr(inference, "fixed_at_bound") %||%
                                          character()) > 0L
      inference
    },
    latents_unsupported_inference = because("no intervals for this model"),
    latents_unsupported_noise = because("no intervals with a noise component"),
    latents_singular_information = because("no intervals: singular information"),
    latents_incomplete_fit = because("no intervals: incomplete fit"),
    latents_bad_inference_data = because("no intervals: data do not reproduce the fit"),
    latents_boundary_fit = because("no intervals: a probability is at its bound"),
    latents_no_converge = because("no intervals: the fit did not converge"))
}
