#' Plot a fitted multilevel latent profile model
#'
#' Draws one view of a fit as a ggplot object, which can be printed, saved
#' with `ggplot2::ggsave()` or styled further with `+ ggplot2::theme()`.
#' Every view shares one profile order (largest first), one profile share (the
#' posterior share every table reports) and one set of Okabe-Ito colours
#' paired with point shapes, so a profile looks the same in every view.
#' Series are labelled directly rather than through a legend.
#'
#' The plots need the ggplot2 package, which latents suggests rather than
#' requires; without it a plot is refused with an error of class
#' `latents_missing_package`.
#'
#' @param x A fitted `multilpa` model.
#' @param what Which view to draw. [plot_views()] lists every value with its
#'   group and a one-line description.
#'
#'   The measurement model: `"profiles"` draws one line per profile across the
#'   continuous indicators, or numeric categorical summaries when there are no
#'   continuous indicators (see `statistic`). `"bars"` draws continuous means as grouped bars that
#'   start at zero. `"heatmap"` draws them as a diverging grid of observed
#'   standard deviations from each indicator's observed mean, the quickest read
#'   when there are many indicators or profiles; for a fit whose indicators are
#'   all categorical it draws the response probabilities instead.
#'   `"raincloud"` shows what the means summarize: for the cases assigned to
#'   each profile, a density, the quartiles and the observations of each
#'   indicator. `"parallel"` draws every case as a line across the indicators,
#'   one panel per profile, faded by the certainty of its assignment.
#'   `"pairs"` draws a scatter-plot matrix of the indicators with each
#'   profile's 95% ellipse from the fitted covariances, the view that shows the
#'   covariance structure. `"responses"` is the categorical counterpart of
#'   `"profiles"`: one line per profile, showing the probability of a chosen
#'   category.
#'
#'   The two-level structure: `"probabilities"` plots profile prevalence within
#'   each group class. `"sequences"` draws one row per group and one column per
#'   position, filled with the assigned profile; it needs a fit made with
#'   `time =`.
#'
#'   Classification: `"sizes"` draws each profile's effective count, the
#'   posterior mass it carries. `"posteriors"` draws the posterior probability
#'   of each case's assigned profile and `"entropy"` each case's entropy
#'   relative to a flat posterior, one strip per profile with its mean marked.
#'   A one-profile fit refuses both with an error of class
#'   `latents_nothing_to_plot`. `"avepp"` draws the average posterior
#'   probability matrix: rows are assigned profiles, columns profiles, and the
#'   diagonal is the avePP usually reported.
#'
#'   `"all"` returns every view the fit has the ingredients for.
#' @param data Optional. The data frame the model was fitted to. A fit carries
#'   the columns it was built from, so intervals are drawn without it; pass it
#'   only to compute them from a different frame.
#' @param scale For the mean views and `"raincloud"`, `"raw"` keeps each
#'   indicator in its input units and `"standardized"` divides its deviation
#'   from the observed mean by its observed standard deviation, the same map
#'   as `as.data.frame(scale = "standardized")`.
#' @param category For `what = "responses"`, which category's probability to
#'   plot: `"last"`, `"first"`, or a single category label or index.
#' @param statistic For `what = "profiles"`, `"mean"` (default), `"median"`
#'   or `"mode"`. With continuous indicators these are the same fitted Gaussian
#'   location. When there are no continuous indicators, summarizes each
#'   categorical indicator's fitted probabilities on its numeric category scale.
#'   Means assume meaningful score spacing; medians are the smallest score with
#'   cumulative probability at least 0.5; tied modes use the smallest score.
#'   Category labels must be distinct finite numeric scores. Categorical
#'   summaries support only `scale = "raw"` and have no confidence intervals;
#'   they do not run parameter inference. Mixed models retain continuous-only
#'   profile plots; use `"responses"` for their categorical indicators.
#' @param labels `TRUE` labels each series at its right end; `FALSE` uses a
#'   legend instead.
#' @param intervals For `"profiles"` and `"responses"`, `TRUE` draws 95%
#'   intervals when the fit has standard errors (a probability's is clipped to
#'   `[0, 1]`). A fit without them is drawn without whiskers and its subtitle
#'   says why.
#' @param cell_labels For `what = "sequences"`, `TRUE` prints the profile
#'   number in each cell while the grid is sparse enough to hold one, so the
#'   profile is not carried by colour alone.
#' @param main,subtitle Title and subtitle. `NULL` uses the view's own.
#' @param ... Nothing further is accepted; an unknown argument raises an error
#'   of class `latents_bad_argument`. Style the returned plot with ggplot2.
#' @return A ggplot object; for `what = "all"`, a `latents_plots` list of them,
#'   named by view, that draws every one when printed.
#' @details `"profiles"` and `"bars"` draw 95% Wald intervals when the fit has
#'   standard errors. Profile numbers are arbitrary, so two fits must have
#'   their labels aligned before their plots are compared. With
#'   `scale = "standardized"` the standard deviations are the observed
#'   indicator standard deviations, not the within-profile ones, so the values
#'   are comparable across indicators but are not effect sizes.
#' @section Errors: `latents_missing_package` without ggplot2;
#'   `latents_no_continuous`, `latents_no_categorical`, `latents_no_time` and
#'   `latents_nothing_to_plot` when the fit lacks what a view needs;
#'   `latents_bad_argument` for an argument the method does not use.
#' @references Zappia, L. and Oshlack, A. (2018). Clustering trees: a
#'   visualization for evaluating clusterings at multiple resolutions.
#'   *GigaScience*, 7(7), giy083.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(7)
#'   example_data <- data.frame(
#'     school = rep(seq_len(12), each = 10),
#'     score_a = rnorm(120), score_b = rnorm(120)
#'   )
#'   fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2)
#'   plot(fit)
#'   plot(fit, scale = "standardized")
#'   plot(fit, what = "bars")
#'   plot(fit, what = "pairs")
#'   plot(fit, what = "avepp")
#' }
#' @seealso [plot_views()] for the catalogue of views.
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
#' @param errors_available `FALSE` for a family whose standard errors are not
#'   wired into its plots, which then draws point estimates.
#' @noRd
.multilpa_draw_view <- function(x, what, data, scale, category, labels,
                                intervals, cell_labels, main, subtitle,
                                errors_available = TRUE, statistic = "mean") {
  if (what %in% c("entropy", "posteriors", "avepp")) {
    .multilpa_refuse_noise(x, sprintf("plot(what = \"%s\")", what),
                           class = "latents_nothing_to_plot")
  }
  with_errors <- isTRUE(intervals) && isTRUE(errors_available)
  switch(what,
    profiles = .gg_view_profiles(
      x, scale, labels,
      if (with_errors) .multilpa_mean_error_matrix(x, data), main, subtitle,
      statistic = statistic, intervals = with_errors),
    bars = .gg_view_bars(
      x, scale,
      if (isTRUE(errors_available)) .multilpa_mean_error_matrix(x, data),
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
.multilpa_require_case_posteriors <- function(x, what) {
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
  if (!identical(x$centering %||% "none", "none")) {
    stop(errorCondition(paste(
      "This view draws the indicators as supplied, and this fit was",
      "centred, so its means are on a different scale."),
      class = "latents_bad_argument", call = NULL))
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

#' Plot a class-enumeration grid
#'
#' `what = "enumeration"` draws information criteria against the number of
#' profiles, one panel per criterion and one line per covariance model and
#' number of group classes. A candidate that failed to converge is marked with
#' a cross on its panel's floor, so a gap in a line reads as a failure and not
#' as a missing candidate. `what = "tree"` draws how profiles split as more are
#' added: each row is one solution divided into its profiles by posterior
#' share, and bands carry the posterior mass shared by a profile and one in the
#' next solution, coloured by lineage. A band that forks is a profile split in
#' two; bands that merge or cross are cases reshuffled.
#'
#' @param x A `multilpa_enumeration` result from [enumerate_classes()],
#'   [enumerate_lpa()] or [enumerate_lca()].
#' @param what `"enumeration"` (the default) or `"tree"`.
#' @param criterion For `"enumeration"`, one or more criterion columns of
#'   `as.data.frame(x)`, such as `"bic_individual"` or `"sabic_groups"`, or
#'   `"all"` for every information criterion. The default draws AIC, BIC
#'   under both sample-size conventions and ICL counted over individuals; a
#'   criterion identical at both levels, as in a single-level grid, is drawn
#'   once. [plot_enumeration()] draws the same view with these arguments
#'   listed on its own.
#' @param combine `TRUE` draws several criteria as panels of one plot; `FALSE`
#'   returns one plot per criterion.
#' @param labels `TRUE` labels each series at its right end.
#' @param mark_minimum `TRUE` rings each criterion's lowest value among
#'   converged candidates. This marks an extremum; it does not select a model.
#' @param main,subtitle Title and subtitle. `NULL` uses the view's own.
#' @param ... Nothing further is accepted; an unknown argument raises an error
#'   of class `latents_bad_argument`.
#' @return A ggplot object, or with `combine = FALSE` and several criteria a
#'   `latents_plots` list of them, one per criterion.
#' @section Errors: `latents_unknown_criterion` for a criterion that is not a
#'   column of the grid; `latents_nothing_to_plot` when no converged candidate
#'   has a finite value, or when a tree has fewer than two numbers of profiles
#'   for every model.
#' @references Zappia, L. and Oshlack, A. (2018). Clustering trees: a
#'   visualization for evaluating clusterings at multiple resolutions.
#'   *GigaScience*, 7(7), giy083.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(7)
#'   example_data <- data.frame(score_a = rnorm(120), score_b = rnorm(120))
#'   candidates <- enumerate_lpa(example_data, c("score_a", "score_b"),
#'                               n_profiles = 1:3, n_starts = 2)
#'   plot(candidates)
#'   plot(candidates, what = "tree")
#' }
#' @export
plot.multilpa_enumeration <- function(x, what = c("enumeration", "tree"),
                                      criterion = c("aic", "bic_groups",
                                                    "bic_individual",
                                                    "icl_individual"),
                                      combine = TRUE, labels = TRUE,
                                      mark_minimum = TRUE, main = NULL,
                                      subtitle = NULL, ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"))
  .multilpa_reject_extra_arguments(
    list(...), "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  what <- match.arg(what)
  if (identical(what, "tree")) {
    .gg_require()
    return(.gg_view_tree(x, main, subtitle))
  }
  plot_enumeration(x, criterion = criterion, combine = combine,
                   labels = labels, mark_minimum = mark_minimum,
                   main = main, subtitle = subtitle)
}

#' Plot information criteria across an enumeration grid
#'
#' Draws information criteria against the number of profiles, one panel per
#' criterion and one line per covariance model and number of group classes. It
#' is the view `plot(x)` draws for an enumeration, as a function of its own so
#' that its arguments are listed (and completed by an editor) without going
#' through `what =`. A candidate that failed to converge is marked with a cross
#' on its panel's floor, so a gap in a line reads as a failure and not as a
#' missing candidate.
#'
#' @param x A `multilpa_enumeration` result from [enumerate_classes()],
#'   [enumerate_lpa()] or [enumerate_lca()].
#' @param criterion One or more criterion columns of `as.data.frame(x)`, such
#'   as `"bic_individual"` or `"sabic_groups"`, or `"all"` for every
#'   information criterion in the grid. The default draws AIC, BIC under both
#'   sample-size conventions and ICL counted over individuals. A criterion
#'   identical at both levels, as in a single-level grid, is drawn once.
#' @param combine `TRUE` draws several criteria as panels of one plot; `FALSE`
#'   returns one plot per criterion.
#' @param labels `TRUE` labels each series at its right end.
#' @param mark_minimum `TRUE` rings each criterion's lowest value among
#'   converged candidates. This marks an extremum; it does not select a model.
#' @param main,subtitle Title and subtitle. `NULL` uses the view's own.
#' @return A ggplot object, one panel per criterion; with `combine = FALSE`
#'   and several criteria, a `latents_plots` list of ggplot objects, one per
#'   criterion, named by its criterion column.
#' @section Errors: `latents_unknown_criterion` for a criterion that is not an
#'   information criterion of the grid; `latents_nothing_to_plot` when no
#'   converged candidate has a finite value.
#' @seealso [plot.multilpa_enumeration()] for the profile tree,
#'   [summary.multilpa_enumeration()] for the candidate each criterion prefers.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(7)
#'   example_data <- data.frame(score_a = rnorm(120), score_b = rnorm(120))
#'   candidates <- enumerate_lpa(example_data, c("score_a", "score_b"),
#'                               n_profiles = 1:3, n_starts = 2)
#'   plot_enumeration(candidates)
#'   plot_enumeration(candidates, criterion = "all")
#' }
#' @export
plot_enumeration <- function(x, criterion = c("aic", "bic_groups",
                                              "bic_individual",
                                              "icl_individual"),
                             combine = TRUE, labels = TRUE,
                             mark_minimum = TRUE, main = NULL,
                             subtitle = NULL) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"),
            "`criterion` must name one or more columns, or be \"all\"" =
              is.character(criterion) && length(criterion) >= 1L &&
              !anyNA(criterion),
            "`combine` must be TRUE or FALSE" = isTRUE(combine) || isFALSE(combine),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`mark_minimum` must be TRUE or FALSE" =
              isTRUE(mark_minimum) || isFALSE(mark_minimum))
  .gg_require()
  grid <- as.data.frame(x)
  if (identical(criterion, "all")) criterion <- .multilpa_enumeration_criteria()
  unknown <- setdiff(criterion, .multilpa_enumeration_criteria())
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "%s is not an information criterion in the enumeration grid.",
      paste(sprintf("`%s`", unknown), collapse = ", ")),
      class = "latents_unknown_criterion", call = NULL))
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
  criteria <- .multilpa_distinct_criteria(grid, criterion)
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
#' @param grid The enumeration table.
#' @param criterion Requested criterion columns.
#' @return The criteria to draw, named by their display titles.
#' @noRd
.multilpa_distinct_criteria <- function(grid, criterion) {
  base <- sub("_(groups|individual)$", "", criterion)
  level <- ifelse(grepl("_groups$", criterion), "groups",
                  ifelse(grepl("_individual$", criterion), "individuals", ""))
  same_at_both <- vapply(base, function(stem) {
    both <- paste0(stem, c("_groups", "_individual"))
    all(both %in% names(grid)) &&
      isTRUE(all.equal(grid[[both[1L]]], grid[[both[2L]]]))
  }, logical(1))
  keep <- !duplicated(ifelse(same_at_both, base, criterion))
  titles <- ifelse(same_at_both | !nzchar(level), toupper(base),
                   sprintf("%s (%s)", toupper(base), level))
  stats::setNames(criterion[keep], titles[keep])
}

#' Plot a covariate model
#'
#' Draws the parts of a covariate fit that are still fixed quantities: the
#' measurement model, the assignments and the classification diagnostics, as
#' [plot.multilpa()] draws them. Profile prevalence cannot be drawn: a
#' covariate model has no single prevalence vector, because prevalence varies
#' with each unit's covariates, so asking for it is refused rather than
#' answered with an average that no unit has.
#'
#' @param x A fitted `multilpa_covariates` model.
#' @param what The view; every view of [plot.multilpa()] except
#'   `"probabilities"`.
#' @param data,scale,category,labels,intervals,cell_labels,main,subtitle,statistic,...
#'   As in [plot.multilpa()].
#' @return A ggplot object; for `what = "all"`, a `latents_plots` list.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(5)
#'   school <- rep(seq_len(16), each = 8)
#'   high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#'   x <- rnorm(128)
#'   profile <- ifelse(
#'     runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#'   )
#'   example_data <- data.frame(
#'     school = school, x = x,
#'     y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'     y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#'   )
#'   fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                   n_group_classes = 2, profile_covariates = "x",
#'                   n_starts = 2)
#'   plot(fit)
#' }
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
  if (identical(what, "probabilities")) {
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

#' The plots this package can draw
#'
#' Returns the catalogue of `what =` values accepted by the `plot()` methods,
#' with the group each belongs to, so the available views can be listed rather
#' than recalled from a help page.
#'
#' @return A base `data.frame` with one row per plot type and the columns
#'   `type` (the value to pass as `what`), `group` (`"measurement"`,
#'   `"structure"`, `"diagnostics"`, `"selection"` or `"every"`), and
#'   `description`.
#' @examples
#' plot_views()
#' @export
plot_views <- function() {
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
      "Information criteria across a candidate grid (plot an enumeration)",
      "How profiles split as more are added (plot an enumeration)",
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
