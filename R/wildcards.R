# Wild-card verbs: one function per picture or table, each taking any object
# the package returns. The function's name is the choice; there is no `what`.
# They read the object and call the existing method that draws or computes
# that view, so a view looks and computes exactly as it does through plot(),
# plot_diagnostics(), plot_enumeration(), diagnostics() and get_results().

#' The kind of object a wild-card verb was given
#'
#' @param x Any object.
#' @return One string: a fit kind of `.latents_plot_kind()`, `"grid"`,
#'   `"comparison"`, `"bootstrap_lrt"`, `"bootstrap"`, `"sensitivity"`, or
#'   `NA` for an object of no kind the package returns.
#' @noRd
.latents_wild_kind <- function(x) {
  if (inherits(x, .latents_enumeration_classes)) return("grid")
  results <- c(latents_comparison = "comparison",
               multilpa_bootstrap_lrt = "bootstrap_lrt",
               latents_bootstrap = "bootstrap",
               latents_sensitivity = "sensitivity")
  matched <- results[intersect(class(x), names(results))]
  if (length(matched) > 0L) return(unname(matched[1L]))
  .latents_plot_kind(x)
}

#' Which existing view each wild-card plot verb draws, per kind of object
#'
#' Fits of profiles, classes and transitions (`profile_model`) draw through
#' `plot(x, what = )`, or through `plot_diagnostics()` for the four
#' classification views, which then add the group level of a fit with group
#' classes. Regression and growth fits draw through `plot(x, what = )`, grids
#' through `plot_enumeration(x, what = )`, and comparisons through `plot(x)`.
#' @return A named list, one element per verb: a named character vector,
#'   kind -> view.
#' @noRd
.latents_wild_views <- function() {
  regression <- c("growth_mixture", "trajectory_regression", "mixture_regression")
  # Every fitted family classifies units; those whose own plot() has no
  # classification view draw it from their posteriors (.latents_wild_classes()).
  classified <- c("profile_model", regression, "group_class_family",
                  "cross_level", "general_transitions")
  each <- function(kinds, view) stats::setNames(rep(view, length(kinds)), kinds)
  list(
    plot_bars = c(profile_model = "bars"),
    plot_heatmap = c(profile_model = "heatmap", grid = "heatmap"),
    plot_raincloud = c(profile_model = "raincloud"),
    plot_parallel = c(profile_model = "parallel"),
    plot_pairs = c(profile_model = "pairs"),
    plot_sizes = c(each(classified, "sizes"), grid = "smallest"),
    plot_entropy = c(each(classified, "entropy"), grid = "entropy"),
    plot_posteriors = each(classified, "posteriors"),
    plot_avepp = c(each(classified, "avepp"), grid = "avepp"),
    plot_trajectories = c(growth_mixture = "trajectories",
                          trajectory_regression = "trajectories",
                          mixture_regression = "fitted"),
    plot_individuals = each(regression[1:2], "individuals"),
    plot_classification = each(regression[1:2], "classification"),
    plot_coefficients = each(regression, "coefficients"),
    plot_criteria = c(grid = "all", comparison = NA_character_),
    plot_deltas = c(grid = "deltas"),
    plot_trust = c(grid = "trust"))
}

#' Draw one wild-card view of any object
#'
#' @param x Any object; a diagnostics result is drawn from its fit.
#' @param verb The calling verb, the key of `.latents_wild_views()`.
#' @param expression The caller's `substitute(x)`, for messages.
#' @param main,subtitle Titles, passed to the drawing method.
#' @param ... Further arguments of a view the fit's own plot() draws, such as
#'   `persons` for the individuals view; passed only to that plot() call.
#' @return A ggplot object, or a `latents_plots` list of the individual and
#'   group levels for a classification view of a fit with group classes.
#' @noRd
.latents_wild_plot <- function(x, verb, expression, main, subtitle, ...) {
  if (inherits(x, "multilpa_diagnostics")) x <- x$fit
  kind <- .latents_wild_kind(x)
  if (is.na(kind)) .latents_wild_refuse_object(x, verb, expression)
  views <- .latents_wild_views()[[verb]]
  if (!kind %in% names(views)) .latents_wild_refuse_view(x, kind, verb)
  view <- views[[kind]]
  # Every refusal is this verb's own, naming the verbs that do draw this
  # object; an existing method's refusal would advise a `what =` instead.
  if (!.latents_wild_applies(x, kind, view)) .latents_wild_refuse_view(x, kind, verb)
  .gg_require()
  if (identical(kind, "grid")) {
    return(plot_enumeration(x, what = view, main = main, subtitle = subtitle))
  }
  if (is.na(view)) return(plot(x, main = main, subtitle = subtitle))
  classification <- c("sizes", "posteriors", "entropy", "avepp")
  if (identical(kind, "profile_model") && view %in% classification) {
    individuals <- plot_diagnostics(x, what = view, main = main, subtitle = subtitle)
    # The group level is a second classification only when there are two or
    # more group classes; a noise fit draws the sizes alone, as
    # plot_diagnostics(what = "all") does.
    groups_too <- .gg_has_group_classes(x) &&
      (identical(view, "sizes") || !isTRUE(x$noise))
    if (!groups_too) return(individuals)
    groups <- plot_diagnostics(x, what = paste0("group_", view), main = main,
                               subtitle = subtitle)
    return(.gg_plots(list(individuals = individuals, groups = groups)))
  }
  # The other families' own "posteriors" views are not this view: regressions
  # draw every row's raw posterior (zero-weight rows included, groups left
  # out), and trajectory and growth fits alias it to the stacked membership
  # chart of plot_classification(). Every classification view of these
  # families is drawn from the eligible, weighted posteriors instead.
  if (view %in% classification) {
    return(.latents_wild_classes(x, view, verb, main, subtitle))
  }
  # A view the fit's own plot() does not offer (the case-level views of a
  # transition fit) is refused here, by class, rather than by match.arg().
  if (!view %in% .latents_plot_offered_views(x)) .latents_wild_refuse_view(x, kind, verb)
  plot(x, what = view, main = main, subtitle = subtitle, ...)
}

#' A classification view drawn from any fit's posteriors
#'
#' For the families whose own plot() has no such view. Each classifying level
#' with two or more classes is drawn (sizes also with one); two levels give a
#' `latents_plots` list with `individuals` and `groups`.
#' @param x A fitted model of a family other than profiles, classes and
#'   homogeneous transitions.
#' @param view `"sizes"`, `"posteriors"`, `"entropy"` or `"avepp"`.
#' @param verb The calling verb, for messages.
#' @param main,subtitle Titles.
#' @return A ggplot object or a `latents_plots` list.
#' @noRd
.latents_wild_classes <- function(x, view, verb, main, subtitle) {
  levels <- .latents_wild_class_levels(x, view)
  units <- lapply(stats::setNames(names(levels), names(levels)),
                  \(level) .latents_wild_unit(x, level))
  units <- Filter(Negate(is.null), units)
  if (length(units) == 0L) {
    stop(errorCondition(sprintf(
      "%s(): %s has no classified units to draw.", verb, .latents_describe(x)),
      class = "latents_nothing_to_plot", call = NULL))
  }
  draw <- switch(view,
    sizes = .latents_wild_sizes,
    avepp = .latents_wild_avepp,
    posteriors = \(unit, main, subtitle) .latents_wild_strip(unit, "certainty", main, subtitle),
    entropy = \(unit, main, subtitle) .latents_wild_strip(unit, "entropy", main, subtitle))
  plots <- lapply(units, draw, main = main, subtitle = subtitle)
  if (length(plots) == 1L) plots[[1L]] else .gg_plots(plots)
}

#' One classifying level of any fit, in the form the drawers read
#' @param x A fitted model.
#' @param level `"individuals"` or `"groups"`.
#' @return A list: `posteriors` (eligible rows), `weights` (`NULL` when
#'   unweighted), `key` (one row per class: `class`, `rank`, `share`,
#'   `count`, `label`, `colour`, `shape`), `class_word`, `unit_word`, `total`;
#'   `NULL` when the level has no posteriors.
#' @noRd
.latents_wild_unit <- function(x, level) {
  pieces <- .latents_classification_inputs(x, level)
  posteriors <- pieces$posteriors
  if (is.null(posteriors) || identical(pieces$weights, NA) || nrow(posteriors) == 0L) {
    return(NULL)
  }
  weights <- pieces$weights %||% rep(1, nrow(posteriors))
  counts <- colSums(posteriors * weights)
  n_classes <- length(counts)
  share <- counts / sum(counts)
  ordering <- order(-share, seq_len(n_classes))
  rank <- integer(n_classes)
  rank[ordering] <- seq_len(n_classes)
  groups <- identical(level, "groups")
  class_word <- if (groups) "group class" else "class"
  label <- sprintf("%s %d (%s)", .gg_capitalize(class_word), seq_len(n_classes),
                   .gg_percent(share))
  key <- data.frame(class = seq_len(n_classes), rank = rank, share = share,
                    count = unname(counts),
                    label = factor(label, levels = label[ordering]),
                    colour = rep_len(.gg_okabe_ito, n_classes)[rank],
                    shape = rep_len(.gg_shapes, n_classes)[rank],
                    stringsAsFactors = FALSE)
  list(posteriors = posteriors, weights = pieces$weights, key = key,
       class_word = class_word, unit_word = if (groups) "group" else "unit",
       level = level, total = sum(counts))
}

#' Effective class sizes of one level, largest first
#' @noRd
.latents_wild_sizes <- function(unit, main, subtitle) {
  key <- unit$key
  key$class_label <- factor(as.character(key$label), levels = rev(levels(key$label)))
  key$bar_label <- sprintf("%.1f  (%.1f%%)", key$count, 100 * key$share)
  ggplot2::ggplot(key, ggplot2::aes(count, class_label, fill = label)) +
    ggplot2::geom_col(width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = bar_label), hjust = -0.08,
                       size = 3.4, colour = "grey20") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.24))) +
    .gg_scales(key, "fill") +
    .gg_titles(main, subtitle, sprintf("%s sizes", .gg_capitalize(unit$class_word)),
               sprintf(paste("Effective number of %ss (sum of posterior",
                             "probabilities); %s in total"), unit$unit_word,
                       format(round(unit$total), big.mark = ","))) +
    ggplot2::labs(x = sprintf("Effective %ss", unit$unit_word), y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Each unit's certainty or relative entropy by assigned class, with means
#' @noRd
.latents_wild_strip <- function(unit, value, main, subtitle) {
  posteriors <- unit$posteriors
  key <- unit$key
  k <- ncol(posteriors)
  assigned <- max.col(posteriors, ties.method = "first")
  # -sum p log p relative to log(K): 0 is certain, 1 is a flat posterior.
  entropy <- -rowSums(posteriors * log(pmax(posteriors, .Machine$double.xmin)))
  values <- if (identical(value, "certainty")) {
    pmin(posteriors[cbind(seq_len(nrow(posteriors)), assigned)], 1)
  } else pmax(entropy / log(k), 0)
  weights <- unit$weights %||% rep(1, length(values))
  y0 <- k + 1L - key$rank
  cases <- data.frame(value = values, class = assigned, y0 = y0[assigned],
                      class_label = key$label[assigned])
  # A deterministic spread within each class's band, so points do not overlap
  # in one line and the plot is identical on every draw.
  within <- stats::ave(cases$value, cases$class, FUN = \(v) rank(v, ties.method = "first"))
  cases$y <- cases$y0 + 0.5 * ((within * 0.618) %% 1 - 0.5)
  means <- data.frame(class = seq_len(k), y0 = y0,
                      value = vapply(seq_len(k), \(j) {
                        chosen <- assigned == j
                        if (!any(chosen)) NA_real_ else
                          sum(values[chosen] * weights[chosen]) / sum(weights[chosen])
                      }, numeric(1)))
  means <- means[!is.na(means$value), , drop = FALSE]
  means$label <- sprintf("mean %.2f", means$value)
  certainty <- identical(value, "certainty")
  reference <- if (certainty) 1 / k else 1
  ggplot2::ggplot(cases, ggplot2::aes(value, y, colour = class_label)) +
    ggplot2::geom_vline(xintercept = reference, linetype = "22",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_point(size = 1.1, alpha = 0.55) +
    ggplot2::geom_segment(data = means,
                          ggplot2::aes(x = value, xend = value,
                                       y = y0 - 0.32, yend = y0 + 0.32),
                          colour = "grey10", linewidth = 0.8, inherit.aes = FALSE) +
    ggplot2::geom_text(data = means,
                       ggplot2::aes(x = value, y = y0 + 0.34, label = label),
                       colour = "grey10", size = 3, vjust = 0,
                       hjust = if (certainty) 1.05 else -0.05, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(limits = if (certainty) c(1 / k, 1) else c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = y0, labels = as.character(key$label),
                                limits = c(0.5, k + 0.5),
                                expand = ggplot2::expansion(0)) +
    .gg_scales(key, "colour") +
    ggplot2::coord_cartesian(clip = "off") +
    .gg_titles(main, subtitle,
               if (certainty) "Classification certainty" else "Classification uncertainty",
               if (certainty) {
                 sprintf("Posterior probability of each %s's assigned %s; dashed line: chance",
                         unit$unit_word, unit$class_word)
               } else {
                 sprintf(paste("Each %s's entropy relative to a flat posterior",
                               "(0 = certain, 1 = flat); relative entropy = %.3f"),
                         unit$unit_word,
                         .multilpa_relative_entropy(posteriors, unit$weights))
               }) +
    ggplot2::labs(x = if (certainty) sprintf("Posterior probability of the assigned %s",
                                             unit$class_word) else
                    sprintf("Relative %s entropy", unit$unit_word), y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Average posterior matrix of one level: assigned class by posterior class
#' @noRd
.latents_wild_avepp <- function(unit, main, subtitle) {
  posteriors <- unit$posteriors
  key <- unit$key
  averages <- .multilpa_average_posterior_matrix(posteriors, unit$weights)
  k <- ncol(averages)
  assigned <- tabulate(max.col(posteriors, ties.method = "first"), k)
  display <- order(key$rank)
  title_word <- .gg_capitalize(unit$class_word)
  cells <- expand.grid(column = seq_len(k), row = seq_len(k))
  cells$average_posterior <- averages[cbind(cells$row, cells$column)]
  assigned_levels <- sprintf("Assigned %d (n = %d)", seq_len(k), assigned)
  cells$assigned_label <- factor(assigned_levels[cells$row],
                                 levels = rev(assigned_levels[display]))
  cells$posterior_label <- factor(sprintf("%s %d", title_word, cells$column),
                                  levels = sprintf("%s %d", title_word, display))
  cells$text_colour <- .gg_ink(.gg_probability_fill(cells$average_posterior))
  cells$label <- ifelse(is.na(cells$average_posterior), "",
                        sprintf("%.2f", cells$average_posterior))
  diagonal <- cells[cells$row == cells$column, , drop = FALSE]
  ggplot2::ggplot(cells, ggplot2::aes(posterior_label, assigned_label,
                                      fill = average_posterior)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_tile(data = diagonal, fill = NA, colour = "grey15", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = text_colour), size = 3.4) +
    ggplot2::scale_colour_identity() +
    .gg_probability_scale("Average\nposterior") +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    .gg_titles(main, subtitle,
               sprintf("Average posterior probability of the %ss", unit$class_word),
               sprintf(paste("Mean posterior the %ss assigned to each %s place on",
                             "each %s; the outlined diagonal is classification",
                             "accuracy"), unit$unit_word, unit$class_word,
                       unit$class_word)) +
    ggplot2::labs(x = NULL, y = NULL) +
    .gg_theme() +
    .gg_tile_theme()
}

#' The wild-card plot verbs that apply to this object
#' @param x The object (a diagnostics result already replaced by its fit).
#' @param kind Its kind of `.latents_wild_kind()`.
#' @return A character vector of verb names, ending with `plot_all`.
#' @noRd
.latents_wild_verbs_for <- function(x, kind) {
  views <- .latents_wild_views()
  mapped <- names(views)[vapply(views, \(v) kind %in% names(v), logical(1))]
  applies <- vapply(mapped, \(verb) .latents_wild_applies(x, kind, views[[verb]][[kind]]),
                    logical(1))
  c(.latents_wild_named_verbs(x, kind), mapped[applies], "plot_all")
}

#' Can this object draw this view?
#' @param x The object.
#' @param kind Its kind.
#' @param view The view a wild-card verb maps to for this kind.
#' @return A single logical.
#' @noRd
.latents_wild_applies <- function(x, kind, view) {
  if (is.na(view)) return(TRUE)
  if (identical(kind, "grid")) {
    views <- plot_views(x)
    return(view %in% views$type[views$supported %in% TRUE])
  }
  classification <- view %in% c("sizes", "posteriors", "entropy", "avepp")
  if (identical(kind, "profile_model")) {
    # A noise component leaves the profile posteriors short of one, so only
    # the sizes are drawn for it (as plot_diagnostics(what = "all") does).
    if (classification && isTRUE(x$noise)) return(identical(view, "sizes"))
    return(view %in% .multilpa_supported_views(x))
  }
  if (classification) return(length(.latents_wild_class_levels(x, view)) > 0L)
  view %in% .latents_plot_offered_views(x)
}

#' The classifying levels a classification view of a non-profile fit draws
#'
#' One class has sizes and a one-by-one AvePP matrix; its certainty and
#' entropy are constant by construction and are not drawn.
#' @param x A fitted model.
#' @param view The classification view.
#' @return A named integer vector (level -> number of classes).
#' @noRd
.latents_wild_class_levels <- function(x, view) {
  if (!inherits(x, .latents_model_classes)) return(integer())
  levels <- .latents_classifying_levels(x)
  levels[!is.na(levels) & (levels > 1L | view %in% c("sizes", "avepp"))]
}

#' The existing named plot functions that draw this fit
#'
#' [plot_profiles()], [plot_responses()], [plot_sequences()],
#' [plot_prevalence()] and [plot_transitions()] take any fit too; a refusal
#' names the ones whose view this fit has.
#' @param x The object.
#' @param kind Its kind.
#' @return A character vector of function names.
#' @noRd
.latents_wild_named_verbs <- function(x, kind) {
  if (!inherits(x, .latents_model_classes) && !identical(kind, "pooled")) {
    return(character())
  }
  registry <- .latents_named_plot_registry()
  rows <- registry[registry$kind %in% .latents_plot_kind(x) &
                     registry$what %in% c("", "lines") &
                     registry$verb != "plot_diagnostics", , drop = FALSE]
  available <- if (identical(kind, "profile_model")) .multilpa_supported_views(x) else
    .latents_plot_offered_views(x)
  unique(rows$verb[is.na(rows$view) | rows$view %in% available])
}

#' Refuse a view the object has no counterpart of
#' @noRd
.latents_wild_refuse_view <- function(x, kind, verb) {
  stop(errorCondition(sprintf(
    "%s() has nothing to draw for %s%s. For this object: %s.", verb,
    .latents_describe(x),
    if (isTRUE(x$noise)) " with a noise component" else "",
    paste(sprintf("%s()", .latents_wild_verbs_for(x, kind)), collapse = ", ")),
    class = c("latents_nothing_to_plot", "latents_bad_argument"), call = NULL))
}

#' Refuse an object the package did not return
#' @noRd
.latents_wild_refuse_object <- function(x, verb, expression) {
  label <- .latents_argument_label(expression, "x")
  stop(errorCondition(sprintf(paste(
    "%s() takes a fitted model, an enumeration grid, a comparison, a bootstrap",
    "test, a sensitivity analysis or diagnostics, but `%s` is %s."),
    verb, label, .latents_describe(x)),
    class = c("latents_not_a_fit", "latents_bad_argument"), call = NULL))
}

#' Plot how many units each class holds, or the size of the smallest class of
#' every candidate in a grid
#'
#' `plot_sizes()` plots how many units belong to each class of a fitted model,
#' as bars, or how large the smallest class of every candidate of an enumeration
#' grid is, as a line.
#'
#' @details
#' **A fitted model.** One horizontal bar per class, the largest at the top.
#' The length of a bar is the effective number of units in the class: the sum
#' of the units' posterior probabilities of that class. Each bar is labelled
#' with that number and its share of the total, and each class has its own
#' colour.
#'
#' A fit classifies at one or two levels. The individual level sorts the
#' cases (rows) into profiles or classes. The group level sorts the groups
#' (the persons, for example) into group classes. A fit with two or more
#' group classes besides its profiles is drawn at both levels, as two plots.
#' An additive fit and a mixture regression with `class_level = "group"`
#' classify the groups only and give one plot. A fit with a noise component
#' draws the sizes of its profiles; the noise cases are left out of the bars
#' and a caption says so.
#'
#' **An enumeration grid.** The x-axis is the number of classes, the y-axis
#' the estimated share of the smallest class of each candidate. There is one
#' line per number of group classes, and per covariance structure when the
#' grid has several. A dashed line marks 5%. A grid with group classes has a
#' second panel for the smallest group class.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. This covers fits with
#'     covariates, growth mixtures, and the additive and cross-level families
#'     of [multilpa()];
#'   * an enumeration grid from [enumerate_lpa()], [enumerate_lca()],
#'     [enumerate_classes()] or [enumerate_regressions()]. A grid of
#'     group-class families (`enumerate_classes(family = )`) has no such
#'     view and is refused;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text. When two plots are returned, both get the
#'   same title and subtitle.
#' @return A ggplot object. A fit that classifies at both levels returns a
#'   `latents_plots` list of two ggplot objects named `individuals` and
#'   `groups`; printing the list draws both.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a grid of group-class families, a grid in
#'     which no converged candidate has two or more classes, pooled
#'     imputations, a comparison, a bootstrap result or a sensitivity analysis.
#'     For an object of the wrong kind, the message names the plot functions
#'     that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [get_classification()] for the sizes as a table;
#'   [plot_posteriors()], [plot_entropy()], [plot_avepp()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_sizes(fit)
#'
#' # Two levels: profiles of the rows and group classes of the students
#' two_level <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
#'                       "student", n_profiles = 2, n_group_classes = 2,
#'                       n_starts = 2, seed = 1)
#' plot_sizes(two_level)
#'
#' # A grid: the smallest class of every candidate
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_sizes(grid)
#' @export
plot_sizes <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_sizes", substitute(x), main, subtitle)
}

#' Plot how uncertain the classification is, by unit for a fitted model or by
#' candidate for a grid
#'
#' `plot_entropy()` plots how uncertain the classification is: for a fitted
#' model, the entropy of each unit, by the class the unit is assigned to; for an
#' enumeration grid, the relative entropy of every candidate.
#'
#' @details
#' **A fitted model.** One point per classified unit (a case or a group). The
#' x-axis is the unit's entropy relative to a flat posterior: 0 means the
#' unit belongs to one class with certainty, 1 means it is equally likely to
#' belong to every class (dashed line). The rows are the assigned (most
#' likely) classes, the largest at the top, each in its own colour; the
#' points are spread within a row so that they do not overlap. A short
#' vertical bar marks the mean of each class. The subtitle gives the relative
#' entropy of the whole classification, which is one minus the mean of the
#' points: there, higher means sharper.
#'
#' A fit with two or more group classes besides its profiles is drawn at
#' both levels, as two plots. A level with one class has no entropy to draw,
#' so a fit with a single profile or class is refused. A fit with a noise
#' component is refused too.
#'
#' **An enumeration grid.** The x-axis is the number of classes (the number
#' of group classes for a grid of group-class families), the y-axis the
#' relative entropy of each candidate: higher is sharper. There is one line
#' per number of group classes, and per covariance structure when the grid
#' has several. A grid with group classes has one panel for the individuals
#' and one for the groups. A candidate with one class has no point.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. This covers fits with
#'     covariates, growth mixtures, and the additive and cross-level families
#'     of [multilpa()];
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()];
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text. When two plots are returned, both get the
#'   same title and subtitle.
#' @return A ggplot object. A fit that classifies at both levels returns a
#'   `latents_plots` list of two ggplot objects named `individuals` and
#'   `groups`; printing the list draws both.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit with one class or with a noise
#'     component, a grid in which no converged candidate has two or more
#'     classes, pooled imputations, a comparison, a bootstrap result or a
#'     sensitivity analysis. For an object of the wrong kind, the message names
#'     the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [get_entropy()] for the relative entropy as a table;
#'   [plot_posteriors()], [plot_avepp()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_entropy(fit)
#'
#' # A grid: the relative entropy of every candidate
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_entropy(grid)
#' @export
plot_entropy <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_entropy", substitute(x), main, subtitle)
}

#' Plot the posterior probability of each unit's assigned class, by class
#'
#' `plot_posteriors()` plots how certain the classification of each unit is: the
#' posterior probability of the class the unit is assigned to, by class.
#'
#' @details
#' One point per classified unit (a case or a group). The x-axis is the
#' posterior probability of the unit's assigned (most likely) class. It runs
#' from chance, one divided by the number of classes (dashed line), to 1,
#' which is certainty. The rows are the assigned classes, the largest at the
#' top, each in its own colour; the points are spread within a row so that
#' they do not overlap. A short vertical bar marks the mean of each class,
#' which is the class's average posterior probability.
#'
#' A fit with two or more group classes besides its profiles is drawn at
#' both levels, as two plots. A level with one class has nothing to draw, so
#' a fit with a single profile or class is refused. A fit with a noise
#' component is refused too. Enumeration grids have no such view; use
#' [plot_avepp()] or [plot_entropy()] for a grid.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. This covers fits with
#'     covariates, growth mixtures, and the additive and cross-level families
#'     of [multilpa()];
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text. When two plots are returned, both get the
#'   same title and subtitle.
#' @return A ggplot object. A fit that classifies at both levels returns a
#'   `latents_plots` list of two ggplot objects named `individuals` and
#'   `groups`; printing the list draws both.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit with one class or with a noise
#'     component, an enumeration grid, pooled imputations, a comparison, a
#'     bootstrap result or a sensitivity analysis. For an object of the wrong
#'     kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [get_classification()] for the average posterior probability of
#'   each class as a table; [plot_entropy()], [plot_avepp()];
#'   [plot_classification()] for every class probability of each person of a
#'   growth model.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_posteriors(fit)
#' @export
plot_posteriors <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_posteriors", substitute(x), main, subtitle)
}

#' Plot the average posterior probability matrix of a fit, or the lowest average
#' posterior probability of every candidate
#'
#' `plot_avepp()` plots the average posterior probability (AvePP) matrix of a
#' fitted model: how much of their probability the units assigned to each class
#' place on every class. For an enumeration grid it plots the lowest average
#' posterior probability of every candidate.
#'
#' @details
#' **A fitted model.** A square of tiles. Each row is the class the units are
#' assigned to (their most likely class), labelled with the number of units
#' assigned. Each column is a class. A tile holds the mean posterior
#' probability of the column's class among the units of the row; the number
#' is printed in the tile and the fill is darker for higher values. The
#' outlined diagonal is the average posterior probability of each class
#' among its own units: values near 1 mean that the class is well separated.
#' A fit with one class draws a single tile.
#'
#' A fit with two or more group classes besides its profiles is drawn at
#' both levels, as two plots. A fit with a noise component is refused.
#'
#' **An enumeration grid.** The x-axis is the number of classes (the number
#' of group classes for a grid of group-class families), the y-axis the
#' lowest average posterior probability of each candidate, that is, the
#' smallest value on its diagonal. There is one line per number of group
#' classes, and per covariance structure when the grid has several. A dashed
#' line marks the cutoff 0.8;
#' `plot_enumeration(x, what = "avepp", cutoff = )` draws the same view with
#' another cutoff. A grid with group classes has one panel for the
#' individuals and one for the groups.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. This covers fits with
#'     covariates, growth mixtures, and the additive and cross-level families
#'     of [multilpa()];
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()];
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text. When two plots are returned, both get the
#'   same title and subtitle.
#' @return A ggplot object. A fit that classifies at both levels returns a
#'   `latents_plots` list of two ggplot objects named `individuals` and
#'   `groups`; printing the list draws both.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit with a noise component, a grid in
#'     which no converged candidate has two or more classes, pooled
#'     imputations, a comparison, a bootstrap result or a sensitivity analysis.
#'     For an object of the wrong kind, the message names the plot functions
#'     that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [get_avepp()] for the matrix as a table; [plot_entropy()],
#'   [plot_posteriors()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_avepp(fit)
#'
#' # A grid: the lowest average posterior probability of every candidate
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_avepp(grid)
#' @export
plot_avepp <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_avepp", substitute(x), main, subtitle)
}

#' Plot profile means as tiles of profiles by indicators, or a grid's criterion
#' as tiles of candidates
#'
#' `plot_heatmap()` plots a grid of tiles: for a profile, class or transition
#' fit, one row per profile and one column per indicator; for an enumeration
#' grid, one tile per candidate, showing how far its information criterion is
#' above the lowest.
#'
#' @details
#' **A fit with continuous indicators.** Each tile is the profile's mean of
#' the indicator minus the indicator's observed mean, in observed standard
#' deviations (weighted ones for a fit with sampling weights). The value is
#' printed in the tile. Tiles below the mean are red, tiles above it blue,
#' and a tile at the mean is white. Categorical, ordinal and count
#' indicators of such a fit are not drawn; a caption names them.
#'
#' **A fit with categorical or ordinal indicators only.** The tiles are
#' response probabilities, darker for higher values. When every indicator
#' has two categories, the tile is the probability of the last category.
#' When an indicator has more, there is one panel per profile, with the
#' indicators in rows and the categories in columns.
#'
#' **An enumeration grid.** The x-axis is the number of classes and the
#' y-axis the number of group classes. Each tile is one candidate, labelled
#' with how far its criterion is above the lowest in the grid (`+12`); the
#' lowest candidate is labelled `lowest`, and a candidate that was not
#' fitted or did not converge is marked `x`. The fill is darker the closer
#' the candidate is to the lowest. The criterion is `bic_groups` (BIC with
#' the number of groups as the sample size) on a profile or class grid and
#' `bic` on a transition or regression grid. A grid with several covariance
#' structures has one panel per structure.
#' `plot_selection(x, what = "heatmap", using = )` draws the same view with
#' another criterion.
#'
#' @param x One of:
#'   * a profile or class fit from [lpa()], [lca()], [multilpa()] or
#'     [multilca()], with or without covariates;
#'   * a transition fit from [lta()] of class `multilpa_transitions`. The
#'     general transition fits (class `multilpa_lta`) have no heat map and
#'     are refused; `class(x)` tells which one [lta()] returned;
#'   * an enumeration grid from [enumerate_lpa()], [enumerate_lca()],
#'     [enumerate_classes()] or [enumerate_regressions()]. A grid of
#'     group-class families (`enumerate_classes(family = )`) is refused;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: an additive or cross-level fit, a mixture
#'     regression, a growth mixture, a general transition fit or a grid of
#'     group-class families. For an object of the wrong kind, the message names
#'     the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_profiles()], [plot_bars()] and [plot_responses()] for
#'   other views of the same means and probabilities; [plot_selection()]
#'   and [plot_enumeration()] for the other views of a grid.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_heatmap(fit)
#'
#' # Categorical indicators: response probabilities
#' classes <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
#'                n_classes = 2, n_starts = 2, seed = 1)
#' plot_heatmap(classes)
#'
#' # A grid: BIC of every candidate above the lowest
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_heatmap(grid)
#' @export
plot_heatmap <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_heatmap", substitute(x), main, subtitle)
}

#' Plot the mean of each continuous indicator in each profile as bars
#'
#' `plot_bars()` plots the mean of each continuous indicator in each profile as
#' bars that start at zero.
#'
#' @details
#' The x-axis lists the continuous indicators and the y-axis is the
#' estimated mean, in each indicator's own units. Within an indicator there
#' is one bar per profile, side by side. A profile has its own fill colour
#' and its own marker at the end of its bar, shown in the legend. When the
#' fit has standard errors, each bar carries its 95% interval; when it has
#' none, the subtitle says why. Categorical, ordinal and count indicators
#' are not drawn; a caption names them.
#'
#' `plot_profiles(x, what = "bars", scale = "standardized")` draws the same
#' bars in standard deviations, and `intervals = FALSE` there leaves the
#' intervals out.
#'
#' @param x One of:
#'   * a profile fit from [lpa()] or [multilpa()], with or without
#'     covariates, that has at least one continuous indicator;
#'   * a transition fit from [lta()] of class `multilpa_transitions` with at
#'     least one continuous indicator. The general transition fits (class
#'     `multilpa_lta`) are refused;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit without continuous indicators (a
#'     latent class fit), an additive or cross-level fit, a mixture regression,
#'     a growth mixture, a general transition fit or an enumeration grid. For
#'     an object of the wrong kind, the message names the plot functions that
#'     apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_profiles()] for the same means as lines and for the
#'   options of this view; [plot_heatmap()], [plot_raincloud()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_bars(fit)
#' @export
plot_bars <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_bars", substitute(x), main, subtitle)
}

#' Plot the distribution of each continuous indicator in each profile as a
#' density, its quartiles and the points
#'
#' `plot_raincloud()` plots the observed values of each continuous indicator,
#' split by profile: a density, a line over the middle half of the values (first
#' to third quartile) with a point at the median, and every case as a point.
#'
#' @details
#' There is one panel per continuous indicator, each with its own x-axis in
#' the indicator's units. The rows of a panel are the profiles, the largest
#' at the top, each in its own colour. For each profile the panel shows the
#' density of the indicator as a filled curve, the interquartile range as a
#' thick line with a white point at the median, and each case as a small
#' point below. A case counts for its most likely profile.
#'
#' The view draws the stored data of the fit. A fit with person or grand
#' centering or scaling is refused, because its stored data are not on the
#' scale of the model. For a transition fit, a case is one answered occasion
#' of a unit, counted for its most likely state.
#'
#' `plot_profiles(x, what = "raincloud", scale = "standardized")` draws the
#' same view in standard deviations.
#'
#' @param x One of:
#'   * a profile fit from [lpa()] or [multilpa()], with or without
#'     covariates, that has at least one continuous indicator and was fitted
#'     without centering or scaling;
#'   * a transition fit from [lta()] with homogeneous transitions (class
#'     `multilpa_transitions`) that has at least one continuous indicator;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit without continuous indicators (a
#'     latent class fit), a centred or scaled fit, a general transition fit
#'     (class `multilpa_lta`), an additive or cross-level fit, a mixture
#'     regression, a growth mixture or an enumeration grid. For an object of
#'     the wrong kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_parallel()] and [plot_pairs()] for the other views of the
#'   cases; [plot_profiles()] for the profile means.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_raincloud(fit)
#' @export
plot_raincloud <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_raincloud", substitute(x), main, subtitle)
}

#' Plot every case as a line across the continuous indicators, one panel per
#' profile
#'
#' `plot_parallel()` plots each case as one line that joins its values on the
#' continuous indicators, with one panel per profile (a parallel-coordinates
#' plot).
#'
#' @details
#' The x-axis lists the continuous indicators and the y-axis is the
#' standardized value (z), so that indicators on different scales can share
#' one axis. Each panel holds the cases whose most likely profile it is,
#' drawn in the profile's colour; a line is fainter the less certain the
#' case's assignment. The cases of the other profiles are drawn in grey
#' behind. The dark line with markers is the mean of the panel's profile.
#'
#' The view needs two or more continuous indicators and draws the stored
#' data of the fit. A fit with person or grand centering or scaling is
#' refused, because its stored data are not on the scale of the model.
#' For a transition fit, a case is one answered occasion of a unit, counted
#' for its most likely state.
#'
#' @param x One of:
#'   * a profile fit from [lpa()] or [multilpa()], with or without
#'     covariates, that has two or more continuous indicators and was fitted
#'     without centering or scaling;
#'   * a transition fit from [lta()] with homogeneous transitions (class
#'     `multilpa_transitions`) that has two or more continuous indicators;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit with fewer than two continuous
#'     indicators, a centred or scaled fit, a general transition fit (class
#'     `multilpa_lta`), an additive or cross-level fit, a mixture regression,
#'     a growth mixture or an enumeration grid. For an object of the wrong
#'     kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_raincloud()] and [plot_pairs()] for the other views of the
#'   cases; [plot_profiles()] for the profile means.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_parallel(fit)
#' @export
plot_parallel <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_parallel", substitute(x), main, subtitle)
}

#' Plot every pair of continuous indicators against each other, with the cases
#' and an ellipse for each profile
#'
#' `plot_pairs()` plots every pair of continuous indicators against each other,
#' with the cases marked by profile and an ellipse for each profile.
#'
#' @details
#' The plot is a matrix of panels with one row and one column per
#' continuous indicator. A panel off the diagonal is a scatter plot of two
#' indicators. Each case is a point with the colour and the marker of its
#' most likely profile, drawn larger the less certain its assignment. Each
#' profile has an ellipse that holds 95% of its fitted distribution, drawn
#' around the profile's mean, which is marked. A panel on the diagonal shows
#' the density of that indicator among the cases of each profile, scaled by
#' the size of the profile.
#'
#' The view needs two or more continuous indicators and draws the stored
#' data of the fit. A fit with person or grand centering or scaling is
#' refused, because its stored data are not on the scale of the model.
#' For a transition fit, a case is one answered occasion of a unit, counted
#' for its most likely state.
#'
#' @param x One of:
#'   * a profile fit from [lpa()] or [multilpa()], with or without
#'     covariates, that has two or more continuous indicators and was fitted
#'     without centering or scaling;
#'   * a transition fit from [lta()] with homogeneous transitions (class
#'     `multilpa_transitions`) that has two or more continuous indicators;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fit with fewer than two continuous
#'     indicators, a centred or scaled fit, a general transition fit (class
#'     `multilpa_lta`), an additive or cross-level fit, a mixture regression,
#'     a growth mixture or an enumeration grid. For an object of the wrong
#'     kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_raincloud()] and [plot_parallel()] for the other views of
#'   the cases; [plot_profiles()] for the profile means.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_pairs(fit)
#' @export
plot_pairs <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_pairs", substitute(x), main, subtitle)
}

#' Plot the curve of the outcome in each class of a mixture regression model
#'
#' `plot_trajectories()` plots the curve of the outcome in each class of a model
#' fitted with [mixture_regression()]: over time when the classes are classes of
#' persons, over a predictor otherwise.
#'
#' @details
#' **Classes of persons.** This is a growth mixture (a fit with `random`) or
#' a fit with `class_level = "group"`. The x-axis is the time variable: the
#' first numeric variable among the random-effect terms and then the
#' predictors. The y-axis is the outcome. For each class the plot shows:
#'
#' * a thick line, the mean trajectory of the class, with the class named at
#'   its right end;
#' * a dark band, the 95% confidence band of that mean;
#' * for a growth mixture, a light band in which 80% of the trajectories of
#'   the class's persons lie;
#' * points with a range, the observed mean of the class at each time, with
#'   its interval (drawn when there are at most 30 distinct times);
#' * thin lines, the observed trajectories of up to 120 persons, coloured by
#'   their most likely class.
#'
#' `plot(x, what = "trajectories")` draws the same view and takes `time`
#' (another time variable), `facet = TRUE` (one panel per class),
#' `max_persons`, `level` and, for a growth mixture, `spread`.
#'
#' **Classes of observations.** For any other mixture regression the x-axis
#' is the first numeric predictor of the formula and the y-axis the outcome.
#' Each row of the data is a point with the colour and marker of its most
#' likely class, and each class has one fitted line, drawn with the other
#' predictors held at their means (factors at their first level).
#' `plot(x, what = "fitted", predictor = )` draws the same view over another
#' numeric predictor.
#'
#' @param x A fit from [mixture_regression()] of any kind, growth mixtures
#'   included, or a [diagnostics()] result, which is drawn from the fit it
#'   checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a profile, class, transition, additive or
#'     cross-level fit, or an enumeration grid. For an object of the wrong
#'     kind, the message names the plot functions that apply to it;
#'   * `latents_bad_argument` when a mixture regression that classifies
#'     observations has no numeric predictor in its formula;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_individuals()] for single persons, [plot_coefficients()]
#'   for the coefficients behind the curves, [plot_classification()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_trajectories(fit)
#'
#' # Classes of observations: the outcome against the first numeric predictor
#' regression <- mixture_regression(score ~ hours, data = study_hours,
#'                                  n_classes = 2, n_starts = 2, seed = 1)
#' plot_trajectories(regression)
#' @export
plot_trajectories <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_trajectories", substitute(x), main, subtitle)
}

#' Plot individual persons' observations and predicted curves beside the
#' trajectory of their class
#'
#' `plot_individuals()` plots one panel per person of a growth model, with the
#' person's observations, the person's own predicted curve and the trajectory of
#' the person's class.
#'
#' @details
#' The x-axis is the time variable and the y-axis the outcome, as in
#' [plot_trajectories()]. In each panel the points are the person's
#' observations, the solid line is the curve predicted for that person and
#' the dashed line is the mean trajectory of the person's most likely class.
#' All three are drawn in the colour of that class. The panel title gives
#' the name of the identifier column with the person's identifier, then the
#' person's class and the posterior probability of that class.
#'
#' @param x One of:
#'   * a growth mixture: a fit from [mixture_regression()] with `random`;
#'   * a fit from [mixture_regression()] with `class_level = "group"`, which
#'     classifies persons;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param persons Which persons to draw:
#'   * one number (default `12`): how many persons. They are taken from the
#'     classes in turn, and within a class they are spread from its most
#'     certain to its least certain member;
#'   * two or more values, or text: the identifiers of the persons to draw,
#'     as they appear in the data, for example `c(3, 17, 42)`. A single
#'     numeric identifier must be given as text, `persons = "1001"`,
#'     because a single number is read as a count.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object with one panel per person.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a mixture regression that classifies
#'     observations (including one with `n_group_classes` of two or more), a
#'     fit of any other family or an enumeration grid. For an object of the
#'     wrong kind, the message names the plot functions that apply to it;
#'   * `latents_bad_argument` when `persons` is a number below 1 or contains a
#'     missing value, or when none of the identifiers is in the fitted data;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_trajectories()] for the classes, [plot_classification()]
#'   for the class probabilities of every person.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_individuals(fit)
#' plot_individuals(fit, persons = 4)
#' plot_individuals(fit, persons = c(3, 17, 42))
#' @export
plot_individuals <- function(x, persons = 12L, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_individuals", substitute(x), main, subtitle,
                     persons = persons)
}

#' Plot each person's posterior probability of every class in a growth model as
#' a stacked bar
#'
#' `plot_classification()` plots one stacked bar per person of a growth model,
#' split by the person's posterior probability of each class.
#'
#' @details
#' The x-axis holds the persons and the y-axis the posterior probability,
#' from 0 to 1. Each bar is one person and is divided into one segment per
#' class, filled in the class's colour, so that a bar of a single colour is
#' a person classified with certainty. There is one panel per assigned (most
#' likely) class, the largest first; within a panel the persons are ordered
#' from the most to the least certain. On the bar of the least certain
#' person of each panel the segments are labelled with their class, where
#' they are wide enough. The subtitle gives, for each class, the average
#' probability of the assigned class.
#'
#' [plot_posteriors()] shows only the probability of the assigned class, and
#' takes fits of every family.
#'
#' @param x One of:
#'   * a growth mixture: a fit from [mixture_regression()] with `random`;
#'   * a fit from [mixture_regression()] with `class_level = "group"`, which
#'     classifies persons;
#'   * a [diagnostics()] result, which is drawn from the fit it checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a mixture regression that classifies
#'     observations (including one with `n_group_classes` of two or more), a
#'     fit of any other family or an enumeration grid. For an object of the
#'     wrong kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_posteriors()] for the certainty of the assigned class
#'   alone; [plot_individuals()]; [get_classification()] for the table.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_classification(fit)
#' @export
plot_classification <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_classification", substitute(x), main, subtitle)
}

#' Plot the regression coefficients of each class of a mixture regression, with
#' 95% confidence intervals
#'
#' `plot_coefficients()` plots the regression coefficients of each class of a
#' model fitted with [mixture_regression()], each with its 95% confidence
#' interval.
#'
#' @details
#' **Classes of persons** (a growth mixture, or a fit with
#' `class_level = "group"`). There is one panel per term of the model, each
#' with its own x-axis. The rows of a panel are the classes, and a row
#' `Shared by all classes` holds a term that is common to them. A point is
#' the estimate and the horizontal bar its interval; each class has its own
#' colour and marker. The estimate, the interval and the Wald p-value are
#' printed beside each row, and a dashed line marks zero.
#'
#' **Classes of observations** (any other mixture regression). There is one
#' panel. The rows are the terms of the model; within a term the classes
#' are drawn one under the other, each with its own colour and marker. A
#' point is the estimate and the horizontal line its interval, and a dotted
#' line marks zero.
#'
#' `plot(x, what = "coefficients", level = )` draws the same view with
#' another confidence level.
#'
#' @param x A fit from [mixture_regression()] of any kind, growth mixtures
#'   included, or a [diagnostics()] result, which is drawn from the fit it
#'   checked.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a profile, class, transition, additive or
#'     cross-level fit, or an enumeration grid. For an object of the wrong
#'     kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_trajectories()] for the curves the coefficients give;
#'   [get_results()] for the coefficients as a table.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_coefficients(fit)
#' @export
plot_coefficients <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_coefficients", substitute(x), main, subtitle)
}

#' Plot every information criterion of a grid, or of a model comparison, against
#' the number of classes
#'
#' `plot_criteria()` plots every information criterion of an enumeration grid
#' against the number of classes, or the criteria of the models compared by
#' [compare_models()].
#'
#' @details
#' **An enumeration grid.** There is one panel per criterion, each with its
#' own y-axis; lower values are preferred. The x-axis is the number of
#' classes (the number of group classes for a grid of group-class families,
#' or when the grid has one number of profiles and several numbers of group
#' classes). There is one line per number of group classes, and per
#' covariance structure when the grid has several. In each panel the lowest
#' value among the converged candidates is circled, and a candidate that did
#' not converge is marked with a cross at the bottom of the panel. The
#' circle marks a minimum; it does not choose a model.
#'
#' The criteria depend on the kind of grid:
#'
#' * a profile or class grid: AIC, KIC, BIC, SABIC, CAIC, AWE, ICL and CLC.
#'   BIC, SABIC, CAIC, AWE and ICL exist in two versions, with the number of
#'   groups or the number of individuals as the sample size; the two are
#'   drawn as separate panels where they differ and as one panel where they
#'   are the same;
#' * a transition grid or a grid of group-class families: AIC, BIC and BIC
#'   with the individuals as the sample size;
#' * a regression grid: AIC, BIC, BIC with the rows as the sample size,
#'   SABIC and ICL.
#'
#' On these grids a criterion that takes the same values as one already
#' drawn is not repeated.
#'
#' [plot_enumeration()] draws the criteria you name, and takes the options
#' of this view (`labels`, `mark_minimum`, `combine`).
#'
#' **A comparison.** The x-axis lists the compared models and the y-axis is
#' the value of the criterion; lower is better. AIC, BIC and SABIC are drawn
#' as one line each, with their own colour and marker.
#'
#' @param x One of:
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()];
#'   * a comparison from [compare_models()].
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fitted model, a diagnostics result, pooled
#'     imputations, a bootstrap result or a sensitivity analysis. For an object
#'     of the wrong kind, the message names the plot functions that apply to
#'     it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [plot_enumeration()] to choose the criteria; [get_candidates()]
#'   and [get_best()] for the values as tables; [plot_deltas()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_criteria(grid)
#'
#' # A comparison of two fitted models
#' one <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 1,
#'            n_starts = 2, seed = 1)
#' two <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_criteria(compare_models(one, two))
#' @export
plot_criteria <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_criteria", substitute(x), main, subtitle)
}

#' Plot how much the information criterion falls each time a class or a group
#' class is added
#'
#' `plot_deltas()` plots how much the information criterion of an enumeration
#' grid falls each time one class, or one group class, is added.
#'
#' @details
#' The plot has up to two panels. In the left panel the x-axis is the number
#' of classes after the step, and each point is the decrease in the
#' criterion from the candidate with one class fewer and the same number of
#' group classes. There is one line per number of group classes, named at
#' its right end. In the right panel the x-axis is the number of group
#' classes after the step, and each point is the decrease from the candidate
#' with one group class fewer and the same number of classes; a grid with
#' one number of group classes has no such panel. A point above the zero
#' line means that the larger model has the lower criterion. A small
#' decrease after a large one is where further classes stop paying.
#'
#' The criterion is `bic_groups` (BIC with the number of groups as the
#' sample size) on a profile or class grid and `bic` on a transition or
#' regression grid. `plot_selection(x, what = "deltas", using = )` draws the
#' same view with another criterion.
#'
#' @param x An enumeration grid from [enumerate_lpa()], [enumerate_lca()],
#'   [enumerate_classes()] or [enumerate_regressions()]. A grid of
#'   group-class families (`enumerate_classes(family = )`) has no such view
#'   and is refused.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fitted model, a comparison or another
#'     result, a grid of group-class families, or a grid in which no converged
#'     candidate has a converged neighbour one step smaller. For an object of
#'     the wrong kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [get_deltas()] and [get_knees()] for the numbers of a profile or
#'   class grid; [plot_selection()] for the other selection views.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_deltas(grid)
#' @export
plot_deltas <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_deltas", substitute(x), main, subtitle)
}

#' Plot whether each candidate of a grid converged, lies on a boundary and
#' replicated its best log likelihood
#'
#' `plot_trust()` plots whether the estimate of each candidate of an enumeration
#' grid can be relied on: did it converge, is it on a boundary, and did more
#' than one start reach its best log likelihood.
#'
#' @details
#' The x-axis is the number of classes and the y-axis the number of group
#' classes, so each tile is one candidate. A grid with several covariance
#' structures has one panel per structure. The fill and the text of a tile
#' give the status of the candidate:
#'
#' * `replicated`: converged, and two or more starts reached the best log
#'   likelihood;
#' * `once`: converged, but only one start reached the best log likelihood,
#'   so a higher maximum may exist;
#' * `n/a`: converged, with no replication count stored;
#' * `boundary`: converged with an estimate on the boundary of the
#'   parameter space;
#' * `not conv.`: the best start did not converge;
#' * `not fitted`: the candidate has no log likelihood.
#'
#' Under the status the tile gives the number of starts that reached the
#' best log likelihood, for example `2x`. Replication does not prove that
#' the global maximum was found.
#'
#' @param x An enumeration grid from [enumerate_lpa()], [enumerate_lca()],
#'   [enumerate_classes()] or [enumerate_regressions()]. A grid of
#'   group-class families (`enumerate_classes(family = )`) has no such view
#'   and is refused.
#' @param main,subtitle Title and subtitle of the plot. `NULL` (default)
#'   keeps the view's own text.
#' @return A ggplot object.
#'
#'   Errors:
#'   * `latents_nothing_to_plot` when the package returned `x` but this
#'     function has no plot of it: a fitted model, a comparison or another
#'     result, or a grid of group-class families. For an object of the wrong
#'     kind, the message names the plot functions that apply to it;
#'   * `latents_not_a_fit` when `x` is not a fit, grid or other result of the
#'     package, for example a data frame;
#'   * `latents_missing_package` when ggplot2 is not installed.
#' @seealso [get_replication()] for the counts of starts as a table;
#'   [get_candidates()] for `converged` and `boundary`; [plot_selection()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_trust(grid)
#' @export
plot_trust <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_trust", substitute(x), main, subtitle)
}

#' Plot every view an object has and return the plots as one list
#'
#' `plot_all()` plots every view that an object of the package has, and returns
#' the plots as one list. Printing the list shows each plot in turn.
#'
#' @details
#' Which plots the list holds depends on the object. The names of the list
#' are the names of the views:
#'
#' * a profile, class or transition fit of [lpa()], [lca()], [multilpa()],
#'   [multilca()] or [lta()] (class `multilpa_transitions`): the views its
#'   data allow among `profiles`, `bars`, `heatmap`, `raincloud`,
#'   `parallel`, `pairs`, `responses`, `probabilities`, `sequences`,
#'   `transitions`, `sizes`, `entropy`, `posteriors` and `avepp`. A fit with
#'   two or more group classes adds `group_sizes`, `group_posteriors`,
#'   `group_entropy` and `group_avepp`;
#' * a fit of another family: the views of its own `plot()`, then the four
#'   classification views `sizes`, `posteriors`, `entropy` and `avepp`. A
#'   mixture regression has `fitted` and `coefficients`; one that classifies
#'   persons also `trajectories`, `individuals` and `classification`; a
#'   growth mixture has those three, `coefficients` and `random`; a general
#'   transition fit `transitions`, `profiles` and `sequences`; an additive
#'   fit `means`, `variances` and `intercepts`; a cross-level fit
#'   `profiles`, `group_means` and `composition`. When such a fit classifies
#'   at two levels, each classification view is drawn twice, as
#'   `sizes_individuals`, `sizes_groups` and so on;
#' * pooled imputations: the pooled estimates of each level the pooled table
#'   has, among `profile`, `group` and `measurement`;
#' * an enumeration grid: `criteria` (every criterion, as [plot_criteria()]),
#'   `entropy`, `avepp`, `deltas`, `heatmap`, `smallest`, `trust` and
#'   `tree`. A grid of group-class families has `criteria`, `entropy`,
#'   `avepp` and `tree`. A grid with several trees has one plot per tree,
#'   named `tree_` and the tree's name;
#' * a comparison, a bootstrap test, a bootstrap table of
#'   [parameter_inference()] or a sensitivity analysis: its one plot, named
#'   `plot`;
#' * a [diagnostics()] result: the plots of the fit it checked.
#'
#' Only the views the object can draw are returned. A view that was tried
#' and refused is left out and named in a message
#' (`Not drawn for this object: ...`).
#'
#' @param x Any object returned by the package: a fitted model, pooled
#'   imputations, an enumeration grid, a comparison, a bootstrap result, a
#'   sensitivity analysis or a diagnostics result.
#' @return A `latents_plots` list of ggplot objects, named by view. Printing
#'   it draws every plot; in an interactive session R waits for a key press
#'   before each new plot.
#'
#'   Errors: `latents_not_a_fit` when `x` is not a fit, grid or other result
#'   of the package, for example a data frame; `latents_nothing_to_plot`
#'   when none of the object's views can be drawn;
#'   `latents_missing_package` when ggplot2 is not installed.
#' @seealso [overview()] prints the summary and the checks as well; the
#'   functions named after single views, such as [plot_profiles()],
#'   [plot_sizes()] and [plot_criteria()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_all(fit)
#' @export
plot_all <- function(x) {
  if (inherits(x, "multilpa_diagnostics")) x <- x$fit
  kind <- .latents_wild_kind(x)
  if (is.na(kind)) .latents_wild_refuse_object(x, "plot_all", substitute(x))
  .gg_require()
  draws <- .latents_wild_all_draws(x, kind)
  # A view can still refuse for a reason only its method knows (a panel with
  # nothing in it); the others are drawn and the refused ones named, as
  # report() and plot(what = "all") do.
  refused <- \(condition) NULL
  plots <- lapply(draws, \(draw) tryCatch(draw(),
    latents_no_plot = refused, latents_nothing_to_plot = refused,
    latents_no_time = refused, latents_no_categorical = refused,
    latents_no_continuous = refused, latents_no_indicator_data = refused,
    latents_unsupported_view = refused, latents_no_group_classes = refused))
  drawn <- !vapply(plots, is.null, logical(1))
  if (!all(drawn)) {
    message(sprintf("Not drawn for this object: %s.",
                    paste(names(draws)[!drawn], collapse = ", ")))
  }
  if (!any(drawn)) {
    stop(errorCondition("plot_all() has no view this object can draw.",
                        class = "latents_nothing_to_plot", call = NULL))
  }
  # A view that is itself a list of plots (a tree per group-class count)
  # contributes each of its plots, named after the view.
  plots <- plots[drawn]
  flat <- lapply(names(plots), \(name) {
    p <- plots[[name]]
    if (!inherits(p, "latents_plots")) return(stats::setNames(list(p), name))
    stats::setNames(unclass(p), paste(name, names(p) %||% seq_along(p), sep = "_"))
  })
  .gg_plots(do.call(c, flat))
}

#' Every view of an object, each as a function that draws it
#' @param x Any object of a known kind.
#' @param kind Its kind.
#' @return A named list of functions of no arguments.
#' @noRd
.latents_wild_all_draws <- function(x, kind) {
  if (identical(kind, "grid")) {
    views <- plot_views(x)
    # Every criterion is one panel plot ("all"); "diagnostics" repeats the
    # entropy and AvePP views and "selection" repeats the rest.
    keep <- views$supported %in% TRUE & (!views$group %in% "criteria" |
                                           views$type %in% "all") &
      !views$type %in% c("diagnostics", "selection")
    types <- views$type[keep]
    names(types) <- ifelse(types == "all", "criteria", types)
    return(lapply(types, \(view) {
      \() if (identical(view, "tree")) plot_tree(x) else plot_enumeration(x, what = view)
    }))
  }
  if (kind %in% c("comparison", "bootstrap_lrt", "bootstrap", "sensitivity")) {
    return(list(plot = \() plot(x)))
  }
  if (identical(kind, "profile_model")) {
    views <- .multilpa_supported_views(x)
    draws <- lapply(stats::setNames(views, views), \(view) {
      \() plot(x, what = view)
    })
    if (.gg_has_group_classes(x) && !isTRUE(x$noise)) {
      group <- paste0("group_", c("sizes", "posteriors", "entropy", "avepp"))
      draws <- c(draws, lapply(stats::setNames(group, group), \(view) {
          \() plot_diagnostics(x, what = view)
      }))
    }
    return(draws)
  }
  views <- .latents_plot_offered_views(x)
  # A pooled fit draws only the levels its pooled table has (no group view
  # with one group class), so no view is tried that is bound to refuse.
  if (identical(kind, "pooled")) views <- intersect(views, unique(get_results(x)$level))
  classification <- c("sizes", "posteriors", "entropy", "avepp")
  # The family's own "posteriors" view is not the certainty view (see
  # .latents_wild_plot()); the classification views come from the posteriors.
  if (inherits(x, .latents_model_classes)) views <- setdiff(views, classification)
  draws <- if (length(views) == 0L) list(plot = \() plot(x)) else
    lapply(stats::setNames(views, views), \(view) \() plot(x, what = view))
  if (!inherits(x, .latents_model_classes)) return(draws)
  c(draws, lapply(stats::setNames(classification, classification), \(view)
    \() .latents_wild_classes(x, view, "plot_all", NULL, NULL)))
}

#' Read one diagnostic table from any fit, grid or diagnostics result
#' @param x The object.
#' @param table The table's name in a diagnostics result.
#' @param verb The calling verb, for messages.
#' @param expression The caller's `substitute(x)`.
#' @return A data frame.
#' @noRd
.latents_wild_table <- function(x, table, verb, expression) {
  if (inherits(x, "multilpa_diagnostics")) {
    if (is.null(x[[table]])) .latents_wild_refuse_table(x$fit, table, verb)
    return(get_results(x, table))
  }
  if (inherits(x, .latents_model_classes)) {
    quality <- diagnostics(x, plots = FALSE)
    if (is.null(quality[[table]])) .latents_wild_refuse_table(x, table, verb)
    return(get_results(quality, table))
  }
  if (!inherits(x, .latents_enumeration_classes)) {
    label <- .latents_argument_label(expression, "x")
    stop(errorCondition(sprintf(paste(
      "%s() takes a fitted model, an enumeration grid or a diagnostics result,",
      "but `%s` is %s."), verb, label, .latents_describe(x)),
      class = c("latents_not_a_fit", "latents_bad_argument"), call = NULL))
  }
  fits <- x$fits
  names(fits) <- .latents_wild_candidate_names(x)
  # A failed candidate is stored as NULL (profile grids) or as the condition
  # that stopped it (regression grids); only fitted models have tables.
  fits <- Filter(\(fit) inherits(fit, .latents_model_classes), fits)
  tables <- lapply(names(fits), \(candidate) {
    quality <- diagnostics(fits[[candidate]], plots = FALSE)
    # Read the table only where this candidate's family defines it, so a
    # grid mixing candidates with and without it stacks the ones that have it.
    if (is.null(quality[[table]])) return(NULL)
    rows <- get_results(quality, table)
    if (nrow(rows) == 0L) return(NULL)
    cbind(data.frame(candidate = candidate, stringsAsFactors = FALSE), rows)
  })
  tables <- Filter(Negate(is.null), tables)
  if (length(tables) == 0L) {
    stop(errorCondition(sprintf(
      "%s(): no candidate of this grid has a `%s` table.", verb, table),
      class = "latents_bad_argument", call = NULL))
  }
  result <- do.call(rbind, tables)
  row.names(result) <- NULL
  result
}

#' Refuse a table a fit's family does not define, naming the ones it has
#' @param x A fitted model.
#' @param table The table's name in a diagnostics result.
#' @param verb The calling verb.
#' @return Does not return. A fit with a noise component raises
#'   `latents_unsupported_noise`, any other `latents_bad_argument`.
#' @noRd
.latents_wild_refuse_table <- function(x, table, verb) {
  getters <- c(checks = "get_checks", entropy = "get_entropy",
               classification = "get_classification",
               average_posteriors = "get_avepp", residuals = "get_residuals",
               replication = "get_replication")
  quality <- diagnostics(x, plots = FALSE)
  has <- getters[vapply(names(getters), \(name) !is.null(quality[[name]]), logical(1))]
  noise <- isTRUE(x$noise)
  stop(errorCondition(sprintf(
    "%s() has no table for %s%s. For this object: %s.", verb, .latents_describe(x),
    if (noise) " with a noise component" else "",
    paste(sprintf("%s()", has), collapse = ", ")),
    class = c(if (noise) "latents_unsupported_noise", "latents_bad_argument"),
    call = NULL))
}

#' Candidate names of a grid, in the order of its stored fits
#' @param x An enumeration grid.
#' @return A character vector, one name per element of `x$fits`.
#' @noRd
.latents_wild_candidate_names <- function(x) {
  if (!is.null(names(x$fits))) return(names(x$fits))
  # A grid of mixture regressions stores its fits unnamed, one per row of its
  # table, and names them by their class counts.
  grid <- x$table
  groups <- grid$n_group_classes %||% rep(1L, nrow(grid))
  sprintf("classes%d_groups%d", as.integer(grid$n_classes), as.integer(groups))
}

#' Get every diagnostic check of a fitted model, with its value, rule and status
#'
#' `get_checks()` returns the checks that [diagnostics()] runs on a fitted
#' model, as one table: each check with its value, the rule it is judged by and
#' its status. The checks cover the estimation (convergence, replication of the
#' best log likelihood, failed starts, estimates at a bound), the class sizes,
#' the quality of the classification (relative entropy, average posterior
#' probability, odds of correct classification), local independence (the
#' bivariate residuals) and the information criteria.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. The table is the one
#'     [diagnostics()] computes with its defaults;
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#'     The tables of the fitted candidates are stacked, with a first column
#'     `candidate` that names the candidate of each row. A candidate whose
#'     estimation failed has no rows;
#'   * the result of [diagnostics()]. The table is the one the result holds.
#' @return A data frame with one row per check (per candidate and check for
#'   a grid) and the columns:
#'   * `domain`: the group of the check, one of `"Estimation"`,
#'     `"Class size"`, `"Classification"`, `"Local independence"` and
#'     `"Information"`;
#'   * `check`: what is checked;
#'   * `level`: `"model"` for a check of the whole fit, `"individuals"` or
#'     `"groups"` for a check of one classification level;
#'   * `value`: the result as a number;
#'   * `display`: the result as formatted text;
#'   * `criterion`: the rule the status follows; empty when there is none;
#'   * `status`: `"ok"` (passed), `"info"` (reported, not judged),
#'     `"check"` (look at it), `"fail"` (failed) or `"n/a"` (not available
#'     for this fit);
#'   * `detail`: a note, for example the pair of indicators with the
#'     smallest adjusted p-value;
#'   * `reference`: the source of the check.
#'
#'   Errors: `latents_not_a_fit` when `x` is none of these, for example a
#'   comparison, pooled imputations or a data frame.
#' @seealso [diagnostics()] for the checks and their rules; [get_entropy()],
#'   [get_classification()], [get_avepp()], [get_replication()] and
#'   [get_residuals()] for the tables behind them; [overview()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_checks(fit)
#' @export
get_checks <- function(x) {
  .latents_wild_table(x, "checks", "get_checks", substitute(x))
}

#' Get the relative entropy, from 0 to 1, of how sharply a model separates its
#' classes
#'
#' `get_entropy()` returns one number per classification level that says how
#' sharply the model separates its classes: the relative entropy, from 0 (every
#' unit is equally likely to be in every class) to 1 (every unit belongs to one
#' class with certainty).
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. The table is the one
#'     [diagnostics()] computes with its defaults;
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#'     The tables of the fitted candidates are stacked, with a first column
#'     `candidate` that names the candidate of each row. A candidate whose
#'     estimation failed has no rows;
#'   * the result of [diagnostics()]. The table is the one the result holds.
#' @return A data frame with one row per classification level (per
#'   candidate and level for a grid) and the columns:
#'   * `level`: `"individuals"` (the cases sorted into profiles or classes)
#'     or `"groups"` (the groups sorted into group classes);
#'   * `n_classes`: the number of classes at that level;
#'   * `n_units`: the number of classified units;
#'   * `entropy_sum`: the sum of the units' entropies;
#'   * `relative_entropy`: one minus the mean entropy of a unit divided by
#'     its largest possible value, the logarithm of `n_classes`. `NA` for a
#'     level with one class.
#'
#'   A profile, class or transition fit always has both rows; without group
#'   classes its `"groups"` row has one class and `NA`. A fit of another
#'   family has a row for each level it classifies.
#'
#'   Errors: `latents_not_a_fit` when `x` is none of these, for example a
#'   comparison, pooled imputations or a data frame.
#' @seealso [plot_entropy()] draws the entropy of each unit;
#'   [get_classification()], [get_avepp()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_entropy(fit)
#' @export
get_entropy <- function(x) {
  .latents_wild_table(x, "entropy", "get_entropy", substitute(x))
}

#' Get how many units each class holds and how certain their assignment is
#'
#' `get_classification()` returns, for every class, how many units it holds,
#' counted two ways, and how certain the assignment of its units is.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. The table is the one
#'     [diagnostics()] computes with its defaults;
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#'     The tables of the fitted candidates are stacked, with a first column
#'     `candidate` that names the candidate of each row. A candidate whose
#'     estimation failed has no rows;
#'   * the result of [diagnostics()]. The table is the one the result holds.
#' @return A data frame with one row per level and class (per candidate,
#'   level and class for a grid) and the columns:
#'   * `level`: `"individuals"` or `"groups"`;
#'   * `class`: the number of the class; `0` is the noise component of a fit
#'     that has one;
#'   * `n_modal`: the number of units whose most likely class this is;
#'   * `proportion_modal`: `n_modal` as a share of the units;
#'   * `estimated_n`: the sum of the units' posterior probabilities of the
#'     class;
#'   * `estimated_proportion`: `estimated_n` as a share of the units;
#'   * `average_posterior`: the mean posterior probability of the class
#'     among the units assigned to it;
#'   * `odds_correct_classification`: the odds of `average_posterior`
#'     divided by the odds of `estimated_proportion`.
#'
#'   A fit without group classes has rows for the individuals only; a fit
#'   that classifies groups only has rows for the groups only.
#'
#'   Errors: `latents_not_a_fit` when `x` is none of these, for example a
#'   comparison, pooled imputations or a data frame.
#' @seealso [plot_sizes()] and [plot_posteriors()] draw it; [get_avepp()]
#'   for the full matrix of average posterior probabilities.
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_classification(fit)
#' @export
get_classification <- function(x) {
  .latents_wild_table(x, "classification", "get_classification", substitute(x))
}

#' Get the average posterior probability of every class among the units assigned
#' to each class
#'
#' `get_avepp()` returns the average posterior probability (AvePP) matrix in
#' long form: for the units assigned to each class, their mean posterior
#' probability of every class.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. The table is the one
#'     [diagnostics()] computes with its defaults;
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#'     The tables of the fitted candidates are stacked, with a first column
#'     `candidate` that names the candidate of each row. A candidate whose
#'     estimation failed has no rows;
#'   * the result of [diagnostics()]. The table is the one the result holds.
#' @return A data frame with one row per level, assigned class and class
#'   (and per candidate for a grid) and the columns:
#'   * `level`: `"individuals"` or `"groups"`;
#'   * `assigned_class`: the class the units are assigned to (their most
#'     likely class);
#'   * `class`: the class whose probability is averaged;
#'   * `n_assigned`: the number of units assigned to `assigned_class`;
#'   * `average_posterior`: the mean posterior probability of `class` among
#'     those units.
#'
#'   The rows in which `assigned_class` equals `class` are the diagonal of
#'   the matrix: values near 1 mean that the class is well separated.
#'
#'   Errors: `latents_not_a_fit` when `x` is none of these, for example a
#'   comparison, pooled imputations or a data frame.
#' @seealso [plot_avepp()] draws the matrix; [get_classification()],
#'   [get_entropy()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_avepp(fit)
#' @export
get_avepp <- function(x) {
  .latents_wild_table(x, "average_posteriors", "get_avepp", substitute(x))
}

#' Get the bivariate residuals: the association between pairs of indicators that
#' each class leaves unexplained
#'
#' `get_residuals()` returns the bivariate residuals, a check of local
#' independence: for each pair of indicators within each class, the association
#' that the model leaves unexplained.
#'
#' @details
#' The statistic depends on the two indicators (column `kind`):
#'
#' * `"gaussian"`, two continuous indicators: the residual correlation,
#'   with a normal test;
#' * `"categorical"`, two categorical or ordinal indicators: a chi-square
#'   test of the observed two-way table against the one the model expects;
#' * `"mixed"` or `"count"`, a continuous with a categorical indicator or a
#'   pair with a count indicator: a correlation, without a test.
#'
#' The p-values are adjusted across the pairs by the Benjamini-Hochberg
#' method. A fit with sampling weights reports the associations without
#' tests. `?get_results` describes the statistics in full.
#'
#' Only profile, class and transition fits have these residuals: fits of
#' [lpa()], [lca()], [multilpa()] and [multilca()], with or without
#' covariates, and transition fits of [lta()] of class
#' `multilpa_transitions`.
#'
#' @param x One of:
#'   * a profile, class or transition fit as described above. The table is
#'     the one [diagnostics()] computes with its defaults, which assesses
#'     each class separately;
#'   * an enumeration grid of such fits. The tables of the fitted
#'     candidates are stacked, with a first column `candidate`;
#'   * the result of [diagnostics()]. The table is the one the result
#'     holds, so `diagnostics(fit, by = "overall")` gives the residuals
#'     pooled over the classes.
#' @return A data frame with one row per class and pair of indicators (and
#'   per candidate for a grid) and the columns:
#'   * `profile`: the class, as `profile_1`, `profile_2`, ...; `overall`
#'     for residuals pooled over the classes;
#'   * `indicator_1`, `indicator_2`: the pair;
#'   * `kind`: the kind of pair, as listed above;
#'   * `observed`, `expected`: the observed association and the one the
#'     model implies;
#'   * `residual`: the residual association;
#'   * `effective_n`: the effective number of units behind it;
#'   * `statistic`, `df`, `p_value`: the test, `NA` where there is none;
#'   * `p_adjusted`: the p-value adjusted across the pairs.
#'
#'   Errors: `latents_bad_argument` when the fit is of a family without
#'   bivariate residuals (a mixture regression, a growth mixture, an
#'   additive or cross-level fit, or a general transition fit of class
#'   `multilpa_lta`), and when no candidate of a grid has them; the message
#'   names the tables the object does have. `latents_unsupported_noise`
#'   when the fit has a noise component. `latents_not_a_fit` when `x` is
#'   not a fit, a grid or a diagnostics result.
#' @seealso [get_results()] for the statistics; [diagnostics()] and
#'   [get_checks()] for the check that reads this table.
#' @examples
#' fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
#'            n_classes = 2, n_starts = 2, seed = 1)
#' get_residuals(fit)
#' @export
get_residuals <- function(x) {
  .latents_wild_table(x, "residuals", "get_residuals", substitute(x))
}

#' Get how many random starts converged, failed and reached the best log
#' likelihood
#'
#' `get_replication()` returns how many of a model's random starts ran,
#' converged and failed, and how many arrived at the best log likelihood. A best
#' solution that only one start reached may not be the highest one.
#'
#' @param x One of:
#'   * a fitted model of any family: from [lpa()], [lca()], [multilpa()],
#'     [multilca()], [lta()] or [mixture_regression()]. The table is the one
#'     [diagnostics()] computes with its defaults;
#'   * an enumeration grid of any kind, from [enumerate_lpa()],
#'     [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#'     The tables of the fitted candidates are stacked, with a first column
#'     `candidate` that names the candidate of each row. A candidate whose
#'     estimation failed has no rows;
#'   * the result of [diagnostics()]. The table is the one the result holds.
#' @return A data frame with one row per fit (one row per candidate for a
#'   grid) and the columns:
#'   * `n_starts`: the number of starts that were run to the end;
#'   * `n_screened`: the number of starts that were only screened in a
#'     short run; 0 when the fit screened none;
#'   * `n_converged`: the number of starts that converged;
#'   * `n_failed`: the number of starts that failed or gave no log
#'     likelihood;
#'   * `n_best_replicated`: the number of starts whose log likelihood is
#'     within `tolerance` of the best;
#'   * `share`: `n_best_replicated` divided by `n_starts`;
#'   * `tolerance`: the difference in log likelihood up to which two starts
#'     count as the same solution;
#'   * `best_start`: the number of the start that gave the fit;
#'   * `log_likelihood`: the log likelihood of the fit.
#'
#'   A count the fit does not record is `NA`.
#'
#'   Errors: `latents_not_a_fit` when `x` is none of these, for example a
#'   comparison, pooled imputations or a data frame.
#' @seealso [plot_trust()] draws the replication of every candidate of a
#'   grid; [get_checks()]; [sensitivity()] refits a model from other seeds.
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_replication(fit)
#'
#' # One row per candidate of a grid
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_replication(grid)
#' @export
get_replication <- function(x) {
  .latents_wild_table(x, "replication", "get_replication", substitute(x))
}

#' Print and plot everything about any object of the package: its summary, its
#' checks and every plot
#'
#' @description
#' `overview()` prints and plots everything about any object the package
#' returns: the summary, the diagnostic checks where the object has them, and
#' every plot the object can draw. What is shown depends on the object:
#'
#' * a profile, class or transition fit of [lpa()], [lca()], [multilpa()],
#'   [multilca()] or [lta()] (class `multilpa_transitions`): the output of
#'   [report()], that is, the summary, the descriptives, the checks and
#'   every plot view;
#' * a fit of another family (a mixture regression, a growth mixture, a
#'   general transition fit, an additive or a cross-level fit): its
#'   summary, the checks of [diagnostics()], and the plots of [plot_all()];
#' * pooled imputations: the summary and the plots of [plot_all()];
#' * an enumeration grid: its summary and the plots of [plot_all()]
#'   (criteria, entropy, AvePP, criterion changes, heat map, smallest
#'   class, trust and tree);
#' * a comparison, a bootstrap test, a bootstrap table of
#'   [parameter_inference()] or a sensitivity analysis: the object or its
#'   summary printed, and its plot;
#' * a [diagnostics()] result: the result printed, and the plots of the fit
#'   it checked.
#'
#' @param x Any object returned by the package: a fitted model, pooled
#'   imputations, an enumeration grid, a comparison, a bootstrap result, a
#'   sensitivity analysis or a diagnostics result.
#' @param plots Whether to draw the plots:
#'   * `TRUE` (default): prints, then draws every plot. The plots need the
#'     ggplot2 package; without it only the printed part is shown;
#'   * `FALSE`: prints only.
#' @return `x`, invisibly. The function is called for what it prints and
#'   draws.
#'
#'   Errors: `latents_not_a_fit` when `x` is not a fit, grid or other
#'   result of the package, for example a data frame; `latents_bad_argument`
#'   when `plots` is not `TRUE` or `FALSE`.
#' @seealso [report()] for the fits it covers; [plot_all()] for the plots
#'   alone; [get_checks()] for the checks as a table.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 2, seed = 1)
#' overview(fit, plots = FALSE)
#' @export
overview <- function(x, plots = TRUE) {
  .latents_require("`plots` must be TRUE or FALSE" = isTRUE(plots) || isFALSE(plots))
  if (.multilpa_any_fit(x)) {
    report(x, plots = plots)
    return(invisible(x))
  }
  kind <- if (inherits(x, "multilpa_diagnostics")) "diagnostics" else .latents_wild_kind(x)
  if (is.na(kind)) .latents_wild_refuse_object(x, "overview", substitute(x))
  if (identical(kind, "diagnostics")) {
    print(x)
    if (isTRUE(plots) && .gg_available()) print(plot_all(x$fit))
    return(invisible(x))
  }
  .latents_wild_print_summary(x)
  if (inherits(x, .latents_model_classes)) {
    cat("\n")
    print(diagnostics(x, plots = FALSE))
  }
  if (isTRUE(plots) && .gg_available()) print(plot_all(x))
  invisible(x)
}

#' Print an object's summary, or the object when its summary has no printer
#' @param x Any object of a known kind.
#' @return `NULL`, invisibly.
#' @noRd
.latents_wild_print_summary <- function(x) {
  if (is.data.frame(x)) {
    print(x)
    return(invisible(NULL))
  }
  summarised <- summary(x)
  # A summary returned as a bare list would print as a dump of its parts;
  # the object's own print method reads it instead.
  if (identical(class(summarised), "list")) print(x) else print(summarised)
  invisible(NULL)
}

#' Read one table of any grid, refusing a grid without it
#' @param x The object.
#' @param table The table's name in `get_results()`.
#' @param verb The calling verb.
#' @param expression The caller's `substitute(x)`.
#' @return A data frame.
#' @noRd
.latents_wild_grid_table <- function(x, table, verb, expression) {
  if (!inherits(x, .latents_enumeration_classes)) {
    label <- .latents_argument_label(expression, "x")
    stop(errorCondition(sprintf(
      "%s() takes an enumeration grid, but `%s` is %s.", verb, label,
      .latents_describe(x)),
      class = c("latents_not_a_fit", "latents_bad_argument"), call = NULL))
  }
  getters <- c(candidates = "get_candidates", fit = "get_candidates",
               best = "get_best", deltas = "get_deltas", knees = "get_knees",
               tree = "get_tree")
  # The grid's catalogue of tables, read without building them: get_results(x,
  # "all") builds every table, and a tree that cannot be built (one usable
  # candidate) would refuse the candidates table with it.
  available <- if (inherits(x, "multilpa_enumeration")) {
    names(.multilpa_catalogue(x))
  } else if (inherits(x, "latents_regression_enumeration")) {
    c("fit", "tree")
  } else {
    method <- utils::getS3method("get_results", intersect(
      class(x), c("latents_family_enumeration", "latents_transition_enumeration"))[1L])
    setdiff(eval(formals(method)$what), "all")
  }
  # A regression grid calls its candidates table "fit".
  if (identical(table, "candidates") && !"candidates" %in% available) table <- "fit"
  if (!table %in% available) {
    has <- unique(getters[intersect(names(getters), available)])
    stop(errorCondition(sprintf(
      "%s() has no table for %s. For this grid: %s.", verb, .latents_describe(x),
      paste(sprintf("%s()", has), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  get_results(x, table)
}

#' Get the table of every candidate model in an enumeration grid, with its
#' criteria and diagnostics
#'
#' `get_candidates()` returns the main table of an enumeration grid: one row per
#' candidate model, with its numbers of classes, its log likelihood, its
#' information criteria, its classification quality and whether its estimation
#' can be relied on.
#'
#' @details
#' The table is `get_results(x, "candidates")` without the columns that are
#' empty in every row. On a class grid, for example, the covariance
#' structure `model` and the checks of Gaussian indicators are left out,
#' and the columns `warnings` and `error` appear only when a candidate has
#' a warning or failed.
#'
#' @param x An enumeration grid of any kind, from [enumerate_lpa()],
#'   [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#' @return A data frame with one row per candidate. The columns depend on
#'   the kind of grid:
#'   * a profile or class grid: `candidate` (its name), `n_profiles`,
#'     `n_group_classes`, `model` (the covariance structure),
#'     `log_likelihood`, `n_parameters`, the criteria (`aic`, `kic`, and
#'     `bic`, `sabic`, `caic`, `awe`, `icl` and `clc` with the suffix
#'     `_groups` or `_individual` for the sample size they use),
#'     `profile_entropy` and `group_entropy` (relative entropy at each
#'     level), `avepp_min`, `avepp_max`, `avepp_mean` and `avepp_below`
#'     (the average posterior probabilities of the classes and how many
#'     fall below the cutoff 0.8; with the prefix `group_` for the group
#'     classes), `converged`, `boundary`, `min_standardized_eigenvalue`
#'     and `min_gaussian_mass` (checks of Gaussian indicators) and
#'     `n_best_replicated` (the number of starts that reached the best log
#'     likelihood);
#'   * a transition grid: `candidate`, `n_profiles`, `n_group_classes`,
#'     `log_likelihood`, `n_parameters`, `aic`, `bic`, `bic_individual`,
#'     `entropy`, `converged`, `boundary` and `n_best_replicated`;
#'   * a grid of group-class families: `candidate`, `model`, `family`,
#'     `between_variance`, `n_group_classes`, `log_likelihood`,
#'     `n_parameters`, `aic`, `bic`, `bic_individual`, `entropy`,
#'     `smallest_class`, `min_effective_groups`, `converged`, `boundary`
#'     and `n_best_replicated`;
#'   * a regression grid: `n_classes`, `n_group_classes`, `log_likelihood`,
#'     `n_parameters`, `aic`, `bic`, `bic_rows`, `sabic`, `icl`, `entropy`,
#'     `smallest_share`, `converged`, `n_best_replicated` and `best_bic`
#'     (`TRUE` for the candidate with the lowest BIC). It has no
#'     `candidate` or `boundary` column.
#'
#'   Errors: `latents_not_a_fit` when `x` is not an enumeration grid.
#' @seealso [get_best()], [get_deltas()] and [get_knees()] for what the
#'   criteria say; [plot_criteria()] draws them; [get_results()] for the
#'   full table.
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_candidates(grid)
#' @export
get_candidates <- function(x) {
  table <- .latents_wild_grid_table(x, "candidates", "get_candidates", substitute(x))
  # A column empty in every row says nothing about these candidates: the
  # covariance structure and Gaussian checks of a class grid, the error of a
  # grid where every candidate fitted. get_results() keeps them for 1.1.9 code.
  empty <- vapply(table, \(column) all(is.na(column) |
    (is.character(column) & !nzchar(trimws(column)))), logical(1))
  table[, !empty, drop = FALSE]
}

#' Get the converged candidate with the lowest value of each information
#' criterion
#'
#' `get_best()` returns, for each information criterion of an enumeration grid,
#' the converged candidate with the lowest value. Several criteria may name the
#' same candidate, and no criterion chooses a model on its own.
#'
#' @param x An enumeration grid from [enumerate_lpa()], [enumerate_lca()]
#'   or [enumerate_classes()]. A regression grid from
#'   [enumerate_regressions()] has no such table and is refused; its
#'   [get_candidates()] table marks the lowest BIC in `best_bic`.
#' @return A data frame.
#'
#'   On a profile or class grid it has one row per criterion and sample
#'   size, with the columns:
#'   * `criterion`: `aic`, `kic`, `bic`, `sabic`, `caic`, `awe`, `icl` or
#'     `clc`;
#'   * `convention`: `"groups"` or `"individuals"`, the sample size the
#'     criterion uses; `NA` for AIC and KIC, which use none;
#'   * `candidate`: the name of the converged candidate with the lowest
#'     value;
#'   * `n_profiles`, `n_group_classes`, `model`: that candidate's numbers
#'     of classes and covariance structure;
#'   * `value`: its value of the criterion.
#'
#'   On a transition grid or a grid of group-class families it is one row:
#'   the row of `get_results(x, "candidates")` with the lowest BIC.
#'
#'   Errors: `latents_bad_argument` for a regression grid, naming the
#'   tables it has; `latents_not_a_fit` when `x` is not an enumeration
#'   grid.
#' @seealso [get_candidates()] for every candidate and criterion;
#'   [get_deltas()], [get_knees()]; [plot_criteria()].
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_best(grid)
#' @export
get_best <- function(x) {
  .latents_wild_grid_table(x, "best", "get_best", substitute(x))
}

#' Get how much each information criterion improves as a class or a group class
#' is added
#'
#' `get_deltas()` returns, for every candidate of a profile or class grid, how
#' much lower its information criterion is than that of the candidate with one
#' profile fewer, and than that of the candidate with one group class fewer. The
#' criterion is `bic_groups`, BIC with the number of groups as the sample size.
#'
#' @param x A profile or class grid from [enumerate_lpa()],
#'   [enumerate_lca()] or [enumerate_classes()]. Transition grids, grids of
#'   group-class families and regression grids have no such table and are
#'   refused.
#' @return A data frame with one row per candidate and the columns:
#'   * `candidate`, `n_profiles`, `n_group_classes`, `model`: the candidate,
#'     its numbers of classes and its covariance structure;
#'   * `criterion`: the name of the criterion, `bic_groups`;
#'   * `value`: the candidate's value of the criterion;
#'   * `delta_profile`: the criterion of the candidate with one profile
#'     fewer (and the same group classes and structure) minus `value`. A
#'     positive number means that this candidate is lower;
#'   * `from_profile`: the name of that smaller candidate;
#'   * `delta_group_class`, `from_group_class`: the same for the candidate
#'     with one group class fewer.
#'
#'   A delta is `NA` where there is no such smaller candidate, or where
#'   either of the two did not converge.
#'
#'   Errors: `latents_bad_argument` for a transition, family or regression
#'   grid, naming the tables it has; `latents_not_a_fit` when `x` is not an
#'   enumeration grid.
#' @seealso [plot_deltas()] draws the changes; [get_knees()], [get_best()].
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_deltas(grid)
#' @export
get_deltas <- function(x) {
  .latents_wild_grid_table(x, "deltas", "get_deltas", substitute(x))
}

#' Get the candidate after which adding classes improves each information
#' criterion much less
#'
#' `get_knees()` returns the knee of every line of a profile or class grid: the
#' candidate after which adding classes lowers the criterion much less than
#' before, and whether the knee stays the same when the largest candidate is
#' left out. The criterion is `bic_groups`, BIC with the number of groups as the
#' sample size.
#'
#' @details
#' Both axes of a line are rescaled to run from 0 to 1. The knee is the
#' candidate that lies furthest below the straight line from the first
#' candidate to the last. A line needs three candidates to have a knee and
#' four for its stability to be assessed, and a line whose criterion does
#' not fall from its first to its last candidate has none.
#'
#' @param x A profile or class grid from [enumerate_lpa()],
#'   [enumerate_lca()] or [enumerate_classes()]. Transition grids, grids of
#'   group-class families and regression grids have no such table and are
#'   refused.
#' @return A data frame with one row per line and the columns:
#'   * `direction`: `"profiles"` (along the number of profiles, at a fixed
#'     number of group classes), `"group_classes"` (along the number of
#'     group classes, at a fixed number of profiles) or `"parameters"`
#'     (every candidate, ordered by its number of parameters);
#'   * `line`: what is held fixed along the line;
#'   * `criterion`: the name of the criterion, `bic_groups`;
#'   * `from`, `to`: the first and the last candidate of the line;
#'   * `n_points`: the number of candidates on the line;
#'   * `knee_candidate`: the candidate at the knee; `NA` when the line has
#'     none;
#'   * `distance`: how far the knee lies below the straight line, in the
#'     rescaled units;
#'   * `knee_without_largest`: the knee when the largest candidate is left
#'     out;
#'   * `stable`: `TRUE` when the two knees are the same candidate; `NA`
#'     when the line has fewer than four candidates;
#'   * `note`: why a knee or its stability is missing.
#'
#'   Errors: `latents_bad_argument` for a transition, family or regression
#'   grid, naming the tables it has; `latents_not_a_fit` when `x` is not an
#'   enumeration grid.
#' @seealso [get_deltas()] for the changes behind the lines; [get_best()];
#'   [plot_deltas()].
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_knees(grid)
#' @export
get_knees <- function(x) {
  .latents_wild_grid_table(x, "knees", "get_knees", substitute(x))
}

#' Get how the classes of each candidate split into the classes of the next, as
#' a table
#'
#' `get_tree()` returns the flows that [plot_tree()] plots, as a table: for each
#' pair of consecutive candidates of a tree, how much of each class of the
#' smaller candidate goes to each class of the larger one.
#'
#' @param x An enumeration grid of any kind, from [enumerate_lpa()],
#'   [enumerate_lca()], [enumerate_classes()] or [enumerate_regressions()].
#' @return A data frame with one row per flow, that is, per pair of a class
#'   of one candidate and a class of the next, and the columns:
#'   * `tree`: the name of the tree. A grid can hold several trees, for
#'     example one per covariance structure;
#'   * `level`: `"individuals"` when the cases are divided into profiles or
#'     classes, `"groups"` when the groups are divided into group classes;
#'   * `from_candidate`, `from_class`, `from_share`: the smaller candidate,
#'     one of its classes and that class's share of the units;
#'   * `to_candidate`, `to_class`, `to_share`: the same for the next
#'     candidate;
#'   * `mass`: the share of all units that goes from `from_class` to
#'     `to_class`;
#'   * `share_of_from`: `mass` as a share of `from_class`;
#'   * `share_of_to`: `mass` as a share of `to_class`.
#'
#'   Errors: `latents_no_tree` when the grid has fewer than two candidates
#'   that can form a tree; `latents_not_a_fit` when `x` is not an
#'   enumeration grid.
#' @seealso [plot_tree()] draws the flows and chooses the trees.
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_tree(grid)
#' @export
get_tree <- function(x) {
  .latents_wild_grid_table(x, "tree", "get_tree", substitute(x))
}
