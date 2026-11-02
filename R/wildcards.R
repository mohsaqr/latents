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
#' @return A ggplot object, or a `latents_plots` list of the individual and
#'   group levels for a classification view of a fit with group classes.
#' @noRd
.latents_wild_plot <- function(x, verb, expression, main, subtitle) {
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
  plot(x, what = view, main = main, subtitle = subtitle)
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

#' Plots by name, for any object
#'
#' @description
#' Each function draws one picture and takes any object the package returns:
#' a fit of any family, an enumeration grid, a model comparison, a bootstrap
#' test, a sensitivity analysis or a diagnostics result. The function's name
#' says what is drawn; there is no `what` to choose. What a function draws
#' depends on the object:
#'
#' | Function | Fits | Grids |
#' |---|---|---|
#' | `plot_sizes()` | effective class sizes | share of the smallest class |
#' | `plot_entropy()` | each unit's relative entropy, by class | relative entropy of every candidate |
#' | `plot_posteriors()` | posterior of each unit's assigned class | |
#' | `plot_avepp()` | average posterior matrix | lowest average posterior of every candidate |
#' | `plot_heatmap()` | profile means, standardized (profile and class fits) | criterion over class counts |
#' | `plot_bars()`, `plot_raincloud()`, `plot_parallel()`, `plot_pairs()` | the profile views of [plot_profiles()] (profile fits) | |
#' | `plot_trajectories()` | class curves of the outcome (regression and growth fits) | |
#' | `plot_individuals()`, `plot_classification()` | persons and their classes (growth and trajectory fits) | |
#' | `plot_coefficients()` | class coefficients with intervals (regression and growth fits) | |
#' | `plot_criteria()` | (a comparison: its criteria) | every information criterion |
#' | `plot_deltas()`, `plot_trust()` | | change in the criterion; convergence and replication |
#' | `plot_all()` | every view the object can draw | every view |
#'
#' The four classification views (`plot_sizes()`, `plot_entropy()`,
#' `plot_posteriors()`, `plot_avepp()`) apply to a fit of every family. For
#' profile, class and homogeneous transition fits they are the views of
#' [plot_diagnostics()]. The other families (mixture regression, growth,
#' additive, cross-level, general transition) draw them from the eligible,
#' weighted posteriors that [diagnostics()] summarises, at each level they
#' classify, so the plot carries the numbers of [get_avepp()],
#' [get_classification()] and [get_entropy()]; units of weight zero are left
#' out. A level with one class has sizes and a one-by-one AvePP matrix, but
#' no certainty or entropy to draw.
#'
#' A diagnostics result is drawn from the fit it was computed on. A
#' classification view of a fit that classifies at two levels (two or more
#' group classes) draws both: the result is then a `latents_plots` list with
#' elements `individuals` and `groups`.
#'
#' The existing [plot_profiles()], [plot_responses()], [plot_sequences()],
#' [plot_prevalence()], [plot_transitions()] and [plot_tree()] also take any
#' fit (the tree, any grid) and complete the set.
#'
#' @param x Any object returned by the package (see the table).
#' @param main,subtitle Optional title and subtitle. `NULL` keeps the default.
#' @return A ggplot object; for the classification views of a fit with group
#'   classes, and for `plot_all()`, a `latents_plots` list of ggplot objects
#'   named by view, which draws every plot when printed.
#'
#'   A function with nothing to draw for the object raises
#'   `latents_nothing_to_plot` (also `latents_bad_argument`), naming the
#'   functions that do apply. An object the package did not return raises
#'   `latents_not_a_fit`. A view whose ingredients the fit lacks raises the
#'   condition of the underlying method, for example
#'   `latents_no_continuous` for `plot_bars()` of a categorical fit.
#'   `plot_all()` leaves such views out and names them in a message.
#' @seealso [overview()] prints and draws everything at once; [get_checks()]
#'   and its siblings return the diagnostic tables.
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_avepp(fit)
#' plot_entropy(fit)
#' plot_bars(fit)
#'
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"), n_profiles = 1:3,
#'                       model = "EEI", n_starts = 2, seed = 1)
#' plot_criteria(grid)
#' plot_avepp(grid)
#' @name wildcard-plots
NULL

#' @rdname wildcard-plots
#' @export
plot_sizes <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_sizes", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_entropy <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_entropy", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_posteriors <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_posteriors", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_avepp <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_avepp", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_heatmap <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_heatmap", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_bars <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_bars", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_raincloud <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_raincloud", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_parallel <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_parallel", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_pairs <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_pairs", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_trajectories <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_trajectories", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_individuals <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_individuals", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_classification <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_classification", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_coefficients <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_coefficients", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_criteria <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_criteria", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_deltas <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_deltas", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
#' @export
plot_trust <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_trust", substitute(x), main, subtitle)
}

#' @rdname wildcard-plots
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

#' Diagnostic tables by name, for any fit, grid or diagnostics result
#'
#' @description
#' Each function returns one diagnostic table as a data frame and takes a fit
#' of any family, an enumeration grid or a [diagnostics()] result. The name is
#' the table; there is no `what` to choose.
#'
#' * `get_checks()`: every check with its status, criterion and reference.
#' * `get_entropy()`: relative entropy, one row per classification level.
#' * `get_classification()`: modal and estimated class sizes, average
#'   posteriors and odds of correct classification, one row per level and
#'   class.
#' * `get_avepp()`: the average posterior matrix, one row per level, assigned
#'   class and posterior class.
#' * `get_residuals()`: bivariate residuals, with Benjamini-Hochberg adjusted
#'   p-values, for families that define them.
#' * `get_replication()`: how many starts reached the best log likelihood.
#'
#' For a fit, the table is the one [diagnostics()] computes (with its
#' defaults), read as `get_results(diagnostics(x), "<table>")`. For a grid,
#' the table of every candidate is stacked, with a first column `candidate`;
#' a candidate whose estimation failed has no row. For a diagnostics result,
#' its own table is returned.
#'
#' @param x A fitted model of any family, an enumeration grid, or the result
#'   of [diagnostics()].
#' @return A data frame with the columns documented under [diagnostics()] and
#'   [get_results()] for that table; for a grid, `candidate` comes first.
#'
#'   A family without the table (bivariate residuals of a regression, growth,
#'   additive or cross-level fit, or of a fit with a noise component) raises
#'   `latents_bad_argument`; a grid none of whose candidates has it does too.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [diagnostics()], the wild-card plots ([plot_avepp()] and others),
#'   and [overview()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_checks(fit)
#' get_entropy(fit)
#' get_avepp(fit)
#'
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"), n_profiles = 1:3,
#'                       model = "EEI", n_starts = 2, seed = 1)
#' get_entropy(grid)
#' @name wildcard-tables
NULL

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

#' @rdname wildcard-tables
#' @export
get_checks <- function(x) {
  .latents_wild_table(x, "checks", "get_checks", substitute(x))
}

#' @rdname wildcard-tables
#' @export
get_entropy <- function(x) {
  .latents_wild_table(x, "entropy", "get_entropy", substitute(x))
}

#' @rdname wildcard-tables
#' @export
get_classification <- function(x) {
  .latents_wild_table(x, "classification", "get_classification", substitute(x))
}

#' @rdname wildcard-tables
#' @export
get_avepp <- function(x) {
  .latents_wild_table(x, "average_posteriors", "get_avepp", substitute(x))
}

#' @rdname wildcard-tables
#' @export
get_residuals <- function(x) {
  .latents_wild_table(x, "residuals", "get_residuals", substitute(x))
}

#' @rdname wildcard-tables
#' @export
get_replication <- function(x) {
  .latents_wild_table(x, "replication", "get_replication", substitute(x))
}

#' Everything about any object, printed and drawn
#'
#' @description
#' `overview()` takes any object the package returns and shows all of it: the
#' summary, the diagnostic checks and every plot the object can draw.
#'
#' * A profile, class or transition fit: [report()] (summary, descriptives,
#'   checks, every plot view).
#' * Any other fit (mixture regression, growth mixture, general transition,
#'   additive, cross-level, pooled): its summary, its [diagnostics()] checks
#'   where defined, and [plot_all()].
#' * An enumeration grid: its summary and [plot_all()] (criteria, entropy,
#'   AvePP, smallest class, BIC changes, trust, tree).
#' * A comparison, bootstrap test, sensitivity analysis or diagnostics result:
#'   the object printed and its plot.
#'
#' @param x Any object returned by the package.
#' @param plots `TRUE` (the default) draws the plots; `FALSE` prints only.
#' @return `x`, invisibly. Called for its printed and drawn output.
#'
#'   An object the package did not return raises `latents_not_a_fit`.
#' @seealso [report()], [plot_all()], [get_checks()].
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

#' Grid tables by name, for any enumeration grid
#'
#' @description
#' Each function returns one table of an enumeration grid as a data frame; the
#' name is the table, with no `what` to choose.
#'
#' * `get_candidates()`: one row per candidate, every criterion and quality
#'   column. Columns empty in every row are left out (for a class grid, the
#'   covariance `model` and the Gaussian checks; `error` when every candidate
#'   fitted); `get_results(x, "candidates")` keeps them.
#' * `get_best()`: the candidate with the lowest value of each criterion.
#' * `get_deltas()`: the change in the criterion with each added class and
#'   group class.
#' * `get_knees()`: where each line of the criterion bends, and whether the
#'   bend holds when the largest candidate is left out.
#' * `get_tree()`: how the classes of each candidate split into those of the
#'   next.
#'
#' [get_classification()], [get_entropy()], [get_avepp()] and [get_checks()]
#' stack the tables of every candidate of a grid as well.
#'
#' @param x An enumeration grid from [enumerate_classes()], [enumerate_lpa()],
#'   [enumerate_lca()] or [enumerate_regressions()].
#' @return A data frame, the table [get_results()] returns under that name
#'   (`"fit"` for the candidates of a regression grid), with `get_candidates()`
#'   leaving out the columns that are empty in every row.
#'
#'   A grid without the table (deltas, knees or tree of a family or transition
#'   grid) raises `latents_bad_argument`, naming the functions it has; any
#'   other object raises `latents_not_a_fit`.
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"), n_profiles = 1:3,
#'                       model = "EEI", n_starts = 2, seed = 1)
#' get_best(grid)
#' get_knees(grid)
#' @name wildcard-grid-tables
NULL

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

#' @rdname wildcard-grid-tables
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

#' @rdname wildcard-grid-tables
#' @export
get_best <- function(x) {
  .latents_wild_grid_table(x, "best", "get_best", substitute(x))
}

#' @rdname wildcard-grid-tables
#' @export
get_deltas <- function(x) {
  .latents_wild_grid_table(x, "deltas", "get_deltas", substitute(x))
}

#' @rdname wildcard-grid-tables
#' @export
get_knees <- function(x) {
  .latents_wild_grid_table(x, "knees", "get_knees", substitute(x))
}

#' @rdname wildcard-grid-tables
#' @export
get_tree <- function(x) {
  .latents_wild_grid_table(x, "tree", "get_tree", substitute(x))
}
