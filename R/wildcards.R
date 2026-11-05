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

#' Class sizes
#'
#' For a fit of any family, the effective number of units in each class
#' with its share, at each level the fit classifies (the cases, and the
#' groups when there are two or more group classes). For an enumeration
#' grid, the share of the smallest class of every candidate, with 5%
#' dashed.
#'
#' @param x A fitted model of any family, an enumeration grid or a
#'   [diagnostics()] result (drawn from its fit).
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object; a `latents_plots` list with elements
#'   `individuals` and `groups` for a fit that classifies at both levels.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [get_classification()] for the table; [plot_posteriors()],
#'   [plot_entropy()], [plot_avepp()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_sizes(fit)
#' @export
plot_sizes <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_sizes", substitute(x), main, subtitle)
}

#' Entropy of the classification
#'
#' For a fit of any family, each unit's entropy relative to a flat
#' posterior, by assigned class, at each level the fit classifies. For an
#' enumeration grid, the relative entropy of every candidate. A level with
#' one class has no entropy to draw.
#'
#' @param x A fitted model of any family, an enumeration grid or a
#'   [diagnostics()] result (drawn from its fit).
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object; a `latents_plots` list with elements
#'   `individuals` and `groups` for a fit that classifies at both levels.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [get_entropy()] for the table; [plot_posteriors()], [plot_avepp()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_entropy(fit)
#' @export
plot_entropy <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_entropy", substitute(x), main, subtitle)
}

#' Certainty of each unit's classification
#'
#' The posterior probability of each unit's assigned class, by class, at
#' each level a fit of any family classifies. A level with one class has
#' nothing to draw.
#'
#' @param x A fitted model of any family or a [diagnostics()] result
#'   (drawn from its fit).
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object; a `latents_plots` list with elements
#'   `individuals` and `groups` for a fit that classifies at both levels.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [get_classification()] for the table; [plot_entropy()], [plot_avepp()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_posteriors(fit)
#' @export
plot_posteriors <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_posteriors", substitute(x), main, subtitle)
}

#' Average posterior probabilities
#'
#' For a fit of any family, the average posterior probability matrix
#' (assigned class by posterior class), at each level the fit classifies.
#' For an enumeration grid, the lowest average posterior of every
#' candidate, with the 0.8 cutoff dashed.
#'
#' @param x A fitted model of any family, an enumeration grid or a
#'   [diagnostics()] result (drawn from its fit).
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object; a `latents_plots` list with elements
#'   `individuals` and `groups` for a fit that classifies at both levels.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [get_avepp()] for the table; [plot_entropy()], [plot_posteriors()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_avepp(fit)
#' @export
plot_avepp <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_avepp", substitute(x), main, subtitle)
}

#' Profile means, or a criterion, as a heat map
#'
#' For a profile or class fit, the profile means in observed standard
#' deviations from the observed mean, or the response probabilities of
#' categorical indicators. For an enumeration grid, the information
#' criterion over profiles and group classes, relative to its lowest
#' value (BIC over groups on a profile or class grid, BIC on the others).
#'
#' @param x A profile or class fit, or an enumeration grid.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_profiles()], [plot_bars()]; [plot_enumeration()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_heatmap(fit)
#' @export
plot_heatmap <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_heatmap", substitute(x), main, subtitle)
}

#' Profile means as bars
#'
#' The profile means of a profile or class fit as grouped bars from zero,
#' with 95% intervals when the fit has standard errors.
#'
#' @param x A profile or class fit with continuous indicators.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_profiles()], [plot_heatmap()], [plot_raincloud()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_bars(fit)
#' @export
plot_bars <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_bars", substitute(x), main, subtitle)
}

#' Distribution of each indicator by profile
#'
#' Each continuous indicator's distribution by assigned profile: density,
#' box and the observations.
#'
#' @param x A profile or class fit with continuous indicators.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_parallel()], [plot_pairs()], [plot_profiles()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_raincloud(fit)
#' @export
plot_raincloud <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_raincloud", substitute(x), main, subtitle)
}

#' Every case as a line across the indicators
#'
#' Each case's values on the continuous indicators as one line, one panel
#' per assigned profile.
#'
#' @param x A profile or class fit with continuous indicators.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_raincloud()], [plot_pairs()], [plot_profiles()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_parallel(fit)
#' @export
plot_parallel <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_parallel", substitute(x), main, subtitle)
}

#' Scatter-plot matrix of the indicators by profile
#'
#' The continuous indicators against each other, cases coloured by
#' assigned profile, with each profile's 95% covariance ellipse.
#'
#' @param x A profile or class fit with continuous indicators.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_raincloud()], [plot_parallel()], [plot_profiles()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' plot_pairs(fit)
#' @export
plot_pairs <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_pairs", substitute(x), main, subtitle)
}

#' Class trajectories of the outcome
#'
#' For a growth mixture or a regression that classifies persons, each
#' class's mean trajectory over time with its band, the observed
#' trajectories behind. For another mixture regression, the outcome
#' against the first numeric predictor, one line per class. The fit's own
#' `plot(x, what = "trajectories")` also takes `time`, `facet`,
#' `max_persons` and `spread`, and `plot(x, what = "fitted")` takes
#' `predictor` and `level`.
#'
#' @param x A growth mixture, or a mixture regression from
#'   [mixture_regression()].
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_individuals()], [plot_coefficients()], [plot_classification()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_trajectories(fit)
#' @export
plot_trajectories <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_trajectories", substitute(x), main, subtitle)
}

#' Trajectories of individual persons
#'
#' The observed and predicted trajectory of each of twelve persons of a
#' growth mixture or a regression that classifies persons, coloured by
#' class. The fit's own `plot(x, what = "individuals")` takes `persons`,
#' a number or identifiers.
#'
#' @param x A growth mixture, or a mixture regression that classifies
#'   persons (`class_level = "group"`).
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_trajectories()], [plot_classification()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_individuals(fit)
#' @export
plot_individuals <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_individuals", substitute(x), main, subtitle)
}

#' Class probabilities of every person
#'
#' One stacked bar per person of a growth mixture or a regression that
#' classifies persons: the probability of each class, persons ordered by
#' assigned class and certainty.
#'
#' @param x A growth mixture, or a mixture regression that classifies
#'   persons (`class_level = "group"`).
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_posteriors()] for the certainty of the assigned class alone;
#'   [plot_individuals()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_classification(fit)
#' @export
plot_classification <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_classification", substitute(x), main, subtitle)
}

#' Class coefficients with intervals
#'
#' Each class's regression coefficients of a mixture regression or growth
#' mixture with their intervals. The fit's own
#' `plot(x, what = "coefficients")` takes `level` (default 0.95).
#'
#' @param x A growth mixture, or a mixture regression from
#'   [mixture_regression()].
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_trajectories()]; [get_results()] for the coefficient table.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' fit <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
#'                           n_classes = 2, id = "student", class_level = "group",
#'                           random = "wave", n_starts = 2, seed = 1)
#' plot_coefficients(fit)
#' @export
plot_coefficients <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_coefficients", substitute(x), main, subtitle)
}

#' Every information criterion
#'
#' For an enumeration grid, every information criterion against the
#' number of classes, as panels of one plot with each criterion's lowest
#' value circled. For a comparison from [compare_models()], the criteria
#' of each compared model.
#'
#' @param x An enumeration grid, or a comparison from [compare_models()].
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [plot_enumeration()] to choose the criteria; [get_candidates()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_criteria(grid)
#' @export
plot_criteria <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_criteria", substitute(x), main, subtitle)
}

#' Change in the criterion with each added class
#'
#' The decrease in the information criterion of an enumeration grid with
#' each added profile and each added group class (BIC over groups on a
#' profile or class grid, BIC on the others).
#'
#' @param x An enumeration grid.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [get_deltas()] and [get_knees()] for the tables; [plot_selection()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_deltas(grid)
#' @export
plot_deltas <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_deltas", substitute(x), main, subtitle)
}

#' Convergence, boundary and replication of every candidate
#'
#' Whether each candidate of an enumeration grid converged, reached a
#' boundary and replicated its best log likelihood across starts.
#'
#' @param x An enumeration grid.
#' @param main,subtitle Title and subtitle; `NULL` (default) keeps the
#'   view's own.
#' @return A ggplot object.
#'   An object with nothing to draw raises `latents_nothing_to_plot`,
#'   naming the plot functions that apply; an object the package did not
#'   return raises `latents_not_a_fit`.
#' @seealso [get_replication()] for the table; [plot_selection()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' plot_trust(grid)
#' @export
plot_trust <- function(x, main = NULL, subtitle = NULL) {
  .latents_wild_plot(x, "plot_trust", substitute(x), main, subtitle)
}

#' Every plot of an object
#'
#' Every view the object can draw, as a list that draws each plot when
#' printed: for a fit, its views and the classification views at each level;
#' for an enumeration grid, the criteria, entropy, AvePP, smallest class,
#' criterion changes, trust and tree; for a comparison, bootstrap test or
#' sensitivity analysis, its one plot; for a diagnostics result, the views
#' of the fit it checked. A view the
#' object lacks the ingredients for is left out and named in a message.
#'
#' @param x Any object returned by the package.
#' @return A `latents_plots` list of ggplot objects named by view. An object
#'   the package did not return raises `latents_not_a_fit`.
#' @seealso [overview()] prints the summary and checks as well.
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

#' Every diagnostic check of a fit
#'
#' Each check [diagnostics()] runs, with its status, criterion and
#' reference: estimation, class sizes, classification, local independence
#' and information. For a fit, the table is the one [diagnostics()] computes with its
#' defaults; for a grid, the tables of every candidate are stacked with
#' a first column `candidate` (a candidate whose estimation failed has no
#' row); for a diagnostics result, its own table.
#'
#' @param x A fitted model of any family, an enumeration grid, or the
#'   result of [diagnostics()].
#' @return A data frame, one row per check (and candidate): `domain, check, level, value, display, criterion, status, detail, reference`.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [diagnostics()] for the checks and their criteria; [get_entropy()],
#'   [get_classification()], [get_avepp()], [get_replication()],
#'   [get_residuals()] for the tables behind them; [overview()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_checks(fit)
#' @export
get_checks <- function(x) {
  .latents_wild_table(x, "checks", "get_checks", substitute(x))
}

#' Relative entropy of the classification
#'
#' The entropy of the posterior class probabilities relative to a flat
#' posterior, at each level the fit classifies. For a fit, the table is the one [diagnostics()] computes with its
#' defaults; for a grid, the tables of every candidate are stacked with
#' a first column `candidate` (a candidate whose estimation failed has no
#' row); for a diagnostics result, its own table.
#'
#' @param x A fitted model of any family, an enumeration grid, or the
#'   result of [diagnostics()].
#' @return A data frame, one row per classification level (and candidate): `level, n_classes, n_units, entropy_sum, relative_entropy`.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [plot_entropy()] draws it; [get_classification()], [get_avepp()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_entropy(fit)
#' @export
get_entropy <- function(x) {
  .latents_wild_table(x, "entropy", "get_entropy", substitute(x))
}

#' Class sizes and classification certainty
#'
#' Modal and estimated class sizes, average posterior probability and
#' odds of correct classification, at each level the fit classifies. For a fit, the table is the one [diagnostics()] computes with its
#' defaults; for a grid, the tables of every candidate are stacked with
#' a first column `candidate` (a candidate whose estimation failed has no
#' row); for a diagnostics result, its own table.
#'
#' @param x A fitted model of any family, an enumeration grid, or the
#'   result of [diagnostics()].
#' @return A data frame, one row per level and class (and candidate): `level, class, n_modal, proportion_modal, estimated_n, estimated_proportion, average_posterior, odds_correct_classification`.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [plot_sizes()] and [plot_posteriors()] draw it; [get_avepp()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_classification(fit)
#' @export
get_classification <- function(x) {
  .latents_wild_table(x, "classification", "get_classification", substitute(x))
}

#' Average posterior probability matrix
#'
#' The average posterior probability of each class among the units
#' assigned to each class, at each level the fit classifies. For a fit, the table is the one [diagnostics()] computes with its
#' defaults; for a grid, the tables of every candidate are stacked with
#' a first column `candidate` (a candidate whose estimation failed has no
#' row); for a diagnostics result, its own table.
#'
#' @param x A fitted model of any family, an enumeration grid, or the
#'   result of [diagnostics()].
#' @return A data frame, one row per level, assigned class and class (and candidate): `level, assigned_class, class, n_assigned, average_posterior`.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [plot_avepp()] draws it; [get_classification()], [get_entropy()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_avepp(fit)
#' @export
get_avepp <- function(x) {
  .latents_wild_table(x, "average_posteriors", "get_avepp", substitute(x))
}

#' Bivariate residuals of the indicators
#'
#' Each pair of indicators within each class, checked for association the
#' model leaves unexplained: a residual correlation with a normal test for
#' two Gaussian indicators, a chi-square of observed against expected
#' counts for two categorical ones, and a descriptive correlation without a
#' test for a mixed pair or one with a count indicator, with
#' Benjamini-Hochberg adjusted p-values, for the families that define
#' them. For a fit, the table is the one [diagnostics()] computes with its
#' defaults; for a grid, the tables of every candidate are stacked with
#' a first column `candidate` (a candidate whose estimation failed has no
#' row); for a diagnostics result, its own table.
#'
#' @param x A fitted model of any family, an enumeration grid, or the
#'   result of [diagnostics()].
#' @return A data frame, one row per class and indicator pair (and candidate): `profile, indicator_1, indicator_2, kind, observed, expected, residual, effective_n, statistic, df, p_value, p_adjusted`.
#'   A family without bivariate residuals (regression, growth, additive
#'   and cross-level fits, a general transition fit with occasion-dependent
#'   structure, or a fit with a noise component) raises
#'   `latents_bad_argument`; a grid none of whose candidates has them
#'   does too.  Any other object raises `latents_not_a_fit`.
#' @seealso [diagnostics()] for the test; [get_checks()].
#' @examples
#' fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
#'            n_classes = 2, n_starts = 2, seed = 1)
#' get_residuals(fit)
#' @export
get_residuals <- function(x) {
  .latents_wild_table(x, "residuals", "get_residuals", substitute(x))
}

#' How many starts reached the best log likelihood
#'
#' The number of starts, how many were screened, converged and failed,
#' and how many replicated the best log likelihood within tolerance. For a fit, the table is the one [diagnostics()] computes with its
#' defaults; for a grid, the tables of every candidate are stacked with
#' a first column `candidate` (a candidate whose estimation failed has no
#' row); for a diagnostics result, its own table.
#'
#' @param x A fitted model of any family, an enumeration grid, or the
#'   result of [diagnostics()].
#' @return A data frame, one row per fit (and candidate): `n_starts, n_screened, n_converged, n_failed, n_best_replicated, share, tolerance, best_start, log_likelihood`.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [plot_trust()] draws it for a grid; [get_checks()].
#' @examples
#' fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
#'            n_starts = 2, seed = 1)
#' get_replication(fit)
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

#' Every candidate of a grid
#'
#' One row per candidate with its class counts, log likelihood, every
#' information criterion and quality measure, and its convergence,
#' boundary and replication status. Columns empty in every row are left
#' out (for a class grid, the covariance `model` and the Gaussian checks;
#' `error` when every candidate fitted); `get_results(x, "candidates")`
#' keeps them.
#'
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()], [enumerate_lca()] or [enumerate_regressions()].
#' @return A data frame, one row per candidate. On a profile or class grid: candidate, n_profiles, n_group_classes, model, log_likelihood, n_parameters, the criteria (`aic`, `bic_groups`, ...), `profile_entropy`, the AvePP summaries, `converged`, `boundary`, `n_best_replicated`, and the Gaussian checks where they apply.
#'   A transition or family grid has one `entropy`, `bic` and `bic_individual`, and no AvePP summaries; a regression grid has `n_classes`, `bic`, `bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share` and `best_bic`, without `candidate` or `boundary`.
#'   Any other object raises `latents_not_a_fit`.
#' @seealso [get_best()], [get_deltas()], [plot_criteria()]; [get_results()] for
#'   the columns.
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

#' The best candidate under each criterion
#'
#' The candidate with the lowest (or, by the criterion's convention,
#' highest) value of each information criterion. Several criteria may
#' name the same candidate; none selects a model on its own.
#'
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()] or [enumerate_lca()].
#' @return A data frame. On a profile or class grid, one row per criterion and convention: `criterion, convention, candidate, n_profiles, n_group_classes, model, value`.
#'   On a transition or family grid, the one row of the candidates table with the lowest BIC.
#'   A regression grid raises `latents_bad_argument`, naming the tables it has. Any other object raises `latents_not_a_fit`.
#' @seealso [get_candidates()], [get_deltas()], [get_knees()].
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_best(grid)
#' @export
get_best <- function(x) {
  .latents_wild_grid_table(x, "best", "get_best", substitute(x))
}

#' Change in the criterion with each added class
#'
#' The decrease in the information criterion from the candidate with one
#' profile fewer, and from the one with one group class fewer (BIC over
#' groups).
#'
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()] or [enumerate_lca()].
#' @return A data frame, one row per candidate: `candidate, n_profiles, n_group_classes, model, criterion, value, delta_profile, from_profile, delta_group_class, from_group_class`.
#'   A regression, transition or family grid has no deltas and raises
#'   `latents_bad_argument`, naming the tables it has.  Any other object raises `latents_not_a_fit`.
#' @seealso [plot_deltas()] draws it; [get_knees()].
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_deltas(grid)
#' @export
get_deltas <- function(x) {
  .latents_wild_grid_table(x, "deltas", "get_deltas", substitute(x))
}

#' Where each criterion line bends
#'
#' The knee of each line of the criterion over the number of classes,
#' and whether it holds when the largest candidate is left out.
#'
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()] or [enumerate_lca()].
#' @return A data frame, one row per line and direction: `direction, line, criterion, from, to, n_points, knee_candidate, distance, knee_without_largest, stable, note`.
#'   A regression, transition or family grid has no knees and raises
#'   `latents_bad_argument`, naming the tables it has.  Any other object raises `latents_not_a_fit`.
#' @seealso [get_deltas()], [get_best()].
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_knees(grid)
#' @export
get_knees <- function(x) {
  .latents_wild_grid_table(x, "knees", "get_knees", substitute(x))
}

#' How the classes of each candidate split into the next
#'
#' The flows [plot_tree()] draws: for each pair of consecutive candidates
#' of a tree, how much of each class goes to each class of the next.
#'
#' @param x An enumeration grid from [enumerate_classes()],
#'   [enumerate_lpa()], [enumerate_lca()] or [enumerate_regressions()].
#' @return A data frame, one row per flow: `tree, level, from_candidate, from_class, from_share, to_candidate, to_class, to_share, mass, share_of_from, share_of_to`.
#'   A grid without a tree raises `latents_bad_argument`, naming the
#'   tables it has.  Any other object raises `latents_not_a_fit`.
#' @seealso [plot_tree()] draws it.
#' @examples
#' grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
#'                       n_profiles = 1:3, model = "EEI",
#'                       n_starts = 2, seed = 1)
#' get_tree(grid)
#' @export
get_tree <- function(x) {
  .latents_wild_grid_table(x, "tree", "get_tree", substitute(x))
}
