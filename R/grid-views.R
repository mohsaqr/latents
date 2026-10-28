# One contract for every enumeration grid: which views it can draw, and the
# table those views read.
#
# The package makes four kinds of grid. A profile or class grid
# (enumerate_classes(), enumerate_lpa(), enumerate_lca()) crosses numbers of
# profiles with numbers of group classes; a transition grid
# (enumerate_classes(time = )) does the same for lta() fits; a group-class
# family grid (enumerate_classes(family = )) crosses families with numbers of
# group classes and has no profiles; a regression grid (enumerate_regressions())
# crosses numbers of classes with numbers of group classes. Every grid verb
# (plot(), plot_enumeration(), plot_selection(), plot_tree(), plot_views())
# asks `.latents_grid_views()` what the grid offers, and the drawing code reads
# the grid through `.latents_grid_table()`, so a grid kind is described once.

#' The kind of an enumeration grid
#' @param x Any object.
#' @return `"profiles"`, `"transitions"`, `"families"` or `"regressions"`, or
#'   `NA_character_` for an object that is not a grid.
#' @noRd
.latents_grid_kind <- function(x) {
  kinds <- c(multilpa_enumeration = "profiles",
             latents_transition_enumeration = "transitions",
             latents_family_enumeration = "families",
             latents_regression_enumeration = "regressions")
  at <- inherits(x, names(kinds), which = TRUE)
  if (!any(at > 0L)) return(NA_character_)
  unname(kinds[which(at > 0L)[1L]])
}

#' The views that are not criteria, in the order they are listed
#' @return A character vector of `what` values.
#' @noRd
.latents_grid_view_names <- function() {
  c("entropy", "avepp", "diagnostics", .latents_selection_views(), "selection",
    "tree")
}

#' What a grid can draw
#'
#' The registry every grid verb consults. A view is offered when the grid's
#' table or fits carry what it draws: the criterion columns of the table; the
#' classification-quality views from the fits' posteriors (or the table's
#' entropy); the model-selection views, which lay the candidates out by number
#' of classes and number of group classes; and the tree, which divides each
#' named candidate's units by their posteriors. A view that is not offered has
#' a reason, which the refusal gives.
#'
#' @param x An enumeration grid.
#' @param verb The calling function's name, for the refusal of a non-grid.
#' @return A list: `kind`; `criteria`, the criterion columns, named by their
#'   display titles; `default`, the criteria drawn when `what` is not given;
#'   `using`, the criterion the selection views read by default; `axis`, the
#'   table column on the x-axis (`"n_profiles"` or `"n_group_classes"`);
#'   `word`, what the classes on that axis are called; `views`, the offered
#'   views other than criteria; `refused`, a named character vector giving the
#'   reason for each view not offered.
#' @noRd
.latents_grid_views <- function(x, verb = "This function") {
  kind <- .latents_grid_kind(x)
  if (is.na(kind)) .latents_refuse_non_grid(x, verb)
  no_profiles <- paste(
    "it lays the candidates out by number of profiles and number of group",
    "classes, and a family grid's candidates have group classes only",
    "(compare them in as.data.frame(x))")
  spec <- switch(
    kind,
    profiles = list(
      criteria = .multilpa_enumeration_criteria(),
      default = c("aic", "bic_groups", "bic_individual", "icl_individual"),
      using = "bic_groups", axis = "n_profiles", word = .latents_class_word(x),
      refused = character()),
    transitions = list(
      criteria = c(AIC = "aic", BIC = "bic", "BIC (individuals)" = "bic_individual"),
      default = c("aic", "bic", "bic_individual"), using = "bic",
      axis = "n_profiles", word = "profiles", refused = character()),
    families = list(
      criteria = c(AIC = "aic", BIC = "bic", "BIC (individuals)" = "bic_individual"),
      default = c("aic", "bic", "bic_individual"), using = "bic",
      axis = "n_group_classes", word = "group classes",
      refused = c(
        deltas = no_profiles, heatmap = no_profiles, smallest = no_profiles,
        trust = no_profiles, selection = no_profiles)),
    regressions = list(
      criteria = c(AIC = "aic", BIC = "bic", "BIC (rows)" = "bic_rows",
                   SABIC = "sabic", ICL = "icl"),
      default = c("aic", "bic", "sabic", "icl"), using = "bic",
      axis = "n_profiles", word = "classes", refused = character()))
  spec$kind <- kind
  spec$views <- setdiff(.latents_grid_view_names(), names(spec$refused))
  spec
}

#' The views of one grid, as `plot_views(x)` returns them
#'
#' @param x An enumeration grid.
#' @return A base `data.frame`, one row per `what` value: `type`, `group`
#'   (`"criteria"`, `"quality"`, `"selection"` or `"tree"`), `supported`,
#'   `description` and `reason` (why a view is not offered; `NA` when it is).
#' @noRd
.latents_grid_view_table <- function(x) {
  views <- .latents_grid_views(x, "plot_views")
  criteria <- unname(views$criteria)
  titles <- names(views$criteria) %||% toupper(criteria)
  if (identical(views$kind, "profiles")) {
    titles <- names(.multilpa_distinct_criteria(data.frame(), criteria))
  }
  named <- .latents_grid_view_names()
  described <- c(
    entropy = "Relative entropy of every candidate, at each level with two or more classes",
    avepp = "The lowest average posterior probability of every candidate, with the cutoff dashed",
    diagnostics = "Entropy and lowest AvePP together",
    deltas = "The decrease in the criterion in `using` with each added class and each added group class",
    heatmap = "The criterion in `using` over classes and group classes, relative to the lowest",
    smallest = "The share of the smallest class of every candidate, with 5% dashed",
    trust = "Convergence, boundary and replication of every candidate",
    selection = "Every selection view the grid offers, and diagnostics, as a list",
    tree = "How classes split as more are added, one plot per tree (see plot_tree())")
  data.frame(
    type = c(criteria, "all", named),
    group = c(rep("criteria", length(criteria) + 1L),
              ifelse(named %in% c("entropy", "avepp", "diagnostics"), "quality",
                     ifelse(named == "tree", "tree", "selection"))),
    supported = c(rep(TRUE, length(criteria) + 1L), named %in% views$views),
    description = c(sprintf("%s across the candidates; lower is preferred", titles),
                    "Every criterion of the grid, as panels of one plot",
                    unname(described[named])),
    reason = c(rep(NA_character_, length(criteria) + 1L),
               unname(views$refused[named])),
    stringsAsFactors = FALSE)
}

#' Refuse an object that is not an enumeration grid
#' @param x The object passed to a grid verb.
#' @param verb The calling function's name, for the message.
#' @return Does not return; raises `latents_bad_argument`.
#' @noRd
.latents_refuse_non_grid <- function(x, verb) {
  stop(errorCondition(sprintf(paste(
    "%s() takes an enumeration grid from enumerate_classes(), enumerate_lpa(),",
    "enumerate_lca() or enumerate_regressions(), not a `%s` object."),
    verb, class(x)[1L]), class = "latents_bad_argument", call = NULL))
}

#' Everything a grid's `what` takes, as text for a message
#' @param views The output of `.latents_grid_views()`.
#' @return One string.
#' @noRd
.latents_grid_offer_text <- function(views) {
  paste(sprintf("\"%s\"", c(unname(views$criteria), "all", views$views)),
        collapse = ", ")
}

#' Refuse a view a grid does not offer
#'
#' Criterion names are checked separately: an unknown name is an unknown
#' criterion. This refuses the named views (quality, selection, tree) that the
#' grid's kind cannot draw, with the reason and what the grid does offer.
#'
#' @param x An enumeration grid.
#' @param what The requested views.
#' @param verb The calling function's name, for the message.
#' @return `NULL`, invisibly, when every requested view is offered.
#' @noRd
.latents_require_grid_views <- function(x, what, verb) {
  views <- .latents_grid_views(x, verb)
  refused <- intersect(what, names(views$refused))
  if (length(refused) == 0L) return(invisible(NULL))
  stop(errorCondition(sprintf(
    "%s() cannot draw `what = \"%s\"` for a `%s` grid: %s. This grid offers %s.",
    verb, refused[1L], class(x)[1L], views$refused[[refused[1L]]],
    .latents_grid_offer_text(views)),
    class = c("latents_unsupported_view", "latents_bad_argument"), call = NULL))
}

#' A grid's candidates in the layout every grid view reads
#'
#' One row per candidate, with `candidate`, `n_group_classes`, `model`,
#' `log_likelihood`, `n_parameters`, the criterion columns, `converged` (never
#' `NA`), `boundary` and `n_best_replicated`; and `n_profiles`, the number of
#' profiles or classes, except on a family grid, whose candidates have none.
#' A regression grid's `n_classes` is its `n_profiles`, its candidates are
#' named `classes2` (with `_groups2` when the grid has group classes), and it
#' records no boundary (`NA`).
#'
#' @param x An enumeration grid.
#' @return A base `data.frame`.
#' @noRd
.latents_grid_table <- function(x) {
  kind <- .latents_grid_kind(x)
  if (!identical(kind, "regressions")) {
    table <- .latents_named_table(x)
    table$converged <- table$converged %in% TRUE
    return(table)
  }
  table <- .subset2(x, "table")
  single_level <- all(table$n_group_classes %in% 1L)
  candidate <- paste0("classes", table$n_classes,
                      if (single_level) "" else paste0("_groups", table$n_group_classes))
  data.frame(
    candidate = candidate, n_profiles = as.integer(table$n_classes),
    n_group_classes = as.integer(table$n_group_classes), model = NA_character_,
    table[setdiff(names(table), c("n_classes", "n_group_classes", "converged"))],
    converged = table$converged %in% TRUE, boundary = NA,
    stringsAsFactors = FALSE)
}

#' The posteriors of one regression-grid fit at one level
#' @param fit A cell of a regression grid: a fit, or the condition it raised.
#' @param level `"individuals"` (the units classified into classes) or
#'   `"groups"` (the groups classified into group classes).
#' @return A posterior matrix, or `NULL` when the cell has none at that level.
#' @noRd
.latents_regression_posteriors <- function(fit, level) {
  if (inherits(fit, "condition") || is.null(fit)) return(NULL)
  expectation <- fit$expectation
  if (inherits(fit, "latents_growth_mixture")) {
    if (identical(level, "groups")) return(expectation$structure$rho)
    return(expectation$posterior)
  }
  if (identical(level, "groups")) {
    return(if (identical(fit$spec$nesting, "two-level")) expectation$rho)
  }
  if (identical(fit$spec$nesting, "group")) expectation$group_tau else
    expectation$tau
}

#' A grid's table with the classification quality of every candidate
#'
#' Adds `profile_entropy` and `group_entropy` (relative entropy at each
#' level) and `avepp_min` and `group_avepp_min` (the lowest average posterior
#' probability at each level) to `.latents_grid_table()`; a profile or class
#' grid gives `.multilpa_enumeration_grid()`, which carries the full AvePP
#' summary. A level a candidate does not have, or has with one class, is `NA`.
#'
#' @param x An enumeration grid.
#' @param cutoff The AvePP cutoff.
#' @return A base `data.frame`, one row per candidate.
#' @noRd
.latents_grid_quality_table <- function(x, cutoff = 0.8) {
  kind <- .latents_grid_kind(x)
  if (identical(kind, "profiles")) return(.multilpa_enumeration_grid(x, cutoff))
  stopifnot("`cutoff` must be a single number between 0 and 1" =
              is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff) &&
              cutoff > 0 && cutoff <= 1)
  grid <- .latents_grid_table(x)
  fits <- unname(.subset2(x, "fits"))
  pieces_of <- function(level) {
    lapply(fits, function(fit) {
      absent <- list(posteriors = NULL, weights = NULL)
      if (is.null(fit) || inherits(fit, "condition")) return(absent)
      actual <- level
      if (identical(kind, "regressions") &&
          inherits(fit, "latents_mixture_regression") &&
          identical(fit$spec$nesting, "group")) {
        if (identical(level, "groups")) return(absent)
        actual <- "groups"
      }
      if (!actual %in% names(.latents_classifying_levels(fit))) return(absent)
      .latents_classification_inputs(fit, actual)
    })
  }
  individual_pieces <- pieces_of("individuals")
  group_pieces <- pieces_of("groups")
  posteriors_of <- function(level) lapply(if (identical(level, "groups"))
    group_pieces else individual_pieces, `[[`, "posteriors")
  level_weights <- function(level) lapply(if (identical(level, "groups"))
    group_pieces else individual_pieces, `[[`, "weights")
  entropy_of <- function(posteriors, weights) {
    unlist(Map(function(p, w) {
      if (is.null(p)) NA_real_ else .multilpa_relative_entropy(p, w)
    }, posteriors, weights), use.names = FALSE) %||% numeric()
  }
  lowest_avepp <- function(posteriors, weights) {
    unlist(Map(function(p, w) {
      .multilpa_avepp_summary(p, cutoff, w)[["min"]]
    }, posteriors, weights), use.names = FALSE) %||% numeric()
  }
  individuals <- posteriors_of("individuals")
  groups <- posteriors_of("groups")
  group_weights <- level_weights("groups")
  # The tables already carry the entropy their own fits report: a family
  # grid's is its group classes', the others' their individuals'.
  grid$profile_entropy <- switch(kind, families = NA_real_, grid$entropy)
  grid$group_entropy <- switch(kind, families = grid$entropy,
                               entropy_of(groups, group_weights))
  grid$avepp_min <- lowest_avepp(individuals, level_weights("individuals"))
  grid$group_avepp_min <- lowest_avepp(groups, group_weights)
  grid
}

#' The criteria to draw, named by their display titles
#'
#' A profile grid pairs each criterion's group and individual conventions and
#' keeps a pair once where a single-level grid makes them equal
#' (`.multilpa_distinct_criteria()`). The other grids name each column; a
#' column equal to one requested before it (a regression grid's `bic_rows` on
#' single-level data, where rows are the units) is drawn once.
#'
#' @param x An enumeration grid.
#' @param grid `.latents_grid_table(x)`.
#' @param criterion The requested criterion columns.
#' @return The columns to draw, named by title.
#' @noRd
.latents_grid_criteria_titles <- function(x, grid, criterion) {
  views <- .latents_grid_views(x)
  if (identical(views$kind, "profiles")) {
    return(.multilpa_distinct_criteria(grid, criterion))
  }
  repeats <- vapply(seq_along(criterion), function(index) {
    any(vapply(criterion[seq_len(index - 1L)], function(earlier) {
      isTRUE(all.equal(grid[[earlier]], grid[[criterion[index]]]))
    }, logical(1)))
  }, logical(1))
  kept <- criterion[!repeats]
  stats::setNames(kept, names(views$criteria)[match(kept, views$criteria)])
}

#' The share of the smallest class of every candidate, at each level
#'
#' A profile, class or transition grid takes each fit's effective class
#' counts over the units it classifies; a regression grid reads the
#' `smallest_share` its table carries for its classes, and has no group-level
#' share. A candidate with fewer than two classes at a level has none (`NA`).
#'
#' @param x An enumeration grid offering the `"smallest"` view.
#' @return A list of two numeric vectors, `individuals` and `groups`, one value
#'   per candidate.
#' @noRd
.latents_grid_smallest_shares <- function(x) {
  if (identical(.latents_grid_kind(x), "regressions")) {
    table <- .subset2(x, "table")
    return(list(individuals = ifelse(table$n_classes >= 2L, table$smallest_share,
                                     NA_real_),
                groups = rep(NA_real_, nrow(table))))
  }
  counts_field <- c(individuals = "effective_profile_counts",
                    groups = "effective_group_counts")
  share <- function(fit, level) {
    # The rows and groups the fit counts (see `drop_empty`).
    posteriors <- if (identical(level, "individuals"))
      .multilpa_classified_posteriors(fit) else .multilpa_classified_group_posteriors(fit)
    if (is.null(posteriors)) return(NA_real_)
    weights <- .multilpa_level_weights(fit, level, nrow(posteriors)) %||% 1
    counts <- fit[[counts_field[[level]]]] %||%
      colSums(posteriors[, seq_len(if (identical(level, "individuals"))
        ncol(fit$subject_posteriors) else ncol(posteriors)), drop = FALSE] * weights)
    if (length(counts) < 2L) return(NA_real_)
    min(counts) / sum(posteriors * weights)
  }
  fits <- .subset2(x, "fits")
  lapply(stats::setNames(names(counts_field), names(counts_field)), function(level) {
    vapply(fits, function(fit) {
      if (is.null(fit)) NA_real_ else share(fit, level)
    }, numeric(1L), USE.NAMES = FALSE)
  })
}

#' Draw a transition, family or regression grid through the shared views
#'
#' The `plot()` methods of these grids keep their 0.9.21 call: `what` (or its
#' former name `criterion`) names a criterion, matched as `match.arg()`
#' matched `"bic"`, `"aic"` and `"bic_individual"`, so a unique abbreviation
#' such as `"bic_i"` still works. Without either, the grid's summary of its
#' default criteria is drawn, as [plot_enumeration()] draws it. Any other value is
#' passed to [plot_enumeration()], which draws every view the grid offers.
#' Arguments of [plot_enumeration()] and, for the tree, of [plot_tree()] given
#' through `...` are passed on; any other is ignored with a
#' `latents_unused_argument` warning, as 0.9.21 ignored it.
#'
#' @param x A transition or group-class family grid.
#' @param what The caller's `what`.
#' @param what_given Whether the caller supplied `what`.
#' @param criterion The caller's `criterion`, `NULL` when not given.
#' @param main,subtitle Titles.
#' @param extra `list(...)` as the method received it.
#' @return A ggplot object, or a `latents_plots` list.
#' @noRd
.latents_plot_former_grid <- function(x, what, what_given, criterion, main,
                                      subtitle, extra) {
  views <- .latents_grid_views(x, "plot")
  # Without `what` or `criterion` the grid's summary is drawn, as
  # plot_enumeration() draws it: every default criterion as a panel.
  what <- if (!what_given && is.null(criterion)) views$default else
    .latents_grid_criterion(.latents_former_argument(
      what, what_given, criterion, c("what", "criterion")), views)
  supplied <- names(extra) %||% rep("", length(extra))
  passed <- c("combine", "labels", "mark_minimum", "cutoff", "using",
              if (identical(what, "tree")) .latents_tree_arguments())
  forwarded <- nzchar(supplied) & supplied %in% passed
  .latents_warn_unused_arguments(extra[!forwarded], "plot()",
                                 .latents_grid_plot_advice())
  do.call(plot_enumeration, c(list(x, what = what, main = main,
                                   subtitle = subtitle), extra[forwarded]))
}
