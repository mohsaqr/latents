# Comparing fitted models on the same data, in one call.

utils::globalVariables(c("model", "criterion", "value"))

#' Compare fitted models on the same data
#'
#' One row per model with its log likelihood, number of parameters,
#' information criteria, classification quality and convergence, so
#' alternative specifications (for example trajectory classes alone against
#' a growth mixture model) are compared without assembling tables by hand.
#' `delta_bic` is each model's BIC above the smallest, and `bic_weight` the
#' Schwarz weight, exp(-delta/2) normalized over the models compared (Wagenmakers
#' and Farrell 2004): the approximate posterior probability of each model if
#' one of them is true, with equal prior probabilities.
#'
#' Information criteria are comparable only on the same data: every model
#' must be fitted to the same outcome, persons (or groups) and observations,
#' or the comparison is refused.
#'
#' Candidates of an enumeration grid are compared by naming them, as the
#' grid's `candidate` column shows them:
#' `compare_models(grid, "classes3_groups4", "classes3_groups5")`, or by their
#' class counts, `compare_models(grid, c(3, 4), c(3, 5))` (one number is one
#' group class). `n_profiles` and `n_group_classes` select candidates instead:
#' `compare_models(grid, n_profiles = 3)` compares every candidate with 3
#' profiles; `compare_models(grid)` compares them all. Latent class models are
#' compared the same way; their `n_profiles` column counts classes. Fits from [multilpa()], [lpa()], [lca()], [lta()] or [get_candidate()]
#' can also be passed directly.
#'
#' @param ... Either an enumeration grid from [enumerate_classes()] followed
#'   by the candidates to compare (none for all of them, or a selection
#'   through `n_profiles` and `n_group_classes`); or two or more fits from [multilpa()],
#'   [lpa()], [lca()], [lta()] or [get_candidate()]; or two or more fits from
#'   [mixture_regression()] (with or without `random`). Name fits to label
#'   the rows; unnamed fits are labelled by their expression, candidates by
#'   their name.
#' @param n_profiles,n_group_classes With a grid, the class counts of the
#'   candidates to compare, instead of naming them; `NULL` takes every count.
#'   Candidates whose fit failed are left out, with a message.
#' @param cutoff For profile and class models, the average posterior
#'   probability below which a class is counted in `avepp_below` and
#'   `group_avepp_below`. Not used for mixture regressions.
#' @return A base `data.frame` of class `latents_comparison` with one row per
#'   model. For mixture regressions: `model`, `n_classes`, `random`,
#'   `n_parameters`, `log_likelihood`, `aic`, `bic`, `sabic`, `icl`,
#'   `entropy`, `smallest_share`, `delta_bic`, `bic_weight`, `converged`. For
#'   profile and class models: `model`, `n_profiles`, `n_group_classes`,
#'   `structure` (the covariance structure, `NA` without continuous
#'   indicators), `n_parameters`, `log_likelihood`, every information
#'   criterion of [enumerate_classes()]'s candidates table, `profile_entropy`,
#'   `group_entropy`, `avepp_min`, `avepp_max`, `avepp_mean`, `avepp_below`
#'   and the same four prefixed `group_`, `smallest_share` (the smallest estimated profile
#'   share, relative to all observation mass including any noise component), `converged`, `n_best_replicated`, and `delta_bic` and
#'   `bic_weight` on `bic_groups` (the BIC that counts groups as the
#'   independent units; on a single-level model it is the ordinary BIC). The
#'   numbers are the ones the grid's candidates table reports. The printout
#'   marks the row with the lowest BIC; that is a fact about the criterion,
#'   not a choice. Raises `latents_bad_argument` for fits it cannot compare,
#'   for a group-class family grid (compare those in their own table), for
#'   candidates named together with `n_profiles` or `n_group_classes`, and for
#'   `n_profiles` or `n_group_classes` without a grid;
#'   `latents_unknown_candidate` for a name or counts not in the grid;
#'   `latents_failed_candidate` for a named candidate that could not be
#'   fitted; and `latents_incomparable_models` for fits to different data. A
#'   `cutoff` outside (0, 1] fails the argument checks.
#' @references
#' Wagenmakers, E.-J., & Farrell, S. (2004). AIC model selection using Akaike
#' weights. *Psychonomic Bulletin & Review*, 11, 192--196.
#' @examples
#' # Fixed trajectories against growth curves with a random intercept,
#' # for the first 60 students
#' few_students <- subset(growth_scores, student <= 60)
#' trajectories <- mixture_regression(score ~ wave, few_students, n_classes = 2,
#'                                    id = "student", class_level = "group",
#'                                    seed = 1)
#' growth <- mixture_regression(score ~ wave, few_students, n_classes = 2,
#'                              id = "student", class_level = "group",
#'                              random = "intercept", seed = 1)
#' compare_models(trajectories = trajectories, growth = growth)
#'
#' # Candidates of an enumeration grid, by name, by class counts, or selected
#' grid <- enumerate_classes(subset(course_engagement, student <= 40),
#'                           c("browse", "lectures"), "student",
#'                           n_profiles = 2:3, n_group_classes = 1:2,
#'                           model = "VVI", n_starts = 2, seed = 1)
#' compare_models(grid, "profiles2_groups2", "profiles3_groups2")
#' compare_models(grid, c(2, 2), c(3, 2))
#' compare_models(grid, n_group_classes = 2)
#' @export
compare_models <- function(..., n_profiles = NULL, n_group_classes = NULL,
                           cutoff = 0.8) {
  fits <- list(...)
  expressions <- vapply(as.list(substitute(list(...)))[-1L], deparse1, character(1))
  labels <- names(fits) %||% rep("", length(fits))
  labels[!nzchar(labels)] <- expressions[!nzchar(labels)]
  if (length(fits) >= 1L && inherits(fits[[1L]], "latents_family_enumeration")) {
    stop(errorCondition(paste(
      "Group-class family grids are compared in their own table:",
      "as.data.frame(grid)."), class = "latents_bad_argument", call = NULL))
  }
  if (length(fits) >= 1L && .latents_is_grid(fits[[1L]])) {
    return(.multilpa_compare_candidates(fits[[1L]], fits[-1L], labels[-1L],
                                        n_profiles, n_group_classes, cutoff))
  }
  if (!is.null(n_profiles) || !is.null(n_group_classes)) {
    stop(errorCondition(paste(
      "`n_profiles` and `n_group_classes` select candidates of an enumeration",
      "grid; pass the grid first."), class = "latents_bad_argument", call = NULL))
  }
  if (length(fits) >= 1L && all(vapply(fits, .latents_comparison_measurement_fit, logical(1)))) {
    return(.multilpa_compare_fits(fits, labels, cutoff))
  }
  supported <- vapply(fits, function(fit) {
    inherits(fit, "latents_mixture_regression") || inherits(fit, "latents_growth_mixture")
  }, logical(1))
  if (length(fits) < 2L || !all(supported)) {
    stop(errorCondition(paste(
      "`compare_models()` takes an enumeration grid and two or more candidate",
      "keys, or two or more fits of one kind: from multilpa(), lpa(), lca() or",
      "get_candidate(), or from mixture_regression()."),
      class = "latents_bad_argument", call = NULL))
  }
  if (anyDuplicated(labels)) {
    stop(errorCondition("Every model needs a distinct name.",
                        class = "latents_bad_argument", call = NULL))
  }
  signature <- function(fit) {
    spec <- fit$spec
    list(y = spec$y, trials = spec$trials, family = spec$family,
         categories = spec$category_levels, groups = spec$group_index,
         group_levels = spec$group_levels, clusters = spec$cluster_index,
         cluster_levels = spec$cluster_levels,
         weights = spec$sampling_weights, rows = spec$kept_rows,
         response = spec$response_name)
  }
  reference <- signature(fits[[1L]])
  same <- vapply(fits, function(fit) isTRUE(all.equal(signature(fit), reference)),
                 logical(1))
  if (!all(same)) {
    stop(errorCondition(sprintf(paste(
      "%s %s fitted to different data from %s: information criteria compare",
      "models only on the same outcome and observations."),
      paste(labels[!same], collapse = ", "), if (sum(!same) > 1L) "are" else "is",
      labels[1L]), class = "latents_incomparable_models", call = NULL))
  }
  rows <- lapply(seq_along(fits), function(i) {
    fit_row <- get_results(fits[[i]], "fit")
    data.frame(model = labels[i], n_classes = as.integer(fit_row$n_classes),
               random = .growth_random_description(fits[[i]]$spec),
               n_parameters = as.integer(fit_row$n_parameters),
               log_likelihood = fit_row$log_likelihood, aic = fit_row$aic,
               bic = fit_row$bic, sabic = fit_row$sabic, icl = fit_row$icl,
               entropy = fit_row$entropy, smallest_share = fit_row$smallest_share,
               converged = isTRUE(fit_row$converged), stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  table$delta_bic <- table$bic - min(table$bic)
  relative <- exp(-table$delta_bic / 2)
  table$bic_weight <- relative / sum(relative)
  table <- table[c("model", "n_classes", "random", "n_parameters", "log_likelihood",
                   "aic", "bic", "sabic", "icl", "entropy", "smallest_share",
                   "delta_bic", "bic_weight", "converged")]
  column <- .latents_column
  table <- .latents_table(table, title = "Model comparison", display = list(
    column("Model", "model", "text"), column("Classes", "n_classes", "integer"),
    column("Random", "random", "text"), column("k", "n_parameters", "integer"),
    column("LogLik", "log_likelihood"), column("BIC", "bic"),
    column(if (isTRUE(l10n_info()$`UTF-8`)) "\u0394BIC" else "dBIC", "delta_bic"), column("Weight", "bic_weight", "share"),
    column("Entropy", "entropy")))
  # Marks the lowest BIC, which is a fact about the criterion, not a choice.
  attr(table, "marks") <- ifelse(seq_len(nrow(table)) == which.min(table$delta_bic),
                                 "  <- lowest BIC", "")
  attr(table, "models") <- stats::setNames(fits, labels)
  class(table) <- c("latents_comparison", class(table))
  table
}

#' @rdname latents-plots
#' @export
plot.latents_comparison <- function(x, main = NULL, subtitle = NULL, ...) {
  .latents_warn_unused_arguments(
    list(...), "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  .gg_require()
  table <- as.data.frame(x)
  # Profile and class models report each criterion under the groups
  # convention, mixture regressions under one name.
  criteria <- c(aic = "AIC", bic = "BIC", bic_groups = "BIC", sabic = "SABIC",
                sabic_groups = "SABIC")
  criteria <- criteria[names(criteria) %in% names(table)]
  long <- do.call(rbind, lapply(names(criteria), function(column) {
    data.frame(model = factor(table$model, levels = table$model),
               criterion = factor(criteria[[column]], levels = unname(criteria)),
               value = table[[column]], best = table[[column]] == min(table[[column]]),
               stringsAsFactors = FALSE)
  }))
  ggplot2::ggplot(long, ggplot2::aes(model, value, group = criterion,
                                     colour = criterion, shape = criterion)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_point(data = long[long$best, , drop = FALSE], size = 5.2, shape = 21,
                        fill = NA, colour = "black", stroke = 0.6) +
    ggplot2::scale_colour_manual(values = .gg_okabe_ito[c(1L, 4L, 3L)]) +
    ggplot2::scale_shape_manual(values = c(16L, 17L, 15L)) +
    ggplot2::labs(x = NULL, y = "Information criterion (lower is better)") +
    .gg_theme() + .gg_top_legend() +
    .gg_titles(main, subtitle, "Model comparison",
               "Each criterion's best model is circled")
}

#' Compare candidates of an enumeration grid
#'
#' Candidates are named one by one (names or class counts), or chosen by
#' class counts with `n_profiles` and `n_group_classes`; naming none compares
#' every candidate. A named candidate whose fit failed is an error; one that
#' a selection takes in is left out, with a message saying so.
#'
#' @param grid A `multilpa_enumeration` or `latents_transition_enumeration`.
#' @param chosen The candidates as `compare_models()` received them.
#' @param labels Their names, where the caller named them.
#' @param n_profiles,n_group_classes Class counts to select, or `NULL`.
#' @param cutoff The AvePP cutoff.
#' @return The comparison table.
#' @noRd
.multilpa_compare_candidates <- function(grid, chosen, labels, n_profiles,
                                         n_group_classes, cutoff) {
  names <- .latents_candidate_names(grid)
  table <- .subset2(grid, "table")
  if (length(chosen) > 0L) {
    if (!is.null(n_profiles) || !is.null(n_group_classes)) {
      stop(errorCondition(paste(
        "Name the candidates, or select them with `n_profiles` and",
        "`n_group_classes`, not both."),
        class = "latents_bad_argument", call = NULL))
    }
    rows <- vapply(chosen, function(candidate) .latents_candidate_row(grid, candidate),
                   integer(1L))
    fits <- lapply(rows, function(at) .latents_candidate_at(grid, at))
    given <- nzchar(names(chosen) %||% rep("", length(chosen)))
    labels <- ifelse(given, labels, names[rows])
  } else {
    counts <- function(value) {
      is.null(value) || (is.numeric(value) && length(value) >= 1L &&
                           all(is.finite(value)) && all(value >= 1) &&
                           all(value == floor(value)))
    }
    if (!counts(n_profiles) || !counts(n_group_classes)) {
      stop(errorCondition(paste(
        "`n_profiles` and `n_group_classes` must be NULL or positive whole",
        "numbers."), class = "latents_bad_argument", call = NULL))
    }
    rows <- which((is.null(n_profiles) | table$n_profiles %in% n_profiles) &
                    (is.null(n_group_classes) |
                       table$n_group_classes %in% n_group_classes))
    if (length(rows) == 0L) {
      stop(errorCondition(sprintf(
        "No candidate matches that selection. The grid has %s.",
        paste(sprintf("\"%s\"", names), collapse = ", ")),
        class = "latents_unknown_candidate", call = NULL))
    }
    failed <- rows[vapply(rows, function(at) is.null(.subset2(grid, "fits")[[at]]),
                          logical(1L))]
    if (length(failed) > 0L) {
      message(sprintf("Left out %s, which could not be fitted.",
                      paste(sprintf("\"%s\"", names[failed]), collapse = ", ")))
    }
    rows <- setdiff(rows, failed)
    fits <- lapply(rows, function(at) .subset2(grid, "fits")[[at]])
    labels <- names[rows]
  }
  if (length(fits) < 2L) {
    stop(errorCondition(sprintf(paste(
      "A comparison needs two or more fitted candidates; this one has %d.",
      "The grid has %s."), length(fits),
      paste(sprintf("\"%s\"", names), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  .multilpa_compare_fits(fits, labels, cutoff)
}

#' Compare profile and class models fitted to the same data
#'
#' @param fits Fitted models of class `multilpa`.
#' @param labels Their row labels.
#' @param cutoff The AvePP cutoff.
#' @return The comparison table.
#' @noRd
.multilpa_compare_fits <- function(fits, labels, cutoff) {
  if (length(fits) < 2L) {
    stop(errorCondition("`compare_models()` needs two or more fits.",
                        class = "latents_bad_argument", call = NULL))
  }
  if (anyDuplicated(labels)) {
    stop(errorCondition("Every model needs a distinct name.",
                        class = "latents_bad_argument", call = NULL))
  }
  stopifnot("`cutoff` must be a single number between 0 and 1" =
              is.numeric(cutoff) && length(cutoff) == 1L && is.finite(cutoff) &&
              cutoff > 0 && cutoff <= 1)
  # The data a likelihood is computed on: the same units in the same groups
  # with the same indicator values.
  signature <- .latents_comparison_signature
  reference <- signature(fits[[1L]])
  same <- vapply(fits, function(fit) isTRUE(all.equal(signature(fit), reference)),
                 logical(1))
  if (!all(same)) {
    stop(errorCondition(sprintf(paste(
      "%s %s fitted to different data from %s: information criteria compare",
      "models only on the same indicators and observations."),
      paste(labels[!same], collapse = ", "), if (sum(!same) > 1L) "are" else "is",
      labels[1L]), class = "latents_incomparable_models", call = NULL))
  }
  statistics <- c("min", "max", "mean", "below")
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    level_summary <- function(posteriors, prefix) {
      values <- .multilpa_avepp_summary(posteriors, cutoff)
      columns <- as.data.frame(as.list(stats::setNames(values,
                                                       paste0(prefix, statistics))))
      columns[[paste0(prefix, "below")]] <- as.integer(values[["below"]])
      columns
    }
    diagnostic <- .latents_comparison_diagnostic_view(fit)
    counts <- fit$effective_profile_counts %||%
      colSums(fit$subject_posteriors *
        (.multilpa_level_weights(diagnostic, "individuals", nrow(fit$subject_posteriors)) %||% 1))
    # Noise is outside the Gaussian profiles, but still occupies sample mass.
    # Profile shares are relative to all classified observations, not to the
    # Gaussian portion alone.
    total <- sum(.multilpa_class_posteriors(fit) *
      (.multilpa_level_weights(diagnostic, "individuals", nrow(fit$subject_posteriors)) %||% 1))
    shares <- counts / total
    cbind(
      data.frame(model = labels[i], n_profiles = as.integer(fit$n_profiles),
                 n_group_classes = as.integer(fit$n_group_classes %||% 1L),
                 structure = if (length(.multilpa_continuous_names(fit)) == 0L)
                   NA_character_ else .latents_comparison_structure(fit),
                 n_parameters = as.integer(fit$n_parameters),
                 log_likelihood = fit$log_likelihood, stringsAsFactors = FALSE),
      .multilpa_enumeration_indices(diagnostic),
      data.frame(profile_entropy = .multilpa_relative_entropy(.multilpa_class_posteriors(fit)),
                 group_entropy = if (is.null(fit$group_posteriors)) NA_real_ else
                   .multilpa_relative_entropy(fit$group_posteriors)),
      level_summary(.multilpa_class_posteriors(fit), "avepp_"),
      level_summary(fit$group_posteriors, "group_avepp_"),
      data.frame(smallest_share = if (length(shares) == 0L) NA_real_ else min(shares),
                 converged = isTRUE(fit$converged),
                 n_best_replicated = fit$n_best_replicated %||% NA_integer_))
  })
  table <- do.call(rbind, rows)
  table$delta_bic <- table$bic_groups - min(table$bic_groups)
  relative <- exp(-table$delta_bic / 2)
  table$bic_weight <- relative / sum(relative)
  column <- .latents_column
  delta <- if (isTRUE(l10n_info()$`UTF-8`)) "\u0394BIC" else "dBIC"
  # The printed view names latent classes as classes, and leaves out the
  # group columns when no model compared has group classes; the table itself
  # keeps every column.
  latent_class <- all(vapply(fits, function(fit)
    length(.multilpa_continuous_names(fit)) == 0L, logical(1)))
  grouped <- any(table$n_group_classes > 1L)
  table <- .latents_table(table, title = "Model comparison", display = c(
    list(column("Model", "model", "text"),
         column(if (latent_class) "Classes" else "Profiles", "n_profiles",
                "integer")),
    if (grouped) list(column("Group classes", "n_group_classes", "integer")),
    list(column("k", "n_parameters", "integer"), column("LogLik", "log_likelihood"),
         column("BIC", "bic_groups"), column(delta, "delta_bic"),
         column("Weight", "bic_weight", "share"),
         column("ICL", "icl_individual"), column("Entropy", "profile_entropy"),
         column("Min AvePP", "avepp_min"),
         column(sprintf("< %g", cutoff), "avepp_below", "integer")),
    if (grouped) list(column("Group entropy", "group_entropy"),
                      column("Group min AvePP", "group_avepp_min"))))
  # Marks the lowest BIC, which is a fact about the criterion, not a choice.
  attr(table, "marks") <- ifelse(seq_len(nrow(table)) == which.min(table$delta_bic),
                                 "  <- lowest BIC", "")
  attr(table, "models") <- stats::setNames(fits, labels)
  class(table) <- c("latents_comparison", class(table))
  table
}

#' The model behind one row of a comparison
#'
#' @param x A `latents_comparison`.
#' @param label The row's `model` label.
#' @return The fitted model.
#' @noRd
.latents_compared_model <- function(x, label) {
  if (missing(label)) {
    stop(errorCondition(paste(
      "Name the row: get_candidate(comparison, \"<model>\"), as the",
      "comparison's `model` column shows it."),
      class = "latents_bad_argument", call = NULL))
  }
  models <- attr(x, "models")
  if (is.null(models)) {
    stop(errorCondition(paste(
      "This comparison carries no models: it was made by an earlier version,",
      "or rows were taken out of it. Compare the models again."),
      class = "latents_incomplete_fit", call = NULL))
  }
  # Rows taken out of the comparison take their models with them.
  present <- intersect(names(models), x$model)
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !label %in% present) {
    stop(errorCondition(sprintf(
      "%s is not a model of this comparison. It has %s.", deparse1(label),
      paste(sprintf("\"%s\"", present), collapse = ", ")),
      class = "latents_unknown_candidate", call = NULL))
  }
  models[[label]]
}

#' Per-class classification quality of every compared model
#'
#' @param x A `latents_comparison` of profile or class models.
#' @return One row per model, level and class: `model`, then a fit's
#'   `"classification"` columns and that level's `relative_entropy`.
#' @noRd
.latents_comparison_classification <- function(x) {
  models <- attr(x, "models")
  if (is.null(models)) {
    stop(errorCondition(paste(
      "This comparison carries no models: it was made by an earlier version,",
      "or rows were taken out of it. Compare the models again."),
      class = "latents_incomplete_fit", call = NULL))
  }
  # Subsetting and reordering data frames keeps attributes. Read only the
  # models whose rows remain, in their current order.
  models <- models[x$model[x$model %in% names(models)]]
  supported <- vapply(models, .latents_comparison_measurement_fit, logical(1))
  if (!all(supported)) {
    stop(errorCondition(paste(
      "Per-class classification is reported for profile, class and transition",
      "models; a mixture regression reports its own:",
      "get_results(fit, \"classification\")."),
      class = "latents_bad_argument", call = NULL))
  }
  rows <- lapply(names(models), function(label) {
    quality <- .multilpa_informative_classification(
      .latents_comparison_diagnostic_view(models[[label]]))
    cbind(data.frame(model = rep(label, nrow(quality))), quality)
  })
  if (length(rows) == 0L) {
    empty <- .multilpa_enumeration_classification(
      structure(list(fits = list()), class = "multilpa_enumeration"))
    empty[c("candidate", "n_profiles", "n_group_classes", "model")] <- NULL
    return(cbind(data.frame(model = character()), empty))
  }
  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result
}

#' Measurement families that share a model-comparison table
#' @noRd
.latents_comparison_measurement_fit <- function(x) {
  inherits(x, c("multilpa", "multilpa_covariates", "multilpa_transitions", "multilpa_lta"))
}

#' Posterior and likelihood fields for the shared diagnostic formulas
#'
#' General transitions store these fields under the same meaning as the
#' homogeneous engine. This private read-only view contains only diagnostic
#' inputs; it is never returned as a fitted model or used for inference.
#' @noRd
.latents_comparison_diagnostic_view <- function(fit) {
  if (!inherits(fit, "multilpa_lta")) return(fit)
  fields <- c("n_parameters", "log_likelihood", "n_groups", "n_observations",
              "n_informative", "subject_posteriors", "group_posteriors",
              "sampling_weights", "weights", "group_index")
  structure(fit[intersect(fields, names(fit))], class = "multilpa_transitions")
}

#' Covariance code even for families that store only the older switches
#' @noRd
.latents_comparison_structure <- function(fit) {
  fit$covariance_structure %||% .multilpa_resolve_structure(
    fit$variance_model %||% "varying", fit$covariance_model %||% "diagonal")
}

#' Original observations, their coding and the transformation the model used
#' @noRd
.latents_comparison_signature <- function(fit) {
  general <- inherits(fit, "multilpa_lta")
  raw <- if (general) sweep(fit$x, 2L, fit$centers, "+") else
    .multilpa_uncentered_indicators(fit)
  # Engines use different storage types for an absent Gaussian block.
  if (!is.null(raw) && ncol(raw) == 0L) raw <- NULL
  list(vars = fit$vars, id = if (isTRUE(fit$single_level)) NULL else fit$id,
       groups = as.character(fit$group_values %||% fit$group_ids),
       group_index = fit$group_index, indicators = raw,
       categorical = if (general) fit$codes else fit$categorical_data,
       categorical_levels = fit$categorical_levels,
       extra = .latents_draw_extra_columns(fit$extra_data),
       weights = unname(fit$sampling_weights), time = fit$time_values,
       centering = fit$centering %||% "none")
}
