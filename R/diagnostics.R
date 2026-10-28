#' Classification entropy sum
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @param weights `NULL`, or one sampling weight per row.
#' @return Sum of `-p log p` over every cell, the classification entropy;
#'   with `weights`, each row's entropy multiplied by its weight.
#' @noRd
.multilpa_entropy_sum <- function(probabilities, weights = NULL) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities),
            "`weights` must be NULL or one finite weight per row" =
              is.null(weights) || (is.numeric(weights) &&
                                     length(weights) == nrow(probabilities) &&
                                     all(is.finite(weights))))
  cells <- ifelse(probabilities > 0, probabilities * log(probabilities), 0)
  if (is.null(weights)) return(-sum(cells))
  -sum(weights * rowSums(cells))
}

#' Relative classification entropy
#'
#' `1 - EN / (n log K)` (Ramaswamy et al. 1993). Under sampling weights `EN`
#' is the weighted entropy sum and `n` the sum of the weights, so the value is
#' that of the population the weights describe (review 2026-10-26 round 2,
#' S-4).
#'
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @param weights `NULL`, or one sampling weight per row.
#' @return Entropy on the zero-to-one scale, or `NA_real_` for a single class.
#' @noRd
.multilpa_relative_entropy <- function(probabilities, weights = NULL) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  if (ncol(probabilities) == 1L) return(NA_real_)
  total <- if (is.null(weights)) nrow(probabilities) else sum(weights)
  1 - .multilpa_entropy_sum(probabilities, weights) /
    (total * log(ncol(probabilities)))
}

#' The formula behind each criterion, keyed by its name
#'
#' `EN` in the entropy-charged criteria is the soft entropy of the
#' complete-data classification: of the group classes and the profiles
#' together, the same under both sample-size conventions. A transition model
#' keeps one entropy per level (its profiles form a Markov chain within a
#' person), and its definitions say so.
#'
#' @param criterion Character vector of criterion names.
#' @param transition `TRUE` for a transition model.
#' @return Character vector of formulas, the same length as `criterion`.
#' @noRd
.multilpa_criterion_definitions <- function(criterion, transition = FALSE) {
  entropy <- if (isTRUE(transition)) {
    paste("; EN = soft classification entropy of this convention's level",
          "(group classes for groups, occasions' profiles for individuals)")
  } else {
    paste("; EN = soft entropy of the complete-data classification (group",
          "classes and profiles), the same under both conventions",
          "(Biernacki, Celeux & Govaert 2000)")
  }
  definitions <- c(deviance = "-2L",
                   aic = "-2L + 2q",
                   kic = "-2L + 3(q + 1) (Cavanaugh 1999, as tidyLPA computes it)",
                   bic = "-2L + q log(n)",
                   sabic = "-2L + q log((n + 2) / 24)",
                   caic = "-2L + q (log(n) + 1)",
                   awe = paste0("-2(L - EN) + 2q (1.5 + log(n))", entropy),
                   icl = paste0("-2L + q log(n) + 2 EN", entropy),
                   clc = paste0("-2L + 2 EN", entropy))
  stopifnot("Every criterion must have a definition" =
              all(criterion %in% names(definitions)))
  unname(definitions[criterion])
}

#' Every likelihood-penalty criterion a fit supports
#'
#' The formulas, the two sample-size conventions, the entropy each criterion
#' penalizes and the references are documented on `?get_results`, which is where
#' a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @param definitions `TRUE` adds the formula and the reference per criterion,
#'   which only the long form has a row for.
#' @param format `"wide"` for the one-row reporting shape, `"long"` for one row
#'   per criterion and sample-size convention.
#' @return A base `data.frame`, one row in the wide form and one row per
#'   criterion and convention in the long form.
#' @noRd
.multilpa_information_criteria <- function(x, definitions = FALSE,
                                 format = c("wide", "long")) {
  .latents_require_fit(x)
  .latents_require("`definitions` must be TRUE or FALSE" =
              isTRUE(definitions) || isFALSE(definitions))
  format <- .latents_match_arg(format, c("wide", "long"), "format")
  if (isTRUE(definitions) && identical(format, "wide")) {
    stop(errorCondition(
      "`definitions` describes one criterion per row, which the wide form has not. Use format = \"long\".",
      class = "latents_bad_argument", call = NULL))
  }
  q <- x$n_parameters
  deviance <- -2 * x$log_likelihood
  conventions <- data.frame(
    convention = c("groups", "individuals"),
    n = c(.multilpa_group_sample_size(x),
          x$n_informative %||% x$n_observations))
  # The entropy charge of ICL, AWE and CLC is that of the same likelihood the
  # deviance is: under sampling weights the weighted sum, which equals the
  # entropy of the data with every unit repeated by its weight (review
  # 2026-10-26 round 2, EN-1). It is the entropy of the complete-data
  # classification, group classes and profiles together, so the two
  # conventions differ only in `n` (round 2, S-1). A transition model keeps
  # one entropy per level: its profiles are a Markov chain within a person,
  # whose path entropy is not the sum of the occasions' entropies.
  transition <- inherits(x, "multilpa_transitions")
  level_entropy <- function(level, posteriors) {
    if (is.null(posteriors)) return(NA_real_)
    .multilpa_entropy_sum(posteriors,
                          .multilpa_level_weights(x, level, nrow(posteriors)))
  }
  entropies <- if (transition) {
    c(groups = level_entropy("groups", .multilpa_classified_group_posteriors(x)),
      individuals = level_entropy("individuals", .multilpa_classified_posteriors(x)))
  } else {
    complete <- .multilpa_complete_entropy(x)
    c(groups = complete, individuals = complete)
  }
  scale_free <- data.frame(
    criterion = c("deviance", "aic", "kic"),
    convention = NA_character_, n = NA_integer_,
    penalty = c(0, 2 * q, 3 * (q + 1)))
  scaled <- do.call(rbind, lapply(seq_len(nrow(conventions)), function(row) {
    convention <- conventions$convention[row]
    n <- conventions$n[row]
    entropy_sum <- entropies[[convention]]
    penalties <- c(bic = q * log(n),
                   sabic = q * log((n + 2) / 24),
                   caic = q * (log(n) + 1),
                   awe = 2 * q * (1.5 + log(n)) + 2 * entropy_sum,
                   icl = q * log(n) + 2 * entropy_sum,
                   # The only criterion here that does not use the sample size:
                   # the same under both conventions, except for a transition
                   # model, whose convention picks the level's entropy.
                   clc = 2 * entropy_sum)
    data.frame(criterion = names(penalties), convention = convention, n = n,
               penalty = unname(penalties))
  }))
  result <- rbind(scale_free, scaled)
  # Every row is the same deviance plus its own penalty, so the `value` column
  # means one thing throughout: smaller is a better-supported model. The penalty
  # itself is not a column: nobody reports it, and it is the difference between
  # two numbers the table already carries.
  result$value <- deviance + result$penalty
  result <- result[c("criterion", "convention", "n", "value")]
  row.names(result) <- NULL
  if (identical(format, "wide")) return(.multilpa_criteria_wide(result, x))
  if (isTRUE(definitions)) {
    result$definition <- .multilpa_criterion_definitions(result$criterion,
                                                         transition)
  }
  result
}

#' One criterion per column, the shape a model-comparison table is reported in
#'
#' The same pivot the enumeration grid uses, so a single fit and a row of
#' `as.data.frame(enumerate_classes(...))` cannot name the same quantity
#' differently. `.multilpa_enumeration_criteria()` owns the column order.
#'
#' @param indices The long table, one row per criterion and convention.
#' @param fit The model the indices came from, for the likelihood and count.
#' @return A one-row base `data.frame`.
#' @noRd
.multilpa_criteria_wide <- function(indices, fit) {
  labels <- .multilpa_criterion_labels(indices)
  wanted <- .multilpa_enumeration_criteria()
  missing_labels <- setdiff(wanted, labels)
  if (length(missing_labels) > 0L) {
    stop(errorCondition(sprintf(
      "The information criteria no longer report %s.",
      paste(missing_labels, collapse = ", ")),
      class = "latents_unknown_criterion", call = NULL))
  }
  values <- stats::setNames(indices$value, labels)
  cbind(data.frame(log_likelihood = fit$log_likelihood,
                   n_parameters = fit$n_parameters),
        as.data.frame(as.list(values[wanted])))
}

#' The wide column name each long row maps to
#'
#' A criterion that uses no sample size carries `NA` as its convention and keeps
#' its bare name; the rest take a `_groups` or `_individual` suffix.
#'
#' @param indices The long table.
#' @return A character vector, one label per row of `indices`.
#' @noRd
.multilpa_criterion_labels <- function(indices) {
  no_convention <- is.na(indices$convention) | indices$convention %in% "none"
  ifelse(no_convention, indices$criterion,
         paste(indices$criterion,
               sub("individuals", "individual", indices$convention), sep = "_"))
}

#' The posterior matrices of the levels a diagnostic was asked for
#' @param object A fitted model of this package.
#' @param level One of `"individuals"`, `"groups"` or `"both"`.
#' @return A named list of posterior matrices, in reporting order.
#' @noRd
.multilpa_posterior_levels <- function(object, level) {
  posteriors <- list(individuals = .multilpa_classified_posteriors(object),
                     groups = .multilpa_classified_group_posteriors(object))
  posteriors[.multilpa_classification_levels(object, level)]
}

#' Which levels a request for one, the other or both resolves to
#'
#' `"both"` on a model with no discrete group classes is the individual level
#' rather than a refusal: the caller asked for every level the fit has, and it
#' has one. Naming `"groups"` explicitly is a different request and refuses.
#'
#' @param object A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A character vector of level names, in reporting order.
#' @noRd
.multilpa_classification_levels <- function(object, level) {
  levels_wanted <- if (identical(level, "both")) c("individuals", "groups") else level
  if (!is.null(object$group_posteriors)) return(levels_wanted)
  if (identical(level, "groups")) {
    stop(errorCondition("This model has no discrete group classes.",
                        class = "latents_no_group_classes", call = NULL))
  }
  "individuals"
}

#' Per-class classification quality at one or both levels
#'
#' The columns, the odds-of-correct-classification formula and the reference
#' are documented on `?get_results`, which is where a caller reaches this table
#' from.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A base `data.frame`, one row per level and class.
#' @noRd
.multilpa_classification_table <- function(x, level = c("individuals",
                                                        "groups", "both")) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(x))
  level <- .latents_match_arg(level, c("individuals", "groups", "both"), "level")
  posteriors <- .multilpa_posterior_levels(x, level)
  result <- do.call(rbind, lapply(names(posteriors), function(which_level) {
    .multilpa_classification_summary(
      posteriors[[which_level]], which_level,
      .multilpa_level_weights(x, which_level, nrow(posteriors[[which_level]])))
  }))
  row.names(result) <- NULL
  result
}

#' The sampling weight of every row of one level's posterior matrix
#'
#' Weights belong to the independent units: the groups' own posterior rows
#' carry them directly, and an individual row carries its unit's weight.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"` or `"groups"`.
#' @param n_rows Rows of that level's posterior matrix.
#' @return `NULL` for an unweighted fit, otherwise one weight per row.
#' @noRd
.multilpa_level_weights <- function(x, level, n_rows) {
  if (!.latents_is_weighted(x)) return(NULL)
  units <- unname(x$sampling_weights)
  weights <- if (identical(level, "groups") || length(units) == n_rows &&
                 is.null(x$group_index)) units else units[x$group_index]
  # Classification statistics take the informative rows and groups only;
  # their weights are the same subset.
  if (length(weights) != n_rows) {
    weights <- weights[if (identical(level, "groups"))
      .multilpa_informative_groups(x) else .multilpa_informative_rows(x)]
  }
  stopifnot("every posterior row must carry its unit's weight" =
              length(weights) == n_rows && !anyNA(weights))
  weights
}

#' Mean posterior of every class within each modal assignment
#'
#' Documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A base `data.frame`, one row per level and ordered pair of classes.
#' @noRd
.multilpa_average_posteriors <- function(x, level = c("individuals", "groups", "both")) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(x))
  level <- .latents_match_arg(level, c("individuals", "groups", "both"), "level")
  posteriors <- .multilpa_posterior_levels(x, level)
  result <- do.call(rbind, lapply(names(posteriors), function(which_level) {
    .multilpa_average_posterior_rows(
      posteriors[[which_level]], which_level,
      .multilpa_level_weights(x, which_level, nrow(posteriors[[which_level]])))
  }))
  row.names(result) <- NULL
  result
}

#' Mean posterior of every class within each modal assignment group
#'
#' Under sampling weights the mean is weighted, so the average posterior is
#' the population's, as the class shares it is compared with are (review
#' 2026-10-26 round 2, S-3).
#'
#' @param probabilities Posterior probability matrix.
#' @param weights `NULL`, or one sampling weight per row.
#' @return A square matrix, rows indexed by assigned class, columns by class.
#' @noRd
.multilpa_average_posterior_matrix <- function(probabilities, weights = NULL) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities),
            "`weights` must be NULL or one weight per row" =
              is.null(weights) || length(weights) == nrow(probabilities))
  n_classes <- ncol(probabilities)
  assigned <- max.col(probabilities, ties.method = "first")
  row_weights <- weights %||% rep(1, nrow(probabilities))
  t(vapply(seq_len(n_classes), function(class) {
    rows <- assigned == class
    total <- sum(row_weights[rows])
    if (!any(rows) || total <= 0) rep(NA_real_, n_classes) else
      colSums(probabilities[rows, , drop = FALSE] * row_weights[rows]) / total
  }, numeric(n_classes)))
}

#' One row per class of classification diagnostics for one level
#'
#' Modal counts are counts of the rows supplied. Under sampling weights every
#' other column estimates the population: the class sizes and shares are the
#' weighted posterior totals and the average posteriors are weighted means, so
#' the odds of correct classification compare the two on the same population
#' (review 2026-10-26 round 2, S-3).
#'
#' @param probabilities Posterior probability matrix.
#' @param level Label recorded in the `level` column.
#' @param weights `NULL`, or one sampling weight per row.
#' @return A tidy `data.frame`, one row per class.
#' @noRd
.multilpa_classification_summary <- function(probabilities, level, weights = NULL) {
  average <- .multilpa_average_posterior_matrix(probabilities, weights)
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  n_modal <- tabulate(assigned, nbins = n_classes)
  estimated_n <- if (is.null(weights)) colSums(probabilities) else
    colSums(probabilities * weights)
  estimated_proportion <- estimated_n /
    if (is.null(weights)) nrow(probabilities) else sum(weights)
  diagonal <- average[cbind(classes, classes)]
  data.frame(level = level, class = .multilpa_class_labels(probabilities),
             n_modal = n_modal,
             proportion_modal = n_modal / nrow(probabilities),
             estimated_n = estimated_n,
             estimated_proportion = estimated_proportion,
             average_posterior = diagonal,
             odds_correct_classification =
               .multilpa_odds_correct(diagonal, estimated_proportion))
}

#' One row per assigned class and class of average posteriors for one level
#' @param probabilities Posterior probability matrix.
#' @param level Label recorded in the `level` column.
#' @param weights `NULL`, or one sampling weight per row; `n_assigned` stays a
#'   count.
#' @return A tidy `data.frame`, one row per assigned class and class.
#' @noRd
.multilpa_average_posterior_rows <- function(probabilities, level, weights = NULL) {
  average <- .multilpa_average_posterior_matrix(probabilities, weights)
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  labels <- .multilpa_class_labels(probabilities)
  data.frame(level = level,
             assigned_class = rep(labels, each = n_classes),
             class = rep(labels, times = n_classes),
             n_assigned = rep(tabulate(assigned, nbins = n_classes),
                              each = n_classes),
             average_posterior = as.vector(t(average)))
}

#' The class numbers a posterior matrix reports
#'
#' Its columns in order, except that a final `noise` column --- a fit's noise
#' component --- is class 0, the number mclust gives it.
#'
#' @param probabilities Posterior probability matrix.
#' @return An integer vector, one label per column.
#' @noRd
.multilpa_class_labels <- function(probabilities) {
  labels <- seq_len(ncol(probabilities))
  names <- colnames(probabilities)
  if (!is.null(names) && identical(names[length(names)], "noise")) {
    labels[length(labels)] <- 0L
  }
  labels
}

#' Odds of correct classification
#' @param average Average posterior probability in the assigned class.
#' @param proportion Model-estimated class proportion.
#' @return Odds ratio; `Inf` for a class whose average posterior is one to
#'   within rounding (perfect classification), and `NA_real_` where the ratio
#'   is undefined (no modal members, or a class proportion of zero or one).
#' @noRd
.multilpa_odds_correct <- function(average, proportion) {
  stopifnot("`average` and `proportion` must have equal length" =
              length(average) == length(proportion))
  sized <- !is.na(average) & !is.na(proportion) & proportion > 0 & proportion < 1
  # Within rounding of one, 1 - average is floating-point noise, and dividing
  # by it would report odds such as 7e15 that measure only that noise.
  perfect <- sized & average >= 1 - sqrt(.Machine$double.eps)
  defined <- sized & !perfect & average > 0
  odds <- rep(NA_real_, length(average))
  odds[perfect] <- Inf
  odds[defined] <- (average[defined] / (1 - average[defined])) /
    (proportion[defined] / (1 - proportion[defined]))
  odds
}

#' Classification entropy at each level a fit has
#'
#' Documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @return A base `data.frame`, one row per level, with `level`, `n_classes`,
#'   `n_units`, `entropy_sum` and `relative_entropy`.
#' @noRd
.multilpa_entropy_table <- function(x) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(x))
  posteriors <- list(individuals = .multilpa_classified_posteriors(x),
                     groups = .multilpa_classified_group_posteriors(x))
  posteriors <- Filter(Negate(is.null), posteriors)
  result <- do.call(rbind, lapply(names(posteriors), function(level) {
    probabilities <- posteriors[[level]]
    # Under sampling weights both are the population's (round 2, S-4);
    # `n_units` stays a count.
    weights <- .multilpa_level_weights(x, level, nrow(probabilities))
    data.frame(level = level, n_classes = ncol(probabilities),
               n_units = nrow(probabilities),
               entropy_sum = .multilpa_entropy_sum(probabilities, weights),
               relative_entropy = .multilpa_relative_entropy(probabilities, weights))
  }))
  row.names(result) <- NULL
  result
}
