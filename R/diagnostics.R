#' Classification entropy sum
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @return Sum of `-p log p` over every cell, the classification entropy.
#' @noRd
.multilpa_entropy_sum <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  positive <- probabilities > 0
  -sum(probabilities[positive] * log(probabilities[positive]))
}

#' Relative classification entropy
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @return Entropy on the zero-to-one scale, or `NA_real_` for a single class.
#' @noRd
.multilpa_relative_entropy <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  if (ncol(probabilities) == 1L) return(NA_real_)
  1 - .multilpa_entropy_sum(probabilities) /
    (nrow(probabilities) * log(ncol(probabilities)))
}

#' The formula behind each criterion, keyed by its name
#' @param criterion Character vector of criterion names.
#' @return Character vector of formulas, the same length as `criterion`.
#' @noRd
.multilpa_criterion_definitions <- function(criterion) {
  definitions <- c(deviance = "-2L",
                   aic = "-2L + 2q",
                   kic = "-2L + 3(q + 1)",
                   bic = "-2L + q log(n)",
                   sabic = "-2L + q log((n + 2) / 24)",
                   caic = "-2L + q (log(n) + 1)",
                   awe = "-2(L - EN) + 2q (1.5 + log(n))",
                   icl = "-2L + q log(n) + 2 EN",
                   clc = "-2L + 2 EN")
  stopifnot("Every criterion must have a definition" =
              all(criterion %in% names(definitions)))
  unname(definitions[criterion])
}

#' Every likelihood-penalty criterion a fit supports
#'
#' The formulas, the two sample-size conventions, the entropy each criterion
#' penalizes and the references are documented on `?get_data`, which is where
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
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x),
            "`definitions` must be TRUE or FALSE" =
              isTRUE(definitions) || isFALSE(definitions))
  format <- match.arg(format)
  if (isTRUE(definitions) && identical(format, "wide")) {
    stop(errorCondition(
      "`definitions` describes one criterion per row, which the wide form has not. Use format = \"long\".",
      class = "multilpa_bad_argument", call = NULL))
  }
  q <- x$n_parameters
  deviance <- -2 * x$log_likelihood
  conventions <- data.frame(
    convention = c("groups", "individuals"),
    n = c(x$n_groups, x$n_informative %||% x$n_observations))
  entropies <- c(groups = if (is.null(x$group_posteriors)) NA_real_ else
                   .multilpa_entropy_sum(x$group_posteriors),
                 individuals = .multilpa_entropy_sum(x$subject_posteriors))
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
                   # The only criterion here that does not use the sample size.
                   # Its convention still matters, because it decides which
                   # level's classification uncertainty is being penalized.
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
    result$definition <- .multilpa_criterion_definitions(result$criterion)
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
      class = "multilpa_unknown_criterion", call = NULL))
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
  posteriors <- list(individuals = object$subject_posteriors,
                     groups = object$group_posteriors)
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
                        class = "multilpa_no_group_classes", call = NULL))
  }
  "individuals"
}

#' Per-class classification quality at one or both levels
#'
#' The columns, the odds-of-correct-classification formula and the reference
#' are documented on `?get_data`, which is where a caller reaches this table
#' from.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A base `data.frame`, one row per level and class.
#' @noRd
.multilpa_classification_table <- function(x, level = c("individuals",
                                                        "groups", "both")) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(x))
  level <- match.arg(level)
  posteriors <- .multilpa_posterior_levels(x, level)
  result <- do.call(rbind, lapply(names(posteriors), function(which_level) {
    .multilpa_classification_summary(posteriors[[which_level]], which_level)
  }))
  row.names(result) <- NULL
  result
}

#' Mean posterior of every class within each modal assignment
#'
#' Documented on `?get_data`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A base `data.frame`, one row per level and ordered pair of classes.
#' @noRd
.multilpa_average_posteriors <- function(x, level = c("individuals", "groups", "both")) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(x))
  level <- match.arg(level)
  posteriors <- .multilpa_posterior_levels(x, level)
  result <- do.call(rbind, lapply(names(posteriors), function(which_level) {
    .multilpa_average_posterior_rows(posteriors[[which_level]], which_level)
  }))
  row.names(result) <- NULL
  result
}

#' Mean posterior of every class within each modal assignment group
#' @param probabilities Posterior probability matrix.
#' @return A square matrix, rows indexed by assigned class, columns by class.
#' @noRd
.multilpa_average_posterior_matrix <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  n_classes <- ncol(probabilities)
  assigned <- max.col(probabilities, ties.method = "first")
  t(vapply(seq_len(n_classes), function(class) {
    rows <- assigned == class
    if (!any(rows)) rep(NA_real_, n_classes) else
      colMeans(probabilities[rows, , drop = FALSE])
  }, numeric(n_classes)))
}

#' One row per class of classification diagnostics for one level
#' @param probabilities Posterior probability matrix.
#' @param level Label recorded in the `level` column.
#' @return A tidy `data.frame`, one row per class.
#' @noRd
.multilpa_classification_summary <- function(probabilities, level) {
  average <- .multilpa_average_posterior_matrix(probabilities)
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  n_modal <- tabulate(assigned, nbins = n_classes)
  estimated_n <- colSums(probabilities)
  estimated_proportion <- estimated_n / nrow(probabilities)
  diagonal <- average[cbind(classes, classes)]
  data.frame(level = level, class = classes, n_modal = n_modal,
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
#' @return A tidy `data.frame`, one row per assigned class and class.
#' @noRd
.multilpa_average_posterior_rows <- function(probabilities, level) {
  average <- .multilpa_average_posterior_matrix(probabilities)
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  data.frame(level = level,
             assigned_class = rep(classes, each = n_classes),
             class = rep(classes, times = n_classes),
             n_assigned = rep(tabulate(assigned, nbins = n_classes),
                              each = n_classes),
             average_posterior = as.vector(t(average)))
}

#' Odds of correct classification
#' @param average Average posterior probability in the assigned class.
#' @param proportion Model-estimated class proportion.
#' @return Odds ratio, or `NA_real_` where the ratio is undefined.
#' @noRd
.multilpa_odds_correct <- function(average, proportion) {
  stopifnot("`average` and `proportion` must have equal length" =
              length(average) == length(proportion))
  defined <- !is.na(average) & average > 0 & average < 1 &
    proportion > 0 & proportion < 1
  odds <- rep(NA_real_, length(average))
  odds[defined] <- (average[defined] / (1 - average[defined])) /
    (proportion[defined] / (1 - proportion[defined]))
  odds
}

#' Classification entropy at each level a fit has
#'
#' Documented on `?get_data`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @return A base `data.frame`, one row per level, with `level`, `n_classes`,
#'   `n_units`, `entropy_sum` and `relative_entropy`.
#' @noRd
.multilpa_entropy_table <- function(x) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(x))
  posteriors <- list(individuals = x$subject_posteriors,
                     groups = x$group_posteriors)
  posteriors <- Filter(Negate(is.null), posteriors)
  result <- do.call(rbind, lapply(names(posteriors), function(level) {
    probabilities <- posteriors[[level]]
    data.frame(level = level, n_classes = ncol(probabilities),
               n_units = nrow(probabilities),
               entropy_sum = .multilpa_entropy_sum(probabilities),
               relative_entropy = .multilpa_relative_entropy(probabilities))
  }))
  row.names(result) <- NULL
  result
}
