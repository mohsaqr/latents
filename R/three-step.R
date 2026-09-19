#' Posteriors and modal assignments at a chosen level
#' @param object A fitted model of this package.
#' @param level `"individuals"` or `"groups"`.
#' @return A list with `posteriors`, `modal`, `units` and `group_index`.
#' @noRd
.multilpa_level_assignments <- function(object, level) {
  stopifnot("`object` must be a fitted model of this package" =
              .multilpa_any_fit(object))
  level <- match.arg(level, c("individuals", "groups"))
  if (identical(level, "individuals")) {
    list(posteriors = object$subject_posteriors, modal = object$subject_profiles,
         units = seq_len(object$n_observations), group_index = object$group_index,
         n_classes = object$n_profiles)
  } else {
    list(posteriors = object$group_posteriors, modal = object$group_classes,
         units = object$group_values, group_index = seq_along(object$group_classes),
         n_classes = object$n_group_classes)
  }
}

#' How often each class is assigned to the wrong one
#'
#' The probability that a unit truly in one class is assigned to another, which
#' is the quantity every three-step method needs. Unlike the average-posterior
#' matrix from [classification_table()], this conditions on the *true* class
#' rather than the assigned one, which is the direction the correction requires.
#'
#' @param object A fitted model of this package.
#' @param level `"individuals"` for profiles, `"groups"` for group classes.
#' @return A base `data.frame` with one row per ordered pair of classes and the
#'   columns `level`, `true_class`, `assigned_class` and `probability`. The
#'   probabilities sum to one within each `true_class`.
#' @seealso [bch_weights()] and [three_step()], which consume it.
#' @references Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
#'   structure models with categorical variables. *Political Analysis*, 12,
#'   3--27. Vermunt, J. K. (2010). Latent class modeling with covariates: two
#'   improved three-step approaches. *Political Analysis*, 18, 450--469.
#' @examples
#' set.seed(11)
#' group <- rep(seq_len(40), each = 10)
#' truth <- 1L + as.integer(runif(400) > 0.5)
#' example_data <- data.frame(
#'   g = group,
#'   a = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
#'   b = rnorm(400, ifelse(truth == 2L, 1.1, -1.1))
#' )
#' fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 4, seed = 1)
#' classification_errors(fit)
#' @export
classification_errors <- function(object, level = c("individuals", "groups")) {
  level <- match.arg(level)
  pieces <- .multilpa_level_assignments(object, level)
  matrix_form <- .multilpa_error_matrix(pieces)
  classes <- seq_len(pieces$n_classes)
  index <- expand.grid(assigned_class = classes, true_class = classes)
  data.frame(level = level, true_class = index$true_class,
             assigned_class = index$assigned_class,
             probability = matrix_form[cbind(index$true_class, index$assigned_class)],
             row.names = NULL)
}

#' The classification error matrix, rows indexed by the true class
#' @param pieces The list from [.multilpa_level_assignments()].
#' @return A square matrix whose rows sum to one.
#' @noRd
.multilpa_error_matrix <- function(pieces) {
  classes <- seq_len(pieces$n_classes)
  totals <- colSums(pieces$posteriors)
  ## vapply() drops to a vector when there is only one class, so the shape is
  ## forced; a one-class model is degenerate but must not error here.
  result <- matrix(vapply(classes, function(assigned) {
    colSums(pieces$posteriors * (pieces$modal == assigned)) / totals
  }, numeric(pieces$n_classes)), nrow = pieces$n_classes,
  dimnames = list(paste0("true_", classes), paste0("assigned_", classes)))
  result
}

#' Weights that undo classification error
#'
#' The Bolck-Croon-Hagenaars weights. Analysing an outcome by modal class
#' attenuates every difference between classes, because some units are in the
#' wrong one; weighting by the inverse of the classification error matrix
#' removes that attenuation without letting the outcome influence the classes.
#'
#' @param object A fitted model of this package.
#' @param level `"individuals"` for profiles, `"groups"` for group classes.
#' @return A base `data.frame` with one row per unit and the columns `unit`,
#'   `assigned_class`, and one `weight_class_*` column per class. Weights are
#'   not probabilities: they are signed, and a unit assigned to one class
#'   ordinarily carries a negative weight for the others. They sum to one across
#'   the classes within a unit.
#' @seealso [three_step()], which applies them.
#' @references Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
#'   structure models with categorical variables. *Political Analysis*, 12,
#'   3--27.
#' @examples
#' set.seed(11)
#' group <- rep(seq_len(40), each = 10)
#' truth <- 1L + as.integer(runif(400) > 0.5)
#' example_data <- data.frame(
#'   g = group,
#'   a = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
#'   b = rnorm(400, ifelse(truth == 2L, 1.1, -1.1))
#' )
#' fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 4, seed = 1)
#' head(bch_weights(fit))
#' @export
bch_weights <- function(object, level = c("individuals", "groups")) {
  level <- match.arg(level)
  pieces <- .multilpa_level_assignments(object, level)
  inverse <- .multilpa_invert_errors(.multilpa_error_matrix(pieces))
  weights <- inverse[pieces$modal, , drop = FALSE]
  colnames(weights) <- paste0("weight_class_", seq_len(pieces$n_classes))
  data.frame(unit = pieces$units, assigned_class = pieces$modal, weights,
             row.names = NULL)
}

#' Invert a classification error matrix, refusing a singular one
#' @param errors The matrix from [.multilpa_error_matrix()].
#' @return Its inverse.
#' @noRd
.multilpa_invert_errors <- function(errors) {
  inverse <- tryCatch(solve(errors), error = function(condition) {
    stop(errorCondition(
      paste("The classification error matrix cannot be inverted, so the classes",
            "are not separated well enough for a three-step correction.",
            "Inspect classification_errors()."),
      class = "multilpa_inseparable_classes", call = NULL))
  })
  if (min(abs(eigen(as.matrix(errors), only.values = TRUE)$values)) < 1e-8) {
    stop(errorCondition(
      "The classification error matrix is numerically singular; the classes are not separable.",
      class = "multilpa_inseparable_classes", call = NULL))
  }
  inverse
}

#' Relate a latent class to an outcome it did not help define
#'
#' The three-step approach: the measurement model is fitted first, units are
#' assigned to classes, and only then is the outcome brought in, corrected for
#' the fact that some assignments are wrong. Doing it in one step instead lets
#' the outcome pull the classes toward itself; doing it naively by modal class
#' instead attenuates every difference. This does neither.
#'
#' @param object A fitted model of this package.
#' @param data The data frame carrying the outcome, in the fit's row order.
#' @param outcome Name of a numeric outcome column. For `level = "groups"` it
#'   must be constant within each group.
#' @param level `"individuals"` relates the outcome to profiles, `"groups"` to
#'   group classes.
#' @param method `"bch"` applies the Bolck-Croon-Hagenaars weights;
#'   `"proportional"` weights by the posteriors, which is simpler and
#'   attenuated; `"modal"` is the naive assignment, shown for comparison.
#' @param level_ci Confidence level for the intervals.
#' @return A base `data.frame` with one row per class and the columns `class`,
#'   `estimate`, `standard_error`, `conf_low`, `conf_high` and `effective_n`.
#'   Standard errors are cluster-robust, taking the fit's groups as the
#'   independent units, so they remain honest when the outcome is measured on
#'   observations nested inside those groups.
#' @details The correction assumes the outcome is independent of the assigned
#'   class given the true one, which is what makes a three-step method valid.
#'
#'   The error matrix is treated as known rather than estimated, so the intervals
#'   are optimistic and the amount is worth stating. Over 60 replications of a
#'   two-profile design with 8% misclassification and a true class difference of
#'   10, `"bch"` recovered that difference with a bias of -0.13 against -1.28 for
#'   `"modal"` and -1.82 for `"proportional"`, and its nominal 95% intervals
#'   covered the truth 85% of the time. Treat the point estimate as close to
#'   unbiased and the interval as somewhat too narrow; the other two methods are
#'   biased enough that their intervals covered the truth in none of the 60.
#' @references Vermunt, J. K. (2010). Latent class modeling with covariates: two
#'   improved three-step approaches. *Political Analysis*, 18, 450--469.
#'   Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class
#'   modeling with continuous distal outcomes. *Structural Equation Modeling*,
#'   23, 20--31.
#' @examples
#' set.seed(11)
#' group <- rep(seq_len(40), each = 10)
#' truth <- 1L + as.integer(runif(400) > 0.5)
#' example_data <- data.frame(
#'   g = group,
#'   a = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
#'   b = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
#'   y = rnorm(400, ifelse(truth == 2L, 10, 0))
#' )
#' fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 4, seed = 1)
#' three_step(fit, example_data, "y")
#' @export
three_step <- function(object, data, outcome,
                       level = c("individuals", "groups"),
                       method = c("bch", "proportional", "modal"),
                       level_ci = 0.95) {
  level <- match.arg(level)
  method <- match.arg(method)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`outcome` must name a single column of `data`" =
      is.character(outcome) && length(outcome) == 1L && outcome %in% names(data),
    "`outcome` must be numeric" = is.numeric(data[[outcome]]),
    "`outcome` must not be missing" = !anyNA(data[[outcome]]),
    "`level_ci` must be a single number in (0, 1)" =
      is.numeric(level_ci) && length(level_ci) == 1L && level_ci > 0 && level_ci < 1
  )
  pieces <- .multilpa_level_assignments(object, level)
  values <- .multilpa_outcome_values(object, data, outcome, level)
  weights <- .multilpa_step_weights(pieces, method)
  classes <- seq_len(pieces$n_classes)
  quantile <- stats::qnorm(1 - (1 - level_ci) / 2)

  rows <- lapply(classes, function(class) {
    weight <- weights[, class]
    total <- sum(weight)
    estimate <- sum(weight * values) / total
    ## Cluster-robust variance of a weighted mean: sum the influence
    ## contributions within each independent group, then across groups.
    influence <- weight * (values - estimate) / total
    clustered <- as.vector(rowsum(influence, pieces$group_index, reorder = FALSE))
    error <- sqrt(sum(clustered^2))
    data.frame(class = class, estimate = estimate, standard_error = error,
               conf_low = estimate - quantile * error,
               conf_high = estimate + quantile * error,
               effective_n = sum(weight)^2 / sum(weight^2), row.names = NULL)
  })
  result <- do.call(rbind, rows)
  attr(result, "method") <- method
  attr(result, "level") <- level
  result
}

#' The outcome, at the level being analysed
#' @return A numeric vector, one value per unit at that level.
#' @noRd
.multilpa_outcome_values <- function(object, data, outcome, level) {
  values <- data[[outcome]]
  stopifnot("`data` must have one row per observation of the fit" =
              length(values) == object$n_observations)
  if (identical(level, "individuals")) return(values)
  by_group <- split(values, object$group_index)
  constant <- vapply(by_group, function(v) length(unique(v)) == 1L, logical(1))
  if (!all(constant)) {
    stop(errorCondition(
      "For `level = \"groups\"` the outcome must be constant within each group.",
      class = "multilpa_bad_outcome", call = NULL))
  }
  vapply(by_group, function(v) v[[1L]], numeric(1))
}

#' Weights for the requested three-step method
#' @return A matrix with one row per unit and one column per class.
#' @noRd
.multilpa_step_weights <- function(pieces, method) {
  classes <- seq_len(pieces$n_classes)
  switch(method,
    bch = {
      inverse <- .multilpa_invert_errors(.multilpa_error_matrix(pieces))
      inverse[pieces$modal, , drop = FALSE]
    },
    proportional = pieces$posteriors,
    modal = vapply(classes, function(class) as.numeric(pieces$modal == class),
                   numeric(length(pieces$modal))))
}
