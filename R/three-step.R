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

#' Class-membership priors from a multinomial logit
#' @param gamma Coefficient matrix, predictors by non-reference classes.
#' @param design The model matrix.
#' @return A matrix of class probabilities, one row per unit.
#' @noRd
.multilpa_logit_prior <- function(gamma, design) {
  eta <- cbind(design %*% gamma, 0)
  scaled <- exp(eta - apply(eta, 1L, max))
  scaled / rowSums(scaled)
}

#' Covariates predicting class membership, corrected for misclassification
#'
#' The R3STEP approach. A multinomial logit of class membership on covariates,
#' fitted after the measurement model rather than alongside it, but treating the
#' assigned class as an error-prone indicator of the true one with the error
#' rates held fixed at what step one found. Regressing the modal class directly
#' instead attenuates every coefficient, because some units are in the wrong
#' class and the covariate cannot explain why.
#'
#' @param object A fitted model of this package.
#' @param data The data frame carrying the covariates, in the fit's row order.
#' @param covariates Character vector of numeric covariate columns. For
#'   `level = "groups"` each must be constant within a group.
#' @param level `"individuals"` predicts profile membership, `"groups"`
#'   predicts group-class membership.
#' @param level_ci Confidence level for the intervals.
#' @param vcov_type `"observed"` uses the observed information; `"robust"` uses
#'   the sandwich clustered on the fit's groups, which is the honest choice when
#'   the covariates are measured on observations nested inside them.
#' @return A base `data.frame` with one row per non-reference class and term,
#'   and the columns `level`, `outcome`, `term`, `estimate`, `standard_error`,
#'   `statistic`, `p_value`, `conf_low` and `conf_high`. Coefficients are log
#'   odds against the final class, which is the reference, matching
#'   [fit_covariates()].
#' @details The error matrix is held fixed rather than estimated jointly, which
#'   is what makes this a three-step method and what keeps the covariates from
#'   reshaping the classes.
#'
#'   Over 60 replications of a two-profile design with a true log-odds slope of
#'   1.2 and intercept -0.3, this recovered the slope with a bias of -0.002 and
#'   the intercept with a bias of +0.006, where a logistic regression on the
#'   modal class was biased by -0.163 on the slope. The correction costs some
#'   precision: the standard deviation of the slope across replications was
#'   0.127 against 0.100 for the naive fit. The nominal 95% intervals covered
#'   the true slope in 95% of those replications.
#'
#'   That coverage was measured where the profiles separate well. The error
#'   matrix is still treated as known rather than estimated, so intervals should
#'   be expected to run narrow where classification is poorer; check
#'   [classification_errors()] before relying on them.
#' @references Vermunt, J. K. (2010). Latent class modeling with covariates: two
#'   improved three-step approaches. *Political Analysis*, 18, 450--469.
#'   Asparouhov, T., & Muthen, B. (2014). Auxiliary variables in mixture
#'   modeling: three-step approaches using Mplus. *Structural Equation
#'   Modeling*, 21, 329--341.
#' @seealso [three_step()] for a distal outcome, and [fit_covariates()] for the
#'   one-step alternative that estimates everything jointly.
#' @examples
#' set.seed(21)
#' g <- rep(seq_len(50), each = 12)
#' x <- rnorm(600)
#' truth <- 1L + as.integer(runif(600) < plogis(-0.3 + 1.2 * x))
#' example_data <- data.frame(
#'   g = g, x = x,
#'   a = rnorm(600, ifelse(truth == 2L, 1.2, -1.2)),
#'   b = rnorm(600, ifelse(truth == 2L, 1.2, -1.2))
#' )
#' fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 4, seed = 1)
#' r3step(fit, example_data, "x")
#' @export
r3step <- function(object, data, covariates,
                   level = c("individuals", "groups"), level_ci = 0.95,
                   vcov_type = c("observed", "robust")) {
  level <- match.arg(level)
  vcov_type <- match.arg(vcov_type)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`covariates` must name columns of `data`" =
      is.character(covariates) && length(covariates) >= 1L &&
      all(covariates %in% names(data)),
    "`covariates` must be numeric" =
      all(vapply(data[covariates], is.numeric, logical(1))),
    "`covariates` must not be missing" = !anyNA(data[covariates]),
    "`level_ci` must be a single number in (0, 1)" =
      is.numeric(level_ci) && length(level_ci) == 1L && level_ci > 0 && level_ci < 1
  )
  pieces <- .multilpa_level_assignments(object, level)
  if (pieces$n_classes < 2L) {
    stop(errorCondition(
      "A single class has no membership to predict.",
      class = "multilpa_inseparable_classes", call = NULL))
  }
  design <- .multilpa_r3step_design(object, data, covariates, level)
  errors <- .multilpa_error_matrix(pieces)
  ## The likelihood of the assigned class, given the covariates and an error
  ## matrix fixed at its step-one value.
  by_assigned <- t(errors)[pieces$modal, , drop = FALSE]
  n_free <- pieces$n_classes - 1L
  shape <- function(theta) matrix(theta, ncol(design), n_free)

  objective <- function(theta) {
    prior <- .multilpa_logit_prior(shape(theta), design)
    -sum(log(pmax(rowSums(prior * by_assigned), 1e-300)))
  }
  score_rows <- function(theta) {
    prior <- .multilpa_logit_prior(shape(theta), design)
    joint <- prior * by_assigned
    posterior <- joint / rowSums(joint)
    residual <- (posterior - prior)[, seq_len(n_free), drop = FALSE]
    do.call(cbind, lapply(seq_len(n_free), function(k) design * residual[, k]))
  }
  gradient <- function(theta) -colSums(score_rows(theta))

  fitted <- stats::optim(rep(0, ncol(design) * n_free), objective, gradient,
                         method = "BFGS", control = list(maxit = 500L))
  if (!identical(fitted$convergence, 0L)) {
    stop(errorCondition(
      sprintf("The membership regression did not converge (optim code %d).",
              fitted$convergence),
      class = "multilpa_no_converge", call = NULL))
  }
  information <- .multilpa_observed_hessian(
    function(step) objective(fitted$par + step),
    function(step) gradient(fitted$par + step),
    rep(1, length(fitted$par)), 1e-4)
  covariance <- information$inverse
  if (identical(vcov_type, "robust")) {
    scores <- rowsum(score_rows(fitted$par), pieces$group_index, reorder = FALSE)
    covariance <- covariance %*% crossprod(scores) %*% covariance
  }
  .multilpa_r3step_frame(fitted$par, covariance, colnames(design), n_free,
                         pieces$n_classes, level, level_ci, vcov_type)
}

#' The covariate design at the level being analysed
#' @return A model matrix with an intercept.
#' @noRd
.multilpa_r3step_design <- function(object, data, covariates, level) {
  stopifnot("`data` must have one row per observation of the fit" =
              nrow(data) == object$n_observations)
  values <- data[covariates]
  if (identical(level, "groups")) {
    constant <- vapply(values, function(column) {
      all(vapply(split(column, object$group_index),
                 function(v) length(unique(v)) == 1L, logical(1)))
    }, logical(1))
    if (!all(constant)) {
      stop(errorCondition(
        "For `level = \"groups\"` every covariate must be constant within a group.",
        class = "multilpa_bad_outcome", call = NULL))
    }
    values <- values[!duplicated(object$group_index), , drop = FALSE]
  }
  design <- cbind(`(Intercept)` = 1, as.matrix(values))
  if (qr(design)$rank < ncol(design)) {
    stop(errorCondition(
      "The covariate design is rank deficient; remove constant or collinear predictors.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  design
}

#' Assemble the R3STEP coefficient table
#' @return One row per non-reference class and term.
#' @noRd
.multilpa_r3step_frame <- function(estimates, covariance, terms, n_free,
                                   n_classes, level, level_ci, vcov_type) {
  errors <- sqrt(pmax(diag(covariance), 0))
  statistic <- estimates / errors
  quantile <- stats::qnorm(1 - (1 - level_ci) / 2)
  labels <- expand.grid(term = terms, outcome = paste0("class_", seq_len(n_free)),
                        stringsAsFactors = FALSE)
  result <- data.frame(
    level = level, outcome = labels$outcome, term = labels$term,
    estimate = estimates, standard_error = errors, statistic = statistic,
    p_value = 2 * stats::pnorm(-abs(statistic)),
    conf_low = estimates - quantile * errors,
    conf_high = estimates + quantile * errors,
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "reference_class") <- n_classes
  attr(result, "vcov_type") <- vcov_type
  result
}
