#' Posteriors and modal assignments at a chosen level
#' @param object A fitted model of this package.
#' @param level `"individuals"` or `"groups"`.
#' @return A list with `posteriors`, `modal`, `units` and `group_index`.
#' @noRd
.multilpa_level_assignments <- function(object, level) {
  stopifnot("`object` must be a fitted model of this package" =
              .multilpa_any_fit(object))
  if (inherits(object, "multilpa_random_intercept")) {
    stop("Three-step methods require a discrete group-class model.")
  }
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
#' @param x A fitted model of this package.
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
classification_errors <- function(x, level = c("individuals", "groups")) {
  level <- match.arg(level)
  pieces <- .multilpa_level_assignments(x, level)
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
#' @param x A fitted model of this package.
#' @param level `"individuals"` for profiles, `"groups"` for group classes.
#' @return A base `data.frame` with one row per unit and class, and the columns
#'   `level`, `unit`, `assigned_class`, `class` and `weight`, ordered by unit and
#'   then by class. This is the expanded layout the method is defined on: one
#'   weighted record per unit per class. Weights are not probabilities: they are
#'   signed, and a unit assigned to one class ordinarily carries a negative
#'   weight for the others. Within a unit they sum to one across the classes.
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
bch_weights <- function(x, level = c("individuals", "groups")) {
  level <- match.arg(level)
  pieces <- .multilpa_level_assignments(x, level)
  inverse <- .multilpa_invert_errors(.multilpa_error_matrix(pieces))
  weights <- inverse[pieces$modal, , drop = FALSE]
  classes <- seq_len(pieces$n_classes)
  ## t() first, so that transposing the matrix into a vector walks a unit's
  ## classes before moving to the next unit and the frame reads unit by unit.
  data.frame(level = level,
             unit = rep(pieces$units, each = pieces$n_classes),
             assigned_class = rep(pieces$modal, each = pieces$n_classes),
             class = rep(classes, times = length(pieces$modal)),
             weight = as.vector(t(weights)),
             row.names = NULL, stringsAsFactors = FALSE)
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
#' @param x A fitted model of this package.
#' @param data The data frame carrying the outcome, in the fit's row order.
#' @param outcome Name of a numeric outcome column. For `level = "groups"` it
#'   must be constant within each group.
#' @param level `"individuals"` relates the outcome to profiles, `"groups"` to
#'   group classes.
#' @param method `"bch"` applies the Bolck-Croon-Hagenaars weights;
#'   `"proportional"` weights by the posteriors, which is simpler and
#'   attenuated; `"modal"` is the naive assignment, shown for comparison.
#' @param ci_level Confidence level for the intervals.
#' @param contrast `"none"` (the default) reports one outcome mean per class.
#'   `"pairs"` reports the difference between every pair of classes instead,
#'   which is the quantity a three-step analysis is usually run to test, with a
#'   statistic and a p-value against a difference of zero.
#' @param adjust Multiplicity correction applied across the pairwise
#'   differences, passed to [stats::p.adjust()]. `"BH"` by default; `"none"`
#'   leaves the p-values uncorrected. Ignored for `contrast = "none"`, which
#'   tests nothing.
#' @param vcov_type `"cluster"`, the default, sums the influence contributions
#'   within the fit's groups and takes those groups as the independent units,
#'   which is the honest choice when the outcome is measured on observations
#'   nested inside them. `"independent"` treats every observation as its own
#'   independent unit instead. That ignores the nesting the model was fitted to
#'   and is only defensible when there is no nesting left to ignore, such as a
#'   fit with a single group; it is never substituted silently, and the choice
#'   is recorded in the result's `vcov_type` attribute.
#' @return A base `data.frame`. For `contrast = "none"` it has one row per class
#'   and the columns `level`, `method`, `class`, `estimate`, `standard_error`,
#'   `conf_low`, `conf_high` and `effective_n`. It carries no `statistic` or
#'   `p_value`, unlike [r3step()]: an outcome mean has no meaningful null value
#'   to be tested against, and a p-value for "this class's mean is zero" would
#'   answer a question nobody asked. The comparison that does have a null is the
#'   difference between two classes, and that is what `contrast = "pairs"`
#'   returns: one row per pair with the columns `level`, `method`, `class`,
#'   `reference_class`, `estimate` (the mean of `class` minus the mean of
#'   `reference_class`), `standard_error`, `statistic`, `p_value`,
#'   `p_value_adjusted`, `conf_low` and `conf_high`, matching [r3step()].
#'
#'   Standard errors are cluster-robust in both shapes under the default
#'   `vcov_type = "cluster"`, taking the fit's groups as the independent units,
#'   so they remain honest when the outcome is measured on observations nested
#'   inside those groups; `vcov_type = "independent"` gives the unclustered
#'   variance instead. The variance that was used is recorded in the result's
#'   `vcov_type` attribute, and for `contrast = "pairs"` the correction applied
#'   to `p_value_adjusted` is recorded in its `adjust` attribute. A
#'   `contrast = "none"` table tests nothing, so it carries no `adjust`
#'   attribute.
#' @details The correction assumes the outcome is independent of the assigned
#'   class given the true one, which is what makes a three-step method valid.
#'
#'   A cluster-robust variance is the sum of one outer product per independent
#'   unit, and those contributions sum to zero at the estimate, so it has rank
#'   at most one less than the number of units. With a single group it is
#'   exactly zero and with as many groups as classes it is singular. Both are
#'   refused with `multilpa_too_few_groups` rather than reported as a very
#'   small standard error; `vcov_type = "independent"` is the labelled way to
#'   ask for the unclustered variance instead.
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
#' three_step(fit, example_data, "y", contrast = "pairs")
#' @export
three_step <- function(x, data, outcome,
                       level = c("individuals", "groups"),
                       method = c("bch", "proportional", "modal"),
                       ci_level = 0.95, contrast = c("none", "pairs"),
                       adjust = c("BH", "holm", "hochberg", "hommel",
                                    "bonferroni", "BY", "none"),
                       vcov_type = c("cluster", "independent")) {
  level <- match.arg(level)
  method <- match.arg(method)
  contrast <- match.arg(contrast)
  adjust <- match.arg(adjust)
  vcov_type <- match.arg(vcov_type)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`outcome` must name a single column of `data`" =
      is.character(outcome) && length(outcome) == 1L && outcome %in% names(data),
    "`outcome` must be numeric" = is.numeric(data[[outcome]]),
    "`outcome` must not be missing" = !anyNA(data[[outcome]]),
    "`outcome` must be finite" = all(is.finite(data[[outcome]])),
    "`ci_level` must be a single number in (0, 1)" =
      is.numeric(ci_level) && length(ci_level) == 1L && ci_level > 0 && ci_level < 1
  )
  pieces <- .multilpa_level_assignments(x, level)
  values <- .multilpa_outcome_values(x, data, outcome, level)
  weights <- .multilpa_step_weights(pieces, method)
  classes <- seq_len(pieces$n_classes)
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)

  totals <- colSums(weights)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop(errorCondition("A class has no positive total outcome weight.",
                        class = "multilpa_inseparable_classes", call = NULL))
  }
  estimates <- as.vector(crossprod(weights, values)) / totals
  ## Cluster-robust variance of a weighted mean: the per-unit influence
  ## contributions are summed within each independent group, then across
  ## groups. Holding them as a matrix lets a difference between two classes be
  ## given a standard error from the same quantities, rather than assuming the
  ## two means are independent when they share every unit.
  influence <- sweep(weights * outer(values, estimates, "-"), 2L, totals, "/")
  units <- if (identical(vcov_type, "cluster")) {
    pieces$group_index
  } else seq_len(nrow(influence))
  clustered <- rowsum(influence, units, reorder = FALSE)
  .multilpa_require_clusters(
    nrow(clustered), pieces$n_classes,
    what = "A cluster-robust standard error for the class outcome means",
    unit = if (identical(vcov_type, "cluster")) "independent groups" else "observations",
    alternative = if (identical(vcov_type, "cluster")) {
      paste("`vcov_type = \"independent\"` treats every observation as its own",
            "unit instead, which ignores the nesting and must be reported as",
            "having done so.")
    })

  if (identical(contrast, "pairs")) {
    return(.multilpa_step_pairs(estimates, clustered, level, method, quantile,
                                adjust, vcov_type))
  }
  errors <- sqrt(colSums(clustered^2))
  result <- data.frame(level = level, method = method, class = classes,
                       estimate = estimates, standard_error = errors,
                       conf_low = estimates - quantile * errors,
                       conf_high = estimates + quantile * errors,
                       effective_n = totals^2 / colSums(weights^2),
                       row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "vcov_type") <- vcov_type
  result
}

#' Pairwise differences between class outcome means
#'
#' @param estimates The per-class weighted means.
#' @param clustered Group-summed influence contributions, one column per class.
#' @param level,method The call's level and method, carried into the result.
#' @param quantile The normal quantile for the interval.
#' @param adjust A [stats::p.adjust()] method applied across the pairs.
#' @param vcov_type The variance requested, carried into the result's attribute.
#' @return One row per unordered pair of classes.
#' @noRd
.multilpa_step_pairs <- function(estimates, clustered, level, method, quantile,
                                 adjust, vcov_type) {
  n_classes <- length(estimates)
  if (n_classes < 2L) {
    stop(errorCondition(
      "A single class has no other class to be compared with.",
      class = "multilpa_inseparable_classes", call = NULL))
  }
  pairs <- utils::combn(n_classes, 2L)
  difference <- estimates[pairs[2L, ]] - estimates[pairs[1L, ]]
  errors <- sqrt(colSums((clustered[, pairs[2L, ], drop = FALSE] -
                            clustered[, pairs[1L, ], drop = FALSE])^2))
  statistic <- difference / errors
  raw <- 2 * stats::pnorm(-abs(statistic))
  result <- data.frame(
    level = level, method = method, class = pairs[2L, ],
    reference_class = pairs[1L, ], estimate = difference,
    standard_error = errors, statistic = statistic, p_value = raw,
    p_value_adjusted = stats::p.adjust(raw, method = adjust),
    conf_low = difference - quantile * errors,
    conf_high = difference + quantile * errors,
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "adjust") <- adjust
  attr(result, "vcov_type") <- vcov_type
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
#' @param x A fitted model of this package.
#' @param data The data frame carrying the covariates, in the fit's row order.
#' @param covariates Character vector of numeric covariate columns. For
#'   `level = "groups"` each must be constant within a group.
#' @param level `"individuals"` predicts profile membership, `"groups"`
#'   predicts group-class membership.
#' @param ci_level Confidence level for the intervals.
#' @param vcov_type `"observed"` uses the observed information; `"robust"` uses
#'   the sandwich clustered on the fit's groups, which is the honest choice when
#'   the covariates are measured on observations nested inside them. `"robust"`
#'   needs more independent groups than the regression has coefficients, because
#'   the group score contributions sum to zero at the estimate and so span at
#'   most one dimension fewer than there are groups; with too few it is refused
#'   with `multilpa_too_few_groups` rather than reporting a variance that is
#'   singular, or, with a single group, numerically zero. `"observed"` remains
#'   available there and does not allow for the nesting.
#' @param adjust Multiplicity correction applied across the covariate terms,
#'   passed to [stats::p.adjust()]. `"BH"` by default; `"none"` leaves the
#'   p-values uncorrected. The family is every covariate term of every
#'   non-reference class, which is the set of tests this call computes; the
#'   intercepts are not part of it.
#' @return A base `data.frame` with one row per non-reference class and term,
#'   and the columns `level`, `outcome`, `term`, `estimate`, `standard_error`,
#'   `statistic`, `p_value`, `p_value_adjusted`, `conf_low` and `conf_high`.
#'   `p_value` is uncorrected and `p_value_adjusted` carries the correction named
#'   by `adjust`, which is also recorded in the result's `adjust` attribute;
#'   it is `NA` on the intercept rows, which are not part of the tested family.
#'   Coefficients are log odds against the final class, which is the reference,
#'   matching [fit_covariates()]; that class is named in the result's
#'   `reference_class` attribute, and the variance that was used in its
#'   `vcov_type` attribute.
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
r3step <- function(x, data, covariates,
                   level = c("individuals", "groups"), ci_level = 0.95,
                   vcov_type = c("observed", "robust"),
                   adjust = c("BH", "holm", "hochberg", "hommel",
                                "bonferroni", "BY", "none")) {
  level <- match.arg(level)
  vcov_type <- match.arg(vcov_type)
  adjust <- match.arg(adjust)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`covariates` must name columns of `data`" =
      is.character(covariates) && length(covariates) >= 1L &&
      all(covariates %in% names(data)),
    "`covariates` must be numeric" =
      all(vapply(data[covariates], is.numeric, logical(1))),
    "`covariates` must not be missing" = !anyNA(data[covariates]),
    "`covariates` must be finite" =
      all(vapply(data[covariates], function(value) all(is.finite(value)), logical(1))),
    "`ci_level` must be a single number in (0, 1)" =
      is.numeric(ci_level) && length(ci_level) == 1L && ci_level > 0 && ci_level < 1
  )
  pieces <- .multilpa_level_assignments(x, level)
  if (pieces$n_classes < 2L) {
    stop(errorCondition(
      "A single class has no membership to predict.",
      class = "multilpa_inseparable_classes", call = NULL))
  }
  design <- .multilpa_r3step_design(x, data, covariates, level)
  errors <- .multilpa_error_matrix(pieces)
  ## The likelihood of the assigned class, given the covariates and an error
  ## matrix fixed at its step-one value.
  by_assigned <- t(errors)[pieces$modal, , drop = FALSE]
  n_free <- pieces$n_classes - 1L
  shape <- function(theta) matrix(theta, ncol(design), n_free)
  ## Refuse a sandwich the groups cannot support before paying for the fit.
  if (identical(vcov_type, "robust")) {
    .multilpa_require_clusters(
      length(unique(pieces$group_index)), ncol(design) * n_free,
      what = "A cluster-robust covariance for the membership regression",
      alternative = paste("`vcov_type = \"observed\"` uses the model-based",
                          "information instead, which does not allow for the",
                          "nesting."))
  }

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
    covariance <- covariance %*% .multilpa_cross_product(scores) %*% covariance
  }
  .multilpa_r3step_frame(fitted$par, covariance, colnames(design), n_free,
                         pieces$n_classes, level, ci_level, vcov_type, adjust)
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
                                   n_classes, level, ci_level, vcov_type,
                                   adjust) {
  errors <- sqrt(pmax(diag(covariance), 0))
  statistic <- estimates / errors
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  labels <- expand.grid(term = terms, outcome = paste0("class_", seq_len(n_free)),
                        stringsAsFactors = FALSE)
  raw <- 2 * stats::pnorm(-abs(statistic))
  ## The intercepts are estimated, not tested: correcting across them would
  ## enlarge the family with hypotheses nobody asked about and make every
  ## covariate look less significant than the data say.
  tested <- labels$term != "(Intercept)"
  adjusted <- rep(NA_real_, length(raw))
  adjusted[tested] <- stats::p.adjust(raw[tested], method = adjust)
  result <- data.frame(
    level = level, outcome = labels$outcome, term = labels$term,
    estimate = estimates, standard_error = errors, statistic = statistic,
    p_value = raw, p_value_adjusted = adjusted,
    conf_low = estimates - quantile * errors,
    conf_high = estimates + quantile * errors,
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "reference_class") <- n_classes
  attr(result, "vcov_type") <- vcov_type
  attr(result, "adjust") <- adjust
  result
}
