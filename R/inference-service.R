# The inference service: covariance matrices from scores, the delta method,
# and Wald tables, for every model in the package.
#
# Every engine supplies what is particular to its likelihood -- the per-unit
# analytic scores, the packing of its parameters, the map to the reported
# scale -- and asks this file for the rest. The arithmetic here is exactly
# the arithmetic each engine used before the consolidation, so the numbers
# did not move. Where engines genuinely differ, the difference is an
# argument naming the engine's variant rather than a second implementation:
#
#   observed information  "five_point" (transitions over occasions, the
#                         additive families) or "central" (regression and
#                         growth mixtures) differences of the analytic score,
#                         or the profile models' scaled optimHess();
#   step rule             h = step * max(1, |theta|) ("max1"),
#                         h = step * (1 + |theta|) ("plus1"), or h = step
#                         ("absolute");
#   positive definiteness  "min_le_eps_max" (min <= 1e-10 max),
#                         "ratio_le_eps" (min / max <= 1e-10) or
#                         "min_le_eps_maxabs" (min <= 1e-10 max |ev|, which
#                         warns and returns NA rather than refusing);
#   critical value        qnorm((1 + level) / 2) ("upper") or
#                         qnorm(1 - (1 - level) / 2) ("one_minus"), which
#                         differ in the last bit at some levels;
#   delta method          J V J' ("quadratic") or rowSums((J V) * J)
#                         ("rowsums").
#
# None of these is a defect to be unified: each is its engine's published
# numbers, and changing one would move them in the last bits.

#' Observed information by differences of an analytic score
#'
#' Differentiates the total score numerically, each column with its own
#' step, and returns minus the symmetrized Jacobian. The analytic score makes
#' the truncation error that of the first derivative, not of a second
#' difference of the likelihood.
#'
#' @param total_score A function of a packed vector returning the summed
#'   score.
#' @param theta The packed estimate (named).
#' @param step Relative step.
#' @param scheme `"central"` (two points) or `"five_point"` (fourth order).
#' @param h_rule `"max1"`, `"plus1"` or `"absolute"`.
#' @return A symmetric matrix named by `theta`.
#' @noRd
.inference_score_information <- function(total_score, theta, step,
                                         scheme = c("central", "five_point"),
                                         h_rule = c("plus1", "max1", "absolute")) {
  scheme <- match.arg(scheme)
  h_rule <- match.arg(h_rule)
  width <- function(value) {
    switch(h_rule, plus1 = step * (1 + abs(value)),
           max1 = step * max(1, abs(value)), absolute = step)
  }
  if (identical(scheme, "central")) {
    columns <- lapply(seq_along(theta), function(j) {
      h <- width(theta[j])
      up <- theta
      down <- theta
      up[j] <- up[j] + h
      down[j] <- down[j] - h
      (total_score(up) - total_score(down)) / (2 * h)
    })
    hessian <- do.call(cbind, columns)
    information <- -(hessian + t(hessian)) / 2
  } else {
    columns <- vapply(seq_along(theta), function(k) {
      h <- width(theta[k])
      at <- function(multiple) {
        point <- theta
        point[k] <- point[k] + multiple * h
        total_score(point)
      }
      (-at(2) + 8 * at(1) - 8 * at(-1) + at(-2)) / (12 * h)
    }, numeric(length(theta)))
    information <- -(columns + t(columns)) / 2
  }
  dimnames(information) <- list(names(theta), names(theta))
  information
}

#' Refuse, or flag, an information matrix that is not positive definite
#'
#' @param information A symmetric matrix.
#' @param rule `"min_le_eps_max"`: refuse when the smallest eigenvalue is at
#'   most 1e-10 times the largest; `"ratio_le_eps"`: refuse when their ratio
#'   is at most 1e-10 (or not finite), and when the smallest is not positive;
#'   `"min_le_eps_maxabs"`: flag when the smallest is at most 1e-10 times the
#'   largest absolute eigenvalue.
#' @param message The refusal's message; `"%.3g"` in it receives the ratio.
#' @return `TRUE` when the matrix may be inverted. A failed `"min_le_eps_max"`
#'   or `"ratio_le_eps"` check raises `latents_singular_information`; a failed
#'   `"min_le_eps_maxabs"` check returns `FALSE`, for the caller to warn.
#' @noRd
.inference_positive_definite <- function(information,
                                         rule = c("min_le_eps_max", "ratio_le_eps",
                                                  "min_le_eps_maxabs"),
                                         message = NULL) {
  rule <- match.arg(rule)
  eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  refuse <- function(text) {
    stop(errorCondition(text, class = "latents_singular_information", call = NULL))
  }
  switch(rule,
    min_le_eps_max = {
      if (!all(is.finite(eigenvalues)) ||
          min(eigenvalues) <= 1e-10 * max(eigenvalues)) {
        refuse(message %||% "The observed information is not positive definite.")
      }
      TRUE
    },
    ratio_le_eps = {
      condition_ratio <- min(eigenvalues) / max(eigenvalues)
      # A negative definite matrix has a positive ratio; its smallest
      # eigenvalue is what shows it.
      if (!is.finite(condition_ratio) || condition_ratio <= 1e-10 ||
          min(eigenvalues) <= 0) {
        refuse(sprintf(message %||% paste(
          "The observed information is not positive definite (eigenvalue",
          "ratio %.3g)."), condition_ratio))
      }
      TRUE
    },
    min_le_eps_maxabs = all(is.finite(eigenvalues)) &&
      min(eigenvalues) > 1e-10 * max(abs(eigenvalues)))
}

#' Inverse of a positive definite information matrix
#'
#' The Cholesky inverse, which is symmetric by construction, after the
#' engine's positive-definiteness check.
#'
#' @inheritParams .inference_positive_definite
#' @return The inverse, without dimnames (the caller names it).
#' @noRd
.inference_invert <- function(information, rule = "min_le_eps_max", message = NULL) {
  .inference_positive_definite(information, rule, message)
  chol2inv(chol(information))
}

#' Refuse a robust or OPG covariance the unit scores cannot support
#'
#' The cross-product of the unit scores has the rank of the scores; with fewer
#' independent units than parameters it is singular, and a sandwich built on it
#' is not a covariance.
#'
#' @param scores Units by parameters.
#' @param message The refusal; `%d` placeholders receive the number of
#'   units, the rank and the number of parameters, in that order, when
#'   `formatted` is `TRUE`.
#' @param formatted Whether `message` is an `sprintf()` format.
#' @param extra Further `sprintf()` arguments after the three counts.
#' @return `NULL`, invisibly; raises `latents_too_few_groups`.
#' @noRd
.inference_require_rank <- function(scores, message, formatted = FALSE,
                                    extra = list()) {
  rank <- qr(scores)$rank
  if (rank < ncol(scores)) {
    text <- if (formatted) do.call(sprintf, c(list(message, nrow(scores), rank,
                                                    ncol(scores)), extra)) else message
    stop(errorCondition(text, class = "latents_too_few_groups", call = NULL))
  }
  invisible(NULL)
}

#' Robust (sandwich) covariance
#'
#' `bread %*% meat %*% bread`, left to right, with the meat the cross-product
#' of the unit scores (clustered scores where units are clustered).
#'
#' @param bread The inverse information.
#' @param meat `crossprod()` of the unit scores.
#' @return The sandwich.
#' @noRd
.inference_sandwich <- function(bread, meat) bread %*% meat %*% bread

#' Delta-method covariance of derived quantities
#'
#' `J V J'`, the covariance of the derived quantities to first order, for the
#' engines that report the matrix rather than its diagonal alone.
#'
#' @param jacobian Derived quantities by packed parameters.
#' @param covariance Covariance of the packed parameters.
#' @return The covariance matrix of the derived quantities.
#' @noRd
.inference_delta_covariance <- function(jacobian, covariance) {
  jacobian %*% covariance %*% t(jacobian)
}

#' Delta-method standard errors of derived quantities
#'
#' @param jacobian Derived quantities by packed parameters.
#' @param covariance Covariance of the packed parameters.
#' @param form `"quadratic"`, the diagonal of `J V J'`, or `"rowsums"`,
#'   `rowSums((J V) * J)`; the same numbers up to the last bits, and each
#'   engine keeps the one it has always reported.
#' @return Standard errors, `sqrt(max(variance, 0))`.
#' @noRd
.inference_delta_se <- function(jacobian, covariance,
                                form = c("quadratic", "rowsums")) {
  form <- match.arg(form)
  variance <- switch(form,
    quadratic = diag(.inference_delta_covariance(jacobian, covariance)),
    rowsums = rowSums((jacobian %*% covariance) * jacobian))
  sqrt(pmax(variance, 0))
}

#' Numerical Jacobian of a transformation by central differences
#'
#' @param transform A function of the packed vector.
#' @param theta The packed estimate.
#' @param step Relative step: `h = step * (1 + |theta_j|)`.
#' @return A matrix, derived quantities by packed parameters.
#' @noRd
.inference_numeric_jacobian <- function(transform, theta, step = 1e-6) {
  estimate <- transform(theta)
  jacobian <- do.call(cbind, lapply(seq_along(theta), function(j) {
    h <- step * (1 + abs(theta[j]))
    up <- theta
    down <- theta
    up[j] <- up[j] + h
    down[j] <- down[j] - h
    (transform(up) - transform(down)) / (2 * h)
  }))
  matrix(jacobian, length(estimate))
}

#' Derivative of a probability vector with respect to its free logits
#'
#' A softmax `p` has derivative `diag(p) - p p'` with respect to all its
#' logits; the reference category's logit is fixed at zero, so its column is
#' dropped -- the last category's for the profile models, the first class's
#' for the additive families.
#'
#' @param probabilities A probability vector.
#' @param reference `"last"` or `"first"`: whose column is dropped.
#' @return A `length(p)` by `length(p) - 1` matrix.
#' @noRd
.inference_simplex_jacobian <- function(probabilities, reference = c("last", "first")) {
  reference <- match.arg(reference)
  block <- diag(probabilities, nrow = length(probabilities)) - tcrossprod(probabilities)
  if (identical(reference, "first")) block[, -1L, drop = FALSE] else
    block[, seq_len(length(probabilities) - 1L), drop = FALSE]
}

#' The two-sided normal critical value
#' @param level Confidence level.
#' @param form `"upper"`, `qnorm((1 + level) / 2)`, or `"one_minus"`,
#'   `qnorm(1 - (1 - level) / 2)`: equal in exact arithmetic, and each engine
#'   keeps the spelling its published intervals used.
#' @return A single number.
#' @noRd
.inference_critical <- function(level, form = c("one_minus", "upper")) {
  form <- match.arg(form)
  switch(form, upper = stats::qnorm((1 + level) / 2),
         one_minus = stats::qnorm(1 - (1 - level) / 2))
}

#' Wald columns from estimates and standard errors
#'
#' The symmetric interval, the z statistic and its two-sided p value.
#'
#' @param estimate Estimates.
#' @param std_error Their standard errors.
#' @param level Confidence level.
#' @param form The critical value's form, see `.inference_critical()`.
#' @return A data frame: `estimate`, `std_error`, `statistic`, `p_value`,
#'   `conf_low`, `conf_high`.
#' @noRd
.inference_wald <- function(estimate, std_error, level, form = "one_minus") {
  z <- .inference_critical(level, form)
  statistic <- estimate / std_error
  data.frame(estimate = estimate, std_error = std_error,
             statistic = statistic,
             p_value = 2 * stats::pnorm(-abs(statistic)),
             conf_low = estimate - z * std_error,
             conf_high = estimate + z * std_error,
             row.names = NULL)
}


#' Differentiate the observed information and check that it is usable
#'
#' Curvature is tested after scaling, so that indicators on different units do
#' not by themselves make an identified model look singular.
#'
#' @return A list with the `scaled` and `natural` Hessians, the scaled
#'   `condition_ratio`, and the `inverse` of the scaled Hessian.
#' @noRd
.multilpa_observed_hessian <- function(scaled_objective, scaled_gradient,
                                     parameter_scale, step) {
  n_parameters <- length(parameter_scale)
  scaled <- tryCatch(
    stats::optimHess(rep(0, n_parameters), scaled_objective, gr = scaled_gradient,
                     control = list(ndeps = rep(step, n_parameters))),
    error = function(error) {
      stop(sprintf("Observed Hessian failed: %s", conditionMessage(error)))
    })
  natural <- scaled / tcrossprod(parameter_scale)
  eigenvalues <- eigen(natural, symmetric = TRUE, only.values = TRUE)$values
  scaled_eigenvalues <- eigen(scaled, symmetric = TRUE, only.values = TRUE)$values
  condition_ratio <- min(scaled_eigenvalues) / max(scaled_eigenvalues)
  if (any(!is.finite(eigenvalues)) || !is.finite(condition_ratio) ||
      min(scaled_eigenvalues) <= 0 || condition_ratio <= 1e-10) {
    stop(errorCondition(
      "Observed information is not positive definite or is numerically singular; Wald inference is unavailable.",
      class = "latents_singular_information", call = NULL))
  }
  inverse <- tryCatch(solve(scaled), error = function(error) {
    stop(sprintf("Observed information could not be inverted: %s",
                 conditionMessage(error)))
  })
  list(scaled = scaled, natural = natural, condition_ratio = condition_ratio,
       inverse = inverse)
}

#' Invert the outer product of the group scores
#'
#' The outer-product-of-gradients (BHHH) estimate of the information: the
#' cross-product of the per-group scores, which estimates the information by the
#' variance of the score rather than by the curvature of the likelihood. It is
#' what `glca` reports.
#'
#' @param cross The (scaled) cross-product of the group scores.
#' @return Its inverse, or `latents_singular_information` when it has none.
#' @noRd
.multilpa_opg_inverse <- function(cross) {
  factor <- tryCatch(chol(cross), error = function(error) NULL)
  if (is.null(factor) ||
      min(diag(factor))^2 < 1e-12 * max(diag(factor))^2) {
    stop(errorCondition(paste(
      "The outer product of the group scores is singular; OPG inference is",
      "unavailable. `vcov_type = \"observed\"` uses the curvature instead."),
      class = "latents_singular_information", call = NULL))
  }
  chol2inv(factor)
}

#' Robust cross-product matrix from per-group scores
#' @param scores Groups-by-parameters score matrix.
#' @return The summed outer product of the group scores.
#' @noRd
.multilpa_cross_product <- function(scores) {
  stopifnot("`scores` must be a numeric matrix" =
              is.matrix(scores) && is.numeric(scores))
  if (any(!is.finite(scores))) {
    stop(errorCondition("Per-group scores are not finite; robust inference is unavailable.",
                        class = "latents_bad_scores", call = NULL))
  }
  ## A single row cannot support a sandwich: at the estimate the rows sum to
  ## zero, so one row is itself the zero vector up to rounding and the meat is
  ## numerical dust rather than a variance.
  .multilpa_require_clusters(nrow(scores), 1L, "A robust cross-product matrix")
  crossprod(scores)
}

#' Refuse a cluster-robust variance too few independent units can support
#'
#' The meat of a cluster-robust sandwich is the sum of one outer product per
#' independent unit, formed from that unit's summed estimating-equation
#' contributions. Those contributions sum to zero across units at the estimate,
#' because that sum *is* the stationarity condition the estimate solves. The
#' meat therefore has rank at most `n_units - 1`, whatever the data: with one
#' unit it is exactly zero, and with `n_units <= n_quantities` it is singular,
#' so some contrast among the reported quantities is handed a variance of
#' exactly zero. Neither is a small standard error; both are refused.
#'
#' The rule this enforces, `n_units > n_quantities`, is the same one
#' [parameter_inference()] applies to a robust model vcov through
#' `.multilpa_check_regularity()`, so every robust path in the package asks for
#' as many independent units as the covariance it reports has dimensions, plus
#' one for the constraint that the contributions sum to zero.
#'
#' @param n_units Number of independent units the contributions were summed
#'   within: groups for a clustered variance, observations for an unclustered
#'   one.
#' @param n_quantities Number of quantities whose covariance is being reported.
#' @param what A noun phrase naming that covariance, opening the message.
#' @param unit What one independent unit is called in the message.
#' @param alternative A sentence naming a legitimate non-clustered alternative,
#'   or `NULL` when there is none.
#' @return `invisible(TRUE)` when the rule holds. Otherwise raises
#'   `latents_too_few_groups`.
#' @references Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide
#'   to cluster-robust inference. *Journal of Human Resources*, 50, 317--372.
#' @noRd
.multilpa_require_clusters <- function(n_units, n_quantities, what,
                                       unit = "independent groups",
                                       alternative = NULL) {
  stopifnot(
    "`n_units` must be a single non-negative whole number" =
      length(n_units) == 1L && is.numeric(n_units) && is.finite(n_units) &&
      n_units >= 0,
    "`n_quantities` must be a single positive whole number" =
      length(n_quantities) == 1L && is.numeric(n_quantities) &&
      is.finite(n_quantities) && n_quantities >= 1,
    "`what` must be a single string" = is.character(what) && length(what) == 1L,
    "`unit` must be a single string" = is.character(unit) && length(unit) == 1L,
    "`alternative` must be a single string or NULL" =
      is.null(alternative) || (is.character(alternative) &&
                                 length(alternative) == 1L))
  if (n_units > n_quantities) return(invisible(TRUE))
  plural <- if (n_quantities == 1L) "y" else "ies"
  reason <- if (n_units <= 1L) {
    paste("The contributions of a single unit are the estimating equation the",
          "estimate solves, so they sum to exactly zero and the variance would",
          "be floating-point noise rather than an estimate.")
  } else {
    sprintf(paste("The contributions sum to zero at the estimate, so %d %s leave at",
                  "most %d independent contributions for %d quantit%s: the",
                  "covariance is singular and some contrast among them would be",
                  "given a variance of exactly zero."),
            n_units, unit, n_units - 1L, n_quantities, plural)
  }
  stop(errorCondition(
    paste(c(sprintf("%s needs more %s than the %d quantit%s it covers; this fit has %d.",
                    what, unit, n_quantities, plural, n_units),
            reason, alternative), collapse = " "),
    class = "latents_too_few_groups", call = NULL))
}

#' Wald intervals on the scale where each parameter is unbounded
#'
#' A symmetric interval around a probability of 0.98 with a standard error of
#' 0.015 runs past one, and one around a small variance runs below zero:
#' neither is a set of values the parameter can take. Each interval is formed
#' where the parameter is unbounded and mapped back -- the logit of a
#' probability, the log of a variance -- by the delta method, with the same
#' standard error. It stays inside the parameter's range and becomes
#' asymmetric near a bound, as the sampling distribution does. Every other
#' parameter keeps the symmetric interval.
#' @param estimate,standard_error Numeric vectors.
#' @param kind `"probability"`, `"positive"` or anything else, per row.
#' @param critical The normal quantile for the level.
#' @return A list with `low` and `high`.
#' @noRd
.multilpa_wald_bounds <- function(estimate, standard_error, kind, critical) {
  low <- estimate - critical * standard_error
  high <- estimate + critical * standard_error
  probability <- kind == "probability" & is.finite(standard_error) &
    estimate > 0 & estimate < 1
  if (any(probability)) {
    centre <- stats::qlogis(estimate[probability])
    spread <- critical * standard_error[probability] /
      (estimate[probability] * (1 - estimate[probability]))
    low[probability] <- stats::plogis(centre - spread)
    high[probability] <- stats::plogis(centre + spread)
  }
  positive <- kind == "positive" & is.finite(standard_error) & estimate > 0
  if (any(positive)) {
    spread <- critical * standard_error[positive] / estimate[positive]
    low[positive] <- estimate[positive] * exp(-spread)
    high[positive] <- estimate[positive] * exp(spread)
  }
  list(low = low, high = high)
}

#' Which scale a parameter's interval is formed on
#' @param parameter,term The tidy label columns.
#' @param continuous The continuous indicator names.
#' @return `"probability"`, `"positive"` or `"real"`, per row.
#' @noRd
.multilpa_interval_kind <- function(parameter, term, continuous) {
  on_diagonal <- parameter == "covariance" &
    term %in% paste(continuous, continuous, sep = ":")
  ifelse(parameter %in% c("probability", "response", "initial_probability",
                          "transition_probability"), "probability",
         ifelse(parameter %in% c("variance", "count_mean", "count_dispersion") |
                  on_diagonal,
                "positive", "real"))
}

#' Add the adjusted p-value column an inference table reports
#'
#' The family is the tests the table actually carries: a bounded parameter has
#' no test, is `NA` in `p_value`, and must not inflate the correction for the
#' parameters that do.
#'
#' @param result An inference table carrying `statistic` and `p_value`.
#' @param adjust One of `.multilpa_p_adjust_methods`.
#' @return `result` with `p_adjusted` inserted after `p_value` and the
#'   `adjust` attribute recorded.
#' @noRd
.multilpa_adjust_p <- function(result, adjust) {
  stopifnot("`result` must carry a `p_value` column" = "p_value" %in% names(result))
  result$p_adjusted <- stats::p.adjust(result$p_value, method = adjust)
  ordered <- c("level", "outcome", "term", "parameter", "estimate",
               "standard_error", "statistic", "p_value", "p_adjusted",
               "conf_low", "conf_high")
  stopifnot("the inference table must carry exactly the documented columns" =
              setequal(names(result), ordered))
  result <- result[, ordered, drop = FALSE]
  attr(result, "adjust") <- adjust
  result
}
