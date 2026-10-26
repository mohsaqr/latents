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
#' @param h_rule `"max1"`, `"plus1"`, `"absolute"` or `"scaled"`. The first two
#'   scale the step with the parameter's value, so they depend on its origin:
#'   an intercept of -9000, which a predictor shifted by 2000 produces, was
#'   differenced with a step of 0.09. `"scaled"` takes `step` times the
#'   coordinate's own `scale`, which a caller sets to a standard-error scale
#'   and which moves with any change of units or origin, but never more than
#'   `step * (1 + |theta|)`: a nearly unidentified coordinate has a huge
#'   standard-error scale, and a log-scale threshold gap at -8 with a scale of
#'   815 was otherwise differenced with a step of 0.8, across a visible part
#'   of its curvature (6% off a numerical Hessian).
#' @param scale For `h_rule = "scaled"`: one positive scale per coordinate;
#'   a non-finite or non-positive entry (a coordinate whose unit scores are
#'   all zero) falls back to the fixed width `1e-5 * (1 + |theta|)`, not to
#'   `step` times it: `step` is a fraction of a standard error and has no
#'   meaning without one.
#' @return A symmetric matrix named by `theta`.
#' @noRd
.inference_score_information <- function(total_score, theta, step,
                                         scheme = c("central", "five_point"),
                                         h_rule = c("plus1", "max1", "absolute",
                                                    "scaled"),
                                         scale = NULL) {
  scheme <- match.arg(scheme)
  h_rule <- match.arg(h_rule)
  stopifnot("`scale` must give one value per coordinate for h_rule = \"scaled\"" =
              !identical(h_rule, "scaled") ||
                (is.numeric(scale) && length(scale) == length(theta)))
  width_at <- function(j) {
    value <- theta[[j]]
    if (identical(h_rule, "scaled")) {
      usable <- is.finite(scale[[j]]) && scale[[j]] > 0
      return(if (usable) step * min(scale[[j]], 1 + abs(value)) else
        1e-5 * (1 + abs(value)))
    }
    switch(h_rule, plus1 = step * (1 + abs(value)),
           max1 = step * max(1, abs(value)), absolute = step)
  }
  if (identical(scheme, "central")) {
    columns <- lapply(seq_along(theta), function(j) {
      h <- width_at(j)
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
      h <- width_at(k)
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

#' Independent groups contributing to the fitted likelihood
#'
#' Empty measurement units do not supply independent estimating equations.
#' Preserve the explicit legacy empty-unit convention, but a unit with zero
#' sampling weight never contributes under either convention. Score matrices
#' keep all their rows for alignment with fitted groups; only validity counts
#' and the finite-sample correction use this selection.
#' @param object A fit with `n_groups`, `group_index` and optional `informative`
#'   and `sampling_weights` entries.
#' @param observed_only Require observed measurement even under the legacy
#'   empty-unit counting convention. Used for validity and rank checks; CR1
#'   retains the explicitly requested legacy convention.
#' @return A logical vector in fitted-group order.
#' @noRd
.multilpa_inference_groups <- function(object, observed_only = FALSE) {
  stopifnot(is.list(object), length(object$n_groups) == 1L,
            isTRUE(observed_only) || isFALSE(observed_only))
  selected <- rep(TRUE, object$n_groups)
  if ((observed_only || !isFALSE(object$drop_empty)) && !is.null(object$informative)) {
    selected[] <- FALSE
    selected[seq_len(object$n_groups) %in%
               object$group_index[object$informative]] <- TRUE
  }
  selected & (object$sampling_weights %||% rep(1, object$n_groups)) > 0
}

#' Small-sample factor of a cluster-robust variance
#'
#' `G / (G - 1)` for `G` independent units: the CR1 correction, which is the
#' default of `sandwich::vcovCL()` (Cameron and Miller, 2015). Without it the
#' variance (CR0) is biased downward when there are few units.
#'
#' @param n_units Number of independent units (groups, or rows when each row
#'   is its own unit).
#' @return A single number above 1.
#' @noRd
.inference_cluster_factor <- function(n_units) {
  stopifnot("`n_units` must be a single count of at least two units" =
              is.numeric(n_units) && length(n_units) == 1L && n_units >= 2)
  n_units / (n_units - 1)
}

#' Robust (sandwich) covariance
#'
#' `bread %*% meat %*% bread`, left to right, with the meat the cross-product
#' of the unit scores (clustered scores where units are clustered), scaled by
#' the CR1 factor `G / (G - 1)`.
#'
#' @param bread The inverse information.
#' @param meat `crossprod()` of the unit scores.
#' @param n_units Number of units whose scores form `meat`.
#' @return The sandwich.
#' @noRd
.inference_sandwich <- function(bread, meat, n_units) {
  bread %*% (meat * .inference_cluster_factor(n_units)) %*% bread
}

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
#' @param scale `NULL`, or one positive scale per coordinate (a standard
#'   error). With it, `h = scaled_step * min(scale_j, 1 + |theta_j|)`; a
#'   non-finite or non-positive entry uses the unscaled step.
#' @param scaled_step The fraction of the scale stepped. The truncation error
#'   of a central difference is about `h^2 f''' / 6`, the rounding error about
#'   `eps |f| / h`: at a hundred-thousandth of a standard error both sit near
#'   1e-11 for an exponential transform. A thousandth let truncation reach
#'   7e-7 (exp at a standard error of 2); a millionth let rounding reach 3e-9.
#' @return A matrix, derived quantities by packed parameters.
#' @noRd
.inference_numeric_jacobian <- function(transform, theta, step = 1e-6, scale = NULL,
                                        scaled_step = 1e-5) {
  stopifnot("`scale` must be NULL or one value per coordinate" =
              is.null(scale) || (is.numeric(scale) && length(scale) == length(theta)))
  estimate <- transform(theta)
  jacobian <- do.call(cbind, lapply(seq_along(theta), function(j) {
    # With a scale (a standard error), the step moves with the coordinate's
    # units and origin, capped as in .inference_score_information(); without
    # one, or where it is unusable, it is relative to 1 + |theta|.
    usable <- !is.null(scale) && is.finite(scale[[j]]) && scale[[j]] > 0
    h <- if (usable) scaled_step * min(scale[[j]], 1 + abs(theta[[j]])) else
      step * (1 + abs(theta[j]))
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


#' Results computed in this session, reused for identical inputs
#' @noRd
.latents_cache <- new.env(parent = emptyenv())

#' Reuse an expensive deterministic result whose inputs are identical
#'
#' Entries are compared with `identical()` on the whole key, so any change to
#' a fit, its data or an argument computes the result anew; the same fit
#' object passed again is recognised at once. Only values are stored: a
#' computation that fails is not, and fails again the next time. The most
#' recently used `size` entries are kept.
#' @param store The name of the store in `.latents_cache`.
#' @param key A list of everything the result depends on.
#' @param compute A function of no arguments that computes the result.
#' @param size The number of entries kept.
#' @return The result.
#' @noRd
.latents_memo <- function(store, key, compute, size = 8L) {
  entries <- .latents_cache[[store]] %||% list()
  hit <- Position(\(entry) identical(entry$key, key), entries)
  if (!is.na(hit)) {
    .latents_cache[[store]] <- c(entries[hit], entries[-hit])
    return(entries[[hit]]$value)
  }
  value <- compute()
  .latents_cache[[store]] <- utils::head(c(list(list(key = key, value = value)),
                                           entries), size)
  value
}

#' Differentiate the observed information and check that it is usable
#'
#' Curvature is tested after scaling, so that indicators on different units do
#' not by themselves make an identified model look singular. The Hessian
#' takes two score evaluations over all the data per parameter, which makes
#' it the slowest step of Wald inference; with a `cache_key` it is computed
#' once per session for those inputs, so a model plotted with intervals again
#' and again (or tabulated) does not differentiate its likelihood each time.
#'
#' @param cache_key `NULL`, or a list of everything the objective and the
#'   gradient depend on (see `.latents_memo()`).
#' @return A list with the `scaled` and `natural` Hessians, the scaled
#'   `condition_ratio`, and the `inverse` of the scaled Hessian.
#' @noRd
.multilpa_observed_hessian <- function(scaled_objective, scaled_gradient,
                                     parameter_scale, step, cache_key = NULL) {
  if (!is.null(cache_key)) {
    return(.latents_memo("observed_information",
                         c(cache_key, list(scale = parameter_scale, step = step)),
                         \() .multilpa_observed_hessian(
                           scaled_objective, scaled_gradient, parameter_scale, step)))
  }
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
#' @param continuous The continuous indicator names, as the data spell them;
#'   the diagonal term `a:a` is built by `.latents_term()`, the helper every
#'   label producer uses.
#' @return `"probability"`, `"positive"` or `"real"`, per row.
#' @noRd
.multilpa_interval_kind <- function(parameter, term, continuous) {
  on_diagonal <- parameter == "covariance" &
    term %in% .latents_term(continuous, continuous)
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

#' Covariance types a fit on person-transformed indicators supports
#'
#' Person centring (and person scaling) replaces each unit's rows by their
#' deviations from that unit's own mean (or divides them by its own standard
#' deviation). The rows a unit contributes then sum to zero and depend on the
#' unit's latent mix, so the likelihood the model maximizes treats as
#' independent draws values that are not, and the inverse observed
#' information understates the variance: in a Monte Carlo study (300 data
#' sets, 40 units of 6) its standard errors of the profile means were 0.78 to
#' 0.93 of the empirical spread (review 2026-10-26 round 2, INF2-2). Each
#' unit's transformed rows are a function of that unit's data alone, so the
#' cluster-robust sandwich over units is valid (0.91 to 1.01 in the same
#' study); the outer product of the scores is not a variance either. A
#' person-transformed fit therefore defaults to `"robust"` and refuses the
#' other types.
#'
#' @param object A fitted model with `centering` and `scaling`.
#' @param vcov_type The requested type, already matched.
#' @param defaulted Whether the caller left `vcov_type` at its default.
#' @return The type to use; raises `latents_unsupported_inference` for
#'   `"observed"` or `"opg"` asked of a person-transformed fit.
#' @noRd
.multilpa_transformed_vcov <- function(object, vcov_type, defaulted) {
  person <- identical(object$centering %||% "none", "person") ||
    identical(object$scaling %||% "none", "person")
  if (!person || identical(vcov_type, "robust")) return(vcov_type)
  if (defaulted) return("robust")
  stop(errorCondition(sprintf(paste(
    "`vcov_type = \"%s\"` treats a unit's person-centred or person-scaled rows",
    "as independent observations, which they are not, and understates the",
    "standard errors. Use `vcov_type = \"robust\"` (the default for these",
    "fits) or `method = \"bootstrap\"`."), vcov_type),
    class = "latents_unsupported_inference", call = NULL))
}

#' The grand transform of the indicators, as estimated quantities
#'
#' `centering = "grand"` subtracts the sample mean of each indicator and
#' `scaling = "grand"` divides by its sample standard deviation. Both are
#' estimates, computed from the same data the model is fitted to, so a
#' reported mean `mu_k - xbar` varies with `xbar` as well as with `mu_k`.
#' Each is written as an estimating equation summed over units, over the
#' observed cells of the units that carry weight (the rows that set them):
#' `sum_i (x_ij - m_j) = 0` for the mean and
#' `sum_i ((x_ij - m_j)^2 - s_j^2 (n_j - 1) / n_j) = 0` for the standard
#' deviation (denominator `n_j - 1`, as `.multilpa_indicator_scales()` uses).
#' A person transform is a function of each unit's own data and needs no
#' equation; it is held at its stored value.
#'
#' @param object A fitted `multilpa` model with stored `indicator_data`.
#' @return `NULL` when the fit has no grand transform; otherwise a list with
#'   `eta` (the means, then the standard deviations, of the indicators the
#'   transform estimated), `steps` (a differencing width per entry),
#'   `transform` (a function of `eta` returning the transformed indicator
#'   matrix), `contributions` (units by entries of `eta`, at the estimate) and
#'   `slope` (the diagonal of minus the derivative of the summed equations).
#' @noRd
.multilpa_grand_transform <- function(object) {
  centering <- object$centering %||% "none"
  scaling <- object$scaling %||% "none"
  grand_mean <- identical(centering, "grand") || identical(scaling, "grand")
  if (!grand_mean || is.null(object$indicator_data) ||
      ncol(object$indicator_data) == 0L) return(NULL)
  raw <- .multilpa_uncentered_indicators(object)
  n_groups <- object$n_groups
  carrying <- if (is.null(object$sampling_weights)) rep(TRUE, nrow(raw)) else
    (object$sampling_weights > 0)[object$group_index]
  observed <- !is.na(raw) & carrying
  counts <- colSums(observed)
  means <- if (identical(centering, "grand")) object$centering_offsets[1L, ] else
    colSums(ifelse(observed, raw, 0)) / counts
  deviations <- ifelse(observed, sweep(raw, 2L, means, "-"), 0)
  spread <- sqrt(colSums(deviations^2) / pmax(counts - 1, 1))
  scaled <- identical(scaling, "grand")
  n_columns <- ncol(raw)
  eta <- c(means, if (scaled) object$scaling_factors[1L, ])
  offsets_for <- function(eta) {
    switch(centering,
      grand = matrix(eta[seq_len(n_columns)], nrow(raw), n_columns, byrow = TRUE),
      person = object$centering_offsets[object$group_index, , drop = FALSE],
      none = matrix(0, nrow(raw), n_columns))
  }
  scales_for <- function(eta) {
    switch(scaling,
      grand = matrix(eta[n_columns + seq_len(n_columns)], nrow(raw), n_columns,
                     byrow = TRUE),
      person = object$scaling_factors[object$group_index, , drop = FALSE],
      none = matrix(1, nrow(raw), n_columns))
  }
  transform <- function(eta) {
    value <- (raw - offsets_for(eta)) / scales_for(eta)
    dimnames(value) <- dimnames(object$indicator_data)
    value
  }
  mean_terms <- rowsum(deviations, object$group_index, reorder = TRUE)
  contributions <- if (!scaled) mean_terms else {
    standard <- eta[n_columns + seq_len(n_columns)]
    square_terms <- deviations^2 -
      observed * matrix(standard^2 * (counts - 1) / counts, nrow(raw), n_columns,
                        byrow = TRUE)
    cbind(mean_terms, rowsum(square_terms, object$group_index, reorder = TRUE))
  }
  contributions <- matrix(contributions, n_groups, length(eta))
  list(eta = unname(eta),
       steps = 1e-4 * unname(c(spread, if (scaled) spread)),
       transform = transform, contributions = contributions,
       slope = unname(c(counts, if (scaled) 2 * eta[n_columns + seq_len(n_columns)] *
                                              (counts - 1))))
}

#' Add the variability of an estimated grand transform to a covariance
#'
#' The model's parameters and the grand transform are estimated jointly, as a
#' stacked system of estimating equations: the model's score, and the mean and
#' standard-deviation equations of `.multilpa_grand_transform()`. To first
#' order `theta_hat - theta = bread S + L E`, with `S` the summed unit scores,
#' `E` the summed transform equations, and `L = bread C diag(1 / slope)`, where
#' `C` is the derivative of the summed score with respect to the transform
#' (taken by central differences of the analytic score on re-transformed data).
#' The covariance of the model alone, `covariance`, is kept for the `S` term --
#' the inverse information, the sandwich or the OPG inverse, as asked -- and the
#' terms in `E` are added from the unit contributions:
#' `covariance + c (bread B_SE L' + L B_ES bread) + c L B_EE L'`, with
#' `c = G / (G - 1)` the CR1 factor when the model's covariance is the robust
#' sandwich, which carries it too, and `c = 1` otherwise (Stefanski and Boos,
#' 2002). The factor must match: a transform the estimates absorb exactly,
#' such as the mean and variance of one z-scored indicator in one profile,
#' then gets a variance of zero, as it should, rather than `(1 - c)` times its
#' model variance, which is negative. Without the added terms a grand-centred mean had the standard error
#' of the uncentred mean, 0.69 to 0.95 of its empirical spread in a Monte Carlo
#' study (review 2026-10-26 round 2, INF2-1).
#'
#' @param object A fitted `multilpa` model.
#' @param covariance The covariance of the free unconstrained coordinates.
#' @param bread The inverse information those coordinates were inverted with
#'   (the OPG inverse for `vcov_type = "opg"`).
#' @param free Positions of the free coordinates in the unconstrained vector.
#' @param codes The categorical codes of the fit.
#' @param robust Whether `covariance` is the CR1 sandwich.
#' @param score_of,group_scores_of `NULL` for a profile model, or, for a
#'   family with its own parameter vector (the covariate fits), functions of a
#'   transformed indicator matrix returning the summed score and the
#'   units-by-coordinates score matrix in the coordinates of `covariance`.
#' @return `covariance`, with the transform's contribution added when the fit
#'   has a grand transform. A variance that comes out negative beyond rounding
#'   (more than 1e-8 of the model's own) raises
#'   `latents_singular_information`; one within rounding of zero is zero.
#' @references Stefanski, L. A., & Boos, D. D. (2002). The calculus of
#'   M-estimation. *The American Statistician*, 56, 29--38.
#' @noRd
.multilpa_transform_covariance <- function(object, covariance, bread, free,
                                           codes = NULL, robust = FALSE,
                                           score_of = NULL, group_scores_of = NULL) {
  nuisance <- .multilpa_grand_transform(object)
  if (is.null(nuisance)) return(covariance)
  if (is.null(score_of) || is.null(group_scores_of)) {
    theta <- unname(.multilpa_coefficients(object, "unconstrained"))
  }
  score_of <- score_of %||% function(transformed) {
    -.multilpa_score(theta, transformed, object, codes)[free]
  }
  group_scores_of <- group_scores_of %||% function(transformed) {
    .multilpa_group_scores(theta, transformed, object, codes)[, free, drop = FALSE]
  }
  total_score <- function(eta) score_of(nuisance$transform(eta))
  derivative <- do.call(cbind, lapply(seq_along(nuisance$eta), function(j) {
    h <- nuisance$steps[[j]]
    up <- nuisance$eta
    down <- nuisance$eta
    up[j] <- up[j] + h
    down[j] <- down[j] - h
    (total_score(up) - total_score(down)) / (2 * h)
  }))
  scores <- group_scores_of(object$indicator_data)
  influence <- bread %*% sweep(matrix(derivative, nrow(bread)), 2L,
                               nuisance$slope, "/")
  factor <- if (isTRUE(robust)) {
    .inference_cluster_factor(sum(.multilpa_inference_groups(object)))
  } else 1
  cross <- bread %*% crossprod(scores, nuisance$contributions) %*% t(influence)
  added <- factor * (cross + t(cross) +
                       influence %*% crossprod(nuisance$contributions) %*% t(influence))
  result <- covariance + added
  result <- (result + t(result)) / 2
  # The terms added for the transform mix the model's information with
  # empirical moments; a variance negative beyond rounding would be read as a
  # zero standard error, so it is refused rather than floored. One within
  # rounding of zero is an estimate the transform fixes exactly.
  variances <- diag(result)
  rounding <- 1e-8 * pmax(abs(diag(covariance)), .Machine$double.eps)
  if (any(!is.finite(variances)) || any(variances < -rounding)) {
    stop(errorCondition(paste(
      "The covariance with the estimated centring or scaling added is not",
      "positive; use `vcov_type = \"robust\"` or `method = \"bootstrap\"`."),
      class = "latents_singular_information", call = NULL))
  }
  diag(result) <- pmax(variances, 0)
  dimnames(result) <- dimnames(covariance)
  result
}
