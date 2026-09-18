#' Parameterize a fitted mixture
#' @param object A fitted `multilpa` model.
#' @param scale Natural estimates or unconstrained fitting coordinates.
#' @return A named numeric parameter vector.
#' @noRd
.multilpa_coefficients <- function(object, scale) {
  stopifnot(inherits(object, "multilpa"), scale %in% c("natural", "unconstrained"))
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  n_indicators <- length(object$indicators)
  mean_names <- as.vector(t(outer(seq_len(n_profiles), object$indicators,
    function(profile, indicator) sprintf("mean[%s,%s]", profile, indicator))))
  variance_names <- if (object$variance_model == "equal") {
    sprintf("variance[shared,%s]", object$indicators)
  } else sub("^mean", "variance", mean_names)
  variance_values <- if (object$variance_model == "equal") object$variances[1L, ] else
    as.vector(t(object$variances))
  profile_indices <- if (scale == "natural") seq_len(n_profiles) else seq_len(n_profiles - 1L)
  group_indices <- if (scale == "natural") seq_len(n_types) else seq_len(n_types - 1L)
  profile_names <- as.vector(t(outer(seq_len(n_types), profile_indices,
    function(group, profile) sprintf("profile_probability[%s,%s]", group, profile))))
  group_names <- sprintf("group_probability[%s]", group_indices)
  profile_values <- as.vector(t(object$profile_probabilities[, profile_indices, drop = FALSE]))
  group_values <- object$group_probabilities[group_indices]
  if (scale == "unconstrained") {
    variance_values <- log(variance_values)
    variance_names <- sub("^variance", "log_variance", variance_names)
    profile_values <- as.vector(t(log(object$profile_probabilities[, profile_indices, drop = FALSE]) -
      log(object$profile_probabilities[, n_profiles])))
    group_values <- log(group_values) - log(object$group_probabilities[n_types])
    profile_names <- sub("^profile_probability", "profile_logit", profile_names)
    group_names <- sub("^group_probability", "group_logit", group_names)
  }
  if (identical(object$covariance_model, "full")) {
    covariance_coordinates <- .multilpa_covariance_coordinates(object, scale)
    variance_values <- unname(covariance_coordinates)
    variance_names <- names(covariance_coordinates)
  }
  stats::setNames(c(as.vector(t(object$means)), variance_values,
                   profile_values, group_values),
                 c(mean_names, variance_names, profile_names, group_names))
}

#' Encode unrestricted Gaussian covariance matrices
#' @param object A fitted full-covariance model.
#' @param scale Natural covariance or log-Cholesky coordinates.
#' @return A named vector of nonredundant covariance parameters.
#' @noRd
.multilpa_covariance_coordinates <- function(object, scale) {
  stopifnot(inherits(object, "multilpa"), identical(object$covariance_model, "full"),
            scale %in% c("natural", "unconstrained"))
  dimension <- length(object$indicators)
  profiles <- if (object$variance_model == "equal") 1L else seq_len(object$n_profiles)
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  positions <- which(lower, arr.ind = TRUE)
  unlist(lapply(profiles, function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance <- matrix(object$covariances[, , profile], dimension, dimension)
    values <- if (scale == "natural") covariance else t(chol(covariance))
    if (scale == "unconstrained") diag(values) <- log(diag(values))
    prefix <- if (scale == "natural") "covariance" else "cholesky"
    labels <- sprintf("%s[%s,%s,%s]", prefix,
      if (object$variance_model == "equal") "shared" else profile,
      object$indicators[positions[, 1L]], object$indicators[positions[, 2L]])
    if (scale == "unconstrained") {
      labels[positions[, 1L] == positions[, 2L]] <- sub("^cholesky", "log_cholesky",
        labels[positions[, 1L] == positions[, 2L]])
    }
    stats::setNames(values[lower], labels)
  }), use.names = TRUE)
}

#' Decode unconstrained mixture coordinates
#' @param theta Parameter vector with log variances and baseline-category logits.
#' @param object A fitted model defining the parameter dimensions.
#' @return The four parameter blocks used by the expectation step.
#' @noRd
.multilpa_decode <- function(theta, object) {
  stopifnot(is.numeric(theta), inherits(object, "multilpa"),
            length(theta) == object$n_parameters)
  theta <- unname(theta)
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  n_indicators <- length(object$indicators)
  n_means <- n_profiles * n_indicators
  n_covariance <- if (identical(object$covariance_model, "full"))
    n_indicators * (n_indicators + 1L) / 2L else n_indicators
  n_variances <- n_covariance * if (object$variance_model == "equal") 1L else n_profiles
  n_logits <- n_types * (n_profiles - 1L)
  means <- matrix(theta[seq_len(n_means)], n_profiles, byrow = TRUE)
  covariances <- NULL
  if (identical(object$covariance_model, "full")) {
    lower <- lower.tri(matrix(0, n_indicators, n_indicators), diag = TRUE)
    covariance_list <- lapply(seq_len(n_profiles), function(profile) {
      stopifnot(is.numeric(profile), length(profile) == 1L)
      block <- if (object$variance_model == "equal") 1L else profile
      factor <- matrix(0, n_indicators, n_indicators)
      factor[lower] <- theta[n_means + (block - 1L) * n_covariance + seq_len(n_covariance)]
      diag(factor) <- exp(diag(factor))
      tcrossprod(factor)
    })
    covariances <- array(unlist(covariance_list, use.names = FALSE),
      c(n_indicators, n_indicators, n_profiles))
    variances <- t(matrix(vapply(covariance_list, diag, numeric(n_indicators)),
      n_indicators, n_profiles))
  } else {
    variances <- matrix(exp(theta[n_means + seq_len(n_variances)]),
      n_profiles, n_indicators, byrow = TRUE)
  }
  profile_logits <- cbind(matrix(theta[n_means + n_variances + seq_len(n_logits)],
    n_types, n_profiles - 1L, byrow = TRUE), 0)
  profile_probabilities <- exp(sweep(profile_logits, 1L,
    .multilpa_log_sum_exp(profile_logits), "-"))
  group_logits <- c(theta[n_means + n_variances + n_logits + seq_len(n_types - 1L)], 0)
  group_probabilities <- exp(group_logits - max(group_logits))
  parameters <- list(means = means, variances = variances,
    profile_probabilities = profile_probabilities,
    group_probabilities = group_probabilities / sum(group_probabilities))
  if (!is.null(covariances)) parameters$covariances <- covariances
  parameters
}

#' Differentiate the observed mixture likelihood
#' @param theta Unconstrained parameter vector.
#' @param x Centered numeric indicator matrix, possibly containing NA.
#' @param object Fitted model defining dimensions and groups.
#' @return The gradient of the negative observed-data log likelihood.
#' @noRd
.multilpa_score <- function(theta, x, object) {
  stopifnot(is.numeric(theta), is.matrix(x), inherits(object, "multilpa"))
  parameters <- .multilpa_decode(theta, object)
  expectation <- .multilpa_expectation(x, object$group_index, parameters)
  if (identical(object$covariance_model, "full")) {
    measurement <- .multilpa_full_measurement_score(parameters, expectation, object)
    mean_score <- measurement$means
    variance_score <- measurement$covariances
  } else {
  observed <- !is.na(x)
  measurement_scores <- lapply(seq_len(object$n_profiles), function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
    residuals[!observed] <- 0
    weights <- expectation$subject_posteriors[, profile]
    list(means = colSums(sweep(residuals, 2L, parameters$variances[profile, ], "/") * weights),
      variances = 0.5 * colSums((sweep(residuals^2, 2L,
        parameters$variances[profile, ], "/") - observed) * weights))
  })
  mean_score <- unlist(lapply(measurement_scores, `[[`, "means"), use.names = FALSE)
  variance_score <- if (object$variance_model == "equal") {
    Reduce(`+`, lapply(measurement_scores, `[[`, "variances"))
  } else unlist(lapply(measurement_scores, `[[`, "variances"), use.names = FALSE)
  }
  profile_score <- unlist(lapply(seq_len(object$n_group_classes), function(group) {
    stopifnot(is.numeric(group), length(group) == 1L)
    counts <- colSums(expectation$joint[[group]])
    (counts - sum(counts) * parameters$profile_probabilities[group, ])[
      seq_len(object$n_profiles - 1L)]
  }), use.names = FALSE)
  group_score <- (colSums(expectation$group_posteriors) -
    object$n_groups * parameters$group_probabilities)[seq_len(object$n_group_classes - 1L)]
  -c(mean_score, variance_score, profile_score, group_score)
}

#' Gaussian scores in log-Cholesky coordinates
#' @param parameters Decoded Gaussian and mixing parameters.
#' @param expectation Nested expectation step including Gaussian moments.
#' @param object A fitted full-covariance model.
#' @return Mean and log-Cholesky score vectors.
#' @noRd
.multilpa_full_measurement_score <- function(parameters, expectation, object) {
  stopifnot(is.list(parameters), is.list(expectation), inherits(object, "multilpa"),
            identical(object$covariance_model, "full"))
  dimension <- length(object$indicators)
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  component_scores <- lapply(seq_len(object$n_profiles), function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance <- matrix(parameters$covariances[, , profile], dimension, dimension)
    precision <- chol2inv(chol(covariance))
    moments <- expectation$gaussian_moments[[profile]]
    weights <- expectation$subject_posteriors[, profile]
    residuals <- sweep(moments$expected, 2L, parameters$means[profile, ], "-")
    scatter <- crossprod(residuals, residuals * weights) +
      Reduce(`+`, lapply(moments$adjustments, function(adjustment) {
        stopifnot(is.list(adjustment))
        adjustment$covariance * sum(weights[adjustment$rows])
      }))
    list(means = as.vector(precision %*% colSums(residuals * weights)),
      covariance = 0.5 * (precision %*% scatter %*% precision - sum(weights) * precision))
  })
  profiles <- if (object$variance_model == "equal") 1L else seq_len(object$n_profiles)
  covariance_scores <- lapply(profiles, function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance_score <- if (object$variance_model == "equal") {
      Reduce(`+`, lapply(component_scores, `[[`, "covariance"))
    } else component_scores[[profile]]$covariance
    factor <- t(chol(matrix(parameters$covariances[, , profile], dimension, dimension)))
    factor_score <- 2 * covariance_score %*% factor
    diag(factor_score) <- diag(factor_score) * diag(factor)
    factor_score[lower]
  })
  list(means = unlist(lapply(component_scores, `[[`, "means"), use.names = FALSE),
       covariances = unlist(covariance_scores, use.names = FALSE))
}

#' Delta-method Jacobian for mixture parameters
#' @param object A fitted model.
#' @return Jacobian mapping unconstrained to natural coefficients.
#' @noRd
.multilpa_inference_jacobian <- function(object) {
  stopifnot(inherits(object, "multilpa"))
  natural <- .multilpa_coefficients(object, "natural")
  theta <- .multilpa_coefficients(object, "unconstrained")
  n_means <- length(object$means)
  full_covariance <- identical(object$covariance_model, "full")
  dimension <- length(object$indicators)
  n_covariance <- if (full_covariance) dimension * (dimension + 1L) / 2L else dimension
  n_variances <- n_covariance * if (object$variance_model == "equal") 1L else object$n_profiles
  n_measurement <- n_means + n_variances
  jacobian <- matrix(0, length(natural), length(theta), dimnames = list(names(natural), names(theta)))
  diag(jacobian)[seq_len(n_means)] <- 1
  if (full_covariance) {
    lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
    positions <- which(lower, arr.ind = TRUE)
    invisible(lapply(seq_len(n_variances / n_covariance), function(profile) {
      stopifnot(is.numeric(profile), length(profile) == 1L)
      factor <- t(chol(matrix(object$covariances[, , profile], dimension, dimension)))
      block <- vapply(seq_len(n_covariance), function(index) {
        stopifnot(is.numeric(index), length(index) == 1L)
        derivative <- matrix(0, dimension, dimension)
        row <- positions[index, 1L]
        column <- positions[index, 2L]
        derivative[row, column] <- if (row == column) factor[row, column] else 1
        covariance_derivative <- derivative %*% t(factor) + factor %*% t(derivative)
        covariance_derivative[lower]
      }, numeric(n_covariance))
      indices <- n_means + (profile - 1L) * n_covariance + seq_len(n_covariance)
      jacobian[indices, indices] <<- block
    }))
  } else {
    diag(jacobian)[n_means + seq_len(n_variances)] <- natural[n_means + seq_len(n_variances)]
  }
  invisible(lapply(seq_len(object$n_group_classes), function(group) {
    stopifnot(is.numeric(group), length(group) == 1L)
    probabilities <- object$profile_probabilities[group, ]
    block <- diag(probabilities, nrow = length(probabilities)) - tcrossprod(probabilities)
    natural_rows <- n_measurement + (group - 1L) * object$n_profiles + seq_len(object$n_profiles)
    free_columns <- n_measurement + (group - 1L) * (object$n_profiles - 1L) + seq_len(object$n_profiles - 1L)
    jacobian[natural_rows, free_columns] <<- block[, seq_len(object$n_profiles - 1L), drop = FALSE]
  }))
  probabilities <- object$group_probabilities
  block <- diag(probabilities, nrow = length(probabilities)) - tcrossprod(probabilities)
  natural_rows <- n_measurement + object$n_group_classes * object$n_profiles + seq_len(object$n_group_classes)
  free_columns <- n_measurement + object$n_group_classes * (object$n_profiles - 1L) + seq_len(object$n_group_classes - 1L)
  jacobian[natural_rows, free_columns] <- block[, seq_len(object$n_group_classes - 1L), drop = FALSE]
  jacobian
}

#' Observed-information inference for multilevel LPA
#'
#' Computes ordinary maximum-likelihood standard errors from the inverse
#' observed Hessian and delta-method Wald confidence intervals. These are local,
#' asymptotic estimates conditional on the selected numbers of classes. They are
#' not robust sandwich errors and do not account for model selection. Mixture
#' boundary solutions and unidentified Hessians do not admit this calculation.
#'
#' @param object A converged `multilpa` fit with inactive variance bounds.
#' @param data The original fitting data in the original row order. Stored
#'   indicator data and group identifiers are checked exactly. For older fits
#'   without stored indicators, the likelihood provides a weaker consistency check.
#' @param level Confidence level strictly between zero and one.
#' @param step Positive finite-difference step for the analytic score derivative.
#' @param vcov_type `"observed"` inverts the observed information.
#'   `"robust"` returns the Huber-White sandwich `A^-1 B A^-1`, where `B`
#'   accumulates the outer product of the per-group scores. Groups are the
#'   independent units, so robust errors relax the assumption that the Gaussian
#'   within-group model is correctly specified. They do not relax the assumption
#'   that groups are independent.
#' @return An `multilpa_inference` list containing natural `estimates`,
#'   `standard_errors`, `covariance`, `confidence_intervals`, unconstrained
#'   estimates and covariance, observed Hessian, score, and its eigenvalue
#'   condition ratio. All class probabilities are returned; their sum
#'   constraints make the natural covariance singular. Wald intervals on the
#'   natural scale are not clipped to parameter bounds. `vcov_type` records
#'   which covariance was returned. When it is `"robust"`, `group_scores` holds
#'   the groups-by-parameters score matrix and `scaling_correction` the MLR
#'   scaling correction factor `tr(A^-1 B) / q`, used for scaled likelihood-ratio
#'   difference tests. Use [as.data.frame()] for a tidy estimate table.
#' @examples
#' set.seed(42)
#' dat <- data.frame(group = rep(1:10, each = 10), y = rnorm(100))
#' fit <- multilpa(dat, "y", "group", 1, 1, n_starts = 1)
#' parameter_inference(fit, dat)$standard_errors
#' @export
parameter_inference <- function(object, data, level = 0.95, step = 1e-4,
                             vcov_type = c("observed", "robust")) {
  stopifnot(inherits(object, "multilpa"), is.data.frame(data),
            is.numeric(level), length(level) == 1L, is.finite(level), level > 0, level < 1,
            is.numeric(step), length(step) == 1L, is.finite(step), step > 0)
  vcov_type <- match.arg(vcov_type)
  .multilpa_check_regularity(object, vcov_type)
  theta <- .multilpa_coefficients(object, "unconstrained")
  prepared <- .multilpa_inference_matrix(object, data)
  x <- prepared$x
  centered_object <- object
  centered_object$means <- sweep(object$means, 2L, prepared$centers, "-")
  centered_theta <- .multilpa_coefficients(centered_object, "unconstrained")
  objective <- function(parameters) {
    stopifnot(is.numeric(parameters))
    -.multilpa_expectation(x, object$group_index,
      .multilpa_decode(parameters, object))$log_likelihood
  }
  score <- function(parameters) {
    stopifnot(is.numeric(parameters))
    .multilpa_score(parameters, x, object)
  }
  fitted_likelihood <- -objective(centered_theta)
  if (abs(fitted_likelihood - object$log_likelihood) > 1e-8 * (1 + abs(object$log_likelihood))) {
    stop("data do not reproduce the fitted log likelihood; supply the original fitting data.")
  }
  parameter_scale <- c(as.vector(t(sqrt(object$variances))),
                       rep(1, length(theta) - length(object$means)))
  scaled_objective <- function(displacement) {
    stopifnot(is.numeric(displacement))
    objective(centered_theta + displacement * parameter_scale)
  }
  scaled_gradient <- function(displacement) {
    stopifnot(is.numeric(displacement))
    score(centered_theta + displacement * parameter_scale) * parameter_scale
  }
  information <- .multilpa_observed_hessian(scaled_objective, scaled_gradient,
                                          parameter_scale, step)
  scaled_hessian <- information$scaled
  hessian <- information$natural
  condition_ratio <- information$condition_ratio
  scaled_inverse <- information$inverse
  group_scores <- NULL
  scaling_correction <- NA_real_
  if (identical(vcov_type, "robust")) {
    group_scores <- .multilpa_group_scores(centered_theta, x, object)
    scaled_cross <- .multilpa_cross_product(sweep(group_scores, 2L, parameter_scale, "*"))
    scaling_correction <- sum(diag(scaled_inverse %*% scaled_cross)) / length(theta)
    scaled_inverse <- scaled_inverse %*% scaled_cross %*% scaled_inverse
  }
  covariance_unconstrained <- scaled_inverse * tcrossprod(parameter_scale)
  dimnames(hessian) <- dimnames(covariance_unconstrained) <- list(names(theta), names(theta))
  jacobian <- .multilpa_inference_jacobian(object)
  covariance <- jacobian %*% covariance_unconstrained %*% t(jacobian)
  standard_errors <- sqrt(pmax(diag(covariance), 0))
  estimates <- .multilpa_coefficients(object, "natural")
  critical <- stats::qnorm((1 + level) / 2)
  intervals <- cbind(estimates - critical * standard_errors, estimates + critical * standard_errors)
  colnames(intervals) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%")
  gradient <- stats::setNames(score(centered_theta), names(theta))
  scaled_score <- max(abs(gradient * parameter_scale))
  if (scaled_score > 0.01) warning("The fitted likelihood has a non-negligible score; refit with a tighter tolerance before using Wald inference.", call. = FALSE)
  result <- list(estimates = estimates, standard_errors = standard_errors, covariance = covariance,
    confidence_intervals = intervals, coefficients_unconstrained = theta,
    covariance_unconstrained = covariance_unconstrained, hessian = hessian,
    gradient = gradient, scaled_score = scaled_score,
    condition_ratio = condition_ratio, level = level, step = step,
    vcov_type = vcov_type, group_scores = group_scores,
    scaling_correction = scaling_correction)
  class(result) <- "multilpa_inference"
  result
}

#' Extract multilevel LPA coefficients
#' @param object A fitted `multilpa` model.
#' @param scale Natural coefficients or unconstrained log variances (diagonal),
#'   log-Cholesky coordinates (full covariance), and baseline-category logits.
#' @param ... Reserved for generic compatibility.
#' @return A named vector. Natural coefficients include all mixing probabilities;
#'   unconstrained coefficients exclude their reference categories.
#' @examples
#' # After fitting: coef(fit)
#' @export
#' @importFrom stats coef
coef.multilpa <- function(object, scale = c("natural", "unconstrained"), ...) {
  stopifnot(inherits(object, "multilpa"))
  .multilpa_coefficients(object, match.arg(scale))
}

#' Extract multilevel LPA covariance estimates
#' @param object A fitted `multilpa` model.
#' @param data Original fitting data, required unless `object$inference` is stored.
#' @param scale Natural or unconstrained parameter scale.
#' @param ... Additional arguments passed to [parameter_inference()].
#' @return The observed-information covariance matrix. On the natural scale,
#'   probability sum constraints make this matrix singular by construction.
#' @examples
#' # After fitting: vcov(fit, data = original_data)
#' @export
#' @importFrom stats vcov
vcov.multilpa <- function(object, data = NULL, scale = c("natural", "unconstrained"), ...) {
  stopifnot(inherits(object, "multilpa"))
  scale <- match.arg(scale)
  information <- if (is.null(data)) object$inference else parameter_inference(object, data, ...)
  if (is.null(information)) stop("Supply original data or store parameter_inference() in object$inference.")
  if (scale == "natural") information$covariance else information$covariance_unconstrained
}

#' Wald confidence intervals for multilevel LPA coefficients
#' @param object A fitted `multilpa` model.
#' @param parm Optional coefficient names or indices; defaults to all coefficients.
#' @param level Confidence level strictly between zero and one.
#' @param data Original fitting data, required unless `object$inference` is stored.
#' @param ... Additional arguments passed to [parameter_inference()].
#' @return A two-column matrix of natural-scale Wald intervals. Bounds are not
#'   clipped to the probability or variance parameter space.
#' @examples
#' # After fitting: confint(fit, data = original_data)
#' @export
#' @importFrom stats confint
confint.multilpa <- function(object, parm, level = 0.95, data = NULL, ...) {
  stopifnot(inherits(object, "multilpa"), is.numeric(level), length(level) == 1L,
            is.finite(level), level > 0, level < 1)
  estimates <- coef.multilpa(object)
  covariance <- vcov.multilpa(object, data = data, ...)
  if (missing(parm)) parm <- names(estimates)
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm)) || any(parm != floor(parm)) ||
        any(parm < 1) || any(parm > length(estimates))) stop("Invalid coefficient indices in parm.")
    parm <- names(estimates)[parm]
  }
  if (!is.character(parm) || anyNA(parm) || !all(parm %in% names(estimates))) {
    stop("parm must identify existing coefficient names or indices.")
  }
  standard_errors <- sqrt(pmax(diag(covariance)[parm], 0))
  critical <- stats::qnorm((1 + level) / 2)
  intervals <- cbind(estimates[parm] - critical * standard_errors,
                     estimates[parm] + critical * standard_errors)
  colnames(intervals) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%")
  intervals
}

#' Check that a fit admits Wald inference at all
#'
#' Boundary solutions, unconverged fits and zero mixing probabilities each put
#' the estimate where the asymptotic normal approximation does not hold, so
#' they are refused rather than reported with a number.
#'
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_check_regularity <- function(object, vcov_type) {
  if (!is.null(object$response_probabilities)) {
    stop(errorCondition("Standard errors are not yet available for categorical indicators; the score functions cover Gaussian measurement only.",
                        class = "multilpa_unsupported_inference", call = NULL))
  }
  if (identical(vcov_type, "robust") && object$n_groups <= object$n_parameters) {
    stop(errorCondition(sprintf(
      "Robust inference needs more groups than parameters; this fit has %d groups and %d parameters.",
      object$n_groups, object$n_parameters),
      class = "multilpa_too_few_groups", call = NULL))
  }
  if (!isTRUE(object$converged)) stop("Inference requires a converged fit.")
  if (isTRUE(object$boundary)) {
    stop("Wald inference is unavailable for a bound-active fit.")
  }
  if (any(object$profile_probabilities <= 0) || any(object$group_probabilities <= 0)) {
    stop("Wald inference requires strictly positive mixing probabilities.")
  }
  invisible(NULL)
}

#' Recover and centre the indicator matrix the fit was built on
#'
#' The Hessian is evaluated at the fitted estimates, so it is only meaningful
#' if `data` really is the data that produced them. Stored indicators are
#' compared exactly; older fits without them fall back to the likelihood check
#' the caller performs next.
#'
#' @return A list with the centred matrix `x` and the `centers` removed from it.
#' @noRd
.multilpa_inference_matrix <- function(object, data) {
  if (nrow(data) != object$n_observations ||
      !all(c(object$indicators, object$cluster) %in% names(data)) ||
      anyDuplicated(names(data)) ||
      !all(vapply(data[, object$indicators, drop = FALSE], is.numeric, logical(1)))) {
    stop("data must contain the original numeric indicators and cluster column.")
  }
  group_index <- match(data[[object$cluster]], object$group_values)
  if (!identical(group_index, object$group_index)) {
    stop("data must retain the original group identifiers and row order.")
  }
  x <- as.matrix(data[, object$indicators, drop = FALSE])
  if (any(is.infinite(x)) || any(is.nan(x)) ||
      (anyNA(x) && !identical(object$missing, "fiml"))) {
    stop("data contain unsupported missing or non-finite indicators.")
  }
  if (!is.null(object$indicator_data) && !identical(x, object$indicator_data)) {
    stop("data must reproduce the original indicator data, including row order and names.")
  }
  centers <- colMeans(x, na.rm = TRUE)
  list(x = sweep(x, 2L, centers, "-"), centers = centers)
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
      condition_ratio <= 1e-10) {
    stop("Observed information is not positive definite or is numerically singular; Wald inference is unavailable.")
  }
  inverse <- tryCatch(solve(scaled), error = function(error) {
    stop(sprintf("Observed information could not be inverted: %s",
                 conditionMessage(error)))
  })
  list(scaled = scaled, natural = natural, condition_ratio = condition_ratio,
       inverse = inverse)
}
