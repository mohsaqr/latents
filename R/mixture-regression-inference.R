# Standard errors for mixture_regression() fits.
#
# The free parameters are packed into one unconstrained vector: the class
# regression coefficients, the shared coefficients, the log residual standard
# deviations (Gaussian), and the free columns of every multinomial-logit block
# (the first class and the first group class are the references).
#
# Scores are analytic. By Fisher's identity the score of a unit's observed-data
# log likelihood is the posterior expectation of its complete-data score, so
# each block's score is its complete-data score weighted by the posteriors the
# E-step already holds. The observed information is the numerical Jacobian of
# the total analytic score, which is accurate to the step's truncation error
# rather than the square of it, as a second difference of the likelihood
# would be.

#' Pack the free parameters into a named vector
#' @param spec The model specification.
#' @param params The parameter list.
#' @return A named numeric vector.
#' @noRd
.mixture_pack <- function(spec, params) {
  class_names <- paste0("class_", seq_len(spec$n_classes))
  regression <- stats::setNames(
    as.vector(params$beta),
    paste0("coefficient.", rep(class_names, each = ncol(spec$x)), ".",
           rep(colnames(spec$x), spec$n_classes)))
  common <- if (ncol(spec$z) == 0L) numeric() else stats::setNames(
    params$common, paste0("coefficient.common.", colnames(spec$z)))
  dispersion <- switch(spec$family,
    gaussian = if (identical(spec$variance, "equal")) {
      c(log_sigma.all = 0.5 * log(params$sigma2[1L]))
    } else stats::setNames(0.5 * log(params$sigma2), paste0("log_sigma.", class_names)),
    negative_binomial = if (identical(spec$variance, "equal")) {
      c(log_dispersion.all = log(params$sigma2[1L]))
    } else stats::setNames(log(params$sigma2), paste0("log_dispersion.", class_names)),
    ordinal = .ordinal_pack(spec, params),
    numeric())
  c(regression, common, dispersion, .mixture_structure_pack(spec, params))
}

#' Unpack a free-parameter vector into a parameter list
#' @param spec The model specification.
#' @param theta Packed vector.
#' @param template A parameter list supplying the shapes.
#' @return A parameter list.
#' @noRd
.mixture_unpack <- function(spec, theta, template) {
  k <- spec$n_classes
  p <- ncol(spec$x)
  q <- ncol(spec$z)
  position <- 0L
  take <- function(count) {
    values <- theta[position + seq_len(count)]
    position <<- position + count
    unname(values)
  }
  params <- template
  params$beta[] <- take(p * k)
  params$common[] <- take(q)
  if (identical(spec$family, "negative_binomial")) {
    params$sigma2 <- if (identical(spec$variance, "equal")) rep(exp(take(1L)), k) else
      exp(take(k))
  }
  if (identical(spec$family, "gaussian")) {
    params$sigma2 <- if (identical(spec$variance, "equal")) {
      rep(exp(2 * take(1L)), k)
    } else exp(2 * take(k))
  }
  if (identical(spec$family, "ordinal")) {
    m <- spec$n_categories - 1L
    params$thresholds[] <- .ordinal_natural(take(m * k), m)
  }
  .mixture_structure_unpack(spec, theta[position + seq_len(length(theta) - position)],
                            params)
}

#' Per-unit analytic scores
#'
#' @param spec The model specification.
#' @param params The parameter list.
#' @param expectation The E-step at `params`.
#' @return A matrix with one row per top-level unit (row, or group) and one
#'   column per packed parameter.
#' @noRd
.mixture_unit_scores <- function(spec, params, expectation) {
  tau <- expectation$tau
  eta <- .mixture_linear_predictors(spec, params)
  ordinal <- if (identical(spec$family, "ordinal")) {
    .ordinal_derivatives(spec, params, eta)
  }
  # The ordinal derivative in eta is -(d/da + d/db) over the category bounds.
  derivative <- if (is.null(ordinal)) {
    .mixture_family_derivatives(spec$family, spec$y, spec$trials, eta,
                                params$sigma2)$gradient
  } else -(ordinal$upper + ordinal$lower)
  weighted <- tau * derivative
  regression <- do.call(cbind, lapply(seq_len(spec$n_classes), function(k) {
    spec$x * weighted[, k]
  }))
  common <- spec$z * rowSums(weighted)
  dispersion <- if (identical(spec$family, "negative_binomial")) {
    standardized <- tau * .latents_negative_binomial_dispersion_derivatives(
      spec$y, exp(eta), params$sigma2)$score
    if (identical(spec$variance, "equal")) {
      matrix(rowSums(standardized), ncol = 1L)
    } else standardized
  } else if (!is.null(ordinal)) {
    .ordinal_threshold_scores(spec, params, tau, ordinal)
  } else if (!identical(spec$family, "gaussian")) {
    matrix(0, spec$n, 0L)
  } else {
    residual2 <- (spec$y - eta)^2
    standardized <- tau * (sweep(residual2, 2L, params$sigma2, "/") - 1)
    if (identical(spec$variance, "equal")) {
      matrix(rowSums(standardized), ncol = 1L)
    } else standardized
  }
  row_scores <- cbind(regression, common, dispersion)
  to_units <- function(scores) {
    if (identical(spec$nesting, "observation")) scores else
      rowsum(scores, spec$group_index, reorder = TRUE)
  }
  mixing <- .mixture_structure_scores(spec, expectation)
  scores <- cbind(to_units(row_scores), mixing)
  colnames(scores) <- names(.mixture_pack(spec, params))
  scores
}

#' Total score at a packed parameter vector
#' @noRd
.mixture_total_score <- function(spec, theta, template) {
  params <- .mixture_unpack(spec, theta, template)
  colSums(.mixture_unit_scores(spec, params,
                              .mixture_expectation(spec, params)))
}

#' Observed information by differentiating the analytic score
#' @param spec The model specification.
#' @param params The parameter list at the estimate.
#' @param step Relative step for central differences.
#' @return A symmetric matrix.
#' @noRd
.mixture_observed_information <- function(spec, params, step = 1e-5) {
  # Central differences with h = step * (1 + |theta|), the shared service's
  # "central"/"plus1" variant.
  .inference_score_information(function(point) .mixture_total_score(spec, point, params),
                               .mixture_pack(spec, params), step,
                               scheme = "central", h_rule = "plus1")
}

#' Covariance of the packed estimates
#'
#' @param fit A `latents_mixture_regression` fit.
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`.
#' @return A list with `theta`, `vcov`, `vcov_type`, `clusters` and the
#'   `information` matrix.
#' @noRd
.mixture_inference <- function(fit, vcov_type = c("observed", "robust", "opg")) {
  vcov_type <- match.arg(vcov_type)
  spec <- fit$spec
  params <- fit$params
  theta <- .mixture_pack(spec, params)
  scores <- .mixture_unit_scores(spec, params, .mixture_expectation(spec, params))
  clustered_on <- if (identical(spec$nesting, "observation") &&
                      !is.null(spec$id)) spec$id else if (
    !identical(spec$nesting, "observation")) spec$id else "row"
  if (identical(spec$nesting, "observation") && !is.null(spec$id)) {
    scores_for_meat <- rowsum(scores, spec$group_index, reorder = TRUE)
  } else scores_for_meat <- scores
  meat <- crossprod(scores_for_meat)
  information <- if (identical(vcov_type, "opg")) crossprod(scores) else
    .mixture_observed_information(spec, params)
  # A negative-binomial dispersion at the Poisson limit is a boundary
  # estimate: the likelihood is flat there, so it is held, its errors are
  # NA, and the others are conditional on it.
  free <- !.mixture_boundary_coordinates(spec, theta)
  if (identical(vcov_type, "robust")) {
    .multilpa_require_clusters(nrow(scores_for_meat), sum(free),
                               "A robust regression covariance",
                               unit = "independent units")
    .inference_require_rank(scores_for_meat[, free, drop = FALSE],
                             "The unit scores are rank deficient; robust inference is unavailable.")
  }
  inverse <- .mixture_invert(information[free, free, drop = FALSE])
  free_vcov <- if (identical(vcov_type, "robust")) {
    .inference_sandwich(inverse, meat[free, free, drop = FALSE])
  } else inverse
  vcov <- matrix(NA_real_, length(theta), length(theta))
  vcov[free, free] <- free_vcov
  dimnames(vcov) <- list(names(theta), names(theta))
  list(theta = theta, vcov = vcov, vcov_type = vcov_type,
       clustered_on = if (identical(vcov_type, "robust")) clustered_on else NA,
       information = information,
       score_norm = max(abs(colSums(scores))))
}

#' Coordinates of a regression mixture on a boundary
#'
#' A negative-binomial log dispersion at its floor (the Poisson limit).
#' @param spec The model specification.
#' @param theta The packed estimates.
#' @return A logical vector over `theta`.
#' @noRd
.mixture_boundary_coordinates <- function(spec, theta) {
  startsWith(names(theta), "log_dispersion.") &
    theta <= log(.latents_min_dispersion) + 1e-6
}

#' Invert an information matrix, refusing a singular one
#' @noRd
.mixture_invert <- function(information) {
  # The "min <= 1e-10 max |ev|" rule flags rather than refuses: a mixture fit
  # reports NA standard errors for an unidentified estimate.
  if (!.inference_positive_definite(information, "min_le_eps_maxabs")) {
    warning(warningCondition(paste(
      "The information matrix is singular or not positive definite: some",
      "parameters are not identified at this estimate (an empty class, a",
      "boundary solution or separation). Their standard errors are NA."),
      class = "latents_singular_information", call = NULL))
    out <- matrix(NA_real_, nrow(information), ncol(information))
    return(out)
  }
  chol2inv(chol(information))
}

#' Delta-method standard errors for a transformation of the estimates
#'
#' @param theta Packed estimates.
#' @param vcov Their covariance.
#' @param transform Function from a packed vector to the derived vector.
#' @return A list of `estimate` and `std_error`.
#' @noRd
.mixture_delta <- function(theta, vcov, transform, step = 1e-6) {
  # The mixture engines have always used the rowSums((J V) * J) form.
  jacobian <- .inference_numeric_jacobian(transform, theta, step)
  if (anyNA(vcov)) {
    errors <- vapply(seq_len(nrow(jacobian)), function(i) {
      used <- which(jacobian[i, ] != 0)
      if (!length(used)) return(0)
      .inference_delta_se(jacobian[i, used, drop = FALSE],
                          vcov[used, used, drop = FALSE], "rowsums")
    }, numeric(1))
    return(list(estimate = transform(theta), std_error = errors))
  }
  list(estimate = transform(theta),
       std_error = .inference_delta_se(jacobian, vcov, "rowsums"))
}
