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
      c(log_sigma.all = 0.5 * log(params$sigma2[[1L]]))
    } else stats::setNames(0.5 * log(params$sigma2), paste0("log_sigma.", class_names)),
    negative_binomial = if (identical(spec$variance, "equal")) {
      c(log_dispersion.all = log(params$sigma2[[1L]]))
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
#'
#' Central differences with each coordinate's step a thousandth of its own
#' standard-error scale, `1 / sqrt(sum of squared unit scores)`. That scale
#' moves with the parameter's units and origin. A step of `1e-5 (1 + |theta|)`
#' did not: with a predictor shifted by 2000 the intercept (about -9000) was
#' differenced with a step of 0.09, and its standard error moved by 0.2%
#' (review 2026-10-26, MR-F6).
#'
#' @param spec The model specification.
#' @param params The parameter list at the estimate.
#' @param scores The unit scores at the estimate.
#' @param step Step as a fraction of each coordinate's standard-error scale.
#' @return A symmetric matrix.
#' @noRd
.mixture_observed_information <- function(spec, params, scores, step = 1e-3) {
  .inference_score_information(function(point) .mixture_total_score(spec, point, params),
                               .mixture_pack(spec, params), step,
                               scheme = "central", h_rule = "scaled",
                               scale = 1 / sqrt(colSums(scores^2)))
}

#' The centred and scaled chart the information is differenced in
#'
#' Every non-constant design column of the regression (`x`, `z`) and of the
#' membership designs (`w`, `v`) is centred on its mean and divided by its
#' standard deviation where an intercept absorbs the shift (a class intercept,
#' an ordinal class's thresholds, a membership intercept, or the group-class
#' intercepts of the two-level class logits), and only scaled otherwise. The
#' likelihood is the same function in either chart, and the packed vectors
#' are linked exactly: phi = A theta and theta = J phi, J = A^-1.
#'
#' In the supplied chart an intercept and the slope of a predictor whose mean
#' lies far from zero relative to its spread are nearly collinear (a
#' correlation of 1 - (sd / mean)^2), so differencing errors in the information
#' are amplified by (mean / sd)^2: standard errors lost accuracy from a shift of
#' 1e4 standard deviations and were all NA from 5e4 (review 2026-10-26, round
#' 2, N6). The pattern is that of `.lta_chart_jacobian()`.
#'
#' @param spec The model specification.
#' @return A list with the chart's `spec`, `to_chart` (A) and `to_raw` (J),
#'   both square over the packed vector.
#' @noRd
.mixture_chart <- function(spec) {
  part <- function(design, absorbed) {
    if (ncol(design) == 0L) return(list(centre = numeric(), scale = numeric()))
    constant <- apply(design, 2L, function(column) max(column) == min(column))
    centre <- if (absorbed) colMeans(design) else numeric(ncol(design))
    centre[constant] <- 0
    shifted <- sweep(design, 2L, centre, "-")
    scale <- sqrt(colMeans(shifted^2))
    scale[constant | !(scale > 0)] <- 1
    list(centre = centre, scale = scale)
  }
  apply_part <- function(design, chart) {
    if (ncol(design) == 0L) return(design)
    sweep(sweep(design, 2L, chart$centre, "-"), 2L, chart$scale, "/")
  }
  ordinal <- identical(spec$family, "ordinal")
  two_level <- identical(spec$nesting, "two-level")
  k <- spec$n_classes
  p <- ncol(spec$x)
  q <- ncol(spec$z)
  intercept <- match("(Intercept)", colnames(spec$x), nomatch = 0L)
  absorbed <- ordinal || intercept > 0L
  x_part <- part(spec$x, absorbed)
  z_part <- part(spec$z, absorbed)
  w_part <- part(spec$w, two_level || "(Intercept)" %in% colnames(spec$w))
  v_part <- if (two_level) part(spec$v, "(Intercept)" %in% colnames(spec$v)) else NULL
  chart_spec <- spec
  chart_spec$x <- apply_part(spec$x, x_part)
  chart_spec$z <- apply_part(spec$z, z_part)
  chart_spec$w <- apply_part(spec$w, w_part)
  if (two_level) chart_spec$v <- apply_part(spec$v, v_part)

  n_theta <- length(.mixture_pack(spec, .mixture_null_params(spec)))
  # phi = (D + C) theta: D scales the covariate coordinates, C moves their
  # centres into the absorbing coordinate (an intercept or first threshold).
  # C's columns are covariate coordinates, whose own rows of C are empty, so
  # (D + C)^-1 = D^-1 - C D^-1 exactly.
  scale <- rep(1, n_theta)
  centre <- matrix(0, n_theta, n_theta)
  regression <- function(class) (class - 1L) * p + seq_len(p)
  common <- p * k + seq_len(q)
  if (q > 0L) scale[common] <- z_part$scale
  dispersion_start <- p * k + q
  invisible(lapply(seq_len(k), function(class) {
    block <- regression(class)
    scale[block] <<- x_part$scale
    target <- if (ordinal) {
      dispersion_start + (class - 1L) * (spec$n_categories - 1L) + 1L
    } else if (intercept > 0L) block[intercept] else NA_integer_
    if (!is.na(target)) {
      sign <- if (ordinal) -1 else 1
      centre[target, block] <<- sign * x_part$centre
      if (q > 0L) centre[target, common] <<- sign * z_part$centre
    }
  }))
  structure_start <- n_theta - length(.mixture_structure_pack(spec, .mixture_null_params(spec)))
  # Two-level class logits: per free class, one intercept per group class
  # (which absorb the centres) and then the shared covariate slopes.
  place <- function(start, n_blocks, n_intercepts, part) {
    width <- n_intercepts + length(part$centre)
    invisible(lapply(seq_len(n_blocks), function(index) {
      first <- start + (index - 1L) * width
      columns <- first + n_intercepts + seq_along(part$centre)
      scale[columns] <<- part$scale
      if (n_intercepts > 0L) {
        centre[first + seq_len(n_intercepts), columns] <<-
          matrix(part$centre, n_intercepts, length(part$centre), byrow = TRUE)
      }
    }))
    start + n_blocks * width
  }
  # A membership design with its own intercept column, one block per free
  # class (or group class).
  fill_design_blocks <- function(start, n_blocks, design, part) {
    width <- ncol(design)
    own <- match("(Intercept)", colnames(design), nomatch = 0L)
    invisible(lapply(seq_len(n_blocks), function(index) {
      columns <- start + (index - 1L) * width + seq_len(width)
      scale[columns] <<- part$scale
      if (own > 0L) centre[columns[own], columns] <<- part$centre
    }))
    start + n_blocks * width
  }
  if (two_level) {
    next_start <- fill_design_blocks(structure_start, spec$n_group_classes - 1L,
                                     spec$v, v_part)
    place(next_start, k - 1L, spec$n_group_classes, w_part)
  } else {
    fill_design_blocks(structure_start, k - 1L, spec$w, w_part)
  }
  to_chart <- diag(scale, n_theta) + centre
  to_raw <- diag(1 / scale, n_theta) - sweep(centre, 2L, scale, "/")
  list(spec = chart_spec, to_chart = to_chart, to_raw = to_raw)
}

#' Covariance of the packed estimates
#'
#' The scores and the information are computed in the centred and scaled
#' chart of `.mixture_chart()` and mapped back exactly: with theta = J phi the
#' covariance is J V J', the information A' I A and the scores S A.
#'
#' @param fit A `latents_mixture_regression` fit.
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`.
#' @return A list with `theta`, `vcov`, `vcov_type`, `clusters` and the
#'   `information` matrix.
#' @noRd
.mixture_inference <- function(fit, vcov_type = c("observed", "robust", "opg")) {
  vcov_type <- match.arg(vcov_type)
  # Excluded units cannot set the chart's origin or spread, contribute score
  # rows, or affect rank and cluster counts. Keep their data/posteriors on the
  # fitted object, but differentiate the same support used by its fit.
  spec <- .mixture_carrying_spec(fit$spec)
  params <- fit$params
  theta <- .mixture_pack(spec, params)
  chart <- .mixture_chart(spec)
  chart_spec <- chart$spec
  phi <- stats::setNames(drop(chart$to_chart %*% theta), names(theta))
  chart_params <- .mixture_unpack(chart_spec, phi, params)
  scores <- .mixture_unit_scores(chart_spec, chart_params,
                                 .mixture_expectation(chart_spec, chart_params))
  clustered_on <- if (identical(spec$nesting, "observation") &&
                      !is.null(spec$id)) spec$id else if (
    !identical(spec$nesting, "observation")) spec$id else "row"
  if (identical(spec$nesting, "observation") && !is.null(spec$id)) {
    scores_for_meat <- rowsum(scores, spec$group_index, reorder = TRUE)
  } else scores_for_meat <- scores
  n_units <- sum((spec$sampling_weights %||% rep(1, nrow(scores_for_meat))) > 0)
  meat <- crossprod(scores_for_meat)
  information <- if (identical(vcov_type, "opg")) crossprod(scores) else
    .mixture_observed_information(chart_spec, chart_params, scores)
  # A negative-binomial dispersion at the Poisson limit, or an ordinal
  # threshold gap at zero, is a boundary estimate: the likelihood is flat
  # there, so it is held, its errors are NA, and the others are conditional on
  # it. These coordinates are the same in both charts.
  free <- !.mixture_boundary_coordinates(spec, theta)
  if (identical(vcov_type, "robust")) {
    .multilpa_require_clusters(n_units, sum(free),
                               "A robust regression covariance",
                               unit = "independent units")
    .inference_require_rank(scores_for_meat[, free, drop = FALSE],
                             "The unit scores are rank deficient; robust inference is unavailable.")
  }
  inverse <- .mixture_invert(information[free, free, drop = FALSE])
  free_vcov <- if (identical(vcov_type, "robust")) {
    .inference_sandwich(inverse, meat[free, free, drop = FALSE], n_units)
  } else inverse
  to_raw <- chart$to_raw[free, free, drop = FALSE]
  vcov <- matrix(NA_real_, length(theta), length(theta))
  vcov[free, free] <- to_raw %*% free_vcov %*% t(to_raw)
  vcov[free, free] <- (vcov[free, free] + t(vcov[free, free])) / 2
  dimnames(vcov) <- list(names(theta), names(theta))
  raw_information <- crossprod(chart$to_chart, information %*% chart$to_chart)
  dimnames(raw_information) <- list(names(theta), names(theta))
  list(theta = theta, vcov = vcov, vcov_type = vcov_type,
       clustered_on = if (identical(vcov_type, "robust")) clustered_on else NA,
       information = raw_information,
       score_norm = max(abs(colSums(scores %*% chart$to_chart))))
}

#' Coordinates of a regression mixture on a boundary
#'
#' A negative-binomial log dispersion within 1e-6 (relative) of its floor
#' (the Poisson limit), or an ordinal log threshold gap below
#' `.ordinal_min_log_gap` (a category the class never uses: two thresholds
#' coincide, review 2026-10-26, round 2, N5).
#' @param spec The model specification.
#' @param theta The packed estimates.
#' @return A logical vector over `theta`.
#' @noRd
.mixture_boundary_coordinates <- function(spec, theta) {
  (startsWith(names(theta), "log_dispersion.") &
     theta <= log(.latents_min_dispersion) + 1e-6) |
    (startsWith(names(theta), "log_threshold_gap.") &
       theta < .ordinal_min_log_gap)
}

#' Invert an information matrix, refusing a singular one
#'
#' The decision and the inverse are taken on the correlation form `D H D`,
#' with `D = diag(H)^(-1/2)`, and mapped back as `D (D H D)^-1 D`. Positive
#' definiteness is unchanged by that scaling, so the decision no longer depends
#' on the units or the origin of any parameter. On the raw matrix, eigenvalues
#' of parameters in different units were compared with one another: with a
#' predictor shifted by 2000 (a calendar year) the intercept's row is about
#' 2000^2 times the slope's, and every standard error was reported NA even
#' with one class (review 2026-10-26, MR-F5). The origin is removed before
#' this, by differencing in the centred chart of `.mixture_chart()` (N6).
#'
#' @param information A symmetric information matrix.
#' @return Its inverse, or a matrix of `NA` with a
#'   `latents_singular_information` warning when it is singular.
#' @noRd
.mixture_invert <- function(information) {
  diagonal <- diag(information)
  usable <- length(diagonal) > 0L && all(is.finite(diagonal)) && all(diagonal > 0)
  scale <- if (usable) 1 / sqrt(diagonal) else rep(1, length(diagonal))
  conditioned <- information * tcrossprod(scale)
  # The "min <= 1e-10 max |ev|" rule flags rather than refuses: a mixture fit
  # reports NA standard errors for an unidentified estimate. A parameter with
  # no information at all (a zero diagonal) is unidentified outright.
  if (!usable || !.inference_positive_definite(conditioned, "min_le_eps_maxabs")) {
    warning(warningCondition(paste(
      "The information matrix is singular or not positive definite: some",
      "parameters are not identified at this estimate (an empty class, a",
      "boundary solution or separation). Their standard errors are NA."),
      class = "latents_singular_information", call = NULL))
    out <- matrix(NA_real_, nrow(information), ncol(information))
    return(out)
  }
  chol2inv(chol(conditioned)) * tcrossprod(scale)
}

#' Delta-method standard errors for a transformation of the estimates
#'
#' @param theta Packed estimates.
#' @param vcov Their covariance.
#' @param transform Function from a packed vector to the derived vector.
#' @return A list of `estimate` and `std_error`.
#' @noRd
.mixture_delta <- function(theta, vcov, transform, step = 1e-6) {
  # The mixture engines have always used the rowSums((J V) * J) form. Each
  # coordinate is stepped by a hundred-thousandth of its own standard error: a fixed
  # 1e-6 (1 + |theta|) was half of a Cholesky element of 2e-6, which an
  # outcome in units of 1e-5 produces, and a random-effect SD's error came out
  # 0.2% off (review 2026-10-26, G2).
  jacobian <- .inference_numeric_jacobian(transform, theta, step,
                                          scale = sqrt(pmax(diag(vcov), 0)))
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
