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
  # paste0() recycles a zero-length argument to one string, so an empty block
  # (one class, no group-class covariates) must be named explicitly as empty.
  labelled <- function(values, labels) {
    if (length(values) == 0L) numeric() else stats::setNames(values, labels)
  }
  class_names <- paste0("class_", seq_len(spec$n_classes))
  regression <- stats::setNames(
    as.vector(params$beta),
    paste0("coefficient.", rep(class_names, each = ncol(spec$x)), ".",
           rep(colnames(spec$x), spec$n_classes)))
  common <- if (ncol(spec$z) == 0L) numeric() else stats::setNames(
    params$common, paste0("coefficient.common.", colnames(spec$z)))
  dispersion <- if (!identical(spec$family, "gaussian")) numeric() else if (
    identical(spec$variance, "equal")) c(log_sigma.all = 0.5 * log(params$sigma2[1L])) else
      stats::setNames(0.5 * log(params$sigma2), paste0("log_sigma.", class_names))
  free_classes <- class_names[-1L]
  mixing <- if (identical(spec$nesting, "two-level")) {
    group_names <- paste0("group_class_", seq_len(spec$n_group_classes))
    delta <- params$delta[, -1L, drop = FALSE]
    logits <- params$class_logits[, -1L, drop = FALSE]
    logit_rows <- c(paste0("(Intercept):", group_names), colnames(spec$w))
    c(labelled(as.vector(delta), paste0(
      "group_membership.", rep(group_names[-1L], each = ncol(spec$v)), ".",
      rep(colnames(spec$v), spec$n_group_classes - 1L))),
      labelled(as.vector(logits), paste0(
        "membership.", rep(free_classes, each = length(logit_rows)), ".",
        rep(logit_rows, length(free_classes)))))
  } else {
    gamma <- params$gamma[, -1L, drop = FALSE]
    labelled(as.vector(gamma), paste0(
      "membership.", rep(free_classes, each = ncol(spec$w)), ".",
      rep(colnames(spec$w), length(free_classes))))
  }
  c(regression, common, dispersion, mixing)
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
  if (identical(spec$family, "gaussian")) {
    params$sigma2 <- if (identical(spec$variance, "equal")) {
      rep(exp(2 * take(1L)), k)
    } else exp(2 * take(k))
  }
  if (identical(spec$nesting, "two-level")) {
    delta_free <- take(ncol(spec$v) * (spec$n_group_classes - 1L))
    params$delta <- cbind(0, matrix(delta_free, ncol(spec$v)))
    logits_free <- take(nrow(template$class_logits) * (k - 1L))
    params$class_logits <- cbind(0, matrix(logits_free,
                                           nrow(template$class_logits)))
  } else {
    gamma_free <- take(ncol(spec$w) * (k - 1L))
    params$gamma <- cbind(0, matrix(gamma_free, ncol(spec$w)))
  }
  params
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
  derivative <- .mixture_family_derivatives(spec$family, spec$y, spec$trials,
                                           eta, params$sigma2)$gradient
  weighted <- tau * derivative
  regression <- do.call(cbind, lapply(seq_len(spec$n_classes), function(k) {
    spec$x * weighted[, k]
  }))
  common <- spec$z * rowSums(weighted)
  dispersion <- if (!identical(spec$family, "gaussian")) {
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
  mixing <- switch(spec$nesting,
    observation = {
      residual <- tau - rowSums(tau) * exp(expectation$log_prior)
      do.call(cbind, lapply(seq.int(2L, length.out = spec$n_classes - 1L),
                            function(k) spec$w * residual[, k]))
    },
    group = {
      residual <- expectation$group_tau -
        rowSums(expectation$group_tau) * exp(expectation$log_prior)
      do.call(cbind, lapply(seq.int(2L, length.out = spec$n_classes - 1L),
                            function(k) spec$w * residual[, k]))
    },
    "two-level" = {
      rho <- expectation$rho
      group_residual <- rho - rowSums(rho) * exp(expectation$log_eta)
      delta_scores <- do.call(cbind, lapply(
        seq.int(2L, length.out = spec$n_group_classes - 1L),
        function(h) spec$v * group_residual[, h]))
      # Row-level class-logit scores, summed to groups below.
      residual_by_group_class <- lapply(
        seq_len(spec$n_group_classes), function(h) {
          expectation$tau_by_group_class[[h]] -
            rho[spec$group_index, h] *
            exp(expectation$log_prior_by_group_class[[h]])
        })
      logit_scores <- do.call(cbind, lapply(
        seq.int(2L, length.out = spec$n_classes - 1L), function(k) {
          intercepts <- vapply(residual_by_group_class, function(r) r[, k],
                               numeric(spec$n))
          intercepts <- matrix(intercepts, spec$n)
          slopes <- spec$w * Reduce(`+`, lapply(residual_by_group_class,
                                                function(r) r[, k]))
          cbind(intercepts, slopes)
        }))
      cbind(delta_scores %||% matrix(0, spec$n_groups, 0L),
            to_units(logit_scores %||% matrix(0, spec$n, 0L)))
    })
  mixing <- mixing %||% matrix(0, if (identical(spec$nesting, "observation"))
    spec$n else spec$n_groups, 0L)
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
  theta <- .mixture_pack(spec, params)
  columns <- lapply(seq_along(theta), function(j) {
    h <- step * (1 + abs(theta[j]))
    up <- theta
    down <- theta
    up[j] <- up[j] + h
    down[j] <- down[j] - h
    (.mixture_total_score(spec, up, params) -
       .mixture_total_score(spec, down, params)) / (2 * h)
  })
  hessian <- do.call(cbind, columns)
  information <- -(hessian + t(hessian)) / 2
  dimnames(information) <- list(names(theta), names(theta))
  information
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
  inverse <- .mixture_invert(information)
  vcov <- if (identical(vcov_type, "robust")) {
    inverse %*% meat %*% inverse
  } else inverse
  dimnames(vcov) <- list(names(theta), names(theta))
  list(theta = theta, vcov = vcov, vcov_type = vcov_type,
       clustered_on = if (identical(vcov_type, "robust")) clustered_on else NA,
       information = information,
       score_norm = max(abs(colSums(scores))))
}

#' Invert an information matrix, refusing a singular one
#' @noRd
.mixture_invert <- function(information) {
  eigen_values <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eigen_values)) || min(eigen_values) <=
      1e-10 * max(abs(eigen_values))) {
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
  estimate <- transform(theta)
  jacobian <- do.call(cbind, lapply(seq_along(theta), function(j) {
    h <- step * (1 + abs(theta[j]))
    up <- theta
    down <- theta
    up[j] <- up[j] + h
    down[j] <- down[j] - h
    (transform(up) - transform(down)) / (2 * h)
  }))
  jacobian <- matrix(jacobian, length(estimate))
  variance <- rowSums((jacobian %*% vcov) * jacobian)
  list(estimate = estimate, std_error = sqrt(pmax(variance, 0)))
}

#' Wald table from estimates and standard errors
#' @noRd
.mixture_wald <- function(estimate, std_error, level) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  statistic <- estimate / std_error
  data.frame(estimate = estimate, std_error = std_error,
             statistic = statistic,
             p_value = 2 * stats::pnorm(-abs(statistic)),
             conf_low = estimate - z * std_error,
             conf_high = estimate + z * std_error,
             row.names = NULL)
}
