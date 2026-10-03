# Estimation scale, analytic scores and the quasi-Newton finish of the growth
# mixture model. The estimation vector is: class-specific coefficients (by
# class), common coefficients, log residual standard deviations (one per
# class, or one), the log-Cholesky factor of each random-effect covariance
# (log diagonal, then lower off-diagonal by column; one matrix per class, or
# one), log class scales of the proportional covariance (all but the last
# class), and the membership logits of classes 2..K.

#' Names of the estimation vector of a growth mixture model
#' @noRd
.growth_names <- function(spec) {
  classes <- paste0("class_", seq_len(spec$n_classes))
  random_names <- colnames(spec$random_design)
  n_random <- length(random_names)
  cholesky <- function(owner) {
    diagonal <- sprintf("random.log_cholesky.%s.%s:%s", owner, random_names, random_names)
    if (isTRUE(spec$random_diagonal) || n_random == 1L) return(diagonal)
    pairs <- which(lower.tri(diag(n_random)), arr.ind = TRUE)
    c(diagonal, sprintf("random.cholesky.%s.%s:%s", owner,
                        random_names[pairs[, 1L]], random_names[pairs[, 2L]]))
  }
  matrices <- switch(spec$random_covariance, varying = classes, "shared")
  c(as.vector(outer(colnames(spec$x), classes,
                    function(term, class) sprintf("coefficient.%s.%s", class, term))),
    if (ncol(spec$z) > 0L) sprintf("coefficient.common.%s", colnames(spec$z)),
    if (identical(spec$variance, "equal")) "log_sigma.shared" else
      sprintf("log_sigma.%s", classes),
    unlist(lapply(matrices, cholesky)),
    if (identical(spec$random_covariance, "proportional"))
      sprintf("random.log_scale.%s", classes[-spec$n_classes]),
    if (spec$n_classes > 1L)
      as.vector(outer(colnames(spec$w), classes[-1L],
                      function(term, class) sprintf("membership.%s.%s", class, term))))
}

#' Pack a growth-mixture parameter list into the estimation vector
#' @noRd
.growth_pack <- function(spec, params) {
  n_random <- ncol(spec$random_design)
  chart <- function(covariance) {
    factor <- .growth_factor(covariance)
    diagonal <- log(diag(factor))
    if (isTRUE(spec$random_diagonal) || n_random == 1L) return(diagonal)
    c(diagonal, factor[lower.tri(factor)])
  }
  sigma <- 0.5 * log(params$sigma2)
  theta <- c(as.vector(params$beta), params$common,
             if (identical(spec$variance, "equal")) sigma[1L] else sigma,
             unlist(lapply(params$random_covariance, chart)),
             if (identical(spec$random_covariance, "proportional"))
               log(params$random_scale[-spec$n_classes]),
             if (spec$n_classes > 1L) as.vector(params$gamma[, -1L, drop = FALSE]))
  stats::setNames(theta, .growth_names(spec))
}

#' Unpack an estimation vector into a growth-mixture parameter list
#' @param template A parameter list with the right shapes and names.
#' @noRd
.growth_unpack <- function(spec, theta, template) {
  n_classes <- spec$n_classes
  n_random <- ncol(spec$random_design)
  position <- 0L
  take <- function(count) {
    values <- theta[position + seq_len(count)]
    position <<- position + count
    unname(values)
  }
  params <- template
  params$beta[] <- take(length(params$beta))
  params$common[] <- take(length(params$common))
  sigma <- take(if (identical(spec$variance, "equal")) 1L else n_classes)
  params$sigma2 <- rep_len(exp(2 * sigma), n_classes)
  n_matrices <- if (identical(spec$random_covariance, "varying")) n_classes else 1L
  per_matrix <- if (isTRUE(spec$random_diagonal) || n_random == 1L) n_random else
    n_random * (n_random + 1L) / 2L
  params$random_covariance <- lapply(seq_len(n_matrices), function(m) {
    values <- take(per_matrix)
    factor <- diag(exp(values[seq_len(n_random)]), n_random)
    if (per_matrix > n_random) factor[lower.tri(factor)] <- values[-seq_len(n_random)]
    tcrossprod(factor)
  })
  params$random_scale <- if (identical(spec$random_covariance, "proportional"))
    c(exp(take(n_classes - 1L)), 1) else rep(1, n_classes)
  if (n_classes > 1L) {
    params$gamma[, -1L] <- take(ncol(spec$w) * (n_classes - 1L))
  }
  stopifnot("the estimation vector must be consumed exactly" =
              position == length(theta))
  params
}

#' Per-person scores of the growth mixture log likelihood
#'
#' Analytic derivatives of each person's log likelihood on the estimation
#' scale. Within class k, with r = y - X beta, V^-1 r = (r - Z m) / sigma^2
#' (m the posterior mean of b), so
#'   d/d beta = X'(r - Z m) / sigma^2,
#'   d/d sigma^2 = (|r - Z m|^2 / sigma^2 - n + tr(Z'Z S) / sigma^2) / (2 sigma^2),
#'   d/d G = (q q' - Z'V^-1 Z) / 2, q = Z'V^-1 r,
#' with Z'V^-1 Z = (Z'Z - Z'Z S Z'Z / sigma^2) / sigma^2 and S the posterior
#' covariance. The mixture weights each class's derivative by the person's
#' posterior.
#'
#' @return A persons-by-parameters matrix (unweighted by sampling weights).
#' @noRd
.growth_scores <- function(spec, stats, params, expectation = NULL) {
  expectation <- expectation %||% .growth_expectation(spec, stats, params)
  n_classes <- spec$n_classes
  n_specific <- ncol(spec$x)
  n_common <- ncol(spec$z)
  n_columns <- n_specific + n_common
  n_random <- ncol(spec$random_design)
  n_persons <- length(stats$n)
  posterior <- expectation$posterior
  covariances <- .growth_covariances(spec, params)
  zz_flat <- matrix(stats$zz, n_random^2, n_persons)
  # Per class: derivatives of log f_k for every person.
  per_class <- lapply(seq_len(n_classes), function(k) {
    sigma2 <- params$sigma2[k]
    coefficients <- .growth_coefficients(params, k)
    residuals <- .growth_residuals(stats, coefficients)
    means <- expectation$means[[k]]
    pattern_covariance <- matrix(expectation$covariances[[k]], n_random^2)[, stats$pattern,
                                                                        drop = FALSE]
    # X'r = X'y - X'X beta, X'Z m and Z'Z m, persons by columns.
    by_column <- function(array3, rows, vector_or_matrix) {
      Reduce(`+`, lapply(seq_along(rows), function(b) {
        slice <- matrix(array3[, b, ], dim(array3)[1L], n_persons)
        weight <- if (is.matrix(vector_or_matrix)) vector_or_matrix[, b] else
          vector_or_matrix[b]
        t(slice) * weight
      }))
    }
    xr <- stats$xy - by_column(stats$xx, seq_len(n_columns), coefficients)
    xzm <- by_column(stats$xz, seq_len(n_random), means)
    zzm <- by_column(stats$zz, seq_len(n_random), means)
    d_beta <- (xr - xzm) / sigma2
    residual_ss <- residuals$rr - 2 * rowSums(means * residuals$zr) + rowSums(means * zzm)
    trace_term <- colSums(zz_flat * pattern_covariance)
    d_sigma2 <- (residual_ss / sigma2 - stats$n + trace_term / sigma2) / (2 * sigma2)
    # q = Z'V^-1 r = (Z'r - Z'Z m) / sigma^2; Z'V^-1 Z depends on the pattern only.
    q <- (residuals$zr - zzm) / sigma2
    zvz <- vapply(seq_len(dim(stats$pattern_zz)[3L]), function(g) {
      zz <- matrix(stats$pattern_zz[, , g], n_random)
      s <- matrix(expectation$covariances[[k]][, , g], n_random)
      as.vector((zz - zz %*% s %*% zz / sigma2) / sigma2)
    }, numeric(n_random^2)) |> matrix(n_random^2)
    outer_q <- t(q[, rep(seq_len(n_random), n_random), drop = FALSE] *
                   q[, rep(seq_len(n_random), each = n_random), drop = FALSE])
    d_g <- (outer_q - zvz[, stats$pattern, drop = FALSE]) / 2
    list(beta = d_beta, sigma2 = d_sigma2, g = d_g, covariance = covariances[[k]])
  })
  # Coefficients.
  beta_scores <- do.call(cbind, lapply(seq_len(n_classes), function(k) {
    per_class[[k]]$beta[, seq_len(n_specific), drop = FALSE] * posterior[, k]
  }))
  common_scores <- if (n_common > 0L) {
    Reduce(`+`, lapply(seq_len(n_classes), function(k) {
      per_class[[k]]$beta[, n_specific + seq_len(n_common), drop = FALSE] * posterior[, k]
    }))
  }
  # log sigma: d/d log sigma = 2 sigma^2 d/d sigma^2.
  sigma_by_class <- vapply(seq_len(n_classes), function(k) {
    2 * params$sigma2[k] * per_class[[k]]$sigma2 * posterior[, k]
  }, numeric(n_persons)) |> matrix(n_persons, n_classes)
  sigma_scores <- if (identical(spec$variance, "equal")) {
    matrix(rowSums(sigma_by_class), n_persons, 1L)
  } else sigma_by_class
  # Random-effect covariance, through its chart.
  chart_scores <- function(d_g_flat, covariance) {
    factor <- .growth_factor(covariance)
    t(apply(d_g_flat, 2L, function(column) {
      d <- matrix(column, n_random)
      d <- (d + t(d)) / 2
      d_factor <- 2 * d %*% factor
      diagonal <- diag(d_factor) * diag(factor)
      if (isTRUE(spec$random_diagonal) || n_random == 1L) return(diagonal)
      c(diagonal, d_factor[lower.tri(d_factor)])
    })) |> matrix(n_persons)
  }
  class_g <- lapply(seq_len(n_classes), function(k) {
    per_class[[k]]$g * rep(posterior[, k], each = n_random^2)
  })
  random_scores <- switch(spec$random_covariance,
    varying = do.call(cbind, lapply(seq_len(n_classes), function(k) {
      chart_scores(class_g[[k]], covariances[[k]])
    })),
    equal = chart_scores(Reduce(`+`, class_g), covariances[[1L]]),
    proportional = {
      shared <- params$random_covariance[[1L]]
      scale2 <- params$random_scale^2
      base <- chart_scores(Reduce(`+`, Map(`*`, class_g, scale2)), shared)
      scale_scores <- vapply(seq_len(n_classes - 1L), function(k) {
        2 * colSums(class_g[[k]] * as.vector(covariances[[k]]))
      }, numeric(n_persons)) |> matrix(n_persons)
      cbind(base, scale_scores)
    })
  # Membership logits of classes 2..K.
  prior <- exp(.mixture_log_softmax(spec$w, params$gamma))
  membership_scores <- if (n_classes > 1L) {
    do.call(cbind, lapply(seq_len(n_classes)[-1L], function(k) {
      spec$w * (posterior[, k] - prior[, k])
    }))
  }
  scores <- cbind(beta_scores, common_scores, sigma_scores, random_scores,
                  membership_scores)
  dimnames(scores) <- list(NULL, .growth_names(spec))
  scores
}

#' Quasi-Newton finish of the growth mixture likelihood
#'
#' EM is reliable from far away and slow near the maximum; L-BFGS-B with the
#' analytic gradient finishes from EM's point. Variances keep their floors as
#' bounds. Convergence is judged by the relative change of the log
#' likelihood; L-BFGS-B's line-search stop (code 52) at the maximum is
#' checked by one restart, as for transition fits.
#'
#' @return A list with `params`, `converged` and `log_likelihood`.
#' @noRd
.growth_quasi_newton <- function(spec, stats, params, tol) {
  theta <- .growth_pack(spec, params)
  names_all <- names(theta)
  lower <- rep(-Inf, length(theta))
  lower[startsWith(names_all, "log_sigma.")] <- 0.5 * log(spec$min_variance)
  lower[startsWith(names_all, "random.log_cholesky.")] <- 0.5 * log(spec$min_variance)
  theta <- pmax(theta, lower)
  weights <- spec$sampling_weights %||% rep(1, length(stats$n))
  value <- function(v) {
    -.growth_expectation(spec, stats, .growth_unpack(spec, v, params))$log_likelihood
  }
  gradient <- function(v) {
    current <- .growth_unpack(spec, v, params)
    -colSums(.growth_scores(spec, stats, current) * weights)
  }
  # Polished to near machine precision whatever `tol` is: with the analytic
  # gradient it is cheap, and `tol` judges convergence below.
  control <- list(maxit = 2000L, factr = 10, pgtol = 0)
  start_value <- value(theta)
  found <- stats::optim(theta, value, gradient, method = "L-BFGS-B", lower = lower,
                        control = control)
  if (!is.finite(found$value) || found$value > start_value) {
    return(list(params = params, converged = FALSE, log_likelihood = -start_value))
  }
  converged <- found$convergence == 0L
  if (identical(found$convergence, 52L)) {
    # The line search stops at the maximum when rounding hides any ascent.
    again <- stats::optim(found$par, value, gradient, method = "L-BFGS-B",
                          lower = lower, control = control)
    if (is.finite(again$value) && again$value <= found$value) {
      converged <- again$convergence == 0L ||
        found$value - again$value <= tol * (1 + abs(again$value))
      found <- again
    }
  }
  list(params = .growth_unpack(spec, found$par, params), converged = converged,
       log_likelihood = -found$value)
}
