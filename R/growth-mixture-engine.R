# Engine of the Gaussian growth mixture model: mixture_regression() with
# `random =`. Within class k, person i's outcomes are
#   y_i = X_i beta_k + C_i beta_c + Z_i b_i + e_i,
#   b_i ~ N(0, G_k), e_i ~ N(0, sigma_k^2 I),
# so marginally y_i | k ~ N(X_i beta_k + C_i beta_c, Z_i G_k Z_i' + sigma_k^2 I)
# (Verbeke and Lesaffre 1996; Muthen and Shedden 1999; Proust-Lima, Philipps
# and Liquet 2017). The likelihood is evaluated per person in the q x q
# random-effect space: with G_k = L L' and A = I + L' Z'Z L / sigma^2,
#   log|V| = n log sigma^2 + log|A|,
#   r'V^-1 r = (r'r - u'A^-1 u / sigma^2) / sigma^2,   u = L' Z'r,
# so G is never inverted and a variance on its boundary stays finite. Every
# quantity comes from per-person cross-products computed once.

#' Per-person cross-products of the growth mixture model
#'
#' @param spec The model specification, with `random_design`.
#' @return A list of per-person arrays: `n` (rows), `yy`, `xy` (persons by
#'   regression columns), `xx` (columns by columns by persons), `zy`
#'   (persons by random effects), `zz` (random by random by persons) and
#'   `xz` (columns by random by persons). Regression columns are the
#'   class-specific design followed by the common design.
#' @noRd
.growth_statistics <- function(spec) {
  design <- cbind(spec$x, spec$z)
  random <- spec$random_design
  y <- spec$y - spec$offset
  rows <- split(seq_len(spec$n), spec$group_index)
  n_columns <- ncol(design)
  n_random <- ncol(random)
  cross <- function(a, b, i) crossprod(a[i, , drop = FALSE], b[i, , drop = FALSE])
  list(
    n = lengths(rows, use.names = FALSE),
    yy = vapply(rows, function(i) sum(y[i]^2), numeric(1), USE.NAMES = FALSE),
    xy = t(vapply(rows, function(i) as.vector(crossprod(design[i, , drop = FALSE], y[i])),
                  numeric(n_columns), USE.NAMES = FALSE)) |>
      matrix(length(rows), n_columns),
    zy = t(vapply(rows, function(i) as.vector(crossprod(random[i, , drop = FALSE], y[i])),
                  numeric(n_random), USE.NAMES = FALSE)) |>
      matrix(length(rows), n_random),
    xx = array(vapply(rows, function(i) cross(design, design, i),
                      matrix(0, n_columns, n_columns), USE.NAMES = FALSE),
               c(n_columns, n_columns, length(rows))),
    zz = array(vapply(rows, function(i) cross(random, random, i),
                      matrix(0, n_random, n_random), USE.NAMES = FALSE),
               c(n_random, n_random, length(rows))),
    xz = array(vapply(rows, function(i) cross(design, random, i),
                      matrix(0, n_columns, n_random), USE.NAMES = FALSE),
               c(n_columns, n_random, length(rows)))) |>
    .growth_patterns()
}

#' Group persons by their random-effect design
#'
#' Two persons with the same Z'Z have the same posterior covariance and the
#' same determinant in every class, so those are computed once per pattern.
#'
#' @param stats Per-person cross-products.
#' @return `stats` with `pattern` (each person's pattern index) and
#'   `pattern_zz` (random by random by patterns).
#' @noRd
.growth_patterns <- function(stats) {
  n_random <- dim(stats$zz)[1L]
  flat <- matrix(stats$zz, n_random^2, length(stats$n))
  key <- apply(signif(flat, 12), 2L, paste, collapse = "\r")
  first <- !duplicated(key)
  stats$pattern <- match(key, key[first])
  stats$pattern_zz <- array(flat[, first, drop = FALSE],
                            c(n_random, n_random, sum(first)))
  stats
}

#' Regression coefficients of one class, class-specific then common
#' @noRd
.growth_coefficients <- function(params, k) c(params$beta[, k], params$common)

#' Each person's residual cross-products for one class's coefficients
#'
#' @param stats Output of `.growth_statistics()`.
#' @param coefficients Regression coefficients (class-specific, then common).
#' @return A list with `rr` (per person r'r) and `zr` (persons by random
#'   effects: Z'r), where r = y - X beta.
#' @noRd
.growth_residuals <- function(stats, coefficients) {
  n_persons <- length(stats$n)
  n_columns <- length(coefficients)
  n_random <- dim(stats$zz)[1L]
  outer_coefficients <- as.vector(tcrossprod(coefficients))
  rr <- stats$yy - 2 * as.vector(stats$xy %*% coefficients) +
    colSums(matrix(stats$xx, n_columns^2, n_persons) * outer_coefficients)
  # X'Z is columns x random x persons; beta'X'Z is random x persons.
  xz_beta <- matrix(crossprod(coefficients, matrix(stats$xz, n_columns,
                                                   n_random * n_persons)),
                    n_random, n_persons)
  list(rr = rr, zr = stats$zy - t(xz_beta))
}

#' One class's per-person log densities and random-effect posteriors
#'
#' Persons who share a design (the same Z'Z, as in a balanced or
#' near-balanced study) share the q x q factorization, which is done once
#' per pattern; everything per person is then a vector operation.
#'
#' @param stats Output of `.growth_statistics()`.
#' @param coefficients The class's regression coefficients (class-specific,
#'   then common).
#' @param factor Lower Cholesky factor `L` of the class's random-effect
#'   covariance.
#' @param sigma2 The class's residual variance.
#' @return A list with `log_density` (one per person), `mean` (persons by
#'   random effects: posterior means of b_i), `covariance` (random by random
#'   by patterns: posterior covariances) and the per-person `zr` and `rr`.
#' @noRd
.growth_class_terms <- function(stats, coefficients, factor, sigma2) {
  n_random <- ncol(factor)
  residuals <- .growth_residuals(stats, coefficients)
  per_pattern <- lapply(seq_len(dim(stats$pattern_zz)[3L]), function(g) {
    zz <- matrix(stats$pattern_zz[, , g], n_random, n_random)
    a_factor <- chol(diag(n_random) + crossprod(factor, zz %*% factor) / sigma2)
    a_inverse <- chol2inv(a_factor)
    list(log_det = 2 * sum(log(diag(a_factor))), a_inverse = a_inverse,
         covariance = factor %*% a_inverse %*% t(factor))
  })
  pattern <- stats$pattern
  u <- residuals$zr %*% factor
  # A^-1 u, person by person through each person's pattern.
  solved <- matrix(0, nrow(u), n_random)
  invisible(lapply(seq_along(per_pattern), function(g) {
    members <- which(pattern == g)
    solved[members, ] <<- u[members, , drop = FALSE] %*% per_pattern[[g]]$a_inverse
  }))
  log_det <- stats$n * log(sigma2) +
    vapply(per_pattern, `[[`, numeric(1), "log_det")[pattern]
  quadratic <- (residuals$rr - rowSums(u * solved) / sigma2) / sigma2
  list(log_density = -0.5 * (stats$n * log(2 * pi) + log_det + quadratic),
       mean = tcrossprod(solved, factor) / sigma2,
       covariance = array(vapply(per_pattern, `[[`, matrix(0, n_random, n_random),
                                 "covariance"),
                          c(n_random, n_random, length(per_pattern))),
       zr = residuals$zr, rr = residuals$rr)
}

#' Random-effect covariance of every class from the parameter list
#' @return A list of q x q matrices, one per class.
#' @noRd
.growth_covariances <- function(spec, params) {
  lapply(seq_len(spec$n_classes), function(k) {
    base <- params$random_covariance[[if (identical(spec$random_covariance, "varying")) k else 1L]]
    if (identical(spec$random_covariance, "proportional")) {
      base * params$random_scale[k]^2
    } else base
  })
}

#' A Cholesky factor that tolerates a covariance on its boundary
#' @param covariance A symmetric positive semi-definite matrix.
#' @return A lower-triangular `L` with `L %*% t(L)` equal to `covariance` up
#'   to a jitter of 1e-12 times its scale on a singular direction.
#' @noRd
.growth_factor <- function(covariance) {
  jitter <- 1e-12 * max(1, max(abs(diag(covariance))))
  t(chol(covariance + diag(jitter, nrow(covariance))))
}

#' E-step of the growth mixture model
#'
#' @param spec The model specification.
#' @param stats Per-person cross-products.
#' @param params The parameter list.
#' @return A list with `log_likelihood` (sampling-weighted when the fit is),
#'   `tau` (persons by classes, weighted posteriors used by the M-step),
#'   `posterior` (per-person posteriors, rows summing to one), `log_density`
#'   (persons by classes), and the per-class random-effect posteriors
#'   `means` and `covariances` (lists by class).
#' @noRd
.growth_expectation <- function(spec, stats, params) {
  covariances <- .growth_covariances(spec, params)
  terms <- lapply(seq_len(spec$n_classes), function(k) {
    .growth_class_terms(stats, .growth_coefficients(params, k),
                        .growth_factor(covariances[[k]]), params$sigma2[k])
  })
  log_density <- vapply(terms, `[[`, numeric(length(stats$n)), "log_density") |>
    matrix(length(stats$n), spec$n_classes)
  structure <- .growth_structure(spec)
  persons <- .mixture_structure_expectation(structure, params, log_density)
  weighted <- if (is.null(structure$sampling_weights)) persons else
    .mixture_weigh(structure, persons)
  # `structure` keeps the (weighted) structure E-step the membership M-step
  # and scores read: posteriors joint with each group class for two levels.
  list(log_likelihood = weighted$log_likelihood,
       unit_log_likelihood = persons$unit_log_likelihood,
       tau = weighted$tau, posterior = persons$tau, structure = weighted,
       log_density = log_density,
       means = lapply(terms, `[[`, "mean"),
       covariances = lapply(terms, `[[`, "covariance"))
}

#' M-step of the growth mixture model (one conditional maximization each)
#'
#' Regression coefficients given the variances, then residual variances,
#' random-effect covariances and class shares: each step maximizes the
#' expected complete-data log likelihood with the random effects as missing
#' data (Laird and Ware 1982), so the likelihood never decreases.
#'
#' @return The updated parameter list.
#' @noRd
.growth_maximization <- function(spec, stats, params, expectation) {
  n_classes <- spec$n_classes
  n_specific <- ncol(spec$x)
  n_common <- ncol(spec$z)
  n_columns <- n_specific + n_common
  n_random <- ncol(spec$random_design)
  n_persons <- length(stats$n)
  tau <- expectation$tau
  xx_flat <- matrix(stats$xx, n_columns^2, n_persons)
  zz_flat <- matrix(stats$zz, n_random^2, n_persons)
  # X_i'Z_i m_i for every person (columns by persons).
  xz_times <- function(means) {
    weighted <- array(matrix(stats$xz, n_columns, n_random * n_persons) *
                        rep(as.vector(t(means)), each = n_columns),
                      c(n_columns, n_random, n_persons))
    matrix(colSums(matrix(aperm(weighted, c(2L, 1L, 3L)), n_random)),
           n_columns, n_persons)
  }
  # Each person's m_i m_i', flattened like Z'Z (random^2 by persons).
  outer_means <- function(means) {
    t(means[, rep(seq_len(n_random), n_random), drop = FALSE] *
        means[, rep(seq_len(n_random), each = n_random), drop = FALSE])
  }
  # Weighted normal equations of every class, then the joint solve for the
  # class-specific and shared coefficients.
  blocks <- lapply(seq_len(n_classes), function(k) {
    weight <- tau[, k] / params$sigma2[k]
    list(gram = matrix(xx_flat %*% weight, n_columns, n_columns),
         right = as.vector(crossprod(stats$xy, weight) -
                             xz_times(expectation$means[[k]]) %*% weight))
  })
  specific <- seq_len(n_specific)
  shared <- n_specific + seq_len(n_common)
  size <- n_classes * n_specific + n_common
  common_at <- n_classes * n_specific + seq_len(n_common)
  system <- matrix(0, size, size)
  right <- numeric(size)
  invisible(lapply(seq_len(n_classes), function(k) {
    at <- (k - 1L) * n_specific + specific
    system[at, at] <<- blocks[[k]]$gram[specific, specific]
    right[at] <<- blocks[[k]]$right[specific]
    if (n_common > 0L) {
      system[at, common_at] <<- blocks[[k]]$gram[specific, shared]
      system[common_at, at] <<- blocks[[k]]$gram[shared, specific]
      system[common_at, common_at] <<- system[common_at, common_at] +
        blocks[[k]]$gram[shared, shared]
      right[common_at] <<- right[common_at] + blocks[[k]]$right[shared]
    }
  }))
  solution <- .mixture_solve(system, right)
  params$beta[] <- solution[seq_len(n_classes * n_specific)]
  if (n_common > 0L) params$common[] <- solution[common_at]

  # Expected residual sums of squares, per class, at the new coefficients:
  # E||y - X beta - Z b||^2 = r'r - 2 m'Z'r + m'Z'Z m + tr(Z'Z S).
  residual <- vapply(seq_len(n_classes), function(k) {
    residuals <- .growth_residuals(stats, .growth_coefficients(params, k))
    means <- expectation$means[[k]]
    covariance_flat <- matrix(expectation$covariances[[k]], n_random^2)
    per_person <- residuals$rr - 2 * rowSums(means * residuals$zr) +
      colSums(zz_flat * outer_means(means)) +
      colSums(zz_flat * covariance_flat[, stats$pattern, drop = FALSE])
    sum(tau[, k] * per_person)
  }, numeric(1))
  rows <- colSums(tau * stats$n)
  params$sigma2 <- if (identical(spec$variance, "equal")) {
    rep(sum(residual) / sum(rows), n_classes)
  } else residual / rows
  params$sigma2 <- pmax(params$sigma2, spec$min_variance)

  # Expected second moments of the random effects, per class.
  moments <- lapply(seq_len(n_classes), function(k) {
    means <- expectation$means[[k]]
    by_pattern <- as.vector(rowsum(tau[, k], stats$pattern, reorder = TRUE))
    crossprod(means * tau[, k], means) +
      matrix(matrix(expectation$covariances[[k]], n_random^2) %*% by_pattern,
             n_random, n_random)
  })
  sizes <- pmax(colSums(tau), .Machine$double.xmin)
  restrict <- function(covariance) {
    covariance <- (covariance + t(covariance)) / 2
    if (isTRUE(spec$random_diagonal)) covariance <- diag(diag(covariance), n_random)
    covariance
  }
  params$random_covariance <- switch(spec$random_covariance,
    varying = lapply(seq_len(n_classes), function(k) restrict(moments[[k]] / sizes[k])),
    equal = list(restrict(Reduce(`+`, moments) / sum(sizes))),
    proportional = {
      # Alternate the closed forms once: the shared matrix given the scales,
      # then each class's scale given the matrix (the last class's is 1).
      scale2 <- params$random_scale^2
      shared_matrix <- restrict(Reduce(`+`, Map(`/`, moments, scale2)) / sum(sizes))
      inverse <- .mixture_solve(shared_matrix, diag(n_random))
      updated <- vapply(seq_len(n_classes), function(k) {
        sum(inverse * moments[[k]]) / (n_random * sizes[k])
      }, numeric(1))
      # Re-express with the last class as the unit of scale; every class's
      # covariance s_k^2 G is unchanged by moving that factor into G.
      params$random_scale <- sqrt(pmax(updated / updated[n_classes], 1e-12))
      list(shared_matrix * updated[n_classes])
    })
  params$random_covariance <- lapply(params$random_covariance, function(covariance) {
    diag(covariance) <- pmax(diag(covariance), spec$min_variance)
    covariance
  })
  .mixture_update_mixing(.growth_structure(spec), params, expectation$structure)
}

#' The latent structure of a growth mixture: classes over persons
#'
#' Each person's random effects are integrated out of their density, so the
#' persons are the units the classes belong to. Without clusters every person
#' has their own class, with membership logits of person covariates -- the
#' `"observation"` structure with persons for rows. With `cluster`, persons
#' are nested in clusters (students in schools) and a group class of the
#' cluster shifts how probable each trajectory class is -- the `"two-level"`
#' structure, whose independent units are the clusters: the multilevel growth
#' mixture model (Asparouhov and Muthen 2008; Vermunt 2003).
#'
#' @param spec The growth specification.
#' @return A structure list (see kernel-structures.R).
#' @noRd
.growth_structure <- function(spec) {
  if (is.null(spec$cluster)) {
    return(list(nesting = "observation", n_classes = spec$n_classes, w = spec$w,
                intercept_only = spec$intercept_only,
                sampling_weights = spec$sampling_weights,
                row_weights = spec$sampling_weights))
  }
  list(nesting = "two-level", n_classes = spec$n_classes, n = spec$n_groups,
       w = spec$cluster_w, n_group_classes = spec$n_cluster_classes,
       v = spec$cluster_v, group_intercept_only = spec$cluster_intercept_only,
       group_index = spec$cluster_index, n_groups = spec$n_clusters)
}

#' Names of a structure's free membership coefficients
#' @noRd
.growth_structure_names <- function(structure) {
  k <- structure$n_classes
  zero <- if (identical(structure$nesting, "two-level")) {
    list(delta = matrix(0, ncol(structure$v), structure$n_group_classes),
         class_logits = matrix(0, structure$n_group_classes + ncol(structure$w), k))
  } else list(gamma = matrix(0, ncol(structure$w), k))
  names(.mixture_structure_pack(structure, zero))
}

#' EM for the growth mixture model from one start
#' @return A list with `params`, `expectation`, `converged`, `iterations`,
#'   `history` and `degenerate`.
#' @noRd
.growth_em <- function(spec, stats, params, max_iter, tol) {
  fit <- .latents_em(params, function(point) .growth_expectation(spec, stats, point),
                     function(point, expectation) {
                       .growth_maximization(spec, stats, point, expectation)
                     },
                     max_iter, tol, .mixture_em_settings())
  params <- fit$state
  expectation <- fit$expectation
  sizes <- colSums(expectation$posterior)
  degenerate <- any(sizes < ncol(spec$x) + 1) ||
    any(params$sigma2 <= spec$min_variance * (1 + 1e-8))
  list(params = params, expectation = expectation, converged = fit$converged,
       iterations = fit$iterations, history = fit$history, degenerate = degenerate)
}
