# An independent reference for covariances on chosen pairs of continuous
# indicators, written from the model's formula. Nothing here calls latents:
#
#   f(x) = sum_k pi_k N(x; mu_k, Sigma_k),
#
# Sigma_k block diagonal: one 2 x 2 block per pair, a variance for every other
# indicator, every other covariance exactly zero ("equal": one Sigma shared
# by the classes). A row with missing values contributes the marginal density
# of its observed entries (FIML). A floor `min_variance` bounds every
# eigenvalue of every block from below.
#
# The maximizer is an EM of its own (conditional moments for missing values,
# masked M-step, eigenvalue clipping per block at a floor) from several
# starts of its own, polished by BFGS on a chart that cannot leave the
# model: Sigma_block = floor * I + L L', L lower triangular with log
# diagonal. Standard errors come from numDeriv::hessian of the log
# likelihood in the raw coordinates (class proportions 2..K, means,
# variances, pair covariances).

#' Indicator blocks: each pair, then every unpaired indicator alone
ref_blocks <- function(vars, pairs) {
  split <- lapply(strsplit(pairs, ":", fixed = TRUE), \(p) match(p, vars))
  paired <- unlist(split)
  c(split, as.list(setdiff(seq_along(vars), paired)))
}

#' The free entries of a covariance matrix (diagonal and the pairs)
ref_pattern <- function(vars, pairs) {
  pattern <- diag(length(vars)) > 0
  invisible(lapply(ref_blocks(vars, pairs), function(block) {
    pattern[block, block] <<- TRUE
  }))
  dimnames(pattern) <- list(vars, vars)
  pattern
}

#' Log likelihood of a mixture at `par` (proportions, means k x d, covariances
#' d x d x k), each row on its observed entries; optionally per row.
ref_loglik <- function(par, x, per_row = FALSE, weights = NULL) {
  k <- length(par$proportions)
  observed <- !is.na(x)
  key <- as.vector(observed %*% 2^(seq_len(ncol(x)) - 1L))
  rows <- split(seq_len(nrow(x)), key)
  values <- numeric(nrow(x))
  invisible(lapply(rows, function(index) {
    o <- which(observed[index[1L], ])
    if (length(o) == 0L) return(NULL)
    logs <- vapply(seq_len(k), function(j) {
      log(par$proportions[j]) +
        mvtnorm::dmvnorm(x[index, o, drop = FALSE], par$means[j, o],
                         par$covariances[o, o, j, drop = TRUE] |> as.matrix(), log = TRUE)
    }, numeric(length(index)))
    logs <- matrix(logs, length(index), k)
    top <- apply(logs, 1L, max)
    values[index] <<- top + log(rowSums(exp(logs - top)))
    NULL
  }))
  if (per_row) return(values)
  sum(values * (weights %||% 1))
}

#' Posterior class probabilities at `par`
ref_posteriors <- function(par, x) {
  k <- length(par$proportions)
  observed <- !is.na(x)
  key <- as.vector(observed %*% 2^(seq_len(ncol(x)) - 1L))
  rows <- split(seq_len(nrow(x)), key)
  out <- matrix(NA_real_, nrow(x), k)
  invisible(lapply(rows, function(index) {
    o <- which(observed[index[1L], ])
    logs <- vapply(seq_len(k), function(j) {
      if (length(o) == 0L) return(rep(log(par$proportions[j]), length(index)))
      log(par$proportions[j]) +
        mvtnorm::dmvnorm(x[index, o, drop = FALSE], par$means[j, o],
                         as.matrix(par$covariances[o, o, j]), log = TRUE)
    }, numeric(length(index)))
    logs <- matrix(logs, length(index), k)
    top <- apply(logs, 1L, max)
    out[index, ] <<- exp(logs - top) / rowSums(exp(logs - top))
    NULL
  }))
  out
}

#' Clip the eigenvalues of each block at the floor (the M-step under an
#' eigenvalue bound: Sigma = V max(Lambda, floor) V' maximizes
#' -log|Sigma| - tr(S Sigma^-1) over {Sigma : Lambda >= floor}).
ref_clip <- function(sigma, blocks, floor) {
  invisible(lapply(blocks, function(block) {
    e <- eigen(sigma[block, block, drop = FALSE], symmetric = TRUE)
    sigma[block, block] <<- e$vectors %*% diag(pmax(e$values, floor), length(block)) %*%
      t(e$vectors)
  }))
  sigma
}

#' One EM step: conditional moments of missing entries, masked M-step
ref_em_step <- function(par, x, pattern, blocks, equal, floor, weights = NULL) {
  k <- length(par$proportions)
  d <- ncol(x)
  w <- weights %||% rep(1, nrow(x))
  tau <- ref_posteriors(par, x) * w
  observed <- !is.na(x)
  key <- as.vector(observed %*% 2^(seq_len(ncol(x)) - 1L))
  rows <- split(seq_len(nrow(x)), key)
  # Expected x and the conditional covariance of the missing entries, per class.
  moments <- lapply(seq_len(k), function(j) {
    expected <- x
    extra <- matrix(0, d, d)
    invisible(lapply(rows, function(index) {
      o <- which(observed[index[1L], ])
      m <- which(!observed[index[1L], ])
      if (length(m) == 0L) return(NULL)
      sigma <- par$covariances[, , j]
      if (length(o) == 0L) {
        expected[index, ] <<- matrix(par$means[j, ], length(index), d, byrow = TRUE)
        extra <<- extra + sigma * sum(tau[index, j])
        return(NULL)
      }
      gain <- sigma[m, o, drop = FALSE] %*% solve(sigma[o, o, drop = FALSE])
      centred <- sweep(x[index, o, drop = FALSE], 2L, par$means[j, o])
      expected[index, m] <<- sweep(centred %*% t(gain), 2L, par$means[j, m], "+")
      conditional <- sigma[m, m, drop = FALSE] - gain %*% sigma[o, m, drop = FALSE]
      extra[m, m] <<- extra[m, m] + conditional * sum(tau[index, j])
      NULL
    }))
    list(expected = expected, extra = extra)
  })
  sizes <- colSums(tau)
  means <- t(vapply(seq_len(k), function(j) {
    colSums(moments[[j]]$expected * tau[, j]) / sizes[j]
  }, numeric(d)))
  scatter <- lapply(seq_len(k), function(j) {
    centred <- sweep(moments[[j]]$expected, 2L, means[j, ])
    crossprod(centred, centred * tau[, j]) + moments[[j]]$extra
  })
  covariances <- if (equal) {
    shared <- Reduce(`+`, scatter) / sum(sizes)
    shared[!pattern] <- 0
    shared <- ref_clip(shared, blocks, floor)
    array(shared, c(d, d, k))
  } else {
    array(vapply(seq_len(k), function(j) {
      s <- scatter[[j]] / sizes[j]
      s[!pattern] <- 0
      ref_clip(s, blocks, floor)
    }, numeric(d * d)), c(d, d, k))
  }
  list(proportions = sizes / sum(sizes), means = means, covariances = covariances)
}

#' EM to a relative change of 1e-12 (or `max_iter` steps)
ref_em <- function(par, x, pattern, blocks, equal, floor, max_iter = 5000L, weights = NULL) {
  previous <- ref_loglik(par, x, weights = weights)
  iteration <- 0L
  change <- Inf
  # An iterative fit: each step needs the previous one.
  while (iteration < max_iter && change > 1e-12) {
    par <- ref_em_step(par, x, pattern, blocks, equal, floor, weights)
    current <- ref_loglik(par, x, weights = weights)
    change <- abs(current - previous) / (1 + abs(current))
    previous <- current
    iteration <- iteration + 1L
  }
  list(par = par, log_likelihood = previous, iterations = iteration)
}

#' A start from a partition: class proportions, means and patterned scatter
ref_start_from_posteriors <- function(tau, x, pattern, blocks, equal, floor) {
  k <- ncol(tau)
  d <- ncol(x)
  filled <- x
  filled[is.na(filled)] <- matrix(colMeans(x, na.rm = TRUE), nrow(x), d, byrow = TRUE)[is.na(x)]
  par <- list(proportions = rep(1 / k, k), means = matrix(0, k, d),
              covariances = array(diag(apply(filled, 2L, stats::var), d), c(d, d, k)))
  # One M-step on the filled data with these memberships.
  sizes <- colSums(tau)
  means <- t(vapply(seq_len(k), \(j) colSums(filled * tau[, j]) / sizes[j], numeric(d)))
  scatter <- lapply(seq_len(k), function(j) {
    centred <- sweep(filled, 2L, means[j, ])
    crossprod(centred, centred * tau[, j]) / sizes[j]
  })
  covariances <- array(vapply(seq_len(k), function(j) {
    s <- if (equal) Reduce(`+`, lapply(seq_len(k), \(i) scatter[[i]] * sizes[i])) / sum(sizes) else
      scatter[[j]]
    s[!pattern] <- 0
    ref_clip(s + diag(1e-6 * diag(s), d), blocks, max(floor, 1e-8 * max(diag(s))))
  }, numeric(d * d)), c(d, d, k))
  list(proportions = sizes / sum(sizes), means = means, covariances = covariances)
}

# -- charts ---------------------------------------------------------------------

#' Raw coordinates: proportions 2..K, means, then per covariance matrix
#' (one if equal) each indicator's variance and each pair's covariance.
ref_pack_raw <- function(par, vars, pairs, equal) {
  k <- length(par$proportions)
  split <- lapply(strsplit(pairs, ":", fixed = TRUE), \(p) match(p, vars))
  matrices <- if (equal) 1L else seq_len(k)
  covariance <- unlist(lapply(matrices, function(j) {
    sigma <- par$covariances[, , j]
    c(diag(sigma), vapply(split, \(p) sigma[p[1L], p[2L]], numeric(1)))
  }))
  c(par$proportions[-1L], as.vector(t(par$means)), covariance)
}

ref_unpack_raw <- function(theta, vars, pairs, k, equal) {
  d <- length(vars)
  split <- lapply(strsplit(pairs, ":", fixed = TRUE), \(p) match(p, vars))
  proportions <- c(1 - sum(theta[seq_len(k - 1L)]), theta[seq_len(k - 1L)])
  means <- matrix(theta[k - 1L + seq_len(k * d)], k, d, byrow = TRUE)
  rest <- theta[-seq_len(k - 1L + k * d)]
  per <- d + length(split)
  build <- function(values) {
    sigma <- diag(values[seq_len(d)], d)
    invisible(Map(function(p, v) {
      sigma[p[1L], p[2L]] <<- v
      sigma[p[2L], p[1L]] <<- v
    }, split, values[d + seq_along(split)]))
    sigma
  }
  covariances <- if (equal) array(build(rest), c(d, d, k)) else
    array(vapply(seq_len(k), \(j) build(rest[(j - 1L) * per + seq_len(per)]),
                 numeric(d * d)), c(d, d, k))
  list(proportions = proportions, means = means, covariances = covariances)
}

#' Names of the raw coordinates, in the order of ref_pack_raw()
ref_raw_names <- function(vars, pairs, k, equal) {
  classes <- if (equal) "all" else as.character(seq_len(k))
  c(if (k > 1L) sprintf("probability:%d", 2:k),
    sprintf("mean:%d:%s", rep(seq_len(k), each = length(vars)), vars),
    unlist(lapply(classes, \(j) c(sprintf("variance:%s:%s", j, vars),
                                  sprintf("covariance:%s:%s", j, pairs)))))
}

#' Floor-respecting chart: logits, means, and per block log-Cholesky of
#' Sigma_block - floor * I.
ref_pack_chart <- function(par, blocks, equal, floor) {
  k <- length(par$proportions)
  matrices <- if (equal) 1L else seq_len(k)
  chol_part <- unlist(lapply(matrices, function(j) {
    unlist(lapply(blocks, function(block) {
      excess <- par$covariances[block, block, j, drop = FALSE] |> matrix(length(block))
      excess <- excess - diag(floor, length(block))
      # A block at the floor has a singular excess; nudge it inside.
      e <- eigen(excess, symmetric = TRUE)
      excess <- e$vectors %*% diag(pmax(e$values, 1e-10 * max(1, abs(e$values))),
                                   length(block)) %*% t(e$vectors)
      l <- t(chol(excess))
      c(log(diag(l)), l[lower.tri(l)])
    }))
  }))
  c(log(par$proportions[-1L] / par$proportions[1L]), as.vector(t(par$means)), chol_part)
}

ref_unpack_chart <- function(theta, blocks, d, k, equal, floor) {
  logits <- c(0, theta[seq_len(k - 1L)])
  proportions <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
  means <- matrix(theta[k - 1L + seq_len(k * d)], k, d, byrow = TRUE)
  rest <- theta[-seq_len(k - 1L + k * d)]
  sizes <- vapply(blocks, \(b) length(b) * (length(b) + 1L) / 2L, numeric(1))
  per <- sum(sizes)
  build <- function(values) {
    sigma <- matrix(0, d, d)
    ends <- cumsum(sizes)
    invisible(lapply(seq_along(blocks), function(i) {
      block <- blocks[[i]]
      v <- values[(ends[i] - sizes[i]) + seq_len(sizes[i])]
      m <- length(block)
      l <- diag(exp(v[seq_len(m)]), m)
      l[lower.tri(l)] <- v[-seq_len(m)]
      sigma[block, block] <<- tcrossprod(l) + diag(floor, m)
    }))
    sigma
  }
  covariances <- if (equal) array(build(rest), c(d, d, k)) else
    array(vapply(seq_len(k), \(j) build(rest[(j - 1L) * per + seq_len(per)]),
                 numeric(d * d)), c(d, d, k))
  list(proportions = proportions, means = means, covariances = covariances)
}

#' Maximize: EM from each start, then BFGS on the floor-respecting chart
ref_fit <- function(x, vars, pairs, k, equal = FALSE, floor = 1e-6, starts = NULL,
                    n_random = 6L, seed = 1L, weights = NULL) {
  pattern <- ref_pattern(vars, pairs)
  blocks <- ref_blocks(vars, pairs)
  d <- length(vars)
  complete <- x[stats::complete.cases(x), , drop = FALSE]
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  partitions <- c(starts, lapply(seq_len(n_random), function(i) {
    if (i <= n_random / 2) {
      # k-means on the standardized complete rows, assigned to every row.
      scaled <- scale(complete)
      centres <- stats::kmeans(scaled, k, nstart = 1L)$centers
      filled <- scale(x, attr(scaled, "scaled:center"), attr(scaled, "scaled:scale"))
      filled[is.na(filled)] <- 0
      distance <- vapply(seq_len(k), \(j) rowSums(sweep(filled, 2L, centres[j, ])^2),
                         numeric(nrow(x)))
      diag(k)[max.col(-matrix(distance, nrow(x)), ties.method = "first"), , drop = FALSE]
    } else {
      tau <- matrix(stats::rexp(nrow(x) * k), nrow(x), k)
      tau / rowSums(tau)
    }
  }))
  fits <- lapply(partitions, function(tau) {
    start <- ref_start_from_posteriors(tau, x, pattern, blocks, equal, floor)
    ref_em(start, x, pattern, blocks, equal, floor, weights = weights)
  })
  scores <- vapply(fits, \(f) f$log_likelihood, numeric(1))
  best <- fits[[which.max(scores)]]
  objective <- function(theta) {
    value <- ref_loglik(ref_unpack_chart(theta, blocks, d, k, equal, floor), x,
                        weights = weights)
    if (is.finite(value)) -value else 1e300
  }
  theta <- ref_pack_chart(best$par, blocks, equal, floor)
  polished <- stats::optim(theta, objective, method = "BFGS",
                           control = list(maxit = 2000L, reltol = 1e-15))
  par <- ref_unpack_chart(polished$par, blocks, d, k, equal, floor)
  log_likelihood <- -polished$value
  # BFGS stops where it cannot improve; keep EM's answer if it is better.
  if (log_likelihood < best$log_likelihood) {
    par <- best$par
    log_likelihood <- best$log_likelihood
  }
  list(par = par, log_likelihood = log_likelihood, em = best$log_likelihood,
       start_scores = scores, posteriors = ref_posteriors(par, x))
}

#' Hessian by plain central differences, step h relative to each coordinate
ref_central_hessian <- function(f, theta, h = 1e-4) {
  p <- length(theta)
  steps <- h * pmax(1, abs(theta))
  hessian <- matrix(0, p, p)
  at <- which(upper.tri(hessian, diag = TRUE), arr.ind = TRUE)
  hessian[at] <- apply(at, 1L, function(ij) {
    i <- replace(numeric(p), ij[1L], steps[ij[1L]])
    j <- replace(numeric(p), ij[2L], steps[ij[2L]])
    (f(theta + i + j) - f(theta + i - j) - f(theta - i + j) + f(theta - i - j)) /
      (4 * steps[ij[1L]] * steps[ij[2L]])
  })
  hessian[lower.tri(hessian)] <- t(hessian)[lower.tri(hessian)]
  hessian
}

#' Standard errors from the Hessian of the log likelihood in raw coordinates.
#'
#' The Hessian is a plain central difference (step 1e-4 relative).
#' numDeriv::hessian()'s Richardson extrapolation is computed beside it:
#' on a class whose posteriors are all 0 or 1 it can be off by 10% for that
#' class's covariance block, where central differences at steps 1e-3 and
#' 1e-4 agree to 1e-4 with each other and with the closed form
#' sqrt(2 sigma^2 / n_k) (README).
ref_standard_errors <- function(par, x, vars, pairs, equal, weights = NULL,
                                richardson = FALSE) {
  k <- length(par$proportions)
  theta <- ref_pack_raw(par, vars, pairs, equal)
  f <- function(t) ref_loglik(ref_unpack_raw(t, vars, pairs, k, equal), x, weights = weights)
  hessian <- ref_central_hessian(f, theta)
  gradient <- numDeriv::grad(f, theta)
  covariance <- solve(-hessian)
  names <- ref_raw_names(vars, pairs, k, equal)
  se <- sqrt(diag(covariance))
  richardson_se <- if (richardson) {
    stats::setNames(sqrt(diag(solve(-numDeriv::hessian(f, theta)))), names)
  }
  # Proportion of class 1 is 1 minus the others: its variance is 1'V1.
  first <- if (k > 1L) sqrt(sum(covariance[seq_len(k - 1L), seq_len(k - 1L)])) else 0
  list(theta = stats::setNames(theta, names),
       standard_error = c(stats::setNames(first, "probability:1"), stats::setNames(se, names)),
       gradient = stats::setNames(gradient, names), hessian = hessian,
       richardson_se = richardson_se,
       decrement = as.numeric(crossprod(gradient, covariance %*% gradient)))
}

#' Eigen-decomposition KKT check of a block at an eigenvalue floor.
#'
#' Sigma_block = R(phi) diag(lambda) R(phi)' for a 2 x 2 block (lambda for a
#' single variance). At a constrained maximum: d loglik / d phi = 0,
#' d loglik / d lambda = 0 for an eigenvalue above the floor, and
#' d loglik / d lambda <= 0 for one at it (the multiplier is its negative).
ref_block_kkt <- function(par, x, block, class, floor, equal = FALSE) {
  k <- length(par$proportions)
  classes <- if (equal) seq_len(k) else class
  sigma <- matrix(par$covariances[block, block, class], length(block))
  set_block <- function(value) {
    out <- par
    invisible(lapply(classes, function(j) out$covariances[block, block, j] <<- value))
    out
  }
  if (length(block) == 1L) {
    f <- function(t) ref_loglik(set_block(matrix(t[1L])), x)
    g <- numDeriv::grad(f, sigma[1L, 1L])
    return(data.frame(lambda_1 = sigma[1L, 1L], lambda_2 = NA_real_,
                      gradient_1 = g, gradient_2 = NA_real_, gradient_angle = NA_real_))
  }
  e <- eigen(sigma, symmetric = TRUE)
  phi <- atan2(e$vectors[2L, 1L], e$vectors[1L, 1L])
  f <- function(t) {
    rotation <- matrix(c(cos(t[3L]), sin(t[3L]), -sin(t[3L]), cos(t[3L])), 2L)
    ref_loglik(set_block(rotation %*% diag(t[1:2]) %*% t(rotation)), x)
  }
  g <- numDeriv::grad(f, c(e$values, phi))
  data.frame(lambda_1 = e$values[1L], lambda_2 = e$values[2L],
             gradient_1 = g[1L], gradient_2 = g[2L], gradient_angle = g[3L])
}

#' The best permutation of `mine`'s columns onto `theirs`
ref_relabel <- function(mine, theirs) {
  permutations <- function(k) if (k == 1L) list(1L) else
    do.call(c, lapply(seq_len(k), \(i) lapply(permutations(k - 1L),
                                             \(p) c(i, setdiff(seq_len(k), i)[p]))))
  candidates <- permutations(ncol(mine))
  gaps <- vapply(candidates, \(p) max(abs(mine[, p, drop = FALSE] - theirs)), numeric(1))
  list(order = candidates[[which.min(gaps)]], gap = min(gaps))
}

# -- the study datasets ---------------------------------------------------------

#' The reference study: many cheap simulab datasets across sizes, numbers of
#' classes and indicators, one or two pairs, varying and equal covariances,
#' and a few with missing values. Deterministic in `index`.
ref_study_case <- function(index) {
  seed <- 20261100L + index
  k <- c(2L, 2L, 3L)[(index %% 3L) + 1L]
  d <- c(2L, 3L, 4L, 5L)[(index %% 4L) + 1L]
  n <- c(40L, 80L, 200L, 500L, 2000L)[(index %% 5L) + 1L]
  # Small samples need separation to be estimable at all.
  n <- if (k == 3L && n < 80L) 80L else n
  equal <- index %% 2L == 0L
  pairs <- if (d >= 4L && index %% 7L < 4L) c("y1:y2", "y3:y4") else
    if (d == 2L) "y1:y2" else c("y1:y3", "y2:y3")[(index %% 2L) + 1L]
  missing <- index %% 10L == 3L && d >= 3L
  set.seed(seed)
  separation <- if (n < 100L) 3 else 2
  means <- t(vapply(seq_len(k), \(j) separation * (j - 1) *
                      sign(sin(seq_len(d) * j + index)), numeric(d)))
  sds <- if (equal) matrix(stats::runif(d, 0.6, 1.4), k, d, byrow = TRUE) else
    matrix(stats::runif(k * d, 0.6, 1.4), k, d)
  vars <- paste0("y", seq_len(d))
  values <- if (equal) matrix(stats::runif(length(pairs), -0.8, 0.8), k, length(pairs),
                              byrow = TRUE) else
    matrix(stats::runif(k * length(pairs), -0.8, 0.8), k, length(pairs))
  correlations <- lapply(seq_len(k), function(j) {
    m <- diag(d)
    invisible(Map(function(pair, v) {
      at <- match(strsplit(pair, ":", fixed = TRUE)[[1L]], vars)
      m[at[1L], at[2L]] <<- v
      m[at[2L], at[1L]] <<- v
    }, pairs, values[j, ]))
    m
  })
  proportions <- if (k == 2L) c(0.55, 0.45) else c(0.4, 0.35, 0.25)
  simulated <- simulab::simulate_lpa(n = n, means = means, sds = sds,
                                     proportions = proportions,
                                     correlations = correlations, seed = seed)
  frame <- as.data.frame(simulated)
  data <- stats::setNames(frame[paste0("indicator_", seq_len(d))], vars)
  if (missing) {
    data[vars[-1L]] <- simulab::inject_missingness(data[vars[-1L]], mechanism = "MCAR",
                                                   proportion = 0.1, seed = seed + 1L)
  }
  list(index = index, seed = seed, data = data, vars = vars, pairs = pairs, k = k,
       equal = equal, missing = missing, n = n,
       truth = list(means = means, sds = sds, correlations = correlations,
                    proportions = proportions),
       membership = as.integer(sub("Profile ", "", frame$profile)))
}

# -- the floor solver of the constrained diagonal family (commit a40d5c7) -------------
#
# Not on the pair path (a pair fit is a patterned full covariance whose
# floor is per-block eigenvalue clipping), but the same commit's floor
# update: the Gaussian M-step over variances sigma_kj = lambda_k a_kj with
# volume lambda and determinant-one shape a, profiles k and indicators j,
#   minimize sum_kj (w_k log sigma_kj + s_kj / sigma_kj)  s.t. sigma_kj >= f,
# VEI (shape shared, volumes free) or EVI (volume shared, shapes free).
# The reference builds its own chart, log lambda and log a with the last
# log shape minus the sum of the others, and solves the convex problem with
# stats::constrOptim, whose barrier keeps it strictly feasible.

#' The chart's design: rows (k, j) profile-major, columns the coordinates
ref_floor_design <- function(k, d, structure) {
  volume <- if (structure == "EVI") matrix(1, k * d, 1L) else
    kronecker(diag(k), matrix(1, d, 1L))
  shape_one <- rbind(diag(d - 1L), rep(-1, d - 1L))
  shape <- if (structure == "VEI") kronecker(matrix(1, k, 1L), shape_one) else
    kronecker(diag(k), shape_one)
  cbind(volume, shape)
}

#' The objective (per unit weight) at variances sigma (k x d)
ref_floor_objective <- function(sigma, scatter, weights) {
  sum(weights * log(sigma) + scatter / sigma) / sum(weights)
}

#' Independent constrained optimum: constrOptim from a feasible interior point
ref_floor_solve <- function(scatter, weights, structure, floor) {
  k <- nrow(scatter)
  d <- ncol(scatter)
  design <- ref_floor_design(k, d, structure)
  counts <- rep(weights, each = d)
  values <- as.vector(t(scatter))
  objective <- function(theta) {
    logs <- as.vector(design %*% theta)
    sum(counts * logs + values * exp(-logs)) / sum(weights)
  }
  gradient <- function(theta) {
    logs <- as.vector(design %*% theta)
    as.vector(crossprod(design, counts - values * exp(-logs))) / sum(weights)
  }
  n_volume <- if (structure == "EVI") 1L else k
  start <- c(rep(log(max(scatter / weights, floor)) + 1, n_volume),
             numeric(ncol(design) - n_volume))
  # A decreasing barrier weight, each solve starting where the last ended. A
  # tiny weight lets the barrier round a slack to zero, after which
  # constrOptim cannot evaluate its start; the last finished solve is kept
  # (its objective is within about (rows x weight) of the optimum).
  solve_at <- function(theta, mu) {
    stats::constrOptim(theta, objective, gradient, ui = design,
                       ci = rep(log(floor), k * d), mu = mu,
                       control = list(reltol = 1e-14, maxit = 5000L),
                       outer.iterations = 200L, outer.eps = 1e-14)$par
  }
  reached <- Reduce(function(state, mu) {
    if (!is.null(state$stopped)) return(state)
    tryCatch(list(par = solve_at(state$par, mu), mu = mu),
             error = function(e) {
               rounded <- grepl("not finite|not in the interior", conditionMessage(e))
               if (!rounded) stop(e)
               c(state, list(stopped = mu))
             })
  }, c(1e-4, 1e-6, 1e-8, 1e-10, 1e-12), list(par = start, mu = NA_real_))
  logs <- as.vector(design %*% reached$par)
  list(sigma = matrix(exp(logs), k, d, byrow = TRUE), value = objective(reached$par),
       mu = reached$mu)
}

#' KKT residuals of a candidate: stationarity on the chart with nonnegative
#' multipliers on the rows at the floor
ref_floor_kkt <- function(sigma, scatter, weights, structure, floor) {
  k <- nrow(scatter)
  d <- ncol(scatter)
  design <- ref_floor_design(k, d, structure)
  logs <- log(as.vector(t(sigma)))
  counts <- rep(weights, each = d)
  values <- as.vector(t(scatter))
  gradient <- as.vector(crossprod(design, counts - values * exp(-logs))) / sum(weights)
  active <- which(abs(logs - log(floor)) <= 1e-9)
  # gradient = design[active, ]' mu with mu >= 0. Active rows can be linearly
  # dependent (a profile with every variance at the floor), so the
  # multipliers are not unique: nonnegative least squares (nnls) finds one
  # nonnegative set when any exists, and the residual is the stationarity.
  multipliers <- if (length(active)) {
    nnls::nnls(t(design[active, , drop = FALSE]), gradient)$x
  } else numeric()
  residual <- gradient - if (length(active))
    as.vector(t(design[active, , drop = FALSE]) %*% multipliers) else 0
  list(feasible = min(logs - log(floor)), stationarity = max(abs(residual)),
       min_multiplier = if (length(active)) min(multipliers) else 0,
       n_active = length(active),
       # The family's own constraints: shared shape (VEI) or volume (EVI) and
       # determinant-one shapes.
       family = if (structure == "VEI") {
         shapes <- sigma / exp(rowMeans(log(sigma)))
         max(abs(sweep(shapes, 2L, shapes[1L, ])))
       } else {
         volumes <- exp(rowMeans(log(sigma)))
         max(abs(volumes - volumes[1L])) / volumes[1L]
       })
}

#' Random floor problems that bind: scatter with tiny entries and a floor
#' between the smallest and the typical unconstrained variance
ref_floor_problem <- function(index) {
  set.seed(20261200L + index)
  k <- 2L + index %% 3L
  d <- 2L + index %% 4L
  weights <- stats::runif(k, 5, 100)
  scatter <- matrix(stats::rexp(k * d), k, d) * weights
  tiny <- sample.int(k * d, 1L + index %% 3L)
  scatter[tiny] <- scatter[tiny] * 10^-stats::runif(length(tiny), 2, 6)
  if (index %% 11L == 0L) scatter[1L, ] <- 0
  ratio <- scatter / weights
  floor <- stats::quantile(ratio[ratio > 0], 0.25) * stats::runif(1L, 0.5, 2)
  list(scatter = scatter, weights = weights, floor = unname(floor),
       structure = c("VEI", "EVI")[index %% 2L + 1L])
}
