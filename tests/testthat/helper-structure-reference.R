# An independent reference for Wald inference under a constrained covariance
# structure. It shares nothing with the package's inference code: its own
# mixture likelihood (written from the definition, with its own Gaussian
# density), its own chart of each structure (Helmert log shapes, Givens-angle
# orientations, a trace-zero matrix logarithm for determinant-one matrices,
# softmax logits against the *first* category) and numerical rather than
# analytic derivatives. Standard errors of natural parameters do not depend on
# the chart at a maximum, so the two must agree.

# The mixture log likelihood of a two-level (or single-level) Gaussian mixture.
reference_mixture_loglik <- function(x, group, means, sigma, profile_probabilities,
                                     group_probabilities) {
  k <- nrow(means)
  log_density <- vapply(seq_len(k), function(profile) {
    factor <- chol(sigma[[profile]])
    residuals <- sweep(x, 2L, means[profile, ], "-")
    solved <- backsolve(factor, t(residuals), transpose = TRUE)
    -0.5 * colSums(solved^2) - sum(log(diag(factor))) - ncol(x) / 2 * log(2 * base::pi)
  }, numeric(nrow(x)))
  per_class <- vapply(seq_along(group_probabilities), function(class) {
    joint <- sweep(log_density, 2L, log(profile_probabilities[class, ]), "+")
    top <- apply(joint, 1L, max)
    rowsum(top + log(rowSums(exp(joint - top))), group, reorder = FALSE)[, 1L]
  }, numeric(length(unique(group))))
  per_class <- matrix(per_class, ncol = length(group_probabilities))
  per_class <- sweep(per_class, 2L, log(group_probabilities), "+")
  top <- apply(per_class, 1L, max)
  sum(top + log(rowSums(exp(per_class - top))))
}

reference_helmert <- function(d) {
  if (d == 1L) return(matrix(0, 1L, 0L))
  basis <- stats::contr.helmert(d)
  sweep(basis, 2L, sqrt(colSums(basis^2)), "/")
}

reference_givens <- function(angles, d) {
  pairs <- which(upper.tri(matrix(0, d, d)), arr.ind = TRUE)
  Reduce(`%*%`, lapply(seq_len(nrow(pairs)), function(pair) {
    rotation <- diag(d)
    a <- pairs[pair, 1L]
    b <- pairs[pair, 2L]
    rotation[a, a] <- rotation[b, b] <- cos(angles[pair])
    rotation[a, b] <- -sin(angles[pair])
    rotation[b, a] <- sin(angles[pair])
    rotation
  }), diag(d))
}

reference_sym_exp <- function(m) {
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  e$vectors %*% (exp(e$values) * t(e$vectors))
}

reference_sym_log <- function(m) {
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  e$vectors %*% (log(e$values) * t(e$vectors))
}

# Build the reference chart of a fit: a starting vector and a map from it to
# (means, covariance list, profile probabilities, group probabilities).
reference_structure_chart <- function(fit) {
  code <- fit$covariance_structure
  k <- fit$n_profiles
  h <- fit$n_group_classes
  d <- ncol(fit$means)
  sigma0 <- lapply(seq_len(k), function(profile) {
    if (is.null(fit$covariances)) diag(unname(fit$variances[profile, ]), d) else
      unname(matrix(fit$covariances[, , profile], d, d))
  })
  log_volume <- vapply(sigma0, function(s) log(det(s)) / d, numeric(1))
  volume_shared <- startsWith(code, "E")
  second <- substr(code, 2L, 2L)
  third <- substr(code, 3L, 3L)
  form <- if (code %in% c("EII", "VII")) "spherical" else if (third == "I") "diagonal" else
    if (code %in% c("EEE", "VEE", "EVV", "VVV")) "logm" else "rotation"
  shape_rows <- if (form == "spherical") 0L else if (second == "E") 1L else k
  shape_of <- if (second == "E") rep(1L, k) else seq_len(k)
  rotation_of <- if (third == "E") rep(1L, k) else seq_len(k)
  rotation_rows <- if (form == "rotation") length(unique(rotation_of)) else 0L
  helmert <- reference_helmert(d)
  anchors <- NULL
  shape0 <- list()
  if (form == "diagonal") {
    t0 <- lapply(seq_len(k), function(p) log(diag(sigma0[[p]])) - log_volume[p])
    shape0 <- if (shape_rows == 1L) list(Reduce(`+`, t0) / k) else t0
    shape0 <- lapply(shape0, function(t) as.vector(crossprod(helmert, t)))
  } else if (form == "rotation") {
    anchors <- if (rotation_rows == 1L) {
      list(eigen(Reduce(`+`, Map(`*`, sigma0, seq_len(k) + 0.37)),
                 symmetric = TRUE)$vectors)
    } else lapply(sigma0, function(s) eigen(s, symmetric = TRUE)$vectors)
    t0 <- lapply(seq_len(k), function(p) {
      a <- anchors[[rotation_of[p]]]
      log(diag(crossprod(a, sigma0[[p]] %*% a))) - log_volume[p]
    })
    shape0 <- if (shape_rows == 1L) list(Reduce(`+`, t0) / k) else t0
    shape0 <- lapply(shape0, function(t) as.vector(crossprod(helmert, t)))
  } else if (form == "logm") {
    unit <- lapply(seq_len(k), function(p) sigma0[[p]] / exp(log_volume[p]))
    base <- if (shape_rows == 1L) list(Reduce(`+`, unit) / k) else unit
    anchors <- lapply(base, reference_sym_log)
    shape0 <- lapply(base, function(b) numeric(d * (d + 1L) / 2L - 1L))
  }
  volume0 <- if (volume_shared) mean(log_volume) else log_volume
  n_rotation <- d * (d - 1L) / 2L
  pieces <- list(means = as.vector(t(unname(fit$means))),
                 volume = volume0,
                 shape = unlist(shape0),
                 rotation = rep(0, rotation_rows * n_rotation),
                 profile = as.vector(t(log(fit$profile_probabilities[, -1L, drop = FALSE]) -
                                         log(fit$profile_probabilities[, 1L]))),
                 group = log(fit$group_probabilities[-1L]) - log(fit$group_probabilities[1L]))
  widths <- lengths(pieces)
  start <- unlist(pieces, use.names = FALSE)
  per_shape <- if (shape_rows == 0L) 0L else length(shape0[[1L]])
  upper <- which(upper.tri(matrix(0, d, d), diag = TRUE))
  decode <- function(psi) {
    at <- cumsum(c(0L, widths))
    block <- function(i) psi[at[i] + seq_len(widths[i])]
    means <- matrix(block(1L), k, d, byrow = TRUE)
    volume <- block(2L)
    shape <- block(3L)
    rotation <- block(4L)
    sigma <- lapply(seq_len(k), function(p) {
      lambda <- exp(volume[if (volume_shared) 1L else p])
      s <- shape[(shape_of[p] - 1L) * per_shape + seq_len(per_shape)]
      switch(form,
        spherical = lambda * diag(d),
        diagonal = lambda * diag(exp(as.vector(helmert %*% s)), d),
        rotation = {
          r <- anchors[[rotation_of[p]]] %*%
            reference_givens(rotation[(rotation_of[p] - 1L) * n_rotation +
                                        seq_len(n_rotation)], d)
          lambda * r %*% (exp(as.vector(helmert %*% s)) * t(r))
        },
        logm = {
          delta <- matrix(0, d, d)
          values <- c(s, 0)
          delta[upper] <- values
          delta <- delta + t(delta) - diag(diag(delta), d)
          # Trace zero: the last diagonal entry is minus the others.
          delta[d, d] <- -sum(diag(delta)[-d])
          lambda * reference_sym_exp(anchors[[shape_of[p]]] + delta)
        })
    })
    softmax_first <- function(logits) {
      v <- c(0, logits)
      exp(v - max(v)) / sum(exp(v - max(v)))
    }
    profile <- matrix(block(5L), h, k - 1L, byrow = TRUE)
    list(means = means, sigma = sigma,
         profile_probabilities = t(apply(profile, 1L, softmax_first)),
         group_probabilities = softmax_first(block(6L)))
  }
  list(start = start, decode = decode, form = form)
}

# Map reference parameters to the fit's natural coefficients, in coef() order.
reference_natural <- function(fit, parameters) {
  d <- ncol(fit$means)
  k <- fit$n_profiles
  ellipsoidal <- !endsWith(fit$covariance_structure, "I")
  shared <- fit$covariance_structure %in% c("EII", "EEI", "EEE")
  rows <- if (shared) 1L else seq_len(k)
  lower <- which(lower.tri(matrix(0, d, d), diag = TRUE))
  spread <- unlist(lapply(rows, function(p) {
    s <- parameters$sigma[[p]]
    if (ellipsoidal) s[lower] else diag(s)
  }))
  c(as.vector(t(parameters$means)), spread,
    as.vector(t(parameters$profile_probabilities)),
    parameters$group_probabilities)
}

# Central-difference Hessian and Jacobian, written out rather than borrowed.
reference_hessian <- function(f, x, step = 2e-4) {
  p <- length(x)
  unit <- diag(step, p)
  entries <- expand.grid(i = seq_len(p), j = seq_len(p))
  entries <- entries[entries$i <= entries$j, ]
  values <- vapply(seq_len(nrow(entries)), function(r) {
    i <- entries$i[r]
    j <- entries$j[r]
    (f(x + unit[, i] + unit[, j]) - f(x + unit[, i] - unit[, j]) -
       f(x - unit[, i] + unit[, j]) + f(x - unit[, i] - unit[, j])) / (4 * step^2)
  }, numeric(1))
  out <- matrix(0, p, p)
  out[cbind(entries$i, entries$j)] <- values
  out[cbind(entries$j, entries$i)] <- values
  out
}

reference_jacobian <- function(f, x, step = 1e-6) {
  p <- length(x)
  unit <- diag(step, p)
  vapply(seq_len(p), function(i) (f(x + unit[, i]) - f(x - unit[, i])) / (2 * step),
         numeric(length(f(x))))
}

# Natural-scale standard errors of a structure fit from the reference chart.
# `hessian` and `jacobian` default to the hand-written differences; the
# equivalence suite passes numDeriv's.
reference_structure_se <- function(fit, data, hessian = reference_hessian,
                                   jacobian = reference_jacobian) {
  chart <- reference_structure_chart(fit)
  x <- as.matrix(data[fit$continuous])
  group <- if (isTRUE(fit$single_level)) seq_len(nrow(data)) else
    match(data[[fit$id]], unique(data[[fit$id]]))
  loglik <- function(psi) {
    p <- chart$decode(psi)
    reference_mixture_loglik(x, group, p$means, p$sigma, p$profile_probabilities,
                             p$group_probabilities)
  }
  information <- -hessian(loglik, chart$start)
  natural <- function(psi) reference_natural(fit, chart$decode(psi))
  map <- jacobian(natural, chart$start)
  covariance <- map %*% solve(information, t(map))
  list(standard_error = sqrt(pmax(diag(covariance), 0)),
       estimate = natural(chart$start),
       log_likelihood = loglik(chart$start))
}
