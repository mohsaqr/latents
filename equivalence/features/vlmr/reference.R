# An independent reference for the Vuong-Lo-Mendell-Rubin test. Nothing here
# calls latents: each model's per-unit log likelihood is written out again in
# its own coordinates (multinomial logits against the LAST class, log
# variances, a log-Cholesky factor, item logits), per-unit scores come from
# numDeriv::jacobian and the Hessian from numDeriv::hessian, Vuong's (1989)
# matrix is assembled from the published formula and its eigenvalues are the
# weights. The tail probability of the weighted chi-square is then computed
# by CompQuadForm's Imhof, Davies and Farebrother algorithms and by Monte
# Carlo. Eigenvalues of W do not depend on the coordinates, so agreeing with
# latents here checks its scores, Hessian, whitening and Imhof together.
#
# A model is described by a list: `model` ("EEI", "VVI", "EEE" or "LCA"),
# `k`, and its estimates `proportions` (k), `means` (k x p), `variances`
# (k x p, rows equal for EEI), `covariance` (p x p, EEE) or `probabilities`
# (k x p, the probability of the item's second category, LCA).

#' Pack estimates into the reference coordinates
ref_pack <- function(estimates) {
  k <- estimates$k
  logits <- if (k > 1L) log(estimates$proportions[-k] / estimates$proportions[k]) else numeric()
  switch(estimates$model,
    EEI = c(logits, as.vector(estimates$means), log(estimates$variances[1L, ])),
    VVI = c(logits, as.vector(estimates$means), as.vector(log(estimates$variances))),
    EEE = {
      root <- t(chol(estimates$covariance))
      diag(root) <- log(diag(root))
      c(logits, as.vector(estimates$means), root[lower.tri(root, diag = TRUE)])
    },
    LCA = c(logits, as.vector(stats::qlogis(estimates$probabilities))))
}

#' The natural size of each reference coordinate: a mean's class standard
#' deviation, one for logits and (log-)variance coordinates
ref_natural_scale <- function(estimates) {
  k <- estimates$k
  if (identical(estimates$model, "LCA")) return(rep(1, length(ref_pack(estimates))))
  sds <- if (identical(estimates$model, "EEE")) {
    matrix(sqrt(diag(estimates$covariance)), k, ncol(estimates$means), byrow = TRUE)
  } else sqrt(estimates$variances)
  rest <- length(ref_pack(estimates)) - (k - 1L) - length(sds)
  if (!identical(estimates$model, "EEE")) return(c(rep(1, k - 1L), as.vector(sds), rep(1, rest)))
  # Log-Cholesky: log diagonals by one, off-diagonals by their row's scale.
  p <- ncol(estimates$means)
  rows <- matrix(sqrt(diag(estimates$covariance)), p, p)
  rows[row(rows) == col(rows)] <- 1
  c(rep(1, k - 1L), as.vector(sds), rows[lower.tri(rows, diag = TRUE)])
}

#' Unpack reference coordinates into proportions and component densities
ref_unpack <- function(theta, model, k, p) {
  logits <- c(theta[seq_len(k - 1L)], 0)
  proportions <- exp(logits - max(logits))
  proportions <- proportions / sum(proportions)
  rest <- theta[seq_along(theta) > k - 1L]
  if (identical(model, "LCA")) {
    return(list(proportions = proportions, probabilities = matrix(stats::plogis(rest), k, p)))
  }
  means <- matrix(rest[seq_len(k * p)], k, p)
  rest <- rest[seq_along(rest) > k * p]
  switch(model,
    EEI = list(proportions = proportions, means = means,
               variances = matrix(exp(rest), k, p, byrow = TRUE)),
    VVI = list(proportions = proportions, means = means, variances = matrix(exp(rest), k, p)),
    EEE = {
      root <- matrix(0, p, p)
      root[lower.tri(root, diag = TRUE)] <- rest
      diag(root) <- exp(diag(root))
      list(proportions = proportions, means = means, covariance = tcrossprod(root))
    })
}

#' Per-unit log likelihoods (unweighted); missing values are skipped
#' (diagonal covariances only), rows with nothing observed contribute zero
ref_unit_loglik <- function(theta, x, model, k) {
  p <- ncol(x)
  parts <- ref_unpack(theta, model, k, p)
  observed <- !is.na(x)
  filled <- x
  filled[!observed] <- 0
  component <- vapply(seq_len(k), function(j) {
    if (identical(model, "LCA")) {
      probability <- matrix(parts$probabilities[j, ], nrow(x), p, byrow = TRUE)
      cell <- ifelse(filled == 2, log(probability), log1p(-probability))
      return(rowSums(cell * observed))
    }
    if (identical(model, "EEE")) {
      stopifnot("EEE reference needs complete rows" = all(observed))
      root <- chol(parts$covariance)
      centred <- sweep(x, 2L, parts$means[j, ])
      z <- forwardsolve(t(root), t(centred))
      return(-0.5 * colSums(z^2) - sum(log(diag(root))) - 0.5 * p * log(2 * pi))
    }
    mean <- matrix(parts$means[j, ], nrow(x), p, byrow = TRUE)
    variance <- matrix(parts$variances[j, ], nrow(x), p, byrow = TRUE)
    cell <- -0.5 * (log(2 * pi * variance) + (filled - mean)^2 / variance)
    rowSums(cell * observed)
  }, numeric(nrow(x)))
  component <- matrix(component, nrow(x), k)
  logged <- sweep(component, 2L, log(parts$proportions), "+")
  top <- apply(logged, 1L, max)
  value <- top + log(rowSums(exp(logged - top)))
  value[rowSums(observed) == 0L] <- 0
  value
}

#' Scores, Hessian and log likelihood of one fitted model
#'
#' With sampling weights the pseudo log likelihood is sum_i w_i l_i, the
#' per-unit scores are w_i s_i and the Hessian is that of the weighted sum.
#' `convention = "frequency"` treats the weights as replication counts, as
#' Latent GOLD's `caseweight` does: the per-unit scores are sqrt(w_i) s_i, so
#' B = sum_i w_i s_i s_i' instead of sum_i w_i^2 s_i s_i'.
ref_model_parts <- function(estimates, x, weights = NULL,
                            convention = c("sampling", "frequency")) {
  convention <- match.arg(convention)
  theta <- ref_pack(estimates)
  w <- weights %||% rep(1, nrow(x))
  # numDeriv steps are relative to each coordinate's size, so coordinates at
  # an estimate near zero (a mean of 0.001) get steps of 1e-7 and roundoff
  # swamps the Hessian of a 4000-row log likelihood, while a mean whose class
  # has a standard deviation of 0.03 is stepped right past its curvature.
  # Differentiating in the affine coordinates u = 1 + (theta - estimate) /
  # scale, with each mean scaled by its class's standard deviation and every
  # logit or log-variance coordinate by one, gives each coordinate a step of
  # its natural size. Vuong's eigenvalues do not depend on the coordinates,
  # so the change is exact.
  scale <- ref_natural_scale(estimates)
  origin <- rep(1, length(theta))
  unit <- function(u) ref_unit_loglik(theta + (u - 1) * scale, x, estimates$model, estimates$k)
  scores <- numDeriv::jacobian(unit, origin) *
    if (identical(convention, "frequency")) sqrt(w) else w
  hessian <- numDeriv::hessian(function(u) sum(w * unit(u)), origin)
  list(theta = theta, scores = scores, hessian = hessian, log_likelihood = sum(w * unit(origin)),
       units = unit(origin),
       # The summed score in the original coordinates.
       score_sum = colSums(scores) / scale)
}

#' Eigenvalues of Vuong's W (Vuong 1989, eq. 3.6; nonnest2::calcLambda's layout)
#'
#' W = [ -B1 A1^-1, -B10 A0^-1 ; B10' A1^-1, B0 A0^-1 ], A the Hessian and B the
#' cross-products of per-unit scores. Returns the real parts, decreasing.
ref_vuong_eigenvalues <- function(alternative, null) {
  b1 <- crossprod(alternative$scores)
  b0 <- crossprod(null$scores)
  b10 <- crossprod(alternative$scores, null$scores)
  a1_inverse <- solve(alternative$hessian)
  a0_inverse <- solve(null$hessian)
  w <- rbind(cbind(-b1 %*% a1_inverse, -b10 %*% a0_inverse),
             cbind(t(b10) %*% a1_inverse, b0 %*% a0_inverse))
  sort(Re(eigen(w, only.values = TRUE)$values), decreasing = TRUE)
}

#' Upper tail of sum lambda_j chi-square(1) by three independent algorithms
ref_tail <- function(q, lambda) {
  # CompQuadForm::imhof warns "Qq + abserr is positive" when a tail near zero
  # comes out slightly negative; that diagnosed note is muffled, nothing else.
  imhof <- withCallingHandlers(
    CompQuadForm::imhof(q, lambda, epsabs = 1e-12, epsrel = 1e-12, limit = 1e5)$Qq,
    warning = function(w) {
      if (grepl("Qq + abserr is positive", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  c(imhof = imhof,
    davies = CompQuadForm::davies(q, lambda, acc = 1e-10, lim = 1e6)$Qq,
    farebrother = if (all(lambda >= 0)) {
      CompQuadForm::farebrother(q, lambda[lambda > 0], eps = 1e-14)$Qq
    } else NA_real_)
}

#' Monte Carlo upper tails at every q, from one set of draws (chunks of 1e5)
ref_monte_carlo <- function(q, lambda, draws = 1e6, seed = 1L) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  chunk <- 1e5
  n_chunks <- ceiling(draws / chunk)
  exceed <- vapply(seq_len(n_chunks), function(index) {
    sums <- crossprod(lambda, matrix(stats::rchisq(chunk * length(lambda), 1), length(lambda)))
    vapply(q, function(value) sum(sums > value), numeric(1))
  }, numeric(length(q)))
  rowSums(matrix(exceed, length(q))) / (chunk * n_chunks)
}

#' The whole reference comparison of two fitted models on the same units
ref_vlmr <- function(alternative, null, x, weights = NULL, n_units = NULL, draws = 1e6,
                     convention = "sampling") {
  a <- ref_model_parts(alternative, x, weights, convention)
  b <- ref_model_parts(null, x, weights, convention)
  lambda <- ref_vuong_eigenvalues(a, b)
  lr <- 2 * (a$log_likelihood - b$log_likelihood)
  added <- length(a$theta) - length(b$theta)
  n_units <- n_units %||% sum(rowSums(!is.na(x)) > 0L)
  adjusted <- lr / (1 + 1 / (added * log(n_units)))
  monte_carlo <- ref_monte_carlo(c(max(lr, 0), max(adjusted, 0)), lambda, draws)
  list(lambda = lambda, lr = lr, adjusted_lr = adjusted, parameters_added = added,
       n_units = n_units, reference_mean = sum(lambda), reference_sd = sqrt(2 * sum(lambda^2)),
       vlmr = c(ref_tail(max(lr, 0), lambda), monte_carlo = monte_carlo[1L]),
       almr = c(ref_tail(max(adjusted, 0), lambda), monte_carlo = monte_carlo[2L]),
       draws = draws, score_sums = c(alternative = max(abs(a$score_sum)),
                                     null = max(abs(b$score_sum))),
       parts = list(alternative = a, null = b))
}
