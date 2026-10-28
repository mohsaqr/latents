pkgload::load_all(quiet = TRUE)
folder <- "validation/release-0.96.1-2026-10-20/covariance-followup/cross-review-v2"
cost <- function(variances, scatter, weights) {
  sum(weights * rowSums(log(variances)) + rowSums(scatter / variances))
}
# Independent scalar reductions in d=2: condition on equal log volume for
# EVI or equal determinant-one log shape for VEI, then optimize one scalar.
reference <- function(scatter, weights, structure, floor) {
  given <- if (structure == "EVI") function(t) {
    volume <- exp(t)
    ratio <- sqrt(scatter[, 1] / scatter[, 2])
    ratio[!is.finite(ratio)] <- ifelse(scatter[!is.finite(ratio), 1] > 0, Inf, 1)
    first <- pmax(floor, pmin(volume^2 / floor, volume * ratio))
    cbind(first, volume^2 / first)
  } else function(t) {
    shape <- c(exp(t), exp(-t))
    volume <- pmax(rowSums(sweep(scatter, 2L, shape, "/")) / (2 * weights),
                    floor / min(shape))
    outer(volume, shape)
  }
  interval <- if (structure == "EVI") {
    c(log(floor), log(max(floor, max(scatter / weights)) * 10))
  } else c(-30, 30)
  optimum <- stats::optimize(function(t) cost(given(t), scatter, weights),
                              interval, tol = 1e-12)
  given(optimum$minimum)
}
set.seed(1907)
scalar <- lapply(seq_len(100), function(index) {
  k <- 2L + index %% 3L
  scatter <- matrix(exp(stats::runif(k * 2, -15, 6)), k, 2)
  if (index %% 4L == 0) scatter[sample(length(scatter), 1L)] <- 0
  weights <- exp(stats::runif(k, -2, 3))
  floor <- exp(stats::runif(1, -5, 1))
  structure <- if (index %% 2) "EVI" else "VEI"
  actual <- .multilpa_bounded_variances(scatter, weights, structure, floor)
  oracle <- reference(scatter, weights, structure, floor)
  excess <- (cost(actual, scatter, weights) - cost(oracle, scatter, weights)) /
    (1 + abs(cost(oracle, scatter, weights)))
  stopifnot(min(actual) >= floor * (1 - 1e-11), excess < 1e-9)
  data.frame(index = index, structure = structure, normalized_excess = excess)
})
scalar <- do.call(rbind, scalar)
print(summary(scalar$normalized_excess), digits = 15)
utils::write.csv(scalar, file.path(folder, "scalar-references.csv"), row.names = FALSE)
# Independently recompute feasibility, gradient stationarity, multiplier
# signs and complementary slackness from final variances, not solver flags.
debug_solver <- .multilpa_bounded_variances
body(debug_solver)[[length(body(debug_solver))]] <- quote(list(
  result = matrix(exp(logs), k, d, byrow = TRUE), active = active,
  iterations = iteration, settled = settled))
set.seed(729)
kkt <- lapply(seq_len(100), function(index) {
  k <- 2L + index %% 3L
  d <- 2L + index %% 4L
  scatter <- matrix(exp(stats::runif(k * d, -12, 6)), k, d)
  if (index %% 3L == 0) scatter[sample(length(scatter), 2L)] <- 0
  weights <- exp(stats::runif(k, -2, 2))
  floor <- exp(stats::runif(1, -4, 1))
  structure <- if (index %% 2) "EVI" else "VEI"
  fit <- debug_solver(scatter, weights, structure, floor)
  variances <- fit$result
  # Explicit independent volume and zero-sum-shape design.
  volume <- if (structure == "EVI") matrix(1, k * d, 1) else
    diag(k)[rep(seq_len(k), each = d), , drop = FALSE]
  shape <- rbind(diag(d - 1), rep(-1, d - 1))
  if (structure == "VEI") shape <- shape[rep(seq_len(d), k), , drop = FALSE] else {
    shape <- kronecker(diag(k), shape)
  }
  design <- cbind(volume, shape)
  slack <- as.vector(t(log(variances))) - log(floor)
  gradient <- as.vector(crossprod(design,
    as.vector(t(matrix(weights, k, d) - scatter / variances)))) / sum(weights)
  active <- fit$active
  if (length(active)) {
    constraints <- design[active, , drop = FALSE]
    lambda <- as.vector(solve(tcrossprod(constraints), constraints %*% gradient))
    residual <- gradient - as.vector(crossprod(constraints, lambda))
    complementarity <- max(abs(lambda * slack[active]))
    minimum_lambda <- min(lambda)
  } else {
    residual <- gradient
    complementarity <- 0
    minimum_lambda <- 0
  }
  tolerance <- 1e-10 * (1 + max(abs(gradient)))
  stopifnot(min(slack) > -1e-10, max(abs(residual)) < tolerance,
              complementarity < tolerance, minimum_lambda > -tolerance)
  data.frame(index = index, structure = structure, iterations = fit$iterations,
    minimum_slack = min(slack), stationarity = max(abs(residual)),
    complementarity = complementarity, minimum_multiplier = minimum_lambda)
})
kkt <- do.call(rbind, kkt)
print(summary(kkt), digits = 15)
utils::write.csv(kkt, file.path(folder, "kkt-certificates.csv"), row.names = FALSE)
# Original exact oracle over all originally problematic neighboring scales.
scales <- 10^seq(-2, -24, by = -1)
oracle <- matrix(c(1, 1, 50.5, 50.5), 2)
sweep <- lapply(scales, function(small) {
  scatter <- rbind(c(small, 1), c(10 * small, 100))
  actual <- .multilpa_bounded_variances(scatter, c(1, 1), "EVI", 1)
  error <- max(abs(actual - oracle))
  stopifnot(error < 1e-10)
  data.frame(small = small, maximum_error = error,
    objective_excess = cost(actual, scatter, c(1, 1)) - cost(oracle, scatter, c(1, 1)))
})
sweep <- do.call(rbind, sweep)
print(sweep, digits = 16)
utils::write.csv(sweep, file.path(folder, "original-oracle-sweep.csv"), row.names = FALSE)
# Degenerate curvature and redundant floor constraints at larger dimensions.
zeros <- expand.grid(k = 2:5, d = 2:6, structure = c("EVI", "VEI"), stringsAsFactors = FALSE)
invisible(lapply(seq_len(nrow(zeros)), function(index) {
  settings <- zeros[index, ]
  result <- .multilpa_bounded_variances(matrix(0, settings$k, settings$d),
    seq_len(settings$k), settings$structure, 0.1)
  stopifnot(max(abs(result - 0.1)) < 1e-12)
}))
cat("Independent checks passed:100scalarreferences,100KKTcertificates,23exactoracle scales,40zero-scatter cases.\n")
