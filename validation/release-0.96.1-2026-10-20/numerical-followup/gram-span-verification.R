pkgload::load_all(quiet = TRUE)
# Equal score-Gram information is not an assertion that the spaces nest.
# Rotated orthonormal subspaces have exact weights +/- sin(angle), plus
# one unit weight for the additional alternative coordinate.
set.seed(917)
basis <- qr.Q(qr(matrix(stats::rnorm(400 * 9), 400, 9)))
results <- lapply(c(0.2, 1e-6, 1e-10), function(angle) {
  alternative_scores <- basis[, 1:5]
  null_scores <- cos(angle) * basis[, 1:4] + sin(angle) * basis[, 6:9]
  alternative <- list(scores = alternative_scores, information = crossprod(alternative_scores))
  null <- list(scores = null_scores, information = crossprod(null_scores))
  lambda <- sort(.latents_vuong_weights(alternative, null))
  expected <- sort(c(rep(-sin(angle), 4), rep(sin(angle), 4), 1))
  error <- max(abs(lambda - expected))
  stopifnot(error < 1e-13)
  data.frame(angle = angle, max_weight_error = error, retained_weights = sum(lambda != 0))
})
print(do.call(rbind, results), digits = 17)
# Independently known ill-conditioned genuine negatives must be preserved or
# refused. Their nonzero probability cannot be silently replaced by one.
set.seed(83)
scores <- qr.Q(qr(matrix(stats::rnorm(400 * 5), 400, 5)))
results <- lapply(list(c(1e-4, 1e-8), c(1e-3, 1e-10)), function(case) {
  transform <- diag(4)
  transform[1, 2] <- 1
  transform[2, 2] <- case[1]
  delta <- case[2]
  alternative <- list(scores = scores, information = diag(5))
  null <- list(scores = scores[, 1:4] %*% transform * sqrt(1 + delta),
               information = crossprod(transform))
  result <- tryCatch(.latents_vuong_weights(alternative, null),
                     latents_inaccurate_probability = identity)
  stopifnot(inherits(result, "latents_inaccurate_probability"))
  data.frame(separation = case[1], delta = delta,
              expected_probability = stats::pf(4 * delta, 1, 4, lower.tail = FALSE),
              disposition = class(result)[1L])
})
print(do.call(rbind, results), digits = 17)
