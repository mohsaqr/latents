pkgload::load_all(quiet = TRUE)
probe <- function(separation) {
  set.seed(83)
  scores <- matrix(rnorm(400 * 5), 400, 5)
  null_scores <- cbind(scores[, 1], scores[, 1] + separation * scores[, 2], scores[, 3:4])
  null_information <- crossprod(null_scores)
  a <- list(scores = scores, information = crossprod(scores))
  b <- list(scores = null_scores, information = null_information)
  lambda <- .latents_vuong_weights(a, b)
  scale <- 1 / sqrt(diag(null_information))
  information <- null_information * tcrossprod(scale)
  root <- chol(information)
  inverse <- backsolve(root, diag(ncol(information)))
  gamma <- (ncol(root) + 1) * .Machine$double.eps /
    (1 - (ncol(root) + 1) * .Machine$double.eps)
  error <- gamma * crossprod(abs(root))
  propagated <- crossprod(abs(inverse), error %*% abs(inverse))
  data.frame(separation = separation, condition = kappa(information, exact = TRUE),
    max_null_weight = max(abs(lambda[abs(lambda - 1) > 0.1])),
    propagated_bound = norm(propagated, "2"), nonzero_weights = sum(lambda != 0))
}
print(do.call(rbind, lapply(c(1, .1, .01, .001, .0001), probe)), digits = 14)
