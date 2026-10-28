pkgload::load_all(quiet = TRUE)
set.seed(8)
scores <- matrix(stats::rnorm(400 * 5), 400, 5)
jacobian <- matrix(stats::rnorm(5 * 4), 5, 4)
null_scores <- scores %*% jacobian
alternative <- list(scores = scores, information = crossprod(scores))
null <- list(scores = null_scores, information = crossprod(null_scores))
# Reproduce the previous implementation from its stated matrix formula.
both <- crossprod(scores, null_scores)
previous_matrix <- rbind(
  cbind(crossprod(scores) %*% solve(alternative$information),
        both %*% solve(null$information)),
  cbind(-t(both) %*% solve(alternative$information),
        -crossprod(null_scores) %*% solve(null$information)))
previous <- Re(eigen(previous_matrix, only.values = TRUE)$values)
current <- .latents_vuong_weights(alternative, null)
cat("Previous weights:\n"); print(previous, digits = 17)
cat("Symmetric weights:\n"); print(current, digits = 17)
cat("Previous p at zero:", .latents_imhof_upper(0, previous), "\n")
cat("Correct p at zero:", .latents_imhof_upper(0, current), "\n")
cat("Known exact p at zero: 1\n")
stopifnot(sum(current != 0) == 1, abs(current[current != 0] - 1) < 1e-12)
