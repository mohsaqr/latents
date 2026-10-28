pkgload::load_all(quiet = TRUE)
scatter <- rbind(c(1e-5, 1), c(1e-4, 100))
reference <- matrix(c(1, 1, 50.5, 50.5), 2)
result <- .multilpa_bounded_variances(scatter, c(1, 1), "EVI", 1)
# The reference has both first dimensions at their floor and their common
# second variance is the pooled second scatter / total membership = 50.5.
# Positive directional derivatives confirm these active constraints.
print(list(actual = result, reference = reference,
  maximum_error = max(abs(result - reference)),
  objective_excess = sum(log(result) + scatter / result) -
    sum(log(reference) + scatter / reference)))
debug_solver <- .multilpa_bounded_variances
body(debug_solver)[[length(body(debug_solver))]] <- quote(list(
  active = active, log_slack = logs - lower, gradient = gradient,
  active_multipliers = solve(tcrossprod(design[active, , drop = FALSE]),
    design[active, , drop = FALSE] %*% gradient), iterations = iteration,
  settled = settled))
print(debug_solver(scatter, c(1, 1), "EVI", 1), digits = 17)
