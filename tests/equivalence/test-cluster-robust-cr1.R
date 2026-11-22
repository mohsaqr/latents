# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Cluster-robust variances against sandwich::vcovCL(), whose default
# adjustment is the CR1 factor G / (G - 1).

sandwich_note <- "Package side: the internal sandwich (.inference_sandwich) behind every vcov_type = \"robust\" of latents, here that of mixture_regression()."

test_that("the sandwich equals sandwich::vcovCL() with its default adjustment", {
  skip_without_reference("sandwich")
  set.seed(5)
  groups <- rep(seq_len(12), each = 6)
  x <- cbind(1, stats::rnorm(72))
  y <- drop(x %*% c(1, 0.4)) + rep(stats::rnorm(12), each = 6) + stats::rnorm(72)
  model <- stats::lm(y ~ x[, 2L])
  bread <- solve(crossprod(x))
  scores <- rowsum(x * stats::residuals(model), groups, reorder = FALSE)
  ours <- latents:::.inference_sandwich(bread, crossprod(scores), nrow(scores))
  reference <- sandwich::vcovCL(model, cluster = groups, type = "HC0")
  expect_reference(unname(ours), unname(reference), reference = "sandwich",
                   covers = "mixture_regression()", quantity = "CR1 sandwich of a linear model",
                   tolerance = 1e-10, measure = "relative",
                   note = sandwich_note)
})

test_that("a one-class robust regression reproduces sandwich::vcovCL()", {
  skip_without_reference("sandwich")
  set.seed(11)
  data <- data.frame(id = rep(seq_len(25), each = 8), x = stats::rnorm(200))
  data$y <- 1 + 0.5 * data$x + rep(stats::rnorm(25), each = 8) + stats::rnorm(200)
  fit <- mixture_regression(y ~ x, data = data, n_classes = 1, id = "id",
                            vcov_type = "robust", n_starts = 1, seed = 1)
  ours <- as.data.frame(get_results(fit, "coefficients"))
  reference <- sqrt(diag(sandwich::vcovCL(stats::lm(y ~ x, data), cluster = ~id,
                                          type = "HC0")))
  # 1e-4 relative, the original test's tolerance (its reason was not recorded).
  expect_reference(ours[c("term", "std_error")],
                   data.frame(term = c("(Intercept)", "x"), std_error = unname(reference)),
                   reference = "sandwich",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()",
                              "as.data.frame.latents_table()"),
                   quantity = "robust standard errors, one-class regression",
                   tolerance = 1e-4, measure = "relative", key = "term")
})


# Moved from tests/testthat/test-cluster-robust-cr1.R: comparisons with independent references.

# A hand-built CR1 reference: G / (G - 1) (X'X)^-1 M (X'X)^-1 with M the
# cross-product of per-cluster score sums. The comparison with
# sandwich::vcovCL() itself is at the top of this file.
cr1_reference <- function(x, residuals, groups) {
  scores <- rowsum(x * residuals, groups, reorder = FALSE)
  n_groups <- nrow(scores)
  bread <- solve(crossprod(x))
  n_groups / (n_groups - 1) * bread %*% crossprod(scores) %*% bread
}

test_that("the sandwich is the hand-built CR1 variance", {
  set.seed(5)
  groups <- rep(seq_len(12), each = 6)
  x <- cbind(1, stats::rnorm(72))
  y <- drop(x %*% c(1, 0.4)) + rep(stats::rnorm(12), each = 6) + stats::rnorm(72)
  residuals <- stats::residuals(stats::lm(y ~ x[, 2L]))
  scores <- rowsum(x * residuals, groups, reorder = FALSE)
  ours <- latents:::.inference_sandwich(solve(crossprod(x)), crossprod(scores), nrow(scores))
  expect_reference(unname(ours), unname(cr1_reference(x, residuals, groups)),
                   reference = "independent: hand-built CR1 sandwich, base R",
                   covers = "mixture_regression()", quantity = "CR1 sandwich of a linear model, against the hand-built CR1",
                   tolerance = 1e-12, measure = "relative", note = sandwich_note)
})

test_that("a one-class robust regression gives the hand-built CR1 standard errors", {
  skip_on_cran()
  set.seed(11)
  data <- data.frame(id = rep(seq_len(25), each = 8), x = stats::rnorm(200))
  data$y <- 1 + 0.5 * data$x + rep(stats::rnorm(25), each = 8) + stats::rnorm(200)
  fit <- mixture_regression(y ~ x, data = data, n_classes = 1, id = "id",
                            vcov_type = "robust", n_starts = 1, seed = 1)
  ours <- as.data.frame(get_results(fit, "coefficients"))
  x <- cbind(1, data$x)
  reference <- sqrt(diag(cr1_reference(x, stats::residuals(stats::lm(y ~ x, data)), data$id)))
  # 1e-4 relative, the original test's tolerance (its reason was not recorded).
  expect_reference(ours[c("term", "std_error")],
                   data.frame(term = c("(Intercept)", "x"), std_error = reference),
                   reference = "independent: hand-built CR1 sandwich, base R",
                   covers = c("mixture_regression()", "get_results.latents_mixture_regression()",
                              "as.data.frame.latents_table()"),
                   quantity = "robust standard errors, one-class regression, against the hand-built CR1",
                   tolerance = 1e-4, measure = "relative", key = "term")
})
