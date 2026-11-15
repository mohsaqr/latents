# Comparisons with independent references moved from tests/testthat/test-growth-r3step.R,
# so that they run in CI's equivalence job rather than under R CMD check.

r3step_fixture <- function(n = 240L, seed = 8L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  times <- 0:4
  x <- stats::rnorm(n)
  class <- 1L + (stats::runif(n) < stats::plogis(-0.3 + 1.0 * x))
  out <- data.frame(id = rep(seq_len(n), each = length(times)), time = rep(times, n),
                    x = rep(x, each = length(times)),
                    w = stats::rnorm(n * length(times)))
  out$y <- c(2, 3.5)[class[out$id]] + c(0.6, -0.4)[class[out$id]] * out$time +
    rep(stats::rnorm(n, 0, 0.8), each = length(times)) +
    stats::rnorm(nrow(out), 0, 0.9)
  out
}

r3step_fit <- function(data, ...) {
  mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group",
                     random = ~ 1, n_starts = 2, seed = 1, vcov_type = "none", ...)
}

test_that("r3step on a growth mixture maximizes the fixed-error likelihood", {
  data <- r3step_fixture()
  fit <- r3step_fit(data)
  ours <- r3step(fit, data, "x", vcov_type = "observed")
  expect_identical(ours$outcome, c("class_1", "class_1"))
  expect_identical(ours$level, c("persons", "persons"))
  # Written out from Vermunt (2010): the modal class m_i, error matrix
  # D[k, j] = sum_i p_ik [m_i = j] / sum_i p_ik, and the log likelihood
  # sum_i log sum_k pi_k(x_i) D[k, m_i], logits against the last class.
  posterior <- predict(fit, type = "posterior")
  p <- as.matrix(posterior[c("posterior_class_1", "posterior_class_2")])
  modal <- max.col(p, ties.method = "first")
  errors <- t(vapply(1:2, function(k) {
    vapply(1:2, function(j) sum(p[modal == j, k]), numeric(1)) / sum(p[, k])
  }, numeric(2)))
  x <- data$x[!duplicated(data$id)]
  negative_log_likelihood <- function(theta) {
    eta <- theta[1L] + theta[2L] * x
    prior <- cbind(stats::plogis(eta), 1 - stats::plogis(eta))
    -sum(log(rowSums(prior * t(errors)[modal, ])))
  }
  reference <- stats::optim(c(0, 0), negative_log_likelihood, method = "BFGS",
                            control = list(reltol = 1e-14), hessian = TRUE)
  expect_equal(ours$estimate, reference$par, tolerance = 1e-5)
  expect_equal(ours$standard_error, sqrt(diag(solve(reference$hessian))),
               tolerance = 1e-3)
})

test_that("the stored growth simulation: r3step coverage meets its predeclared gate", {
  # validation/growth-simulation.R, run on a registry fixed before the first
  # replicate (validation/GROWTH_SIMULATION.md). The growth-mixture article
  # reports this coverage; the gate (95% intervals covering in 0.925-0.975 of
  # the datasets, at least 98% available) was set before the run. Under weak
  # separation the coverage falls below it, and the article says so.
  path <- test_path("..", "..", "validation", "growth-simulation-r3step.csv")
  skip_if_not(file.exists(path), "the stored growth simulation is not available")
  rows <- utils::read.csv(path, stringsAsFactors = FALSE)
  registry <- utils::read.csv(test_path("..", "..", "validation",
                                        "growth-simulation-registry.csv"))
  expect_setequal(unique(rows$cell), registry$cell)
  replicates <- max(rows$replicate)
  expect_identical(replicates, 1000L)
  for (cell in registry$cell) {
    own <- rows[rows$cell == cell, ]
    available <- sum(!is.na(own$covered)) / replicates
    coverage <- mean(own$covered, na.rm = TRUE)
    if (identical(cell, "weak")) {
      expect_lt(coverage, 0.925, label = cell)
    } else {
      expect_gte(available, 0.98, label = cell)
      expect_gte(coverage, 0.925, label = cell)
      expect_lte(coverage, 0.975, label = cell)
    }
  }
})
