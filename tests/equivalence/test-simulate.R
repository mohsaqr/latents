# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# simulate() for mixture regressions and growth mixtures. The draws (classes, random effects, residuals) are made
# inside the function with only nsim and seed as arguments, so the evidence is statistical: Monte Carlo summaries
# of the simulated outcomes against the moments the fitted model implies, written out in base R from coef().
#
# Each draw takes a class from the membership probabilities (the prior, not the posterior), then a normal outcome
# around that class's line (with, for a growth mixture, a person-level random intercept), on the fitted design.
# For row r the class-k mean is mu_rk, the class-k variance s_rk^2, so the implied moments are
#   E[y] = sum_k pi_k mu_rk,  E[y^2] = sum_k pi_k (mu_rk^2 + s_rk^2),
#   E[y^4] = sum_k pi_k (mu_rk^4 + 6 mu_rk^2 s_rk^2 + 3 s_rk^4)   (normal moments),
# and the variance of a mean of independent draws follows from them. Tolerances are four standard errors of each
# summary, computed from these implied variances and the number of draws before any run; summaries are pooled only
# over draws that are independent (rows of a regression; persons at one time point of a growth model).

seeds <- 1:5
simulations <- 200L
mc_reference <- "independent: moments implied by the fitted mixture, base R (Monte Carlo summary)"

#' Class shares from intercept-only membership logits (class 1 is the baseline)
class_shares <- function(estimates, k) {
  logits <- c(0, vapply(seq_len(k)[-1L], \(j) estimates[[sprintf("membership.class_%d.(Intercept)", j)]],
                        numeric(1)))
  exp(logits) / sum(exp(logits))
}

#' Mixture moments of one row: a matrix of class means and one of class variances, rows by classes
mixture_moments <- function(means, variances, shares) {
  weight <- matrix(shares, nrow(means), ncol(means), byrow = TRUE)
  first <- rowSums(weight * means)
  second <- rowSums(weight * (means^2 + variances))
  fourth <- rowSums(weight * (means^4 + 6 * means^2 * variances + 3 * variances^2))
  data.frame(first = first, second = second, fourth = fourth)
}

regression_simulation_data <- function() {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(3)
  n <- 300
  x <- stats::runif(n)
  in_second <- stats::rbinom(n, 1, 0.4) == 1
  data.frame(x = x, y = ifelse(in_second, 1 + 2 * x, 4 - x) + stats::rnorm(n, sd = 0.5))
}

growth_simulation_data <- function() {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(11)
  n <- 150
  rising <- stats::rbinom(n, 1, 0.6) == 1
  out <- data.frame(id = rep(seq_len(n), each = 5), time = rep(0:4, n))
  intercept <- ifelse(rising, 2, 6)[out$id]
  slope <- ifelse(rising, 0.5, -0.5)[out$id]
  out$y <- intercept + slope * out$time + stats::rnorm(n, 0, 0.8)[out$id] + stats::rnorm(nrow(out), 0, 0.6)
  out
}

test_that("simulated mixture-regression outcomes have the fitted model's mean, cross moment and second moment", {
  data <- regression_simulation_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 2, n_starts = 5, seed = 1)
  estimates <- coef(fit)
  at <- \(name, j) estimates[[sprintf(name, j)]]
  means <- vapply(1:2, \(j) at("coefficient.class_%d.(Intercept)", j) + at("coefficient.class_%d.x", j) * data$x,
                  numeric(nrow(data)))
  variances <- matrix(vapply(1:2, \(j) exp(2 * at("log_sigma.class_%d", j)), numeric(1)), nrow(data), 2L,
                      byrow = TRUE)
  implied <- mixture_moments(means, variances, class_shares(estimates, 2L))
  n <- nrow(data)
  draws <- n * simulations
  row_variance <- implied$second - implied$first^2
  tolerance_mean <- 4 * sqrt(sum(row_variance) * simulations) / draws
  tolerance_cross <- 4 * sqrt(sum(data$x^2 * row_variance) * simulations) / draws
  tolerance_second <- 4 * sqrt(sum(implied$fourth - implied$second^2) * simulations) / draws
  ours <- do.call(rbind, lapply(seeds, \(s) {
    outcomes <- as.matrix(simulate(fit, nsim = simulations, seed = s))
    data.frame(seed = s, mean = mean(outcomes), cross = mean(data$x * outcomes), second = mean(outcomes^2))
  }))
  theirs <- data.frame(seed = seeds, mean = mean(implied$first), cross = mean(data$x * implied$first),
                       second = mean(implied$second))
  covers <- "simulate.latents_mixture_regression()"
  expect_reference(subset(ours, select = c(seed, mean)), subset(theirs, select = c(seed, mean)),
                   reference = mc_reference, covers = covers, quantity = "Monte Carlo summary: mean outcome",
                   tolerance = tolerance_mean, key = "seed")
  expect_reference(subset(ours, select = c(seed, cross)), subset(theirs, select = c(seed, cross)),
                   reference = mc_reference, covers = covers, quantity = "Monte Carlo summary: mean of x times outcome",
                   tolerance = tolerance_cross, key = "seed")
  expect_reference(subset(ours, select = c(seed, second)), subset(theirs, select = c(seed, second)),
                   reference = mc_reference, covers = covers, quantity = "Monte Carlo summary: mean squared outcome",
                   tolerance = tolerance_second, key = "seed")
})

test_that("simulated growth-mixture outcomes have the fitted model's mean and second moment at each time", {
  data <- growth_simulation_data()
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group", random = ~ 1,
                            n_starts = 5, seed = 1)
  estimates <- coef(fit)
  at <- \(name, j) estimates[[sprintf(name, j)]]
  times <- 0:4
  persons <- length(unique(data$id))
  # A person's rows share the class and the random intercept, so draws are pooled over persons at one time only.
  means <- vapply(1:2, \(j) at("coefficient.class_%d.(Intercept)", j) + at("coefficient.class_%d.time", j) * times,
                  numeric(length(times)))
  variances <- matrix(vapply(1:2, \(j) exp(2 * at("log_sigma.class_%d", j)) +
                               exp(2 * at("random.log_cholesky.class_%d.(Intercept):(Intercept)", j)), numeric(1)),
                      length(times), 2L, byrow = TRUE)
  implied <- transform(mixture_moments(means, variances, class_shares(estimates, 2L)), time = times)
  draws <- persons * simulations
  tolerance_mean <- 4 * sqrt(max(implied$second - implied$first^2) / draws)
  tolerance_second <- 4 * sqrt(max(implied$fourth - implied$second^2) / draws)
  ours <- do.call(rbind, lapply(seeds, \(s) {
    outcomes <- as.matrix(simulate(fit, nsim = simulations, seed = s))
    by_time <- split(seq_len(nrow(data)), data$time)
    data.frame(seed = s, time = as.integer(names(by_time)),
               mean = vapply(by_time, \(rows) mean(outcomes[rows, ]), numeric(1)),
               second = vapply(by_time, \(rows) mean(outcomes[rows, ]^2), numeric(1)))
  }))
  theirs <- merge(data.frame(seed = seeds), data.frame(time = times, mean = implied$first, second = implied$second))
  covers <- "simulate.latents_growth_mixture()"
  expect_reference(subset(ours, select = c(seed, time, mean)), subset(theirs, select = c(seed, time, mean)),
                   reference = mc_reference, covers = covers, quantity = "Monte Carlo summary: mean outcome per time",
                   tolerance = tolerance_mean, key = c("seed", "time"))
  expect_reference(subset(ours, select = c(seed, time, second)), subset(theirs, select = c(seed, time, second)),
                   reference = mc_reference, covers = covers,
                   quantity = "Monte Carlo summary: mean squared outcome per time", tolerance = tolerance_second,
                   key = c("seed", "time"))
})
