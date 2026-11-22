# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Mixture regressions against outside programs:
# - enumerate_regressions() on mpg ~ wt (datasets::mtcars) against mixtools
#   (Benaglia et al., 2009): regmixEM()'s EM continued from each candidate's
#   estimates (so the class labels stay latents') reaches the same maximum,
#   and get_tree()'s flows are computed from mixtools' posteriors. The
#   one-class candidate is stats::lm() (test-accessors-selection.R).
# - simulate() of a two-class mixture regression against flexmix (Leisch,
#   2004; Gruen and Leisch, 2008): every simulated outcome is refitted by
#   flexmix, and the mean of flexmix's estimates over the draws is the fitted
#   model's parameters within four Monte Carlo standard errors.
# flexmix is not compared on the grid: its Gaussian M-step divides the
# weighted residual sum of squares by the weighted count less the number of
# coefficients, so its sigma is not the maximum likelihood one (one class:
# 3.0459 against lm()'s 2.9491 = 3.0459 sqrt(30/32)) and its likelihood sits
# below the maximum (-80.0473 against -80.0147).
# regmixEM() stops at an absolute change of 1e-12 in the likelihood; the grid
# here is fitted with tol = 1e-13 (mixture_regression()'s default, 1e-8,
# leaves the two-class EM 9.4e-7 below the maximum on this slow ridge, with
# posteriors 2e-4 away): 1e-8 relative on the likelihood, 1e-5 absolute on
# posterior masses.

flexmix_control <- list(tol = 1e-12, iter.max = 100000L, minprior = 0)

regression_grid <- quietly(enumerate_regressions(mpg ~ wt, datasets::mtcars, n_classes = 1:3,
                                                 n_starts = 5, seed = 1, tol = 1e-13, max_iter = 100000L))

# mixtools' EM from the estimates of latents' k-class candidate (k >= 2).
continued_mixtools <- function(k) {
  estimates <- coef(get_results(regression_grid, "model", n_classes = k))
  at <- \(name, class) estimates[[sprintf(name, class)]]
  beta <- vapply(seq_len(k), \(j) c(at("coefficient.class_%d.(Intercept)", j), at("coefficient.class_%d.wt", j)),
                 numeric(2))
  logits <- c(0, vapply(seq_len(k)[-1L], \(j) at("membership.class_%d.(Intercept)", j), numeric(1)))
  printed <- utils::capture.output(fitted <- mixtools::regmixEM(
    datasets::mtcars$mpg, datasets::mtcars$wt, lambda = exp(logits) / sum(exp(logits)), beta = beta,
    sigma = vapply(seq_len(k), \(j) exp(at("log_sigma.class_%d", j)), numeric(1)), k = k,
    epsilon = 1e-12, maxit = 100000L))
  fitted
}

test_that("the regression grid's maxima are mixtools' (regmixEM continued from latents' estimates)", {
  skip_without_reference("mixtools")
  expected <- do.call(rbind, lapply(2:3, \(k) data.frame(n_classes = k, log_likelihood = continued_mixtools(k)$loglik)))
  grid <- subset(as.data.frame(regression_grid), n_classes >= 2)
  expect_reference(grid, expected, reference = "mixtools",
                   covers = c("enumerate_regressions()", "as.data.frame.latents_regression_enumeration()"),
                   quantity = "maximised log likelihood per candidate (mixtools continued)", tolerance = 1e-8,
                   measure = "relative", key = "n_classes")
  expect_reference(subset(summary(regression_grid), n_classes >= 2), expected, reference = "mixtools",
                   covers = c("enumerate_regressions()", "summary.latents_regression_enumeration()"),
                   quantity = "maximised log likelihood per candidate, summary() (mixtools continued)",
                   tolerance = 1e-8, measure = "relative", key = "n_classes")
  expect_reference(subset(get_results(regression_grid), n_classes >= 2), expected, reference = "mixtools",
                   covers = c("enumerate_regressions()", "get_results.latents_regression_enumeration()"),
                   quantity = "maximised log likelihood per candidate, get_results() (mixtools continued)",
                   tolerance = 1e-8, measure = "relative", key = "n_classes")
})

test_that("get_tree()'s flows are those of mixtools' posteriors", {
  skip_without_reference("mixtools")
  n <- nrow(datasets::mtcars)
  posteriors <- c(list(matrix(1, n, 1L)), lapply(2:3, \(k) continued_mixtools(k)$posterior))
  # The flow from class a of k classes to class b of k + 1: the mean over rows
  # of p_k(a) p_(k+1)(b); the shares are the mean posteriors.
  expected <- do.call(rbind, lapply(1:2, \(k) {
    flow <- crossprod(posteriors[[k]], posteriors[[k + 1L]]) / n
    cells <- expand.grid(from_class = seq_len(k), to_class = seq_len(k + 1L))
    from_share <- colMeans(posteriors[[k]])[cells$from_class]
    to_share <- colMeans(posteriors[[k + 1L]])[cells$to_class]
    mass <- flow[cbind(cells$from_class, cells$to_class)]
    data.frame(from_candidate = sprintf("classes%d", k), from_class = cells$from_class,
               to_candidate = sprintf("classes%d", k + 1L), to_class = cells$to_class,
               from_share = from_share, to_share = to_share, mass = mass,
               share_of_from = mass / from_share, share_of_to = mass / to_share)
  }))
  expect_reference(get_tree(regression_grid), expected, reference = "mixtools",
                   covers = c("get_tree()", "enumerate_regressions()"),
                   quantity = "flows between consecutive candidates, from mixtools' posteriors", tolerance = 1e-5,
                   key = c("from_candidate", "from_class", "to_candidate", "to_class"))
})

# Three hundred rows, two lines (the design of test-simulate.R).
simulation_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(3)
  n <- 300
  x <- stats::runif(n)
  in_second <- stats::rbinom(n, 1, 0.4) == 1
  data.frame(x = x, y = ifelse(in_second, 1 + 2 * x, 4 - x) + stats::rnorm(n, sd = 0.5))
})
simulation_fit <- mixture_regression(y ~ x, simulation_data, n_classes = 2, n_starts = 5, seed = 1)

# flexmix's estimates of one simulated outcome, in coef()'s names, its classes
# ordered to match the fitted model's coefficients.
flexmix_estimates <- function(y, truth) {
  data <- data.frame(x = simulation_data$x, y = y)
  fitted <- withr::with_seed(1L, flexmix::stepFlexmix(y ~ x, data = data, k = 2L, nrep = 3L,
                                                      control = c(flexmix_control, verbose = 0L), verbose = FALSE))
  parameters <- flexmix::parameters(fitted)
  prior <- flexmix::prior(fitted)
  coefficients <- parameters[c("coef.(Intercept)", "coef.x"), , drop = FALSE]
  wanted <- rbind(truth[c("coefficient.class_1.(Intercept)", "coefficient.class_1.x")],
                  truth[c("coefficient.class_2.(Intercept)", "coefficient.class_2.x")])
  same <- sum((t(coefficients) - wanted)^2)
  swapped <- sum((t(coefficients[, 2:1]) - wanted)^2)
  order <- if (same <= swapped) 1:2 else 2:1
  c("coefficient.class_1.(Intercept)" = parameters[["coef.(Intercept)", order[1L]]],
    "coefficient.class_1.x" = parameters[["coef.x", order[1L]]],
    "coefficient.class_2.(Intercept)" = parameters[["coef.(Intercept)", order[2L]]],
    "coefficient.class_2.x" = parameters[["coef.x", order[2L]]],
    "log_sigma.class_1" = log(parameters[["sigma", order[1L]]]),
    "log_sigma.class_2" = log(parameters[["sigma", order[2L]]]),
    "membership.class_2.(Intercept)" = log(prior[[order[2L]]] / prior[[order[1L]]]))
}

test_that("simulate() of a mixture regression: flexmix recovers the fitted parameters from the draws", {
  skip_without_reference("flexmix")
  truth <- coef(simulation_fit)
  draws <- 100L
  # Four Monte Carlo standard errors of a mean of `draws` estimates, from the
  # fit's own standard errors (its data have the simulated design and size),
  # fixed before any draw. The finite-sample bias of the maximum likelihood
  # estimates (O(1/n), below 0.01 for the log sds) is well inside it.
  bound <- 4 * sqrt(diag(vcov(simulation_fit)))[names(truth)] / sqrt(draws)
  invisible(lapply(c(11L, 12L, 13L), \(seed) {
    simulated <- simulate(simulation_fit, nsim = draws, seed = seed)
    estimates <- vapply(simulated, flexmix_estimates, numeric(length(truth)), truth = truth)
    means <- rowMeans(estimates)
    invisible(lapply(names(truth), \(name) {
      expect_reference(means[[name]], truth[[name]], reference = "flexmix",
                       covers = c("simulate.latents_mixture_regression()", "mixture_regression()"),
                       quantity = sprintf("mean flexmix estimate over %d draws, %s, seed %d", draws, name, seed),
                       tolerance = bound[[name]], measure = "absolute",
                       note = "Four Monte Carlo standard errors from vcov() of the fit, fixed before the draws.")
    }))
  }))
})
