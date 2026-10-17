# Growth mixtures with correlated residuals (mixture_regression(
# residual_correlation =)) against lcmm::hlme(cor = AR(time)) and
# hlme(cor = BM(time)) (Proust-Lima, Philipps and Liquet 2017). hlme's
# process covariances, identified by matching its likelihood: AR,
# s^2 exp(-rho |t_j - t_l|) (cor1 = rho, cor2 = s); BM, s^2 min(t_j, t_l)
# (cor1 = s). As in test-growth-mixture-lcmm.R: (1) the same likelihood at
# hlme's estimates, (2) an optimum at least as high as hlme's.

skip_if_not_installed("lcmm")

correlated_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  set.seed(31)
  n <- 220L
  times <- 0:5
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  out <- data.frame(id = rep(seq_len(n), each = length(times)), time = rep(times, n))
  intercept <- stats::rnorm(n)
  process <- unlist(lapply(seq_len(n), function(i) {
    stats::filter(stats::rnorm(length(times), 0, 0.8), 0.6, method = "recursive")
  }))
  person_class <- class[out$id]
  out$y <- c(1, 4)[person_class] + c(0.5, -0.3)[person_class] * out$time +
    intercept[out$id] + process + stats::rnorm(nrow(out), 0, 0.5)
  out <- out[stats::runif(nrow(out)) > 0.1, ]
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  out
})

fit_correlated_hlme <- function(ng, cor) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  one <- eval(bquote(lcmm::hlme(y ~ time, random = ~1, subject = "id", ng = 1,
                                cor = .(cor), data = correlated_data, verbose = FALSE)))
  if (ng == 1L) return(one)
  hlme <- lcmm::hlme
  model <- bquote(hlme(y ~ time, mixture = ~ time, random = ~1, subject = "id",
                       ng = .(ng), cor = .(cor), data = correlated_data, verbose = FALSE))
  eval(bquote(lcmm::gridsearch(.(model), rep = 10, maxiter = 30, minit = one)))
}

configurations <- list(
  list(label = "AR, one class", ng = 1L, kind = "ar", cor = quote(AR(time))),
  list(label = "AR, two classes", ng = 2L, kind = "ar", cor = quote(AR(time))),
  list(label = "Brownian motion, two classes", ng = 2L, kind = "brownian",
       cor = quote(BM(time))))

invisible(lapply(configurations, function(cfg) {
  test_that(sprintf("correlated growth mixture matches lcmm::hlme: %s", cfg$label), {
    theirs <- fit_correlated_hlme(cfg$ng, cfg$cor)
    expect_equal(theirs$conv, 1)
    ours <- quietly(mixture_regression(y ~ time, correlated_data, n_classes = cfg$ng, id = "id",
                               class_level = "group", random = ~1,
                               random_covariance = "equal", variance = "equal",
                               residual_correlation = cfg$kind, n_starts = 4, seed = 2,
                               vcov_type = "none"), "latents_random_boundary")
    expect_true(ours$converged)
    at_theirs <- latents:::.growth_expectation(ours$spec, ours$stats,
                                               params_from_hlme(ours, theirs))$log_likelihood
    expect_equal(at_theirs, theirs$loglik, tolerance = 1e-8)
    # The optimum, within one likelihood (ours). These data put some fits on
    # a boundary (the residual or the random-intercept variance near zero),
    # where hlme goes below the variance floor latents keeps (min_variance
    # times the outcome's variance); its estimates are then projected onto
    # that floor, a point latents could have reached.
    floor <- ours$spec$min_variance
    feasible <- params_from_hlme(ours, theirs)
    feasible$sigma2 <- pmax(feasible$sigma2, floor)
    feasible$random_covariance <- lapply(feasible$random_covariance, function(g) {
      diag(g) <- pmax(diag(g), floor)
      g
    })
    feasible$correlation[["sd"]] <- max(feasible$correlation[["sd"]], sqrt(floor))
    at_feasible <- latents:::.growth_expectation(ours$spec, ours$stats,
                                                 feasible)$log_likelihood
    expect_gte(ours$log_likelihood, at_feasible - 1e-6)
    expect_lt(abs(ours$log_likelihood - at_theirs), 1e-4)
    expect_equal(ours$n_parameters, length(theirs$best))
  })
}))
