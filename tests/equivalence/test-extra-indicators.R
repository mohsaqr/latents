# Ordinal (adjacent-category logit) and count (Poisson) indicators.

extra_rows <- function(seed = 1, n = 600) {
  set.seed(seed)
  profile <- rbinom(n, 1, 0.4) + 1
  probabilities <- exp(.latents_ordinal_log_probabilities(c(0.5, 0, -1), c(1.2, 0)))
  levels <- c("never", "rarely", "often", "always")
  data.frame(
    y = rnorm(n, c(0, 1.5)[profile]),
    o = factor(vapply(profile, function(h) sample(levels, 1L, prob = probabilities[h, ]),
                      character(1)), levels = levels),
    k = rpois(n, c(1, 5)[profile]))
}

extra_fit <- function(data, ...) {
  quietly(lpa(data, c("y", "o", "k"), 2, ordinal = "o", count = "k", n_starts = 3,
              seed = 1, ...), "latents_single_level")
}

numeric_gradient <- function(objective, theta, step = 1e-6) {
  vapply(seq_along(theta), function(index) {
    up <- theta
    down <- theta
    up[index] <- up[index] + step
    down[index] <- down[index] - step
    (objective(up) - objective(down)) / (2 * step)
  }, numeric(1))
}

test_that("fits agree with Latent GOLD 6.1", {
  skip_on_cran()
  fixture <- readRDS(test_path("fixtures", "latentgold-ordinal.rds"))
  vapply(names(fixture), function(name) {
    case <- fixture[[name]]
    fit <- quietly(multilpa(case$data, case$vars, if (case$two_level) "g" else NULL,
                            case$n_profiles,
                            n_group_classes = if (case$two_level) case$n_group_classes else 1L,
                            ordinal = case$vars[case$types == "ordinal"],
                            count = case$vars[case$types %in% c("count", "negbin")],
                            count_model = if (any(case$types == "negbin"))
                              "negative_binomial" else "poisson",
                            count_dispersion = case$count_dispersion %||% "varying",
                            n_starts = 10, seed = 1, tol = 1e-10, max_iter = 5000),
                   "latents_single_level")
    # Latent GOLD prints four decimals.
    expect_lt(abs(fit$log_likelihood - case$latent_gold_ll), 1e-4)
    expect_equal(fit$n_parameters, case$latent_gold_npar)
    TRUE
  }, logical(1))
})
