# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Ordinal (adjacent-category logit) and count (Poisson) indicators.

extra_rows <- function(seed = 1, n = 600) {
  set.seed(seed)
  profile <- rbinom(n, 1, 0.4) + 1
  probabilities <- exp(latents:::.latents_ordinal_log_probabilities(c(0.5, 0, -1), c(1.2, 0)))
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
  # The fixture holds the data of the seven cases of
  # equivalence/latentgold-ordinal/ and the log likelihood and parameter count
  # of each Latent GOLD 6.1 listing (equivalence/latentgold-ordinal/returned/),
  # read by equivalence/latentgold-ordinal/make-fixture.R. It has no provenance
  # header.
  fixture <- readRDS(reference_data("fixtures/latentgold-ordinal.rds"))
  fits <- lapply(fixture, function(case) {
    quietly(multilpa(case$data, case$vars, if (case$two_level) "g" else NULL,
                     case$n_profiles,
                     n_group_classes = if (case$two_level) case$n_group_classes else 1L,
                     ordinal = case$vars[case$types == "ordinal"],
                     count = case$vars[case$types %in% c("count", "negbin")],
                     count_model = if (any(case$types == "negbin"))
                       "negative_binomial" else "poisson",
                     count_dispersion = case$count_dispersion %||% "varying",
                     n_starts = 10, seed = 1, tol = 1e-10, max_iter = 5000),
            "latents_single_level")
  })
  note <- "Latent GOLD 6.1 listings equivalence/latentgold-ordinal/returned/<case>.lst, stored in tests/equivalence/fixtures/latentgold-ordinal.rds by equivalence/latentgold-ordinal/make-fixture.R."
  covers <- c("multilpa()", "logLik.multilpa()")
  # Latent GOLD prints four decimals.
  expect_reference(vapply(fits, \(fit) as.numeric(logLik(fit)), numeric(1)),
                   vapply(fixture, \(case) case$latent_gold_ll, numeric(1)),
                   reference = "Latent GOLD", covers = covers,
                   quantity = "log likelihood, ordinal and count indicators",
                   tolerance = 1e-4, measure = "absolute", note = note)
  expect_reference(vapply(fits, \(fit) as.numeric(attr(logLik(fit), "df")), numeric(1)),
                   vapply(fixture, \(case) case$latent_gold_npar, numeric(1)),
                   reference = "Latent GOLD", covers = covers,
                   quantity = "number of parameters, ordinal and count indicators",
                   tolerance = 0, note = note)
})


# Adjacent-category ordinal log probabilities from the definition
# log P(k) / P(k - 1) = a_k - a_(k-1) + eta_c (a_0 = 0), one row per class:
# written here so the optimiser's objective shares no code with the package.
hand_ordinal_log_probabilities <- function(intercepts, locations) {
  scores <- outer(locations, seq_len(length(intercepts) + 1L) - 1L) +
    matrix(c(0, intercepts), length(locations), length(intercepts) + 1L, byrow = TRUE)
  scores - log(rowSums(exp(scores)))
}

test_that("the ordinal M-step maximizes its concave objective", {
  set.seed(2)
  counts <- matrix(runif(15, 1, 40), 3L, 5L)
  solved <- latents:::.latents_ordinal_maximize(counts, numeric(4), numeric(3))
  objective <- function(theta) {
    -sum(counts * hand_ordinal_log_probabilities(theta[1:4], c(theta[5:6], 0)))
  }
  reference <- optim(numeric(6), objective, method = "BFGS",
                     control = list(reltol = 1e-15, maxit = 10000))
  ours <- objective(c(solved$intercepts, solved$locations[1:2]))
  expect_lte(ours, reference$value + 1e-10)
  # BFGS stops at its own criterion on a flat concave surface: 1e-5 relative,
  # the original test's tolerance.
  expect_reference(c(solved$intercepts, solved$locations[1:2]), reference$par,
                   reference = "stats",
                   covers = c("lpa()", "multilpa()"),
                   quantity = "adjacent-category ordinal M-step: intercepts and two class locations",
                   tolerance = 1e-5, measure = "relative",
                   key = c(sprintf("intercept %d", 1:4), sprintf("location %d", 1:2)),
                   note = "Package side: the internal ordinal M-step (.latents_ordinal_maximize) of lpa()/multilpa(ordinal =); reference: stats::optim() BFGS maximising the adjacent-category likelihood written in this file (hand_ordinal_log_probabilities()).")
  expect_identical(solved$locations[3L], 0)
})
