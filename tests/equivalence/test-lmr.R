# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-lmr.R,
# so that they run in CI's equivalence job rather than under R CMD check.

# The internal Imhof tail (.latents_imhof_upper) gives the VLMR and adjusted
# LMR p-values of get_results(grid, "lmr").
imhof_covers <- "get_results.multilpa_enumeration()"
imhof_note <- "Package side: the internal weighted chi-square tail (.latents_imhof_upper) behind the VLMR and adjusted LMR p-values of get_results(grid, \"lmr\")."

test_that("the Imhof tail matches exact weighted chi-square probabilities", {
  # One weight, or equal weights, give a scaled chi-square. The original test
  # compared each statistic alone with testthat's relative 1e-8; kept, per number
  # (the smallest probability here is 4e-8).
  statistics <- c(0.05, 0.5, 3, 9, 30)
  labels <- sprintf("q = %g", statistics)
  imhof <- function(weights) {
    vapply(statistics, \(q) latents:::.latents_imhof_upper(q, weights), numeric(1))
  }
  expect_reference(imhof(1), stats::pchisq(statistics, 1, lower.tail = FALSE),
                   reference = "stats", covers = imhof_covers,
                   quantity = "upper tail, one unit weight (chi-square 1)",
                   tolerance = 1e-8, measure = "relative", key = labels, note = imhof_note)
  expect_reference(imhof(rep(1, 3)), stats::pchisq(statistics, 3, lower.tail = FALSE),
                   reference = "stats", covers = imhof_covers,
                   quantity = "upper tail, three unit weights (chi-square 3)",
                   tolerance = 1e-8, measure = "relative", key = labels, note = imhof_note)
  expect_reference(imhof(2.5), stats::pchisq(statistics / 2.5, 1, lower.tail = FALSE),
                   reference = "stats", covers = imhof_covers,
                   quantity = "upper tail, one weight of 2.5 (scaled chi-square 1)",
                   tolerance = 1e-8, measure = "relative", key = labels, note = imhof_note)
  # Two weights, against the convolution of their distributions: P(a X + b Y >
  # q) with y = t^2, so the chi-square(1) density times dy is 2 dnorm(t) dt,
  # smooth at zero. One weight small, so the tail decays slowly; one
  # negative, as an eigenvalue of Vuong's matrix can be.
  two <- function(q, a, b) {
    f <- function(t) stats::pchisq(pmax((q - b * t^2) / a, 0), 1, lower.tail = FALSE) *
      2 * stats::dnorm(t)
    if (b > 0) {
      kink <- sqrt(q / b)
      stats::integrate(f, 0, kink, rel.tol = 1e-12)$value +
        2 * stats::pnorm(kink, lower.tail = FALSE)
    } else stats::integrate(f, 0, Inf, rel.tol = 1e-12)$value
  }
  grid <- expand.grid(q = c(1e-6, 0.5, 5.5, 30, 80), a = c(0.9, 3, 0.01),
                      b = c(0.05, -0.4, -3, 50))
  grid_labels <- with(grid, sprintf("q = %g, weights %g and %g", q, a, b))
  # Absolute: the tails run down to 1e-16 and below, where the Imhof
  # integral's accuracy (about 1e-12) is absolute, not relative. The
  # original test's tolerance, 1e-10.
  expect_reference(mapply(function(q, a, b) latents:::.latents_imhof_upper(q, c(a, b)),
                          grid$q, grid$a, grid$b),
                   mapply(two, grid$q, grid$a, grid$b),
                   reference = "independent: convolution of two scaled chi-squares, stats::integrate()",
                   covers = imhof_covers, quantity = "upper tail, two weights",
                   tolerance = 1e-10, measure = "absolute", key = grid_labels, note = imhof_note)
  # A statistic of zero with weights of both signs has no oscillation to
  # extrapolate: P(a X_m - b Y_k > 0) = P(F(m, k) > b k / (a m)), exactly.
  # Extrapolating partial sums returned 5.6e-7 for 4.2e-7 here.
  cases <- expand.grid(m = c(1, 3), b = c(0.1, 10, 1e4), k = c(1, 3, 6))
  # Absolute, as above: the F tails run down to small values where the
  # accuracy is absolute. The original test's tolerance, 1e-10.
  expect_reference(mapply(function(m, b, k) latents:::.latents_imhof_upper(0, c(rep(1, m), rep(-b, k))),
                          cases$m, cases$b, cases$k),
                   with(cases, stats::pf(b * k / m, m, k, lower.tail = FALSE)),
                   reference = "stats", covers = imhof_covers,
                   quantity = "probability of a positive difference of chi-squares (F tail)",
                   tolerance = 1e-10, measure = "absolute",
                   key = with(cases, sprintf("m = %g, b = %g, k = %g", m, b, k)), note = imhof_note)
  # The accuracy is absolute: about 1e-12, here on a probability of 4e-7.
  expect_reference(latents:::.latents_imhof_upper(0, c(1, -1e4, -1e4, -1e4)),
                   stats::pf(3e4, 1, 3, lower.tail = FALSE),
                   reference = "stats", covers = imhof_covers,
                   quantity = "probability of a positive difference, weights 1 and three of -1e4 (F tail)",
                   tolerance = 1e-12, measure = "absolute", note = imhof_note)
  # The probability does not depend on the weights' units.
  sizes <- c(1e-12, 1e-6, 1e6)
  expect_reference(vapply(sizes, \(size) latents:::.latents_imhof_upper(0, size * c(1, -2)), numeric(1)),
                   rep(stats::pf(2, 1, 1, lower.tail = FALSE), length(sizes)),
                   reference = "stats", covers = imhof_covers,
                   quantity = "probability of a positive difference, weights 1 and -2 in other units (F tail)",
                   tolerance = 1e-10, measure = "relative", key = sprintf("units %g", sizes),
                   note = imhof_note)
  invisible(lapply(sizes, \(size) {
    expect_equal(latents:::.latents_imhof_upper(5 * size, size * c(2, 1, 0.5, -0.4)),
                 latents:::.latents_imhof_upper(5, c(2, 1, 0.5, -0.4)), tolerance = 1e-10)
  }))
})

test_that("VLMR and adjusted LMR agree with Latent GOLD", {
  skip_on_cran()
  # Latent GOLD 6.1, `latent Cluster nominal 2:3` with one variance per
  # indicator shared by the classes, on these simulated data (2 true classes):
  # LR 8.3027, adjusted 8.0439, reference mean 4.1216 and standard deviation
  # 2.8763, p-values 0.0848 (VLMR) and 0.0941 (adjusted).
  set.seed(2024)
  n <- 500
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  data <- data.frame(y1 = stats::rnorm(n, c(0, 1.6)[class]),
                     y2 = stats::rnorm(n, c(0, 1.4)[class]),
                     y3 = stats::rnorm(n, c(0, 1.2)[class]),
                     y4 = stats::rnorm(n, c(0, 1)[class]))
  models <- quietly(enumerate_lpa(data, paste0("y", 1:4), n_profiles = 2:3,
                                  model = "EEI", n_starts = 50, seed = 1, tol = 1e-12))
  test <- get_results(models, "lmr")
  printed <- c(lr = 8.3027, adjusted_lr = 8.0439, reference_mean = 4.1216,
               reference_sd = 2.8763, vlmr_p = 0.0848, almr_p = 0.0941)
  # Latent GOLD prints four decimals: 5e-5 absolute is agreement.
  expect_reference(vapply(names(printed), \(name) test[[name]], numeric(1)), printed,
                   reference = "Latent GOLD",
                   covers = c("enumerate_lpa()", "get_results.multilpa_enumeration()"),
                   quantity = "VLMR and adjusted LMR, two against three classes, EEI",
                   tolerance = 5e-5, measure = "absolute",
                   note = "Latent GOLD 6.1 output typed into the test (the comment above), from the listing kept in the private equivalence repository: latents/local/equivalence/latentgold-extra/vlmr/null2_equal.lst (syntax null2_equal.lgs, data null2.dat).")
  # The adjustment is LR / (1 + 1 / (d log n)), d the added parameters.
  expect_equal(test$adjusted_lr, test$lr / (1 + 1 / (5 * log(500))), tolerance = 1e-12)
})
