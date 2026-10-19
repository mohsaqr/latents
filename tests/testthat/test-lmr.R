# The Vuong-Lo-Mendell-Rubin test (Vuong 1989; Lo, Mendell and Rubin 2001)
# as Latent GOLD implements it (Vermunt 2024).

test_that("the Imhof tail matches exact weighted chi-square probabilities", {
  # One weight, or equal weights, give a scaled chi-square.
  invisible(lapply(c(0.05, 0.5, 3, 9, 30), \(q) {
    expect_equal(.latents_imhof_upper(q, 1), stats::pchisq(q, 1, lower.tail = FALSE),
                 tolerance = 1e-8)
    expect_equal(.latents_imhof_upper(q, rep(1, 3)),
                 stats::pchisq(q, 3, lower.tail = FALSE), tolerance = 1e-8)
    expect_equal(.latents_imhof_upper(q, 2.5), stats::pchisq(q / 2.5, 1, lower.tail = FALSE),
                 tolerance = 1e-8)
  }))
  # Two weights, against the convolution of their distributions; one weight
  # small, so the tail decays slowly, which adaptive quadrature alone fails on.
  two <- function(q, a, b) {
    stats::integrate(function(y) {
      stats::pchisq(pmax((q - b * y) / a, 0), 1, lower.tail = FALSE) * stats::dchisq(y, 1)
    }, 0, Inf, rel.tol = 1e-12, subdivisions = 5000L)$value
  }
  invisible(lapply(list(c(0.9, 0.05), c(0.01, 3)), \(weights) {
    invisible(lapply(c(0.5, 5.5, 9, 30), \(q) {
      expect_equal(.latents_imhof_upper(q, weights), two(q, weights[1], weights[2]),
                   tolerance = 1e-7)
    }))
  }))
  # A negative weight (an eigenvalue of Vuong's matrix can be negative): the
  # same convolution, which holds for either sign of the second weight.
  invisible(lapply(c(0.5, 3, 9), \(q) {
    expect_equal(.latents_imhof_upper(q, c(1.5, -0.4)), two(q, 1.5, -0.4), tolerance = 1e-7)
  }))
})

test_that("for regular nested models the weights are one per added parameter", {
  # Information equal to the variance of the score, and the smaller model's
  # scores the larger one's mapped by the Jacobian of the nesting: Vuong's
  # matrix then has p1 - p0 unit eigenvalues and zeros, so the reference
  # distribution is chi-square on p1 - p0 degrees of freedom.
  set.seed(8)
  scores <- matrix(stats::rnorm(400 * 5), 400, 5)
  jacobian <- matrix(stats::rnorm(5 * 3), 5, 3)
  information <- crossprod(scores)
  alternative <- list(scores = scores, information = information)
  null <- list(scores = scores %*% jacobian,
               information = t(jacobian) %*% information %*% jacobian)
  lambda <- sort(.latents_vuong_weights(alternative, null), decreasing = TRUE)
  expect_equal(lambda[1:2], c(1, 1), tolerance = 1e-8)
  # The rest are zero to the rounding of an unsymmetric eigendecomposition.
  expect_lt(max(abs(lambda[-(1:2)])), 1e-6)
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
  within <- function(actual, printed) {
    expect_true(abs(actual - printed) <= 5e-5, info = format(actual, digits = 10))
  }
  within(test$lr, 8.3027)
  within(test$adjusted_lr, 8.0439)
  within(test$reference_mean, 4.1216)
  within(test$reference_sd, 2.8763)
  within(test$vlmr_p, 0.0848)
  within(test$almr_p, 0.0941)
  # The adjustment is LR / (1 + 1 / (d log n)), d the added parameters.
  expect_equal(test$adjusted_lr, test$lr / (1 + 1 / (5 * log(500))), tolerance = 1e-12)
})

test_that("the comparison table lists each neighbour pair, with its reason when unavailable", {
  skip_on_cran()
  data <- subset(course_engagement, student <= 60)
  grid <- quietly(enumerate_classes(data, c("browse", "lectures", "forum_read"),
                                    "student", n_profiles = 1:3,
                                    n_group_classes = 1:2, model = "VVI",
                                    n_starts = 2, seed = 1, tol = 1e-10))
  test <- get_results(grid, "lmr")
  expect_s3_class(test, "data.frame")
  expect_named(test, c("direction", "null_candidate", "alternative_candidate",
                       "n_profiles", "n_group_classes", "model", "n_units",
                       "parameters_added", "lr", "vlmr_p", "adjusted_lr",
                       "almr_p", "reference_mean", "reference_sd", "note"))
  expect_setequal(unique(test$direction), c("profiles", "group_classes"))
  # Every comparison adds one profile or one group class to its null model.
  table <- as.data.frame(grid)
  null_rows <- match(test$null_candidate, table$candidate)
  alternative_rows <- match(test$alternative_candidate, table$candidate)
  added <- (table$n_profiles[alternative_rows] - table$n_profiles[null_rows]) +
    (table$n_group_classes[alternative_rows] - table$n_group_classes[null_rows])
  expect_true(all(added == 1L))
  computed <- !is.na(test$lr)
  expect_true(all(test$vlmr_p[computed] >= 0 & test$vlmr_p[computed] <= 1))
  expect_true(all(test$adjusted_lr[computed] < test$lr[computed]))
  expect_true(all(test$n_units[computed] == 60L))
  # The table needs every candidate's information, so "all" leaves it out.
  expect_false("lmr" %in% names(get_results(grid, "all")))
  single <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                  n_profiles = 2, model = "VVI", n_starts = 1,
                                  seed = 1))
  expect_error(get_results(single, "lmr"), class = "latents_no_lmr")
})
