# Comparisons with independent references moved from tests/testthat/test-review-1026-kernel-r2.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("K2-1: VEI fits indicators whose scales differ by 1e9, as mclust does", {
  skip_if_not_installed("mclust")
  # A rate (SD about 0.01) and a population (SD about 1e7). The old code
  # refused every start ("EM likelihood decreased beyond numerical roundoff");
  # mclust::me reaches -3803.70304389.
  withr::local_seed(4)
  x <- rbind(cbind(stats::rnorm(150, 0.05, 0.010), stats::rnorm(150, 3e7, 1.0e7)),
             cbind(stats::rnorm(100, 0.09, 0.015), stats::rnorm(100, 6e7, 1.5e7)))
  colnames(x) <- c("rate", "pop")
  z <- cbind(rep(1:0, c(150, 100)), rep(0:1, c(150, 100)))
  # mclust dispatches to meVEI() and mstepVEI() by name in its namespace.
  reference <- do.call(mclust::me, list(
    modelName = "VEI", data = x, z = z,
    control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(1e5, 1e5))),
    envir = asNamespace("mclust"))
  first <- do.call(mclust::mstep, list(modelName = "VEI", data = x, z = z),
                   envir = asNamespace("mclust"))$parameters
  start <- list(means = t(first$mean), profile_probabilities = matrix(first$pro, 1L),
                group_probabilities = 1,
                variances = t(vapply(1:2, function(k) diag(first$variance$sigma[, , k]),
                                     numeric(2))))
  fit <- suppressMessages(suppressWarnings(multilpa(
    as.data.frame(x), colnames(x), id = NULL, n_profiles = 2, model = "VEI",
    start = start, n_starts = 1, tol = 1e-12, max_iter = 20000,
    acceleration = "none")))
  expect_equal(fit$log_likelihood, reference$loglik, tolerance = 1e-10)
  expect_equal(unname(fit$variances),
               unname(t(vapply(1:2, function(k) diag(reference$parameters$variance$sigma[, , k]),
                               numeric(2)))), tolerance = 1e-5)
})

test_that("K2-1, K2-4: VEI, EVE and VVE M-steps are mclust's and unit-free at 1e-9 and 1e9", {
  skip_if_not_installed("mclust")
  # Measured after the fix (iris, soft species partition): VEI equals
  # mclust::mstep to 4.3e-10 relative at units 1e-9, 1 and 1e9; EVE and VVE
  # reach an objective 1.4e-5 and 0.80 below mclust's at 1 and 1e9, and
  # mclust returns no finite covariance at 1e-9. Before the fix VEI at 1e9
  # (bound scaled to the units) and EVE/VVE at 1e-9 were off by 53%.
  x0 <- as.matrix(iris[, 1:4])
  z <- mclust::unmap(as.integer(iris$Species)) * 0.8 + 0.2 / 3
  weights <- colSums(z)
  objective <- function(sigma, scatter) sum(vapply(seq_along(weights), function(k) {
    weights[k] * as.numeric(determinant(sigma[, , k])$modulus) +
      sum(solve(sigma[, , k]) * scatter[[k]])
  }, numeric(1)))
  m_step <- function(model, scatter, floor) {
    if (model == "VEI") {
      v <- .multilpa_structure_variances(.multilpa_scatter_diagonals(scatter), weights,
                                         "VEI", floor)
      return(array(unlist(lapply(seq_along(weights), function(k) diag(v[k, ]))),
                   c(4, 4, 3)))
    }
    .multilpa_structure_covariances(scatter, weights, model, floor)
  }
  scatter_at <- function(units) {
    x <- x0 * units
    centres <- sweep(crossprod(z, x), 1L, weights, "/")
    lapply(seq_along(weights), function(k) {
      r <- sweep(x, 2L, centres[k, ])
      crossprod(r, r * z[, k])
    })
  }
  invisible(lapply(c("VEI", "EVE", "VVE"), function(model) {
    base <- m_step(model, scatter_at(1), 1e-6)
    invisible(lapply(c(1e-9, 1e9), function(units) {
      scaled <- m_step(model, scatter_at(units), 1e-6 * units^2)
      expect_lt(max(abs(scaled / units^2 - base)) / max(abs(base)), 1e-9)
      reference <- do.call(mclust::mstep, list(modelName = model, data = x0 * units,
                                               z = z),
                           envir = asNamespace("mclust"))$parameters
      sigma <- reference$variance$sigma
      if (all(is.finite(sigma))) {
        if (model == "VEI") {
          expect_lt(max(abs(scaled - sigma)) / max(abs(sigma)), 1e-8)
        } else {
          expect_lte(objective(scaled, scatter_at(units)) -
                       objective(sigma, scatter_at(units)), 1e-8)
        }
      }
    }))
  }))
})
