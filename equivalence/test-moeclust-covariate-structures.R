# Membership-covariate models in every covariance structure against
# MoEClust::MoE_clust(gating = ~ z), which fits the mclust structures with a
# multinomial logit on class membership (single level). Both programs search
# with 20 starts and tolerances near 1e-14 on simulated data with two
# profiles, three indicators and a covariate that moves membership.
#
# For each structure, latents' log likelihood may not fall below MoEClust's.
# Where MoEClust's EM rose to its stop, the two agree on the log likelihood,
# the parameter count and the class-specific slope (|beta_2 - beta_1|, which
# does not depend on the reference class). Where MoEClust stopped on a
# decrease (its VVE orientation step is not monotone), latents must be higher
# and its covariances must keep the structure's shared orientation.

skip_if_not_installed("MoEClust")
skip_if_not_installed("withr")

covariate_structure_data <- function() {
  withr::with_seed(41L, {
    n <- 500L
    z <- stats::rnorm(n)
    state <- 1L + stats::rbinom(n, 1L, stats::plogis(-0.3 + 1.1 * z))
    means <- rbind(c(0, 0, 0), c(2, 1.5, 1))
    first <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3L)
    second <- matrix(c(1.6, -0.3, 0.1, -0.3, 0.8, 0.2, 0.1, 0.2, 1.2), 3L)
    noise <- matrix(stats::rnorm(3L * n), n, 3L)
    rows <- t(vapply(seq_len(n), function(i) {
      means[state[i], ] + drop(noise[i, ] %*% chol(if (state[i] == 1L) first else second))
    }, numeric(3)))
    data.frame(y1 = rows[, 1L], y2 = rows[, 2L], y3 = rows[, 3L], z = z)
  })
}

structures <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE",
                "EEV", "VEV", "EVV", "VVV")

test_that("covariate models match MoEClust in every covariance structure", {
  data <- covariate_structure_data()
  invisible(lapply(structures, function(structure) {
    ours <- quietly(lpa(data, c("y1", "y2", "y3"), n_profiles = 2, model = structure,
                        profile_covariates = "z", n_starts = 20, seed = 1, tol = 1e-12,
                        max_iter = 20000))
    printed <- utils::capture.output(theirs <- withr::with_seed(1L, MoEClust::MoE_clust(
      data[c("y1", "y2", "y3")], G = 2, modelNames = structure, gating = ~ z,
      network.data = data,
      control = MoEClust::MoE_control(tol = c(1e-14, 1e-14, 1e-14),
                                      itmax = c(200000L, 200000L, 200000L),
                                      init.z = "random", nstarts = 20L, verbose = FALSE))))
    trace <- as.numeric(theirs$loglik)
    their_likelihood <- utils::tail(trace, 1L)
    expect_gte(ours$log_likelihood, their_likelihood - 1e-6)
    climbed <- length(trace) < 2L || utils::tail(diff(trace), 1L) >= -1e-10
    if (climbed) {
      gating <- stats::coef(theirs$gating)
      their_slope <- abs(if (is.matrix(gating)) gating[, "z"] else gating[["z"]])
      expect_equal(ours$log_likelihood, their_likelihood, tolerance = 1e-6,
                   label = paste(structure, "log likelihood"))
      expect_identical(as.integer(ours$n_parameters), as.integer(theirs$df),
                       label = paste(structure, "parameters"))
      expect_equal(abs(ours$profile_coefficients["z", 1L]), unname(their_slope),
                   tolerance = 1e-4, label = paste(structure, "slope"))
    } else {
      expect_gt(ours$log_likelihood, their_likelihood, label = structure)
      # A shared orientation: the second covariance is diagonal in the
      # first's eigenvectors.
      if (substr(structure, 3L, 3L) == "E") {
        axes <- eigen(ours$covariances[, , 1L], symmetric = TRUE)$vectors
        rotated <- crossprod(axes, ours$covariances[, , 2L] %*% axes)
        expect_lt(max(abs(rotated[upper.tri(rotated)])), 1e-8)
      }
    }
  }))
})
