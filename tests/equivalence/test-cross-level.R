# Cross-level families (Houle et al., 2026, manifest specification). The
# reference likelihood is computed directly from the published formula,
#   L_j = sum_h omega_h f(m_j | h) prod_i sum_k pi_{k|h} f(x_ij | k),
# sharing no code with the package's E-step.

cross_level_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  type <- rep(1:2, each = 40L)
  profile_means <- rbind(c(-1.5, -1), c(1.5, 1))
  do.call(rbind, lapply(seq_along(type), function(j) {
    shares <- if (type[j] == 1L) c(0.8, 0.2) else c(0.25, 0.75)
    k <- sample(1:2, 6L, replace = TRUE, prob = shares)
    data.frame(team = sprintf("t%02d", j), type = type[j],
               y1 = stats::rnorm(6L, profile_means[k, 1L], 0.6),
               y2 = stats::rnorm(6L, profile_means[k, 2L], 0.6))
  }))
})
cross_level_vars <- c("y1", "y2")

# Direct likelihood of the manifest cross-level specification.
cross_level_reference <- function(fit, data) {
  groups <- split(data[, cross_level_vars], factor(data$team, levels = unique(data$team)))
  # Restricted: one prevalence vector for every class, the mean posterior.
  composition <- if (identical(fit$family, "restricted_cross_level")) {
    prevalence <- colMeans(fit$subject_posteriors)
    matrix(prevalence, nrow(fit$composition), length(prevalence), byrow = TRUE)
  } else fit$composition
  sum(vapply(groups, function(x) {
    x <- as.matrix(x)
    m <- colMeans(x)
    per_class <- vapply(seq_along(fit$group_probabilities), function(h) {
      individual <- vapply(seq_len(nrow(fit$profile_means)), function(k) {
        exp(rowSums(stats::dnorm(x, rep(fit$profile_means[k, ], each = nrow(x)),
                                 rep(sqrt(fit$profile_variances[k, ]), each = nrow(x)),
                                 log = TRUE)))
      }, numeric(nrow(x)))
      individual <- matrix(individual, nrow(x))
      log(fit$group_probabilities[h]) +
        sum(stats::dnorm(m, fit$group_means[h, ], sqrt(fit$group_mean_variances[h, ]),
                         log = TRUE)) +
        sum(log(individual %*% composition[h, ]))
    }, numeric(1))
    max(per_class) + log(sum(exp(per_class - max(per_class))))
  }, numeric(1)))
}

test_that("full cross-level likelihood equals the published formula", {
  fit <- multilpa(cross_level_data, cross_level_vars, "team", n_profiles = 2,
                  n_group_classes = 2, family = "full_cross_level", tol = 1e-10,
                  seed = 1)
  expect_true(fit$converged)
  expect_equal(fit$log_likelihood, cross_level_reference(fit, cross_level_data),
               tolerance = 1e-8)
  expect_true(all(diff(fit$log_likelihood_history) >=
                    -1e-10 * (1 + abs(fit$log_likelihood))))
  # 2*2 profile means + 2*2 variances + 2*(2-1) prevalences + 2*2 group means
  # + 2*2 group-mean variances + 1 class weight.
  expect_identical(fit$n_parameters, 19L)
})

test_that("restricted cross-level is the product of its two profile models", {
  fit <- multilpa(cross_level_data, cross_level_vars, "team", n_profiles = 2,
                  n_group_classes = 2, family = "restricted_cross_level", seed = 1)
  expect_equal(fit$log_likelihood, cross_level_reference(fit, cross_level_data),
               tolerance = 1e-6)
  individual <- withCallingHandlers(
    lpa(cross_level_data, cross_level_vars, n_profiles = 2, seed = 1),
    latents_single_level = \(m) invokeRestart("muffleMessage"))
  expect_equal(get_results(fit, "starts")$log_likelihood[1L],
               individual$log_likelihood, tolerance = 1e-6)
  expect_identical(fit$n_parameters, 18L)
})

test_that("the full cross-level EM reaches a stationary point of the formula", {
  skip_if_not_installed("numDeriv")
  fit <- multilpa(cross_level_data, cross_level_vars, "team", n_profiles = 2,
                  n_group_classes = 2, family = "full_cross_level", tol = 1e-12,
                  max_iter = 5000, seed = 1)
  # Perturb the group-mean block (class means, log variances) and the profile
  # means; every derivative of the reference likelihood must vanish.
  theta <- c(as.vector(fit$group_means), log(as.vector(fit$group_mean_variances)),
             as.vector(fit$profile_means))
  at <- function(v) {
    moved <- fit
    moved$group_means[] <- v[1:4]
    moved$group_mean_variances[] <- exp(v[5:8])
    moved$profile_means[] <- v[9:12]
    cross_level_reference(moved, cross_level_data)
  }
  expect_lt(max(abs(numDeriv::grad(at, theta))), 1e-3)
})
