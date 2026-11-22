# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# The multilevel growth mixture model (phase E7 of validation/ENGINE_DESIGN.md):
# the growth block's person densities under the two-level structure, with
# persons nested in clusters and a group class per cluster. References: a
# dense likelihood that shares no code with the engine (full marginal
# covariance per person, a sum over group classes per cluster), numerical
# derivatives, the single-level model it reduces to, and recovery of a
# simulated truth.

multilevel_data <- function(seed = 1, n_schools = 16, per_school = 10, waves = 0:3) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  school_class <- rep(1:2, length.out = n_schools)
  improving <- c(0.85, 0.2)[school_class]
  schools <- lapply(seq_len(n_schools), function(s) {
    students <- lapply(seq_len(per_school), function(i) {
      k <- if (stats::runif(1) < improving[s]) 1L else 2L
      effects <- c(stats::rnorm(1, 0, 2), stats::rnorm(1, 0, 0.4))
      mean <- if (k == 1L) 50 + 3 * waves else 55 + 0 * waves
      data.frame(school = s, student = (s - 1) * per_school + i, wave = waves,
                 score = mean + effects[1L] + effects[2L] * waves +
                   stats::rnorm(length(waves), 0, 1.5),
                 trajectory = k, school_class = school_class[s])
    })
    do.call(rbind, students)
  })
  do.call(rbind, schools)
}

multilevel_fit <- function(data, ...) {
  mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                     class_level = "group", random = "wave",
                     random_covariance = "equal", n_group_classes = 2,
                     cluster = "school", n_starts = 3, seed = 1, ...)
}

# Dense likelihood: each person's marginal N(X beta_k, Z G Z' + sigma_k^2 I),
# mixed over classes within the cluster's group class, mixed over group
# classes per cluster. The fit has no group-membership covariate, so every
# cluster's group-class shares are the softmax of the intercept row of
# `delta`; `random_covariance = "equal"` gives one G for every class.
multilevel_dense <- function(fit, data) {
  params <- fit$params
  group_logits <- params$delta[1L, ]
  group_prior <- exp(group_logits - max(group_logits)) /
    sum(exp(group_logits - max(group_logits)))
  covariance_g <- params$random_covariance[[1L]]
  dense_log_density <- function(rows, k) {
    x <- cbind(1, rows$wave)
    covariance <- x %*% covariance_g %*% t(x) +
      diag(params$sigma2[k], nrow(rows))
    residual <- rows$score - as.vector(x %*% params$beta[, k])
    factor <- chol(covariance)
    -0.5 * (nrow(rows) * log(2 * pi) + 2 * sum(log(diag(factor))) +
              sum(backsolve(factor, residual, transpose = TRUE)^2))
  }
  schools <- split(data, data$school)
  sum(vapply(schools, function(school) {
    persons <- split(school, school$student)
    by_group_class <- vapply(seq_along(group_prior), function(h) {
      shares <- exp(params$class_logits[h, ]) / sum(exp(params$class_logits[h, ]))
      log(group_prior[h]) + sum(vapply(persons, function(rows) {
        log(sum(shares * exp(vapply(seq_len(fit$spec$n_classes), function(k) {
          dense_log_density(rows, k)
        }, numeric(1)))))
      }, numeric(1)))
    }, numeric(1))
    top <- max(by_group_class)
    top + log(sum(exp(by_group_class - top)))
  }, numeric(1)))
}

test_that("the analytic cluster scores are the likelihood's gradient", {
  skip_without_reference("numDeriv")
  data <- multilevel_data()
  fit <- multilevel_fit(data)
  theta <- latents:::.growth_pack(fit$spec, fit$params) + 0.05
  value <- function(v) {
    params <- latents:::.growth_unpack(fit$spec, v, fit$params)
    latents:::.growth_expectation(fit$spec, fit$stats, params)$log_likelihood
  }
  params <- latents:::.growth_unpack(fit$spec, theta, fit$params)
  scores <- latents:::.growth_scores(fit$spec, fit$stats, params)
  # One row per cluster: the clusters are the independent units.
  expect_identical(nrow(scores), fit$spec$n_clusters)
  # Scaled: numDeriv's error is absolute near zero, so a relative measure
  # would judge the rounding noise of a near-zero entry. 1e-6, the original.
  expect_reference(unname(colSums(scores)), numDeriv::grad(value, theta),
                   reference = "numDeriv", covers = "mixture_regression()",
                   quantity = "cluster score sum at a point off the optimum", tolerance = 1e-6,
                   measure = "scaled",
                   note = "Package side: the internal analytic scores (.growth_scores) behind vcov() of the multilevel growth mixture; reference: numDeriv::grad() of the package's likelihood.")
})


test_that("the likelihood equals the dense multilevel likelihood", {
  skip_on_cran()
  data <- multilevel_data()
  fit <- multilevel_fit(data)
  expect_s3_class(fit, "latents_growth_mixture")
  # The same likelihood summed in another order: 1e-10, the original tolerance.
  expect_reference(as.numeric(logLik(fit)), multilevel_dense(fit, data),
                   reference = "independent: dense multilevel growth mixture likelihood, base R",
                   covers = c("mixture_regression()", "logLik.latents_growth_mixture()"),
                   quantity = "log likelihood at the estimate", tolerance = 1e-10,
                   measure = "relative")
  # Off the maximum too.
  moved <- fit
  moved$params$class_logits[1L, 2L] <- moved$params$class_logits[1L, 2L] + 0.3
  moved$params$delta[1L, 2L] <- 0.4
  moved$params$beta[2L, 1L] <- moved$params$beta[2L, 1L] - 0.2
  expectation <- latents:::.growth_expectation(moved$spec, moved$stats, moved$params)
  expect_reference(expectation$log_likelihood, multilevel_dense(moved, data),
                   reference = "independent: dense multilevel growth mixture likelihood, base R",
                   covers = "mixture_regression()",
                   quantity = "log likelihood off the maximum", tolerance = 1e-10,
                   measure = "relative",
                   note = "Package side: the internal likelihood (.growth_expectation) that mixture_regression() maximises, at moved parameters.")
})
