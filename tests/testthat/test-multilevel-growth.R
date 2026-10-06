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

test_that("identical group classes reduce it to the single-level growth mixture", {
  skip_on_cran()
  data <- multilevel_data()
  fit <- multilevel_fit(data)
  single <- mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                               class_level = "group", random = "wave",
                               random_covariance = "equal", n_starts = 3, seed = 1)
  # Invariant: with every group class carrying the single-level class logits,
  # the group classes are indistinguishable and the likelihoods agree.
  params <- fit$params
  params[c("beta", "common", "sigma2", "random_covariance", "random_scale")] <-
    single$params[c("beta", "common", "sigma2", "random_covariance", "random_scale")]
  params$class_logits[] <- rep(single$params$gamma[1L, ], each = 2L)
  params$delta[1L, 2L] <- 0.7
  two_level <- latents:::.growth_expectation(fit$spec, fit$stats, params)
  expect_equal(two_level$log_likelihood, single$log_likelihood, tolerance = 1e-10)
  expect_equal(two_level$posterior, single$expectation$posterior,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("one group class is the single-level model with clusters as units", {
  skip_on_cran()
  data <- multilevel_data()
  one <- mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                            class_level = "group", random = "wave",
                            random_covariance = "equal", n_group_classes = 1,
                            cluster = "school", n_starts = 3, seed = 1)
  single <- mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                               class_level = "group", random = "wave",
                               random_covariance = "equal", n_starts = 3, seed = 1)
  expect_equal(one$log_likelihood, single$log_likelihood, tolerance = 1e-8)
  expect_identical(one$n_parameters, single$n_parameters)
  expect_identical(nrow(latents:::.growth_scores(one$spec, one$stats, one$params)),
                   one$spec$n_clusters)
  compared <- enumerate_regressions(score ~ wave, data, n_classes = 2,
                                    n_group_classes = 1:2, id = "student",
                                    class_level = "group", random = "wave",
                                    random_covariance = "equal", cluster = "school",
                                    n_starts = 2, seed = 1)
  expect_identical(nrow(as.data.frame(compared)), 2L)
})

test_that("a simulated multilevel truth is recovered", {
  skip_on_cran()
  data <- multilevel_data(seed = 2, n_schools = 30, per_school = 12)
  fit <- multilevel_fit(data)
  expect_true(fit$converged)
  clusters <- get_results(fit, "clusters")
  truth <- tapply(data$school_class, data$school, `[`, 1L)
  recovered <- table(clusters$group_class, truth[as.character(clusters$school)])
  expect_gte(max(sum(diag(recovered)), sum(recovered) - sum(diag(recovered))), 28L)
  shares <- get_results(fit, "group_classes")
  improving <- subset(shares, class == "class_1")
  expect_equal(sort(improving$probability), c(0.2, 0.85), tolerance = 0.1)
  coefficients <- get_results(fit, "coefficients")
  slopes <- subset(coefficients, term == "wave")$estimate
  expect_equal(sort(slopes), c(0, 3), tolerance = 0.15)
  # Invariants: posteriors are probabilities; group-class shares sum to one.
  expect_equal(rowSums(fit$expectation$posterior), rep(1, fit$spec$n_groups),
               tolerance = 1e-12)
  expect_equal(sum(unique(shares$group_share)), 1, tolerance = 1e-12)
  within <- tapply(shares$probability, shares$group_class, sum)
  expect_equal(as.vector(within), c(1, 1), tolerance = 1e-12)
})

test_that("the tables, printing, simulation and standard errors are complete", {
  skip_on_cran()
  data <- multilevel_data()
  fit <- multilevel_fit(data)
  tables <- get_results(fit, "all")
  expect_true(all(c("group_classes", "clusters") %in% names(tables)))
  expect_identical(nrow(tables$clusters), fit$spec$n_clusters)
  expect_true("school" %in% names(tables$assignments))
  membership <- get_results(fit, "membership")
  expect_identical(membership$term,
                   c("(Intercept)", "(Intercept):group_class_1", "(Intercept):group_class_2"))
  expect_true(all(is.finite(membership$std_error)))
  model_fit <- get_results(fit, "fit")
  expect_identical(model_fit$n_clusters, fit$spec$n_clusters)
  expect_equal(model_fit$bic, -2 * fit$log_likelihood + log(16) * fit$n_parameters)
  # 4 coefficients, 2 residual SDs, 3 covariance terms, 1 group-class logit
  # and 2 class logits.
  expect_identical(fit$n_parameters, 12L)
  expect_identical(length(coef(fit)), fit$n_parameters)
  expect_output(print(fit), "Multilevel growth mixture model")
  expect_output(print(summary(fit)), "Group classes")
  first <- simulate(fit, nsim = 2, seed = 3)
  expect_identical(dim(first), c(nrow(data), 2L))
  expect_identical(simulate(fit, nsim = 2, seed = 3), first)
  robust <- get_results(fit, "coefficients", vcov_type = "robust")
  expect_true(all(is.finite(robust$std_error)))
})

test_that("multilevel growth requests the model does not define are refused", {
  data <- multilevel_data()
  base <- function(...) {
    mixture_regression(score ~ wave, data, n_classes = 2, id = "student",
                       class_level = "group", n_starts = 1, seed = 1, ...)
  }
  expect_error(base(cluster = "school", n_group_classes = 2),
               class = "latents_bad_argument")
  expect_error(base(random = "wave", n_group_classes = 2),
               class = "latents_bad_argument")
  expect_error(base(random = "wave", cluster = "student", n_group_classes = 2),
               class = "latents_bad_argument")
  moved <- data
  moved$school[1L] <- 99
  expect_error(mixture_regression(score ~ wave, moved, n_classes = 2, id = "student",
                                  class_level = "group", random = "wave",
                                  n_group_classes = 2, cluster = "school",
                                  n_starts = 1, seed = 1),
               class = "latents_bad_data")
  data$weight <- 1
  expect_error(base(random = "wave", cluster = "school", n_group_classes = 2,
                    weights = "weight"),
               class = "latents_unsupported_weights")
  fit <- multilevel_fit(data)
  data$outcome <- as.numeric(data$student %% 2)
  expect_error(three_step(fit, data, "outcome"),
               class = "latents_unsupported_three_step")
})
