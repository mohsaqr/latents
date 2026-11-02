# Audit 2026-11-01, growth-mixture slice: references for the marginal
# likelihood, predictions and simulation, the bootstrap table's identity, and
# regression tests for the plot fixes.

skip_on_cran()

audit_growth_data <- function(n = 60L, seed = 3L, n_classes = 1L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  times <- 0:4
  class <- sample(seq_len(n_classes), n, TRUE)
  out <- data.frame(id = rep(seq_len(n), each = length(times)),
                    time = rep(times, n))
  intercept <- c(2, 6, 10, 14, 18)[class]
  slope <- c(0.5, -0.4, 0.8, -0.8, 0)[class]
  out$y <- intercept[out$id] + slope[out$id] * out$time +
    rep(stats::rnorm(n), each = length(times)) +
    rep(stats::rnorm(n, 0, 0.4), each = length(times)) * out$time +
    stats::rnorm(nrow(out), 0, 0.7)
  out$true_class <- class[out$id]
  out[stats::runif(nrow(out)) > 0.1, ]
}

test_that("a one-class marginal likelihood equals nlme's ML fit and a hand integration", {
  skip_if_not_installed("nlme")
  data <- audit_growth_data()
  fit <- mixture_regression(y ~ time, data, n_classes = 1, id = "id",
                            class_level = "group", random = ~ 1 + time,
                            n_starts = 1, seed = 1, vcov_type = "none")
  reference <- nlme::lme(y ~ time, random = ~ 1 + time | id, data = data, method = "ML",
                         control = nlme::lmeControl(tolerance = 1e-10, msTol = 1e-10))
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)), tolerance = 1e-8)
  expect_equal(unname(fit$params$beta[, 1L]), unname(nlme::fixef(reference)),
               tolerance = 1e-5)
  # The marginal normal y_i ~ N(X beta, Z G Z' + sigma^2 I), integrated by hand.
  params <- fit$params
  by_person <- split(seq_len(nrow(data)), data$id)
  hand <- sum(vapply(by_person, function(rows) {
    design <- cbind(1, data$time[rows])
    covariance <- design %*% params$random_covariance[[1L]] %*% t(design) +
      diag(params$sigma2[1L], length(rows))
    residual <- data$y[rows] - as.vector(design %*% params$beta[, 1L])
    -0.5 * (length(rows) * log(2 * pi) +
              as.numeric(determinant(covariance)$modulus) +
              sum(residual * solve(covariance, residual)))
  }, numeric(1)))
  expect_equal(fit$log_likelihood, hand, tolerance = 1e-10)
})

test_that("predicted class trajectories are the design times the coefficients", {
  data <- audit_growth_data(n = 80L, n_classes = 2L)
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", random = ~ 1, n_starts = 1,
                            seed = 1, vcov_type = "none")
  predicted <- predict(fit, type = "class_response")
  design <- cbind(1, data$time)
  expected <- c(design %*% fit$params$beta[, 1L], design %*% fit$params$beta[, 2L])
  expect_equal(predicted$fitted, expected, tolerance = 1e-12)
  expect_identical(predicted$row, rep(seq_len(nrow(data)), 2L))
})

test_that("data simulated from a fit recover its parameters roughly", {
  data <- audit_growth_data(n = 300L, n_classes = 2L, seed = 9L)
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", random = ~ 1, n_starts = 1,
                            seed = 1, vcov_type = "none")
  simulated <- data
  simulated$y <- simulate(fit, nsim = 1, seed = 11)$sim_1
  refit <- mixture_regression(y ~ time, simulated, n_classes = 2, id = "id",
                              class_level = "group", random = ~ 1, n_starts = 1,
                              seed = 1, vcov_type = "none")
  expect_equal(refit$params$beta, fit$params$beta, tolerance = 0.1)
  expect_equal(refit$params$sigma2, fit$params$sigma2, tolerance = 0.15)
  expect_equal(unlist(refit$params$random_covariance),
               unlist(fit$params$random_covariance), tolerance = 0.3)
})

test_that("a growth bootstrap table carries the identity of its replicates", {
  data <- audit_growth_data(n = 80L, n_classes = 2L)
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", random = ~ 1, n_starts = 1,
                            seed = 1, vcov_type = "none")
  boot <- suppressWarnings(parameter_inference(fit, method = "bootstrap", iter = 3,
                                               n_starts = 1, seed = 2),
                           classes = "latents_bootstrap_dropped")
  identity <- attr(boot, "replicate_identity")
  expect_identical(names(identity), c("block", "class", "term", "estimate"))
  expect_identical(identity$estimate, boot$estimate)
  replicates <- as.data.frame(boot, what = "replicates")
  expect_identical(nrow(replicates), nrow(boot) * attr(boot, "n_valid"))
  expect_error(as.data.frame(boot[2:1, ], what = "replicates"),
               class = "latents_bad_argument")
})

test_that("the classification plot reports the weighted average posterior of the tables", {
  data <- audit_growth_data(n = 120L, n_classes = 2L, seed = 5L)
  person_weight <- round(stats::setNames(seq(0.2, 3, length.out = 120L), 1:120), 2)
  person_weight[1:10] <- 0
  data$w <- unname(person_weight[as.character(data$id)])
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", random = ~ 1, n_starts = 1,
                            seed = 1, weights = "w")
  classes <- as.data.frame(get_results(fit, "classes"))
  # Independently: the weighted mean of each carrying person's largest
  # posterior within their modal class.
  posterior <- fit$expectation$posterior
  modal <- max.col(posterior, ties.method = "first")
  carrying <- fit$spec$sampling_weights > 0
  expected <- vapply(1:2, function(k) {
    members <- modal == k & carrying
    sum(apply(posterior[members, , drop = FALSE], 1L, max) *
          fit$spec$sampling_weights[members]) / sum(fit$spec$sampling_weights[members])
  }, numeric(1))
  expect_equal(classes$average_posterior, expected, tolerance = 1e-12)
  drawn <- plot(fit, what = "classification")
  expect_match(drawn$labels$subtitle,
               sprintf("Class 1 %.2f, Class 2 %.2f", expected[1L], expected[2L]),
               fixed = TRUE)
})

test_that("the individuals view shows every class when persons must be cut", {
  data <- audit_growth_data(n = 100L, n_classes = 5L, seed = 7L)
  fit <- mixture_regression(y ~ time, data, n_classes = 5, id = "id",
                            class_level = "group", n_starts = 1, seed = 1,
                            vcov_type = "none")
  drawn <- plot(fit, what = "individuals", persons = 12)
  shown <- unique(drawn$data[c("person", "class")])
  expect_identical(nrow(shown), 12L)
  expect_setequal(shown$class, paste0("class_", 1:5))
})

test_that("growth plot arguments outside their range are refused with a class", {
  data <- audit_growth_data(n = 60L, n_classes = 1L)
  fit <- mixture_regression(y ~ time, data, n_classes = 1, id = "id",
                            class_level = "group", random = ~ 1, n_starts = 1,
                            seed = 1)
  expect_error(plot(fit, level = 1.5), class = "latents_bad_argument")
  expect_error(plot(fit, level = NA), class = "latents_bad_argument")
  expect_error(plot(fit, spread = NA), class = "latents_bad_argument")
  expect_error(plot(fit, facet = "yes"), class = "latents_bad_argument")
  expect_error(plot(fit, max_persons = NA), class = "latents_bad_argument")
  expect_error(plot(fit, what = "individuals", persons = -3),
               class = "latents_bad_argument")
  expect_error(plot(fit, what = "individuals", persons = NA),
               class = "latents_bad_argument")
  expect_s3_class(plot(fit, what = "individuals", persons = 3), "ggplot")
})
