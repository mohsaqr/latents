# Comparisons with independent references moved from tests/testthat/test-local-dependence.R,
# so that they run in CI's equivalence job rather than under R CMD check.

dependent_items <- function(seed = 5, n = 1500) {
  set.seed(seed)
  class <- sample(1:2, n, TRUE)
  # Two items driven by one latent response within each class, so they stay
  # associated given the class; two more that are independent given it.
  latent <- stats::rnorm(n)
  first <- ifelse(class == 1, latent > -0.5, latent > 0.8)
  second <- ifelse(class == 1, latent + stats::rnorm(n, 0, 0.6) > -0.5,
                   latent + stats::rnorm(n, 0, 0.6) > 0.8)
  yes_no <- function(value) factor(value, c(FALSE, TRUE), c("no", "yes"))
  data.frame(a = yes_no(first), b = yes_no(second),
             c = factor(stats::rbinom(n, 1, c(0.8, 0.2)[class]), 0:1, c("no", "yes")),
             e = factor(stats::rbinom(n, 1, c(0.75, 0.3)[class]), 0:1, c("no", "yes")))
}

joint_fit <- function(data, pairs = "a:b") {
  quietly(lca(data, c("a", "b", "c", "e"), 2, n_starts = 10, seed = 1, tol = 1e-12,
              local_dependence = pairs))
}

correlated_pair <- function(seed = 9, n = 1200) {
  set.seed(seed)
  class <- sample(1:2, n, TRUE)
  correlation <- c(0.6, -0.5)[class]
  first <- stats::rnorm(n)
  second <- correlation * first + sqrt(1 - correlation^2) * stats::rnorm(n)
  data.frame(y1 = c(0, 2.5)[class] + first, y2 = c(0, 2)[class] + second,
             y3 = stats::rnorm(n, c(0, 1.5)[class]), y4 = stats::rnorm(n, c(0, 1)[class]))
}

test_that("a jointly modelled pair agrees with Latent GOLD's class-specific association", {
  skip_on_cran()
  # Latent GOLD 6.1 on these data, `a <-> b | Cluster`: log likelihood
  # -3614.1388 with 11 parameters, class sizes 0.5093 and 0.4907, and these
  # marginal probabilities of `a` and `b` being "yes" in each class.
  fit <- joint_fit(dependent_items())
  within <- function(actual, printed) {
    expect_true(all(abs(actual - printed) <= 5e-5),
                info = paste(format(actual, digits = 8), collapse = " "))
  }
  within(fit$log_likelihood, -3614.1388)
  expect_equal(fit$n_parameters, 11)
  sizes <- colMeans(fit$subject_posteriors)
  order <- order(-sizes)
  within(sizes[order], c(0.5093, 0.4907))
  table <- get_results(fit, "local_dependence")
  yes_a <- aggregate(probability ~ profile, subset(table, category_1 == "yes"), sum)
  yes_b <- aggregate(probability ~ profile, subset(table, category_2 == "yes"), sum)
  within(yes_a$probability[order], c(0.6914, 0.2184))
  within(yes_b$probability[order], c(0.6827, 0.2750))
})

test_that("a covariance for a continuous pair agrees with Latent GOLD", {
  skip_on_cran()
  # Latent GOLD 6.1, class-specific variances and `y1 <-> y2 | Cluster`:
  # log likelihood -7272.5017 with 19 parameters, the pair's covariances
  # 0.5968 (s.e. 0.0580) and -0.4362 (s.e. 0.0442).
  fit <- quietly(lpa(correlated_pair(), paste0("y", 1:4), 2, n_starts = 10, seed = 1,
                     tol = 1e-12, local_dependence = "y1:y2"))
  within <- function(actual, printed) {
    expect_true(all(abs(actual - printed) <= 5e-5),
                info = paste(format(actual, digits = 8), collapse = " "))
  }
  within(fit$log_likelihood, -7272.5017)
  expect_equal(fit$n_parameters, 19)
  table <- get_results(fit, "local_dependence")
  within(sort(table$covariance), c(-0.4362, 0.5968))
  inference <- parameter_inference(fit)
  pair_rows <- subset(inference, parameter == "covariance" & term == "y2:y1")
  within(sort(pair_rows$standard_error), c(0.0442, 0.0580))
  # The other covariances are zero and carry no coordinate.
  expect_true(all(fit$covariances[3, 1, ] == 0) && all(fit$covariances[4, 2, ] == 0))
  expect_false(any(subset(inference, parameter == "covariance")$term %in%
                     c("y3:y1", "y4:y2")))
})
