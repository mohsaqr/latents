skip_on_cran()

test_that("zero-weight sentinel groups leave additive estimates unchanged", {
  withr::local_seed(271)
  data <- data.frame(
    id = rep(seq_len(30L), each = 4L),
    x = rep(stats::rnorm(30L), each = 4L) + stats::rnorm(120L),
    weight = 1)
  fit <- function(data) multilpa(
    data, "x", "id", family = "additive", n_group_classes = 2,
    weights = "weight", n_starts = 2, seed = 3)
  # This fixture estimates one between-group variance at zero; only that
  # known, classed warning is expected. Other warnings must remain visible.
  expect_warning(reference <- fit(data), class = "latents_boundary")
  excluded <- data.frame(id = rep(31:32, each = 4L),
                          x = 1e6 + rep(c(-1, 1), 4L), weight = 0)
  expect_warning(augmented <- fit(rbind(data, excluded)), class = "latents_boundary")

  expect_true(reference$converged)
  expect_identical(augmented$converged, reference$converged)
  expect_equal(augmented$log_likelihood, reference$log_likelihood,
               tolerance = 1e-10)
  expect_equal(augmented$means, reference$means, tolerance = 1e-8)
  expect_equal(augmented$within_variances, reference$within_variances,
               tolerance = 1e-8)
  expect_equal(augmented$between_variances, reference$between_variances,
               tolerance = 1e-8)
  expect_equal(augmented$group_probabilities, reference$group_probabilities,
               tolerance = 1e-8)
  expect_equal(unname(augmented$group_posteriors[seq_len(30L), ]),
               unname(reference$group_posteriors), tolerance = 1e-8)
  expect_equal(rowSums(augmented$group_posteriors),
               rep(1, nrow(augmented$group_posteriors)),
               ignore_attr = TRUE, tolerance = 1e-14)
})

test_that("one profile plus noise retains its classification checks", {
  withr::local_seed(821)
  data <- data.frame(
    x = c(stats::rnorm(250L), stats::runif(50L, -15, 15)),
    y = c(stats::rnorm(250L), stats::runif(50L, -15, 15)))
  fit <- lpa(data, c("x", "y"), n_profiles = 1,
              noise = TRUE, n_starts = 2, seed = 1)
  quality <- diagnostics(fit, plots = FALSE)
  checks <- get_results(quality, "checks")
  entropy <- get_results(quality, "entropy")
  classification <- get_results(quality, "classification")
  expect_equal(entropy$n_classes[entropy$level == "individuals"], 2)
  expect_setequal(classification$class, c(0L, 1L))
  expect_true(all(classification$n_modal > 0L))
  value <- function(name) checks$value[checks$check == name &
                                       checks$level == "individuals"]
  expect_equal(value("relative entropy"),
               entropy$relative_entropy[entropy$level == "individuals"])
  expect_equal(value("lowest average posterior"),
               min(classification$average_posterior))
  expect_equal(value("smallest class"), min(classification$estimated_proportion))
  expect_false(any(checks$status[checks$domain == "Classification"] == "n/a"))
  expect_identical(checks$status[checks$domain == "Local independence"], "n/a")

  ordinary <- lpa(data, c("x", "y"), n_profiles = 1,
                   n_starts = 1, seed = 1)
  checks <- get_results(diagnostics(ordinary, plots = FALSE), "checks")
  expect_true(all(checks$status[checks$domain == "Classification"] == "n/a"))
})
