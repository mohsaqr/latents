# Comparisons with independent references moved from tests/testthat/test-latents-controls.R,
# so that they run in CI's equivalence job rather than under R CMD check.

# A two-level latent class design with a student-level covariate, far from
# any boundary: two classes, five binary items, a slope of 0.8 on the logit.
.lca_design <- function(n_groups = 60, per = 8, continuous = FALSE) {
  set.seed(2024)
  g <- rep(seq_len(n_groups), each = per)
  cluster <- rep(sample(1:2, n_groups, TRUE), each = per)
  z <- rep(stats::rnorm(n_groups), each = per)
  class <- ifelse(stats::runif(length(g)) <
                    stats::plogis(ifelse(cluster == 1L, 1, -1) + 0.8 * z), 1L, 2L)
  draw <- function(p1, p2) {
    ifelse(stats::runif(length(g)) < ifelse(class == 1L, p1, p2), 1L, 2L)
  }
  data <- data.frame(g = g, z = z, y1 = draw(0.85, 0.2), y2 = draw(0.8, 0.25),
                     y3 = draw(0.75, 0.3), y4 = draw(0.2, 0.7),
                     y5 = draw(0.9, 0.35))
  if (continuous) data$w <- stats::rnorm(length(g), ifelse(class == 1L, 1, -1))
  data
}

.lca_items <- paste0("y", 1:5)

test_that("OPG standard errors reproduce glca on the same data", {
  skip_on_cran()
  # glca 1.4.2, item(y1, ..., y5) ~ z, group = g, nclass = 2, ncluster = 2,
  # on this seeded design (150 groups of 10): log likelihood -4486.5647 and
  # standard errors 0.16569, 0.16736 (the two intercepts) and 0.13636 (z).
  data <- .lca_design(n_groups = 150, per = 10)
  fit <- multilca(data, vars = .lca_items, id = "g", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "z", n_starts = 1,
                  seed = 1)
  expect_equal(fit$log_likelihood, -4486.5647, tolerance = 1e-7)
  opg <- subset(parameter_inference(fit, vcov_type = "opg"), level == "profile")
  expect_equal(sort(opg$standard_error), sort(c(0.16569, 0.16736, 0.13636)),
               tolerance = 2e-3)
  # The observed information is what the curvature says; glca does not use it.
  observed <- subset(parameter_inference(fit), level == "profile" & term == "z")
  expect_lt(observed$standard_error, subset(opg, term == "z")$standard_error)
})
