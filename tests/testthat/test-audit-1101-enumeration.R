skip_on_cran()

# Audit 2026-11-01, slice `enumeration`: regression tests for the fixes in
# R/enumeration.R (bootstrap_lrt() covariate values, get_candidate()
# argument classes).

audit_covariate_data <- function(seed = 5L) {
  set.seed(seed)
  n_groups <- 40L
  size <- 6L
  n <- n_groups * size
  z <- stats::rnorm(n)
  class_of <- 1L + stats::rbinom(n, 1L, stats::plogis(1.2 * z))
  data.frame(g = rep(seq_len(n_groups), each = size), z = z,
             y1 = stats::rnorm(n, c(-1.5, 1.5)[class_of]),
             y2 = stats::rnorm(n, c(-1, 1)[class_of]))
}

audit_covariate_fit <- function(data, n_profiles, n_group_classes = 1L) {
  quietly(multilpa(data, c("y1", "y2"), "g", n_profiles = n_profiles,
                            n_group_classes = n_group_classes,
                            profile_covariates = "z", n_starts = 2L, seed = 1L))
}

test_that("bootstrap_lrt() refuses a shared covariate with other values as incomparable", {
  data <- audit_covariate_data()
  other <- data
  set.seed(9)
  other$z <- stats::rnorm(nrow(data))
  null_model <- audit_covariate_fit(data, 2L)
  alternative <- audit_covariate_fit(other, 3L)
  # Without `data` the alternative's stored columns are used; with either
  # frame supplied the pair is the same pair of different data sets.
  expect_error(bootstrap_lrt(null_model, alternative, iter = 2L, n_starts = 1L, seed = 1L),
               class = "latents_incomparable_models")
  expect_error(bootstrap_lrt(null_model, alternative, data = data, iter = 2L,
                             n_starts = 1L, seed = 1L),
               class = "latents_incomparable_models")
  expect_error(.multilpa_check_covariate_nesting(null_model, alternative),
               class = "latents_incomparable_models")
})

test_that("the covariate value check passes nested pairs at either level", {
  data <- audit_covariate_data()
  # A profile step and a group-class step (whose stacked designs differ in
  # rows) with the same covariate values are not refused.
  expect_null(.multilpa_check_covariate_nesting(audit_covariate_fit(data, 2L),
                                                audit_covariate_fit(data, 3L)))
  expect_null(.multilpa_check_covariate_nesting(audit_covariate_fit(data, 2L, 1L),
                                                audit_covariate_fit(data, 2L, 2L)))
})

test_that("get_candidate() refuses a fractional count or a vector model with latents_bad_argument", {
  data <- data.frame(g = rep(seq_len(10), each = 10),
                     y = c(stats::rnorm(50, -2), stats::rnorm(50, 2)))
  grid <- quietly(enumerate_classes(data, "y", "g", n_profiles = 1:2,
                                             n_group_classes = 1, model = "EEI",
                                             n_starts = 1, seed = 1))
  expect_error(get_candidate(grid, 2.5), class = "latents_bad_argument")
  expect_error(get_candidate(grid, 2, n_group_classes = 0), class = "latents_bad_argument")
  expect_error(get_candidate(grid, 2, model = c("EEI", "VVI")),
               class = "latents_bad_argument")
  expect_identical(get_candidate(grid, 2, model = "EEI")$n_profiles, 2L)
  expect_error(candidate_fit(grid, 1.5), class = "latents_bad_argument")
})

test_that("a family grid checks its count instead of reporting it truncated", {
  set.seed(8)
  school <- rep(seq_len(40), each = 6)
  shift <- rep(sample(1:2, 40, replace = TRUE), each = 6)
  family_data <- data.frame(school = school,
                            y1 = stats::rnorm(240, c(-1, 1)[shift]),
                            y2 = stats::rnorm(240, c(0, 1)[shift]))
  grid <- quietly(enumerate_classes(family_data, c("y1", "y2"), "school",
                                             family = "additive", n_group_classes = 1:2,
                                             between_variance = "equal",
                                             n_starts = 1, seed = 1))
  expect_error(get_candidate(grid, n_group_classes = 2.5, model = "additive_equal"),
               class = "latents_bad_argument")
  expect_error(get_candidate(grid, n_group_classes = 2, model = c("a", "b")),
               class = "latents_bad_argument")
  expect_s3_class(get_candidate(grid, n_group_classes = 2, model = "additive_equal"),
                  "multilpa_additive")
})
