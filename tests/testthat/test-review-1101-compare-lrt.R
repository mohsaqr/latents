# Review 2026-11-01, compare_models(test = "lrt"): F1 shared predictors must
# keep their values, F3 ordinal category order is part of the coding, F6 the
# slope vocabulary is the fitter's, F4 a singular covariance restriction gets
# no chi-square reference.

skip_on_cran()

lrt_table <- function(...) as.data.frame(compare_models(..., test = "lrt"))

# Two Gaussian indicators, membership driven by z; w and v are unrelated.
membership_rows <- function(n = 300L, seed = 300L) {
  set.seed(seed)
  z <- stats::rnorm(n)
  state <- 1L + stats::rbinom(n, 1L, stats::plogis(2 * z))
  data.frame(a = stats::rnorm(n, c(-1.5, 1.5)[state]),
             b = stats::rnorm(n, c(-1, 1)[state]),
             z = z, w = stats::rnorm(n), v = stats::rnorm(n))
}
membership_fit <- function(data, covariates, ...) {
  quietly(lpa(data, c("a", "b"), 2, model = "EEI", profile_covariates = covariates,
              n_starts = 2, seed = 1, ...))
}

test_that("F1: a shared membership predictor with other values is not nested", {
  original <- membership_rows()
  shuffled <- original
  shuffled$z <- withr::with_seed(1, sample(original$z))
  small <- membership_fit(shuffled, "z")
  large <- membership_fit(original, c("z", "w"))
  expect_error(compare_models(small, large, test = "lrt"), class = "latents_bad_nesting")
})

test_that("F1: adding a predictor to unchanged or recoded shared columns is nested", {
  original <- membership_rows()
  large <- membership_fit(original, c("z", "w"))
  small <- membership_fit(original, "z")
  tested <- lrt_table(small, large)
  expect_identical(tested$lr_df[2L], 1L)
  expect_equal(tested$lr_statistic[2L],
               2 * (large$log_likelihood - small$log_likelihood), tolerance = 1e-12)
  # An affine recoding spans the same design: the null model, its likelihood
  # and so the statistic are unchanged.
  recoded <- original
  recoded$z <- 10 * original$z + 3
  rescaled <- membership_fit(recoded, "z", tol = 1e-10)
  refit <- membership_fit(original, "z", tol = 1e-10)
  expect_equal(rescaled$log_likelihood, refit$log_likelihood, tolerance = 1e-6)
  expect_identical(lrt_table(rescaled, large)$lr_df[2L], 1L)
})

test_that("F1: group-class membership predictors keep their values", {
  set.seed(42)
  groups <- 60L
  size <- 6L
  g <- rep(seq_len(groups), each = size)
  x <- stats::rnorm(groups)
  group_class <- 1L + stats::rbinom(groups, 1L, stats::plogis(1.5 * x))
  state <- 1L + stats::rbinom(groups * size, 1L, c(0.2, 0.8)[group_class[g]])
  data <- data.frame(g = g, x = x[g], u = stats::rnorm(groups)[g],
                     a = stats::rnorm(groups * size, c(-2, 2)[state]),
                     b = stats::rnorm(groups * size, c(-1.5, 1.5)[state]))
  shuffled <- data
  shuffled$x <- withr::with_seed(2, sample(x))[g]
  fit <- function(frame, covariates) {
    quietly(multilpa(frame, c("a", "b"), "g", 2, 2, group_covariates = covariates,
                     model = "EEI", n_starts = 2, seed = 1))
  }
  large <- fit(data, c("x", "u"))
  expect_identical(lrt_table(fit(data, "x"), large)$lr_df[2L], 1L)
  expect_error(compare_models(fit(shuffled, "x"), large, test = "lrt"),
               class = "latents_bad_nesting")
})

test_that("F1: a direct-effect predictor with other values is not nested", {
  original <- membership_rows()
  original$a <- original$a + 0.5 * original$w
  shuffled <- original
  shuffled$w <- withr::with_seed(3, sample(original$w))
  fit <- function(frame, effects) {
    quietly(lpa(frame, c("a", "b"), 2, model = "EEI", direct_effects = effects,
                n_starts = 2, seed = 1))
  }
  large <- fit(original, c("a ~ w", "b ~ v"))
  expect_identical(lrt_table(fit(original, "a ~ w"), large)$lr_df[2L], 1L)
  expect_error(compare_models(fit(shuffled, "a ~ w"), large, test = "lrt"),
               class = "latents_bad_nesting")
})

test_that("F1: a transition covariate with other values is not nested", {
  set.seed(5)
  persons <- 200L
  x <- stats::rnorm(persons)
  first <- stats::rbinom(persons, 1L, 0.5) + 1L
  second <- ifelse(stats::runif(persons) < stats::plogis(1 + x), first, 3L - first)
  occasion <- function(time, state) {
    data.frame(person = seq_len(persons), occasion = time,
               y = stats::rnorm(persons, c(-2, 2)[state]), x = x,
               v = withr::with_seed(9, stats::rnorm(persons)))
  }
  panel <- rbind(occasion(1L, first), occasion(2L, second))
  shuffled <- panel
  shuffled$x <- rep(withr::with_seed(4, sample(x)), 2L)
  fit <- function(frame, covariates) {
    quietly(lta(frame, "y", "person", n_profiles = 2, time = "occasion",
                transition_covariates = covariates, n_starts = 2, seed = 1))
  }
  large <- fit(panel, c("x", "v"))
  small <- fit(panel, "x")
  # One added covariate, one slope per off-diagonal transition logit.
  expect_identical(lrt_table(small, large)$lr_df[2L],
                   as.integer(large$n_parameters - small$n_parameters))
  expect_error(compare_models(fit(shuffled, "x"), large, test = "lrt"),
               class = "latents_bad_nesting")
})

ordinal_rows <- function(n = 300L) {
  set.seed(931)
  z <- stats::rnorm(n)
  state <- 1L + stats::rbinom(n, 1L, stats::plogis(3 * z))
  labels <- c("low", "middle", "high")
  draw <- function(k) {
    value <- exp(c(0, 0.1, -0.5) + 0:2 * c(-1, 1)[k])
    sample(labels, 1L, prob = value / sum(value))
  }
  u <- vapply(state, draw, character(1))
  t <- vapply(state, draw, character(1))
  data.frame(a = stats::rnorm(n, c(-1.5, 1.5)[state]), z = z,
             u = ordered(u, levels = labels), t = ordered(t, levels = labels))
}

test_that("F3: the order of ordinal categories is part of the coding", {
  data <- ordinal_rows()
  permuted <- data
  # The same labels and responses, the middle and last categories swapped:
  # a nonmonotonic reordering of the scores 0, 1, 2.
  permuted$u <- ordered(data$u, levels = c("low", "high", "middle"))
  fit <- function(frame, vars, ...) {
    quietly(lpa(frame, vars, 2, ordinal = intersect(vars, c("u", "t")),
                n_starts = 2, seed = 1, ...))
  }
  # Mixed Gaussian and ordinal indicators.
  plain <- fit(data, c("a", "u"))
  added <- fit(data, c("a", "u"), profile_covariates = "z")
  expect_identical(lrt_table(plain, added)$lr_df[2L], 1L)
  expect_error(compare_models(plain, fit(permuted, c("a", "u"), profile_covariates = "z"),
                              test = "lrt"),
               class = "latents_bad_nesting")
  # Ordinal indicators only.
  ordinal_plain <- fit(data, c("u", "t"))
  ordinal_added <- fit(data, c("u", "t"), profile_covariates = "z")
  expect_identical(lrt_table(ordinal_plain, ordinal_added)$lr_df[2L], 1L)
  expect_error(compare_models(ordinal_plain,
                              fit(permuted, c("u", "t"), profile_covariates = "z"),
                              test = "lrt"),
               class = "latents_bad_nesting")
})

test_that("F6: slopes shared across group classes nest in slopes by group class", {
  set.seed(917)
  g <- rep(1:60, each = 8)
  z <- stats::rnorm(length(g))
  group <- rep(rep(1:2, each = 30), each = 8)
  state <- 1L + stats::rbinom(length(g), 1L, stats::plogis(c(-1.8, 1.8)[group] + 0.5 * z))
  data <- data.frame(g, z, w = stats::rnorm(length(g)),
                     a = stats::rnorm(length(g), c(-2, 2)[state]),
                     b = stats::rnorm(length(g), c(-1.5, 1.5)[state]))
  fit <- function(slopes, covariates = "z") {
    quietly(multilpa(data, c("a", "b"), "g", 2, 2, profile_covariates = covariates,
                     profile_slopes = slopes, model = "EEI", n_starts = 3, seed = 1,
                     tol = 1e-9, max_iter = 5000))
  }
  shared <- fit("shared")
  separate <- fit("group_class")
  tested <- lrt_table(shared, separate)
  expect_identical(tested$lr_df[2L], as.integer(separate$n_parameters - shared$n_parameters))
  expect_true(is.finite(tested$lr_p_value[2L]))
  expect_error(compare_models(separate, shared, test = "lrt"), class = "latents_bad_nesting")
  # Two fits with slopes by group class, one restricted in its covariates.
  expect_identical(lrt_table(separate, fit("group_class", c("z", "w")))$lr_df[2L], 2L)
  # A value outside the vocabulary is refused by class, not a missing value in `if`.
  unknown <- shared
  unknown$profile_slopes <- "varying"
  expect_error(compare_models(unknown, separate, test = "lrt"), class = "latents_bad_nesting")
})

test_that("F4: exactly the spherical nulls within EVE, VVE, EEV and VEV are singular", {
  structures <- .multilpa_structures()
  stand_in <- function(structure, profiles = 2L, continuous = c("a", "b")) {
    structure(list(n_profiles = profiles, continuous = continuous,
                   covariance_structure = structure), class = "multilpa")
  }
  pairs <- expand.grid(null = structures, alternative = structures,
                       stringsAsFactors = FALSE)
  pairs$nested <- mapply(function(inner, outer) {
    all(vapply(seq_len(3L), function(position) {
      rank <- function(code) match(substr(code, position, position), c("I", "E", "V"))
      rank(inner) <= rank(outer)
    }, logical(1)))
  }, pairs$null, pairs$alternative)
  nested <- subset(pairs, nested & null != alternative)
  nested$singular <- mapply(function(inner, outer) {
    .latents_lrt_singular_covariance(stand_in(inner), stand_in(outer))
  }, nested$null, nested$alternative)
  expected <- c("EII:EVE", "EII:VVE", "EII:EEV", "EII:VEV", "VII:VVE", "VII:VEV")
  expect_setequal(with(subset(nested, singular), paste(null, alternative, sep = ":")),
                  expected)
  # One profile or one continuous indicator leaves no orientation to lose.
  expect_false(.latents_lrt_singular_covariance(stand_in("EII", 1L), stand_in("EVE", 1L)))
  expect_false(.latents_lrt_singular_covariance(stand_in("EII", continuous = "a"),
                                                stand_in("EVE", continuous = "a")))
})

test_that("F4: a singular covariance pair is refused, a regular pair still tested", {
  set.seed(91471)
  state <- rep(1:2, each = 150L)
  data <- data.frame(a = stats::rnorm(300L, c(-6, 6)[state]), b = stats::rnorm(300L))
  fit <- function(model) quietly(lpa(data, c("a", "b"), 2, model = model,
                                     n_starts = 1, seed = 1))
  spherical <- fit("EII")
  expect_error(compare_models(spherical, fit("EVE"), test = "lrt"),
               class = "latents_singular_restriction")
  expect_error(compare_models(spherical, fit("EEV"), test = "lrt"),
               class = "latents_bad_nesting")
  diagonal <- fit("EEI")
  varying <- fit("VVI")
  tested <- lrt_table(diagonal, varying)
  expect_identical(tested$lr_df[2L], as.integer(varying$n_parameters - diagonal$n_parameters))
  expect_equal(tested$lr_p_value[2L],
               stats::pchisq(tested$lr_statistic[2L], tested$lr_df[2L], lower.tail = FALSE),
               tolerance = 1e-12)
  # A non-spherical null fixes the orientation: EEI within EVE stays regular.
  expect_true(is.finite(lrt_table(diagonal, fit("EVE"))$lr_p_value[2L]))
  expect_true(is.finite(lrt_table(spherical, fit("EEE"))$lr_p_value[2L]))
})
