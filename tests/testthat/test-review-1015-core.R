# Regressions from the 2026-10-15 review of the core fitting files
# (multilpa.R, lpa.R, categorical.R, centering.R, prior.R). Each test failed
# before its repair.

skip_on_cran()

.review_core_data <- function(scale = 1) {
  skip_if_not_installed("withr")
  data <- withr::with_seed(71, data.frame(
    g = rep(seq_len(30), each = 6),
    a = rnorm(180) + rep(c(0, 2), length.out = 180),
    b = rnorm(180)))
  data$a <- data$a * scale
  data$b <- data$b * scale
  data
}

test_that("lpa() and lca() record the call the caller wrote, so update() refits", {
  items <- c("cognitive_strategies", "intrinsic_value")
  fit <- quietly(lpa(srl, items, n_profiles = 2, n_starts = 1, seed = 1))
  expect_identical(as.character(fit$call[[1L]]), "lpa")
  expect_identical(fit$call$data, as.name("srl"))
  # Before the repair the stored call read `multilpa(data = data, ...)`, which
  # update() evaluated against whatever `data` meant in the caller's frame.
  refit <- quietly(update(fit, n_profiles = 3))
  expect_identical(refit$n_profiles, 3L)
  expect_identical(refit$n_observations, nrow(srl))

  activities <- c("time_with_friends", "on_social_media")
  # Two classes on two binary items are not identified (Goodman 1974).
  classes <- quietly(lca(student_esm, activities, n_classes = 2, n_starts = 1,
                         seed = 1), classes = c(.multilpa_expected_warnings, "latents_unidentified"))
  expect_identical(as.character(classes$call[[1L]]), "lca")
  expect_identical(classes$call$n_classes, 2)
  # The direct two-level verb still reports its own call.
  direct <- quietly(multilca(subset(student_esm, day <= 1), activities,
                             id = "student", n_profiles = 2,
                             n_group_classes = 1, n_starts = 1, seed = 1),
                    classes = c(.multilpa_expected_warnings, "latents_unidentified"))
  expect_identical(as.character(direct$call[[1L]]), "multilca")
})

test_that("lca() refuses n_profiles with a classed error", {
  expect_error(lca(student_esm, c("time_with_friends", "on_social_media"),
                   n_classes = 2, n_profiles = 3),
               class = "latents_bad_argument")
})

test_that("a start that is not a list is refused as a bad start", {
  data <- .review_core_data()
  expect_error(multilpa(data, c("a", "b"), "g", n_profiles = 2,
                        n_group_classes = 1, start = 5, n_starts = 1),
               class = "latents_bad_start")
  expect_error(multilpa(data, c("a", "b"), "g", n_profiles = 2,
                        n_group_classes = 1, start = 5, fixed = "means",
                        n_starts = 1),
               class = "latents_bad_start")
  expect_error(multilpa(data, c("a", "b"), "g", n_profiles = 2,
                        n_group_classes = 1, start = "x", noise = TRUE,
                        n_starts = 1),
               class = "latents_bad_start")
})

test_that("a single-level fit takes only exactly one group class", {
  data <- .review_core_data()
  expect_error(multilpa(data, c("a", "b"), id = NULL, n_profiles = 2,
                        n_group_classes = 1.5, n_starts = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, c("a", "b"), id = NULL, n_profiles = 2,
                        n_group_classes = "1", n_starts = 1),
               class = "latents_bad_argument")
  fit <- quietly(multilpa(data, c("a", "b"), id = NULL, n_profiles = 2,
                          n_group_classes = 1, n_starts = 1, seed = 1))
  expect_identical(fit$n_group_classes, 1L)
})

test_that("prior_control() refuses a matrix scale that is not positive definite", {
  expect_error(prior_control(scale = matrix(c(1, 2, 2, 1), 2)),
               class = "latents_bad_argument")
  expect_error(prior_control(scale = matrix(1:6, 2)),
               class = "latents_bad_argument")
  expect_error(prior_control(scale = matrix(c(2, 1, 0, 2), 2)),
               class = "latents_bad_argument")
  expect_error(prior_control(scale = c(1, 2)), class = "latents_bad_argument")
  expect_s3_class(prior_control(scale = diag(2)), "latents_prior")
  expect_s3_class(prior_control(scale = 0.5), "latents_prior")
})

test_that("categorical values that print alike are refused, not merged in labels", {
  frame <- data.frame(a = rep(c(0.3, 0.1 + 0.2, 1), 20),
                      b = rep(c(1, 2, 2, 1, 3, 3), 10))
  expect_error(.multilpa_encode_categorical(frame),
               class = "latents_bad_categorical")
  # Distinct labels are unaffected.
  encoded <- .multilpa_encode_categorical(transform(frame, a = round(a, 6)))
  expect_identical(encoded$levels$a, c("0.3", "1"))
})

test_that("person centring is invariant to the order group indices appear in", {
  x <- matrix(c(1, 2, 10, 14, 100, 106), ncol = 1, dimnames = list(NULL, "a"))
  forward <- .multilpa_center_indicators(x, c(1L, 1L, 2L, 2L, 3L, 3L), 3L,
                                         "person")
  backward <- .multilpa_center_indicators(x, c(3L, 3L, 2L, 2L, 1L, 1L), 3L,
                                          "person")
  expect_equal(as.vector(forward$x), c(-0.5, 0.5, -2, 2, -3, 3))
  expect_equal(as.vector(backward$x), c(-0.5, 0.5, -2, 2, -3, 3))
  expect_equal(as.vector(backward$offsets), c(103, 12, 1.5))
})

test_that("centring does not refuse indicators for being recorded in small units", {
  for (centering in c("grand", "person")) {
    unit <- quietly(multilpa(.review_core_data(), c("a", "b"), "g",
                             n_profiles = 2, n_group_classes = 1,
                             centering = centering, n_starts = 1, seed = 1,
                             min_variance = 1e-30, tol = 1e-12))
    small <- quietly(multilpa(.review_core_data(1e-9), c("a", "b"), "g",
                              n_profiles = 2, n_group_classes = 1,
                              centering = centering, n_starts = 1, seed = 1,
                              min_variance = 1e-30, tol = 1e-12))
    # A change of units moves the log likelihood by the Jacobian only. One
    # deterministic start and a tight tolerance, because the relative
    # convergence rule itself reads |log L|, which the units shift.
    expect_equal(logLik(small)[[1L]],
                 logLik(unit)[[1L]] - 2 * 180 * log(1e-9),
                 tolerance = 1e-6)
  }
  # Constant within every unit is still refused, at any scale.
  flat <- .review_core_data(1e-9)
  flat$b <- rep(seq_len(30), each = 6) * 1e-9
  expect_error(multilpa(flat, c("a", "b"), "g", n_profiles = 2,
                        n_group_classes = 1, centering = "person",
                        n_starts = 1),
               class = "latents_bad_data")
})

test_that("a weighted fit flags small classes on its weighted effective counts", {
  skip_if_not_installed("withr")
  data <- withr::with_seed(2, data.frame(
    a = c(rnorm(57), rnorm(3, 8)), b = c(rnorm(57), rnorm(3, 8))))
  # Three outlying rows carry little weight: three sampled rows, but a
  # weighted effective membership far below one.
  data$w <- rep(c(1, 0.05), c(57, 3))
  warned <- character()
  fit <- withCallingHandlers(
    lpa(data, c("a", "b"), n_profiles = 2, weights = "w", n_starts = 3,
        seed = 1, variance_model = "equal"),
    warning = function(w) {
      warned <<- c(warned, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  expect_lt(min(fit$effective_profile_counts), 1)
  expect_true(fit$small_classes)
  expect_true("latents_small_classes" %in% warned)
  # The flag reads the counts the fit reports, weighted or not.
  expect_identical(fit$small_classes,
                   any(fit$effective_profile_counts < 1) ||
                     any(fit$effective_group_counts < 1))
})

test_that("capped VEI M-steps are counted and reported once per fit", {
  latents:::.multilpa_reset_structure_log()
  on.exit(latents:::.multilpa_reset_structure_log(), add = TRUE)
  diagonals <- rbind(c(4, 1, 1), c(1, 1, 9))
  # Before the repair each capped M-step warned on its own.
  expect_silent(.multilpa_structure_variances(diagonals, c(10, 20), "VEI",
                                              1e-10, max_iter = 1L))
  expect_silent(.multilpa_structure_variances(diagonals, c(10, 20), "VEI",
                                              1e-10, max_iter = 1L))
  expect_warning(latents:::.multilpa_report_structure_log(),
                 "2 M-step\\(s\\) of the VEI", class = "latents_no_converge")
})
