# Two-stage bootstrap for fit_staged(): every resample refits the measurement
# and then the group-class structure, so the second-stage errors carry the
# first stage's uncertainty instead of treating the measurement as known.

skip_on_cran()

.staged_fixture <- function(n_groups = 50L, size = 8L, seed = 11L) {
  set.seed(seed)
  n <- n_groups * size
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size) + 1L
  profile <- 2L - stats::rbinom(n, 1L, c(0.8, 0.25)[group_class])
  data.frame(g = rep(seq_len(n_groups), each = size),
             y1 = stats::rnorm(n, c(0, 2)[profile]),
             y2 = stats::rnorm(n, c(0, 1.8)[profile]))
}

.staged_fit <- function(data) {
  quietly(fit_staged(data, c("y1", "y2"), "g", 2L, 2L, n_starts = 3, seed = 1))
}

.staged_bootstrap <- function(fit, iter = 40L, seed = 2L) {
  quietly(parameter_inference(fit, method = "bootstrap", iter = iter,
                              n_starts = 2L, seed = seed),
          classes = c(.multilpa_expected_warnings, "latents_bootstrap_dropped"))
}

test_that("a staged bootstrap reports every parameter, the held measurement included", {
  fit <- .staged_fit(.staged_fixture())
  wald <- parameter_inference(fit)
  bootstrap <- .staged_bootstrap(fit)
  expect_identical(attr(bootstrap, "stages_resampled"), 2L)
  expect_identical(nrow(bootstrap), length(.multilpa_coefficients(fit, "natural")))
  # The conditional Wald table has no rows for the held measurement.
  expect_gt(nrow(bootstrap), nrow(wald))
  measurement <- subset(bootstrap, level == "measurement")
  expect_setequal(unique(measurement$parameter), c("mean", "variance"))
  expect_true(all(is.finite(measurement$standard_error) &
                    measurement$standard_error > 0))
  expect_true(all(measurement$conf_low < measurement$estimate &
                    measurement$estimate < measurement$conf_high))
})

test_that("each resample re-estimates the measurement rather than holding it", {
  fit <- .staged_fit(.staged_fixture())
  bootstrap <- .staged_bootstrap(fit)
  replicates <- attr(bootstrap, "replicates")
  means <- grepl("mean", colnames(replicates))
  # A refit that held the first stage would give every replicate the same
  # means, exactly.
  expect_true(all(apply(replicates[, means, drop = FALSE], 2L, stats::sd) > 0.01))
})

test_that("propagating the first stage widens the second-stage errors", {
  fit <- .staged_fit(.staged_fixture(n_groups = 80L))
  wald <- parameter_inference(fit)
  bootstrap <- .staged_bootstrap(fit, iter = 150L)
  key <- function(table) paste(table$level, table$outcome, table$term)
  profile <- subset(wald, level == "profile")
  resampled <- bootstrap$standard_error[match(key(profile), key(bootstrap))]
  # A bootstrap that held the first stage reproduces these conditional errors
  # (measured: 0.0271 against 0.0279, 0.0524 against 0.0528); refitting the
  # measurement adds its variability, 14% and 32% on this fixture. The margin
  # is set between the two, so holding the first stage fails this test.
  expect_gt(mean(resampled / profile$standard_error), 1.1)
})

test_that("a staged bootstrap is reproducible from its seed", {
  fit <- .staged_fit(.staged_fixture())
  first <- .staged_bootstrap(fit, iter = 10L, seed = 5L)
  again <- .staged_bootstrap(fit, iter = 10L, seed = 5L)
  expect_identical(first$standard_error, again$standard_error)
})
