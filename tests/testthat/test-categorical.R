categorical_fixture <- function(seed = 99L, n_groups = 40L, per_group = 15L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per_group)
  group_class <- rep(c(1L, 2L), length.out = n_groups)
  prevalence <- rbind(c(0.80, 0.20), c(0.25, 0.75))
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = prevalence[group_class[j], ])
  }, integer(1))
  probabilities <- rbind(c(0.85, 0.80, 0.75, 0.20, 0.15),
                         c(0.15, 0.20, 0.25, 0.80, 0.85))
  responses <- vapply(seq_len(5L), function(indicator) {
    1L + stats::rbinom(length(group), 1L, probabilities[profile, indicator])
  }, integer(length(group)))
  colnames(responses) <- paste0("v", 1:5)
  data.frame(school = group, responses)
}

test_that("categorical encoding is deterministic and type agnostic", {
  frame <- data.frame(
    f = factor(c("lo", "hi", "mid", "lo"), levels = c("lo", "mid", "hi")),
    ch = c("y", "n", "y", NA),
    b = c(0, 1, 1, 0),
    lg = c(TRUE, FALSE, TRUE, TRUE), stringsAsFactors = FALSE)
  encoded <- .multilpa_encode_categorical(frame)
  expect_identical(unname(encoded$n_categories), c(3L, 2L, 2L, 2L))
  # Factor levels keep their declared order, not their order of appearance.
  expect_identical(encoded$levels$f, c("lo", "mid", "hi"))
  expect_identical(encoded$codes[, "f"], c(1L, 3L, 2L, 1L))
  expect_true(is.na(encoded$codes[4L, "ch"]))
  # Codes recover the original values exactly.
  expect_identical(encoded$levels$f[encoded$codes[, "f"]], as.character(frame$f))
  # Row order cannot change the coding.
  shuffled <- .multilpa_encode_categorical(frame[c(3L, 1L, 4L, 2L), , drop = FALSE])
  expect_identical(shuffled$levels, encoded$levels)
})

test_that("categorical encoding rejects unusable indicators by condition class", {
  expect_error(.multilpa_encode_categorical(data.frame(a = c(1, 1, 1))),
               class = "multilpa_bad_categorical")
  expect_error(.multilpa_encode_categorical(data.frame(a = c(1, NA, NA))),
               class = "multilpa_bad_categorical")
  expect_error(.multilpa_encode_categorical(data.frame(a = c(1, Inf, 2))),
               class = "multilpa_bad_categorical")
})

test_that("the categorical M-step is the posterior-weighted response share", {
  codes <- matrix(c(1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L), 4L, 2L)
  posteriors <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), 4L, 2L)
  probabilities <- .multilpa_categorical_maximize(codes, posteriors, c(2L, 2L), 1e-10)
  # With hard assignments the estimate is the plain within-profile share.
  expect_equal(probabilities[[1L]][1L, ], c(0.5, 0.5))
  expect_equal(probabilities[[2L]][1L, ], c(1, 0), tolerance = 1e-9)
  expect_equal(probabilities[[2L]][2L, ], c(0, 1), tolerance = 1e-9)
  expect_true(all(vapply(probabilities, function(block) {
    max(abs(rowSums(block) - 1))
  }, numeric(1)) < 1e-12))
  # The probability bound is respected and is a genuine constraint.
  bounded <- .multilpa_categorical_maximize(codes, posteriors, c(2L, 2L), 0.01)
  expect_gte(min(bounded[[2L]]), 0.01 * (1 - 1e-9))
  expect_error(.multilpa_categorical_maximize(codes,
    matrix(c(1, 1, 0, 0, 0, 0, 0, 0), 4L, 2L), c(2L, 2L), 1e-10),
    class = "multilpa_empty_profile")
})

test_that("missing categorical cells contribute nothing to the log density", {
  codes <- matrix(c(1L, 2L, NA, 2L, 1L, 1L), 3L, 2L)
  probabilities <- list(rbind(c(0.7, 0.3), c(0.4, 0.6)),
                        rbind(c(0.8, 0.2), c(0.5, 0.5)))
  density <- .multilpa_categorical_log_density(codes, probabilities)
  expect_identical(dim(density), c(3L, 2L))
  # Row three observes only the second indicator, so its density is that term.
  expect_equal(density[3L, ], log(c(0.8, 0.5)))
  expect_equal(density[1L, ], log(c(0.7 * 0.2, 0.4 * 0.5)))
  expect_true(all(is.finite(density)))
})

test_that("thresholds invert to the cumulative response probabilities", {
  probabilities <- rbind(c(0.2, 0.3, 0.5), c(0.6, 0.1, 0.3))
  thresholds <- .multilpa_categorical_thresholds(probabilities)
  expect_identical(dim(thresholds), c(2L, 2L))
  expect_equal(stats::plogis(thresholds[, 1L]), probabilities[, 1L])
  expect_equal(stats::plogis(thresholds[, 2L]),
               probabilities[, 1L] + probabilities[, 2L])
  # A binary indicator yields one threshold, the logit of the first category.
  binary <- .multilpa_categorical_thresholds(rbind(c(0.25, 0.75)))
  expect_equal(as.vector(binary), stats::qlogis(0.25))
})

test_that("two-level latent class analysis recovers its generating structure", {
  dat <- categorical_fixture()
  vars <- paste0("v", 1:5)
  fit <- multilpa(dat, vars, "school", n_profiles = 2, n_group_classes = 2,
                    categorical = vars, n_starts = 20, seed = 3)
  expect_identical(fit$measurement_model, "categorical")
  expect_true(fit$converged)
  # (H - 1) + H (K - 1) + K sum(C - 1), with no Gaussian parameters.
  expect_equal(fit$n_parameters, 1 + 2 * 1 + 2 * 5)
  expect_identical(dim(fit$means), c(2L, 0L))
  expect_length(fit$response_probabilities, 5L)
  expect_identical(names(fit$response_probabilities), vars)
  expect_true(all(vapply(fit$response_probabilities, function(block) {
    max(abs(rowSums(block) - 1))
  }, numeric(1)) < 1e-10))
  # The two profiles must separate on every indicator, as generated.
  separation <- vapply(fit$response_probabilities, function(block) {
    abs(block[1L, 2L] - block[2L, 2L])
  }, numeric(1))
  expect_true(all(separation > 0.4))
  # Group classes must differ in profile prevalence, which is the multilevel part.
  expect_gt(max(abs(fit$profile_probabilities[1L, ] -
                      fit$profile_probabilities[2L, ])), 0.4)
  expect_equal(unname(rowSums(fit$profile_probabilities)), rep(1, 2))
})

test_that("categorical fits reproduce genuine Mplus two-level results", {
  reference <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-categorical.rds"))
  vars <- paste0("u", 1:5)
  fit <- multilpa(reference$data, vars, "clus", n_profiles = 2,
                    n_group_classes = 2, categorical = vars,
                    n_starts = 40, seed = 20260918, tol = 1e-13, max_iter = 20000)
  expect_true(fit$converged)
  expect_equal(fit$n_parameters, as.numeric(reference$n_parameters))
  thresholds <- vapply(fit$response_probabilities, function(block) {
    as.vector(.multilpa_categorical_thresholds(block))
  }, numeric(2))
  profile_order <- order(thresholds[, 1L])[rank(reference$mplus_thresholds[, 1L])]
  group_order <- order(fit$profile_probabilities[, profile_order[1L]])[
    rank(reference$mplus_profile_probabilities[, 1L])]
  expect_lt(max(abs(thresholds[profile_order, ] - reference$mplus_thresholds)), 1e-4)
  expect_lt(max(abs(fit$profile_probabilities[group_order, profile_order] -
                      reference$mplus_profile_probabilities)), 1e-5)
  expect_lt(max(abs(fit$group_probabilities[group_order] -
                      reference$mplus_group_probabilities)), 1e-5)
  expect_lt(abs(fit$log_likelihood - reference$mplus_log_likelihood), 1e-3)
})

test_that("ordinal and mixed-mode measurement fit and count parameters correctly", {
  set.seed(21)
  dat <- categorical_fixture(seed = 21L, n_groups = 30L, per_group = 12L)
  ordinal <- dat
  # Collapse two binary indicators into one three-category indicator.
  ordinal$v1 <- ordinal$v1 + ordinal$v2 - 1L
  ordinal$v2 <- NULL
  vars <- c("v1", "v3", "v4", "v5")
  fit <- multilpa(ordinal, vars, "school", 2, 2, categorical = vars,
                    n_starts = 12, seed = 4)
  expect_identical(ncol(fit$response_probabilities$v1), 3L)
  expect_equal(fit$n_parameters, 1 + 2 + 2 * (2 + 1 + 1 + 1))
  mixed <- dat
  mixed$score <- stats::rnorm(nrow(dat))
  mixed_fit <- multilpa(mixed, c("score", paste0("v", 1:5)), "school", 2, 2,
                          categorical = paste0("v", 1:5), n_starts = 12, seed = 4)
  expect_identical(mixed_fit$measurement_model, "mixed")
  expect_identical(mixed_fit$continuous, "score")
  expect_identical(mixed_fit$categorical, paste0("v", 1:5))
  # One Gaussian mean and variance per profile, plus the categorical block.
  expect_equal(mixed_fit$n_parameters, 1 + 2 + (2 + 2) + 2 * 5)
  expect_identical(dim(mixed_fit$means), c(2L, 1L))
  expect_identical(nrow(as.data.frame(mixed_fit)), 2L)
  expect_identical(nrow(as.data.frame(mixed_fit, what = "responses")), 20L)
})

test_that("categorical indicators support observed-data maximum likelihood", {
  dat <- categorical_fixture(seed = 7L, n_groups = 30L, per_group = 14L)
  vars <- paste0("v", 1:5)
  incomplete <- dat
  set.seed(12)
  incomplete$v1[sample(nrow(dat), 60L)] <- NA
  incomplete$v4[sample(nrow(dat), 40L)] <- NA
  expect_error(multilpa(incomplete, vars, "school", 2, 2,
                          categorical = vars, n_starts = 2, seed = 1),
               "missing or non-finite")
  fit <- multilpa(incomplete, vars, "school", 2, 2,
                    categorical = vars, missing = "fiml",
                    n_starts = 12, seed = 1)
  expect_true(fit$converged)
  expect_equal(unname(fit$n_observed_by_indicator[["v1"]]), nrow(dat) - 60)
  expect_equal(unname(fit$n_observed_by_indicator[["v2"]]), nrow(dat))
  # Every individual still receives posterior probabilities.
  expect_identical(nrow(fit$subject_posteriors), nrow(dat))
  expect_equal(unname(rowSums(fit$subject_posteriors)), rep(1, nrow(dat)))
})

test_that("the responses accessor reports probabilities and thresholds", {
  dat <- categorical_fixture(seed = 15L, n_groups = 25L, per_group = 12L)
  vars <- paste0("v", 1:5)
  fit <- multilpa(dat, vars, "school", 2, 2, categorical = vars,
                    n_starts = 10, seed = 6)
  responses <- as.data.frame(fit, what = "responses")
  expect_identical(names(responses),
    c("profile", "indicator", "category", "probability", "threshold"))
  expect_identical(nrow(responses), 2L * 5L * 2L)
  totals <- tapply(responses$probability,
                   list(responses$profile, responses$indicator), sum)
  expect_equal(as.vector(totals), rep(1, 10))
  # The final category has cumulative probability one, so no finite threshold.
  final <- subset(responses, responses$category == "2")
  expect_true(all(is.na(final$threshold)))
  first <- subset(responses, responses$category == "1")
  expect_equal(stats::plogis(first$threshold), first$probability)
  # A Gaussian-only fit returns an empty table of the same shape, not an error.
  gaussian <- multilpa(data.frame(g = rep(1:10, each = 8), y = stats::rnorm(80)),
                         "y", "g", 2, 1, n_starts = 3, seed = 1)
  expect_identical(nrow(as.data.frame(gaussian, what = "responses")), 0L)
})

test_that("unsupported categorical combinations are refused by condition class", {
  dat <- categorical_fixture(seed = 31L, n_groups = 25L, per_group = 12L)
  vars <- paste0("v", 1:5)
  fit <- multilpa(dat, vars, "school", 2, 2, categorical = vars,
                    n_starts = 8, seed = 2)
  smaller <- multilpa(dat, vars, "school", 1, 1, categorical = vars,
                        n_starts = 3, seed = 2)
  # Categorical measurement now carries analytic scores, so inference works.
  categorical_inference <- parameter_inference(fit, dat)
  expect_true(all(categorical_inference$standard_error >= 0))
  expect_true("response" %in% categorical_inference$parameter)
  # The parametric bootstrap now simulates categorical indicators too. Replicate
  # validity depends on the fit converging, which this small model does not
  # always manage, so only the machinery is asserted here.
  # The bootstrap compares models differing by one class at one level. One
  # profile with two group classes is unidentifiable, so the step is taken at
  # the group level instead: two profiles, one class against two.
  one_class <- multilpa(dat, vars, "school", 2, 1,
                        categorical = vars, n_starts = 3, seed = 2)
  bootstrap <- quietly(
    bootstrap_lrt(one_class, fit, dat, iter = 3, n_starts = 3, seed = 1))
  expect_gt(bootstrap$statistic, 0)
  expect_equal(bootstrap$iter, 3L)
  expect_true(all(c("replicate", "statistic", "valid") %in%
                    names(bootstrap$replicates)))
  expect_equal(nrow(bootstrap$replicates), 3L)
  expect_error(multilpa(dat, vars, "school", 2, 2, categorical = "absent",
                          n_starts = 2), class = "multilpa_bad_categorical")
  expect_error(multilpa(dat, vars, "school", 2, 2,
                          categorical = c("v1", "v1"), n_starts = 2),
               class = "multilpa_bad_categorical")
  # Categorical starts are supported; a start that names nothing the model uses
  # is rejected on its contents rather than refused outright.
  expect_error(multilpa(dat, vars, "school", 2, 2, categorical = vars,
                          n_starts = 1, start = list(a = 1)),
               "response_probabilities")
})

test_that("categorical models plot their response probabilities", {
  dat <- categorical_fixture(seed = 41L, n_groups = 25L, per_group = 12L)
  vars <- paste0("v", 1:5)
  fit <- multilpa(dat, vars, "school", 2, 2, categorical = vars,
                    n_starts = 8, seed = 2)
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 900, height = 600)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE, after = FALSE)
  expect_identical(plot(fit, what = "responses"), fit)
  expect_identical(plot(fit, what = "responses", category = "first"), fit)
  expect_identical(plot(fit, what = "probabilities"), fit)
  expect_error(plot(fit, what = "responses", category = "nope"),
               class = "multilpa_unknown_category")
  # A model with no continuous indicators cannot draw a profile-means plot.
  expect_error(plot(fit, what = "profiles"), class = "multilpa_no_continuous")
  gaussian <- multilpa(data.frame(g = rep(1:10, each = 8), y = stats::rnorm(80)),
                         "y", "g", 2, 1, n_starts = 3, seed = 1)
  expect_error(plot(gaussian, what = "responses"), class = "multilpa_no_categorical")
})
