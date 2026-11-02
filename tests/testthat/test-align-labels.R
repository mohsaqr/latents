# align_labels(): renumbering one fit's profiles and group classes to match
# another's. The exact check refits the same model from the reference's own
# estimates listed in another order: EM then reaches the same maximum under
# permuted labels, independently of the relabelling code, and the aligned
# refit must equal the refit from the unpermuted start in every table.

skip_on_cran()

permuted_start <- function(fit, profile_order, group_order = seq_len(fit$n_group_classes)) {
  start <- unclass(starting_values(fit))
  if (ncol(start$means) > 0L) start$means <- start$means[profile_order, , drop = FALSE]
  if (!is.null(start$covariances)) {
    start$covariances <- start$covariances[, , profile_order, drop = FALSE]
  }
  if (!is.null(start$variances)) start$variances <- start$variances[profile_order, , drop = FALSE]
  if (!is.null(start$response_probabilities)) {
    start$response_probabilities <- lapply(start$response_probabilities, function(block) {
      block[profile_order, , drop = FALSE]
    })
  }
  start$profile_probabilities <- start$profile_probabilities[group_order, profile_order,
                                                             drop = FALSE]
  start$group_probabilities <- start$group_probabilities[group_order]
  start
}

expect_same_tables <- function(aligned, control) {
  tables <- names(.multilpa_catalogue(control))
  invisible(lapply(tables, function(table) {
    # A table that does not apply to the model is refused; both fits must
    # refuse it with the same condition.
    read <- function(x) tryCatch(get_results(x, table), error = function(e) class(e)[1L])
    expect_equal(read(aligned), read(control), tolerance = 1e-6, ignore_attr = TRUE,
                 label = table)
  }))
}

expect_aligned_refit <- function(refit, reference, profile_order, group_order = NULL) {
  groups <- seq_len(reference$n_group_classes)
  control <- refit(permuted_start(reference, seq_len(reference$n_profiles), groups))
  swapped <- refit(permuted_start(reference, profile_order, group_order %||% groups))
  expect_equal(swapped$log_likelihood, control$log_likelihood, tolerance = 1e-10)
  expect_false(identical(swapped$subject_profiles, control$subject_profiles))
  aligned <- align_labels(swapped, control)
  expect_identical(aligned$subject_profiles, control$subject_profiles)
  expect_identical(aligned$group_classes, control$group_classes)
  expect_same_tables(aligned, control)
}

test_that("a refit with permuted labels is aligned back in every table", {
  vars <- c("eruptions", "waiting")
  reference <- lpa(datasets::faithful, vars, n_profiles = 3, model = "VVV", n_starts = 3,
                   seed = 1)
  refit <- function(start) {
    lpa(datasets::faithful, vars, n_profiles = 3, model = "VVV", n_starts = 1, start = start)
  }
  expect_aligned_refit(refit, reference, c(3L, 1L, 2L))
})

test_that("latent classes and the group classes of a two-level model are aligned", {
  items <- c("time_with_friends", "on_social_media", "sports", "walking")
  first_week <- subset(student_esm, day <= 3)
  # A tight tolerance leaves no score large enough to warn under Wald inference.
  classes <- lca(first_week, items, n_classes = 3, n_starts = 3, seed = 1, tol = 1e-10)
  expect_aligned_refit(function(start) {
    lca(first_week, items, n_classes = 3, n_starts = 1, start = start, tol = 1e-10)
  }, classes, c(2L, 3L, 1L))
  feelings <- c("happy", "relaxed", "worried")
  two_level <- function(...) {
    quietly(multilpa(student_esm, feelings, "student", n_profiles = 3, n_group_classes = 2,
                     ...), "latents_boundary")
  }
  expect_aligned_refit(function(start) two_level(n_starts = 1, start = start),
                       two_level(n_starts = 3, seed = 1), c(3L, 1L, 2L), c(2L, 1L))
})

test_that("aligning changes labels only", {
  first_days <- subset(student_esm, day <= 3)
  last_days <- subset(student_esm, day > 3)
  feelings <- c("happy", "exhausted")
  early <- lpa(first_days, feelings, n_profiles = 2, n_starts = 2, seed = 1)
  late <- lpa(last_days, feelings, n_profiles = 2, n_starts = 2, seed = 1)
  aligned <- align_labels(late, early)
  # Each aligned profile is the late profile nearest the early profile of
  # the same label.
  nearest <- vapply(seq_len(2L), function(k) {
    which.min(colSums((t(late$means) - early$means[k, ])^2))
  }, integer(1))
  expect_identical(nearest, 2:1)
  expect_equal(aligned$means, late$means[nearest, ], ignore_attr = TRUE)
  expect_identical(aligned$log_likelihood, late$log_likelihood)
  expect_identical(table(aligned$subject_profiles), table(late$subject_profiles)[nearest],
                   ignore_attr = TRUE)
  expect_equal(rowSums(aligned$subject_posteriors), rep(1, nrow(late$subject_posteriors)),
               ignore_attr = TRUE)
  # A fit aligned to itself is unchanged.
  expect_identical(align_labels(early, early), early)
})

test_that("a covariate fit and a full-covariance fit align to simpler models", {
  # Three separated profiles whose shares move with z.
  set.seed(3)
  z <- stats::rnorm(450)
  membership <- vapply(z, function(value) {
    sample.int(3L, 1L, prob = exp(c(0.8, 0, -0.8) * value))
  }, integer(1))
  centres <- rbind(c(-3, 0), c(0, 3), c(3, -3))
  data <- data.frame(y1 = centres[membership, 1L] + stats::rnorm(450),
                     y2 = centres[membership, 2L] + stats::rnorm(450), z = z)
  fit <- function(...) lpa(data, c("y1", "y2"), n_profiles = 3, n_starts = 2, seed = 1, ...)
  plain <- fit(model = "VVI")
  covariate <- fit(model = "VVI", profile_covariates = "z")
  scrambled <- .multilpa_cov_permute(covariate, c(3L, 1L, 2L), 1L)
  aligned <- align_labels(scrambled, plain)
  expect_s3_class(aligned, "multilpa_covariates")
  expect_equal(aligned$means, covariate$means, tolerance = 1e-12)
  # The membership logits are re-expressed against the restored last profile.
  expect_equal(aligned$profile_coefficients, covariate$profile_coefficients,
               tolerance = 1e-10)
  expect_equal(get_results(aligned, "posteriors"), get_results(covariate, "posteriors"),
               tolerance = 1e-10)
  full <- align_labels(fit(model = "VVV"), plain)
  expect_s3_class(full, "multilpa")
  expect_identical(dim(full$covariances), c(2L, 2L, 3L))
})

test_that("models that cannot be matched are refused", {
  vars <- c("eruptions", "waiting")
  two <- lpa(datasets::faithful, vars, n_profiles = 2, n_starts = 1, seed = 1)
  three <- lpa(datasets::faithful, vars, n_profiles = 3, n_starts = 1, seed = 1)
  one_indicator <- lpa(datasets::faithful, "waiting", n_profiles = 2, n_starts = 1, seed = 1)
  expect_error(align_labels(three, two), class = "latents_bad_argument")
  expect_error(align_labels(one_indicator, two), class = "latents_bad_argument")
  expect_error(align_labels(two, as.data.frame(two)), class = "latents_bad_argument")
  transitions <- lta(course_engagement, c("browse", "lectures"), "student", n_profiles = 2,
                     time = "sequence", n_starts = 1, seed = 1)
  expect_error(align_labels(transitions, transitions), class = "latents_bad_argument")
})
