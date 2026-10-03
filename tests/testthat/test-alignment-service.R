# The alignment and simulation services: matching labels, relabelling exactly,
# seed scoping, categorical draws and the bootstrap replicate loops.

alignment_profile_fit <- function() {
  quietly(multilpa(engagement_small, c("browse", "lectures"), "student", 3, 2,
                   n_starts = 2, seed = 1))
}

alignment_likelihood <- function(model) {
  x <- .multilpa_center_like(model, as.matrix(engagement_small[c("browse", "lectures")]))
  .multilpa_expectation(x, model$group_index, model)$log_likelihood
}

alignment_covariate_rows <- function() {
  set.seed(5)
  n_groups <- 40L
  size <- 8L
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size), z = stats::rnorm(n))
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z))
  frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  frame$y2 <- stats::rnorm(n, c(0, 2)[state])
  frame
}

test_that("matching a relabelled profile fit recovers the inverse permutation", {
  fit <- alignment_profile_fit()
  scramble <- c(3L, 1L, 2L)
  scrambled <- .multilpa_permute_group_classes(
    .multilpa_permute_profiles(fit, scramble), c(2L, 1L))
  # Relabelling is exact: the likelihood does not move.
  expect_equal(alignment_likelihood(scrambled), fit$log_likelihood, tolerance = 1e-10)
  alignment <- .multilpa_pool_align(scrambled, fit)
  expect_identical(alignment$profile_order, order(scramble))
  expect_identical(alignment$group_order, c(2L, 1L))
  expect_identical(alignment$fit$means, fit$means)
  expect_identical(alignment$fit$variances, fit$variances)
  expect_identical(alignment$fit$profile_probabilities, fit$profile_probabilities)
  expect_identical(alignment$fit$subject_profiles, fit$subject_profiles)
  expect_identical(alignment$fit$group_classes, fit$group_classes)
  expect_equal(alignment_likelihood(alignment$fit), fit$log_likelihood, tolerance = 1e-10)
  # The replicate-only form returns the same relabelled fit.
  expect_identical(.multilpa_align_labels(scrambled, fit), alignment$fit)
})

test_that("matching a relabelled covariate fit restores its logits and likelihood", {
  data <- alignment_covariate_rows()
  fit <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L, profile_covariates = "z",
                          n_starts = 2, seed = 1),
                 classes = c(.multilpa_expected_warnings, "latents_extreme_coefficients"))
  scrambled <- .multilpa_cov_permute(fit, c(2L, 1L), 1L)
  expect_false(isTRUE(all.equal(scrambled$profile_coefficients, fit$profile_coefficients)))
  expect_silent(.multilpa_cov_check_likelihood(scrambled, data))
  alignment <- .multilpa_pool_align(scrambled, fit)
  expect_identical(alignment$profile_order, c(2L, 1L))
  expect_equal(alignment$fit$profile_coefficients, fit$profile_coefficients,
               tolerance = 1e-12)
  expect_identical(alignment$fit$means, fit$means)
  expect_silent(.multilpa_cov_check_likelihood(alignment$fit, data))
})

test_that("the measurement signature permutes with the classes", {
  fit <- alignment_profile_fit()
  order <- c(2L, 3L, 1L)
  moved <- .multilpa_permute_profiles(fit, order)
  reference <- .multilpa_profile_signature(fit, fit)
  expect_identical(.multilpa_profile_signature(moved, fit), reference[order, , drop = FALSE])
  # Directly on synthetic blocks, with and without a mixing column.
  block <- list(means = matrix(c(0, 2, 5, 1, -1, 3), 3L),
                variances = matrix(c(1, 4, 2, 1, 1, 9), 3L),
                response_probabilities = list(a = matrix(c(.2, .5, .9, .8, .5, .1), 3L)))
  shuffled <- block
  shuffled$means <- block$means[order, , drop = FALSE]
  shuffled$response_probabilities$a <- block$response_probabilities$a[order, , drop = FALSE]
  base <- .latents_measurement_signature(list(block), list(block), mixing = c(.5, .3, .2))
  expect_identical(.latents_measurement_signature(list(shuffled), list(block),
                                                  mixing = c(.5, .3, .2)[order]),
                   base[order, , drop = FALSE])
  expect_identical(ncol(base), 5L)
  # Means are in units of the reference's within-class standard deviation.
  expect_equal(base[, 1L], block$means[, 1L] / sqrt(mean(block$variances[, 1L])))
  expect_null(.latents_measurement_signature(list(list()), list(list())))
})

test_that("ordinal and count signatures are permutation invariant", {
  set.seed(712)
  n <- 240L
  z <- rep(1:2, each = n / 2)
  data <- data.frame(g = rep(seq_len(n / 3), each = 3),
                     o = ordered(vapply(z, function(k) {
                       sample(1:4, 1L, prob = if (k == 1L) c(.7, .2, .08, .02) else
                         c(.02, .08, .2, .7))
                     }, integer(1))),
                     k = stats::rnbinom(n, mu = c(2, 14)[z], size = 2))
  fit <- quietly(multilpa(data, c("o", "k"), "g", 2L, 1L, ordinal = "o", count = "k",
                          count_model = "negative_binomial", n_starts = 2, seed = 7))
  moved <- .multilpa_permute_profiles(fit, 2:1)
  # The ordinal reference is rebased, so the probabilities agree to rounding.
  expect_equal(.multilpa_profile_signature(moved, fit),
               .multilpa_profile_signature(fit, fit)[2:1, , drop = FALSE],
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(.multilpa_pool_align(moved, fit)$profile_order, 2:1)
})

test_that("rebased logits keep every class probability", {
  softmax <- function(eta) exp(eta) / rowSums(exp(eta))
  full <- cbind(matrix(c(0.4, -1.2, 2, 0.3), 2L), 0)
  order <- c(3L, 1L, 2L)
  last <- .latents_rebase_logits(full, order)
  expect_identical(unname(last[, 3L]), c(0, 0))
  expect_equal(softmax(last), softmax(full)[, order], tolerance = 1e-14)
  first_full <- cbind(0, full[, 1:2])
  first <- .latents_rebase_logits(first_full, order, reference = "first")
  expect_identical(unname(first[, 1L]), c(0, 0))
  expect_equal(softmax(first), softmax(first_full)[, order], tolerance = 1e-14)
  # Rebasing and rebasing back with the inverse is the identity.
  expect_equal(.latents_rebase_logits(last, order(order)), full, tolerance = 1e-14)
})

test_that("mismatched shapes are refused", {
  expect_error(.multilpa_match_order(matrix(0, 2L, 3L), matrix(0, 3L, 3L)),
               "same shape")
  expect_error(.latents_measurement_signature(list(list()), list()),
               "same length")
  expect_error(.latents_rebase_logits(matrix(0, 2L, 3L), 1:2), "permutation")
  expect_error(.latents_rebase_logits(1:3, 1:3), "matrix")
})

test_that("the seed helper restores the caller's random state", {
  withr::local_preserve_seed()
  draw <- function(seed) {
    .latents_local_seed(seed)
    stats::runif(2L)
  }
  set.seed(99)
  before <- .Random.seed
  first <- draw(7)
  expect_identical(.Random.seed, before)
  set.seed(7)
  expect_identical(first, stats::runif(2L))
  # Without a seed, nothing is scoped: the caller's stream advances.
  set.seed(99)
  unseeded <- draw(NULL)
  set.seed(99)
  expect_identical(unseeded, stats::runif(2L))
  # A session without a random state gets none back.
  rm(".Random.seed", envir = globalenv())
  draw(7)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # The state comes back on an error as well.
  set.seed(3)
  before <- .Random.seed
  expect_error((function() {
    .latents_local_seed(1)
    stats::runif(1L)
    stop("boom")
  })(), "boom")
  expect_identical(.Random.seed, before)
  # The expression form does the same.
  set.seed(4)
  before <- .Random.seed
  expect_identical(.mixture_with_seed(7, stats::runif(2L)), first)
  expect_identical(.Random.seed, before)
})

test_that("a seeded fit leaves the caller's stream untouched", {
  withr::local_preserve_seed()
  set.seed(11)
  before <- .Random.seed
  quietly(multilpa(engagement_small, c("browse", "lectures"), "student", 2, 1,
                   n_starts = 2, seed = 5))
  expect_identical(.Random.seed, before)
})

test_that("categorical draws follow the inverse-CDF definition", {
  withr::local_preserve_seed()
  set.seed(1)
  probabilities <- matrix(stats::runif(60L), 20L)
  probabilities <- probabilities / rowSums(probabilities)
  set.seed(2)
  drawn <- .latents_draw_rows(probabilities)
  set.seed(2)
  uniforms <- stats::runif(nrow(probabilities))
  expected <- vapply(seq_len(nrow(probabilities)), function(i) {
    which(cumsum(probabilities[i, ]) >= uniforms[i])[1L]
  }, integer(1))
  expect_identical(drawn, expected)
  # One column still consumes one uniform per row.
  set.seed(2)
  expect_identical(.latents_draw_rows(matrix(1, 4L, 1L)), rep(1L, 4L))
  next_uniform <- stats::runif(1L)
  set.seed(2)
  expect_identical(next_uniform, stats::runif(5L)[5L])
})

test_that("a uniform above a short cumulative total falls in the last category", {
  withr::local_preserve_seed()
  # A seed whose first uniform exceeds the row's total of 0.9999.
  seed <- Find(function(s) {
    set.seed(s)
    stats::runif(1L) > 0.9999
  }, seq_len(100000L))
  short <- matrix(c(0.3, 0.3, 0.3999), 1L)
  set.seed(seed)
  uniform <- stats::runif(1L)
  expect_gt(uniform, sum(short))
  set.seed(seed)
  expect_identical(.latents_draw_rows(short), 3L)
  # The regression mixture's former rule gives the same category.
  expect_identical(as.integer(pmin(1L + sum(uniform > cumsum(short)), ncol(short))), 3L)
})

test_that("the cluster bootstrap records each failure with its reason", {
  withr::local_preserve_seed()
  data <- data.frame(g = rep(1:4, each = 2), y = seq_len(8))
  calls <- 0L
  refit <- function(resampled) {
    calls <<- calls + 1L
    if (calls == 1L) stop("refit refused")
    list(converged = calls != 2L, value = mean(resampled$y), groups = max(resampled$g))
  }
  extract <- function(fit) if (calls == 3L) "could not be matched" else c(fit$value, fit$groups)
  set.seed(8)
  result <- .latents_cluster_bootstrap(data, data$g, 4L, "g", 4L, refit, extract,
                                       c("mean", "groups"))
  expect_identical(result$messages[1:3], c("refit refused", "the replicate did not converge",
                                           "could not be matched"))
  expect_true(is.na(result$messages[4L]))
  expect_true(all(is.na(result$estimates[1:3, ])))
  expect_identical(colnames(result$estimates), c("mean", "groups"))
  # Resampled groups are renumbered, one per draw.
  expect_identical(unname(result$estimates[4L, "groups"]), 4)
  # Draws are one sample of the groups per replicate, in replicate order.
  set.seed(8)
  draws <- lapply(1:4, function(i) sample.int(4L, 4L, replace = TRUE))
  expect_equal(unname(result$estimates[4L, "mean"]),
               mean(unlist(split(data$y, data$g)[draws[[4L]]])))
})

test_that("the likelihood-ratio replicate loop records warnings and errors", {
  replicate <- function(i) {
    if (i == 2L) stop("refit failed")
    warning("start did not converge")
    list(statistic = c(-0.5, NA, 3)[i], valid = TRUE, boundary = FALSE,
         logits_settled = NA, null_replications = 2L, alternative_replications = 1L)
  }
  table <- .latents_blrt_replicates(3L, replicate, muffle_warnings = TRUE)
  expect_identical(table$statistic, c(0, NA, 3))
  expect_identical(table$valid, c(TRUE, FALSE, TRUE))
  expect_identical(table$error, c(NA, "refit failed", NA))
  expect_identical(table$warnings[c(1L, 3L)], rep("start did not converge", 2L))
  expect_warning(.latents_blrt_replicates(1L, replicate, muffle_warnings = FALSE),
                 "start did not converge")
  expect_warning(tally <- .latents_blrt_p_value(table, 1, 3L),
                 class = "latents_failed_replicates")
  expect_true(is.na(tally$p_value))
  expect_identical(tally$n_valid, 2L)
  valid <- table[c(1L, 3L), ]
  tally <- .latents_blrt_p_value(valid, 1, 2L)
  expect_identical(tally$p_value, 2 / 3)
  expect_identical(tally$monte_carlo_se, sqrt(2 / 3 * (1 / 3) / 3))
})
