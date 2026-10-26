# Round-2 review of selection and three-step (2026-10-26), regressions.
# Each test fails on the code before its fix; see
# validation/review-2026-10-26/round-2/selection.md.

r2_items <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  classes <- sample(1:2, n, TRUE)
  item <- function() {
    vapply(classes, function(k) {
      sample(1:4, 1L, prob = if (k == 1L) c(0.5, 0.3, 0.15, 0.05) else
        c(0.05, 0.15, 0.3, 0.5))
    }, integer(1))
  }
  data.frame(i1 = item(), i2 = item(), i3 = item(), i4 = item())
}

test_that("S-2: nominal and ordinal fits of the same items are compared", {
  items <- r2_items()
  vars <- c("i1", "i2", "i3", "i4")
  nominal <- lca(items, vars, n_classes = 2, n_starts = 2, seed = 1)
  ordinal <- lca(items, vars, n_classes = 2, ordinal = vars, n_starts = 2, seed = 1)
  comparison <- compare_models(nominal = nominal, ordinal = ordinal)
  table <- as.data.frame(comparison)
  expect_identical(table$model, c("nominal", "ordinal"))
  expect_equal(table$log_likelihood,
               c(nominal$log_likelihood, ordinal$log_likelihood))
  expect_equal(table$bic_groups, c(
    -2 * nominal$log_likelihood + nominal$n_parameters * log(nrow(items)),
    -2 * ordinal$log_likelihood + ordinal$n_parameters * log(nrow(items))))
  # The adjacent-category ordinal model is nested in the nominal one.
  expect_gte(nominal$log_likelihood, ordinal$log_likelihood - 1e-6)
  expect_lt(ordinal$n_parameters, nominal$n_parameters)
})

test_that("S-2: a changed response is still refused across nominal and ordinal", {
  items <- r2_items()
  vars <- c("i1", "i2", "i3", "i4")
  changed <- items
  changed$i3[[5L]] <- if (changed$i3[[5L]] == 1L) 2L else 1L
  nominal <- lca(items, vars, n_classes = 2, n_starts = 2, seed = 1)
  ordinal <- lca(changed, vars, n_classes = 2, ordinal = vars, n_starts = 2, seed = 1)
  expect_error(compare_models(nominal, ordinal),
               class = "latents_incomparable_models")
})

test_that("S-2: an item modelled as Gaussian is not compared with it as ordinal", {
  items <- r2_items()
  vars <- c("i1", "i2", "i3", "i4")
  gaussian <- lpa(items, vars, n_profiles = 2, n_starts = 2, seed = 1)
  ordinal <- lca(items, vars, n_classes = 2, ordinal = vars, n_starts = 2, seed = 1)
  expect_error(compare_models(gaussian, ordinal),
               class = "latents_incomparable_models")
})

test_that("S-6: a negative variance in the group-class r3step table is NA, with a warning", {
  covariance <- diag(c(0.04, -0.01, 0.0025, 0.09))
  expect_warning(
    table <- .multilpa_r3step_group_frame(
      estimates = c(0.5, 1.2, -0.3, 0.1), covariance = covariance, terms = "x",
      n_free = 1L, n_group_classes = 2L, n_classes = 2L, ci_level = 0.95,
      vcov_type = "observed", adjust = "BH", log_likelihood = -10),
    class = "latents_singular_information")
  expect_true(is.na(table$standard_error[[2L]]))
  expect_true(is.na(table$statistic[[2L]]))
  expect_true(is.na(table$conf_low[[2L]]))
  expect_equal(table$standard_error[-2L], sqrt(c(0.04, 0.0025, 0.09)))
  # The same rule as the pooled table.
  expect_warning(
    pooled <- .multilpa_r3step_frame(
      c(0.5, 1.2), diag(c(0.04, -0.01)), c("(Intercept)", "x"), 1L, 2L,
      "individuals", 0.95, "observed", "BH"),
    class = "latents_singular_information")
  expect_true(is.na(pooled$standard_error[[2L]]))
})

r2_entropy <- function(p) -sum(p[p > 0] * log(p[p > 0]))

r2_two_level <- function(seed = 42L) {
  set.seed(seed)
  n_groups <- 60L; size <- 10L
  g <- rep(seq_len(n_groups), each = size)
  group_class <- sample(1:2, n_groups, TRUE)[g]
  profile <- ifelse(runif(n_groups * size) < ifelse(group_class == 1L, 0.8, 0.25), 1L, 2L)
  data.frame(g = g,
             a = rnorm(n_groups * size, ifelse(profile == 1L, -0.9, 0.9)),
             b = rnorm(n_groups * size, ifelse(profile == 1L, -0.9, 0.9)))
}

test_that("S-1: a single-level fit has one ICL, AWE and CLC under both conventions", {
  data <- r2_two_level()
  fit <- lpa(data, c("a", "b"), n_profiles = 2, n_starts = 3, seed = 1)
  criteria <- get_results(fit, "information_criteria")
  en <- r2_entropy(fit$subject_posteriors)
  n <- nrow(data); q <- fit$n_parameters; deviance <- -2 * fit$log_likelihood
  expect_equal(criteria$icl_groups, deviance + q * log(n) + 2 * en)
  expect_equal(criteria$icl_groups, criteria$icl_individual)
  expect_equal(criteria$awe_groups, criteria$awe_individual)
  expect_equal(criteria$clc_groups, deviance + 2 * en)
  expect_equal(criteria$clc_groups, criteria$clc_individual)
})

test_that("S-1: two-level criteria charge the complete-data entropy of group classes and profiles", {
  data <- r2_two_level()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                  n_starts = 10, seed = 1)
  # Independent: P(X_i = k | C = h, y_i) from the fitted densities.
  density <- vapply(1:2, function(k) {
    stats::dnorm(data$a, fit$means[k, 1], sqrt(fit$variances[k, 1])) *
      stats::dnorm(data$b, fit$means[k, 2], sqrt(fit$variances[k, 2]))
  }, numeric(nrow(data)))
  groups <- fit$group_posteriors
  within <- vapply(1:2, function(h) {
    conditional <- sweep(density, 2L, fit$profile_probabilities[h, ], "*")
    conditional <- conditional / rowSums(conditional)
    as.vector(rowsum(-rowSums(conditional * log(conditional)), fit$group_index,
                     reorder = FALSE))
  }, numeric(nrow(groups)))
  joint <- r2_entropy(groups) + sum(groups * within)
  q <- fit$n_parameters; deviance <- -2 * fit$log_likelihood
  criteria <- get_results(fit, "information_criteria")
  expect_equal(criteria$icl_groups, deviance + q * log(60) + 2 * joint)
  expect_equal(criteria$icl_individual, deviance + q * log(nrow(data)) + 2 * joint)
  expect_equal(criteria$clc_groups, criteria$clc_individual)
  expect_equal(criteria$clc_groups, deviance + 2 * joint)
  expect_equal(criteria$awe_groups, deviance + 2 * joint + 2 * q * (1.5 + log(60)))
})

test_that("S-1: CLC no longer selects the largest single-level candidate by deviance alone", {
  data <- r2_two_level()
  grid <- suppressMessages(enumerate_lpa(data, c("a", "b"), n_profiles = 1:4,
                                         model = "VVI", n_starts = 3, seed = 1))
  best <- get_results(grid, "best")
  clc <- best[best$criterion == "clc", ]
  expect_identical(clc$candidate[clc$convention == "groups"],
                   clc$candidate[clc$convention == "individuals"])
})

test_that("S-1, S-7: the printed definitions name the complete-data soft entropy", {
  data <- r2_two_level()
  fit <- lpa(data, c("a", "b"), n_profiles = 2, n_starts = 2, seed = 1)
  definitions <- get_results(fit, "information_criteria", format = "long",
                             definitions = TRUE)
  icl <- definitions$definition[definitions$criterion == "icl"]
  expect_true(all(grepl("EN", icl)))
  expect_true(all(grepl("complete-data", icl)))
})

r2_weighted <- function() {
  set.seed(7)
  n <- 400L
  classes <- sample(1:2, n, TRUE)
  data <- data.frame(id = seq_len(n),
                     a = rnorm(n, ifelse(classes == 1L, -0.8, 0.8)),
                     b = rnorm(n, ifelse(classes == 1L, -0.8, 0.8)))
  data$w <- ifelse(classes == 1L, 4, 1)
  data
}

test_that("S-3: AvePP and the odds of correct classification are weighted together", {
  fit <- lpa(r2_weighted(), c("a", "b"), n_profiles = 2, weights = "w",
             n_starts = 5, seed = 1)
  posteriors <- fit$subject_posteriors
  w <- unname(fit$sampling_weights)
  modal <- max.col(posteriors, ties.method = "first")
  avepp <- vapply(1:2, function(k) stats::weighted.mean(posteriors[modal == k, k],
                                                       w[modal == k]), numeric(1))
  share <- unname(colSums(posteriors * w)) / sum(w)
  table <- get_results(fit, "classification")
  table <- table[table$level == "individuals", ]
  expect_equal(table$average_posterior, avepp)
  expect_equal(table$odds_correct_classification,
               (avepp / (1 - avepp)) / (share / (1 - share)))
  expect_identical(table$n_modal, tabulate(modal, 2L))
})

test_that("S-4: entropy, the error matrix and BCH weights use the sampling weights", {
  fit <- lpa(r2_weighted(), c("a", "b"), n_profiles = 2, weights = "w",
             n_starts = 5, seed = 1)
  posteriors <- fit$subject_posteriors
  w <- unname(fit$sampling_weights)
  terms <- -rowSums(posteriors * log(posteriors))
  entropy <- get_results(fit, "entropy")
  entropy <- entropy[entropy$level == "individuals", ]
  expect_equal(entropy$entropy_sum, sum(w * terms))
  expect_equal(entropy$relative_entropy, 1 - sum(w * terms) / (sum(w) * log(2)))
  criteria <- get_results(fit, "information_criteria")
  expect_equal(criteria$icl_individual, -2 * fit$log_likelihood +
                 fit$n_parameters * log(400) + 2 * sum(w * terms))
  modal <- max.col(posteriors, ties.method = "first")
  errors <- t(vapply(1:2, function(t) {
    vapply(1:2, function(s) sum(w * posteriors[, t] * (modal == s)), numeric(1)) /
      sum(w * posteriors[, t])
  }, numeric(2)))
  table <- get_results(fit, "classification_errors")
  table <- table[table$level == "individuals", ]
  expect_equal(table$probability, errors[cbind(table$true_class, table$assigned_class)])
  bch <- get_results(fit, "bch_weights", level = "individuals")
  expect_equal(bch$weight, as.vector(t(solve(errors)[modal, ])))
})

test_that("S-1: a two-level fit saved without the row entropies reports NA entropy criteria", {
  data <- r2_two_level()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                  n_starts = 3, seed = 1)
  old <- fit
  old$classification_entropy_rows <- NULL
  criteria <- get_results(old, "information_criteria")
  expect_true(all(is.na(unlist(criteria[c("icl_groups", "icl_individual",
                                          "awe_groups", "awe_individual",
                                          "clc_groups", "clc_individual")]))))
  expect_equal(criteria$bic_groups, get_results(fit, "information_criteria")$bic_groups)
  # With one group class nothing needs storing: the row posteriors suffice.
  single <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                     n_starts = 3, seed = 1)
  stored <- get_results(single, "information_criteria")
  single$classification_entropy_rows <- NULL
  expect_equal(get_results(single, "information_criteria"), stored)
})

test_that("S-1: a transition fit keeps one entropy per level", {
  fit <- lta(subset(course_engagement, student <= 40),
             vars = c("browse", "lectures", "forum_read"), id = "student",
             n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
  criteria <- get_results(fit, "information_criteria")
  individual <- r2_entropy(.multilpa_classified_posteriors(fit))
  expect_equal(criteria$clc_individual, -2 * fit$log_likelihood + 2 * individual)
  definitions <- get_results(fit, "information_criteria", format = "long",
                             definitions = TRUE)
  expect_true(all(grepl("level", definitions$definition[definitions$criterion == "icl"])))
})

test_that("S-3: the average-posterior plot draws the weighted AvePP of the table", {
  skip_if_not_installed("ggplot2")
  fit <- lpa(r2_weighted(), c("a", "b"), n_profiles = 2, weights = "w",
             n_starts = 5, seed = 1)
  drawn <- plot_diagnostics(fit, what = "avepp")$data
  table <- get_results(fit, "average_posteriors")
  table <- table[table$level == "individuals", ]
  expect_equal(drawn$average_posterior[order(drawn$row, drawn$column)],
               table$average_posterior)
})
