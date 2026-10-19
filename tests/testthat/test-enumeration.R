test_that("enumeration retains failures and both BIC conventions", {
  set.seed(12)
  d <- data.frame(g = rep(1:30, each = 10), y = c(rnorm(150, -2), rnorm(150, 2)))
  result <- enumerate_classes(d, "y", "g", n_profiles = 1:2, n_group_classes = 1:2,
                             model = NULL, n_starts = 1, seed = 12)
  grid <- as.data.frame(result)
  expect_equal(nrow(grid), 4L)
  expect_match(grid$error[3], "profile")
  expect_error(get_candidate(result, n_profiles = 1, n_group_classes = 2),
               class = "latents_failed_candidate")
  # Every criterion the grid names must actually carry numbers. `aic` was
  # silently an all-missing column named `NA.`, so a finiteness assertion on
  # `grid$aic` passed vacuously against a column that did not exist.
  expect_true("aic" %in% names(grid))
  expect_true(all(is.finite(grid$aic[c(1, 2, 4)])))
  expect_true(all(c("kic", "clc_groups", "clc_individual") %in% names(grid)))
  one_profile <- get_candidate(result, n_profiles = 1, n_group_classes = 1)
  expect_equal(grid$bic_groups[1], one_profile$bic)
  expect_equal(grid$bic_individual[1], one_profile$bic_individual)
  expect_true(is.na(grid$profile_entropy[1]))
  expect_gt(grid$profile_entropy[2], 0)
  expect_lte(grid$profile_entropy[2], 1)
  expect_error(enumerate_classes(d, "y", "g", n_profiles = 0, n_starts = 1), "n_profiles")
  # `n_profiles` names the grid here and the count everywhere else, so a single
  # value is a one-row grid rather than the error the old `profiles =` spelling
  # had to raise when the same word arrived through `...`.
  single <- enumerate_classes(d, "y", "g", n_profiles = 2, n_group_classes = 1,
                              model = "VVI", n_starts = 1, seed = 1)
  expect_equal(nrow(as.data.frame(single)), 1L)
  # A name multilpa() does not take is refused before any candidate is fitted,
  # so the refusal is visible at once instead of as a grid of failed fits.
  expect_error(enumerate_classes(d, "y", "g", n_profiles = 2, n_group_classes = 1,
                                 nonsense = 1, seed = 1, n_starts = 1),
               class = "latents_bad_argument")
})

test_that("simulation preserves group layout and full covariance moments", {
  set.seed(9)
  noise <- matrix(rnorm(2000), 1000, 2) %*% chol(matrix(c(1, .6, .6, 2), 2))
  d <- data.frame(g = rep(1:100, each = 10), a = noise[, 1], b = noise[, 2])
  model <- multilpa(d, c("a", "b"), "g", 1, 1, covariance_model = "full", n_starts = 1)
  simulated <- .multilpa_simulate(model)
  expect_identical(simulated$g, d$g)
  expect_identical(names(simulated), c("a", "b", "g"))
  expect_equal(cov(simulated[c("a", "b")]), cov(d[c("a", "b")]), tolerance = .2)
  expect_false(anyNA(simulated))
  one <- multilpa(d, "a", "g", 1, 1, covariance_model = "full", n_starts = 1)
  single_indicator <- .multilpa_simulate(one)
  expect_equal(dim(single_indicator), c(1000L, 2L))
  expect_identical(single_indicator$g, d$g)
  expect_equal(mean(single_indicator$a), mean(d$a), tolerance = .1)
})

test_that("simulation preserves categorical types and declared factor order", {
  set.seed(7)
  data <- data.frame(
    g = rep(seq_len(20L), each = 5L),
    y = stats::rnorm(100L),
    category = ordered(rep(c("10", "2"), 50L), levels = c("10", "2")),
    rating = rep(c(10L, 2L), each = 50L))
  fit <- quietly(multilpa(data, c("y", "category", "rating"), "g", 2L, 1L,
                           categorical = c("category", "rating"),
                           n_starts = 1, max_iter = 80L, seed = 1L))
  simulated <- .multilpa_simulate(fit)
  expect_identical(class(simulated$category), c("ordered", "factor"))
  expect_identical(levels(simulated$category), c("10", "2"))
  expect_type(simulated$rating, "integer")
  expect_identical(names(simulated), c("y", "category", "rating", "g"))
  encoded <- .multilpa_encode_categorical(
    simulated[c("category", "rating")])
  expect_identical(encoded$levels$category, fit$categorical_levels$category)
  expect_identical(encoded$levels$rating, fit$categorical_levels$rating)
})

test_that("bootstrap withholds p-values if a replicate cannot be fitted", {
  set.seed(53)
  d <- data.frame(g = rep(1:20, each = 5), y = c(rnorm(50, -2), rnorm(50, 2)))
  small <- multilpa(d, "y", "g", 1, 1, n_starts = 2, seed = 12)
  large <- multilpa(d, "y", "g", 2, 1, n_starts = 2, seed = 12)
  testthat::local_mocked_bindings(multilpa = function(...) stop("test optimization failure"))
  expect_warning(result <- bootstrap_lrt(small, large, d, iter = 2, seed = 1),
                 class = "latents_failed_replicates")
  expect_s3_class(result, "multilpa_bootstrap_lrt")
  ## Methods added in this sweep reach the generic only after the NAMESPACE is
  ## regenerated, so they are exercised by explicit call here. Dispatch itself
  ## is asserted in test-tidy-classes.R.
  test <- get_results(result, "test")
  expect_identical(nrow(test), 1L)
  expect_true(is.na(test$p_value))
  expect_equal(test$n_valid, 0L)
  replicates <- get_results(result, "replicates")
  expect_true(all(replicates$error == "test optimization failure"))
  expect_error(plot.multilpa_bootstrap_lrt(result),
               class = "latents_nothing_to_plot")
})

test_that("bootstrap refits generated data and reports finite simulation correction", {
  set.seed(53)
  d <- data.frame(g = rep(1:30, each = 8), y = rnorm(240))
  small <- multilpa(d, "y", "g", 1, 1, variance_model = "equal", n_starts = 1, seed = 8)
  large <- multilpa(d, "y", "g", 2, 1, variance_model = "equal", n_starts = 1, seed = 8,
                      max_iter = 3000, tol = 1e-7)
  rng <- .Random.seed
  result <- bootstrap_lrt(small, large, d, iter = 3, n_starts = 1,
                                  max_iter = 3000, tol = 1e-7, seed = 42)
  expect_identical(.Random.seed, rng)
  test <- get_results(result, "test")
  replicates <- get_results(result, "replicates")
  expect_identical(nrow(test), 1L)
  expect_equal(test$n_valid, 3L)
  expect_equal(test$statistic, max(0, 2 * (large$log_likelihood - small$log_likelihood)))
  expect_equal(test$p_value, (1 + sum(replicates$statistic >= test$statistic)) / 4)
  expect_equal(nrow(replicates), 3L)
  expect_true(all(is.na(replicates$error)))
  expect_true(is.finite(test$monte_carlo_se))
  expect_output(print.multilpa_bootstrap_lrt(result), "Parametric bootstrap")
  expect_output(print.summary_multilpa_bootstrap_lrt(
    summary.multilpa_bootstrap_lrt(result)), "simulated, not chi-square")
  expect_error(bootstrap_lrt(small, small, d, iter = 2), "one class")
  changed <- d
  changed$y <- d$y + 1
  expect_error(bootstrap_lrt(small, large, changed, iter = 2), "reproduce")
  changed$y[1] <- NA_real_
  # Missing values are refused unless both fits integrated them out.
  expect_error(bootstrap_lrt(small, large, changed, iter = 2),
               class = "latents_bad_data")
})

test_that("full covariance print and summary expose residual matrices", {
  set.seed(3)
  d <- data.frame(g = rep(1:10, each = 10), a = rnorm(100), b = rnorm(100))
  fit <- multilpa(d, c("a", "b"), "g", 1, 1, covariance_model = "full", n_starts = 1)
  expect_output(print(fit), "full residual covariance")
  expect_equal(summary(fit)$covariances, fit$covariances)
  expect_output(print(summary(fit)), "-- covariances", fixed = TRUE)
})

test_that("boundary convergence noise is not mistaken for a reversed likelihood", {
  skip_on_cran()
  vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  null_fit <- multilpa(engagement_small, vars, "student", n_profiles = 2,
                       n_group_classes = 1, n_starts = 1, seed = 1, tol = 1e-8)
  alt_fit <- multilpa(engagement_small, vars, "student", n_profiles = 2,
                      n_group_classes = 2, n_starts = 1, seed = 1, tol = 1e-8)
  # Under the null the alternative converges to the null solution, so every
  # replicate statistic sits at zero plus EM noise. EM stops on a RELATIVE
  # change, so that noise is about `2 * tol * |log_likelihood|`, which on a
  # likelihood of this size is several times the fixed -1e-5 window this check
  # used to apply. The result was
  # that healthy replicates were reported as "Nonconvergence or reversed
  # likelihood" and a perfectly ordinary comparison withheld its p-value.
  # A replicate refit on simulated data can land a small class; that warning is
  # recorded in the replicate table and also propagates, so muffle just it.
  # Eight replicates rather than the twelve this test used on the smaller
  # bundled dataset it was written against: each replicate now refits 1422
  # observations on five indicators. The claim is that EVERY replicate is
  # valid, and eight carry it as well as twelve did.
  result <- quietly(bootstrap_lrt(null_fit, alt_fit, iter = 8, n_starts = 1,
                                  max_iter = 2000, tol = 1e-8, seed = 7))
  test <- get_results(result, "test")
  expect_identical(test$n_valid, 8L)
  expect_false(is.na(test$p_value))

  replicates <- get_results(result, "replicates")
  expect_true(all(replicates$valid))
  expect_true(all(is.na(replicates$error)))
  # A replicate that raised nothing reports NA, not an empty string sitting
  # beside an NA-valued error column. Some replicates legitimately do warn, so
  # the contract is "never empty", not "always absent".
  expect_false(any(!is.na(replicates$warnings) & !nzchar(replicates$warnings)))

  # An unconverged model is refused before any statistic is formed, by class.
  # This guard fires ahead of the reversed-likelihood one, and correctly so: a
  # reversal means the alternative was badly optimised, which is what this
  # catches. That ordering makes `latents_reversed_likelihood` hard to reach
  # from here, so it is not asserted in this test.
  poor_start <- list(means = matrix(rep(c(0, 1), each = length(vars)),
                                    nrow = 2, byrow = TRUE),
                     variances = matrix(1, 2, length(vars)),
                     profile_probabilities = matrix(0.5, 2, 2),
                     group_probabilities = c(0.5, 0.5))
  unconverged <- quietly(multilpa(engagement_small, vars, "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  n_starts = 1, seed = 1, max_iter = 0,
                                  start = poor_start, tol = 1e-8))
  expect_error(bootstrap_lrt(null_fit, unconverged, iter = 2, n_starts = 1,
                             seed = 1),
               class = "latents_no_converge")
})

test_that("covariance models are named with `model`, and a stray argument is refused", {
  grid <- quietly(enumerate_classes(engagement_small, c("browse", "lectures"),
                                    "student", n_profiles = 2, n_group_classes = 1,
                                    model = c("VVI", "EEI"), n_starts = 1,
                                    max_iter = 20, seed = 1))
  table <- as.data.frame(grid)
  expect_identical(sort(table$model), c("EEI", "VVI"))
  expect_s3_class(get_candidate(grid, n_profiles = 2, n_group_classes = 1,
                                model = "EEI"), "multilpa")
  expect_true("model" %in% names(get_results(grid, "best")))
  # The old name would otherwise fail inside every candidate.
  expect_error(enumerate_classes(engagement_small, "browse", "student",
                                 n_profiles = 2, n_group_classes = 1,
                                 structure = "VVI", n_starts = 1),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(engagement_small, "browse", "student",
                                 n_profiles = 2, n_group_classes = 1,
                                 structure = "VVI", n_starts = 1),
               "`model`")
})

test_that("printing a grid shows the criteria and the diagnostics", {
  grid <- quietly(enumerate_classes(engagement_small, c("browse", "lectures"),
                                    "student", n_profiles = 1:2,
                                    n_group_classes = 1:2, n_starts = 1,
                                    max_iter = 20, seed = 1))
  printed <- paste(utils::capture.output(print(grid)), collapse = "\n")
  invisible(lapply(c("aic", "bic_groups", "bic_individual", "icl_individual",
                     "group_entropy", "boundary", "n_best_replicated"),
                   function(column) expect_match(printed, column, fixed = TRUE)))
})

test_that("the enumeration plot draws several criteria as panels or plots", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grid <- quietly(enumerate_classes(engagement_small, c("browse", "lectures"),
                                    "student", n_profiles = 1:3,
                                    n_group_classes = 1:2, n_starts = 1,
                                    max_iter = 20, seed = 1))
  combined <- expect_plot(plot(grid))
  # One panel per distinct criterion.
  expect_identical(levels(combined$data$title),
                   names(.multilpa_distinct_criteria(
                     as.data.frame(grid),
                     c("aic", "bic_groups", "bic_individual", "icl_individual"))))
  expect_plot(plot(grid, what = "bic_groups"))
  expect_plot(plot(grid, what = c("aic", "bic_individual")))
  # Separate plots: one per criterion.
  separate <- expect_plots(plot(grid, what = c("aic", "bic_groups"),
                                combine = FALSE))
  expect_length(separate, 2L)
  expect_plot(plot(grid, what = "tree", combine = TRUE))
  expect_error(plot(grid, combine = NA))
  expect_error(plot(grid, what = c("aic", "not_a_criterion")),
               class = "latents_unknown_criterion")
})

test_that("enumerate_classes() enumerates single-level models with id = NULL", {
  flowers <- stats::setNames(iris[1:4], c("sl", "sw", "pl", "pw"))
  # The two-level verb says what it fitted and which verbs fit it by name.
  expect_message(
    grid <- enumerate_classes(flowers, names(flowers), id = NULL,
                              n_profiles = 1:3, model = c("VVV", "EEE"),
                              n_starts = 1, seed = 1),
    class = "latents_single_level")
  # Leaving `id` out is refused and points to the single-level verbs.
  expect_error(enumerate_classes(flowers, names(flowers), n_profiles = 1:2, n_starts = 1),
               class = "latents_bad_argument")
  table <- as.data.frame(grid)
  expect_identical(nrow(table), 6L)
  expect_true(all(table$n_group_classes == 1L))
  expect_true(all(is.na(table$error)))
  # The single-level notice is given once, not recorded against every row.
  expect_false(any(grepl("single-level", table$warnings)))
  # Each cell is the fit multilpa() makes on its own.
  alone <- withCallingHandlers(
    multilpa(flowers, names(flowers), id = NULL, n_profiles = 3, n_starts = 1,
             seed = 1, volume = "varying", shape = "varying",
             orientation = "varying"),
    latents_single_level = function(w) invokeRestart("muffleMessage"))
  expect_equal(get_candidate(grid, n_profiles = 3, n_group_classes = 1,
                             model = "VVV")$log_likelihood,
               alone$log_likelihood)
  expect_error(enumerate_classes(flowers, names(flowers), id = NULL,
                                 n_profiles = 1:2, n_group_classes = 1:2, n_starts = 1),
               class = "latents_bad_argument")
})

test_that("`model` takes the basic set by default, all 14, or codes", {
  flowers <- stats::setNames(iris[1:4], c("sl", "sw", "pl", "pw"))
  enumerate <- function(...) {
    suppressMessages(enumerate_classes(flowers, names(flowers), id = NULL,
                                       n_profiles = 1:2, n_starts = 1,
                                       seed = 1, ...))
  }
  basic <- as.data.frame(enumerate())
  # Equal or varying variances, crossed with covariances absent or present.
  expect_identical(unique(basic$model), c("EEI", "VVI", "EEE", "VVV"))
  expect_identical(nrow(basic), 8L)
  # The default is the set named "basic", not a separate path.
  expect_equal(as.data.frame(enumerate(model = "basic"))$log_likelihood,
               basic$log_likelihood)
  # Codes still work, and a set combines with codes without duplicates.
  expect_identical(unique(as.data.frame(enumerate(model = "VEV"))$model), "VEV")
  expect_identical(unique(as.data.frame(enumerate(model = c("basic", "EEE",
                                                             "VEV")))$model),
                   c("EEI", "VVI", "EEE", "VVV", "VEV"))
  # NULL crosses nothing and fits what multilpa() fits by default.
  plain <- as.data.frame(enumerate(model = NULL))
  expect_identical(nrow(plain), 2L)
  default_fit <- lpa(flowers, names(flowers), n_profiles = 2, n_starts = 1,
                     seed = 1)
  expect_identical(unique(plain$model), default_fit$covariance_structure)
  expect_error(enumerate(model = "most"), class = "latents_bad_argument")
})

test_that("`model = \"all\"` fits the 14 structures", {
  flowers <- stats::setNames(iris[1:4], c("sl", "sw", "pl", "pw"))
  grid <- suppressMessages(enumerate_classes(flowers, names(flowers),
    id = NULL, n_profiles = 2, model = "all", n_starts = 1, seed = 1))
  expect_setequal(as.data.frame(grid)$model, latents:::.multilpa_structures())
})

test_that("structures the data cannot distinguish are fitted once", {
  set.seed(3)
  one <- data.frame(y = c(rnorm(50), rnorm(50, 3, 2)))
  # One indicator: only equal or varying variance remains.
  grid <- suppressMessages(enumerate_classes(one, "y", id = NULL,
    n_profiles = 2, n_starts = 1, seed = 1))
  expect_identical(as.data.frame(grid)$model, c("EEI", "VVI"))
  grid_all <- suppressMessages(enumerate_classes(one, "y", id = NULL,
    n_profiles = 2, model = "all", n_starts = 1, seed = 1))
  expect_identical(nrow(as.data.frame(grid_all)), 2L)
  # Only categorical indicators: no structure, one candidate per count, and
  # get_candidate() needs no `model`.
  answers <- data.frame(a = factor(sample(c("x", "y"), 100, TRUE)),
                        b = factor(sample(c("x", "y"), 100, TRUE)))
  classes <- suppressMessages(enumerate_classes(answers, c("a", "b"),
    id = NULL, n_profiles = 1:2, categorical = c("a", "b"), n_starts = 1,
    seed = 1))
  expect_identical(nrow(as.data.frame(classes)), 2L)
  expect_s3_class(get_candidate(classes, n_profiles = 2), "multilpa")
  expect_error(enumerate_classes(answers, c("a", "b"), id = NULL,
                                 n_profiles = 1:2, model = "basic",
                                 categorical = c("a", "b"), n_starts = 1),
               class = "latents_bad_argument")
})

test_that("a structure set another way replaces the basic set", {
  flowers <- stats::setNames(iris[1:4], c("sl", "sw", "pl", "pw"))
  grid <- suppressMessages(enumerate_classes(flowers, names(flowers),
    id = NULL, n_profiles = 1:2, covariance_model = "full", n_starts = 1,
    seed = 1))
  expect_identical(nrow(as.data.frame(grid)), 2L)
  expect_error(enumerate_classes(flowers, names(flowers), id = NULL,
                                 n_profiles = 1:2, model = "EEE",
                                 covariance_model = "full", n_starts = 1),
               class = "latents_bad_argument")
})

test_that("the grid compares every candidate's AvePP and entropy", {
  set.seed(3)
  d <- data.frame(g = rep(1:30, each = 10),
                  y = c(rnorm(150, -2), rnorm(150, 2)),
                  z = rnorm(300))
  result <- enumerate_classes(d, c("y", "z"), "g", n_profiles = 1:3,
                              n_group_classes = 1:2, model = "VVI",
                              n_starts = 2, seed = 3)
  grid <- as.data.frame(result)
  summaries <- c(paste0("avepp_", c("min", "max", "mean", "below")),
                 paste0("group_avepp_", c("min", "max", "mean", "below")))
  # Next to the entropies they are read with.
  after <- match("group_entropy", names(grid))
  expect_identical(names(grid)[after + seq_along(summaries)], summaries)
  expect_type(grid$avepp_below, "integer")
  # Calibration against the definition: each class's mean posterior over the
  # units modally assigned to it.
  fit <- get_candidate(result, n_profiles = 3, n_group_classes = 1)
  posteriors <- fit$subject_posteriors
  assigned <- max.col(posteriors, ties.method = "first")
  avepp <- vapply(seq_len(3L), \(k) mean(posteriors[assigned == k, k]),
                  numeric(1L))
  row <- subset(grid, n_profiles == 3 & n_group_classes == 1)
  expect_equal(row$avepp_min, min(avepp))
  expect_equal(row$avepp_max, max(avepp))
  expect_equal(row$avepp_mean, mean(avepp))
  expect_identical(row$avepp_below, sum(avepp < 0.8))
  # The cutoff is an argument; the extremes bracket every class.
  # "Below" is strict: a class whose units all have posterior exactly 1 (this
  # grid's 2 x 2 candidate) is not below a cutoff of 1.
  everything <- get_results(result, "candidates", cutoff = 1)
  fitted <- subset(everything, n_profiles > 1 & !is.na(log_likelihood))
  expect_true(all(fitted$avepp_below <= fitted$n_profiles))
  imperfect <- subset(fitted, avepp_max < 1)
  expect_gt(nrow(imperfect), 0L)
  expect_identical(imperfect$avepp_below, imperfect$n_profiles)
  none <- get_results(result, "candidates", cutoff = min(avepp) / 2)
  expect_identical(subset(none, n_profiles == 3 & n_group_classes == 1)$
                     avepp_below, 0L)
  expect_true(all(with(subset(grid, !is.na(avepp_min)),
                       avepp_min <= avepp_mean &
                         avepp_mean <= avepp_max &
                         avepp_max <= 1)))
  expect_error(get_results(result, "candidates", cutoff = 2), "cutoff")
  expect_error(as.data.frame(result, cutoff = 0.7),
               class = "latents_bad_argument")
  # One class classifies nothing, at either level.
  expect_true(all(is.na(subset(grid, n_profiles == 1)$avepp_min)))
  expect_true(all(is.na(subset(grid, n_group_classes == 1)$group_avepp_below)))

  classification <- get_results(result, "classification")
  expect_s3_class(classification, "data.frame")
  expect_true(all(c("n_profiles", "n_group_classes", "model", "level", "class",
                    "average_posterior", "relative_entropy")
                  %in% names(classification)))
  # Invariants: one row per class of every fitted level with two or more
  # classes; its AvePPs give the candidate's summaries; the entropy is the
  # candidate's entropy.
  expect_false(any(classification$n_profiles == 1 &
                     classification$level == "individuals"))
  expect_false(any(classification$n_group_classes == 1 &
                     classification$level == "groups"))
  three <- subset(classification, n_profiles == 3 & n_group_classes == 1)
  expect_equal(nrow(three), 3L)
  expect_equal(three$average_posterior, avepp)
  expect_equal(unique(three$relative_entropy), row$profile_entropy)
  failed <- subset(grid, is.na(log_likelihood))
  expect_equal(nrow(merge(failed, classification,
                          by = c("n_profiles", "n_group_classes"))), 0L)
  expect_identical(get_results(summary(result), "classification"), classification)

  # Computed from the stored fits when read, so a grid fitted before the
  # columns existed reports them without refitting.
  expect_false(any(summaries %in% names(result$table)))
  # A grid that stored such columns itself is read the same way: the stored
  # copy is replaced by the values its fits give.
  old <- result
  old$table$avepp_min <- -1
  old$table$profile_min_average_posterior <- -1
  expect_identical(as.data.frame(old), grid)
  expect_identical(get_results(summary(old), "candidates"), grid)
  printed <- paste(utils::capture.output(print(old)), collapse = "\n")
  invisible(lapply(summaries, \(column) expect_match(printed, column, fixed = TRUE)))
})

test_that("AvePP summaries handle empty and single classes", {
  # Class 2 has no modally assigned unit: no AvePP, counted below the cutoff.
  posteriors <- cbind(c(0.9, 0.8, 0.7), c(0.1, 0.2, 0.3))
  summary_empty <- .multilpa_avepp_summary(posteriors, cutoff = 0.5)
  expect_equal(unname(summary_empty), c(0.8, 0.8, 0.8, 1))
  expect_true(all(is.na(.multilpa_avepp_summary(NULL, 0.8))))
  expect_true(all(is.na(.multilpa_avepp_summary(matrix(1, 3, 1), 0.8))))
  separated <- cbind(c(0.9, 0.2, 0.4), c(0.1, 0.8, 0.6))
  expect_equal(.multilpa_avepp_summary(separated, cutoff = 0.75),
               c(min = 0.7, max = 0.9, mean = 0.8, below = 1))
})

test_that("a candidate records a readable call, not a copy of the data", {
  set.seed(4)
  frame <- data.frame(g = rep(1:30, each = 10), y = rnorm(300), z = rnorm(300))
  grid <- enumerate_classes(frame, c("y", "z"), "g", n_profiles = 1:2,
                            n_group_classes = 1, model = "VVI", n_starts = 1,
                            seed = 4)
  call <- get_candidate(grid, 2, 1)$call
  expect_identical(call[[1L]], as.name("multilpa"))
  expect_identical(call$data, as.name("frame"))
  expect_identical(call$n_profiles, 2L)
  expect_false(any(vapply(as.list(call), \(part) is.data.frame(part) ||
                            is.function(part), logical(1L))))
  # The grid no longer grows by one data copy per candidate: a candidate's
  # call is a few hundred bytes whatever the data size.
  expect_lt(length(serialize(call, NULL)), 2000)
})

test_that("compare_models() compares candidates of a grid by name or counts", {
  set.seed(5)
  n_groups <- 50
  g <- rep(seq_len(n_groups), each = 10)
  group_class <- rep(sample(1:2, n_groups, TRUE), each = 10)
  profile <- ifelse(runif(length(g)) < ifelse(group_class == 1, 0.8, 0.2), 1, 2)
  d <- data.frame(g = g, y = rnorm(length(g), 3 * profile),
                  z = rnorm(length(g), profile))
  grid <- quietly(enumerate_classes(d, c("y", "z"), "g", n_profiles = 2:3,
                                    n_group_classes = 1:2, model = "VVI",
                                    n_starts = 2, seed = 5))
  compared <- compare_models(grid, "profiles2_groups2", "profiles3_groups2",
                             c(3, 1))
  expect_s3_class(compared, c("latents_comparison", "data.frame"))
  expect_identical(compared$model,
                   c("profiles2_groups2", "profiles3_groups2", "profiles3_groups1"))
  # The numbers are the candidates table's, row for row.
  candidates <- as.data.frame(grid)
  shared <- c("log_likelihood", "n_parameters", "aic", "bic_groups",
              "bic_individual", "icl_individual", "profile_entropy",
              "group_entropy", "avepp_min", "avepp_below",
              "group_avepp_mean")
  matched <- merge(compared, candidates, by.x = "model", by.y = "candidate",
                   suffixes = c("", "_grid"), sort = FALSE)
  invisible(lapply(shared, function(column) {
    expect_equal(matched[[column]], matched[[paste0(column, "_grid")]],
                 info = column)
  }))
  # Schwarz weights: a distribution with the best model at delta zero.
  expect_equal(sum(compared$bic_weight), 1)
  expect_equal(min(compared$delta_bic), 0)
  expect_equal(compared$delta_bic, compared$bic_groups - min(compared$bic_groups))
  # Order of the candidates does not change any row.
  reversed <- compare_models(grid, c(3, 1), "profiles3_groups2", "profiles2_groups2")
  expect_equal(as.data.frame(reversed)[3:1, "bic_weight"], compared$bic_weight)
  # The same as passing the fits themselves, named.
  direct <- compare_models(profiles2_groups2 = grid$profiles2_groups2,
                           profiles3_groups2 = grid$profiles3_groups2,
                           profiles3_groups1 = grid$profiles3_groups1)
  expect_equal(as.data.frame(direct), as.data.frame(compared))
  named <- compare_models(grid, small = "profiles2_groups2", large = c(3, 2))
  expect_identical(named$model, c("small", "large"))
  expect_identical(compare_models(grid, c(2, 2), c(3, 2), cutoff = 1)$avepp_below,
                   c(2L, 3L))
  # Selecting by class counts instead of naming.
  by_profiles <- compare_models(grid, n_profiles = 3)
  expect_identical(by_profiles$model, c("profiles3_groups1", "profiles3_groups2"))
  block <- compare_models(grid, n_profiles = 2:3, n_group_classes = 2)
  expect_identical(block$model, c("profiles2_groups2", "profiles3_groups2"))
  expect_identical(nrow(compare_models(grid)), 4L)

  expect_error(compare_models(grid, "profiles2_groups2"), class = "latents_bad_argument")
  expect_error(compare_models(grid, "profiles2_groups2", "profiles9_groups9"),
               class = "latents_unknown_candidate")
  expect_error(compare_models(grid, c(2, 2), c(2, 2, 1)), class = "latents_bad_argument")
  expect_error(compare_models(grid, c(2, 2), c(5, 1)), class = "latents_unknown_candidate")
  expect_error(compare_models(grid, c(2, 2), n_profiles = 3),
               class = "latents_bad_argument")
  expect_error(compare_models(grid, n_profiles = 7), class = "latents_unknown_candidate")
  expect_error(compare_models(grid, n_profiles = 2, n_group_classes = 1),
               class = "latents_bad_argument")
  expect_error(compare_models(grid$profiles2_groups1, grid$profiles2_groups2,
                              n_profiles = 2), class = "latents_bad_argument")
  other <- quietly(multilpa(d[-1, ], c("y", "z"), "g", 2, 1, n_starts = 1, seed = 1))
  expect_error(compare_models(a = grid$profiles2_groups1, b = other),
               class = "latents_incomparable_models")
  skip_if_not_installed("ggplot2")
  drawn <- expect_plot(plot(compared))
  expect_setequal(levels(drawn$data$criterion), c("AIC", "BIC", "SABIC"))
})

test_that("a selection leaves out candidates that failed, and says so", {
  set.seed(12)
  d <- data.frame(g = rep(1:30, each = 10), y = c(rnorm(150, -2), rnorm(150, 2)))
  grid <- enumerate_classes(d, "y", "g", n_profiles = 1:2, n_group_classes = 1:2,
                            model = NULL, n_starts = 1, seed = 12)
  failed <- subset(as.data.frame(grid), is.na(log_likelihood))$candidate
  expect_identical(failed, "profiles1_groups2")
  expect_message(compared <- compare_models(grid, n_profiles = 1:2),
                 "profiles1_groups2")
  expect_false("profiles1_groups2" %in% compared$model)
  expect_error(compare_models(grid, "profiles1_groups1", "profiles1_groups2"),
               class = "latents_failed_candidate")
  expect_null(grid$profiles1_groups2)
  expect_error(get_candidate(grid, "profiles1_groups2"),
               class = "latents_failed_candidate")
})

test_that("compare_models() compares latent class candidates and fits", {
  set.seed(5)
  n <- 400
  latent <- sample(1:3, n, TRUE)
  item <- function(p) factor(ifelse(runif(n) < p[latent], "yes", "no"))
  d <- data.frame(g = rep(1:40, each = 10), i1 = item(c(.9, .2, .5)),
                  i2 = item(c(.8, .1, .9)), i3 = item(c(.85, .3, .2)),
                  i4 = item(c(.1, .8, .5)))
  items <- paste0("i", 1:4)
  single <- quietly(enumerate_lca(d, items, n_classes = 2:3, n_starts = 2, seed = 1))
  # A single-level latent class grid names candidates by classes alone.
  expect_identical(as.data.frame(single)$candidate, c("classes2", "classes3"))
  by_name <- compare_models(single, "classes2", "classes3")
  expect_equal(as.data.frame(by_name)[-1L],
               as.data.frame(compare_models(single, 2, 3))[-1L])
  # lca() fits directly give the same numbers as the grid's candidates.
  direct <- compare_models(two = single$classes2, three = single$classes3)
  expect_equal(direct$bic_groups, by_name$bic_groups)
  printed <- paste(utils::capture.output(print(by_name)), collapse = "\n")
  expect_match(printed, "Classes", fixed = TRUE)
  expect_false(grepl("Group classes", printed, fixed = TRUE))
  expect_true(all(is.na(by_name$structure)))

  multilevel <- quietly(enumerate_classes(d, items, "g", n_profiles = 2,
                                          n_group_classes = 1:2,
                                          categorical = items, n_starts = 2,
                                          seed = 1))
  expect_identical(names(multilevel$fits), c("classes2_groups1", "classes2_groups2"))
  grouped <- compare_models(multilevel)
  printed <- paste(utils::capture.output(print(grouped)), collapse = "\n")
  expect_match(printed, "Group classes", fixed = TRUE)
  expect_error(compare_models(single, 2, 2.5), class = "latents_bad_argument")
  expect_error(compare_models(single, 2, 0), class = "latents_bad_argument")
})

test_that("one name reaches a candidate everywhere", {
  set.seed(5)
  d <- data.frame(g = rep(1:40, each = 10),
                  y = c(rnorm(200, -2), rnorm(200, 2)), z = rnorm(400))
  grid <- quietly(enumerate_classes(d, c("y", "z"), "g", n_profiles = 1:2,
                                    n_group_classes = 1:2, model = "VVI",
                                    n_starts = 2, seed = 5))
  names <- c("profiles1_groups1", "profiles2_groups1", "profiles1_groups2",
             "profiles2_groups2")
  expect_identical(names(grid$fits), names)
  expect_identical(as.data.frame(grid)$candidate, names)
  expect_identical(names(as.data.frame(grid))[1L], "candidate")
  # Every route returns the identical model.
  model <- get_candidate(grid, 2, 2)
  expect_identical(grid$profiles2_groups2, model)
  expect_identical(grid$fits$profiles2_groups2, model)
  expect_identical(get_candidate(grid, "profiles2_groups2"), model)
  expect_identical(get_candidate(grid, c(2, 2)), model)
  expect_identical(get_candidate(grid, 2, 2, model = "VVI"), model)
  expect_identical(candidate_fit(grid, "profiles2_groups2"), model)
  compared <- compare_models(grid, "profiles2_groups1", "profiles2_groups2")
  expect_identical(get_candidate(compared, "profiles2_groups2"), model)
  # The grid's own components are untouched, and $ still matches as on a list.
  expect_s3_class(grid$table, "data.frame")
  expect_type(grid$fits, "list")
  expect_true(is.call(grid$call))
  expect_identical(grid$tab, grid$table)
  expect_null(grid$profiles9_groups9)
  # An editor completes the components and every candidate.
  expect_identical(utils::.DollarNames(grid, "^profiles2"),
                   c("profiles2_groups1", "profiles2_groups2"))
  expect_true(all(c("table", "fits", "call") %in% utils::.DollarNames(grid, "")))
  # The criteria and classification tables name their candidates.
  picks <- get_results(grid, "best")
  expect_true(all(picks$candidate %in% names))
  classification <- get_results(grid, "classification")
  expect_identical(names(classification)[1L], "candidate")
  expect_true(all(classification$candidate %in% names))

  # A grid saved before names existed is named the same way.
  old <- grid
  old$fits <- unname(old$fits)
  old$table$candidate <- NULL
  expect_identical(as.data.frame(old)$candidate, names)
  expect_identical(old$profiles2_groups2, model)
  expect_identical(get_candidate(old, "profiles2_groups2"), model)

  # Comparison models and per-class rows.
  tables <- get_results(compared, "all")
  expect_named(tables, c("comparison", "classification"))
  expect_null(attr(tables$comparison, "models"))
  expect_identical(unique(tables$classification$model),
                   c("profiles2_groups1", "profiles2_groups2"))
  invisible(lapply(unique(tables$classification$model), function(label) {
    rows <- subset(tables$classification, model == label & level == "individuals")
    expect_equal(min(rows$average_posterior),
                 subset(compared, model == label)$avepp_min)
  }))
  expect_error(get_candidate(compared, "profiles9_groups9"),
               class = "latents_unknown_candidate")
  expect_error(get_candidate(head(compared, 1L), "profiles2_groups2"),
               class = "latents_unknown_candidate")
  stripped <- compared
  attr(stripped, "models") <- NULL
  expect_error(get_candidate(stripped, "profiles2_groups2"),
               class = "latents_incomplete_fit")
  expect_error(get_candidate(grid, "profiles2_groups2", n_group_classes = 2),
               class = "latents_bad_argument")
  expect_error(get_candidate(grid, "two"), class = "latents_unknown_candidate")
})

test_that("names follow the grid: structures, single level, transitions, families", {
  set.seed(6)
  d <- data.frame(g = rep(1:30, each = 10), y = rnorm(300), z = rnorm(300))
  crossed <- quietly(enumerate_lpa(d, c("y", "z"), n_profiles = 1:2,
                                   model = c("VVI", "EEE"), n_starts = 1, seed = 1))
  expect_identical(names(crossed$fits),
                   c("profiles1_VVI", "profiles2_VVI", "profiles1_EEE", "profiles2_EEE"))
  expect_error(get_candidate(crossed, 2), class = "latents_unknown_candidate")
  expect_identical(get_candidate(crossed, "profiles2_EEE"),
                   get_candidate(crossed, 2, model = "EEE"))
  transitions <- quietly(enumerate_classes(
    subset(course_engagement, student <= 30), c("browse", "lectures"), "student",
    time = "sequence", n_profiles = 2:3, n_group_classes = 1, n_starts = 1,
    seed = 1))
  expect_identical(as.data.frame(transitions)$candidate,
                   c("profiles2_groups1", "profiles3_groups1"))
  expect_identical(transitions$profiles3_groups1, get_candidate(transitions, 3, 1))
  ratings <- data.frame(team = rep(1:40, each = 6), climate = rnorm(240))
  families <- quietly(enumerate_classes(ratings, "climate", "team",
                                        family = c("additive", "dispersion"),
                                        n_group_classes = 1:2, n_starts = 1,
                                        seed = 1))
  expect_true("additive_equal_groups2" %in% as.data.frame(families)$candidate)
  expect_identical(families$dispersion_groups2,
                   get_candidate(families, n_group_classes = 2, model = "dispersion"))
  expect_identical(get_candidate(families, "dispersion_groups2"),
                   families$dispersion_groups2)
  expect_error(compare_models(families, "additive_groups1", "additive_groups2"),
               class = "latents_bad_argument")
})

test_that("plot(grid, what = ) draws entropy and AvePP across the candidates", {
  skip_if_not_installed("ggplot2")
  set.seed(5)
  n <- 400
  group <- rep(1:40, each = 10)
  group_class <- rep(sample(1:2, 40, TRUE), each = 10)
  latent <- ifelse(runif(n) < ifelse(group_class == 1, 0.8, 0.2), 1, 2)
  item <- function(p) factor(ifelse(runif(n) < p[latent], "yes", "no"))
  d <- data.frame(g = group, i1 = item(c(.9, .2)), i2 = item(c(.8, .1)),
                  i3 = item(c(.85, .3)), i4 = item(c(.1, .8)))
  items <- paste0("i", 1:4)
  grid <- quietly(enumerate_classes(d, items, "g", n_profiles = 2:3,
                                    n_group_classes = 1:2, categorical = items,
                                    n_starts = 2, seed = 1))
  candidates <- as.data.frame(grid)
  avepp <- expect_plot(plot(grid, what = "avepp"))
  expect_identical(levels(avepp$data$title),
                   c("Lowest AvePP, individuals", "Lowest AvePP, groups"))
  expect_identical(avepp$labels$x, "Number of classes")
  # The drawn values are the candidates table's.
  individuals <- subset(avepp$data, title == "Lowest AvePP, individuals")
  expect_equal(sort(individuals$value), sort(candidates$avepp_min))
  groups <- subset(avepp$data, title == "Lowest AvePP, groups")
  expect_equal(sort(groups$value),
               sort(candidates$group_avepp_min[is.finite(candidates$group_avepp_min)]))
  # The cutoff is drawn and reported.
  cut <- plot_enumeration(grid, what = "avepp", cutoff = 0.9)
  expect_true(any(vapply(cut$layers, \(layer) inherits(layer$geom, "GeomHline"),
                         logical(1))))
  expect_match(cut$labels$subtitle, "0.9", fixed = TRUE)
  entropy <- expect_plot(plot(grid, what = "entropy"))
  expect_equal(sort(subset(entropy$data, title == "Entropy, individuals")$value),
               sort(candidates$profile_entropy))
  both <- expect_plot(plot(grid, what = "diagnostics"))
  expect_length(levels(both$data$title), 4L)
  separate <- plot(grid, what = "diagnostics", combine = FALSE)
  expect_s3_class(separate, "latents_plots")
  expect_named(separate, c("entropy", "avepp"))
  # A single-level grid has no group panels.
  single <- quietly(enumerate_lca(d, items, n_classes = 2:3, n_starts = 2, seed = 1))
  expect_identical(levels(plot(single, what = "diagnostics")$data$title),
                   c("Entropy", "Lowest AvePP"))
  expect_error(plot(grid, what = c("aic", "avepp")), class = "latents_bad_argument")
  expect_error(plot(grid, what = "avepp", cutoff = 2), "cutoff")
  one_class <- quietly(enumerate_lca(d, items, n_classes = 1, n_starts = 1, seed = 1))
  expect_error(plot(one_class, what = "entropy"), class = "latents_nothing_to_plot")
})

test_that("Kneedle reproduces a hand-worked curve and a published value", {
  # Rescaled: x = 0, 1/3, 2/3, 1 and y = 1, 1/9, 1/18, 0. The line is 1 - x,
  # so the points lie 0, 5/9, 5/18, 0 below it, and perpendicular to it
  # those distances shrink by sqrt(2).
  knee <- .latents_kneedle(1:4, c(10, 2, 1.5, 1))
  expect_identical(knee$at, 2L)
  expect_equal(knee$distance, c(0, 5 / 9, 5 / 18, 0) / sqrt(2))
  # A two-level LCA grid's BIC (groups) at 4 group classes, 2 to 8 classes:
  # Kneedle was reported at 4 classes with a perpendicular distance of 0.1915.
  bic <- c(1143453, 1101283, 1083371, 1070010, 1060205, 1051111, 1044007)
  published <- .latents_kneedle(2:8, bic)
  expect_identical((2:8)[published$at], 4L)
  expect_equal(round(published$distance[published$at], 4), 0.1915)
  # A straight line has no knee; nor do fewer than three points.
  expect_identical(.latents_kneedle(1:5, 10:6)$at, NA_integer_)
  expect_identical(.latents_kneedle(1:2, c(5, 1))$at, NA_integer_)
  # The knee does not depend on the criterion's origin or units.
  expect_identical(.latents_kneedle(1:4, 1000 + 7 * c(10, 2, 1.5, 1))$at, 2L)
})

test_that("deltas and knees read the candidates table in both directions", {
  set.seed(5)
  n <- 400
  group <- rep(1:40, each = 10)
  group_class <- rep(sample(1:2, 40, TRUE), each = 10)
  latent <- ifelse(runif(n) < ifelse(group_class == 1, 0.8, 0.2), 1, 2)
  d <- data.frame(g = group, y = rnorm(n, 3 * latent), z = rnorm(n, latent))
  grid <- quietly(enumerate_classes(d, c("y", "z"), "g", n_profiles = 1:4,
                                    n_group_classes = 1:3, model = "VVI",
                                    n_starts = 2, seed = 1))
  candidates <- as.data.frame(grid)
  bic <- function(name) subset(candidates, candidate == name)$bic_groups
  deltas <- get_results(grid, "deltas")
  expect_identical(deltas$candidate, candidates$candidate)
  row <- subset(deltas, candidate == "profiles3_groups2")
  expect_equal(row$delta_profile, bic("profiles2_groups2") - bic("profiles3_groups2"))
  expect_identical(row$from_profile, "profiles2_groups2")
  expect_equal(row$delta_group_class, bic("profiles3_groups1") - bic("profiles3_groups2"))
  expect_identical(row$from_group_class, "profiles3_groups1")
  # No smaller neighbour, or a neighbour that failed: no delta.
  expect_true(is.na(subset(deltas, candidate == "profiles1_groups1")$delta_profile))
  # profiles1_groups2 could not be fitted, so profiles2_groups2 has no
  # profile delta.
  expect_true(is.na(bic("profiles1_groups2")))
  expect_true(is.na(subset(deltas, candidate == "profiles2_groups2")$delta_profile))
  expect_true(all(is.na(subset(deltas, is.na(value))$delta_profile)))
  icl <- get_results(grid, "deltas", using = "icl_individual")
  expect_identical(unique(icl$criterion), "icl_individual")

  knees <- get_results(grid, "knees")
  expect_setequal(unique(knees$direction), c("profiles", "group_classes", "parameters"))
  expect_true(all(knees$knee_candidate %in% c(candidates$candidate, NA)))
  expect_true(all(knees$n_points >= 1L))
  expect_named(get_results(grid, "all"),
               c("candidates", "best", "classification", "deltas", "knees",
                 "tree"))

  skip_if_not_installed("ggplot2")
  views <- c("deltas", "heatmap", "smallest", "trust")
  invisible(lapply(views, function(view) expect_plot(plot(grid, what = view))))
  everything <- plot(grid, what = "selection")
  expect_s3_class(everything, "latents_plots")
  expect_named(everything, c(views, "diagnostics"))
  drawn <- plot(grid, what = "deltas")
  expect_setequal(as.character(unique(drawn$data$panel)),
                  c("Adding a profile", "Adding a group class"))
  expect_equal(sort(subset(drawn$data, panel == "Adding a profile")$delta),
               sort(deltas$delta_profile[is.finite(deltas$delta_profile)]))
  heatmap <- plot(grid, what = "heatmap")
  expect_identical(subset(heatmap$data, tile_label == "lowest")$candidate,
                   subset(get_results(grid, "best"),
                          criterion == "bic" & convention == "groups")$candidate)
  expect_false(any(heatmap$data$tile_label == "best"))
  trust <- plot(grid, what = "trust")
  expect_true("not fitted" %in% trust$data$status)
  expect_error(plot(grid, what = "deltas", using = "nonsense"),
               class = "latents_bad_argument")
  expect_error(plot(grid, what = c("aic", "deltas")), class = "latents_bad_argument")
  expect_error(plot(grid, what = "elbow"), class = "latents_unknown_criterion")
})

test_that("review fixes: alias, transitions, rising lines, empty grids", {
  grid <- quietly(enumerate_lpa(srl, c("self_efficacy", "test_anxiety"),
                                n_profiles = 1:3, model = "VVI", n_starts = 2,
                                seed = 1))
  # The 0.9.21 alias calls the package's function, whatever the caller has
  # named `get_candidate`.
  local({
    get_candidate <- function(...) "not the package's"
    expect_s3_class(candidate_fit(grid, 2), "multilpa")
  })
  # A grid made through a wrapper names the caller's data.
  expect_identical(get_candidate(grid, 2)$call$data, as.name("srl"))

  transitions <- quietly(enumerate_classes(
    subset(course_engagement, student <= 40), c("browse", "lectures"),
    "student", time = "sequence", n_profiles = 2:3, n_group_classes = 1,
    n_starts = 1, seed = 1))
  compared <- compare_models(transitions)
  expect_gt(nrow(get_results(compared, "classification")), 0L)
  skip_if_not_installed("ggplot2")
  expect_error(plot(transitions, what = "bic_groups"), class = "latents_bad_argument")
  # 0.9.21 matched these with match.arg(), so abbreviations still work.
  expect_identical(plot(transitions, what = "bic_i")$labels$y,
                   plot(transitions, what = "bic_individual")$labels$y)
  expect_warning(plot(transitions, colour = "red"), "\"bic\"",
                 class = "latents_unused_argument")

  # A criterion that rises along a line has no knee.
  rising <- .latents_kneedle(1:4, c(100, 101, 150, 300))
  expect_identical(rising$at, NA_integer_)
  expect_true(rising$rising)

  # A grid where nothing could be fitted gives empty tables, not NULL.
  d <- data.frame(g = rep(1:30, each = 10), y = stats::rnorm(300))
  empty <- quietly(enumerate_classes(d, "y", "g", n_profiles = 1,
                                     n_group_classes = 2, model = NULL,
                                     n_starts = 1, seed = 1))
  expect_s3_class(get_results(empty, "classification"), "data.frame")
  expect_identical(nrow(get_results(empty, "classification")), 0L)
  expect_identical(nrow(get_results(empty, "knees")), 0L)
  expect_error(plot(empty, what = "heatmap"), class = "latents_nothing_to_plot")

  # A deltas panel whose lines have one point draws without warnings.
  two <- quietly(enumerate_classes(d, "y", "g", n_profiles = 1:3,
                                   n_group_classes = 1:2, model = "VVI",
                                   n_starts = 1, seed = 1))
  expect_no_warning(ggplot2::ggplot_build(plot(two, what = "deltas")))
  comparison <- compare_models(two, n_profiles = 2:3)
  expect_error(get_candidate(comparison), class = "latents_bad_argument")
  expect_error(compare_models(two, n_profiles = 0), class = "latents_bad_argument")
})
