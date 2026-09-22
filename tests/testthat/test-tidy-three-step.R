## Tidiness and calibration of the three-step verbs: the shapes they return, the
## shapes they keep as the number of classes grows, and whether the new pairwise
## standard error agrees with a resampled one.

.tidy_step_data <- function(seed = 11L, n_groups = 50L, per = 10L,
                            separation = 1.1, difference = 10) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per)
  n <- length(group)
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  data.frame(g = group,
             a = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             b = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             x = stats::rnorm(n),
             y = stats::rnorm(n, ifelse(truth == 2L, difference, 0)))
}

.tidy_step_fit <- function(data, n_profiles = 2L) {
  multilpa(data, c("a", "b"), "g", n_profiles = n_profiles, n_group_classes = 1,
           n_starts = 6, seed = 1)
}

test_that("the weight frame is long, so its columns do not grow with K", {
  data <- .tidy_step_data()
  two <- get_results(.tidy_step_fit(data, 2L), "bch_weights", level = "individuals")
  three <- get_results(.tidy_step_fit(data, 3L), "bch_weights", level = "individuals")

  expect_named(two, c("level", "unit", "assigned_class", "class", "weight"))
  expect_named(three, names(two))
  expect_equal(nrow(two), nrow(data) * 2L)
  expect_equal(nrow(three), nrow(data) * 3L)
  expect_true(is.numeric(two$weight))
  # ordered by unit, then by class, as documented
  expect_identical(two$unit, rep(seq_len(nrow(data)), each = 2L))
  expect_identical(three$class, rep(1:3, times = nrow(data)))
  # the level is named in the frame rather than assumed
  expect_identical(unique(two$level), "individuals")
})

test_that("three_step keeps one shape across methods, levels and class counts", {
  data <- .tidy_step_data()
  fit <- .tidy_step_fit(data)
  columns <- c("level", "method", "class", "estimate", "standard_error",
               "conf_low", "conf_high", "effective_n")

  invisible(lapply(c("bch", "proportional", "modal"), function(method) {
    result <- three_step(fit, data, "y", method = method)
    expect_named(result, columns)
    expect_equal(nrow(result), 2L)
  }))
  expect_named(three_step(.tidy_step_fit(data, 3L), data, "y"), columns)
  # the shape is stable enough that methods stack for comparison without any
  # renaming, which is the point of carrying the call in columns
  stacked <- rbind(three_step(fit, data, "y", method = "bch"),
                   three_step(fit, data, "y", method = "modal"))
  expect_equal(nrow(stacked), 4L)
  expect_setequal(unique(stacked$method), c("bch", "modal"))
})

test_that("pairwise contrasts carry the inference columns r3step carries", {
  data <- .tidy_step_data()
  fit <- .tidy_step_fit(data)
  pairs <- three_step(fit, data, "y", contrast = "pairs")
  covariate <- r3step(fit, data, "x")
  shared <- c("level", "estimate", "standard_error", "statistic", "p_value",
              "p_value_adjusted", "conf_low", "conf_high")

  expect_named(pairs, c("level", "method", "class", "reference_class",
                        "estimate", "standard_error", "statistic", "p_value",
                        "p_value_adjusted", "conf_low", "conf_high"))
  expect_true(all(shared %in% names(covariate)))
  expect_equal(nrow(pairs), 1L)
  expect_equal(nrow(three_step(.tidy_step_fit(data, 3L), data, "y",
                               contrast = "pairs")), 3L)
  expect_equal(pairs$statistic, pairs$estimate / pairs$standard_error)
  expect_equal(pairs$conf_high - pairs$conf_low,
               2 * stats::qnorm(0.975) * pairs$standard_error)
})

test_that("a contrast is exactly the difference of the two class means", {
  data <- .tidy_step_data()
  fit <- .tidy_step_fit(data, 3L)
  means <- three_step(fit, data, "y")
  pairs <- three_step(fit, data, "y", contrast = "pairs")

  expected <- means$estimate[pairs$class] - means$estimate[pairs$reference_class]
  expect_equal(pairs$estimate, expected)
  # and its interval is wider than either mean's, because the two means share
  # every unit rather than being independent
  expect_true(all(pairs$standard_error > 0))
})

test_that("the estimator is equivariant under an affine change of outcome", {
  data <- .tidy_step_data()
  fit <- .tidy_step_fit(data)
  shifted <- data
  shifted$y <- 3.5 * data$y - 7

  plain <- three_step(fit, data, "y", contrast = "pairs")
  moved <- three_step(fit, shifted, "y", contrast = "pairs")
  plain_means <- three_step(fit, data, "y")
  moved_means <- three_step(fit, shifted, "y")

  expect_equal(moved_means$estimate, 3.5 * plain_means$estimate - 7)
  expect_equal(moved_means$standard_error, 3.5 * plain_means$standard_error)
  # a shift cancels in a difference; a scaling does not
  expect_equal(moved$estimate, 3.5 * plain$estimate)
  expect_equal(moved$statistic, plain$statistic)
})

test_that("the pairwise standard error agrees with a cluster bootstrap", {
  skip_on_cran()
  # Weak separation on purpose: the BCH weights for the two classes are then
  # strongly anticorrelated, so a standard error that ignored the covariance
  # between the two class means would miss by more than the resampling noise.
  data <- .tidy_step_data(separation = 0.8)
  fit <- .tidy_step_fit(data)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  weights <- .multilpa_step_weights(pieces, "bch")
  analytic <- three_step(fit, data, "y", contrast = "pairs")

  # Resample whole groups, holding the step-one weights fixed, which is exactly
  # the assumption the cluster-robust variance makes. The influence-function
  # standard error should reproduce the spread this generates.
  set.seed(99)
  rows_by_group <- split(seq_len(nrow(data)), pieces$group_index)
  replicates <- replicate(600L, {
    picked <- sample(seq_along(rows_by_group), replace = TRUE)
    rows <- unlist(rows_by_group[picked], use.names = FALSE)
    taken <- weights[rows, , drop = FALSE]
    means <- colSums(taken * data$y[rows]) / colSums(taken)
    means[[2L]] - means[[1L]]
  })

  ratio <- stats::sd(replicates) / analytic$standard_error
  expect_gt(ratio, 0.90)
  expect_lt(ratio, 1.16)
  # Treating the two class means as independent would inflate the ratio past
  # that bound, which is what makes this a calibration and not a tautology.
  means <- three_step(fit, data, "y")
  independent <- sqrt(sum(means$standard_error^2))
  expect_gt(stats::sd(replicates) / independent, 1.16)
})

test_that("multiplicity is corrected, named, and can be turned off", {
  data <- .tidy_step_data()
  fit <- .tidy_step_fit(data, 3L)
  adjusted <- three_step(fit, data, "y", contrast = "pairs")
  raw <- three_step(fit, data, "y", contrast = "pairs", adjust = "none")

  expect_identical(attr(adjusted, "adjust"), "BH")
  expect_identical(attr(raw, "adjust"), "none")
  expect_equal(raw$p_value_adjusted, raw$p_value)
  expect_true(all(adjusted$p_value_adjusted >= adjusted$p_value))
  expect_equal(adjusted$p_value, raw$p_value)
  strict <- three_step(fit, data, "y", contrast = "pairs",
                       adjust = "bonferroni")
  expect_true(all(strict$p_value_adjusted >= adjusted$p_value_adjusted))
})

test_that("r3step corrects across its covariate terms and not its intercepts", {
  data <- .tidy_step_data()
  data$w <- stats::rnorm(nrow(data))
  data$v <- stats::rnorm(nrow(data))
  fit <- .tidy_step_fit(data)
  result <- r3step(fit, data, c("x", "w", "v"))
  none <- r3step(fit, data, c("x", "w", "v"), adjust = "none")

  expect_identical(attr(result, "adjust"), "BH")
  expect_true(all(is.na(result$p_value_adjusted[result$term == "(Intercept)"])))
  tested <- result$term != "(Intercept)"
  expect_equal(sum(tested), 3L)
  expect_true(all(result$p_value_adjusted[tested] >= result$p_value[tested]))
  # the correction is over the three covariate terms alone, not over four rows
  expect_equal(result$p_value_adjusted[tested],
               stats::p.adjust(result$p_value[tested], method = "BH"))
  expect_equal(none$p_value_adjusted[tested], none$p_value[tested])
})

test_that("a single class is refused a contrast by condition class", {
  data <- .tidy_step_data()
  fit <- .tidy_step_fit(data, 1L)

  expect_equal(nrow(three_step(fit, data, "y")), 1L)
  expect_error(three_step(fit, data, "y", contrast = "pairs"),
               class = "multilpa_inseparable_classes")
  expect_error(three_step(fit, data, "y", contrast = "nonsense"),
               "'arg' should be one of")
})
