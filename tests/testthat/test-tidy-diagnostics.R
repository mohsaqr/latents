# The tidy contract of the diagnostic verbs: one meaning per column, no magic
# strings, no documentation stored as data, one shape per verb.

.tidy_fit <- function() {
  set.seed(11)
  group <- rep(seq_len(30L), each = 8L)
  group_class <- rep(c(1L, 2L), length.out = 30L)
  prevalence <- rbind(c(0.85, 0.15), c(0.2, 0.8))
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = prevalence[group_class[j], ])
  }, integer(1))
  means <- rbind(c(-2, -2), c(2, 2))
  y <- t(vapply(profile, function(k) stats::rnorm(2L, means[k, ], 1), numeric(2)))
  dat <- data.frame(g = group, a = y[, 1L], b = y[, 2L])
  multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
}

test_that("a criterion that uses no sample size says NA, not a third level", {
  criteria <- information_criteria(.tidy_fit(), format = "long")
  scale_free <- subset(criteria, is.na(convention))
  expect_setequal(scale_free$criterion, c("deviance", "aic", "kic"))
  expect_true(all(is.na(scale_free$n)))
  # The regression this guards: "none" made `convention` a three-level
  # variable, so selecting one convention silently dropped aic and kic.
  expect_false(any(criteria$convention %in% "none"))
  expect_identical(nrow(subset(criteria, convention == "individuals")), 6L)
  expect_identical(nrow(subset(criteria, convention == "groups")), 6L)
})

test_that("the value column points one way on every row", {
  fit <- .tidy_fit()
  criteria <- information_criteria(fit, format = "long")
  # Every row is the same deviance plus a non-negative complexity charge, so a
  # smaller value is a better-supported model throughout. The log likelihood,
  # for which the opposite held, is logLik() and is not in the long table.
  expect_false("log_likelihood" %in% criteria$criterion)
  expect_equal(subset(criteria, criterion == "deviance")$value,
               -2 * as.numeric(logLik(fit)))
  expect_true(all(criteria$value >= -2 * as.numeric(logLik(fit)) - 1e-8))
  # The wide form reports the likelihood instead, because that is what a
  # model-comparison table publishes.
  expect_true("log_likelihood" %in% names(information_criteria(fit)))
})

test_that("formulas are documentation, available only on request", {
  fit <- .tidy_fit()
  expect_false("definition" %in% names(information_criteria(fit, format = "long")))
  explained <- information_criteria(fit, format = "long", definitions = TRUE)
  expect_identical(names(explained),
    c("criterion", "convention", "n", "value", "definition"))
  # The column is a pure function of `criterion`, which is why it is optional.
  by_criterion <- tapply(explained$definition, explained$criterion,
                         function(value) length(unique(value)))
  expect_true(all(by_criterion == 1L))
  expect_identical(nrow(explained), nrow(information_criteria(fit, format = "long")))
  expect_error(information_criteria(fit, format = "long", definitions = "yes"),
               "must be TRUE or FALSE")
  # A formula per row needs rows; the wide form has one, so the combination is
  # refused rather than quietly ignored.
  expect_error(information_criteria(fit, definitions = TRUE),
               class = "multilpa_bad_argument")
})

test_that("the wide and long forms are the same numbers", {
  fit <- .tidy_fit()
  wide <- information_criteria(fit)
  long <- information_criteria(fit, format = "long")
  expect_identical(nrow(wide), 1L)
  expect_equal(wide$log_likelihood, fit$log_likelihood)
  expect_equal(wide$bic_groups,
               subset(long, criterion == "bic" & convention == "groups")$value)
  expect_equal(wide$bic_individual,
               subset(long, criterion == "bic" & convention == "individuals")$value)
  expect_equal(wide$aic, subset(long, criterion == "aic")$value)
})

test_that("a single fit and its enumeration grid name the same things", {
  # Both go through one pivot, so this cannot drift. A criterion added to the
  # long table but not to the grid, or renamed in either, fails here.
  fit <- .tidy_fit()
  wide <- information_criteria(fit)
  set.seed(3)
  dat <- data.frame(g = rep(seq_len(20L), each = 6L),
                    a = stats::rnorm(120L), b = stats::rnorm(120L))
  grid <- as.data.frame(enumerate_classes(
    dat, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
    n_starts = 2, seed = 1))
  shared <- intersect(names(wide), names(grid))
  expect_setequal(shared, setdiff(names(wide), character()))
  expect_true(all(c("aic", "kic", "bic_groups", "bic_individual",
                    "log_likelihood", "n_parameters") %in% shared))
})

test_that("classification_table has one shape and average_posteriors the other", {
  fit <- .tidy_fit()
  summary_table <- classification_table(fit, level = "both")
  expect_identical(names(summary_table),
    c("level", "class", "n_modal", "proportion_modal", "estimated_n",
      "estimated_proportion", "average_posterior",
      "odds_correct_classification"))
  cross <- average_posteriors(fit, level = "both")
  expect_identical(names(cross),
    c("level", "assigned_class", "class", "n_assigned", "average_posterior"))
  # One row per level x assigned class x class.
  expect_identical(nrow(cross), 8L)
  # Each assigned class's row is a probability distribution.
  totals <- tapply(cross$average_posterior, interaction(cross$level,
                                                        cross$assigned_class), sum)
  expect_equal(unname(as.vector(totals)), rep(1, 4L))
  # The summary column is the diagonal of the cross-tabulation.
  diagonal <- subset(cross, assigned_class == class)
  expect_equal(diagonal$average_posterior, summary_table$average_posterior)
  # The removed flag is refused by class, not by an "unused argument" message.
  expect_error(classification_table(fit, detail = TRUE),
               class = "multilpa_removed_argument")
})

test_that("average_posteriors is not classification_errors under another name", {
  fit <- .tidy_fit()
  cross <- average_posteriors(fit, level = "individuals")
  errors <- classification_errors(fit, level = "individuals")
  # Same square, opposite conditioning: one averages posteriors within the
  # assigned class, the other distributes the true class over assignments.
  expect_identical(nrow(cross), nrow(errors))
  expect_false(isTRUE(all.equal(sort(cross$average_posterior),
                                sort(errors$probability))))
})

test_that("a class with no modal members reports undefined, not zero", {
  # A hand-built posterior whose second class is never modal.
  probabilities <- cbind(c(0.9, 0.8, 0.7), c(0.1, 0.2, 0.3))
  object <- structure(list(subject_posteriors = probabilities), class = "multilpa")
  table <- classification_table(object)
  empty <- subset(table, class == 2L)
  expect_identical(empty$n_modal, 0L)
  expect_true(is.na(empty$average_posterior))
  expect_true(is.na(empty$odds_correct_classification))
  # Its estimated size is still positive: the model gives it weight, modal
  # assignment does not.
  expect_gt(empty$estimated_n, 0)
  cross <- average_posteriors(object)
  expect_true(all(is.na(subset(cross, assigned_class == 2L)$average_posterior)))
})

test_that("Cramer's V matches an independent computation", {
  # Calibration against stats::chisq.test, not against our own arithmetic.
  counts <- matrix(c(30, 10, 5, 55), 2L, 2L)
  reference <- quietly(stats::chisq.test(counts, correct = FALSE))
  expect_equal(.multilpa_cramers_v(counts),
               sqrt(unname(reference$statistic) / sum(counts)))
  # For a two-by-two table Cramer's V is the absolute indicator correlation.
  long <- data.frame(a = rep(c(0, 0, 1, 1), c(30, 5, 10, 55)),
                     b = rep(c(0, 1, 0, 1), c(30, 5, 10, 55)))
  expect_equal(.multilpa_cramers_v(counts), abs(stats::cor(long$a, long$b)))
  three <- matrix(c(20, 5, 3, 4, 18, 25), 3L, 2L)
  reference_three <- quietly(stats::chisq.test(three, correct = FALSE))
  expect_equal(.multilpa_cramers_v(three),
               sqrt(unname(reference_three$statistic) / (sum(three) * 1L)))
  # A product of its own margins carries no association at all.
  expect_equal(.multilpa_cramers_v(outer(c(30, 70), c(40, 60)) / 100), 0)
  # Undefined rather than zero where there is no table to speak of.
  expect_true(is.na(.multilpa_cramers_v(matrix(0, 2L, 2L))))
  expect_true(is.na(.multilpa_cramers_v(matrix(c(3, 4), 2L, 1L))))
})

.mixed_fit_data <- function() {
  set.seed(21)
  n_groups <- 40L
  school <- rep(seq_len(n_groups), each = 12L)
  profile <- rep(rep(c(1L, 2L), length.out = n_groups), each = 12L)
  n <- length(school)
  shared_continuous <- stats::rnorm(n)
  shared_binary <- stats::runif(n) < 0.5
  chance <- function() stats::runif(n) < ifelse(profile == 2L, 0.8, 0.2)
  data.frame(school = school,
             a = stats::rnorm(n, ifelse(profile == 2L, 2, -2), 0.8) + shared_continuous,
             b = stats::rnorm(n, ifelse(profile == 2L, 2, -2), 0.8) + shared_continuous,
             u1 = as.integer(ifelse(shared_binary, TRUE, chance())),
             u2 = as.integer(ifelse(shared_binary, TRUE, chance())),
             u3 = as.integer(chance()))
}

test_that("a mixed fit puts both kinds of pair on one scale", {
  data <- .mixed_fit_data()
  fit <- multilpa(data, c("a", "b", "u1", "u2", "u3"), "school", n_profiles = 2,
                  n_group_classes = 1, categorical = c("u1", "u2", "u3"),
                  n_starts = 5, seed = 3)
  residuals <- bivariate_residuals(fit, data)
  expect_setequal(unique(residuals$kind), c("gaussian", "categorical"))
  # `residual` is observed minus expected for every kind, not only one.
  expect_equal(residuals$residual, residuals$observed - residuals$expected)
  # Both associations live on the correlation scale, so the ranking across
  # kinds compares like with like.
  expect_true(all(abs(residuals$observed) <= 1 + 1e-12))
  expect_true(all(abs(residuals$expected) <= 1 + 1e-12))
  expect_equal(abs(residuals$residual),
               sort(abs(residuals$residual), decreasing = TRUE))
  # Within a profile the measurement model is a product of margins, so the
  # categorical association it implies is zero, as the Gaussian one is.
  categorical <- subset(residuals, kind == "categorical")
  expect_equal(categorical$expected, rep(0, nrow(categorical)),
               tolerance = 1e-8)
  # Pooling mixes profiles, so the model then implies a nonzero association.
  pooled <- subset(bivariate_residuals(fit, data, by = "overall"),
                   kind == "categorical")
  expect_gt(max(pooled$expected), 0.05)
  # The planted dependence is still the worst pair, over both kinds.
  worst <- head(bivariate_residuals(fit, data, by = "overall"), 1L)
  expect_setequal(c(worst$indicator_1, worst$indicator_2), c("a", "b"))
})

test_that("df is a column that does not apply to a Gaussian pair", {
  data <- .mixed_fit_data()
  fit <- multilpa(data, c("a", "b", "u1", "u2", "u3"), "school", n_profiles = 2,
                  n_group_classes = 1, categorical = c("u1", "u2", "u3"),
                  n_starts = 5, seed = 3)
  residuals <- bivariate_residuals(fit, data, by = "overall")
  expect_type(residuals$df, "integer")
  # A Fisher z is a standard normal deviate and has no degrees of freedom.
  expect_true(all(is.na(subset(residuals, kind == "gaussian")$df)))
  # A chi-square does, and every categorical row reports it.
  expect_false(any(is.na(subset(residuals, kind == "categorical")$df)))
  # `statistic` and `p_value`, unlike `df`, apply to both kinds.
  expect_false(any(is.na(residuals$statistic)))
  expect_false(any(is.na(residuals$p_value)))
  # The empty table declares the same column types as a populated one.
  single <- multilpa(data, "a", "school", n_profiles = 2, n_group_classes = 1,
                     n_starts = 2, seed = 1)
  empty <- bivariate_residuals(single, data)
  expect_identical(names(empty), names(residuals))
  expect_type(empty$df, "integer")
})

test_that("a model with no discrete group classes refuses by condition class", {
  set.seed(174)
  data <- data.frame(g = rep(seq_len(20L), each = 5L),
                     y = rep(stats::rnorm(20L), each = 5L) + stats::rnorm(100L))
  fit <- fit_random_intercept(data, "y", "g", n_profiles = 1, n_starts = 1)
  expect_error(classification_table(fit, level = "groups"),
               class = "multilpa_no_group_classes")
  expect_error(average_posteriors(fit, level = "groups"),
               class = "multilpa_no_group_classes")
  expect_error(bivariate_residuals(fit, data),
               class = "multilpa_no_group_classes")
  # "both" degrades to the level the model actually has.
  expect_identical(unique(classification_table(fit, level = "both")$level),
                   "individuals")
  expect_identical(unique(average_posteriors(fit, level = "both")$level),
                   "individuals")
})
