# Review 2026-11-01, F7 and F10: three_step() with a categorical distal
# outcome.
# F7: a BCH probability outside [0, 1] is reported without an interval or
#     test, and so is every pairwise difference that involves it; the other
#     differences keep their inference.
# F10: the analysed categories are those of the units that enter step three,
#      so a category carried only by an excluded row or group changes nothing.

skip_on_cran()

# The two-class fixture of test-three-step-categorical.R in which both
# probabilities of class 1 fall outside [0, 1].
out_of_range_case <- function() {
  set.seed(10)
  truth <- sample.int(2L, 300L, replace = TRUE, prob = c(0.8, 0.2))
  data <- data.frame(y = stats::rnorm(300L, c(0, 1.2)[truth]),
                     z = stats::rnorm(300L, c(0, 1.2)[truth]))
  data$outcome <- factor(ifelse(stats::runif(300L) < c(0.02, 0.3)[truth], "rare", "common"))
  fit <- quietly(lpa(data, c("y", "z"), n_profiles = 2, n_starts = 2, seed = 1),
                 "latents_boundary")
  list(data = data, fit = fit)
}

# Three classes, three categories: only class 1's share of "c" falls below 0
# (its true share is 0), so some differences involve it and others do not.
partial_case <- function() {
  set.seed(2)
  truth <- sample.int(3L, 400L, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  centres <- rbind(c(0, 0), c(1.3, 1.3), c(3.5, -1))
  data <- data.frame(y = stats::rnorm(400L, centres[truth, 1L]),
                     z = stats::rnorm(400L, centres[truth, 2L]))
  shares <- rbind(c(0.97, 0.03, 0), c(0.5, 0.3, 0.2), c(0.3, 0.3, 0.4))
  data$outcome <- factor(vapply(truth, function(class) {
    sample(c("a", "b", "c"), 1L, prob = shares[class, ])
  }, character(1)), levels = c("a", "b", "c"))
  fit <- lpa(data, c("y", "z"), n_profiles = 3, model = "VVI", n_starts = 2, seed = 1)
  list(data = data, fit = fit)
}

# Inference of the differences that do not involve an out-of-range cell is
# the plain Wald inference of that difference; the adjustment counts only them.
expect_unaffected_inference <- function(pairs, involved, ci_level = 0.95) {
  kept <- subset(pairs, !involved)
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  expect_true(all(is.finite(kept$statistic) & is.finite(kept$p_value) &
                    is.finite(kept$p_adjusted) & is.finite(kept$conf_low) &
                    is.finite(kept$conf_high)))
  expect_equal(kept$statistic, kept$estimate / kept$standard_error, tolerance = 1e-12)
  expect_equal(kept$p_value, 2 * stats::pnorm(-abs(kept$statistic)), tolerance = 1e-12)
  expect_equal(kept$conf_low, kept$estimate - quantile * kept$standard_error,
               tolerance = 1e-12)
  expect_equal(kept$conf_high, kept$estimate + quantile * kept$standard_error,
               tolerance = 1e-12)
  expect_equal(kept$p_adjusted, stats::p.adjust(kept$p_value, method = attr(pairs, "adjust")),
               tolerance = 1e-12)
}

test_that("F7: pairwise differences of out-of-range BCH probabilities have no test", {
  case <- out_of_range_case()
  expect_warning(individual <- three_step(case$fit, case$data, "outcome"),
                 class = "latents_bch_out_of_range")
  expect_true(all(subset(individual, class == 1L)$estimate < 0 |
                    subset(individual, class == 1L)$estimate > 1))
  expect_warning(pairs <- three_step(case$fit, case$data, "outcome", contrast = "pairs"),
                 class = "latents_bch_out_of_range")
  expect_identical(nrow(pairs), 2L)
  # Every difference involves class 1: estimate and standard error stay, the
  # statistic, p-values and interval do not.
  expect_true(all(is.finite(pairs$estimate) & is.finite(pairs$standard_error)))
  expect_true(all(is.na(pairs$statistic) & is.na(pairs$p_value) & is.na(pairs$p_adjusted) &
                    is.na(pairs$conf_low) & is.na(pairs$conf_high)))
  expect_equal(pairs$estimate,
               subset(individual, class == 2L)$estimate -
                 subset(individual, class == 1L)$estimate, tolerance = 1e-12)
  expect_warning(estimated <- three_step(case$fit, case$data, "outcome", contrast = "pairs",
                                         step_one = "estimated"),
                 class = "latents_bch_out_of_range")
  expect_true(all(is.na(estimated$statistic) & is.na(estimated$p_value) &
                    is.na(estimated$conf_low) & is.na(estimated$conf_high)))
  expect_true(all(is.finite(estimated$standard_error)))
})

test_that("F7: only the differences involving an out-of-range cell lose inference", {
  case <- partial_case()
  expect_warning(individual <- three_step(case$fit, case$data, "outcome"),
                 class = "latents_bch_out_of_range")
  invalid <- subset(individual, estimate < 0 | estimate > 1)
  expect_identical(nrow(invalid), 1L)
  # Class 1's other categories are in range, so the cell is partly affected.
  expect_identical(sum(individual$class == invalid$class), 3L)
  is_invalid <- function(class, category) {
    class == invalid$class & category == invalid$category
  }
  lapply(c("fixed", "estimated"), function(step_one) {
    expect_warning(pairs <- three_step(case$fit, case$data, "outcome", contrast = "pairs",
                                       step_one = step_one),
                   class = "latents_bch_out_of_range")
    expect_identical(nrow(pairs), 9L)
    involved <- is_invalid(pairs$class, pairs$category) |
      is_invalid(pairs$reference_class, pairs$category)
    # Two pairs touch class 1 in the out-of-range category; the remaining
    # seven, including class 1's other categories, keep their inference.
    expect_identical(sum(involved), 2L)
    affected <- subset(pairs, involved)
    expect_true(all(is.finite(affected$estimate) & is.finite(affected$standard_error)))
    expect_true(all(is.na(affected$statistic) & is.na(affected$p_value) &
                      is.na(affected$p_adjusted) & is.na(affected$conf_low) &
                      is.na(affected$conf_high)))
    expect_unaffected_inference(pairs, involved)
  })
  # The fixed and estimated step one differ only in the standard errors.
  fixed <- quietly(three_step(case$fit, case$data, "outcome", contrast = "pairs"),
                   "latents_bch_out_of_range")
  estimated <- quietly(three_step(case$fit, case$data, "outcome", contrast = "pairs",
                                  step_one = "estimated"), "latents_bch_out_of_range")
  expect_equal(fixed$estimate, estimated$estimate, tolerance = 1e-12)
  # The individual table and the overall test keep their policy.
  expect_true(is.na(invalid$conf_low) && is.na(invalid$conf_high))
  overall <- quietly(three_step(case$fit, case$data, "outcome", contrast = "overall"),
                     "latents_bch_out_of_range")
  expect_true(is.na(overall$statistic))
})

test_that("F7: in-range probabilities keep all their pairwise inference", {
  case <- partial_case()
  pairs <- three_step(case$fit, case$data, "outcome", method = "proportional",
                      contrast = "pairs")
  expect_unaffected_inference(pairs, rep(FALSE, nrow(pairs)))
})

# A binary outcome on the Old Faithful data, then the same data with one
# appended row whose indicators are all missing (left out of step three by the
# default classification policy) and whose outcome is a category no retained
# unit has. `relabelled` is the extended data with that row given a retained
# category instead: on the same fit, it is the result the unseen category
# must not change.
binary_case <- function() {
  data <- datasets::faithful
  data$outcome <- factor(rep(c("yes", "no"), length.out = nrow(data)))
  fit <- lpa(data, c("eruptions", "waiting"), n_profiles = 2, n_starts = 2, seed = 1)
  extended <- rbind(transform(data, outcome = as.character(outcome)),
                    data.frame(eruptions = NA, waiting = NA, outcome = "excluded"))
  relabelled <- transform(extended, outcome = factor(ifelse(outcome == "excluded", "yes",
                                                            outcome)))
  extended$outcome <- factor(extended$outcome)
  extended_fit <- lpa(extended, c("eruptions", "waiting"), n_profiles = 2, n_starts = 2,
                      seed = 1, missing = "fiml")
  list(data = data, fit = fit, extended = extended, relabelled = relabelled,
       extended_fit = extended_fit)
}

expect_same_three_step <- function(before, after, tolerance = 1e-10) {
  expect_identical(names(after), names(before))
  expect_identical(nrow(after), nrow(before))
  numeric_columns <- names(before)[vapply(before, is.numeric, logical(1))]
  lapply(numeric_columns, function(column) {
    expect_equal(after[[column]], before[[column]], tolerance = tolerance, label = column)
  })
  if ("category" %in% names(before)) expect_identical(after$category, before$category)
}

test_that("F10: a category seen only in an excluded row is not analysed", {
  case <- binary_case()
  expect_false(.multilpa_informative_rows(case$extended_fit)[nrow(case$extended)])
  grid <- expand.grid(method = c("bch", "proportional", "modal", "ml"),
                      contrast = c("none", "pairs", "overall"),
                      step_one = c("fixed", "estimated"), stringsAsFactors = FALSE)
  Map(function(method, contrast, step_one) {
    after <- three_step(case$extended_fit, case$extended, "outcome", method = method,
                        contrast = contrast, step_one = step_one)
    expect_false("excluded" %in% after$category)
    expect_true(all(is.finite(after$standard_error)))
    expect_same_three_step(
      three_step(case$extended_fit, case$relabelled, "outcome", method = method,
                 contrast = contrast, step_one = step_one),
      after)
    # Against the fit without the row. The row's empty contribution makes EM
    # stop one iteration earlier (posteriors differ by about 2.5e-9), so the
    # two fits, not the coding, set this tolerance.
    expect_same_three_step(
      three_step(case$fit, case$data, "outcome", method = method, contrast = contrast,
                 step_one = step_one),
      after, tolerance = 1e-6)
  }, grid$method, grid$contrast, grid$step_one)
  ml <- three_step(case$extended_fit, case$extended, "outcome", method = "ml")
  expect_identical(ml$category, rep(c("no", "yes"), 2L))
  expect_true(all(is.finite(ml$standard_error) & is.finite(ml$conf_low)))
})

test_that("F10: categories keep their factor order among the analysed units", {
  case <- binary_case()
  # "excluded" placed between the retained levels, and first.
  middle <- transform(case$extended,
                      outcome = factor(outcome, levels = c("yes", "excluded", "no")))
  first <- transform(case$extended,
                     outcome = factor(outcome, levels = c("excluded", "no", "yes")))
  reference <- transform(case$relabelled, outcome = factor(outcome, levels = c("yes", "no")))
  after <- three_step(case$extended_fit, middle, "outcome", method = "ml")
  expect_identical(after$category, rep(c("yes", "no"), 2L))
  expect_same_three_step(three_step(case$extended_fit, reference, "outcome", method = "ml"),
                         after)
  expect_same_three_step(three_step(case$extended_fit, case$relabelled, "outcome"),
                         three_step(case$extended_fit, first, "outcome"))
  # A character outcome is coded the same way.
  character_outcome <- transform(case$extended, outcome = as.character(outcome))
  expect_same_three_step(
    three_step(case$extended_fit, case$relabelled, "outcome", method = "ml"),
    three_step(case$extended_fit, character_outcome, "outcome", method = "ml"))
})

test_that("F10: one category among the analysed units is refused", {
  case <- binary_case()
  # Every retained unit says "yes"; only the excluded row says otherwise.
  single <- transform(case$extended,
                      outcome = factor(ifelse(is.na(eruptions), "excluded", "yes")))
  expect_error(three_step(case$extended_fit, single, "outcome"), class = "latents_bad_outcome")
  expect_error(three_step(case$extended_fit, single, "outcome", method = "ml"),
               class = "latents_bad_outcome")
})

# Forty groups of ten in two group classes, then the same with a 41st group
# whose rows have nothing observed, so the whole group is left out.
grouped_case <- function() {
  set.seed(11)
  n_groups <- 40L
  group_class <- 1L + as.integer(stats::runif(n_groups) > 0.5)
  group <- rep(seq_len(n_groups), each = 10L)
  truth <- 1L + as.integer(stats::runif(400L) < c(0.15, 0.85)[group_class[group]])
  data <- data.frame(g = group,
                     a = stats::rnorm(400L, ifelse(truth == 2L, 1.5, -1.5)),
                     b = stats::rnorm(400L, ifelse(truth == 2L, 1.5, -1.5)))
  pass <- stats::runif(n_groups) < c(0.3, 0.7)[group_class]
  data$outcome <- factor(ifelse(pass, "pass", "fail")[group], levels = c("fail", "pass"))
  added <- data.frame(g = n_groups + 1L, a = NA_real_, b = NA_real_,
                      outcome = factor("withdrawn", levels = "withdrawn"))
  extended <- rbind(transform(data, outcome = as.character(outcome)),
                    transform(added, outcome = as.character(outcome)))
  extended$outcome <- factor(extended$outcome, levels = c("fail", "withdrawn", "pass"))
  fit_of <- function(data, ...) {
    multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
             n_starts = 4, seed = 1, ...)
  }
  list(data = data, fit = fit_of(data), extended = extended,
       extended_fit = fit_of(extended, missing = "fiml"))
}

test_that("F10: a wholly excluded group carrying an unseen category changes nothing", {
  case <- grouped_case()
  expect_false(.multilpa_informative_groups(case$extended_fit)[41L])
  # The same fit with the excluded group's outcome relabelled to a retained
  # category isolates the coding of the outcome from the fit.
  relabelled <- transform(case$extended,
                          outcome = factor(ifelse(outcome == "withdrawn", "fail",
                                                  as.character(outcome)),
                                           levels = c("fail", "pass")))
  grid <- expand.grid(method = c("bch", "proportional", "modal", "ml"),
                      level = c("groups", "individuals"),
                      contrast = c("none", "pairs"), stringsAsFactors = FALSE)
  Map(function(method, level, contrast) {
    after <- three_step(case$extended_fit, case$extended, "outcome", level = level,
                        method = method, contrast = contrast)
    expect_false("withdrawn" %in% after$category)
    expect_same_three_step(
      three_step(case$extended_fit, relabelled, "outcome", level = level, method = method,
                 contrast = contrast),
      after)
  }, grid$method, grid$level, grid$contrast)
  ml <- three_step(case$extended_fit, case$extended, "outcome", level = "groups",
                   method = "ml")
  expect_true(all(is.finite(ml$standard_error)))
})
