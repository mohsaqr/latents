# Follow-up review of 2026-10-02: defects found by independent re-inspection
# of the 0.9.7 review (all present before that review).

# Regression tests from past audits; CI runs them on every platform.
skip_on_cran()

followup_lta_data <- function(n = 120L, occasions = 4L, seed = 41L) {
  set.seed(seed)
  transition <- matrix(c(0.85, 0.15, 0.15, 0.85), 2L)
  states <- matrix(NA_integer_, n, occasions)
  states[, 1L] <- sample.int(2L, n, replace = TRUE)
  invisible(lapply(seq_len(occasions)[-1L], function(t) {
    states[, t] <<- vapply(states[, t - 1L], function(from)
      sample.int(2L, 1L, prob = transition[from, ]), integer(1))
  }))
  s <- as.vector(t(states))
  y1 <- c(0, 3)[s] + stats::rnorm(length(s))
  data.frame(id = rep(seq_len(n), each = occasions),
             time = rep(seq_len(occasions), n),
             y1 = y1, y2 = c(0.5, 3.5)[s] + 0.6 * y1 + stats::rnorm(length(s)))
}

test_that("bootstrap_lrt simulates full-covariance transition nulls", {
  skip_on_cran()
  data <- followup_lta_data()
  null <- lta(data, c("y1", "y2"), "id", 2, time = "time", transitions = "occasion",
              model = "VVV", n_starts = 1, seed = 1)
  alternative <- lta(data, c("y1", "y2"), "id", 2, time = "time",
                     transitions = "occasion", order = 2, model = "VVV",
                     n_starts = 1, seed = 1)
  # 0.9.7 failed every replicate (a swallowed refusal) and returned p = NA;
  # 0.9.8 refused up front; full covariances are now simulated.
  result <- bootstrap_lrt(null, alternative, data = data, iter = 4,
                          n_starts = 1, seed = 1)
  test <- get_results(result, "test")
  expect_identical(test$n_valid, 4L)
  expect_true(is.finite(test$p_value))
  homogeneous <- lta(data, c("y1", "y2"), "id", 2, time = "time", model = "VVV",
                     n_starts = 1, seed = 1)
  expect_true(is.finite(get_results(bootstrap_lrt(homogeneous, null, data = data,
                                                  iter = 2, n_starts = 1, seed = 1),
                                    "test")$p_value))
})

test_that("transition parameter_inference handles every generic argument explicitly", {
  data <- followup_lta_data()
  homogeneous <- lta(data, c("y1", "y2"), "id", 2, time = "time", n_starts = 1, seed = 1)
  general <- lta(data, c("y1", "y2"), "id", 2, time = "time", transitions = "occasion",
                 n_starts = 1, seed = 1)
  expect_s3_class(homogeneous, "multilpa_transitions")
  expect_s3_class(general, "multilpa_lta")
  invisible(lapply(list(homogeneous, general), function(fit) {
    # The bootstrap needs the persons to resample: refused, not ignored.
    expect_error(parameter_inference(fit, method = "bootstrap", iter = 3, seed = 1),
                 class = "latents_bad_argument")
    # Nothing is held at a bound in transition inference: on an interior fit
    # "fix" and "error" give the same table.
    expect_identical(parameter_inference(fit, boundary = "fix"), parameter_inference(fit))
  }))
})

test_that("categorical transition response tables keep their standard errors", {
  data <- followup_lta_data()
  set.seed(4)
  data$c1 <- factor(ifelse(data$y1 + stats::rnorm(nrow(data)) > 1.5, "hi", "lo"))
  fit <- lta(data, c("y1", "c1"), "id", 2, time = "time", categorical = "c1",
             n_starts = 1, seed = 1)
  expect_true(isTRUE(fit$converged))
  # 0.9.8 refused boundary = "fix", so this table silently lost every error.
  responses <- get_results(fit, "responses")
  expect_true(all(is.finite(responses$probability_standard_error)))
  wald <- parameter_inference(fit)
  expect_equal(sort(unique(round(responses$probability_standard_error, 10))),
               sort(unique(round(subset(wald, parameter == "response")$standard_error, 10))))
})

test_that("general transition inference applies the requested multiplicity correction", {
  data <- followup_lta_data()
  fit <- lta(data, c("y1", "y2"), "id", 2, time = "time", transitions = "occasion",
             n_starts = 1, seed = 1)
  plain <- parameter_inference(fit)
  expect_identical(plain$p_adjusted, plain$p_value)
  expect_identical(attr(plain, "adjust"), "none")
  corrected <- parameter_inference(fit, adjust = "holm")
  expect_identical(attr(corrected, "adjust"), "holm")
  # The family is the tests the table carries, computed independently here.
  tested <- !is.na(corrected$p_value)
  expect_gt(sum(tested), 1L)
  reference <- stats::p.adjust(corrected$p_value[tested], method = "holm")
  expect_equal(corrected$p_adjusted[tested], reference, tolerance = 1e-14)
  expect_true(all(corrected$p_adjusted[tested] >= corrected$p_value[tested]))
  expect_identical(match("p_adjusted", names(corrected)),
                   match("p_value", names(corrected)) + 1L)
  # Every other column is unchanged by the correction.
  expect_identical(corrected[setdiff(names(corrected), "p_adjusted")],
                   plain[setdiff(names(plain), "p_adjusted")])
})

test_that("near-Poisson negative-binomial dispersion reaches the profile maximum", {
  set.seed(2)
  x1 <- stats::rnorm(2000L)
  y <- stats::rnbinom(2000L, size = 1e9, mu = 500)
  data <- data.frame(x1 = x1, c1 = y)
  fit <- lpa(data, c("x1", "c1"), 1, count = "c1",
             count_model = "negative_binomial", n_starts = 1)
  # One profile: the NB2 mean MLE is the sample mean, so the dispersion MLE is
  # the maximum of the profile likelihood in alpha alone.
  profile <- function(log_alpha) {
    sum(stats::dnbinom(y, size = exp(-log_alpha), mu = mean(y), log = TRUE))
  }
  reference <- stats::optimize(profile, c(log(1e-9), 0), maximum = TRUE, tol = 1e-12)
  expect_true(isTRUE(fit$converged))
  # Before the fix the Newton step stalled at 4.2e-7 against 1.4e-5.
  fitted <- log(as.vector(fit$count_dispersion))
  expect_equal(fitted, reference$maximum, tolerance = 1e-4)
  expect_gte(profile(fitted), reference$objective - 1e-8)
  expect_true(all(is.finite(diag(vcov(fit)))))
})

test_that("the negative-binomial M-step ascends from either side of an indefinite region", {
  set.seed(2)
  stats::rnorm(2000L)
  y <- stats::rnbinom(2000L, size = 1e9, mu = 500)
  posteriors <- matrix(1, length(y), 1L)
  reference <- stats::optimize(function(log_alpha) {
    sum(stats::dnbinom(y, size = exp(-log_alpha), mu = mean(y), log = TRUE))
  }, c(log(1e-9), 0), maximum = TRUE, tol = 1e-12)$maximum
  # Starts below the maximum (where the objective is convex in log alpha) and
  # above it (the moment start, whose first full step overshot) both arrive.
  invisible(lapply(c(1e-8, 1e-7, 0.05, 1), function(start) {
    point <- .latents_negative_binomial_maximize(y, posteriors, mean(y) * 1.01, start)
    expect_equal(log(point$dispersion), reference, tolerance = 1e-4)
  }))
})

test_that("bootstrap_lrt accepts the caller's frame for single-level fits", {
  set.seed(3)
  data <- data.frame(x1 = c(stats::rnorm(60), stats::rnorm(60, 3)),
                     x2 = c(stats::rnorm(60), stats::rnorm(60, 2)))
  one <- lpa(data, c("x1", "x2"), 1, n_starts = 1)
  two <- lpa(data, c("x1", "x2"), 2, n_starts = 1, seed = 1)
  stored <- bootstrap_lrt(one, two, iter = 2, n_starts = 1, seed = 5)
  supplied <- bootstrap_lrt(one, two, data = data, iter = 2, n_starts = 1, seed = 5)
  expect_equal(supplied$statistic, stored$statistic, tolerance = 1e-12)
  expect_equal(supplied$replicates$statistic, stored$replicates$statistic,
               tolerance = 1e-10)
  # Reordered rows still do not reproduce the fit.
  expect_error(bootstrap_lrt(one, two, data = data[rev(seq_len(nrow(data))), ],
                             iter = 2, n_starts = 1, seed = 5),
               class = "latents_bad_inference_data")
})

test_that("weighted classification tables use weighted class totals", {
  set.seed(1)
  n <- 120L
  data <- data.frame(x1 = stats::rnorm(n, rep(c(0, 3), each = 60L)),
                     x2 = stats::rnorm(n, rep(c(0, 2), each = 60L)))
  data$w <- sample(1:4, n, replace = TRUE)
  replicated <- data[rep(seq_len(n), data$w), ]
  weighted <- lpa(data, c("x1", "x2"), 2, weights = "w", seed = 1, n_starts = 1,
                  tol = 1e-12)
  frequency <- lpa(replicated, c("x1", "x2"), 2, seed = 1, n_starts = 1, tol = 1e-12)
  table_w <- get_results(weighted, "classification", level = "individuals")
  table_f <- get_results(frequency, "classification", level = "individuals")
  # Integer weights are frequency weights: shares and OCC match replication.
  expect_equal(table_w$estimated_proportion, table_f$estimated_proportion,
               tolerance = 1e-6)
  # The odds of correct classification compare the (per-unit) average
  # posterior with the model's weighted share, as documented.
  odds <- function(p) p / (1 - p)
  expect_equal(table_w$odds_correct_classification,
               odds(table_w$average_posterior) / odds(table_w$estimated_proportion),
               tolerance = 1e-12)
  expect_equal(table_w$estimated_proportion, as.vector(weighted$profile_probabilities),
               tolerance = 1e-8)
  expect_equal(table_w$estimated_n, unname(weighted$effective_profile_counts),
               tolerance = 1e-12)
  expect_equal(sum(table_w$estimated_proportion), 1, tolerance = 1e-12)
  # Modal counts stay per-unit descriptions of the rows supplied.
  expect_identical(sum(table_w$n_modal), n)
  unweighted <- lpa(data, c("x1", "x2"), 2, seed = 1, n_starts = 1, tol = 1e-12)
  expect_equal(get_results(unweighted, "classification", level = "individuals")$estimated_n,
               unname(colSums(unweighted$subject_posteriors)), tolerance = 1e-12)
})

test_that("restricted cross-level composition and class counts use group weights", {
  set.seed(61)
  n_groups <- 40L
  size <- 6L
  group <- rep(seq_len(n_groups), each = size)
  group_class <- sample(1:2, n_groups, TRUE)
  profile <- ifelse(stats::runif(n_groups * size) <
                      c(0.8, 0.2)[group_class[group]], 1, 2)
  group_mean <- stats::rnorm(n_groups, c(0, 1.5)[group_class], 0.4)
  data <- data.frame(id = group,
                     x1 = stats::rnorm(n_groups * size, c(0, 3)[profile] + group_mean[group]),
                     x2 = stats::rnorm(n_groups * size, c(0, 2)[profile] + group_mean[group]))
  group_weight <- sample(1:3, n_groups, TRUE)
  data$w <- group_weight[data$id]
  copies <- rep(seq_len(n_groups), group_weight)
  replicated <- do.call(rbind, lapply(seq_along(copies), function(j) {
    rows <- data[data$id == copies[j], ]
    rows$id <- j
    rows
  }))
  fit <- function(frame, weights) {
    multilpa(frame, c("x1", "x2"), "id", 2, 2, family = "restricted_cross_level",
             weights = weights, seed = 4, n_starts = 1, tol = 1e-14, max_iter = 20000)
  }
  weighted <- fit(data, "w")
  frequency <- fit(replicated, NULL)
  scale <- n_groups / sum(group_weight)
  classes_w <- get_results(weighted, "group_classes")
  classes_f <- get_results(frequency, "group_classes")
  expect_equal(classes_w$count, classes_f$count * scale, tolerance = 1e-5)
  expect_equal(classes_w$count, classes_w$weight * n_groups, tolerance = 1e-10)
  expect_equal(get_results(weighted, "composition")$share,
               get_results(frequency, "composition")$share, tolerance = 1e-5)
})
