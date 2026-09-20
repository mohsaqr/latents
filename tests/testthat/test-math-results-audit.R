test_that("missing-row information criteria agree across all result accessors", {
  set.seed(183)
  d <- data.frame(g = rep(1:8, each = 6), y = rnorm(48))
  d$y[1:3] <- NA_real_
  fit <- multilpa(d, "y", "g", 1, 1, missing = "fiml", n_starts = 1)
  criteria <- information_criteria(fit)
  individual <- subset(criteria, convention == "individuals")
  expect_equal(individual$n, rep(45, 6))
  # One normal population: independently calculate its maximized likelihood.
  y <- na.omit(d$y)
  ll <- sum(dnorm(y, mean(y), sqrt(mean((y - mean(y))^2)), log = TRUE))
  expect_equal(fit$log_likelihood, ll, tolerance = 1e-8)
  expect_equal(individual$value[individual$criterion == "bic"],
               -2 * ll + 2 * log(length(y)), tolerance = 1e-8)
  expect_equal(individual$value[individual$criterion == "bic"], fit$bic_individual)
  grid <- enumerate_classes(d, "y", "g", n_profiles = 1, n_group_classes = 1,
                            missing = "fiml", n_starts = 1)
  expect_equal(grid$table$bic_individual, grid$fits[[1]]$bic_individual)
  expect_equal(as.data.frame(fit, what = "information_criteria"), criteria)
})

test_that("dropping a univariate covariance preserves one variance per profile", {
  start <- list(means = matrix(c(-1, 1), 2, 1),
                covariances = array(c(2, 3), c(1, 1, 2)),
                profile_probabilities = matrix(c(.4, .6), 1, 2),
                group_probabilities = 1)
  diagonal <- starting_values(start, covariance = "drop")
  expect_equal(diagonal$variances, matrix(c(2, 3), 2, 1))
  d <- data.frame(g = rep(1:3, each = 2), y = c(-2, -1, 0, 1, 2, 3))
  fit <- multilpa(d, "y", "g", 2, 1, start = diagonal, max_iter = 0, n_starts = 1)
  expected <- sum(log(.4 * dnorm(d$y, -1, sqrt(2)) +
                     .6 * dnorm(d$y, 1, sqrt(3))))
  expect_equal(fit$log_likelihood, expected, tolerance = 1e-12)
})

test_that("random-intercept diagnostics distinguish profiles from group classes", {
  set.seed(174)
  d <- data.frame(g = rep(1:20, each = 5),
                  y = rep(rnorm(20), each = 5) + rnorm(100))
  fit <- fit_random_intercept(d, "y", "g", n_profiles = 1, n_starts = 1)
  indices <- information_criteria(fit)
  expect_equal(subset(indices, criterion == "bic" & convention == "groups")$value,
               fit$bic)
  expect_equal(subset(indices, criterion == "bic" & convention == "individuals")$value,
               fit$bic_individual)
  expect_true(all(is.na(subset(indices, convention == "groups" &
                                criterion %in% c("awe", "icl"))$value)))
  expect_equal(classification_table(fit, level = "both")$level, "individuals")
  expect_error(classification_table(fit, level = "groups"), "no discrete group")
  expect_equal(entropy_table(fit)$level, "individuals")
  expect_equal(entropy_table(fit)$entropy_sum, 0)
})

test_that("categorical bootstraps preserve probability constraints and data identity", {
  set.seed(113)
  d <- data.frame(g = rep(1:20, each = 10),
                  a = sample(c("a", "b"), 200, TRUE),
                  b = sample(c("c", "d"), 200, TRUE))
  small <- multilpa(d, c("a", "b"), "g", 1, 1, categorical = c("a", "b"),
                    min_probability = .1, n_starts = 1)
  large <- multilpa(d, c("a", "b"), "g", 2, 1, categorical = c("a", "b"),
                    min_probability = .1, n_starts = 2, seed = 1)
  changed <- d
  changed$a <- ifelse(d$a == "a", "renamed-a", "renamed-b")
  expect_error(bootstrap_lrt(small, large, changed, iter = 2), "categorical levels")
  incompatible <- large
  incompatible$min_probability <- .05
  expect_error(bootstrap_lrt(small, incompatible, d, iter = 2), "same observations")
  bounds <- numeric()
  real_fit <- multilpa
  local_mocked_bindings(multilpa = function(..., min_probability) {
    bounds <<- c(bounds, min_probability)
    real_fit(..., min_probability = min_probability)
  })
  result <- suppressWarnings(bootstrap_lrt(small, large, d, iter = 2,
                                           n_starts = 1, max_iter = 2000, seed = 5))
  expect_equal(bounds, rep(.1, 4))
  expect_equal(nrow(result$replicates), 2L)
})

test_that("categorical bootstrap refuses incomplete observations", {
  d <- data.frame(g = rep(1:10, each = 8), a = rep(c(1, 2, 2, 1), 20),
                  b = rep(c(1, 1, 2, 2), 20))
  d$a[1] <- NA
  small <- multilpa(d, c("a", "b"), "g", 1, 1, categorical = c("a", "b"),
                    missing = "fiml", n_starts = 1)
  large <- multilpa(d, c("a", "b"), "g", 2, 1, categorical = c("a", "b"),
                    missing = "fiml", n_starts = 1, seed = 1)
  expect_error(bootstrap_lrt(small, large, d, iter = 2), "complete finite")
})

test_that("sequence totals count observed groups rather than unused factor levels", {
  d <- data.frame(g = factor(rep(c("b", "a"), each = 3),
                             levels = c("a", "b", "unused")),
                  time = rep(1:3, 2), y = c(-2, -1, -3, 2, 1, 3))
  fit <- multilpa(d, "y", "g", 1, 1, n_starts = 1, time = "time")
  result <- sequence_summary(fit)
  expect_equal(result$groups, 2L)
  expect_equal(result$observations, 6L)
  expect_equal(result$complete, 2L)
  expect_equal(result$mean_length, 3)
  expect_false(anyNA(result))
})

test_that("sequence outputs keep distinct numeric groups with identical printed labels", {
  d <- data.frame(g = rep(c(1, 1 + 1e-15), each = 3),
                  time = rep(1:3, 2), y = c(-2, -1, -3, 2, 1, 3))
  fit <- multilpa(d, "y", "g", 1, 1, n_starts = 1)
  # Attach the already-valid ordering to isolate output from input validation.
  fit$time_values <- d$time
  fit$time <- "time"
  wide <- sequences(fit, "wide")
  expect_equal(dim(wide), c(2L, 5L))
  # The point of this test: two groups that PRINT identically must stay distinct.
  # They now survive as native doubles in a real column, which is stronger than the
  # make.unique()'d row names this replaced - those kept them apart only as strings.
  expect_identical(wide$group, unique(d$g))
  expect_equal(sequence_summary(fit)$groups, 2L)
  expect_equal(sequence_summary(fit)$observations, 6L)
})
