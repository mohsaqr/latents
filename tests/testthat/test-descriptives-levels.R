# The within-group and between-group standard deviations of descriptives().

test_that("sd_within and sd_between equal the hand computation", {
  data <- data.frame(person = c(1, 1, 1, 2, 2, 3, 3, 3, 3),
                     score = c(2, 4, 6, 5, 7, 1, 1, 2, 4))
  described <- descriptives(data, vars = "score", id = "person")
  means <- c(4, 6, 2)
  # Pooled within: sum of squared deviations from each person's own mean,
  # over n - k = 9 - 3 degrees of freedom.
  within <- sqrt((4 + 0 + 4 + 1 + 1 + 1 + 1 + 0 + 4) / 6)
  expect_equal(described$sd_within, within, tolerance = 1e-12)
  expect_equal(described$sd_between, stats::sd(means), tolerance = 1e-12)
})

test_that("the level SDs agree with an analysis of variance and the group means", {
  set.seed(11)
  data <- data.frame(id = rep(1:40, times = sample(2:9, 40, replace = TRUE)))
  data$y <- stats::rnorm(nrow(data), mean = stats::rnorm(40)[data$id])
  data$y[c(3, 50, 77)] <- NA
  described <- descriptives(data, vars = "y", id = "id")
  observed <- subset(data, !is.na(y))
  anova_table <- stats::anova(stats::lm(y ~ factor(id), data = observed))
  expect_equal(described$sd_within, sqrt(anova_table[["Mean Sq"]][2]),
               tolerance = 1e-12)
  person_means <- aggregate(y ~ id, data = observed, FUN = mean)
  expect_equal(described$sd_between, stats::sd(person_means$y), tolerance = 1e-12)
})

test_that("shifting and scaling a variable shifts nothing and scales both SDs", {
  data <- data.frame(id = rep(1:6, each = 4), y = c(1:24) %% 7 + rep(1:6, each = 4))
  base <- descriptives(data, vars = "y", id = "id")
  data$y <- 10 + 3 * data$y
  moved <- descriptives(data, vars = "y", id = "id")
  expect_equal(moved$sd_within, 3 * base$sd_within, tolerance = 1e-12)
  expect_equal(moved$sd_between, 3 * base$sd_between, tolerance = 1e-12)
  expect_equal(moved$icc, base$icc, tolerance = 1e-12)
})

test_that("undefined level SDs are NA, not errors", {
  one_each <- data.frame(id = 1:5, y = c(2, 4, 1, 5, 3))
  described <- descriptives(one_each, vars = "y", id = "id")
  expect_true(is.na(described$sd_within))
  expect_equal(described$sd_between, stats::sd(c(2, 4, 1, 5, 3)))
  one_group <- data.frame(id = 1, y = c(2, 4, 6))
  described <- descriptives(one_group, vars = "y", id = "id")
  expect_true(is.na(described$sd_between))
  expect_equal(described$sd_within, 2)
  text <- data.frame(id = rep(1:3, each = 2), y = letters[1:6])
  described <- descriptives(text, vars = "y", id = "id")
  expect_true(is.na(described$sd_within) && is.na(described$sd_between))
})

test_that("the level SDs are reported per stratum with by", {
  data <- data.frame(id = rep(1:8, each = 3), y = c(1:24) %% 5,
                     cohort = rep(c("a", "b"), each = 12))
  by_cohort <- descriptives(data, vars = "y", id = "id", by = "cohort")
  alone <- descriptives(subset(data, cohort == "a"), vars = "y", id = "id")
  expect_equal(subset(by_cohort, cohort == "a")$sd_within, alone$sd_within)
  expect_equal(subset(by_cohort, cohort == "a")$sd_between, alone$sd_between)
})
