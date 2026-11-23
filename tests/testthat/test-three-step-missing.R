skip_on_cran()

# A distal outcome missing for some units: `missing = "exclude"` leaves those
# units out of step three, as Latent GOLD's `missing excludeall` does.

.missing_data <- function(seed = 11L, n_groups = 60L, per = 8L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per)
  n <- length(group)
  person_class <- rep(1L + as.integer(stats::runif(n_groups) > 0.5), each = per)
  truth <- ifelse(stats::runif(n) < ifelse(person_class == 2L, 0.85, 0.15), 2L, 1L)
  data.frame(g = group,
             a = stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)),
             b = stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)),
             y = stats::rnorm(n, ifelse(truth == 2L, 5, 0)),
             person_y = rep(stats::rnorm(n_groups, 0, 1), each = per) +
               ifelse(person_class == 2L, 3, 0))
}

.missing_fit <- function(data, n_group_classes = 1) {
  multilpa(data, c("a", "b"), "g", n_profiles = 2,
           n_group_classes = n_group_classes, n_starts = 6, seed = 1)
}

# BCH by its definition (Bolck, Croon & Hagenaars, 2004; Vermunt, 2010), from
# the posteriors and modal classes of the units analysed.
.bch_by_hand <- function(posteriors, modal, y) {
  classes <- seq_len(ncol(posteriors))
  errors <- vapply(classes, function(assigned) {
    colSums(posteriors * (modal == assigned)) / colSums(posteriors)
  }, numeric(length(classes)))
  weights <- solve(errors)[modal, , drop = FALSE]
  unname(as.vector(crossprod(weights, y)) / colSums(weights))
}

test_that("the default still refuses a missing outcome", {
  data <- .missing_data()
  fit <- .missing_fit(data)
  data$y[c(3, 40)] <- NA
  expect_error(three_step(fit, data, "y"), class = "latents_bad_argument")
  expect_error(three_step(fit, data, "y", missing = "drop"),
               class = "latents_bad_argument")
})

test_that("an outcome missing everywhere is refused even when excluding", {
  data <- .missing_data()
  fit <- .missing_fit(data)
  data$y <- NA_real_
  expect_error(three_step(fit, data, "y", missing = "exclude"),
               class = "latents_bad_argument")
})

test_that("excluding changes nothing when nothing is missing", {
  data <- .missing_data()
  fit <- .missing_fit(data)
  methods <- c("bch", "proportional", "modal", "ml")
  invisible(lapply(methods, function(method) {
    expect_identical(
      three_step(fit, data, "y", method = method, missing = "exclude"),
      three_step(fit, data, "y", method = method))
  }))
})

test_that("excluded rows leave step three and the BCH means follow the definition", {
  data <- .missing_data()
  fit <- .missing_fit(data)
  dropped <- c(2L, 17L, 33L, 101L, 250L, 251L, 400L)
  data$y[dropped] <- NA
  result <- three_step(fit, data, "y", missing = "exclude")

  kept <- setdiff(seq_len(nrow(data)), dropped)
  expected <- .bch_by_hand(fit$subject_posteriors[kept, , drop = FALSE],
                           fit$subject_profiles[kept], data$y[kept])
  expect_equal(result$estimate, expected, tolerance = 1e-10)
  # The fit is unchanged: step one is not refitted on the reduced data.
  expect_identical(nrow(fit$subject_posteriors), nrow(data))
})

test_that("class differences use only the units with the outcome", {
  data <- .missing_data(seed = 5L)
  fit <- .missing_fit(data)
  dropped <- c(4L, 90L, 91L, 300L)
  data$y[dropped] <- NA
  with_na <- three_step(fit, data, "y", missing = "exclude", contrast = "pairs")
  kept <- setdiff(seq_len(nrow(data)), dropped)
  expected <- .bch_by_hand(fit$subject_posteriors[kept, , drop = FALSE],
                           fit$subject_profiles[kept], data$y[kept])
  expect_equal(with_na$estimate, expected[2] - expected[1], tolerance = 1e-10)
})

test_that("person-level outcomes: groups missing the outcome leave step three", {
  data <- .missing_data(seed = 7L)
  fit <- .missing_fit(data, n_group_classes = 2)
  missing_groups <- c(3L, 18L, 44L)
  data$person_y[data$g %in% missing_groups] <- NA
  result <- three_step(fit, data, "person_y", level = "groups",
                       missing = "exclude")
  expect_identical(nrow(result), 2L)
  expect_true(all(is.finite(result$estimate)))
  expect_true(all(is.finite(result$standard_error)))

  kept_groups <- setdiff(seq_len(60L), missing_groups)
  first_rows <- match(kept_groups, data$g)
  expected <- .bch_by_hand(fit$group_posteriors[kept_groups, , drop = FALSE],
                           fit$group_classes[kept_groups],
                           data$person_y[first_rows])
  expect_equal(result$estimate, expected, tolerance = 1e-10)

  overall <- three_step(fit, data, "person_y", level = "groups",
                        missing = "exclude", contrast = "overall")
  expect_identical(nrow(overall), 1L)
  expect_equal(overall$df, 1)
})

test_that("a group's outcome missing on only some rows is refused", {
  data <- .missing_data(seed = 7L)
  fit <- .missing_fit(data, n_group_classes = 2)
  data$person_y[1L] <- NA
  expect_error(three_step(fit, data, "person_y", level = "groups",
                          missing = "exclude"),
               class = "latents_bad_outcome")
})

test_that("a categorical outcome with missing values is analysed on the rest", {
  data <- .missing_data()
  fit <- .missing_fit(data)
  data$passed <- data$y > 2.5
  data$passed[c(5L, 6L, 200L)] <- NA
  result <- three_step(fit, data, "passed", missing = "exclude")
  totals <- tapply(result$estimate, result$class, sum)
  expect_equal(as.vector(totals), c(1, 1), tolerance = 1e-10)
  # Each class's probability of TRUE is the BCH mean of the indicator over
  # the units that have the outcome, and only those.
  kept <- which(!is.na(data$passed))
  expected <- .bch_by_hand(fit$subject_posteriors[kept, , drop = FALSE],
                           fit$subject_profiles[kept],
                           as.numeric(data$passed[kept]))
  estimates <- as.data.frame(result)
  expect_equal(subset(estimates, category == "TRUE")$estimate, expected,
               tolerance = 1e-10)
  expect_error(three_step(fit, data, "passed"), class = "latents_bad_argument")
})
