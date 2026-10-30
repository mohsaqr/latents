# Review 2026-10-26: rows with nothing observed are not classified, so the
# classification error matrix, the BCH weights and every three-step estimator
# leave them out, as entropy, the average posteriors and the information
# criteria already did. Before the fix the error matrix counted them: a row
# whose posterior is the class prior is assigned to the larger class, and on
# this fixture P(assigned 1 | true 2) was 0.025 rather than about 5e-5.

skip_on_cran()

empty_row_case <- function() {
  data <- datasets::faithful[seq_len(120L), ]
  empty <- c(5L, 9L, 30L)
  data[empty, ] <- NA
  data$unit <- factor(seq_len(nrow(data)))
  fit <- quietly(multilpa(data, vars = c("eruptions", "waiting"), id = "unit",
                          n_profiles = 2, n_group_classes = 1, n_starts = 2,
                          seed = 1, missing = "fiml"))
  list(data = data, empty = empty, fit = fit)
}

# The Vermunt (2010) error matrix written from its definition, on the rows
# given: rows are the true class, columns the modal assignment.
hand_error_matrix <- function(posteriors) {
  modal <- max.col(posteriors, ties.method = "first")
  classes <- seq_len(ncol(posteriors))
  t(vapply(classes, function(true_class) {
    vapply(classes, function(assigned) {
      sum(posteriors[modal == assigned, true_class]) / sum(posteriors[, true_class])
    }, numeric(1))
  }, numeric(length(classes))))
}

test_that("the error matrix and BCH weights use only the classified rows", {
  case <- empty_row_case()
  fit <- case$fit
  informative <- -case$empty
  posteriors <- fit$subject_posteriors[informative, , drop = FALSE]
  expected <- hand_error_matrix(posteriors)
  errors <- get_results(fit, "classification_errors", level = "individuals")
  expect_equal(errors$probability,
               expected[cbind(errors$true_class, errors$assigned_class)],
               tolerance = 1e-12)
  weights <- get_results(fit, "bch_weights", level = "individuals")
  expect_setequal(unique(weights$unit), setdiff(seq_len(120L), case$empty))
  inverse <- solve(expected)
  expect_equal(weights$weight, inverse[cbind(weights$assigned_class, weights$class)],
               tolerance = 1e-10)
  # The same rows as the tables that already followed the rule.
  expect_identical(sum(get_results(fit, "classification")$n_modal[1:2]), 117L)
})

test_that("step three ignores an outcome recorded on a row with nothing observed", {
  case <- empty_row_case()
  withr::local_seed(3)
  data <- case$data
  data$outcome <- stats::rnorm(nrow(data))
  moved <- data
  moved$outcome[case$empty] <- 1e3
  methods <- c("bch", "proportional", "modal", "ml")
  invisible(lapply(methods, function(method) {
    reference <- three_step(case$fit, data, outcome = "outcome", method = method)
    shifted <- three_step(case$fit, moved, outcome = "outcome", method = method)
    expect_equal(as.data.frame(shifted)$estimate, as.data.frame(reference)$estimate,
                 tolerance = 1e-10, info = method)
  }))
  # And the BCH means are the closed form on the 117 classified rows.
  kept <- setdiff(seq_len(120L), case$empty)
  weights <- solve(hand_error_matrix(case$fit$subject_posteriors[kept, ]))[
    case$fit$subject_profiles[kept], ]
  expect_equal(as.data.frame(three_step(case$fit, data, outcome = "outcome",
                                        method = "bch"))$estimate,
               as.vector(crossprod(weights, data$outcome[kept])) / colSums(weights),
               tolerance = 1e-10)
})

test_that("r3step leaves out the empty rows and still runs by group class", {
  case <- empty_row_case()
  withr::local_seed(4)
  data <- case$data
  data$x <- stats::rnorm(nrow(data))
  moved <- data
  moved$x[case$empty] <- 50
  expect_equal(as.data.frame(r3step(case$fit, moved, covariates = "x"))$estimate,
               as.data.frame(r3step(case$fit, data, covariates = "x"))$estimate,
               tolerance = 1e-8)
})

test_that("a group with nothing observed is not a classified group", {
  withr::local_seed(5)
  groups <- 30L
  data <- data.frame(school = rep(seq_len(groups), each = 8L))
  shift <- c(-1.5, 1.5)[1L + (data$school %% 2L)]
  data$y1 <- stats::rnorm(nrow(data), shift)
  data$y2 <- stats::rnorm(nrow(data), -shift)
  data[data$school == 7L, c("y1", "y2")] <- NA
  fit <- quietly(multilpa(data, vars = c("y1", "y2"), id = "school",
                          n_profiles = 2, n_group_classes = 2, n_starts = 3,
                          seed = 1, missing = "fiml"))
  errors <- get_results(fit, "classification_errors", level = "groups")
  kept <- setdiff(seq_len(groups), 7L)
  expected <- hand_error_matrix(fit$group_posteriors[kept, , drop = FALSE])
  expect_equal(errors$probability,
               expected[cbind(errors$true_class, errors$assigned_class)],
               tolerance = 1e-12)
  expect_false(7L %in% get_results(fit, "bch_weights", level = "groups")$unit)
  data$distal <- stats::rnorm(groups)[data$school]
  expect_s3_class(three_step(fit, data, outcome = "distal", level = "groups"),
                  "data.frame")
})

test_that("drop_empty = FALSE keeps every row, as in 0.9.37 and earlier", {
  case <- empty_row_case()
  everything <- quietly(multilpa(case$data, vars = c("eruptions", "waiting"),
                                 id = "unit", n_profiles = 2, n_group_classes = 1,
                                 n_starts = 2, seed = 1, missing = "fiml",
                                 drop_empty = FALSE))
  errors <- get_results(everything, "classification_errors",
                        level = "individuals")
  expected <- hand_error_matrix(everything$subject_posteriors)
  expect_equal(errors$probability,
               expected[cbind(errors$true_class, errors$assigned_class)],
               tolerance = 1e-12)
})
