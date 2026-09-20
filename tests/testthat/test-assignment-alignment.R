# `assignments()` and `bivariate_residuals()` pair every row of the frame they
# are given with the posterior the fit holds for that position. A row count is
# not evidence that the two are in the same order, so the columns the fit
# recognises are compared against it row by row.

.aligned_data <- function(seed = 41L, n_groups = 12L, per = 12L) {
  set.seed(seed)
  n <- n_groups * per
  half <- rep(c(-4, 4), each = n / 2L)
  data.frame(g = rep(seq_len(n_groups), each = per),
             a = half + stats::rnorm(n, sd = 0.3),
             b = half + stats::rnorm(n, sd = 0.3),
             outcome = stats::rnorm(n))
}

.aligned_fit <- function(data) {
  multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
           n_starts = 2, seed = 1)
}

test_that("a conforming frame is joined exactly as before", {
  data <- .aligned_data()
  fit <- .aligned_fit(data)
  supplied <- assignments(fit, data = data)
  carried <- assignments(fit)

  expect_identical(nrow(supplied), nrow(data))
  expect_identical(supplied$profile, fit$subject_profiles)
  expect_identical(supplied$profile, carried$profile)
  # Columns the model never saw come through untouched, which is the point of
  # supplying a frame at all.
  expect_identical(supplied$outcome, data$outcome)
  expect_identical(supplied$a, data$a)
})

test_that("a reordered frame is refused rather than silently misaligned", {
  data <- .aligned_data()
  fit <- .aligned_fit(data)
  reversed <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  expect_error(assignments(fit, data = reversed),
               class = "multilpa_bad_inference_data")

  # A permutation within the fitted group order is caught by the indicators,
  # not by the identifier, so the check is not merely a group-order check.
  set.seed(11)
  shuffled <- data
  shuffled$a <- sample(data$a)
  expect_error(assignments(fit, data = shuffled),
               class = "multilpa_bad_inference_data")
})

test_that("a frame with no column of the fit warns that it cannot be checked", {
  data <- .aligned_data()
  fit <- .aligned_fit(data)
  # Nothing here can be compared with the fit, so the row order is assumed;
  # saying so is the contract, and it is a warning rather than silence.
  expect_warning(assignments(fit, data = data.frame(outcome = data$outcome)),
                 class = "multilpa_unverified_alignment")
})

test_that("the row-count contract is unchanged", {
  data <- .aligned_data()
  fit <- .aligned_fit(data)
  # A wrong number of rows is a data contract failure, not a model-nesting one;
  # `multilpa_bad_nesting` is reserved for comparing two models.
  expect_error(assignments(fit, data = head(data, 10L)),
               class = "multilpa_bad_inference_data")
})

test_that("bivariate residuals accept the fitted frame and refuse another order", {
  data <- .aligned_data()
  fit <- .aligned_fit(data)
  carried <- bivariate_residuals(fit)
  supplied <- bivariate_residuals(fit, data)

  # Supplying the frame the fit was built from changes nothing.
  expect_equal(supplied$residual, carried$residual)
  expect_identical(supplied$profile, carried$profile)

  reversed <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  expect_error(bivariate_residuals(fit, reversed),
               class = "multilpa_bad_inference_data")
  expect_error(bivariate_residuals(fit, reversed, by = "overall"),
               class = "multilpa_bad_inference_data")
})

test_that("starting values carry the labels a categorical block is held by", {
  set.seed(5)
  n <- 240L
  groups <- rep(seq_len(24L), each = 10L)
  latent <- rep(rep(c(1L, 2L), length.out = 24L), each = 10L)
  draw <- function(probability) factor(
    ifelse(stats::runif(n) < ifelse(latent == 1L, probability, 1 - probability),
           "yes", "no"), levels = c("no", "yes"))
  data <- data.frame(g = groups, i1 = draw(0.95), i2 = draw(0.6),
                     i3 = draw(0.75))
  items <- c("i1", "i2", "i3")
  fit <- multilpa(data, items, "g", n_profiles = 2, n_group_classes = 1,
                  categorical = items, n_starts = 3, seed = 1)
  start <- starting_values(fit, what = "measurement")

  # The labels are the only record of which block was estimated for which item.
  expect_identical(names(start$response_probabilities), items)
  expect_identical(colnames(start$response_probabilities[["i1"]]),
                   c("no", "yes"))

  held <- function(order) multilpa(data, order, "g", n_profiles = 2,
                                   n_group_classes = 1, categorical = order,
                                   fixed = "measurement", start = start,
                                   n_starts = 1, seed = 1)
  # Naming the same items in another order is the same held measurement model,
  # so it must reach the same likelihood rather than a quietly different one.
  expect_equal(as.numeric(logLik(held(c("i3", "i1", "i2")))),
               as.numeric(logLik(held(items))), tolerance = 1e-10)

  # A label that names an item this fit does not have is refused, not guessed.
  mislabelled <- start
  names(mislabelled$response_probabilities) <- c("i1", "i2", "elsewhere")
  expect_error(multilpa(data, items, "g", n_profiles = 2, n_group_classes = 1,
                        categorical = items, fixed = "measurement",
                        start = mislabelled, n_starts = 1, seed = 1),
               class = "multilpa_bad_start")

  # The tidy view reports the labels the start carries, matching what the fitted
  # object's own response table reports for the same quantity.
  responses <- as.data.frame(start, what = "responses")
  expect_type(responses$indicator, "character")
  expect_type(responses$category, "character")
  expect_setequal(unique(responses$indicator), items)

  # An unlabelled start is still supported, and has only positions to report.
  positional <- start
  positional$response_probabilities <-
    unname(lapply(positional$response_probabilities, unname))
  bare <- as.data.frame(positional, what = "responses")
  expect_type(bare$indicator, "integer")
  expect_type(bare$category, "integer")
})

test_that("alignment survives a round trip that changes storage, not order", {
  data <- .aligned_data()
  fit <- .aligned_fit(data)
  # The identifier read back as text and the indicators read back from a
  # written file are the same rows in the same order, and must still align.
  as_text <- data
  as_text$g <- as.character(as_text$g)
  expect_identical(assignments(fit, data = as_text)$profile,
                   fit$subject_profiles)

  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE, after = FALSE)
  utils::write.csv(data, path, row.names = FALSE)
  round_tripped <- utils::read.csv(path)
  expect_identical(assignments(fit, data = round_tripped)$profile,
                   fit$subject_profiles)
})
