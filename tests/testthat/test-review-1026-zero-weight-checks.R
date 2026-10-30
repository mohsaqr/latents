# Review 2026-10-26, ML-2: a unit of sampling weight zero is absent from the
# fit, so the data checks must judge the rows of positive weight only. Each
# scenario below is refused on its weighted rows; the zero-weight rows added
# to it would, if counted, make it look acceptable.

# Ten weighted groups of five rows, two clear profiles on `y1`, and `n_zero`
# zero-weight groups appended after them.
skip_on_cran()

.zw_frame <- function(n_positive = 10L, size = 5L, n_zero = 10L,
                      zero_size = 5L) {
  set.seed(1026)
  positive <- data.frame(
    g = rep(seq_len(n_positive), each = size),
    y1 = rep(c(0, 4), length.out = n_positive * size) +
      stats::rnorm(n_positive * size, sd = 0.5),
    w = 1)
  zero <- data.frame(
    g = n_positive + rep(seq_len(n_zero), each = zero_size),
    y1 = stats::rnorm(n_zero * zero_size, sd = 2),
    w = 0)
  rbind(positive, zero)
}

# One scenario per check: the data (zero-weight rows included), the call that
# fits it, and the condition class the weighted rows alone deserve.
.zw_scenarios <- function() {
  constant <- .zw_frame()
  constant$y2 <- ifelse(constant$w > 0, 0.1, stats::rnorm(nrow(constant)))

  unobserved <- .zw_frame()
  unobserved$y2 <- ifelse(unobserved$w > 0, NA, stats::rnorm(nrow(unobserved)))

  one_category <- .zw_frame()
  one_category$c1 <- ifelse(one_category$w > 0, 1L,
                            rep(1:2, length.out = nrow(one_category)))

  no_category <- .zw_frame()
  no_category$c1 <- ifelse(no_category$w > 0, NA_integer_,
                           rep(1:2, length.out = nrow(no_category)))

  # Two distinct indicator rows on the weighted groups; three profiles asked.
  two_rows <- .zw_frame()
  two_rows$y1 <- ifelse(two_rows$w > 0, rep(c(0, 1), length.out = nrow(two_rows)),
                        two_rows$y1)

  # Two weighted groups; three group classes asked.
  two_groups <- .zw_frame(n_positive = 2L, size = 10L)

  # Weighted groups are singletons; the zero-weight groups are not.
  singletons <- .zw_frame(n_positive = 30L, size = 1L)

  list(
    constant = list(data = constant, class = "latents_bad_data",
      fit = function(data) multilpa(data, vars = c("y1", "y2"), id = "g",
        n_profiles = 2, weights = "w", n_starts = 1, seed = 1)),
    constant_single_level = list(data = constant, class = "latents_bad_data",
      fit = function(data) lpa(data, vars = c("y1", "y2"), n_profiles = 2,
        weights = "w", n_starts = 1, seed = 1)),
    unobserved = list(data = unobserved, class = "latents_bad_data",
      fit = function(data) multilpa(data, vars = c("y1", "y2"), id = "g",
        n_profiles = 2, weights = "w", missing = "fiml", n_starts = 1,
        seed = 1)),
    one_category = list(data = one_category, class = "latents_bad_categorical",
      fit = function(data) multilpa(data, vars = c("y1", "c1"), id = "g",
        categorical = "c1", n_profiles = 2, weights = "w", n_starts = 1,
        seed = 1)),
    no_category = list(data = no_category, class = "latents_bad_categorical",
      fit = function(data) multilpa(data, vars = c("y1", "c1"), id = "g",
        categorical = "c1", n_profiles = 2, weights = "w", missing = "fiml",
        n_starts = 1, seed = 1)),
    distinct_rows = list(data = two_rows, class = "latents_unidentified",
      fit = function(data) multilpa(data, vars = "y1", id = "g",
        n_profiles = 3, weights = "w", n_starts = 1, seed = 1, max_iter = 20)),
    group_classes = list(data = two_groups, class = "latents_unidentified",
      fit = function(data) multilpa(data, vars = "y1", id = "g",
        n_profiles = 2, n_group_classes = 3, weights = "w", n_starts = 1,
        seed = 1, max_iter = 20)),
    singleton_groups = list(data = singletons, class = "latents_unidentified",
      fit = function(data) multilpa(data, vars = "y1", id = "g",
        n_profiles = 2, n_group_classes = 2, weights = "w", n_starts = 1,
        seed = 1, max_iter = 20)),
    person_centering = list(data = singletons, class = "latents_bad_data",
      fit = function(data) multilpa(data, vars = "y1", id = "g",
        n_profiles = 2, centering = "person", weights = "w", n_starts = 1,
        seed = 1, max_iter = 20)))
}

# The class of the condition a fit raises, or "fitted" when it returns.
.zw_outcome <- function(fit, data) {
  result <- tryCatch(fit(data), error = identity)
  if (inherits(result, "error")) class(result)[[1L]] else "fitted"
}

test_that("an indicator constant on the weighted rows is refused", {
  scenario <- .zw_scenarios()$constant
  expect_error(scenario$fit(scenario$data), class = "latents_bad_data")
})

test_that("a single-level fit refuses an indicator constant on the weighted rows", {
  scenario <- .zw_scenarios()$constant_single_level
  expect_error(scenario$fit(scenario$data), class = "latents_bad_data")
})

test_that("an indicator observed only in zero-weight rows is refused", {
  scenario <- .zw_scenarios()$unobserved
  expect_error(scenario$fit(scenario$data), class = "latents_bad_data")
})

test_that("a categorical indicator with one weighted category is refused", {
  scenario <- .zw_scenarios()$one_category
  expect_error(scenario$fit(scenario$data), class = "latents_bad_categorical")
})

test_that("a categorical indicator missing on every weighted row is refused", {
  scenario <- .zw_scenarios()$no_category
  expect_error(scenario$fit(scenario$data), class = "latents_bad_categorical")
})

test_that("zero-weight rows do not count as distinct indicator rows", {
  scenario <- .zw_scenarios()$distinct_rows
  expect_error(scenario$fit(scenario$data), class = "latents_unidentified")
})

test_that("zero-weight groups do not count towards group classes", {
  scenario <- .zw_scenarios()$group_classes
  expect_error(scenario$fit(scenario$data), class = "latents_unidentified")
})

test_that("zero-weight groups do not lift the singleton-group refusal", {
  scenario <- .zw_scenarios()$singleton_groups
  expect_error(scenario$fit(scenario$data), class = "latents_unidentified")
})

test_that("zero-weight groups cannot keep a person-centred column alive", {
  scenario <- .zw_scenarios()$person_centering
  expect_error(scenario$fit(scenario$data), class = "latents_bad_data")
})

test_that("each check gives one outcome with the zero-weight rows present or dropped", {
  scenarios <- .zw_scenarios()
  outcomes <- vapply(scenarios, function(scenario) {
    c(present = .zw_outcome(scenario$fit, scenario$data),
      dropped = .zw_outcome(scenario$fit, subset(scenario$data, w > 0)),
      expected = scenario$class)
  }, character(3))
  expect_identical(outcomes["present", ], outcomes["dropped", ])
  expect_identical(outcomes["present", ], outcomes["expected", ])
})

test_that("person scaling ignores a constant zero-weight group", {
  data <- .zw_frame()
  data$y1[data$g == 11] <- 2.5
  present <- multilpa(data, vars = "y1", id = "g", n_profiles = 2,
                      centering = "person", scaling = "person", weights = "w",
                      n_starts = 1, seed = 1)
  dropped <- multilpa(subset(data, w > 0), vars = "y1", id = "g",
                      n_profiles = 2, centering = "person", scaling = "person",
                      weights = "w", n_starts = 1, seed = 1)
  expect_s3_class(present, "multilpa")
  expect_equal(as.numeric(logLik(present)), as.numeric(logLik(dropped)),
               tolerance = 1e-6)
  # A weighted group that is constant is still refused.
  data$y1[data$g == 1] <- 2.5
  expect_error(multilpa(data, vars = "y1", id = "g", n_profiles = 2,
                        centering = "person", scaling = "person",
                        weights = "w", n_starts = 1, seed = 1),
               class = "latents_bad_data")
})

test_that("the categorical check counts weighted values only when weights exclude rows", {
  value <- c(1, 1, 1, 2)
  expect_error(.multilpa_categorical_levels(value, "c",
                                            contributing = c(TRUE, TRUE, TRUE, FALSE)),
               class = "latents_bad_categorical")
  # Every category is still coded, so each row keeps a code.
  expect_identical(.multilpa_categorical_levels(value, "c",
                                                contributing = rep(TRUE, 4L)),
                   c(1, 2))
})

test_that("rows of positive weight are read leniently before the weights are validated", {
  data <- data.frame(w = c(0, 1, 2))
  expect_identical(.multilpa_contributing_rows(data, "w"), c(FALSE, TRUE, TRUE))
  expect_identical(.multilpa_contributing_rows(data, NULL), rep(TRUE, 3L))
  # Weights the resolver will refuse leave every row in, so its refusal wins.
  expect_identical(.multilpa_contributing_rows(data.frame(w = c(-1, 1, 2)), "w"),
                   rep(TRUE, 3L))
  expect_identical(.multilpa_contributing_rows(data.frame(w = c(0, 0, 0)), "w"),
                   rep(TRUE, 3L))
})
