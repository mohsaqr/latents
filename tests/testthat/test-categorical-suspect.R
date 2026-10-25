suspect_data <- function() {
  .latents_local_seed(4101L)
  n <- 120L
  data.frame(g = rep(seq_len(24L), each = 5L), occasion = rep(seq_len(5L), 24L),
             item = sample(1:5, n, replace = TRUE),
             wide = sample(1:14, n, replace = TRUE),
             decimal = round(stats::rnorm(n, 3), 2))
}

suspect_warnings <- function(expr) {
  caught <- list()
  withCallingHandlers(suppressMessages(expr), warning = function(condition) {
    if (inherits(condition, "latents_suspect_categorical")) {
      caught[[length(caught) + 1L]] <<- condition
    }
    invokeRestart("muffleWarning")
  })
  caught
}

test_that("integer codes with at most ten categories are not reported", {
  d <- suspect_data()
  expect_length(suspect_warnings(
    multilca(d, "item", id = "g", n_profiles = 2, n_starts = 1, seed = 1, max_iter = 5)), 0L)
  d$ten <- rep(1:10, length.out = nrow(d))
  expect_length(suspect_warnings(
    multilca(d, c("item", "ten"), id = "g", n_profiles = 2, n_starts = 1, seed = 1,
             max_iter = 5)), 0L)
})

test_that("more than ten categories and non-integer values are reported once", {
  d <- suspect_data()
  wide <- suspect_warnings(
    multilca(d, c("item", "wide"), id = "g", n_profiles = 2, n_starts = 1, seed = 1,
             max_iter = 5))
  expect_length(wide, 1L)
  expect_identical(wide[[1L]]$indicators, "wide")
  expect_match(conditionMessage(wide[[1L]]), "`wide` has 14 categories (1 to 14)", fixed = TRUE)
  expect_match(conditionMessage(wide[[1L]]), "latent profile model", fixed = TRUE)
  both <- suspect_warnings(
    multilca(d, c("decimal", "item", "wide"), id = "g", n_profiles = 2, n_starts = 1,
             seed = 1, max_iter = 5))
  expect_length(both, 1L)
  expect_identical(both[[1L]]$indicators, c("decimal", "wide"))
  expect_match(conditionMessage(both[[1L]]), "2 categorical indicators may be continuous",
               fixed = TRUE)
  expect_match(conditionMessage(both[[1L]]), "`decimal` has non-integer values", fixed = TRUE)
})

test_that("factor indicators are never reported, whatever their level count", {
  d <- suspect_data()
  d$wide <- factor(d$wide)
  expect_length(suspect_warnings(
    multilca(d, "wide", id = "g", n_profiles = 2, n_starts = 1, seed = 1, max_iter = 5)), 0L)
})

test_that("an enumeration reports a suspect indicator once, not per candidate", {
  d <- suspect_data()
  lca <- suspect_warnings(enumerate_classes(d, vars = c("item", "wide"), id = "g",
    type = "lca", n_profiles = 2:3, n_group_classes = 1:2, n_starts = 1, seed = 1,
    max_iter = 5))
  expect_length(lca, 1L)
  lta <- suspect_warnings(enumerate_classes(d, vars = c("item", "wide"), id = "g",
    type = "lta", time = "occasion", categorical = c("item", "wide"),
    n_profiles = 2:3, n_group_classes = 1L, n_starts = 1, seed = 1, max_iter = 5))
  expect_length(lta, 1L)
})

test_that("a direct transition fit reports a suspect categorical indicator", {
  d <- suspect_data()
  caught <- suspect_warnings(lta(d, c("item", "wide"), id = "g", time = "occasion",
    categorical = c("item", "wide"), n_profiles = 2, n_starts = 1, seed = 1,
    max_iter = 5))
  expect_length(caught, 1L)
  expect_s3_class(caught[[1L]], "latents_suspect_categorical")
})
