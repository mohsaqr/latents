# enumerate_classes() names, once per grid, the converged candidates whose
# best log likelihood was reached by one start only (latents_single_start).

single_start_data <- function() {
  withr::local_seed(1)
  data <- data.frame(id = rep(1:30, each = 6), y1 = stats::rnorm(180),
                     y2 = stats::rnorm(180), t = rep(1:6, 30))
  data$y1[1:90] <- data$y1[1:90] + 2
  data
}

# Every expected fixture condition except the one under test.
without_notice <- setdiff(.multilpa_expected_warnings, "latents_single_start")

collect_notices <- function(expr) {
  notices <- character()
  value <- withCallingHandlers(
    quietly(expr, classes = without_notice),
    latents_single_start = function(notice) {
      notices <<- c(notices, conditionMessage(notice))
      invokeRestart("muffleMessage")
    })
  list(value = value, notices = notices)
}

test_that("a profile grid names exactly its single-start candidates, once", {
  run <- collect_notices(enumerate_classes(
    single_start_data(), c("y1", "y2"), "id", n_profiles = 2:3,
    n_group_classes = 1:2, model = "VVV", n_starts = 4, seed = 1))
  grid <- as.data.frame(run$value)
  once <- subset(grid, converged & n_best_replicated %in% 1L)$candidate
  others <- setdiff(grid$candidate, once)
  expect_gt(length(once), 0L)
  expect_gt(length(others), 0L)
  expect_length(run$notices, 1L)
  expect_true(all(vapply(once, grepl, logical(1), x = run$notices, fixed = TRUE)))
  expect_false(any(vapply(others, grepl, logical(1), x = run$notices, fixed = TRUE)))
  expect_match(run$notices, sprintf("^%d of %d candidates", length(once),
                                    sum(grid$converged)))
})

test_that("a grid whose candidates all replicate raises no notice", {
  expect_no_message(
    quietly(enumerate_lpa(single_start_data(), c("y1", "y2"), n_profiles = 1:2,
                          n_starts = 5, seed = 1), classes = without_notice),
    class = "latents_single_start")
})

test_that("a transition grid raises the notice from its own table", {
  run <- collect_notices(enumerate_classes(
    single_start_data(), c("y1", "y2"), "id", n_profiles = 2:3, time = "t",
    n_starts = 3, seed = 1))
  grid <- as.data.frame(run$value)
  once <- subset(grid, converged & n_best_replicated %in% 1L)$candidate
  expect_length(run$notices, as.integer(length(once) > 0L))
  expect_true(all(vapply(once, grepl, logical(1), x = paste(run$notices, ""),
                         fixed = TRUE)))
})

test_that("the rule leaves out unconverged, unfitted and one-start candidates", {
  expect_identical(
    .latents_single_start_candidates(c("a", "b", "c", "d", "e", "f", "g"),
                                     c(TRUE, FALSE, TRUE, NA, TRUE, TRUE, TRUE),
                                     c(1L, 1L, 2L, 1L, NA, 1L, 1L),
                                     c(5L, 5L, 5L, 5L, 5L, 1L, NA)),
    "a")
  expect_identical(.latents_single_start_candidates(character(), logical(),
                                                    integer(), integer()),
                   character())
})

test_that("a grid fitted from one start per candidate raises no notice", {
  expect_no_message(
    quietly(enumerate_classes(single_start_data(), c("y1", "y2"), "id",
                              n_profiles = 2:3, n_group_classes = 1:2,
                              model = "VVV", n_starts = 1, seed = 1),
            classes = without_notice),
    class = "latents_single_start")
})

test_that("the notice is a classed message and lists at most six names", {
  result <- list(table = data.frame(converged = rep(TRUE, 8),
                                    n_best_replicated = c(rep(1L, 7), 3L)),
                 fits = stats::setNames(rep(list(list(n_starts = 4L, log_likelihood = -1)), 8),
                                        sprintf("cand%d", 1:8)))
  expect_message(.latents_single_start_notice(result), class = "latents_single_start")
  text <- tryCatch(.latents_single_start_notice(result),
                   latents_single_start = conditionMessage)
  expect_match(text, "^7 of 8 candidates")
  expect_match(text, "cand6, ...)", fixed = TRUE)
  expect_no_match(text, "cand7", fixed = TRUE)
  result$table$n_best_replicated <- 2L
  expect_no_message(.latents_single_start_notice(result))
  expect_identical(suppressMessages(.latents_single_start_notice(result)), result)
})

test_that("summary() prints the same sentence as the notice", {
  run <- collect_notices(enumerate_classes(
    single_start_data(), c("y1", "y2"), "id", n_profiles = 2:3,
    n_group_classes = 1:2, model = "VVV", n_starts = 4, seed = 1))
  printed <- paste(capture.output(print(summary(run$value))), collapse = " ")
  expected <- sub("\n$", "", run$notices)
  # The print wraps at 76 columns; compare the words.
  squash <- function(text) gsub("[[:space:]!]+", " ", text)
  expect_match(squash(printed), squash(expected), fixed = TRUE)
})
