# Review 2026-10-27: a condition raised by a fitting verb never carries the
# data in its call. A verb reached through do.call() recorded its arguments by
# value, so a rejected input printed, through a host that deparses the call in
# full, every value of the data set.

# More than 1.5 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

condition_data <- function() {
  withr::local_seed(20261027)
  data.frame(id = rep(seq_len(400L), each = 5L), t = rep(1:5, 400L),
             a = factor(sample(3L, 2000L, TRUE)), b = factor(sample(3L, 2000L, TRUE)),
             y = stats::rnorm(2000L))
}

call_bytes <- function(condition) sum(nchar(deparse(conditionCall(condition))))
printed_bytes <- function(condition) {
  nchar(paste(utils::capture.output(print(condition)), collapse = "\n"))
}

test_that("a bounded call keeps language and replaces data by a placeholder", {
  big <- as.list(seq_len(5000L))
  call <- as.call(list(function(x, y, z) NULL, x = quote(data), y = big, z = 1:3))
  bounded <- .latents_bounded_call(call)
  expect_identical(bounded[[1L]], as.symbol("<function>"))
  expect_identical(bounded$x, quote(data))
  expect_identical(bounded$y, as.symbol("<list[5000]>"))
  expect_identical(bounded$z, 1:3)
  expect_lt(sum(nchar(deparse(bounded))), 120L)
  frame <- data.frame(a = seq_len(1000L), b = 1)
  expect_identical(.latents_bounded_call(quote(f(x)))[[2L]], quote(x))
  expect_identical(.latents_bounded_call(as.call(list(as.symbol("f"), frame)))[[2L]],
                   as.symbol("<data.frame[1000 x 2]>"))
  expect_null(.latents_bounded_call(NULL))
})

test_that("the arguments of a call are bound to symbols, evaluated once in the caller's frame", {
  env <- new.env()
  caller <- new.env()
  caller$v <- c("a", "b")
  caller$evaluations <- 0L
  call <- as.call(list(as.symbol("f"), data = data.frame(a = 1:3), vars = quote(v),
                       n = quote({ evaluations <- evaluations + 1L; 2L }), 3L))
  bound <- .latents_bind_values(call, env, caller)
  expect_identical(bound$data, as.symbol(".data"))
  expect_identical(bound$vars, quote(v))
  expect_identical(bound$n, as.symbol(".n"))
  expect_identical(bound[[5L]], as.symbol(".argument_5"))
  expect_identical(get(".data", env), data.frame(a = 1:3))
  expect_identical(get("v", env), c("a", "b"))
  expect_identical(get(".n", env), 2L)
  expect_identical(get(".argument_5", env), 3L)
  # The expression ran once, in the caller's frame.
  expect_identical(caller$evaluations, 1L)
  expect_lt(sum(nchar(deparse(bound))), 80L)
})

test_that("every fitting verb refuses a non data frame with a small classed error", {
  data <- condition_data()
  # Any fit serves as the wrong `data`; this one is not identified (CORE-5).
  expect_warning(grid <- lca(data, c("a", "b"), n_classes = 2, n_starts = 1, seed = 1),
                 class = "latents_unidentified")
  attempts <- list(
    multilpa = function() multilpa(grid, vars = "y", id = "id", n_profiles = 2),
    multilca = function() multilca(grid, vars = c("a", "b"), id = "id", n_profiles = 2,
                                   n_group_classes = 2, time = "t", missing = "fiml"),
    lpa = function() lpa(grid, "y", n_profiles = 2),
    lca = function() lca(grid, c("a", "b"), n_classes = 2),
    lta = function() lta(grid, c("a", "b"), id = "id", n_profiles = 2, time = "t"),
    do_call = function() do.call(multilpa, list(data = grid, vars = "y", id = "id",
                                                n_profiles = 2)))
  invisible(lapply(names(attempts), function(name) {
    condition <- expect_error(attempts[[name]](), class = "latents_bad_data")
    expect_match(conditionMessage(condition), "data frame", label = name)
    expect_lt(call_bytes(condition), 300L, label = name)
    expect_lt(printed_bytes(condition), 600L, label = name)
  }))
})

test_that("a validation error inside the setup names the arguments, not the data", {
  data <- condition_data()
  # `vars` names a column that does not exist: refused inside the setup,
  # after the data frame check, so the condition's call is the setup call.
  condition <- expect_error(
    do.call(multilpa, list(data = data, vars = "nope", id = "id", n_profiles = 2)),
    class = "latents_bad_data")
  expect_lt(call_bytes(condition), 300L)
  condition <- expect_error(
    do.call(multilca, list(data = data, vars = "nope", id = "id", n_profiles = 2)),
    class = "latents_bad_data")
  expect_lt(call_bytes(condition), 300L)
})

test_that("a collected candidate error and a shared data error carry bounded calls", {
  collector <- .latents_task_collector(function(i) {
    do.call(function(data) stop("bad candidate"), list(data = as.list(seq_len(5000L))))
  })
  result <- collector(1L)$value
  expect_s3_class(result, "latents_task_error")
  expect_lt(call_bytes(result$condition), 120L)
  bad <- do.call(function(data) {
    errorCondition("the data are bad", class = "latents_bad_data", call = sys.call())
  }, list(data = as.list(seq_len(5000L))))
  condition <- expect_error(.latents_refuse_shared_data_error(list(bad, bad)),
                            class = "latents_bad_data")
  expect_lt(call_bytes(condition), 120L)
})

test_that("the fits themselves are unchanged by the symbol forwarding", {
  data <- condition_data()
  direct <- multilpa(data, vars = c("a", "b"), id = "id", n_profiles = 2,
                     n_group_classes = 2, categorical = c("a", "b"), n_starts = 1, seed = 1)
  through <- multilca(data, vars = c("a", "b"), id = "id", n_profiles = 2,
                      n_group_classes = 2, n_starts = 1, seed = 1)
  expect_equal(through$log_likelihood, direct$log_likelihood, tolerance = 1e-12)
  expect_identical(through$call[[1L]], as.symbol("multilca"))
  by_value <- do.call(multilpa, list(data = data, vars = c("a", "b"), id = "id",
                                     n_profiles = 2, n_group_classes = 2,
                                     categorical = c("a", "b"), n_starts = 1, seed = 1))
  expect_equal(by_value$log_likelihood, direct$log_likelihood, tolerance = 1e-12)
  transition <- lta(data, c("a", "b"), id = "id", n_profiles = 2, time = "t",
                    categorical = c("a", "b"), n_starts = 1, seed = 1)
  expect_true(transition$converged)
  expect_identical(transition$call[[1L]], as.symbol("lta"))
})
