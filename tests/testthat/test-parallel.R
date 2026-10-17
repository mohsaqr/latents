# Fitting candidates on several cores changes the time an enumeration takes and
# nothing it returns: every candidate sets and restores its own seed.

strip_calls <- function(grid) {
  grid$call <- NULL
  grid$fits <- lapply(grid$fits, function(fit) {
    if (!is.null(fit)) fit$call <- NULL
    fit
  })
  grid
}

collect_conditions <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(expr,
    warning = function(condition) {
      seen <<- c(seen, paste("warning:", conditionMessage(condition)))
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      seen <<- c(seen, paste("message:", conditionMessage(condition)))
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = seen)
}

test_that("a profile grid on several cores is identical to one core", {
  skip_on_cran()
  skip_on_os("windows")
  data <- subset(course_engagement, student <= 40)
  vars <- c("browse", "lectures", "forum_read")
  one <- enumerate_classes(data, vars, "student", n_profiles = 1:3,
                           n_group_classes = 1:2, model = "VVI", n_starts = 2,
                           seed = 1)
  two <- enumerate_classes(data, vars, "student", n_profiles = 1:3,
                           n_group_classes = 1:2, model = "VVI", n_starts = 2,
                           seed = 1, cores = 2)
  expect_identical(strip_calls(two), strip_calls(one))
  expect_identical(as.data.frame(two), as.data.frame(one))
})

test_that("transition and family grids on several cores are identical to one core", {
  skip_on_cran()
  skip_on_os("windows")
  data <- subset(course_engagement, student <= 30)
  vars <- c("browse", "lectures")
  transition <- lapply(1:2, \(cores) enumerate_classes(
    data, vars, "student", time = "sequence", n_profiles = 2:3, n_starts = 2,
    seed = 1, cores = cores))
  expect_identical(strip_calls(transition[[2]]), strip_calls(transition[[1]]))
  family <- lapply(1:2, \(cores) collect_conditions(enumerate_classes(
    data, vars, "student", family = c("additive", "dispersion"),
    n_group_classes = 1:2, seed = 1, cores = cores)))
  expect_identical(strip_calls(family[[2]]$value), strip_calls(family[[1]]$value))
})

test_that("candidates' warnings reach the caller in grid order, as on one core", {
  skip_on_cran()
  skip_on_os("windows")
  vars <- c("cognitive_strategies", "intrinsic_value", "self_efficacy",
            "self_regulation", "test_anxiety")
  runs <- lapply(1:2, \(cores) collect_conditions(enumerate_lpa(
    srl, vars, n_profiles = 6:8, model = "VVI", n_starts = 3, seed = 1,
    cores = cores)))
  expect_true(length(runs[[1]]$conditions) > 0L)
  expect_identical(runs[[2]]$conditions, runs[[1]]$conditions)
  expect_identical(as.data.frame(runs[[2]]$value)$warnings,
                   as.data.frame(runs[[1]]$value)$warnings)
})

test_that("the socket backend returns values and conditions in order", {
  skip_on_cran()
  skip_if_not(requireNamespace("latents", quietly = TRUE))
  fun <- function(i) {
    if (i == 2L) warning(warningCondition(paste("warned", i), class = "test_warning"))
    message(paste("message", i))
    i * 10L
  }
  seen <- list()
  values <- withCallingHandlers(
    .latents_map_candidates(1:3, fun, cores = 2L, backend = "socket"),
    warning = function(condition) {
      seen[[length(seen) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      seen[[length(seen) + 1L]] <<- condition
      invokeRestart("muffleMessage")
    })
  expect_identical(values, list(10L, 20L, 30L))
  expect_identical(vapply(seen, conditionMessage, character(1)),
                   c("message 1\n", "warned 2", "message 2\n", "message 3\n"))
  expect_true(inherits(seen[[2L]], "test_warning"))
})

test_that("cores is validated by class, and one core is lapply", {
  data <- subset(course_engagement, student <= 10)
  invisible(lapply(list(0, 1.5, NA, "2", c(1, 2)), function(cores) {
    expect_error(enumerate_classes(data, "browse", "student", n_profiles = 1,
                                   n_group_classes = 1, cores = cores),
                 class = "latents_bad_argument")
  }))
  expect_identical(.latents_map_candidates(1:3, \(i) i^2, cores = 1L), list(1, 4, 9))
})
