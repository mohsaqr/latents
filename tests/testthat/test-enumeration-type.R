enumeration_type_data <- function() {
  .latents_local_seed(731L)
  data <- data.frame(person = rep(seq_len(12L), each = 6L),
                     time = rep(seq_len(6L), 12L))
  state <- (data$person + data$time) %% 2L
  data$a <- stats::rbinom(nrow(data), 1L, 0.2 + 0.6 * state)
  data$b <- stats::rbinom(nrow(data), 2L, 0.2 + 0.6 * state)
  data$c <- stats::rbinom(nrow(data), 1L, 0.8 - 0.6 * state)
  data$x <- 2 * state + stats::rnorm(nrow(data), sd = 0.5)
  data$y <- -state + stats::rnorm(nrow(data), sd = 0.6)
  data
}

enumeration_type_capture <- function(expr) {
  conditions <- list()
  value <- withCallingHandlers(expr,
    warning = function(condition) {
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = conditions)
}

enumeration_type_without_calls <- function(value) {
  value$call <- NULL
  if (!is.null(value$fits)) {
    value$fits <- lapply(value$fits, function(fit) {
      if (!is.null(fit)) fit$call <- NULL
      fit
    })
  }
  value
}

test_that("enumeration appends type without moving existing formals", {
  expect_identical(names(formals(enumerate_classes)),
    c("data", "vars", "id", "n_profiles", "n_group_classes", "model",
      "seed", "family", "between_variance", "time", "...", "cores", "type"))
  expect_null(formals(enumerate_classes)$type)
})

test_that("explicit enumeration types preserve equivalent inferred fits", {
  data <- enumeration_type_data()
  categorical <- c("a", "b", "c")
  settings <- list(
    lca = list(vars = categorical, categorical = categorical),
    lpa = list(vars = c("x", "y"), model = "VVI"),
    lta = list(vars = categorical, categorical = categorical, time = "time"))
  invisible(lapply(names(settings), function(type) {
    arguments <- c(list(data = data, id = "person", n_profiles = 2:3,
                        n_group_classes = 1L, n_starts = 2L,
                        max_iter = 15L, seed = 34L), settings[[type]])
    inferred <- enumeration_type_capture(do.call(enumerate_classes, arguments))
    forced_arguments <- arguments
    if (type == "lca") forced_arguments$categorical <- NULL
    explicit <- enumeration_type_capture(do.call(enumerate_classes,
      c(forced_arguments, list(type = type))))
    expect_identical(explicit$value$type, type)
    expect_identical(inferred$value$type, type)
    expect_identical(enumeration_type_without_calls(explicit$value),
                     enumeration_type_without_calls(inferred$value))
    expect_identical(explicit$conditions, inferred$conditions)
    expect_identical(as.data.frame(explicit$value), as.data.frame(inferred$value))
    expect_true(all(is.finite(explicit$value$table$log_likelihood)))
    expect_true(all(is.na(explicit$value$table$error)))
    expect_output(print(explicit$value), toupper(type))
    null <- enumeration_type_capture(do.call(enumerate_classes,
      c(arguments, list(type = NULL))))
    expect_identical(enumeration_type_without_calls(null$value),
                     enumeration_type_without_calls(inferred$value))
    expect_identical(null$conditions, inferred$conditions)
  }))
})

test_that("lca type adds a supplied categorical subset like multilca", {
  data <- enumeration_type_data()
  vars <- c("a", "b", "c")
  arguments <- list(data = data, vars = vars, id = "person", n_profiles = 2L,
                    n_group_classes = 1L, n_starts = 2L,
                    max_iter = 15L, seed = 34L)
  explicit <- enumeration_type_capture(do.call(enumerate_classes,
    c(arguments, list(type = "lca", categorical = "a"))))
  direct <- enumeration_type_capture(do.call(multilca,
    c(arguments, list(categorical = "a"))))
  expect_identical(enumeration_type_without_calls(explicit$value$fits[[1L]]),
                   enumeration_type_without_calls(direct$value))
  expect_identical(explicit$conditions, direct$conditions)
  expect_identical(explicit$value$type, "lca")
})

test_that("lca type preserves ordinal and count indicators like multilca", {
  data <- enumeration_type_data()
  arguments <- list(data = data, vars = c("a", "b", "c"), id = "person",
                    n_profiles = 2L, n_group_classes = 1L, ordinal = "b",
                    count = "c", categorical = "a", n_starts = 2L,
                    max_iter = 15L, seed = 34L)
  explicit <- enumeration_type_capture(do.call(enumerate_classes,
    c(arguments, list(type = "lca"))))
  direct <- enumeration_type_capture(do.call(multilca, arguments))
  expect_identical(enumeration_type_without_calls(explicit$value$fits[[1L]]),
                   enumeration_type_without_calls(direct$value))
  expect_identical(explicit$conditions, direct$conditions)
  expect_identical(explicit$value$type, "lca")
})

test_that("lta type preserves Gaussian and mixed measurement", {
  data <- enumeration_type_data()
  settings <- list(list(vars = c("x", "y")),
                    list(vars = c("a", "x"), categorical = "a"))
  invisible(lapply(settings, function(setting) {
    arguments <- c(list(data = data, id = "person", time = "time",
                        n_profiles = 2L, n_group_classes = 1L,
                        n_starts = 1L, max_iter = 10L, seed = 34L), setting)
    inferred <- enumeration_type_capture(do.call(enumerate_classes, arguments))
    explicit <- enumeration_type_capture(do.call(enumerate_classes,
      c(arguments, list(type = "lta"))))
    expect_identical(enumeration_type_without_calls(explicit$value),
                     enumeration_type_without_calls(inferred$value))
    expect_identical(explicit$conditions, inferred$conditions)
    expect_identical(explicit$value$type, "lta")
    expect_true(all(is.finite(explicit$value$table$log_likelihood)))
    expect_true(all(is.na(explicit$value$table$error)))
  }))
})

test_that("enumeration type validation rejects invalid and contradictory choices", {
  data <- enumeration_type_data()
  arguments <- list(data = data, vars = c("x", "y"), id = "person",
                    n_profiles = 2L, n_group_classes = 1L,
                    n_starts = 1L, max_iter = 10L, seed = 34L)
  invalid <- list("other", "l", "", NA_character_, c("lca", "lpa"),
                   character(), 1L, TRUE, list("lca"))
  invisible(lapply(invalid, function(type) {
    expect_error(do.call(enumerate_classes, c(arguments, list(type = type))),
                 class = "latents_bad_argument")
  }))
  conflicts <- list(
    list(type = "lta"),
    list(type = "lca", time = "time"),
    list(type = "lpa", time = "time"),
    list(type = "lpa", categorical = "x"),
    list(type = "lpa", ordinal = "x"),
    list(type = "lpa", count = "x"),
    list(type = "lca", family = "additive"),
    list(type = "lpa", family = "dispersion"),
    list(type = "lta", time = "time", family = "additive_dispersion"))
  invisible(lapply(conflicts, function(conflict) {
    expect_error(do.call(enumerate_classes, c(arguments, conflict)),
                 class = "latents_bad_argument")
    condition <- tryCatch(do.call(enumerate_classes, c(arguments, conflict)),
                           error = identity)
    expect_s3_class(condition, "latents_bad_argument")
    expect_match(conditionMessage(condition), "type")
    expect_match(conditionMessage(condition), "time|categorical|ordinal|count|family")
    expect_match(conditionMessage(condition), "[Uu]se|[Pp]ass|[Ss]et|[Ll]eave|[Rr]emove|[Ss]upply")
  }))
})

test_that("type prefixes and empty categorical declarations follow match.arg style", {
  data <- enumeration_type_data()
  arguments <- list(data = data, vars = c("x", "y"), id = "person",
                    n_profiles = 2L, n_group_classes = 1L, model = "VVI",
                    n_starts = 1L, max_iter = 10L, seed = 34L)
  full <- enumeration_type_capture(do.call(enumerate_classes,
    c(arguments, list(type = "lpa"))))
  # Empty categorical declarations are equivalent to the omitted declaration.
  invisible(lapply(list(NULL, character()), function(categorical) {
    prefix <- enumeration_type_capture(do.call(enumerate_classes,
      c(arguments, list(type = "lp", categorical = categorical))))
    expect_identical(enumeration_type_without_calls(prefix$value),
                     enumeration_type_without_calls(full$value))
    expect_identical(prefix$value$table, full$value$table)
    expect_identical(prefix$value$type, "lpa")
    expect_identical(prefix$conditions, full$conditions)
  }))
})

test_that("NULL retains mixed and group-family enumeration behaviour", {
  data <- enumeration_type_data()
  settings <- list(
    mixed = list(vars = c("a", "x"), n_profiles = 2L,
                  categorical = "a", model = "VVI"),
    family = list(vars = c("x", "y"), family = "additive"))
  invisible(lapply(names(settings), function(type) {
    arguments <- c(list(data = data, id = "person", n_group_classes = 1L,
                        n_starts = 1L, max_iter = 10L, seed = 34L),
                   settings[[type]])
    inferred <- enumeration_type_capture(do.call(enumerate_classes, arguments))
    null <- enumeration_type_capture(do.call(enumerate_classes,
      c(arguments, list(type = NULL))))
    expect_identical(null$value$type, type)
    expect_identical(enumeration_type_without_calls(null$value),
                     enumeration_type_without_calls(inferred$value))
    expect_identical(null$conditions, inferred$conditions)
    expect_true(all(is.finite(null$value$table$log_likelihood)))
  }))
})

test_that("lpa type retains the existing refusal of factor indicators", {
  data <- enumeration_type_data()
  data$x <- factor(data$a)
  arguments <- list(data = data, vars = c("x", "y"), id = "person",
                    n_profiles = 2L, n_group_classes = 1L, model = "VVI",
                    n_starts = 1L, max_iter = 10L, seed = 34L)
  inferred <- enumeration_type_capture(do.call(enumerate_classes, arguments))
  forced <- enumeration_type_capture(do.call(enumerate_classes,
    c(arguments, list(type = "lpa"))))
  expect_identical(enumeration_type_without_calls(forced$value),
                   enumeration_type_without_calls(inferred$value))
  expect_identical(forced$conditions, inferred$conditions)
  expect_null(forced$value$fits[[1L]])
  expect_match(forced$value$table$error, "[Nn]umeric|[Cc]ontinuous|[Ff]actor")
})
