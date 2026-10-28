# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

profile_pool_data <- function() {
  .latents_local_seed(1854L)
  group <- rep(seq_len(12L), each = 6L)
  state <- stats::rbinom(length(group), 1L, 0.25 + 0.5 * (group %% 2L))
  data.frame(group = group,
    a = stats::rbinom(length(group), 1L, 0.15 + 0.65 * state),
    b = stats::rbinom(length(group), 2L, 0.2 + 0.6 * state),
    c = stats::rbinom(length(group), 1L, 0.8 - 0.6 * state),
    x = 2 * state + stats::rnorm(length(group), sd = 0.6),
    y = -1.5 * state + stats::rnorm(length(group), sd = 0.8),
    predictor = stats::rnorm(length(group)))
}

profile_pool_capture <- function(expr) {
  conditions <- list()
  value <- withCallingHandlers(expr,
    warning = function(condition) {
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      # Progress messages report the worker count; they are not part of the
      # qualifications compared between serial and parallel fits.
      if (!inherits(condition, "latents_progress")) {
        conditions[[length(conditions) + 1L]] <<- condition
      }
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = conditions)
}

profile_pool_without_calls <- function(value) {
  value$call <- NULL
  if (!is.null(value$fits)) {
    value$fits <- lapply(value$fits, function(fit) {
      if (!is.null(fit)) fit$call <- NULL
      fit
    })
  }
  value
}

profile_pool_skip_forks <- function() {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(!.latents_fork_safe(), "This BLAS cannot run safely in forked workers")
}

profile_pool_arguments <- function(data, kind) {
  common <- list(data = data, id = "group", n_starts = 2L,
                 max_iter = 15L, seed = 12L)
  if (kind == "categorical") {
    c(common, list(vars = c("a", "b", "c"), categorical = c("a", "b", "c"),
                    model = NULL))
  } else {
    c(common, list(vars = c("x", "y"), model = kind,
                    n_screen = if (kind == "EEE") 4L else NULL))
  }
}

test_that("profile grids match serial and direct categorical and continuous fits", {
  profile_pool_skip_forks()
  data <- profile_pool_data()
  invisible(lapply(c("categorical", "EEE", "VVI"), function(kind) {
    arguments <- profile_pool_arguments(data, kind)
    invisible(lapply(c(FALSE, TRUE), function(single_level) {
      current <- arguments
      if (single_level) current["id"] <- list(NULL)
      grid_arguments <- c(current, list(n_profiles = 2:3))
      # Omitting n_group_classes is meaningful for single-level enumeration.
      if (!single_level) grid_arguments$n_group_classes <- 1:2
      serial <- profile_pool_capture(do.call(enumerate_classes,
        c(grid_arguments, list(cores = 1L))))
      spread <- profile_pool_capture(do.call(enumerate_classes,
        c(grid_arguments, list(cores = 3L))))
      expect_true(all(is.finite(serial$value$table$log_likelihood)))
      expect_identical(profile_pool_without_calls(spread$value),
                       profile_pool_without_calls(serial$value))
      expect_identical(spread$conditions, serial$conditions)
      direct <- lapply(seq_len(nrow(serial$value$table)), function(i) {
        direct_arguments <- c(current,
          list(n_profiles = serial$value$table$n_profiles[i],
               n_group_classes = serial$value$table$n_group_classes[i]))
        fitter <- if (kind == "categorical") multilca else multilpa
        profile_pool_capture(do.call(fitter, direct_arguments))
      })
      expect_identical(lapply(serial$value$fits, profile_pool_without_calls) |>
                         unname(),
                       lapply(direct, function(result) {
                         profile_pool_without_calls(result$value)
                       }))
      if (!single_level) {
        expect_identical(serial$conditions,
                         do.call(c, lapply(direct, `[[`, "conditions")))
      }
    }))
  }))
})

test_that("profile enumeration shares one unnested start pool per search phase", {
  data <- profile_pool_data()
  mapper <- .latents_map_candidates
  active <- 0L
  calls <- list()
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores,
                                                          costs = NULL, ...) {
    active <<- active + 1L
    on.exit(active <<- active - 1L, add = TRUE)
    expect_identical(active, 1L)
    calls[[length(calls) + 1L]] <<- list(tasks = indices, cores = cores,
                                       costs = costs)
    mapper(indices, fun, cores = 1L, costs = costs, ...)
  })
  invisible(lapply(c("categorical", "EEE", "VVI"), function(kind) {
    calls <<- list()
    result <- profile_pool_capture(do.call(enumerate_classes,
      c(profile_pool_arguments(data, kind),
        list(n_profiles = 2:3, n_group_classes = 1:2, cores = 3L))))
    expect_true(all(is.na(result$value$table$error)))
    # The four seeded setups share the pool first, one task per candidate.
    expect_identical(vapply(calls, function(call) length(call$tasks), integer(1)),
                     if (kind == "EEE") c(4L, 16L, 8L) else c(4L, 8L))
    expect_true(all(vapply(calls, `[[`, numeric(1), "cores") <= 3L))
    invisible(lapply(calls, function(call) {
      expect_length(call$costs, length(call$tasks))
      expect_true(all(is.finite(call$costs) & call$costs > 0))
      candidates <- vapply(call$tasks, `[[`, integer(1), "candidate")
      by_candidate <- vapply(split(call$costs, candidates), unique, numeric(1))
      # More profiles or group classes must not receive a smaller cost.
      expect_true(all(diff(by_candidate[c(1L, 2L, 4L)]) > 0))
    }))
  }))
})

test_that("profile start dispatch is longest first without exceeding cores", {
  data <- profile_pool_data()
  mapper <- .latents_map_candidates
  requested <- list()
  dispatched <- list()
  local_mocked_bindings(
    .latents_fork_safe = function(...) TRUE,
    .latents_map_candidates = function(indices, fun, cores, costs = NULL, ...) {
      requested[[length(requested) + 1L]] <<- list(tasks = indices, costs = costs)
      mapper(indices, fun, cores, costs = costs, ...)
    })
  local_mocked_bindings(mclapply = function(X, FUN, mc.cores,
                                            mc.preschedule, ...) {
    dispatched[[length(dispatched) + 1L]] <<- list(tasks = X, cores = mc.cores,
                                                 preschedule = mc.preschedule)
    lapply(X, FUN, ...)
  }, .package = "parallel")
  result <- profile_pool_capture(do.call(enumerate_classes,
    c(profile_pool_arguments(data, "categorical"),
      list(n_profiles = 2:3, n_group_classes = 1:2, cores = 3L))))
  expect_true(all(is.na(result$value$table$error)))
  # One pool for the setups, then one for the starts; each longest first.
  expect_length(requested, 2L)
  expect_length(dispatched, 2L)
  invisible(lapply(seq_along(requested), function(i) {
    expected_order <- order(-requested[[i]]$costs, seq_along(requested[[i]]$costs))
    expect_identical(dispatched[[i]]$tasks, requested[[i]]$tasks[expected_order])
    expect_identical(dispatched[[i]]$cores, 3L)
    expect_false(dispatched[[i]]$preschedule)
  }))
})

test_that("serialized profile worker payloads own caller data and search controls", {
  data <- profile_pool_data()
  worker_names <- c("profile_worker_data", "profile_worker_vars",
                    "profile_worker_categorical", "profile_worker_model",
                    "profile_worker_starts", "profile_worker_screen")
  on.exit(rm(list = intersect(worker_names, ls(globalenv())),
             envir = globalenv()), add = TRUE)
  settings <- c("categorical", "EEE", "VVI")
  expected <- lapply(settings, function(kind) profile_pool_capture(
    do.call(enumerate_classes, c(profile_pool_arguments(data, kind),
      list(n_profiles = 2:3, n_group_classes = 1:2, cores = 1L)))))
  mapper <- .latents_map_candidates
  replayed <- 0L
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    tasks <- unserialize(serialize(indices, NULL))
    shipped <- unserialize(serialize(fun, NULL))
    replayed <<- replayed + 1L
    rm(list = intersect(worker_names, ls(globalenv())), envir = globalenv())
    mapper(tasks, shipped, cores = 1L, ...)
  })
  invisible(lapply(seq_along(settings), function(i) {
    categorical <- settings[i] == "categorical"
    assign("profile_worker_data", data, envir = globalenv())
    assign("profile_worker_vars", if (categorical) c("a", "b", "c") else
      c("x", "y"), envir = globalenv())
    assign("profile_worker_categorical", if (categorical) c("a", "b", "c") else
      character(), envir = globalenv())
    assign("profile_worker_model", if (categorical) NULL else settings[i],
           envir = globalenv())
    assign("profile_worker_starts", 2L, envir = globalenv())
    assign("profile_worker_screen", if (settings[i] == "EEE") 4L else NULL,
           envir = globalenv())
    result <- profile_pool_capture(eval(quote(enumerate_classes(
      profile_worker_data, profile_worker_vars, "group", n_profiles = 2:3,
      n_group_classes = 1:2, categorical = profile_worker_categorical,
      model = profile_worker_model, n_starts = profile_worker_starts,
      n_screen = profile_worker_screen, max_iter = 15L, seed = 12L,
      cores = 3L)), globalenv()))
    expect_identical(profile_pool_without_calls(result$value),
                     profile_pool_without_calls(expected[[i]]$value))
    expect_identical(result$conditions, expected[[i]]$conditions)
  }))
  # Three seeded setup pools, then the categorical and VVI start pools and
  # the EEE screening and final pools.
  expect_identical(replayed, 7L)
})

test_that("membership covariate grids retain candidate parallelism and serial starts", {
  data <- profile_pool_data()
  arguments <- c(profile_pool_arguments(data, "VVI"),
                 list(profile_covariates = "predictor", n_profiles = 2:3,
                      n_group_classes = 1:2))
  expected <- profile_pool_capture(do.call(enumerate_classes,
    c(arguments, list(cores = 1L))))
  mapper <- .latents_map_candidates
  calls <- list()
  active <- 0L
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    # The covariate fitter makes a serial base fit to initialise its logits.
    # A serial mapper creates no workers and may run inside a candidate task.
    if (cores == 1L) return(mapper(indices, fun, cores = 1L, ...))
    active <<- active + 1L
    on.exit(active <<- active - 1L, add = TRUE)
    expect_identical(active, 1L)
    calls[[length(calls) + 1L]] <<- list(tasks = indices, cores = cores)
    mapper(indices, fun, cores = 1L, ...)
  })
  result <- profile_pool_capture(do.call(enumerate_classes,
    c(arguments, list(cores = 3L))))
  expect_identical(profile_pool_without_calls(result$value),
                   profile_pool_without_calls(expected$value))
  expect_identical(result$conditions, expected$conditions)
  expect_true(all(is.na(result$value$table$error)))
  expect_length(calls, 1L)
  expect_length(calls[[1L]]$tasks, 4L)
  expect_identical(calls[[1L]]$cores, 3L)
})

test_that("joined categorical pairs retain exact candidate fits and worker payloads", {
  data <- profile_pool_data()
  arguments <- c(profile_pool_arguments(data, "categorical"),
    list(local_dependence = "a:b", n_profiles = 2:3, n_group_classes = 1:2))
  expected <- profile_pool_capture(do.call(enumerate_classes,
    c(arguments, list(cores = 1L))))
  expect_true(all(is.na(expected$value$table$error)))
  direct <- lapply(seq_len(nrow(expected$value$table)), function(i) {
    current <- arguments
    current$n_profiles <- expected$value$table$n_profiles[i]
    current$n_group_classes <- expected$value$table$n_group_classes[i]
    profile_pool_capture(do.call(multilca, current))
  })
  expect_identical(lapply(expected$value$fits, profile_pool_without_calls) |>
                     unname(),
                   lapply(direct, function(result) {
                     profile_pool_without_calls(result$value)
                   }))
  worker_names <- c("joined_worker_data", "joined_worker_vars", "joined_worker_pairs")
  on.exit(rm(list = intersect(worker_names, ls(globalenv())),
             envir = globalenv()), add = TRUE)
  assign("joined_worker_data", data, envir = globalenv())
  assign("joined_worker_vars", c("a", "b", "c"), envir = globalenv())
  assign("joined_worker_pairs", "a:b", envir = globalenv())
  mapper <- .latents_map_candidates
  parallel_calls <- 0L
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    if (cores == 1L) return(mapper(indices, fun, cores = 1L, ...))
    parallel_calls <<- parallel_calls + 1L
    tasks <- unserialize(serialize(indices, NULL))
    shipped <- unserialize(serialize(fun, NULL))
    expect_length(tasks, 4L)
    expect_true(all(vapply(tasks, function(task) !is.null(task$arguments), logical(1))))
    expect_true(all(vapply(tasks, function(task) (task$arguments$cores %||% 1L) == 1L,
                            logical(1))))
    rm(list = intersect(worker_names, ls(globalenv())), envir = globalenv())
    mapper(tasks, shipped, cores = 1L, ...)
  })
  result <- profile_pool_capture(eval(quote(enumerate_classes(
    joined_worker_data, joined_worker_vars, "group", categorical = joined_worker_vars,
    local_dependence = joined_worker_pairs, n_profiles = 2:3, n_group_classes = 1:2,
    model = NULL, n_starts = 2L, max_iter = 15L, seed = 12L, cores = 3L)),
    globalenv()))
  expect_identical(parallel_calls, 1L)
  expect_identical(profile_pool_without_calls(result$value),
                   profile_pool_without_calls(expected$value))
  expect_identical(result$conditions, expected$conditions)
})

test_that("group-class family enumeration retains its public direct fits", {
  data <- profile_pool_data()
  arguments <- list(data = data, vars = c("x", "y"), id = "group",
                    family = "additive", between_variance = "equal",
                    n_starts = 1L, max_iter = 15L, seed = 12L)
  mapper <- .latents_map_candidates
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    mapper(indices, fun, cores = 1L, ...)
  })
  result <- profile_pool_capture(do.call(enumerate_classes,
    c(arguments, list(n_group_classes = 1:2, cores = 3L))))
  expect_true(all(is.na(result$value$table$error)))
  direct <- lapply(seq_len(2L), function(n_group_classes) profile_pool_capture(
    do.call(multilpa, c(arguments, list(n_group_classes = n_group_classes)))))
  expect_identical(lapply(result$value$fits, profile_pool_without_calls) |> unname(),
                   lapply(direct, function(result) {
                     profile_pool_without_calls(result$value)
                   }))
})

test_that("unseeded profile starts draw in the direct fits' parent grid order", {
  data <- profile_pool_data()
  mapper <- .latents_map_candidates
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    mapper(indices, fun, cores = 1L, ...)
  })
  invisible(lapply(c("categorical", "EEE", "VVI"), function(kind) {
    arguments <- profile_pool_arguments(data, kind)
    arguments["seed"] <- list(NULL)
    enumerated <- function() {
      .latents_local_seed(902L)
      result <- profile_pool_capture(do.call(enumerate_classes,
        c(arguments, list(n_profiles = 2:3, n_group_classes = 1L, cores = 3L))))
      list(result = result, rng = .Random.seed)
    }
    fitted <- function() {
      .latents_local_seed(902L)
      result <- lapply(2:3, function(n_profiles) profile_pool_capture(
        do.call(multilpa, c(arguments,
          list(n_profiles = n_profiles, n_group_classes = 1L)))))
      list(result = result, rng = .Random.seed)
    }
    grid <- enumerated()
    direct <- fitted()
    expect_identical(grid$rng, direct$rng)
    expect_identical(lapply(grid$result$value$fits, profile_pool_without_calls) |>
                       unname(),
                     lapply(direct$result, function(result) {
                       profile_pool_without_calls(result$value)
                     }))
    expect_identical(grid$result$conditions,
                     do.call(c, lapply(direct$result, `[[`, "conditions")))
  }))
})

profile_pool_recovery_grid <- function(data, cores) {
  do.call(enumerate_classes, c(profile_pool_arguments(data, "categorical"),
    list(n_profiles = 2:3, n_group_classes = 1L, cores = cores)))
}

profile_pool_lost <- function(result) {
  Filter(\(condition) inherits(condition, "latents_worker_lost"), result$conditions)
}

profile_pool_kill <- function(parent_pid) {
  stopifnot("only an isolated worker can be killed" = Sys.getpid() != parent_pid)
  tools::pskill(Sys.getpid(), tools::SIGKILL)
}

test_that("a recovered categorical grid start preserves all fitted results", {
  profile_pool_skip_forks()
  data <- profile_pool_data()
  expected <- profile_pool_capture(profile_pool_recovery_grid(data, 1L))
  marker <- tempfile("latents-profile-recovered-")
  on.exit(unlink(marker), add = TRUE)
  parent_pid <- Sys.getpid()
  execute <- .multilpa_execute_start
  local_mocked_bindings(.multilpa_execute_start = function(task) {
    if (task$candidate == 1L && task$start_index == 1L && !file.exists(marker)) {
      writeLines("first attempt", marker)
      profile_pool_kill(parent_pid)
    }
    execute(task)
  })
  result <- profile_pool_capture(profile_pool_recovery_grid(data, 2L))
  expect_identical(profile_pool_without_calls(result$value),
                   profile_pool_without_calls(expected$value))
  lost <- profile_pool_lost(result)
  expect_length(lost, 1L)
  expect_identical(lost[[1L]]$lost, 1L)
  expect_identical(lost[[1L]]$recovered, 1L)
})

test_that("a persistently killed categorical start leaves other starts intact", {
  profile_pool_skip_forks()
  data <- profile_pool_data()
  expected <- profile_pool_capture(profile_pool_recovery_grid(data, 1L))
  parent_pid <- Sys.getpid()
  execute <- .multilpa_execute_start
  local_mocked_bindings(.multilpa_execute_start = function(task) {
    if (task$candidate == 1L && task$start_index == 1L) profile_pool_kill(parent_pid)
    execute(task)
  })
  result <- profile_pool_capture(profile_pool_recovery_grid(data, 2L))
  expect_true(all(is.na(result$value$table$error)))
  starts <- get_results(result$value$fits[[1L]], "starts")
  expect_match(starts$error[1L], "worker process ended", ignore.case = TRUE)
  expect_match(starts$error[1L], "memory", ignore.case = TRUE)
  expect_true(is.na(starts$error[2L]))
  expect_identical(starts$log_likelihood[2L],
                   get_results(expected$value$fits[[1L]], "starts")$log_likelihood[2L])
  expect_identical(result$value$fits[[1L]]$n_failed_starts, 1L)
  expect_identical(result$value$fits[[2L]], expected$value$fits[[2L]])
  expect_length(Filter(\(condition) inherits(condition, "latents_failed_starts"),
                       result$conditions), 1L)
  expect_length(profile_pool_lost(result), 1L)
  expect_identical(profile_pool_lost(result)[[1L]]$remaining, 1L)
})

test_that("losing every categorical start fails only that candidate", {
  profile_pool_skip_forks()
  data <- profile_pool_data()
  expected <- profile_pool_capture(profile_pool_recovery_grid(data, 1L))
  parent_pid <- Sys.getpid()
  execute <- .multilpa_execute_start
  local_mocked_bindings(.multilpa_execute_start = function(task) {
    if (task$candidate == 1L) profile_pool_kill(parent_pid)
    execute(task)
  })
  result <- profile_pool_capture(profile_pool_recovery_grid(data, 2L))
  expect_null(result$value$fits[[1L]])
  expect_match(result$value$table$error[1L], "All 2 starts failed", fixed = TRUE)
  expect_match(result$value$table$error[1L], "worker process ended", ignore.case = TRUE)
  expect_identical(result$value$fits[[2L]], expected$value$fits[[2L]])
  expect_identical(result$value$table[2L, ], expected$value$table[2L, ])
  expect_length(profile_pool_lost(result), 1L)
  expect_identical(profile_pool_lost(result)[[1L]]$remaining, 2L)
})

test_that("direct profile fitting evaluates a supplied cores expression once", {
  data <- profile_pool_data()
  evaluations <- 0L
  result <- profile_pool_capture(multilpa(data, c("x", "y"), "group", 2L, 1L,
    n_starts = 1L, max_iter = 10L, seed = 12L,
    cores = { evaluations <- evaluations + 1L; 1L }))
  expect_identical(evaluations, 1L)
  expect_s3_class(result$value, "multilpa")
})
