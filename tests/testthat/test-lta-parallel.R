# These small panels exercise both transition engines without a large fixture.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

lta_parallel_panel <- function() {
  .latents_local_seed(8401L)
  panel <- data.frame(person = rep(seq_len(12L), each = 5L),
                      beep = rep(seq_len(5L), 12L))
  state <- (panel$person + panel$beep) %% 2L
  panel$a <- 2 * state + stats::rnorm(nrow(panel), sd = 0.6)
  panel$b <- -1.5 * state + stats::rnorm(nrow(panel), sd = 0.7)
  panel
}

lta_parallel_conditions <- function(expr) {
  conditions <- character()
  value <- withCallingHandlers(expr,
    warning = function(condition) {
      conditions <<- c(conditions, conditionMessage(condition))
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      # Progress messages report the worker count; they are not part of the
      # qualifications compared between serial and parallel fits.
      if (!inherits(condition, "latents_progress")) {
        conditions <<- c(conditions, conditionMessage(condition))
      }
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = conditions)
}

lta_parallel_without_calls <- function(value) {
  value$call <- NULL
  if (!is.null(value$fits)) {
    value$fits <- lapply(value$fits, function(fit) {
      if (!is.null(fit)) fit$call <- NULL
      fit
    })
  }
  value
}

test_that("lta appends and validates cores without moving existing arguments", {
  expect_identical(tail(names(formals(lta)), 2L), c("acceleration", "cores"))
  expect_identical(formals(lta)$cores, 1L)
  panel <- lta_parallel_panel()
  invisible(lapply(list(0, -1L, 1.5, NA_real_, Inf, "2", c(1L, 2L)),
    function(cores) {
      expect_error(lta(panel, c("a", "b"), "person", 2L, time = "beep",
                       n_starts = 1L, cores = cores),
                   class = "latents_bad_argument")
    }))
})

test_that("both lta engines have identical forked starts and emEM screening", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(!.latents_fork_safe(), "This BLAS cannot run safely in forked workers")
  panel <- lta_parallel_panel()
  settings <- expand.grid(transitions = c("homogeneous", "occasion"),
                          search = c("kmeans", "emEM"),
                          stringsAsFactors = FALSE)
  invisible(lapply(seq_len(nrow(settings)), function(i) {
    fit <- function(cores) lta_parallel_conditions(lta(
      panel, c("a", "b"), "person", 2L, time = "beep",
      n_group_classes = 2L, n_starts = 3L, max_iter = 25L, seed = 12L,
      model = "VVI", transitions = settings$transitions[i],
      search = settings$search[i],
      n_screen = if (settings$search[i] == "emEM") 6L else NULL,
      cores = cores))
    serial <- fit(1L)
    spread <- fit(2L)
    expect_identical(lta_parallel_without_calls(spread$value),
                     lta_parallel_without_calls(serial$value))
    expect_identical(spread$conditions, serial$conditions)
    expect_identical(get_results(spread$value, "starts"),
                     get_results(serial$value, "starts"))
  }))
})

test_that("transition grids match serial and direct fits for both searches and engines", {
  skip_on_cran()
  skip_if_cores_limited()
  skip_on_os("windows")
  skip_if(!.latents_fork_safe(), "This BLAS cannot run safely in forked workers")
  panel <- lta_parallel_panel()
  settings <- expand.grid(transitions = c("homogeneous", "occasion"),
                          search = c("kmeans", "emEM"),
                          stringsAsFactors = FALSE)
  invisible(lapply(seq_len(nrow(settings)), function(i) {
    fit <- function(cores) lta_parallel_conditions(enumerate_classes(
      panel, c("a", "b"), "person", time = "beep", n_profiles = 2:3,
      n_group_classes = 1:2, model = "VVI", n_starts = 2L, max_iter = 20L,
      seed = 12L, transitions = settings$transitions[i],
      search = settings$search[i],
      n_screen = if (settings$search[i] == "emEM") 4L else NULL, cores = cores))
    serial <- fit(1L)
    expect_true(all(is.finite(serial$value$table$log_likelihood)))
    invisible(lapply(c(2L, 3L), function(cores) {
      spread <- fit(cores)
      expect_identical(lta_parallel_without_calls(spread$value),
                       lta_parallel_without_calls(serial$value))
      expect_identical(spread$conditions, serial$conditions)
      expect_identical(as.data.frame(spread$value), as.data.frame(serial$value))
    }))
    grid <- serial$value$table
    direct <- lapply(seq_len(nrow(grid)), function(j) lta_parallel_conditions(lta(
      panel, c("a", "b"), "person", grid$n_profiles[j], time = "beep",
      n_group_classes = grid$n_group_classes[j], model = "VVI", n_starts = 2L,
      max_iter = 20L, seed = 12L, transitions = settings$transitions[i],
      search = settings$search[i],
      n_screen = if (settings$search[i] == "emEM") 4L else NULL, cores = 1L)))
    expect_identical(lapply(serial$value$fits, lta_parallel_without_calls) |>
                       unname(),
                     lapply(direct, function(result) {
                       lta_parallel_without_calls(result$value)
                     }))
    expect_identical(serial$conditions,
                     unlist(lapply(direct, `[[`, "conditions"),
                            use.names = FALSE) |> as.character())
  }))
})

test_that("serialized transition workers carry caller data and start controls", {
  # Replay the actual task closure in a session lacking the caller's globals.
  # Keep the real mapper's condition envelope, including signal = FALSE.
  worker_names <- c("lta_worker_panel", "lta_worker_vars", "lta_worker_id",
                    "lta_worker_time", "lta_worker_starts", "lta_worker_search",
                    "lta_worker_transitions")
  on.exit(rm(list = intersect(worker_names, ls(globalenv())),
             envir = globalenv()), add = TRUE)
  panel <- lta_parallel_panel()
  settings <- expand.grid(transitions = c("homogeneous", "occasion"),
                          search = c("kmeans", "emEM"),
                          stringsAsFactors = FALSE)
  serial <- lapply(seq_len(nrow(settings)), function(i) lta_parallel_conditions(
    enumerate_classes(panel, c("a", "b"), "person", time = "beep",
      n_profiles = 2:3, n_group_classes = 1:2, n_starts = 2L,
      search = settings$search[i], transitions = settings$transitions[i],
      max_iter = 20L, seed = 12L, cores = 1L)))
  mapper <- .latents_map_candidates
  worker_view <- function(indices, fun, cores, ...) {
    shipped_indices <- unserialize(serialize(indices, NULL))
    shipped <- unserialize(serialize(fun, NULL))
    rm(list = intersect(worker_names, ls(globalenv())), envir = globalenv())
    mapper(shipped_indices, shipped, cores = 1L, ...)
  }
  local_mocked_bindings(.latents_map_candidates = worker_view)
  invisible(lapply(seq_len(nrow(settings)), function(i) {
    assign("lta_worker_panel", panel, envir = globalenv())
    assign("lta_worker_vars", c("a", "b"), envir = globalenv())
    assign("lta_worker_id", "person", envir = globalenv())
    assign("lta_worker_time", "beep", envir = globalenv())
    assign("lta_worker_starts", 2L, envir = globalenv())
    assign("lta_worker_search", settings$search[i], envir = globalenv())
    assign("lta_worker_transitions", settings$transitions[i], envir = globalenv())
    result <- lta_parallel_conditions(eval(quote(enumerate_classes(
      lta_worker_panel, lta_worker_vars, lta_worker_id,
      time = lta_worker_time, n_profiles = 2:3, n_group_classes = 1:2,
      n_starts = lta_worker_starts, search = lta_worker_search,
      transitions = lta_worker_transitions,
      max_iter = 20L, seed = 12L, cores = 2L)), globalenv()))
    expect_true(all(is.na(result$value$table$error)))
    expect_true(all(is.finite(result$value$table$log_likelihood)))
    expect_length(result$value$fits, 4L)
    expect_identical(lta_parallel_without_calls(result$value),
                     lta_parallel_without_calls(serial[[i]]$value))
    expect_identical(result$conditions, serial[[i]]$conditions)
  }))
})

test_that("flat task dispatch is longest first and respects the worker budget", {
  indices <- c(17L, 11L, 8L, 6L, 4L, 3L, 2L)
  costs <- c(1, 24, 3, 24, 12, 2, 48)
  requested <- list()
  local_mocked_bindings(mclapply = function(X, FUN, mc.cores, mc.preschedule, ...) {
    requested[[length(requested) + 1L]] <<- list(indices = X, cores = mc.cores,
                                               preschedule = mc.preschedule)
    lapply(X, FUN, ...)
  }, .package = "parallel")
  invisible(lapply(c(2L, 3L, 4L, 8L, 56L), function(budget) {
    result <- .latents_map_candidates(indices, identity, budget, backend = "fork",
                                      costs = costs)
    expect_identical(unlist(result, use.names = FALSE), indices)
    dispatch <- tail(requested, 1L)[[1L]]
    expect_identical(dispatch$indices, indices[order(-costs, seq_along(costs))])
    expect_identical(dispatch$cores, min(budget, length(indices)))
    expect_false(dispatch$preschedule)
  }))
})

test_that("transition enumeration uses one unnested pool per search phase", {
  panel <- lta_parallel_panel()
  mapper <- .latents_map_candidates
  calls <- list()
  active <- 0L
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores,
                                                          costs = NULL, ...) {
    active <<- active + 1L
    on.exit(active <<- active - 1L, add = TRUE)
    expect_identical(active, 1L)
    calls[[length(calls) + 1L]] <<- list(tasks = length(indices), cores = cores,
                                       costs = costs)
    mapper(indices, fun, cores = 1L, costs = costs, ...)
  })
  invisible(lapply(c("kmeans", "emEM"), function(search) {
    calls <<- list()
    result <- lta_parallel_conditions(enumerate_classes(
      panel, c("a", "b"), "person", time = "beep", n_profiles = 2:3,
      n_group_classes = 1:2, n_starts = 2L, search = search,
      n_screen = if (search == "emEM") 4L else NULL,
      max_iter = 20L, seed = 12L, cores = 3L))
    expect_true(all(is.na(result$value$table$error)))
    # The four seeded setups share the pool first, one task per candidate.
    expect_identical(vapply(calls, `[[`, integer(1), "tasks"),
                     if (search == "emEM") c(4L, 16L, 8L) else c(4L, 8L))
    expect_true(all(vapply(calls, `[[`, numeric(1), "cores") <= 3L))
    costs <- result$value$table$n_group_classes * result$value$table$n_profiles^2
    expected <- if (search == "emEM") c(1L, 4L, 2L) else c(1L, 2L)
    invisible(lapply(seq_along(calls), function(i) {
      expect_identical(sort(calls[[i]]$costs), sort(rep(costs, each = expected[i])))
    }))
  }))
})

test_that("unseeded transition setup draws in the same parent grid order as lta", {
  panel <- lta_parallel_panel()
  mapper <- .latents_map_candidates
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    mapper(indices, fun, cores = 1L, ...)
  })
  invisible(lapply(c("kmeans", "emEM"), function(search) {
    grid <- function() {
      .latents_local_seed(902L)
      result <- lta_parallel_conditions(enumerate_classes(
        panel, c("a", "b"), "person", time = "beep", n_profiles = 2:3,
        n_group_classes = 1L, n_starts = 2L, search = search,
        n_screen = if (search == "emEM") 4L else NULL,
        max_iter = 12L, seed = NULL, cores = 3L))
      list(result = result, rng = .Random.seed)
    }
    direct <- function() {
      .latents_local_seed(902L)
      result <- lapply(2:3, function(n_profiles) lta_parallel_conditions(lta(
        panel, c("a", "b"), "person", n_profiles, time = "beep",
        n_group_classes = 1L, n_starts = 2L, search = search,
        n_screen = if (search == "emEM") 4L else NULL,
        max_iter = 12L, seed = NULL)))
      list(result = result, rng = .Random.seed)
    }
    enumerated <- grid()
    fitted <- direct()
    expect_identical(enumerated$rng, fitted$rng)
    expect_identical(lapply(enumerated$result$value$fits,
                            lta_parallel_without_calls) |> unname(),
                     lapply(fitted$result, function(result) {
                       lta_parallel_without_calls(result$value)
                     }))
    expect_identical(enumerated$result$conditions,
                     unlist(lapply(fitted$result, `[[`, "conditions"),
                            use.names = FALSE) |> as.character())
  }))
})

test_that("longest-first tasks replay classed conditions in input order", {
  seen <- list()
  result <- withCallingHandlers(.latents_map_candidates(
    seq_len(7L), function(i) {
      warning(warningCondition(paste("candidate", i),
                               class = "lta_parallel_test_warning"))
      i
    }, cores = 1L, costs = seq_len(7L)^2),
    warning = function(condition) {
      seen[[length(seen) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    })
  expect_identical(result, as.list(seq_len(7L)))
  expect_identical(vapply(seen, conditionMessage, character(1)),
                   paste("candidate", seq_len(7L)))
  expect_true(all(vapply(seen, inherits, logical(1), "lta_parallel_test_warning")))
})

test_that("covariance qualifications count every screened and continued start once", {
  panel <- lta_parallel_panel()
  original_em <- .multilpa_transition_em
  mapper <- .latents_map_candidates
  local_mocked_bindings(
    .multilpa_transition_em = function(...) {
      .multilpa_warn_structure("VVI", 7L)
      original_em(...)
    },
    .latents_map_candidates = function(indices, fun, cores, ...) {
      mapper(indices, fun, cores = 1L, ...)
    })
  invisible(lapply(c("kmeans", "emEM"), function(search) {
    arguments <- list(data = panel, vars = c("a", "b"), id = "person",
                      time = "beep", n_profiles = 2L, n_group_classes = 1L,
                      model = "VVI", n_starts = 2L, max_iter = 2L,
                      search = search,
                      n_screen = if (search == "emEM") 4L else NULL, seed = 12L)
    direct <- lta_parallel_conditions(do.call(lta, arguments))
    enumerated <- lta_parallel_conditions(do.call(enumerate_classes,
                                                  c(arguments, list(cores = 3L))))
    qualification <- direct$conditions[grepl("M-step(s)", direct$conditions,
                                             fixed = TRUE)]
    expect_length(qualification, 1L)
    expected <- if (search == "emEM") 6L else 2L
    expect_match(qualification, sprintf("%d M-step(s) of the VVI", expected),
                 fixed = TRUE)
    expect_identical(enumerated$conditions, direct$conditions)
    expect_identical(lta_parallel_without_calls(enumerated$value$fits[[1L]]),
                     lta_parallel_without_calls(direct$value))
  }))
})

test_that("covariance qualifications survive failure of every transition start", {
  panel <- lta_parallel_panel()
  mapper <- .latents_map_candidates
  local_mocked_bindings(
    .multilpa_transition_em = function(...) {
      .multilpa_warn_structure("VVI", 7L)
      stop(errorCondition("deliberate failure",
                           class = "lta_parallel_deliberate_failure"))
    },
    .latents_map_candidates = function(indices, fun, cores, ...) {
      mapper(indices, fun, cores = 1L, ...)
    })
  invisible(lapply(c("kmeans", "emEM"), function(search) {
    arguments <- list(data = panel, vars = c("a", "b"), id = "person",
                      time = "beep", n_profiles = 2L, n_group_classes = 1L,
                      model = "VVI", n_starts = 2L, max_iter = 2L,
                      search = search,
                      n_screen = if (search == "emEM") 4L else NULL, seed = 12L)
    direct <- lta_parallel_conditions(tryCatch(do.call(lta, arguments),
                                               error = identity))
    enumerated <- lta_parallel_conditions(do.call(enumerate_classes,
                                                  c(arguments, list(cores = 3L))))
    expect_s3_class(direct$value, "latents_all_starts_failed")
    expect_length(direct$conditions, 1L)
    expected <- if (search == "emEM") 4L else 2L
    expect_match(direct$conditions, sprintf("%d M-step(s) of the VVI", expected),
                 fixed = TRUE)
    expect_identical(enumerated$conditions, direct$conditions)
    expect_identical(enumerated$value$table$error, conditionMessage(direct$value))
    expect_null(enumerated$value$fits[[1L]])
  }))
})

test_that("transition enumeration propagates a lost flat start worker", {
  local_mocked_bindings(.latents_map_candidates = function(...) {
    stop(errorCondition("A flat start worker was lost.",
                         class = "latents_parallel_failure", call = NULL))
  })
  panel <- lta_parallel_panel()
  invisible(lapply(c("homogeneous", "occasion"), function(transitions) {
    expect_error(enumerate_classes(
      panel, c("a", "b"), "person", time = "beep", n_profiles = 2L,
      n_group_classes = 1L, n_starts = 2L, seed = 12L, cores = 2L,
      transitions = transitions), "A flat start worker was lost",
      class = "latents_parallel_failure")
  }))
})

test_that("direct general lta workers serialize both emEM stages without caller globals", {
  panel <- lta_parallel_panel()
  serial <- lta_parallel_conditions(lta(
    panel, c("a", "b"), "person", 2L, time = "beep", n_starts = 2L,
    transitions = "occasion", search = "emEM", n_screen = 4L,
    max_iter = 10L, seed = 12L, cores = 1L))
  globals <- list(lta_general_worker_panel = panel,
                   lta_general_worker_vars = c("a", "b"),
                   lta_general_worker_id = "person",
                   lta_general_worker_time = "beep",
                   lta_general_worker_starts = 2L,
                   lta_general_worker_transitions = "occasion",
                   lta_general_worker_search = "emEM",
                   lta_general_worker_screen = 4L,
                   lta_general_worker_iter = 10L)
  invisible(list2env(globals, envir = globalenv()))
  worker_names <- names(globals)
  on.exit(rm(list = intersect(worker_names, ls(globalenv())),
             envir = globalenv()), add = TRUE)
  calls <- integer()
  mapper <- .latents_map_candidates
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    calls <<- c(calls, cores)
    shipped_indices <- unserialize(serialize(indices, NULL))
    shipped <- unserialize(serialize(fun, NULL))
    rm(list = intersect(worker_names, ls(globalenv())), envir = globalenv())
    mapper(shipped_indices, shipped, cores = 1L, ...)
  })
  spread <- lta_parallel_conditions(eval(quote(lta(
    lta_general_worker_panel, lta_general_worker_vars, lta_general_worker_id,
    2L, time = lta_general_worker_time, n_starts = lta_general_worker_starts,
    transitions = lta_general_worker_transitions, search = lta_general_worker_search,
    n_screen = lta_general_worker_screen, max_iter = lta_general_worker_iter,
    seed = 12L, cores = 2L)), globalenv()))
  expect_identical(calls, c(2L, 2L))
  expect_identical(lta_parallel_without_calls(spread$value),
                   lta_parallel_without_calls(serial$value))
  expect_identical(spread$conditions, serial$conditions)
})
