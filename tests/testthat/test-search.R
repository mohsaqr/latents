# The start search (R/search.R, `.latents_run_screened()`, the pause and
# resume of `.latents_em()`): emEM screening, the kmeans search it extends,
# and the controls a refit inherits. The release identity of
# `search = "kmeans"` (every field of 36 fits, the RNG state on exit) is
# checked against the recorded 0.96.1 panel by equivalence/search/.

search_toy <- function(rate = 0.9) {
  target <- 1 / (1 - rate)
  list(evaluate = function(state) list(log_likelihood = -sum((state - target)^2)),
       maximize = function(state, expectation) rate * state + 1)
}

search_accelerator <- list(flatten = function(state) state,
                           restore = function(state, theta) theta,
                           valid = function(state) all(is.finite(state)))

test_that("a paused and resumed EM run takes exactly the steps of an unpaused one", {
  toy <- search_toy()
  settings <- list(plain = .latents_em_settings(),
                   squarem = .latents_em_settings(accelerate = search_accelerator))
  invisible(lapply(settings, function(setting) {
    straight <- .latents_em(c(0, 3), toy$evaluate, toy$maximize, max_iter = 40,
                            tol = 1e-14, settings = setting)
    invisible(lapply(c(0, 1, 2, 3, 4, 5, 7, 11, 39, 40, 100), function(pause) {
      paused <- .latents_em(c(0, 3), toy$evaluate, toy$maximize, max_iter = 40,
                            tol = 1e-14, settings = setting, pause = pause)
      expect_lte(paused$iterations, 40L)
      resumed <- .latents_em(NULL, toy$evaluate, toy$maximize, max_iter = 40,
                             tol = 1e-14, settings = setting,
                             resume = paused$checkpoint)
      expect_identical(resumed$state, straight$state)
      expect_identical(resumed$history, straight$history)
      expect_identical(resumed$iterations, straight$iterations)
      expect_identical(resumed$converged, straight$converged)
    }))
  }))
})

test_that("a run ended by a non-finite step stays ended, with the E-step before it", {
  evaluate <- function(state) list(log_likelihood = if (state > 2) NaN else -abs(state - 5))
  maximize <- function(state, expectation) state + 1
  setting <- .latents_em_settings(stop_on_nonfinite = TRUE)
  straight <- .latents_em(0, evaluate, maximize, max_iter = 10, tol = 1e-12,
                          settings = setting)
  paused <- .latents_em(0, evaluate, maximize, max_iter = 10, tol = 1e-12,
                        settings = setting, pause = 5)
  resumed <- .latents_em(NULL, evaluate, maximize, max_iter = 10, tol = 1e-12,
                         settings = setting, resume = paused$checkpoint)
  expect_identical(resumed$state, straight$state)
  expect_identical(resumed$expectation, straight$expectation)
  expect_identical(resumed$iterations, straight$iterations)
})

search_runner <- function(screened_values, final_values, fail = integer()) {
  # A run whose screened likelihood is `screened_values[i]` and whose final
  # one is `final_values[i]`; the runner sees only what a real engine returns.
  list(
    prepare = function(start_index) start_index,
    run = function(start_index, values, pause = Inf, resume = NULL) {
      if (start_index %in% fail) stop("start ", start_index, " failed")
      index <- if (is.null(resume)) values else resume$state
      done <- !is.null(resume)
      list(expectation = list(log_likelihood = if (done) final_values[index] else
             screened_values[index]),
           converged = done, iterations = if (done) 30L else 5L,
           checkpoint = list(state = index))
    })
}

test_that("screening runs the protected start and the best screened starts to the end", {
  toy <- search_runner(screened_values = c(-50, -10, -30, -10, -20, -40),
                       final_values = c(-9, -5, -8, -6, -7, -1))
  found <- .latents_run_starts(3, toy$run, "likelihood", prepare = toy$prepare,
                               screen = list(n_pool = 6, iterations = 5, protect = 1L))
  # Start 1 is protected; then 2 and 4 (tied at -10, pool order), not 6.
  expect_identical(found$stage, c("final", "final", "screened", "final",
                                  "screened", "screened"))
  expect_identical(found$best_start, 2L)
  expect_identical(found$screened_scores, c(-50, -10, -30, -10, -20, -40))
  expect_identical(found$scores[c(1, 2, 4)], c(-9, -5, -6))
  # A start that only screened keeps its screened likelihood and is not chosen
  # even when that is higher than every final one.
  high <- search_runner(screened_values = c(-50, -10, -30, -10, -2, -40),
                        final_values = c(-9, -5, -8, -6, -7, -1))
  chosen <- .latents_run_starts(2, high$run, "likelihood", prepare = high$prepare,
                                screen = list(n_pool = 6, iterations = 5,
                                              protect = integer()))
  expect_identical(which(chosen$stage == "final"), c(2L, 5L))
  expect_identical(chosen$best_start, 2L)
})

test_that("failed starts are recorded, and the fit is refused only when all fail", {
  toy <- search_runner(screened_values = c(-50, -10, -30, -10), final_values = c(-9, -5, -8, -6),
                       fail = c(2L, 3L))
  found <- .latents_run_starts(2, toy$run, "likelihood", prepare = toy$prepare,
                               screen = list(n_pool = 4, iterations = 5, protect = 1L))
  expect_identical(found$valid, c(TRUE, FALSE, FALSE, TRUE))
  expect_identical(which(found$stage == "final"), c(1L, 4L))
  none <- search_runner(c(-1, -2), c(-1, -2), fail = 1:2)
  expect_error(.latents_run_starts(1, none$run, "likelihood", prepare = none$prepare,
                                   screen = list(n_pool = 2, iterations = 5,
                                                 protect = 1L)),
               class = "latents_all_starts_failed")
})

test_that("random starts assign every row and leave no cluster empty", {
  set.seed(3)
  input <- matrix(rnorm(40), 20, 2)
  partition <- .multilpa_random_assignments(input, 3, "partition")
  expect_length(partition, 20L)
  expect_identical(sort(unique(partition)), 1:3)
  # A random partition deals the rows out evenly.
  expect_true(all(tabulate(partition, 3) %in% c(6L, 7L)))
  # Repeated rows cannot empty a cluster: each drawn row keeps its own.
  repeated <- rbind(input[1:5, ], input[1:5, ], input[1:5, ], input[1:5, ])
  invisible(lapply(seq_len(25), function(i) {
    centres <- .multilpa_random_assignments(repeated, 4, "centres")
    expect_identical(sort(unique(centres)), 1:4)
  }))
})

test_that("emEM resumes the kmeans search's own starts exactly", {
  items <- names(srl)
  kmeans <- quietly(lpa(srl, items, 4, model = "VVV", seed = 2, search = "kmeans"))
  emem <- quietly(lpa(srl, items, 4, model = "VVV", seed = 2, search = "emEM"))
  kmeans_starts <- get_results(kmeans, "starts")
  emem_starts <- get_results(emem, "starts")
  expect_identical(nrow(emem_starts), 50L)
  expect_identical(sum(emem_starts$stage == "final"), 10L)
  resumed <- subset(emem_starts, stage == "final" & start <= 10)
  expect_gt(nrow(resumed), 0L)
  # Each start the two searches share ends at the same likelihood, bit for
  # bit: screening pauses it at a cycle boundary and the final stage resumes it.
  expect_identical(resumed$log_likelihood, kmeans_starts$log_likelihood[resumed$start])
  expect_identical(emem_starts$kind[1:3], c("ward", "kmeans", "kmeans"))
  expect_setequal(unique(emem_starts$kind[11:50]), c("partition", "centres"))
  expect_gte(emem$log_likelihood, max(kmeans_starts$log_likelihood[resumed$start]))
})

test_that("the kmeans search keeps its table and records its controls", {
  fit <- quietly(lpa(srl, names(srl), 2, seed = 1, search = "kmeans"))
  expect_identical(names(get_results(fit, "starts")),
                   c("start", "log_likelihood", "converged", "iterations", "error",
                     "boundary"))
  expect_identical(fit$search, "kmeans")
  expect_identical(fit$n_starts, 10L)
  expect_null(fit$n_screen)
  expect_identical(fit$select_start, "likelihood")
  expect_identical(fit$boundary_policy, "allow")
  # The kmeans search is the default; emEM is asked for by name.
  default <- quietly(lpa(srl, names(srl), 2, seed = 1))
  expect_identical(default$search, "kmeans")
  expect_identical(get_results(default, "starts"), get_results(fit, "starts"))
})

test_that("emEM is reproducible under a seed and leaves the caller's stream alone", {
  items <- names(srl)
  first <- quietly(lpa(srl, items, 3, seed = 9, search = "emEM", n_starts = 4))
  set.seed(5)
  before <- .Random.seed
  second <- quietly(lpa(srl, items, 3, seed = 9, search = "emEM", n_starts = 4))
  expect_identical(.Random.seed, before)
  expect_identical(second$log_likelihood, first$log_likelihood)
  expect_identical(get_results(second, "starts"), get_results(first, "starts"))
  expect_identical(second$n_screen, 20L)
  smaller <- quietly(lpa(srl, items, 3, seed = 9, search = "emEM", n_starts = 4,
                         n_screen = 8))
  expect_identical(nrow(get_results(smaller, "starts")), 8L)
})

test_that("emEM gives the same fit on one core and on two", {
  skip_on_cran()
  items <- names(srl)
  one <- quietly(lpa(srl, items, 3, seed = 4, search = "emEM", n_starts = 4))
  two <- quietly(lpa(srl, items, 3, seed = 4, search = "emEM", n_starts = 4, cores = 2))
  expect_identical(two$log_likelihood, one$log_likelihood)
  expect_identical(get_results(two, "starts"), get_results(one, "starts"))
})

test_that("a supplied start always runs to the end, and alone when it is the only start", {
  items <- names(srl)
  first <- quietly(lpa(srl, items, 2, seed = 1, n_starts = 2))
  alone <- quietly(lpa(srl, items, 2, start = starting_values(first), n_starts = 1,
                       search = "emEM"))
  expect_identical(nrow(get_results(alone, "starts")), 1L)
  expect_null(alone$n_screen)
  pooled <- quietly(lpa(srl, items, 2, start = starting_values(first), n_starts = 3,
                        search = "emEM", seed = 2))
  starts <- get_results(pooled, "starts")
  expect_identical(starts$kind[1], "supplied")
  expect_identical(starts$stage[1], "final")
})

test_that("invalid search controls are refused by class", {
  items <- names(srl)
  expect_error(lpa(srl, items, 2, search = "emEM", n_screen = 5),
               class = "latents_bad_argument")
  expect_error(lpa(srl, items, 2, search = "emEM", n_screen = 20.5),
               class = "latents_bad_argument")
  expect_error(lpa(srl, items, 2, n_screen = 50, search = "kmeans"),
               class = "latents_bad_argument")
})

test_that("every fit path runs the emEM search and records it", {
  engagement <- c("browse", "lectures", "forum_read")
  fits <- quietly(list(
    covariates = multilpa(course_engagement, engagement, "course", 2, 2,
                          profile_covariates = "previous_grade", seed = 1,
                          n_starts = 2, search = "emEM"),
    additive = multilpa(course_engagement, engagement, "course", n_group_classes = 2,
                        family = "additive", seed = 1, n_starts = 2, search = "emEM"),
    full_cross_level = multilpa(course_engagement, engagement, "course", 2, 2,
                                family = "full_cross_level", seed = 1, n_starts = 2,
                                search = "emEM"),
    transitions = lta(course_engagement, engagement, "student", 2, time = "sequence",
                      seed = 1, n_starts = 2, search = "emEM"),
    general_transitions = lta(course_engagement, engagement, "student", 2,
                              time = "sequence", mover_stayer = TRUE, seed = 1,
                              n_starts = 2, search = "emEM"),
    staged = fit_staged(course_engagement, engagement, "course", n_profiles = 2,
                        n_group_classes = 2, seed = 1, n_starts = 2, search = "emEM")),
    # The indicators are standardised within each course, so every course mean
    # is 0 and two additive classes of them coincide.
    classes = c(.multilpa_expected_warnings, "latents_weak_class"))
  invisible(lapply(names(fits), function(path) {
    fit <- fits[[path]]
    expect_identical(fit$search, "emEM", info = path)
    expect_identical(fit$n_screen, 10L, info = path)
    starts <- get_results(fit, "starts")
    expect_identical(nrow(starts), 10L, info = path)
    expect_identical(sum(starts$stage == "final"), 2L, info = path)
  }))
})

test_that("the bootstrap refits each model with the controls it was fitted with", {
  items <- names(srl)
  one <- quietly(lpa(srl, items, 1, model = "EEE", seed = 1, search = "emEM",
                     n_starts = 3, tol = 1e-7))
  two <- quietly(lpa(srl, items, 2, model = "EEE", seed = 1, search = "emEM",
                     n_starts = 3, tol = 1e-7, acceleration = "none",
                     select_start = "converged"))
  calls <- list()
  real <- multilpa
  local_mocked_bindings(multilpa = function(...) {
    arguments <- list(...)
    calls[[length(calls) + 1L]] <<- arguments[c("n_starts", "n_screen", "search",
                                                 "tol", "max_iter", "acceleration",
                                                 "select_start", "boundary")]
    real(...)
  })
  result <- quietly(bootstrap_lrt(one, two, iter = 2, seed = 3))
  expect_length(calls, 4L)
  null_calls <- calls[c(1, 3)]
  alternative_calls <- calls[c(2, 4)]
  invisible(lapply(null_calls, function(arguments) {
    expect_identical(arguments$acceleration, "squarem")
    expect_identical(arguments$select_start, "likelihood")
  }))
  invisible(lapply(alternative_calls, function(arguments) {
    expect_identical(arguments$search, "emEM")
    expect_identical(arguments$n_starts, 3L)
    expect_identical(arguments$n_screen, 15L)
    expect_identical(arguments$tol, 1e-7)
    expect_identical(arguments$acceleration, "none")
    expect_identical(arguments$select_start, "converged")
    expect_identical(arguments$boundary, "allow")
  }))
  controls <- result$refit_controls
  expect_identical(controls$role, c("null", "alternative"))
  expect_identical(controls$acceleration, c("squarem", "none"))
  # An explicit number of starts applies to both models, with the default pool.
  calls <- list()
  overridden <- quietly(bootstrap_lrt(one, two, iter = 2, seed = 3, n_starts = 2))
  expect_identical(overridden$refit_controls$n_starts, c(2L, 2L))
  expect_identical(unique(vapply(calls, function(arguments) arguments$n_starts,
                                 integer(1))), 2L)
  expect_true(all(vapply(calls, function(arguments) is.null(arguments$n_screen),
                         logical(1))))
})

test_that("a covariate fit records the stopping rule its bootstrap refits inherit", {
  set.seed(5)
  data <- data.frame(g = rep(1:30, each = 8), a = stats::rnorm(240), b = stats::rnorm(240),
                     z = stats::rnorm(240))
  fit <- quietly(multilpa(data, c("a", "b"), "g", 2, 1, profile_covariates = "z",
                          n_starts = 2, seed = 1, max_iter = 777L, tol = 1e-9))
  controls <- .latents_refit_controls(fit)
  expect_identical(controls$max_iter, 777L)
  expect_identical(controls$tol, 1e-9)
})

test_that("a staged fit's bootstrap refits keep its screening pool", {
  set.seed(6)
  data <- data.frame(g = rep(1:30, each = 8), a = stats::rnorm(240), b = stats::rnorm(240))
  staged <- quietly(fit_staged(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                               n_starts = 2, n_screen = 8, seed = 1, search = "emEM"))
  screens <- list()
  real <- fit_staged
  local_mocked_bindings(fit_staged = function(...) {
    screens[[length(screens) + 1L]] <<- list(...)$n_screen
    real(...)
  })
  quietly(.multilpa_bootstrap_replicates(staged, data, iter = 1, n_starts = 2,
                                         max_iter = 200L, tol = 1e-6))
  expect_identical(screens, list(8L))
})

test_that("summaries count replications out of the starts run to the end", {
  fit <- quietly(lpa(srl, names(srl), 2, n_starts = 2, seed = 1, search = "emEM"))
  expect_identical(summary(fit)$n_starts, 2L)
  expect_match(capture.output(print(summary(fit))), "replicated in [0-9]+/2 starts",
               all = FALSE)
  data <- subset(course_engagement, student <= 30)
  transitions <- quietly(lta(data, c("browse", "lectures"), "student", 2, time = "sequence",
                             n_starts = 2, seed = 1, search = "emEM"))
  expect_match(capture.output(print(summary(transitions))), "replicated in [0-9]+/2 starts",
               all = FALSE)
})

test_that("a fit made before the controls were recorded is refitted as before", {
  fit <- quietly(lpa(srl, names(srl), 2, seed = 1, max_iter = 500, tol = 1e-6))
  fit[c("search", "n_starts", "n_screen", "select_start", "boundary_policy")] <- NULL
  controls <- .latents_refit_controls(fit)
  expect_identical(controls[c("search", "n_starts", "max_iter", "tol", "select_start",
                              "boundary")],
                   list(search = "kmeans", n_starts = 10L, max_iter = 1000L,
                        tol = 1e-8, select_start = "likelihood", boundary = "allow"))
})
