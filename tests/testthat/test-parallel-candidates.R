test_that("lost profile and family candidates use the ordinary failed-candidate row", {
  invisible(lapply(c("profiles", "additive"), function(family) {
    original_fit <- multilpa
    fitting <- integer()
    lose_candidate <- FALSE
    loss <- paste("The worker process ended without returning a result,",
                  "likely killed for memory.")
    worker_data <- data.frame(
      person = rep(seq_len(12L), each = 6L),
      y = sin(seq_len(72L)) + rep(seq_len(12L) / 8, each = 6L),
      z = cos(seq_len(72L)) + rep(seq_len(12L) / 10, each = 6L))
    assign("candidate_loss_data", worker_data, envir = globalenv())
    on.exit(rm(list = intersect("candidate_loss_data", ls(envir = globalenv())),
               envir = globalenv()))
    fit_or_error <- function() {
      fitting <<- c(fitting, n_group_classes)
      if (n_group_classes == 2L) {
        stop(errorCondition(loss, class = "candidate_test_failure"))
      }
      argument_names <- names(as.list(match.call())[-1L])
      arguments <- lapply(argument_names, get, envir = environment(),
                          inherits = FALSE)
      names(arguments) <- argument_names
      do.call(original_fit, arguments)
    }
    # Enumeration validates arguments against the fitting function's formals.
    formals(fit_or_error) <- formals(original_fit)
    local_mocked_bindings(
      multilpa = fit_or_error,
      # Exercise profile candidate fallback here; individual flat-start loss
      # is covered by test-profile-pool.R.
      .multilpa_pool_fallback = function(arguments) TRUE,
      .latents_map_candidates = function(indices, fun, cores, on_lost = NULL,
                                        signal = TRUE, ...) {
        if (!lose_candidate) {
          values <- lapply(indices, function(i) .latents_collect_conditions(function() fun(i)))
          return(if (signal) lapply(values, `[[`, "value") else values)
        }
        # Replay the worker payload without the caller's global data name.
        shipped <- unserialize(serialize(list(indices = indices, fun = fun), NULL))
        rm(list = intersect("candidate_loss_data", ls(envir = globalenv())),
           envir = globalenv())
        values <- lapply(shipped$indices, function(i) {
          .latents_collect_conditions(function() {
            candidate <- if (is.list(i)) i$candidate else i
            if (candidate == 2L) on_lost(i, loss) else shipped$fun(i)
          })
        })
        if (signal) lapply(values, `[[`, "value") else values
      })
    request <- if (family == "profiles") quote(enumerate_classes(
      candidate_loss_data, c("y", "z"), "person", n_profiles = 1L,
      n_group_classes = 1:2, model = "VVI", n_starts = 1L, max_iter = 30L,
      seed = 19L, cores = 2L)) else quote(enumerate_classes(
        candidate_loss_data, c("y", "z"), "person", family = "additive",
        n_group_classes = 1:2, between_variance = "varying", n_starts = 1L,
        max_iter = 30L, seed = 19L, cores = 2L))
    ordinary_error <- suppressWarnings(eval(request, envir = globalenv()))
    expect_identical(fitting, 1:2, info = family)
    fitting <- integer()
    lose_candidate <- TRUE
    worker_loss <- suppressWarnings(eval(request, envir = globalenv()))
    expect_identical(fitting, 1L, info = family)
    expect_identical(worker_loss$table, ordinary_error$table, info = family)
    expect_identical(worker_loss$fits, ordinary_error$fits, info = family)
    expect_null(worker_loss$fits[[2L]], info = family)
    expect_false(is.null(worker_loss$fits[[1L]]), info = family)
    expect_match(worker_loss$table$error[[2L]], "worker process ended", fixed = TRUE)
    expect_true(is.finite(worker_loss$table$log_likelihood[[1L]]), info = family)
  }))
})
