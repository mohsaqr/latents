skip_on_cran()

astra_surface_data <- function() {
  withr::local_seed(1027)
  data.frame(id = rep(seq_len(12L), each = 5L), x = stats::rnorm(60),
             y = stats::rnorm(60), a = sample(0:1, 60, TRUE), b = sample(0:1, 60, TRUE))
}

test_that("call bindings cannot overwrite arguments or the forwarding function", {
  data <- astra_surface_data()
  reference <- lpa(data, c("x", "y"), n_profiles = 1, n_starts = 1, seed = 1)
  caller <- list2env(list(.vars = data, .vars.1 = 1L, .multilpa_setup = data),
                     parent = environment())
  collision <- eval(quote(multilpa(.vars, vars = c("x", "y"), id = NULL,
                                   n_profiles = .vars.1, n_starts = 1, seed = 1)), caller)
  reserved <- eval(quote(multilpa(.multilpa_setup, vars = c("x", "y"), id = NULL,
                                  n_profiles = 1, n_starts = 1, seed = 1)), caller)
  expect_equal(collision$log_likelihood, reference$log_likelihood)
  expect_equal(reserved$log_likelihood, reference$log_likelihood)
  categorical <- eval(quote(multilca(.vars, vars = c("a", "b"), id = NULL,
                                     n_profiles = 1, n_starts = 1, seed = 1)), caller)
  expect_equal(categorical$log_likelihood,
                lca(data, c("a", "b"), n_classes = 1, n_starts = 1, seed = 1)$log_likelihood)
})

test_that("multilca forces data, variables, controls and dots exactly once", {
  data <- astra_surface_data()
  counts <- new.env(parent = emptyenv())
  counted <- function(name, value) {
    counts[[name]] <- (counts[[name]] %||% 0L) + 1L
    value
  }
  fitted <- multilca(counted("data", data), vars = counted("vars", c("a", "b")),
                      id = counted("id", NULL), n_profiles = counted("profiles", 1L),
                      n_group_classes = counted("groups", 1L),
                      n_starts = counted("starts", 1L), seed = counted("seed", 1L))
  expect_true(all(unlist(as.list(counts)) == 1L))
  expect_setequal(ls(counts), c("data", "vars", "id", "profiles", "groups", "starts", "seed"))
  expect_equal(fitted$log_likelihood,
                lca(data, c("a", "b"), n_classes = 1, n_starts = 1, seed = 1)$log_likelihood)
  expect_identical(fitted$call[[1L]], quote(multilca))
})

test_that("bounded forwarding retains omission and NULL controls", {
  data <- astra_surface_data()
  fitted <- multilca(data, c("a", "b"), id = NULL, n_profiles = 1, seed = NULL,
                      n_starts = 1)
  expect_identical(fitted$n_group_classes, 1L)
  bad <- tryCatch(multilca(data, "absent", id = NULL, n_profiles = 1), error = identity)
  expect_s3_class(bad, "latents_bad_data")
  expect_null(conditionCall(bad))
})

astra_sensitivity_table <- function(log_likelihood, converged) {
  n <- length(log_likelihood)
  result <- data.frame(seed = seq_len(n), log_likelihood = log_likelihood,
                        converged = converged, iterations = rep(0L, n),
                        optimum = rep(1L, n), best = rep(TRUE, n), agreement = rep(1, n))
  attr(result, "tolerance") <- 1e-3
  class(result) <- c("latents_sensitivity", "data.frame")
  result
}

test_that("sensitivity plots never call unfinished refits maxima", {
  skip_if_not_installed("ggplot2")
  unfinished <- plot(astra_sensitivity_table(c(-10, -10, -11), rep(FALSE, 3)))
  expect_match(unfinished$labels$title, "No seed converged")
  expect_identical(unique(as.character(unfinished$data$status)), "unfinished")
  expect_match(unfinished$labels$caption, "unfinished for seeds 1, 2, 3")
  # An unconverged run may exceed the converged maximum: it must not steal
  # the converged winner's status or enter the replicated-maximum count.
  mixed <- plot(astra_sensitivity_table(c(-8, -10, -10.0001, -12, NA),
                                        c(FALSE, TRUE, TRUE, TRUE, NA)))
  expect_match(mixed$labels$title, "2 of 5 seeds reached the best converged maximum")
  first_panel <- mixed$data[mixed$data$panel == levels(mixed$data$panel)[[1L]], ]
  expect_identical(as.character(first_panel$status),
                   c("unfinished", "best maximum", "best maximum", "a lower maximum"))
  expect_equal(first_panel$value, c(-2, 0, 0.0001, 2), tolerance = 1e-10)
  expect_match(mixed$labels$caption, "failed .* seed 5")
  built <- ggplot2::ggplot_build(mixed)
  expect_identical(built$layout$panel_params[[1L]]$x$get_limits(), as.character(1:5))
})

test_that("zero-iteration sensitivity refits remain unfinished in the plot", {
  skip_if_not_installed("ggplot2")
  fit <- suppressWarnings(lpa(astra_surface_data(), c("x", "y"), n_profiles = 2,
                               n_starts = 1, seed = 1))
  seeds <- suppressWarnings(sensitivity(fit, seeds = 1:3, max_iter = 0, n_starts = 1))
  expect_true(all(!seeds$converged))
  expect_true(all(seeds$iterations == 0L))
  shown <- plot(seeds)
  expect_match(shown$labels$title, "No seed converged")
  expect_true(all(shown$data$status == "unfinished"))
})
