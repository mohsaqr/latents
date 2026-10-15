# End-to-end regressions for the 2026-10-13 plotting/prediction audit.
audit_plot_data <- function() {
  set.seed(1310)
  data <- data.frame(g = rep(seq_len(20L), each = 10L),
                     x = stats::rnorm(200), z = stats::rnorm(200))
  data$x <- data$x + rep(c(-2, 2), each = 100)
  data$w <- rep(c(1, 10), each = 100)
  data
}

test_that("named diagnostics apply the generic noise restrictions", {
  skip_if_not_installed("ggplot2")
  fit <- lpa(audit_plot_data(), c("x", "z"), 2, noise = TRUE,
              n_starts = 1, seed = 1)
  invisible(lapply(c("entropy", "avepp", "posteriors"), function(view) {
    expect_error(plot_diagnostics(fit, what = view),
                 class = "latents_unsupported_noise")
    expect_error(plot(fit, what = view), class = "latents_unsupported_noise")
  }))
  expect_plot(plot_diagnostics(fit, what = "sizes"))
  expect_named(plot_diagnostics(fit), "sizes")
})

test_that("weighted tree shares agree with fitted effective proportions", {
  skip_if_not_installed("ggplot2")
  data <- audit_plot_data()
  grids <- list(
    rows = enumerate_lpa(data, c("x", "z"), 1:2, model = "EEI", weights = "w",
                          n_starts = 1, seed = 1),
    groups = enumerate_classes(data, c("x", "z"), "g", n_profiles = 1:2,
                                n_group_classes = 1, model = "EEI", weights = "w",
                                n_starts = 1, seed = 1)
  )
  invisible(lapply(grids, function(models) {
    tree <- expect_plot(plot_tree(models))$layers[[2L]]$data
    invisible(lapply(1:2, function(k) {
      counts <- get_results(get_candidate(models, n_profiles = k), "counts")
      expected <- counts$effective_proportion[counts$level == "individuals"]
      drawn <- tree$share[tree$n_profiles == k]
      expect_equal(sort(drawn), sort(expected), tolerance = 1e-10)
    }))
  }))
})

test_that("released view abbreviations and named subclass dispatch still work", {
  skip_if_not_installed("ggplot2")
  data <- audit_plot_data()
  models <- enumerate_lpa(data, c("x", "z"), 1:2, model = "EEI",
                          n_starts = 1, seed = 1)
  invisible(lapply(c("e", "enum", "enumeration"), function(view) {
    expect_equal(plot(models, what = view)$data, plot(models)$data)
  }))
  invisible(lapply(c("t", "tr", "tree"), function(view) {
    expect_equal(plot(models, what = view)$layers[[2L]]$data,
                 plot_tree(models)$layers[[2L]]$data)
  }))
  expect_equal(plot(models, what = NULL)$data, plot(models)$data)
  fit <- get_candidate(models, n_profiles = 2)
  class(fit) <- c("downstream_profile", class(fit))
  expect_equal(plot_profiles(fit, intervals = FALSE)$data,
               plot(fit, intervals = FALSE)$data)
  expect_plots(plot(fit, what = "all", intervals = FALSE))
})

test_that("bar plots can skip inference and error bars", {
  skip_if_not_installed("ggplot2")
  fit <- lpa(audit_plot_data(), c("x", "z"), 2, n_starts = 1, seed = 1)
  local_mocked_bindings(.multilpa_mean_error_matrix = function(...) {
    stop("inference should not be evaluated")
  })
  bars <- expect_plot(plot(fit, what = "bars", intervals = FALSE))
  expect_false(any(vapply(bars$layers, function(layer) {
    inherits(layer$geom, "GeomErrorbar")
  }, logical(1))))
})

audit_binomial_fit <- function() {
  set.seed(2301)
  data <- data.frame(x = stats::rnorm(160), n = rep(c(10L, 20L), 80))
  data$success <- stats::rbinom(160, data$n, stats::plogis(.5 + data$x))
  data$failure <- data$n - data$success
  list(data = data, fit = mixture_regression(cbind(success, failure) ~ x,
    data, 1, family = "binomial", n_starts = 1, seed = 1))
}

test_that("binomial count predictions preserve exposure and success-count units", {
  case <- audit_binomial_fit()
  data <- case$data
  fit <- case$fit
  expected <- predict(fit)
  expect_true("trials" %in% names(formals(
    utils::getS3method("predict", "latents_mixture_regression"))))
  expect_equal(predict(fit, data), expected)
  expect_equal(predict(fit, data, type = "class_response"),
               predict(fit, type = "class_response"))
  predictors <- data[c("x", "n")]
  expect_equal(predict(fit, predictors, trials = "n"), expected)
  expect_equal(predict(fit, predictors, trials = data$n), expected)
  per_trial <- predict(fit, predictors, trials = 1)$fitted
  expect_equal(expected$fitted, per_trial * data$n)
  expect_error(predict(fit, predictors), class = "latents_bad_data")
  expect_error(predict(fit, predictors, type = "class_response"),
               class = "latents_bad_data")
  expect_equal(predict(fit, predictors, type = "class_response", trials = 1)$fitted,
               per_trial)
  expect_equal(predict(fit, predictors, trials = 10)$fitted, 10 * per_trial)
  expect_error(predict(fit, predictors, trials = -.1), class = "latents_bad_data")
  expect_error(predict(fit, predictors, trials = c(1, 2)), class = "latents_bad_data")
  expect_error(predict(fit, predictors, trials = NA_real_), class = "latents_bad_data")
  expect_equal(predict(fit, data, type = "posterior", trials = "n"),
               predict(fit, data, type = "posterior"))
  expect_error(predict(fit, data, type = "posterior", trials = 1),
               class = "latents_bad_data")
  expect_equal(nrow(predict(fit, data[FALSE, ])), 0L)
  expect_equal(predict(fit, predictors, trials = 0)$fitted, rep(0, nrow(data)))
})

test_that("Bernoulli predictions do not require or validate unused outcomes", {
  set.seed(2302)
  data <- data.frame(x = stats::rnorm(100))
  data$y <- factor(ifelse(stats::rbinom(100, 1, stats::plogis(data$x)), "yes", "no"))
  fit <- mixture_regression(y ~ x, data, 1, family = "binomial", n_starts = 1, seed = 1)
  new <- data.frame(x = 0, y = factor("unused"))
  expect_equal(predict(fit, new), predict(fit, new["x"]))
  expect_equal(predict(fit, new, trials = 20)$fitted,
               20 * predict(fit, new)$fitted)
})

test_that("mixture prediction errors match the documented data condition", {
  case <- audit_binomial_fit()
  expect_error(predict(case$fit, data.frame(x = NA_real_), trials = 10),
               class = "latents_bad_data")
  expect_error(predict(case$fit, data.frame(z = 1), trials = 10),
               class = "latents_bad_data")
  expect_error(predict(case$fit, data.frame(x = Inf), trials = 10),
               class = "latents_bad_data")
  expect_error(predict(case$fit, data.frame(x = 1), trials = "absent"),
               class = "latents_bad_data")
})

test_that("ported regression plots draw once and remain reusable ggplots", {
  skip_if_not_installed("ggplot2")
  case <- audit_binomial_fit()
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  printed <- 0L
  plot_class <- class(ggplot2::ggplot())[1L]
  original <- utils::getS3method("print", plot_class)
  # Counting the generic dispatch while still printing verifies actual drawing.
  method <- function(x, ...) {
    printed <<- printed + 1L
    original(x, ...)
  }
  base::registerS3method("print", plot_class, method,
                         envir = asNamespace("ggplot2"))
  on.exit(base::registerS3method("print", plot_class, original,
                                 envir = asNamespace("ggplot2")), add = TRUE)
  grid::grid.newpage()
  invisible(lapply(seq_len(2L), function(i) {
    plot(case$fit)
    invisible(NULL)
  }))
  expect_equal(printed, 2L)
  expect_gt(length(grid::grid.ls(print = FALSE)$name), 0L)
  value <- withVisible(plot(case$fit))
  expect_s3_class(value$value, "ggplot")
  expect_false(value$visible)
  expect_equal(printed, 3L)
  expect_plot(value$value + ggplot2::theme(legend.position = "bottom"))
})
