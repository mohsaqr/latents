# Released signatures are deliberately hard-coded, independently of HEAD.
test_that("0.9.21 formal positions are retained", {
  released <- list(
    multilpa = c("data", "vars", "id", "n_profiles", "n_group_classes",
      "profile_covariates", "group_covariates", "variance_model", "n_starts",
      "max_iter", "tol", "min_variance", "seed", "start", "missing",
      "profile_slopes", "covariance_model", "categorical", "min_probability",
      "time", "fixed", "centering", "volume", "shape", "orientation", "model",
      "select_start", "prior", "noise", "acceleration", "family",
      "between_variance", "weights", "ordinal", "count", "count_model",
      "count_dispersion"),
    plot.latents_family_enumeration = c("x", "criterion", "main", "subtitle", "..."),
    plot.latents_transition_enumeration = c("x", "criterion", "main", "subtitle", "..."),
    plot.latents_pooled = c("x", "level", "main", "..."),
    plot.multilpa = c("x", "what", "data", "scale", "category", "labels",
      "intervals", "cell_labels", "main", "subtitle", "..."),
    plot.multilpa_covariates = c("x", "what", "data", "scale", "category", "labels",
      "intervals", "cell_labels", "main", "subtitle", "..."),
    plot.multilpa_enumeration = c("x", "what", "criterion", "combine", "labels",
      "mark_minimum", "main", "subtitle", "..."),
    plot.multilpa_transitions = c("x", "what", "data", "scale", "category", "labels",
      "cell_labels", "main", "subtitle", "..."),
    # 0.9.21 shipped only an erroring stub, predict(object, ...), so no
    # released call can rely on what followed `object`.
    predict.latents_growth_mixture = "object"
  )
  invisible(lapply(names(released), function(name) {
    current <- names(formals(get(name, envir = asNamespace("latents"))))
    expect_identical(current[seq_along(released[[name]])], released[[name]], info = name)
  }))
  expect_identical(tail(names(formals(multilpa)), 1L), "boundary")
  expect_identical(eval(formals(plot.multilpa_enumeration)$what),
                   c("enumeration", "tree"))
  expect_identical(eval(formals(plot.multilpa_enumeration)$criterion),
                   c("aic", "bic_groups", "bic_individual", "icl_individual"))
  invisible(lapply(list(plot.latents_family_enumeration,
                        plot.latents_transition_enumeration), function(method) {
    expect_identical(eval(formals(method)$criterion), c("bic", "aic", "bic_individual"))
  }))
  expect_identical(eval(formals(plot.latents_pooled)$level),
                   c("profile", "group", "measurement"))
})

.api_0921_data <- function() {
  set.seed(921)
  data.frame(g = rep(seq_len(20L), each = 6L),
             a = stats::rnorm(120L, rep(c(-3, 3), each = 60L)),
             b = stats::rnorm(120L, rep(c(-2, 2), each = 60L)))
}

.api_0921_grid <- function() {
  enumerate_classes(.api_0921_data(), c("a", "b"), "g", n_profiles = 1:2,
                    n_group_classes = 1L, model = "VVI", n_starts = 1L, seed = 1L)
}

.api_0921_bootstrap <- function() {
  structure(list(replicates = data.frame(statistic = c(0.2, 0.5, 1, 2, 3)),
                 statistic = 2.5, n_valid = 5L, iter = 5L, p_value = 1 / 3),
            class = "multilpa_bootstrap_lrt")
}

.api_0921_pooled <- function() {
  estimates <- data.frame(level = "profile", outcome = "profile 2",
                          term = "(Intercept)", parameter = "logit", estimate = 0.5,
                          conf_low = 0.2, conf_high = 0.8)
  structure(list(estimates = estimates, imputations = estimates,
                 level = 0.95, m = 2L), class = "latents_pooled")
}

test_that("a prior at released multilpa position 28 is fitted", {
  data <- .api_0921_data()
  prior <- prior_control()
  # Every element is unnamed, including all intervening released defaults.
  arguments <- list(data, c("a", "b"), NULL, 2L, 1L, character(), character(),
                    "varying", 1L, 1000L, 1e-8, 1e-6, 1L, NULL, "error",
                    "shared", "diagonal", character(), 1e-10, NULL, character(),
                    "none", NULL, NULL, NULL, NULL, "likelihood", prior)
  positional <- do.call(multilpa, arguments)
  named <- multilpa(data, c("a", "b"), NULL, 2L, 1L, n_starts = 1L,
                    seed = 1L, prior = prior)
  expect_s3_class(positional, "multilpa")
  expect_identical(positional$prior, prior)
  expect_equal(positional$log_likelihood, named$log_likelihood)
  expect_equal(positional$means, named$means)
})

test_that("enumeration criteria remain available beside best", {
  grid <- .api_0921_grid()
  summary <- summary(grid)
  tables <- get_results(grid, "all")
  expect_s3_class(summary$criteria, "data.frame")
  expect_s3_class(tables$criteria, "data.frame")
  expect_identical(summary$criteria, summary$best)
  expect_identical(tables$criteria, tables$best)
  expect_identical(get_results(summary, "criteria"), tables$criteria)
  # A saved 0.9.21 summary contains criteria, with no best field anywhere.
  summary$best <- NULL
  summary$tables$best <- NULL
  expect_output(print(summary), "distinct candidate")
})

test_that("the explicit released enumeration default vector draws", {
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid <- .api_0921_grid()
  drawn <- plot(grid, what = c("enumeration", "tree"))
  expect_s3_class(drawn, "ggplot")
  expect_equal(ggplot2::ggplot_build(drawn)$data,
                ggplot2::ggplot_build(plot(grid))$data)
  expect_equal(plot(grid, c("enumeration", "tree"), "aic")$data,
                plot(grid, what = "aic")$data)
})

test_that("released plot arguments are accepted and the plots stay ggplot", {
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  bootstrap <- .api_0921_bootstrap()
  # The 0.9.21 base-graphics style constants have nothing to set on a ggplot;
  # they are accepted and reported, never an error.
  expect_warning(styled <- plot(bootstrap, style = list(line_width = 6)),
                 class = "latents_unused_argument")
  expect_s3_class(styled, "ggplot")
  expect_warning(plot(bootstrap, line_width = 4), class = "latents_unused_argument")
  returned <- withVisible(plot(bootstrap))
  expect_false(returned$visible)
  expect_s3_class(returned$value, "ggplot")
  pooled <- .api_0921_pooled()
  positional <- plot(pooled, "profile", "Released positional title")
  expect_s3_class(positional, "ggplot")
  expect_identical(positional$labels$title, "Released positional title")
  expect_equal(ggplot2::ggplot_build(plot(pooled, level = "profile"))$data,
               ggplot2::ggplot_build(plot(pooled, what = "profile"))$data)
})

test_that("family and transition grids take their released positional criterion", {
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  data <- .api_0921_data()
  data$time <- rep(seq_len(6L), 20L)
  family <- quietly(enumerate_classes(data, "a", "g", family = "additive",
                                      n_group_classes = 1:2, n_starts = 1L, seed = 1L))
  transition <- quietly(enumerate_classes(data, "a", "g", time = "time",
                                          n_profiles = 2:3, n_group_classes = 1L,
                                          n_starts = 1L, seed = 1L))
  expect_s3_class(family, "latents_family_enumeration")
  expect_s3_class(transition, "latents_transition_enumeration")
  invisible(lapply(list(family, transition), function(grid) {
    # Without a criterion the grid summary is drawn, every criterion a panel,
    # as plot_enumeration() draws it (0.9.21 drew BIC alone).
    expect_equal(plot(grid)$data, plot_enumeration(grid)$data)
    expect_equal(plot(grid, c("bic", "aic", "bic_individual"))$data,
                  plot(grid, what = "bic")$data)
    expect_equal(plot(grid, "aic", "Released title")$data,
                  plot(grid, what = "aic", main = "Released title")$data)
  }))
})
