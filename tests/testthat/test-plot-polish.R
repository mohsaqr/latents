# The plot and print behaviour added for the comparison articles: intervals
# and observations on the profile plot, reasons instead of errors when there
# are no intervals, the categorical heatmap, the covariate views, one BIC panel
# for a single-level grid, the compact enumeration print, the LCA header and
# the single-class intercept label.

.polish_draw <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

.polish_fit <- function(...) {
  quietly(multilpa(subset(course_engagement, student <= 40),
                   c("browse", "lectures", "forum_read"), "student",
                   n_profiles = 2, n_starts = 2, seed = 1, ...))
}

test_that("the profile plot draws intervals, and the raincloud the data", {
  fit <- .polish_fit(n_group_classes = 1)
  errors <- .multilpa_mean_error_matrix(fit, NULL)
  expect_true(is.matrix(errors))
  expect_identical(dim(errors), dim(fit$means))
  expect_true(all(errors > 0))
  observed <- .multilpa_plot_observations(fit, .multilpa_continuous_names(fit))
  expect_identical(nrow(observed$values), length(observed$profile))
  expect_true(all(observed$profile %in% seq_len(fit$n_profiles)))
  skip_if_not_installed("ggplot2")
  with_whiskers <- expect_plot(plot(fit))
  # The whiskers are the fit's own 95% Wald intervals.
  drawn <- with_whiskers$data
  expect_equal(drawn$upper - drawn$mean,
               stats::qnorm(0.975) * errors[cbind(drawn$profile, drawn$position)])
  without <- expect_plot(plot(fit, intervals = FALSE))
  expect_true(all(is.na(without$data$upper)))
  expect_plot(plot(fit, what = "raincloud"))
  expect_plot(plot(fit, what = "raincloud", scale = "standardized"))
  # Drawing does not touch the random-number stream.
  set.seed(3)
  before <- stats::runif(1)
  set.seed(3)
  ggplot2::ggplot_build(plot(fit, what = "raincloud"))
  expect_identical(stats::runif(1), before)
  categorical <- quietly(multilca(student_esm, c("sports", "walking"), id = NULL,
                                  n_profiles = 2, n_starts = 1, seed = 1))
  expect_error(plot(categorical, what = "raincloud"),
               class = "latents_no_continuous")
})

test_that("a fit without intervals says why instead of failing", {
  fit <- quietly(multilpa(subset(course_engagement, student <= 40),
                          c("browse", "lectures", "forum_read"), "student",
                          n_profiles = 2, n_group_classes = 1, n_starts = 1,
                          seed = 1, max_iter = 2))
  expect_false(fit$converged)
  reason <- .multilpa_mean_error_matrix(fit, NULL)
  expect_identical(reason, "no intervals: the fit did not converge")
  skip_if_not_installed("ggplot2")
  # The plot is drawn without whiskers, and its subtitle says why.
  bars <- expect_plot(plot(fit, what = "bars"))
  expect_match(bars$labels$subtitle, reason, fixed = TRUE)
  profiles <- expect_plot(plot(fit))
  expect_match(profiles$labels$subtitle, reason, fixed = TRUE)
})

test_that("an all-categorical fit gets a response heatmap", {
  binary <- quietly(multilca(student_esm,
                             c("time_with_friends", "on_social_media", "sports"),
                             id = NULL, n_profiles = 2, n_starts = 2, seed = 1))
  skip_if_not_installed("ggplot2")
  responses <- expect_plot(plot(binary, what = "heatmap"))
  expect_true(all(responses$data$probability >= 0 &
                    responses$data$probability <= 1))
  mixed <- .polish_fit(n_group_classes = 1)
  expect_plot(plot(mixed, what = "heatmap"))
})

test_that("covariate fits draw the measurement and classification views", {
  skip_if_not_installed("ggplot2")
  fit <- .polish_fit(n_group_classes = 1, profile_covariates = "previous_grade")
  invisible(lapply(
    c("profiles", "bars", "heatmap", "raincloud", "parallel", "pairs", "sizes",
      "avepp", "entropy", "posteriors"),
    function(view) expect_plot(plot(fit, what = view))))
})

test_that("a single-level grid draws one BIC panel, named for display", {
  grid <- data.frame(bic_groups = c(1, 2), bic_individual = c(1, 2),
                     icl_groups = c(1, 2), icl_individual = c(3, 4), aic = c(5, 6))
  chosen <- .multilpa_distinct_criteria(
    grid, c("aic", "bic_groups", "bic_individual", "icl_individual"))
  expect_identical(unname(chosen), c("aic", "bic_groups", "icl_individual"))
  expect_identical(names(chosen), c("AIC", "BIC", "ICL (individuals)"))
})

test_that("the enumeration print shows a compact candidates table", {
  grid <- quietly(enumerate_classes(iris, c("Sepal.Length", "Petal.Length"),
                                    id = NULL, n_profiles = 1:2,
                                    model = c("EEI", "VVI"), n_starts = 1,
                                    seed = 1))
  compact <- .multilpa_compact_candidates(get_results(grid, "candidates"))
  expect_true(all(c("n_profiles", "model", "log_likelihood", "bic") %in%
                    names(compact)))
  expect_false(any(c("bic_groups", "bic_individual", "n_group_classes") %in%
                     names(compact)))
  expect_output(print(summary(grid)), "returns all of them")
})

test_that("a latent class model is printed as one", {
  fit <- quietly(multilca(student_esm,
                          c("time_with_friends", "on_social_media", "sports"),
                          id = NULL, n_profiles = 2, n_starts = 2, seed = 1))
  output <- utils::capture.output(print(fit))
  expect_match(output[1L], "^Latent class analysis: 2 classes")
  expect_false(any(grepl("residual covariance", output)))
  expect_true(any(grepl("3 categorical indicators", output)))
})

test_that("with one group class the membership intercept is named as one", {
  fit <- .polish_fit(n_group_classes = 1, profile_covariates = "previous_grade")
  terms <- subset(get_results(fit, "coefficients"), level == "profile")$term
  expect_identical(terms, c("(Intercept)", "previous_grade"))
})

test_that("lpa() and lca() are the single-level fits, without the notice", {
  vars <- c("cognitive_strategies", "intrinsic_value", "self_efficacy")
  expect_no_message(fit <- lpa(srl, vars, n_profiles = 2, n_starts = 2, seed = 1))
  same <- suppressMessages(multilpa(srl, vars, id = NULL, n_profiles = 2,
                                    n_starts = 2, seed = 1))
  expect_equal(fit$log_likelihood, same$log_likelihood)
  expect_true(isTRUE(fit$single_level))
  activities <- c("time_with_friends", "on_social_media", "sports")
  expect_no_message(classes <- lca(student_esm, activities, n_classes = 2,
                                   n_starts = 2, seed = 1))
  expect_identical(classes$n_profiles, 2L)
  expect_equal(classes$log_likelihood,
               suppressMessages(multilca(student_esm, activities, id = NULL,
                                         n_profiles = 2, n_starts = 2,
                                         seed = 1))$log_likelihood)
  expect_error(lpa(srl, vars, n_profiles = 2, id = "x"),
               class = "latents_bad_argument")
  expect_error(lca(student_esm, activities, n_classes = 2, n_group_classes = 2),
               class = "latents_bad_argument")
})

test_that("`model` names the covariance structure as enumerate_classes() does", {
  vars <- c("cognitive_strategies", "intrinsic_value", "test_anxiety")
  named <- lpa(srl, vars, n_profiles = 2, model = "EEE", n_starts = 2, seed = 1)
  switched <- lpa(srl, vars, n_profiles = 2, variance_model = "equal",
                  covariance_model = "full", n_starts = 2, seed = 1)
  expect_identical(named$covariance_structure, "EEE")
  expect_equal(named$log_likelihood, switched$log_likelihood)
  expect_identical(lpa(srl, vars, n_profiles = 2, model = "VEV", n_starts = 1,
                       seed = 1)$covariance_structure, "VEV")
  expect_error(lpa(srl, vars, n_profiles = 2, model = "EEE",
                   variance_model = "equal"), class = "latents_bad_argument")
  expect_error(lpa(srl, vars, n_profiles = 2, model = "XYZ"),
               class = "latents_bad_argument")
})

test_that("the measurement tables carry standard errors without being given data", {
  vars <- c("cognitive_strategies", "intrinsic_value", "test_anxiety")
  invisible(lapply(c("EEI", "VVI", "EEE", "VVV"), function(structure) {
    fit <- lpa(srl, vars, n_profiles = 2, model = structure, n_starts = 2, seed = 1)
    profiles <- get_results(fit, "profiles")
    expect_false(anyNA(profiles$mean_standard_error), label = structure)
    expect_false(anyNA(profiles$variance_standard_error), label = structure)
    # The same errors parameter_inference() reports for the means.
    inference <- parameter_inference(fit)
    means <- inference[inference$parameter == "mean", ]
    expect_equal(sort(profiles$mean_standard_error),
                 sort(means$standard_error), label = structure)
  }))
  # A fit without errors still returns its table, with the reason given.
  unconverged <- quietly(lpa(srl, vars, n_profiles = 2, n_starts = 1, seed = 1,
                             max_iter = 2))
  expect_message(table <- get_results(unconverged, "profiles"),
                 class = "latents_no_standard_errors")
  expect_true(all(is.na(table$mean_standard_error)))
  # Printing shows the estimates without computing errors.
  expect_no_message(utils::capture.output(print(unconverged)))
})

test_that("diagnostics() draws its plots by default", {
  fit <- lpa(srl, c("self_efficacy", "test_anxiety"), n_profiles = 2,
             n_starts = 2, seed = 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  before <- grDevices::recordPlot()
  result <- diagnostics(fit)
  expect_s3_class(result, "multilpa_diagnostics")
  expect_false(identical(grDevices::recordPlot(), before))
  expect_s3_class(diagnostics(fit, plots = FALSE), "multilpa_diagnostics")
})
