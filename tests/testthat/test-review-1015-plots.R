# Regressions from the 2026-10-15 plot review. Each test reads the built plot
# (its data, layers or labels), not pixels.

review_lta_data <- function() {
  subset(course_engagement, student <= 40)
}

test_that("response whiskers are the inference table's own bounds", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
             n_classes = 2, n_starts = 1, seed = 1, tol = 1e-10)
  drawn <- expect_plot(plot(fit, what = "responses"))$data
  table <- parameter_inference(fit, boundary = "fix")
  table <- table[table$parameter == "response", , drop = FALSE]
  at <- match(sprintf("profile_%d\r%s:yes", drawn$profile,
                      c("time_with_friends", "on_social_media",
                        "sports")[drawn$position]),
              paste(table$outcome, table$term, sep = "\r"))
  # A probability held at its bound has no row, and so no whisker.
  reported <- !is.na(at)
  # Three of the six on x86 BLAS, four on arm64 macOS.
  expect_gte(sum(reported), 3L)
  expect_true(all(is.na(drawn$lower[!reported])))
  # A probability's interval is built on the logit scale, so it is not the
  # symmetric estimate +/- 1.96 SE the plot used to draw.
  expect_equal(drawn$lower, table$conf_low[at])
  expect_equal(drawn$upper, table$conf_high[at])
  symmetric <- drawn$mean[reported] -
    stats::qnorm(0.975) * table$standard_error[at[reported]]
  expect_gt(max(abs(drawn$lower[reported] - pmax(0, symmetric))), 1e-3)
  # The named function draws the same whiskers.
  expect_equal(plot_responses(fit)$data$lower, drawn$lower)
})

test_that("named plots refuse objects they cannot draw with a classed error", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  general <- quietly(lta(review_lta_data(), c("browse", "lectures"),
                         "student", 2, time = "sequence",
                         transitions = "occasion", n_starts = 1, seed = 1))
  expect_s3_class(general, "multilpa_lta")
  for (verb in list(plot_profiles, plot_responses, plot_prevalence,
                    plot_diagnostics)) {
    expect_error(verb(general), class = "latents_bad_argument")
  }
  expect_error(plot_profiles(general), "plot_transitions")
  # The general fit's transition matrix is still drawn by its own function.
  expect_plot(plot_transitions(general))
  grid <- enumerate_lpa(srl, c("intrinsic_value", "test_anxiety"), 1:2,
                        n_starts = 1, seed = 1)
  expect_error(plot_diagnostics(grid), class = "latents_bad_argument")
  expect_error(plot_profiles(data.frame(a = 1)), class = "latents_bad_argument")
})

test_that("plot_enumeration() and plot_tree() draw a transition grid", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grid <- enumerate_classes(review_lta_data(), c("browse", "lectures"),
                            "student", 2:3, 1, time = "sequence",
                            n_starts = 1, seed = 1)
  expect_s3_class(grid, "latents_transition_enumeration")
  # Refused before the grid contract; every grid kind now draws its views
  # (test-grid-contract.R).
  expect_plot(plot_enumeration(grid))
  expect_plot(plot_tree(grid))
  expect_plot(plot(grid, what = "bic"))
  expect_error(plot_enumeration(grid, what = "bic_groups"),
               class = "latents_bad_argument")
})

test_that("a noise fit's prevalence plot accounts for the noise share", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- lpa(srl, c("intrinsic_value", "self_efficacy", "test_anxiety"), 2,
             noise = TRUE, n_starts = 2, seed = 1)
  drawn <- expect_plot(plot_prevalence(fit))
  noise <- fit$noise_probability
  expect_equal(sum(drawn$data$mean) + noise, 1)
  expect_identical(drawn$labels$caption,
                   sprintf("The noise component takes the remaining %s.",
                           .gg_percent(noise)))
  plain <- lpa(srl, c("intrinsic_value", "test_anxiety"), 2, n_starts = 1,
               seed = 1)
  expect_null(plot_prevalence(plain)$labels$caption)
})

test_that("profile views say which indicators they leave out", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  esm <- transform(student_esm, happy = factor(happy, ordered = TRUE))
  ordinal <- multilpa(esm, c("happy", "exhausted", "relaxed"), "student", 2,
                      1, ordinal = "happy", n_starts = 1, seed = 1)
  caption <- "Not drawn: ordinal happy (get_results(fit, \"ordinal\"))."
  for (view in c("profiles", "bars", "heatmap")) {
    drawn <- expect_plot(plot(ordinal, what = view))
    expect_identical(drawn$labels$caption, caption)
  }
  expect_identical(levels(plot(ordinal)$data$indicator),
                   c("exhausted", "relaxed"))
  mixed <- multilpa(course_engagement,
                    c("browse", "lectures", "engagement"), "student", 2, 1,
                    categorical = "engagement", n_starts = 1, seed = 1)
  expect_identical(
    plot(mixed)$labels$caption,
    "Not drawn: categorical engagement (what = \"responses\").")
  plain <- lpa(srl, c("intrinsic_value", "test_anxiety"), 2, n_starts = 1,
               seed = 1)
  expect_null(plot(plain)$labels$caption)
})

test_that("covariate fits refuse an abbreviated prevalence with the reason", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- multilpa(course_engagement, c("browse", "lectures"), "student", 2, 2,
                  profile_covariates = "previous_grade", n_starts = 1,
                  seed = 1)
  expect_s3_class(fit, "multilpa_covariates")
  expect_error(plot(fit, what = "prob"), class = "latents_nothing_to_plot")
  expect_error(plot(fit, what = "probabilities"),
               class = "latents_nothing_to_plot")
  # "p" stays ambiguous, as match.arg() makes it, and is now a classed refusal.
  expect_error(plot(fit, what = "p"), class = "latents_bad_argument")
})

test_that("a centred fit's case-level views are refused as nothing to plot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fit <- multilpa(course_engagement, c("browse", "lectures"), "student", 2, 2,
                  centering = "person", n_starts = 1, seed = 1)
  expect_error(plot(fit, what = "raincloud"), class = "latents_nothing_to_plot")
  # The former class is kept for code that caught it.
  expect_error(plot(fit, what = "parallel"), class = "latents_bad_argument")
})
