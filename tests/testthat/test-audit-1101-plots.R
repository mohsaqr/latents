skip_on_cran()
skip_if_not_installed("ggplot2")

# Audit 2026-11-01, plots slice. Each block guards one fix in R/plot.R,
# R/plot-gg.R or R/plot-named.R.

audit_plot_data <- function() {
  data <- course_engagement
  # Row weights for a single-level fit: some zero, some three.
  data$w <- 1
  data$w[seq(1L, nrow(data), by = 7L)] <- 0
  data$w[seq(2L, nrow(data), by = 5L)] <- 3
  # Group weights for a two-level fit: a third of the students weigh zero.
  data$gw <- c(0, 1, 3)[(data$student %% 3L) + 1L]
  data
}

audit_vars <- c("browse", "lectures", "forum_read")

audit_weighted_lpa <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- lpa(audit_plot_data(), audit_vars, n_profiles = 3, n_starts = 2,
                  seed = 1, weights = "w")
    }
    fit
  }
})

audit_weighted_multilpa <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- multilpa(audit_plot_data(), audit_vars, "student", n_profiles = 3,
                       n_group_classes = 2, n_starts = 2, seed = 1,
                       weights = "gw")
    }
    fit
  }
})

# Weighted mean of each unit's largest posterior within its modal class, and
# weighted relative entropy, from the posteriors and weights alone.
audit_weighted_classification <- function(posteriors, weights) {
  assigned <- max.col(posteriors, ties.method = "first")
  certainty <- posteriors[cbind(seq_len(nrow(posteriors)), assigned)]
  average <- vapply(seq_len(ncol(posteriors)), \(k) {
    rows <- assigned == k
    sum(weights[rows] * certainty[rows]) / sum(weights[rows])
  }, numeric(1))
  cells <- ifelse(posteriors > 0, posteriors * log(posteriors), 0)
  entropy <- 1 - sum(weights * -rowSums(cells)) /
    (sum(weights) * log(ncol(posteriors)))
  list(average = average, entropy = entropy)
}

strip_means <- function(plot) {
  means <- plot$layers[[5L]]$data
  means$value[order(means$profile)]
}

test_that("weighted strip means and entropy use the sampling weights (individuals)", {
  fit <- audit_weighted_lpa()
  data <- audit_plot_data()
  assignments <- get_results(fit, "assignments")
  posterior_columns <- grep("^posterior_profile_", names(assignments), value = TRUE)
  carrying <- data$w > 0
  expected <- audit_weighted_classification(
    as.matrix(assignments[carrying, posterior_columns]),
    data$w[carrying] / mean(data$w))

  posteriors_plot <- plot_diagnostics(fit, what = "posteriors")
  expect_equal(strip_means(posteriors_plot), expected$average, tolerance = 1e-10)
  table <- get_results(fit, "classification")
  expect_equal(strip_means(posteriors_plot),
               table$average_posterior[table$level == "individuals"],
               tolerance = 1e-12)
  # The label on each mean marker shows the same number.
  expect_setequal(posteriors_plot$layers[[5L]]$data$label,
                  sprintf("mean %.2f", expected$average))

  entropy_plot <- plot_diagnostics(fit, what = "entropy")
  expect_match(entropy_plot$labels$subtitle,
               sprintf("1 - mean = %.3f", expected$entropy), fixed = TRUE)
  entropy_table <- get_results(fit, "entropy")
  expect_match(entropy_plot$labels$subtitle, sprintf(
    "1 - mean = %.3f",
    entropy_table$relative_entropy[entropy_table$level == "individuals"]),
    fixed = TRUE)
  # Weighting matters here: the unweighted values differ in the drawn digits.
  unweighted <- audit_weighted_classification(
    as.matrix(assignments[carrying, posterior_columns]), rep(1, sum(carrying)))
  expect_false(identical(sprintf("%.3f", unweighted$entropy),
                         sprintf("%.3f", expected$entropy)))
})

test_that("weighted strip means and entropy use the sampling weights (groups)", {
  fit <- audit_weighted_multilpa()
  table <- get_results(fit, "classification")
  for (level in c("individuals", "groups")) {
    prefix <- if (identical(level, "groups")) "group_" else ""
    plot <- plot_diagnostics(fit, what = paste0(prefix, "posteriors"))
    expect_equal(strip_means(plot),
                 table$average_posterior[table$level == level],
                 tolerance = 1e-12)
    entropy <- get_results(fit, "entropy")
    expect_match(
      plot_diagnostics(fit, what = paste0(prefix, "entropy"))$labels$subtitle,
      sprintf("1 - mean = %.3f", entropy$relative_entropy[entropy$level == level]),
      fixed = TRUE)
  }
})

test_that("unweighted strip means are the plain means, as before", {
  fit <- lpa(course_engagement, audit_vars, n_profiles = 2, n_starts = 2,
             seed = 1)
  posteriors <- fit$subject_posteriors
  assigned <- max.col(posteriors, ties.method = "first")
  certainty <- posteriors[cbind(seq_len(nrow(posteriors)), assigned)]
  expect_identical(strip_means(plot_diagnostics(fit, what = "posteriors")),
                   as.vector(tapply(certainty, assigned, mean)))
})

test_that("the named plot verbs and plot() refuse bad arguments by class", {
  fit <- lpa(course_engagement, audit_vars, n_profiles = 2, n_starts = 2,
             seed = 1)
  expect_error(plot_profiles(fit, what = "lyne"), class = "latents_bad_argument")
  expect_error(plot_profiles(fit, statistic = "average"),
               class = "latents_bad_argument")
  expect_error(plot_profiles(fit, scale = "z"), class = "latents_bad_argument")
  expect_error(plot_profiles(fit, intervals = NA), class = "latents_bad_argument")
  expect_error(plot_diagnostics(fit, what = "spread"),
               class = "latents_bad_argument")
  expect_error(plot(fit, labels = NA), class = "latents_bad_argument")
  expect_error(plot(fit, what = "responses", category = c(1, 2)),
               class = "latents_bad_argument")
  # Abbreviations keep working, as match.arg() took them.
  expect_identical(
    ggplot2::layer_data(plot_profiles(fit, what = "bar", intervals = FALSE), 2L),
    ggplot2::layer_data(plot(fit, what = "bars", intervals = FALSE), 2L))
  expect_s3_class(plot_diagnostics(fit, what = "avep"), "ggplot")

  grid <- enumerate_lpa(course_engagement, audit_vars, n_profiles = 1:3,
                        n_starts = 2, seed = 1)
  expect_error(plot_selection(grid, what = "delta_bic"),
               class = "latents_bad_argument")
  expect_error(plot_selection(grid, cutoff = 2), class = "latents_bad_argument")
  expect_error(plot_enumeration(grid, what = "avepp", cutoff = 0),
               class = "latents_bad_argument")
  expect_error(plot_enumeration(grid, combine = NA),
               class = "latents_bad_argument")
  expect_error(plot_enumeration(grid, labels = "yes"),
               class = "latents_bad_argument")
})

test_that("a fractional category index is refused, a whole one is used", {
  data <- course_engagement
  data$post <- factor(ifelse(data$forum_post > 0, "yes", "no"))
  data$present <- factor(ifelse(data$attendance > 0, "hi", "lo"))
  data$grade <- cut(data$previous_grade, c(-Inf, -0.5, 0.5, Inf),
                    labels = c("low", "mid", "high"))
  fit <- lca(data, c("post", "present", "grade"), n_classes = 2, n_starts = 2,
             seed = 1)
  expect_error(plot_responses(fit, category = 1.6, intervals = FALSE),
               class = "latents_unknown_category")
  second <- plot_responses(fit, category = 2, intervals = FALSE)$data
  expected <- vapply(seq_len(nrow(second)), \(row) {
    fit$response_probabilities[[second$position[row]]][second$profile[row], 2L]
  }, numeric(1))
  expect_equal(second$mean, expected)
})

test_that("pairs marks each profile by shape and each ellipse by its mean", {
  fit <- lpa(course_engagement, audit_vars, n_profiles = 3, n_starts = 2,
             seed = 1)
  plot <- expect_plot(plot_profiles(fit, what = "pairs"))
  key <- latents:::.gg_profile_key(fit)
  # The cases: one shape per profile, the shape the other views give it.
  cases <- plot$layers[[4L]]$data
  built <- ggplot2::layer_data(plot, 4L)
  shapes <- tapply(built$shape, as.character(cases$profile_label),
                   \(shape) unique(shape))
  expect_identical(as.vector(shapes[as.character(key$label)]), key$shape)
  # The centre markers sit at the fitted means.
  centres <- plot$layers[[6L]]$data
  profile <- match(as.character(centres$profile_label), as.character(key$label))
  expect_equal(centres$x,
               fit$means[cbind(profile, match(as.character(centres$xvar),
                                              colnames(fit$means)))])
  expect_equal(centres$y,
               fit$means[cbind(profile, match(as.character(centres$yvar),
                                              colnames(fit$means)))])
  expect_equal(nrow(centres), 3 * choose(length(audit_vars), 2L))
})
