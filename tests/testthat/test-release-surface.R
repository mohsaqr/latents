# Release review 2026-10-27, cluster "surface": the user-facing contract of
# plot(), summary(), report(), descriptives(), get_candidate(),
# starting_values(), diagnostics() and bootstrap_lrt(). Each block names the
# finding it guards; each fails on the code before the fix.

surface_lpa <- function() {
  lpa(srl, names(srl), n_profiles = 2, n_starts = 2, seed = 1)
}

# API-3 ---------------------------------------------------------------------

test_that("plot() of an LCA on labelled categories draws its response view by default", {
  skip_if_not_installed("ggplot2")
  fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
             n_classes = 2, n_starts = 2, seed = 1)
  drawn <- plot(fit)
  expect_s3_class(drawn, "ggplot")
  expect_identical(.multilpa_default_view(fit), "responses")
  expect_identical(
    ggplot2::ggplot_build(drawn)$plot$labels$title,
    ggplot2::ggplot_build(plot(fit, what = "responses"))$plot$labels$title)
  # Asked for by name, the profile summaries still refuse, as documented.
  expect_error(plot(fit, what = "profiles"), class = "latents_bad_argument")
  expect_error(plot_profiles(fit), class = "latents_bad_argument")
})

test_that("a fit with numeric categories or a continuous indicator keeps the profile default", {
  skip_if_not_installed("ggplot2")
  numeric_classes <- lca(student_esm, c("happy", "relaxed"), n_classes = 2,
                         n_starts = 2, seed = 1)
  expect_identical(.multilpa_default_view(numeric_classes), "profiles")
  # report() and plot(what = "all") still draw its profile summaries.
  expect_true("profiles" %in% .multilpa_supported_views(numeric_classes))
  # Two classes of two binary items are not identified (CORE-5).
  expect_warning(labelled <- lca(student_esm, c("time_with_friends", "sports"),
                                 n_classes = 2, n_starts = 1, seed = 1),
                 class = "latents_unidentified")
  expect_false("profiles" %in% .multilpa_supported_views(labelled))
  expect_identical(.multilpa_default_view(surface_lpa()), "profiles")
})

# API-5 ---------------------------------------------------------------------

test_that("an unknown or malformed plot view is a classed refusal", {
  skip_if_not_installed("ggplot2")
  fit <- surface_lpa()
  expect_error(plot(fit, what = "nope"), class = "latents_bad_argument")
  expect_error(plot(fit, what = NA_character_), class = "latents_bad_argument")
  expect_error(plot(fit, what = c("profiles", "bars")), class = "latents_bad_argument")
  expect_error(plot(fit, scale = "log"), class = "latents_bad_argument")
  expect_error(plot(fit, what = "profiles", statistic = "max"),
               class = "latents_bad_argument")
  # The released abbreviations still name their view.
  expect_s3_class(plot(fit, what = "bar"), "ggplot")
})

test_that("the classed matcher keeps match.arg()'s contract", {
  choices <- c("profiles", "probabilities", "bars")
  expect_identical(.latents_match_arg(NULL, choices), "profiles")
  expect_identical(.latents_match_arg(choices, choices), "profiles")
  expect_identical(.latents_match_arg("bars", choices), "bars")
  expect_identical(.latents_match_arg("ba", choices), "bars")
  # Ambiguous, unknown, missing, several or non-character: refused by class.
  expect_error(.latents_match_arg("pro", choices), class = "latents_bad_argument")
  expect_error(.latents_match_arg("nope", choices), class = "latents_bad_argument")
  expect_error(.latents_match_arg(NA_character_, choices), class = "latents_bad_argument")
  expect_error(.latents_match_arg(c("bars", "profiles"), choices),
               class = "latents_bad_argument")
  expect_error(.latents_match_arg(1, choices), class = "latents_bad_argument")
  condition <- tryCatch(.latents_match_arg("nope", choices), error = identity)
  expect_null(conditionCall(condition))
})

test_that("get_results() table arguments matched by name are classed", {
  fit <- surface_lpa()
  expect_error(get_results(fit, "classification", level = "nope"),
               class = "latents_bad_argument")
  expect_error(diagnostics(fit, by = "nope", plots = FALSE),
               class = "latents_bad_argument")
})

test_that("get_candidate() refuses a grid of mixture regressions by class", {
  grid <- enumerate_regressions(score ~ hours, data = study_hours,
                                n_classes = 1:2, n_starts = 1, seed = 1)
  expect_error(get_candidate(grid, 2), class = "latents_bad_argument")
  condition <- tryCatch(get_candidate(grid, 2), error = identity)
  expect_match(conditionMessage(condition), "get_results", fixed = TRUE)
})

test_that("starting_values() and descriptives() refuse a non-fit by class", {
  expect_error(starting_values(1), class = "latents_bad_argument")
  expect_error(starting_values(surface_lpa(), covariance = "nope"),
               class = "latents_bad_argument")
  expect_error(descriptives(1:3), class = "latents_not_a_fit")
})

# API-4 ---------------------------------------------------------------------

test_that("summary() of a profile fit names the arguments it does not use", {
  fit <- surface_lpa()
  # 0.9.21 ignored them without a word; the summary is still built, and the
  # warning names what was not used.
  expect_warning(robust <- summary(fit, vcov_type = "robust"),
                 class = "latents_unused_argument")
  expect_identical(robust, summary(fit))
  expect_warning(summary(fit, level = 0.5), class = "latents_unused_argument")
  expect_warning(summary(fit, foo = 1), "`foo`", class = "latents_unused_argument")
  expect_silent(summary(fit))
})

test_that("summary() of a grid and of a bootstrap test names unused arguments", {
  grid <- enumerate_lpa(srl, names(srl), n_profiles = 1:2, n_starts = 1,
                        seed = 1, model = "EEI")
  expect_warning(summary(grid, levl = 0.5), class = "latents_unused_argument")
  expect_s3_class(suppressWarnings(summary(grid, levl = 0.5)),
                  "summary_multilpa_enumeration")
  ce <- subset(course_engagement, student <= 30)
  null_fit <- multilpa(ce, "browse", "student", n_profiles = 1, n_group_classes = 1,
                       n_starts = 1, seed = 1)
  alternative <- multilpa(ce, "browse", "student", n_profiles = 2,
                          n_group_classes = 1, n_starts = 1, seed = 1)
  test <- bootstrap_lrt(null_fit, alternative, iter = 2, n_starts = 1, seed = 1)
  expect_warning(summary(test, rows = 5), class = "latents_unused_argument")
})

# API-6 ---------------------------------------------------------------------

test_that("summary() of diagnostics and of starting values is the object, not a component dump", {
  fit <- surface_lpa()
  quality <- diagnostics(fit, plots = FALSE)
  expect_s3_class(summary(quality), "multilpa_diagnostics")
  expect_identical(get_results(summary(quality), "checks"),
                   get_results(quality, "checks"))
  printed <- capture.output(print(summary(quality)))
  expect_false(any(grepl("Length +Class +Mode", printed)))
  expect_error(summary(quality, foo = 1), class = "latents_bad_argument")

  start <- starting_values(fit)
  expect_s3_class(summary(start), "multilpa_start")
  expect_identical(get_results(summary(start)), get_results(start))
  expect_false(any(grepl("Length +Class +Mode", capture.output(print(summary(start))))))
})

# API-7 ---------------------------------------------------------------------

test_that("report() and descriptives() name the families they take", {
  regression <- mixture_regression(score ~ hours, data = study_hours,
                                   n_classes = 2, n_starts = 1, seed = 1)
  condition <- tryCatch(report(regression, plots = FALSE), error = identity)
  expect_s3_class(condition, "latents_not_a_fit")
  # The refusal no longer says the fit must be what it already is.
  expect_false(grepl("fitted model of this package", conditionMessage(condition),
                     fixed = TRUE))
  expect_match(conditionMessage(condition), "diagnostics()", fixed = TRUE)
  expect_error(descriptives(regression), class = "latents_not_a_fit")
  expect_error(descriptives(regression), class = "latents_bad_argument")
})

# API-9 ---------------------------------------------------------------------

test_that("diagnostics() of a fit without discrete group classes reports individuals only", {
  fit <- surface_lpa()
  quality <- diagnostics(fit, plots = FALSE)
  expect_setequal(unique(get_results(quality, "classification")$level), "individuals")
  expect_setequal(unique(get_results(quality, "average_posteriors")$level),
                  "individuals")
  expect_false(.multilpa_has_group_classes(fit))
  ce <- subset(course_engagement, student <= 50)
  vars <- c("browse", "lectures", "forum_read")
  one_class <- multilpa(ce, vars, "student", n_profiles = 2, n_group_classes = 1,
                        n_starts = 1, seed = 1)
  expect_false(.multilpa_has_group_classes(one_class))
  expect_setequal(unique(get_results(diagnostics(one_class, plots = FALSE),
                                     "classification")$level), "individuals")
  two_classes <- multilpa(ce, vars, "student", n_profiles = 2, n_group_classes = 2,
                          n_starts = 1, seed = 1)
  expect_true(.multilpa_has_group_classes(two_classes))
  expect_setequal(unique(get_results(diagnostics(two_classes, plots = FALSE),
                                     "classification")$level),
                  c("individuals", "groups"))
})

# API-10 --------------------------------------------------------------------

test_that("the bars view keys its profiles by shape as well as fill", {
  skip_if_not_installed("ggplot2")
  drawn <- plot(surface_lpa(), what = "bars")
  mapped <- unlist(lapply(drawn$layers, function(layer) names(layer$mapping)))
  expect_true("shape" %in% mapped)
  built <- ggplot2::ggplot_build(drawn)
  shapes <- unlist(lapply(built$data, function(layer) layer$shape))
  expect_identical(length(unique(shapes)), 2L)
})

test_that("the trust view states each status in its tile text", {
  status <- c("best replicated", "best found once", "not converged",
              "not converged", "on a boundary", "replication unavailable",
              "not fitted")
  labels <- .latents_trust_label(status, c(3L, 1L, 1L, 3L, 2L, NA, NA))
  expect_identical(labels, c("replicated\n3x", "once\n1x", "not conv.\n1x",
                             "not conv.\n3x", "boundary\n2x", "n/a", "not fitted"))
  # Statuses that share a count no longer share a label.
  expect_false(identical(labels[[2L]], labels[[3L]]))
  skip_if_not_installed("ggplot2")
  grid <- suppressWarnings(enumerate_lpa(srl, names(srl), n_profiles = 1:3,
                                         model = "VVI", n_starts = 3, seed = 1,
                                         max_iter = 4))
  built <- ggplot2::ggplot_build(plot_selection(grid, what = "trust"))
  text <- unlist(lapply(built$data, function(layer) layer$label))
  table <- as.data.frame(grid)
  unconverged <- !(table$converged %in% TRUE) & !is.na(table$log_likelihood)
  expect_true(all(startsWith(text[unconverged], "not conv.")))
})

# INF-6 (bootstrap_lrt) -------------------------------------------------------

test_that("bootstrap_lrt() refuses invalid scalars by class with named messages", {
  ce <- subset(course_engagement, student <= 40)
  vars <- c("browse", "lectures")
  null_fit <- multilpa(ce, vars, "student", n_profiles = 1, n_group_classes = 1,
                       n_starts = 1, seed = 1)
  alternative <- multilpa(ce, vars, "student", n_profiles = 2, n_group_classes = 1,
                          n_starts = 1, seed = 1)
  expect_error(bootstrap_lrt(null_fit, alternative, iter = 1),
               class = "latents_bad_argument")
  condition <- tryCatch(bootstrap_lrt(null_fit, alternative, iter = 2, n_starts = 0),
                        error = identity)
  expect_s3_class(condition, "latents_bad_argument")
  expect_match(conditionMessage(condition), "`n_starts`", fixed = TRUE)
  expect_error(bootstrap_lrt(null_fit, alternative, iter = 2, max_iter = 0),
               class = "latents_bad_argument")
  expect_error(bootstrap_lrt(null_fit, alternative, iter = 2, tol = -1),
               class = "latents_bad_argument")
})

# Handoffs applied after the merge -------------------------------------------

test_that("get_results() refuses a `what` that is not one table name, by class (API-5)", {
  fit <- surface_lpa()
  expect_error(get_results(fit, NA_character_), class = "latents_bad_argument")
  expect_error(get_results(fit, c("profiles", "posteriors")), class = "latents_bad_argument")
  expect_error(get_results(fit, 1), class = "latents_bad_argument")
  expect_error(get_results(fit, "posteriors", format = "nope"), class = "latents_bad_argument")
  expect_error(get_results(fit, "residuals", by = "nope"), class = "latents_bad_argument")
  # Abbreviated table arguments still match, as match.arg() matched them.
  expect_identical(get_results(fit, "posteriors", format = "w"),
                   get_results(fit, "posteriors", format = "wide"))
  regression <- mixture_regression(score ~ hours, study_hours, n_classes = 2,
                                   n_starts = 1, seed = 1)
  expect_error(get_results(regression, NA_character_), class = "latents_bad_argument")
})

test_that("a fit without discrete group classes defaults its classification tables to individuals (API-9)", {
  single <- surface_lpa()
  classification <- get_results(single, "classification")
  expect_identical(unique(classification$level), "individuals")
  expect_identical(classification, get_results(single, "classification", level = "individuals"))
  expect_identical(unique(get_results(single, "bch_weights")$level), "individuals")
  # The group posteriors stay available by name, but leave "all".
  expect_s3_class(get_results(single, "group_posteriors"), "data.frame")
  expect_false("group_posteriors" %in% names(get_results(single, "all")))
  # A two-level fit with group classes keeps both levels.
  two_level <- multilpa(subset(course_engagement, student <= 60),
                        c("browse", "lectures"), "student", n_profiles = 2,
                        n_group_classes = 2, n_starts = 2, seed = 1)
  expect_setequal(unique(get_results(two_level, "classification")$level),
                  c("individuals", "groups"))
  expect_true("group_posteriors" %in% names(get_results(two_level, "all")))
})

test_that("lca() refuses a non-character `categorical` by class (CORE-7)", {
  expect_error(lca(student_esm, c("sports", "walking"), 2, categorical = 1),
               class = "latents_bad_argument")
})

test_that("the summaries that took `...` silently now name what they ignore (API-4)", {
  d <- subset(course_engagement, student <= 40)
  covariate <- quietly(multilpa(d, c("browse", "lectures"), "student", n_profiles = 2,
                                n_group_classes = 1, profile_covariates = "forum_read",
                                n_starts = 1, seed = 1))
  expect_warning(summary(covariate, levl = 0.5), class = "latents_unused_argument")
  transitions <- quietly(lta(engagement_small, c("browse", "lectures"), "student",
                             n_profiles = 2, time = "sequence", n_starts = 1, seed = 1))
  expect_warning(summary(transitions, levl = 0.5), class = "latents_unused_argument")
})

test_that("bootstrap_lrt() refuses a `data` that is not a data frame by class", {
  ce <- subset(course_engagement, student <= 40)
  vars <- c("browse", "lectures")
  null_fit <- multilpa(ce, vars, "student", n_profiles = 1, n_group_classes = 1,
                       n_starts = 1, seed = 1)
  alternative <- multilpa(ce, vars, "student", n_profiles = 2, n_group_classes = 1,
                          n_starts = 1, seed = 1)
  expect_error(bootstrap_lrt(null_fit, alternative, data = as.list(ce), iter = 2),
               class = "latents_bad_data")
})
