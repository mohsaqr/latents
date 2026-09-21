# `get_data()` is the one verb every table now comes from, so what is tested
# here is the contract rather than any one table: that the catalogue and
# `what = "all"` cannot disagree, that a refused argument is named instead of
# dropped, and that a summary serves back exactly what its fit produced.

activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

two_level <- function(...) {
  multilpa(course_engagement, vars = activity, id = "student", n_profiles = 2,
           n_group_classes = 2, n_starts = 2, seed = 1, ...)
}

test_that("every table in a fit's catalogue is a data frame", {
  fit <- two_level()
  catalogue <- setdiff(names(multilpa:::.multilpa_catalogue(fit)), "all")
  expect_gt(length(catalogue), 15L)
  absent <- c("multilpa_no_time", "multilpa_no_group_classes")
  for (what in catalogue) {
    table <- tryCatch(get_data(fit, what), error = function(condition) {
      if (inherits(condition, absent)) NULL else stop(condition)
    })
    if (is.null(table)) next
    expect_s3_class(table, "data.frame")
    expect_false(is.null(names(table)))
  }
})

test_that("`all` is built from the same definitions as a single `what`", {
  fit <- two_level(time = "sequence")
  every <- get_data(fit, "all")
  expect_type(every, "list")
  expect_named(every)
  # Not merely equal: the same definition produced both, so a drift in a
  # default between the two paths would show up here as a difference.
  for (what in names(every)) {
    expect_identical(every[[what]], get_data(fit, what), info = what)
  }
})

test_that("`all` omits only the tables this fit cannot produce", {
  with_time <- get_data(two_level(time = "sequence"), "all")
  without <- get_data(two_level(), "all")
  expect_true(all(c("sequences", "sequence_summary") %in% names(with_time)))
  expect_false(any(c("sequences", "sequence_summary") %in% names(without)))
  expect_identical(setdiff(names(with_time), names(without)),
                   c("sequences", "sequence_summary"))
})

test_that("a table a fit cannot produce still refuses by class when named", {
  fit <- two_level()
  expect_error(get_data(fit, "sequences"), class = "multilpa_no_time")
  expect_error(get_data(fit, "sequence_summary"), class = "multilpa_no_time")
  intercepts <- fit_random_intercept(course_engagement, activity, "student",
                                     n_profiles = 2, n_starts = 1, seed = 1)
  expect_error(get_data(intercepts, "residuals"),
               class = "multilpa_no_group_classes")
  # Three-step tables refused an intercept fit with an unclassed error before
  # 0.11.2, which `all` could not tell apart from a defect.
  expect_error(get_data(intercepts, "classification_errors"),
               class = "multilpa_no_group_classes")
  expect_error(get_data(intercepts, "bch_weights"),
               class = "multilpa_no_group_classes")
  expect_false("residuals" %in% names(get_data(intercepts, "all")))
})

test_that("an unknown table names the ones this object has", {
  fit <- two_level()
  expect_error(get_data(fit, "nonsense"), class = "multilpa_bad_argument")
  expect_error(get_data(fit, "nonsense"), "profile_probabilities")
  expect_error(get_data(fit, c("entropy", "profiles")),
               "must be a single table name")
})

test_that("an argument a table does not take is named, not dropped", {
  fit <- two_level()
  expect_error(get_data(fit, "entropy", level = "groups"),
               class = "multilpa_bad_argument")
  expect_error(get_data(fit, "profiles", nonsense = 1),
               class = "multilpa_bad_argument")
  expect_error(get_data(fit, "profiles", "standardized"),
               class = "multilpa_bad_argument")
  expect_error(get_data(fit, "all", data = course_engagement),
               class = "multilpa_bad_argument")
})

test_that("get_data has no method for an unrelated object", {
  expect_error(get_data(data.frame(a = 1)), class = "multilpa_bad_argument")
})

test_that("the classification tables default to every level the fit has", {
  fit <- two_level()
  intercepts <- fit_random_intercept(course_engagement, activity, "student",
                                     n_profiles = 2, n_starts = 1, seed = 1)
  for (what in c("classification", "average_posteriors",
                 "classification_errors", "bch_weights")) {
    both <- get_data(fit, what)
    expect_setequal(unique(both$level), c("individuals", "groups"))
    expect_identical(both, get_data(fit, what, level = "both"), info = what)
    one <- get_data(fit, what, level = "individuals")
    expect_identical(unique(one$level), "individuals")
    expect_lt(nrow(one), nrow(both))
  }
  # A family with one level reports that level rather than refusing.
  expect_identical(unique(get_data(intercepts, "classification")$level),
                   "individuals")
})

test_that("a level a fit has not is refused rather than silently dropped", {
  intercepts <- fit_random_intercept(course_engagement, activity, "student",
                                     n_profiles = 2, n_starts = 1, seed = 1)
  expect_error(get_data(intercepts, "classification", level = "groups"),
               class = "multilpa_no_group_classes")
})

test_that("as.data.frame coerces to the primary table and takes nothing else", {
  fit <- two_level()
  expect_identical(as.data.frame(fit), get_data(fit, "profiles"))
  expect_error(as.data.frame(fit, what = "entropy"),
               class = "multilpa_bad_argument")
  expect_error(as.data.frame(fit, scale = "standardized"),
               class = "multilpa_bad_argument")
  moves <- fit_transitions(course_engagement, activity, "student",
                           n_profiles = 2, n_group_classes = 2,
                           time = "sequence", n_starts = 2, seed = 1)
  # The primary table is the measurement model for every fitted family, so a
  # transition fit coerces to its profiles and not to its transition matrix.
  expect_identical(as.data.frame(moves), get_data(moves, "profiles"))
  expect_false(identical(as.data.frame(moves), get_data(moves, "transitions")))
})

test_that("a summary serves back exactly the tables its fit produced", {
  fit <- two_level()
  summarised <- summary(fit)
  from_fit <- get_data(fit, "all")
  expect_identical(get_data(summarised, "all"), from_fit)
  for (what in names(from_fit)) {
    expect_identical(get_data(summarised, what), from_fit[[what]], info = what)
  }
  expect_identical(as.data.frame(summarised), get_data(fit, "profiles"))
})

test_that("the summary print shows every table and honours `rows`", {
  fit <- two_level()
  short <- capture.output(print(summary(fit), rows = 2))
  long <- capture.output(print(summary(fit), rows = 8))
  expect_gt(length(long), length(short))
  for (what in names(get_data(fit, "all"))) {
    expect_true(any(grepl(paste0("-- ", what, " "), short, fixed = TRUE)),
                info = what)
  }
  expect_true(any(grepl("more rows", short, fixed = TRUE)))
  expect_error(print(summary(fit), rows = -1), "non-negative")
})

test_that("the enumeration and bootstrap objects carry their own tables", {
  candidates <- enumerate_classes(course_engagement, activity, "student",
                                  n_profiles = 1:2, n_group_classes = 1,
                                  n_starts = 2, seed = 1)
  expect_identical(as.data.frame(candidates), get_data(candidates, "candidates"))
  expect_s3_class(get_data(candidates, "criteria"), "data.frame")
  expect_named(get_data(candidates, "all"), c("candidates", "criteria"))
  expect_identical(get_data(summary(candidates), "criteria"),
                   get_data(candidates, "criteria"))
})

test_that("a starting-value set and a diagnostics object have catalogues", {
  fit <- two_level()
  start <- starting_values(fit)
  expect_identical(as.data.frame(start), get_data(start, "profiles"))
  expect_named(get_data(start, "all"),
               c("profiles", "responses", "profile_probabilities",
                 "covariances"))
  quality <- diagnostics(fit)
  expect_named(get_data(quality, "all"),
               c("entropy", "classification", "average_posteriors",
                 "residuals"))
  expect_identical(get_data(quality, "entropy"), get_data(fit, "entropy"))
  # The gathered table is named for what it holds, not for the individual
  # posteriors it used to collide with.
  expect_error(get_data(quality, "posteriors"), class = "multilpa_bad_argument")
})

test_that("the tables the summaries used to own are on the fit", {
  fit <- two_level()
  counts <- get_data(fit, "counts")
  expect_named(counts, c("level", "class", "effective_count",
                         "effective_proportion"))
  expect_equal(sum(counts$effective_count[counts$level == "individuals"]),
               fit$n_observations)
  expect_equal(as.vector(tapply(counts$effective_proportion, counts$level, sum)),
               c(1, 1))
  covariances <- get_data(fit, "covariances")
  expect_equal(nrow(covariances),
               fit$n_profiles * length(activity)^2)
  # The diagonal parameterization states exact zeros off the diagonal rather
  # than leaving them missing.
  off_diagonal <- covariances$covariance[
    covariances$indicator != covariances$indicator_2]
  expect_true(all(off_diagonal == 0))
  model <- get_data(fit, "model")
  expect_equal(nrow(model), 1L)
  expect_equal(model$log_likelihood, fit$log_likelihood)
  expect_equal(model$n_profiles, fit$n_profiles)
  expect_equal(model$bic_groups, fit$bic)
})

test_that("`model` reports the columns of the family, not a padded union", {
  intercepts <- fit_random_intercept(course_engagement, activity, "student",
                                     n_profiles = 2, n_starts = 1, seed = 1)
  covariates <- multilpa(course_engagement, activity, "student", n_profiles = 2,
                         n_group_classes = 2,
                         profile_covariates = "previous_grade",
                         n_starts = 2, seed = 1)
  expect_true("quadrature_check_passed" %in% names(get_data(intercepts, "model")))
  expect_true("n_profile_covariates" %in% names(get_data(covariates, "model")))
  expect_false("quadrature_check_passed" %in% names(get_data(two_level(), "model")))
  # The one table whose columns are the same for every family, which is what
  # makes it the one to compare fits across families with.
  expect_identical(names(get_data(intercepts, "information_criteria")),
                   names(get_data(two_level(), "information_criteria")))
  # Both families spell the two conventions the same way.
  for (fit in list(intercepts, covariates, two_level())) {
    expect_true(all(c("bic_groups", "bic_individual") %in%
                      names(get_data(fit, "model"))))
  }
})

test_that("truth cross-tabulates against the level each column describes", {
  fit <- two_level()
  recovery <- get_data(fit, "assignments", data = course_engagement,
                       truth = c("engagement", "student_type"))
  expect_named(recovery, c("assignment", "class", "truth", "value", "n",
                           "proportion"))
  # `engagement` varies within a student, so it is an enrolment-level truth;
  # `student_type` does not, so it is a student-level one.
  expect_identical(unique(recovery$assignment[recovery$truth == "engagement"]),
                   "profile")
  expect_identical(unique(recovery$assignment[recovery$truth == "student_type"]),
                   "group_class")
  # Every enrolment is counted once per truth column.
  expect_equal(sum(recovery$n[recovery$truth == "engagement"]),
               nrow(course_engagement))
  # The proportions are shares within a truth value, so they sum to one.
  shares <- as.vector(tapply(recovery$proportion,
                             paste(recovery$truth, recovery$value), sum))
  expect_equal(shares, rep(1, length(shares)))
})

test_that("truth reproduces the cross-tabulation it replaces", {
  fit <- two_level()
  joined <- get_data(fit, "assignments", data = course_engagement)
  expected <- table(joined$profile, joined$engagement)
  recovery <- get_data(fit, "assignments", data = course_engagement,
                       truth = "engagement")
  taken <- recovery$n[recovery$class == 1L & recovery$value == "engaged"]
  expect_equal(taken, as.integer(expected["1", "engaged"]))
  expect_equal(sum(recovery$n), sum(expected))
})

test_that("truth refuses a column it cannot use", {
  fit <- two_level()
  expect_error(get_data(fit, "assignments", data = course_engagement,
                        truth = "absent"), class = "multilpa_bad_data")
  expect_error(get_data(fit, "assignments", data = course_engagement,
                        truth = "profile"), class = "multilpa_bad_data")
  expect_error(get_data(fit, "assignments", data = course_engagement,
                        truth = c("engagement", "engagement")),
               "must not be duplicated")
})

test_that("plot draws every supported view under what = \"all\"", {
  fit <- two_level()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE, after = FALSE)
  drawn <- expect_invisible(plot(fit, what = "all"))
  expect_s3_class(drawn, "multilpa")
  # Evaluated where the package namespace is *not* on the search path, which is
  # every ordinary user session: "all" re-issues the caller's call, and a call
  # naming the unexported method would be unresolvable there.
  clean <- new.env(parent = globalenv())
  assign("fit", fit, envir = clean)
  expect_s3_class(evalq(plot(fit, what = "all"), clean), "multilpa")
  # "all" is a request, not a view, so it must not appear among the views it
  # walks: it would otherwise call itself.
  expect_false("all" %in% multilpa:::.multilpa_supported_views(fit))
})

test_that("a fit prints its means one row per profile, not one per cell", {
  # Three profiles over five indicators is fifteen console lines for a table
  # that has three rows. The long table is still what `get_data()` returns.
  indicators <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  fit <- quietly(multilpa(course_engagement, indicators, "student",
                          n_profiles = 3L, n_group_classes = 1L, n_starts = 3L,
                          seed = 1L, max_iter = 5000))
  wide <- multilpa:::.multilpa_wide_means(fit)
  expect_identical(nrow(wide), 3L)
  expect_identical(names(wide), c("profile", indicators))

  # Reshaped from the long table, so the two cannot report different numbers.
  long <- get_data(fit, "profiles")
  expect_equal(wide$forum_post, long$mean[long$indicator == "forum_post"])
  expect_identical(nrow(long), 15L)

  printed <- utils::capture.output(print(fit))
  expect_true(any(grepl("attendance", printed, fixed = TRUE)))
  # How big each profile is, beside what it looks like, matched on the label.
  expect_true(any(grepl("proportion", printed, fixed = TRUE)))
  sizes <- get_data(fit, "counts")
  individuals <- sizes[sizes$level == "individuals", , drop = FALSE]
  expect_true(any(grepl(format(individuals$effective_count[1L], digits = 7),
                        printed, fixed = TRUE)))
  # The spread columns belong to the long table, and are pointed at, not shown.
  expect_false(any(grepl("standard_deviation", printed, fixed = TRUE)))
  expect_true(any(grepl("get_data(x, \"profiles\")", printed, fixed = TRUE)))
})
