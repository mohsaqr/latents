# `get_results()` is the one verb every table now comes from, so what is tested
# here is the contract rather than any one table: that the catalogue and
# `what = "all"` cannot disagree, that a refused argument is named instead of
# dropped, and that a summary serves back exactly what its fit produced.

activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

two_level <- function(...) {
  multilpa(engagement_small, vars = activity, id = "student", n_profiles = 2,
           n_group_classes = 2, n_starts = 2, seed = 1, ...)
}

test_that("every table in a fit's catalogue is a data frame", {
  fit <- two_level()
  catalogue <- setdiff(names(multilpa:::.multilpa_catalogue(fit)), "all")
  expect_gt(length(catalogue), 15L)
  absent <- c("multilpa_no_time", "multilpa_no_group_classes")
  for (what in catalogue) {
    table <- tryCatch(get_results(fit, what), error = function(condition) {
      if (inherits(condition, absent)) NULL else stop(condition)
    })
    if (is.null(table)) next
    expect_s3_class(table, "data.frame")
    expect_false(is.null(names(table)))
  }
})

test_that("`all` is built from the same definitions as a single `what`", {
  fit <- two_level(time = "sequence")
  every <- get_results(fit, "all")
  expect_type(every, "list")
  expect_named(every)
  # Not merely equal: the same definition produced both, so a drift in a
  # default between the two paths would show up here as a difference.
  for (what in names(every)) {
    expect_identical(every[[what]], get_results(fit, what), info = what)
  }
})

test_that("`all` omits only the tables this fit cannot produce", {
  with_time <- get_results(two_level(time = "sequence"), "all")
  without <- get_results(two_level(), "all")
  expect_true(all(c("sequences", "sequence_summary") %in% names(with_time)))
  expect_false(any(c("sequences", "sequence_summary") %in% names(without)))
  expect_identical(setdiff(names(with_time), names(without)),
                   c("sequences", "sequence_summary"))
})

test_that("a table a fit cannot produce still refuses by class when named", {
  fit <- two_level()
  expect_error(get_results(fit, "sequences"), class = "multilpa_no_time")
  expect_error(get_results(fit, "sequence_summary"), class = "multilpa_no_time")
})

test_that("an unknown table names the ones this object has", {
  fit <- two_level()
  expect_error(get_results(fit, "nonsense"), class = "multilpa_bad_argument")
  expect_error(get_results(fit, "nonsense"), "profile_probabilities")
  expect_error(get_results(fit, c("entropy", "profiles")),
               "must be a single table name")
})

test_that("an argument a table does not take is named, not dropped", {
  fit <- two_level()
  expect_error(get_results(fit, "entropy", level = "groups"),
               class = "multilpa_bad_argument")
  expect_error(get_results(fit, "profiles", nonsense = 1),
               class = "multilpa_bad_argument")
  expect_error(get_results(fit, "profiles", "standardized"),
               class = "multilpa_bad_argument")
  expect_error(get_results(fit, "all", data = engagement_small),
               class = "multilpa_bad_argument")
})

test_that("get_results has no method for an unrelated object", {
  expect_error(get_results(data.frame(a = 1)), class = "multilpa_bad_argument")
})

test_that("the classification tables default to every level the fit has", {
  fit <- two_level()
  for (what in c("classification", "average_posteriors",
                 "classification_errors", "bch_weights")) {
    both <- get_results(fit, what)
    expect_setequal(unique(both$level), c("individuals", "groups"))
    expect_identical(both, get_results(fit, what, level = "both"), info = what)
    one <- get_results(fit, what, level = "individuals")
    expect_identical(unique(one$level), "individuals")
    expect_lt(nrow(one), nrow(both))
  }
})

test_that("as.data.frame coerces to the primary table and takes nothing else", {
  fit <- two_level()
  expect_identical(as.data.frame(fit), get_results(fit, "profiles"))
  expect_error(as.data.frame(fit, what = "entropy"),
               class = "multilpa_bad_argument")
  expect_error(as.data.frame(fit, scale = "standardized"),
               class = "multilpa_bad_argument")
  moves <- lta(engagement_small, activity, "student",
                           n_profiles = 2, n_group_classes = 2,
                           time = "sequence", n_starts = 2, seed = 1)
  # The primary table is the measurement model for every fitted family, so a
  # transition fit coerces to its profiles and not to its transition matrix.
  expect_identical(as.data.frame(moves), get_results(moves, "profiles"))
  expect_false(identical(as.data.frame(moves), get_results(moves, "transitions")))
})

test_that("a summary serves back exactly the tables its fit produced", {
  fit <- two_level()
  summarised <- summary(fit)
  from_fit <- get_results(fit, "all")
  expect_identical(get_results(summarised, "all"), from_fit)
  for (what in names(from_fit)) {
    expect_identical(get_results(summarised, what), from_fit[[what]], info = what)
  }
  expect_identical(as.data.frame(summarised), get_results(fit, "profiles"))
})

test_that("the summary print shows every table and honours `rows`", {
  fit <- two_level()
  short <- capture.output(print(summary(fit), rows = 2))
  long <- capture.output(print(summary(fit), rows = 8))
  expect_gt(length(long), length(short))
  for (what in names(get_results(fit, "all"))) {
    expect_true(any(grepl(paste0("-- ", what, " "), short, fixed = TRUE)),
                info = what)
  }
  expect_true(any(grepl("more rows", short, fixed = TRUE)))
  expect_error(print(summary(fit), rows = -1), "non-negative")
})

test_that("the enumeration and bootstrap objects carry their own tables", {
  candidates <- enumerate_classes(engagement_small, activity, "student",
                                  n_profiles = 1:2, n_group_classes = 1,
                                  n_starts = 2, seed = 1)
  expect_identical(as.data.frame(candidates), get_results(candidates, "candidates"))
  expect_s3_class(get_results(candidates, "criteria"), "data.frame")
  expect_named(get_results(candidates, "all"), c("candidates", "criteria"))
  expect_identical(get_results(summary(candidates), "criteria"),
                   get_results(candidates, "criteria"))
})

test_that("a starting-value set and a diagnostics object have catalogues", {
  fit <- two_level()
  start <- starting_values(fit)
  expect_identical(as.data.frame(start), get_results(start, "profiles"))
  expect_named(get_results(start, "all"),
               c("profiles", "responses", "profile_probabilities",
                 "covariances"))
  quality <- diagnostics(fit)
  expect_named(get_results(quality, "all"),
               c("entropy", "classification", "average_posteriors",
                 "residuals"))
  expect_identical(get_results(quality, "entropy"), get_results(fit, "entropy"))
  # The gathered table is named for what it holds, not for the individual
  # posteriors it used to collide with.
  expect_error(get_results(quality, "posteriors"), class = "multilpa_bad_argument")
})

test_that("the tables the summaries used to own are on the fit", {
  fit <- two_level()
  counts <- get_results(fit, "counts")
  expect_named(counts, c("level", "class", "effective_count",
                         "effective_proportion"))
  expect_equal(sum(counts$effective_count[counts$level == "individuals"]),
               fit$n_observations)
  expect_equal(as.vector(tapply(counts$effective_proportion, counts$level, sum)),
               c(1, 1))
  covariances <- get_results(fit, "covariances")
  expect_equal(nrow(covariances),
               fit$n_profiles * length(activity)^2)
  # The diagonal parameterization states exact zeros off the diagonal rather
  # than leaving them missing.
  off_diagonal <- covariances$covariance[
    covariances$indicator != covariances$indicator_2]
  expect_true(all(off_diagonal == 0))
  model <- get_results(fit, "model")
  expect_equal(nrow(model), 1L)
  expect_equal(model$log_likelihood, fit$log_likelihood)
  expect_equal(model$n_profiles, fit$n_profiles)
  expect_equal(model$bic_groups, fit$bic)
})

test_that("`model` reports the columns of the family, not a padded union", {
  covariates <- multilpa(engagement_small, activity, "student", n_profiles = 2,
                         n_group_classes = 2,
                         profile_covariates = "previous_grade",
                         n_starts = 2, seed = 1)
  expect_true("n_profile_covariates" %in% names(get_results(covariates, "model")))
  expect_false("n_profile_covariates" %in% names(get_results(two_level(), "model")))
  # The one table whose columns are the same for every family, which is what
  # makes it the one to compare fits across families with.
  expect_identical(names(get_results(covariates, "information_criteria")),
                   names(get_results(two_level(), "information_criteria")))
  # Both families spell the two conventions the same way.
  for (fit in list(covariates, two_level())) {
    expect_true(all(c("bic_groups", "bic_individual") %in%
                      names(get_results(fit, "model"))))
  }
})

test_that("truth cross-tabulates against the level each column describes", {
  fit <- two_level()
  recovery <- get_results(fit, "assignments", data = engagement_small,
                       truth = c("engagement", "student_type"))
  expect_named(recovery, c("assignment", "class", "truth", "value", "n",
                           "proportion"))
  # `engagement` varies within a student, so it is an enrolment-level truth;
  # `student_type` does not, so it is a student-level one.
  expect_identical(unique(recovery$assignment[recovery$truth == "engagement"]),
                   "profile")
  expect_identical(unique(recovery$assignment[recovery$truth == "student_type"]),
                   "group_class")
  # Observation-level truth counts enrolments; group-level truth counts students.
  expect_equal(sum(recovery$n[recovery$truth == "engagement"]),
               nrow(engagement_small))
  expect_equal(sum(recovery$n[recovery$truth == "student_type"]),
               length(unique(engagement_small$student)))
  # The proportions are shares within a truth value, so they sum to one.
  shares <- as.vector(tapply(recovery$proportion,
                             paste(recovery$truth, recovery$value), sum))
  expect_equal(shares, rep(1, length(shares)))
})

test_that("group recovery gives each group one vote and omits unused truth levels", {
  fit <- two_level()
  data <- engagement_small
  data$known_group <- factor(data$student_type,
                             levels = c(as.character(unique(data$student_type)),
                                        "unused"))
  joined <- get_results(fit, "assignments", data = data)
  groups <- joined[!duplicated(joined$student), ]
  recovery <- get_results(fit, "assignments", data = data,
                          truth = "known_group")
  expected <- table(groups$group_class, groups$known_group)
  expect_identical(unique(recovery$assignment), "group_class")
  expect_equal(sum(recovery$n), nrow(groups))
  expect_false(any(as.character(recovery$value) == "unused"))
  expect_true(all(is.finite(recovery$proportion)))
  for (class in seq_len(fit$n_group_classes)) {
    for (value in setdiff(colnames(expected), "unused")) {
      observed <- recovery$n[recovery$class == class & recovery$value == value]
      expect_equal(observed, as.integer(expected[as.character(class), value]))
    }
  }
})

test_that("truth reproduces the cross-tabulation it replaces", {
  fit <- two_level()
  joined <- get_results(fit, "assignments", data = engagement_small)
  expected <- table(joined$profile, joined$engagement)
  recovery <- get_results(fit, "assignments", data = engagement_small,
                       truth = "engagement")
  taken <- recovery$n[recovery$class == 1L & recovery$value == "engaged"]
  expect_equal(taken, as.integer(expected["1", "engaged"]))
  expect_equal(sum(recovery$n), sum(expected))
})

test_that("truth refuses a column it cannot use", {
  fit <- two_level()
  expect_error(get_results(fit, "assignments", data = engagement_small,
                        truth = "absent"), class = "multilpa_bad_data")
  expect_error(get_results(fit, "assignments", data = engagement_small,
                        truth = "profile"), class = "multilpa_bad_data")
  expect_error(get_results(fit, "assignments", data = engagement_small,
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
  # that has three rows. The long table is still what `get_results()` returns.
  indicators <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  fit <- quietly(multilpa(engagement_small, indicators, "student",
                          n_profiles = 3L, n_group_classes = 1L, n_starts = 3L,
                          seed = 1L, max_iter = 5000))
  wide <- multilpa:::.multilpa_wide_means(fit)
  expect_identical(nrow(wide), 3L)
  expect_identical(names(wide), c("profile", indicators))

  # Reshaped from the long table, so the two cannot report different numbers.
  long <- get_results(fit, "profiles")
  expect_equal(wide$forum_post, long$mean[long$indicator == "forum_post"])
  expect_identical(nrow(long), 15L)

  printed <- utils::capture.output(print(fit))
  expect_true(any(grepl("attendance", printed, fixed = TRUE)))
  # How big each profile is, beside what it looks like, matched on the label.
  expect_true(any(grepl("proportion", printed, fixed = TRUE)))
  sizes <- get_results(fit, "counts")
  individuals <- sizes[sizes$level == "individuals", , drop = FALSE]
  # The count is printed to the column's precision; its value truncated to two
  # decimals is a prefix of the printed number whatever its size.
  shown <- formatC(trunc(individuals$effective_count[1L] * 100) / 100,
                   format = "f", digits = 2)
  expect_true(any(grepl(shown, printed, fixed = TRUE)))
  # The spread columns belong to the long table, and are pointed at, not shown.
  expect_false(any(grepl("standard_deviation", printed, fixed = TRUE)))
  expect_true(any(grepl("get_results(x, \"profiles\")", printed, fixed = TRUE)))
})
