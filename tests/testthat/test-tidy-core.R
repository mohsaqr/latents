# The public surface of the two entry points: their argument contracts, the
# classes their failures carry, and the tidy shape of what comes back.
#
# Error paths are asserted by CLASS, never by message text: the class is the
# contract documented on ?"multilpa-conditions", the message is not.

.core_fixture <- function(n_groups = 20L, per_group = 10L, seed = 5L) {
  set.seed(seed)
  n <- n_groups * per_group
  data <- data.frame(
    school = rep(seq_len(n_groups), each = per_group),
    kind = rep(c("mixed", "high"), each = n / 2L))
  data$high <- stats::rbinom(n, 1L, ifelse(data$kind == "high", 0.85, 0.15))
  data$a <- stats::rnorm(n, mean = 2 * data$high)
  data$b <- stats::rnorm(n, mean = 2 * data$high)
  data$c <- ifelse(stats::rbinom(n, 1L, ifelse(data$high == 1L, 0.8, 0.2)) == 1L,
                   "yes", "no")
  data
}

test_that("multilpa() gives every match.arg argument a vector default", {
  formals_used <- formals(multilpa)
  # Each of these is consumed by match.arg(), which requires the full set of
  # permitted values as the default.
  expect_identical(eval(formals_used$variance_model), c("varying", "equal"))
  expect_identical(eval(formals_used$missing), c("error", "fiml"))
  expect_identical(eval(formals_used$covariance_model), c("diagonal", "full"))
  staged_formals <- formals(fit_staged)
  expect_identical(eval(staged_formals$variance_model), c("varying", "equal"))
  expect_identical(eval(staged_formals$missing), c("error", "fiml"))
  expect_identical(eval(staged_formals$covariance_model), c("diagonal", "full"))
})

test_that("match.arg rejects an unlisted choice for every such argument", {
  data <- .core_fixture()
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        variance_model = "free"))
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        missing = "drop"))
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        covariance_model = "banded"))
})

test_that("a broken data contract is refused by class, not by message", {
  data <- .core_fixture()
  # An absent indicator, a duplicated indicator, and a group column that is
  # also an indicator all break the same contract.
  expect_error(multilpa(data, c("a", "absent"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
  expect_error(multilpa(data, c("a", "a"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
  expect_error(multilpa(data, c("a", "school"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
  # A factor is not a continuous indicator.
  factored <- transform(data, a = factor(sample(c("p", "q"), nrow(data), TRUE)))
  expect_error(multilpa(factored, c("a", "b"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
  # A constant indicator identifies nothing.
  flat <- transform(data, a = 1)
  expect_error(multilpa(flat, c("a", "b"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
  # Missing values are rejected under the default missing = "error".
  gappy <- data
  gappy$a <- replace(gappy$a, 1L, NA_real_)
  expect_error(multilpa(gappy, c("a", "b"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
  # A missing group identifier.
  ungrouped <- data
  ungrouped$school <- replace(ungrouped$school, 1L, NA_integer_)
  expect_error(multilpa(ungrouped, c("a", "b"), "school", 2L, 1L, n_starts = 1),
               class = "multilpa_bad_data")
})

test_that("a scalar argument outside its contract is refused by class", {
  data <- .core_fixture()
  expect_error(multilpa(data, c("a", "b"), "school", 0L, 1L, n_starts = 1),
               class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("a", "b"), "school", 2.5, 1L, n_starts = 1),
               class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 0L),
               class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        max_iter = -1L), class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        tol = 0), class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        min_variance = -1), class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        seed = -1), class = "multilpa_bad_argument")
})

test_that("an unidentifiable request is refused by class", {
  data <- .core_fixture()
  # More group classes than groups.
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 100L, n_starts = 1),
               class = "multilpa_unidentified")
  # Several group classes with a single profile.
  expect_error(multilpa(data, c("a", "b"), "school", 1L, 2L, n_starts = 1),
               class = "multilpa_unidentified")
  # Only singleton groups.
  singletons <- transform(data, school = seq_len(nrow(data)))
  expect_error(multilpa(singletons, c("a", "b"), "school", 2L, 2L, n_starts = 1),
               class = "multilpa_unidentified")
  # More profiles than distinct observed indicator rows.
  tiny <- data.frame(school = rep(1:4, each = 2),
                     a = rep(c(0, 1), 4), b = rep(c(0, 1), 4))
  expect_error(multilpa(tiny, c("a", "b"), "school", 5L, 1L, n_starts = 1),
               class = "multilpa_unidentified")
})

test_that("the data contract is checked before `time` is", {
  data <- .core_fixture()
  data$occasion <- rep(seq_len(10L), times = 20L)
  # `time` is validated against the group column, so a group column that does
  # not exist has to be reported as the data-contract failure it is, not as an
  # opaque failure inside split() during the time check.
  expect_error(multilpa(data, c("a", "b"), "absent", 2L, 1L, n_starts = 1,
                        time = "occasion"), class = "multilpa_bad_data")
})

test_that("`time` may not name an indicator, as its contract says", {
  data <- .core_fixture()
  # The contract message has always promised this; only the `id` column was
  # actually excluded, so `time = "a"` used to be accepted silently.
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        time = "a"),
               "must not be one of `vars` or the `id` column")
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        time = "school"),
               "must not be one of `vars` or the `id` column")
})

test_that("a repeated or missing position is refused by class", {
  data <- .core_fixture()
  data$occasion <- rep(seq_len(10L), times = 20L)
  repeated <- data
  repeated$occasion <- replace(repeated$occasion, 2L, 1L)
  expect_error(multilpa(repeated, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        time = "occasion"), class = "multilpa_bad_time")
  gappy <- data
  gappy$occasion <- replace(gappy$occasion, 1L, NA_integer_)
  expect_error(multilpa(gappy, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        time = "occasion"), class = "multilpa_bad_time")
})

test_that("a malformed start is refused by class", {
  data <- .core_fixture()
  fit <- quietly(multilpa(data, c("a", "b"), "school", 2L, 1L,
                                   n_starts = 2, seed = 1))
  start <- starting_values(fit)
  wrong_shape <- start
  wrong_shape$means <- matrix(0, 3L, 2L)
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        start = wrong_shape), class = "multilpa_bad_start")
  missing_block <- start
  missing_block$variances <- NULL
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 1L, n_starts = 1,
                        start = missing_block), class = "multilpa_bad_start")
  negative <- start
  negative$group_probabilities <- c(2, -1)
  expect_error(multilpa(data, c("a", "b"), "school", 2L, 2L, n_starts = 1,
                        start = negative), class = "multilpa_bad_start")
})

test_that("a qualified fit warns by class and never tells the caller to use $", {
  data <- .core_fixture(n_groups = 8L, per_group = 4L, seed = 9L)
  seen <- list()
  fit <- withCallingHandlers(
    multilpa(data, c("a", "b"), "school", 3L, 2L, n_starts = 2, seed = 4,
             max_iter = 1L),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(fit, "multilpa")
  # This fixture is chosen to qualify the fit; if it stops warning the test is
  # no longer checking anything and must be told so.
  expect_gt(length(seen), 0L)
  messages <- vapply(seen, conditionMessage, character(1))
  classes <- vapply(seen, function(w) class(w)[[1L]], character(1))
  # Every warning multilpa() raises is catchable by its own class ...
  expect_true(all(startsWith(classes, "multilpa_")))
  # ... and none of them directs the caller at a `$` on the result.
  expect_false(any(grepl("$", messages, fixed = TRUE)))
})

test_that("an unconverged fit raises a catchable multilpa_unconverged warning", {
  data <- .core_fixture()
  expect_warning(
    multilpa(data, c("a", "b"), "school", 2L, 2L, n_starts = 1, seed = 1,
             max_iter = 1L),
    class = "multilpa_unconverged")
})

test_that("the starts diagnostics are reachable as a tidy frame, not by $", {
  data <- .core_fixture()
  fit <- quietly(multilpa(data, c("a", "b"), "school", 2L, 2L,
                                   n_starts = 3, seed = 1))
  starts <- get_results(fit, "starts")
  expect_s3_class(starts, "data.frame")
  expect_identical(nrow(starts), 3L)
  expect_true(all(c("start", "log_likelihood", "converged", "iterations",
                    "error", "boundary") %in% names(starts)))
})

test_that("multilpa() returns a tidy frame for each `what`", {
  data <- .core_fixture()
  fit <- quietly(multilpa(data, c("a", "b", "c"), "school", 2L, 2L,
                                   categorical = "c", n_starts = 3, seed = 1))
  profiles <- as.data.frame(fit)
  expect_s3_class(profiles, "data.frame")
  # One row per profile x continuous indicator.
  expect_identical(nrow(profiles), 4L)
  weights <- get_results(fit, "profile_probabilities")
  expect_s3_class(weights, "data.frame")
  # One row per group class x profile, never a K-wide matrix.
  expect_identical(nrow(weights), 4L)
  expect_true(all(abs(rowSums(matrix(weights$probability, 2L, 2L,
                                     byrow = TRUE)) - 1) < 1e-8))
})

test_that("fit_staged() reaches the staged workflow in one call", {
  data <- .core_fixture()
  staged <- fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                       n_group_classes = 2L, n_starts = 3, seed = 1)
  expect_s3_class(staged, "multilpa")
  expect_true(isTRUE(staged$staged))
  stages <- get_results(staged, "stages")
  expect_s3_class(stages, "data.frame")
  expect_identical(nrow(stages), 2L)
  expect_identical(stages$stage, c("measurement", "membership"))
  # The held measurement is exactly the first stage's, not a nearby refit.
  expect_equal(staged$means, staged$stage_one$means)
  expect_equal(staged$variances, staged$stage_one$variances)
})

test_that("fit_staged() refuses a staging request it cannot meet", {
  data <- .core_fixture()
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                          n_group_classes = 1L, n_starts = 2, seed = 1))
  # "3" >= 2L compares as strings and used to pass silently.
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                          n_group_classes = "3", n_starts = 2, seed = 1))
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                          n_group_classes = 2.5, n_starts = 2, seed = 1))
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                          n_group_classes = 2L, n_starts = 2, seed = 1,
                          measurement = "not a fit"))
})

test_that("the categorical and utility helpers stay internal", {
  exported <- getNamespaceExports("multilpa")
  # These return bare matrices and lists by design, which is only acceptable
  # because nothing outside the package can call them.
  internal <- c(".multilpa_encode_categorical", ".multilpa_categorical_log_density",
                ".multilpa_categorical_maximize", ".multilpa_bound_probabilities",
                ".multilpa_categorical_initialize", ".multilpa_categorical_design",
                ".multilpa_categorical_parameters", ".multilpa_categorical_thresholds",
                ".multilpa_any_fit")
  expect_true(all(vapply(internal, exists, logical(1),
                         envir = asNamespace("multilpa"), inherits = FALSE)))
  expect_length(intersect(internal, exported), 0L)
  # Nothing in the namespace may be exported under a private `.multilpa_` name.
  expect_length(grep("^\\.multilpa_", exported, value = TRUE), 0L)
})

test_that("every condition class raised in R/ is documented, and every documented class is raised", {
  skip_on_cran()
  # A built package does not ship R/, so this is a source-tree check only.
  source_dir <- file.path(testthat::test_path("..", ".."), "R")
  skip_if_not(dir.exists(source_dir), "R/ sources are not available")
  files <- list.files(source_dir, pattern = "[.]R$", full.names = TRUE)

  raised <- unlist(lapply(files, function(path) {
    lines <- readLines(path, warn = FALSE)
    text <- paste(ifelse(grepl("^\\s*[^#]", lines), lines, ""), collapse = "\n")
    positions <- gregexpr('class\\s*=\\s*(c\\()?\\s*"(multilpa[A-Za-z0-9_]*)"',
                          text)[[1L]]
    if (positions[1L] == -1L) return(character(0))
    found <- regmatches(text, gregexpr(
      'class\\s*=\\s*(c\\()?\\s*"(multilpa[A-Za-z0-9_]*)"', text))[[1L]]
    last_at <- function(pattern, before) {
      hits <- gregexpr(pattern, substr(text, 1L, before), perl = TRUE)[[1L]]
      if (hits[1L] == -1L) 0L else max(hits)
    }
    # A `class =` literal belongs to the constructor that opened most recently
    # before it; `structure(x, class = "multilpa_covariates")` is an S3 object
    # class, not a condition.
    is_condition <- vapply(positions, function(s) {
      last_at("(error|warning|simple|message)Condition\\(", s) >
        last_at("(structure\\(|class\\s*<-)", s)
    }, logical(1))
    sub('.*"(multilpa[A-Za-z0-9_]*)".*', "\\1", found[is_condition])
  }))

  doc_text <- paste(readLines(file.path(source_dir, "conditions.R"),
                              warn = FALSE), collapse = "\n")
  items <- regmatches(doc_text, gregexpr("\\\\item\\{[^}]*\\}", doc_text))[[1L]]
  documented <- unique(gsub("`", "", unlist(regmatches(items,
    gregexpr("`multilpa_[A-Za-z0-9_]+`", items)))))

  expect_setequal(unique(raised), documented)
})

test_that("every warning the package raises carries a class", {
  skip_on_cran()
  source_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(source_dir), "package sources not available")
  # `?"multilpa-conditions"` states the classes are the contract and the
  # messages are not. A bare `warning()` breaks that promise silently: a caller
  # muffling an expected qualification can only match message text, which
  # changes between versions, so it ends up muffling every warning instead.
  # A bare `stop()` remains allowed for one-off internal guards, which no
  # caller is expected to catch; only warnings are enforced here.
  offenders <- unlist(lapply(list.files(source_dir, pattern = "\\.R$",
                                        full.names = TRUE), function(path) {
    lines <- readLines(path, warn = FALSE)
    bare <- grepl("(^|[^a-zA-Z._])warning\\s*\\(\\s*[\"']", lines) &
      !grepl("^\\s*#", lines)
    if (!any(bare)) return(NULL)
    sprintf("%s:%d", basename(path), which(bare))
  }), use.names = FALSE)
  expect_identical(offenders, NULL)
})
