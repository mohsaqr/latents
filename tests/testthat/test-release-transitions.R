# Regressions for the pre-release review of the transition models (cluster
# "transitions"): TRANS-1, -3, -4, -5, -7, -8, -9, INF-5 and the transition
# parts of API-2, API-4, API-5 and API-6. Each test fails on the code before
# the fix.

# Two well separated profiles over up to four occasions, unequal lengths.
release_panel <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(14)
  moves <- matrix(c(0.8, 0.2, 0.3, 0.7), 2L, byrow = TRUE)
  do.call(rbind, lapply(seq_len(50L), function(i) {
    length_i <- sample(2:4, 1L)
    state <- Reduce(function(previous, step) sample(1:2, 1L, prob = moves[previous, ]),
                    seq_len(length_i - 1L), sample(1:2, 1L), accumulate = TRUE)
    data.frame(id = i, time = seq_len(length_i),
               y1 = stats::rnorm(length_i, c(0, 2)[state]),
               y2 = stats::rnorm(length_i, c(0, 1.5)[state], 0.8),
               state = state, wt = 1)
  }))
})

# Three units of weight zero observed at six occasions, beyond every weighted
# unit's last occasion: absent from the likelihood by the package's rule.
release_zero_weight <- local({
  extra <- do.call(rbind, lapply(1001:1003, function(i) {
    data.frame(id = i, time = seq_len(6L), y1 = (i - 1002) + seq_len(6L) / 10,
               y2 = -(i - 1002) - seq_len(6L) / 7, state = 1L, wt = 0)
  }))
  rbind(release_panel, extra)
})

# A third profile, 40 units away, occupied only at the last occasion of every
# third group (and by single-occasion groups): no move ever leaves it.
release_empty_origin <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(31)
  panel <- do.call(rbind, lapply(seq_len(45L), function(i) {
    state <- c(sample(1:2, 3L, replace = TRUE), if (i %% 3L == 0L) 3L else sample(1:2, 1L))
    data.frame(id = i, time = 1:4, y1 = stats::rnorm(4L, c(0, 40, 80)[state]),
               y2 = stats::rnorm(4L, c(0, 40, 80)[state]), one = 1, w = stats::rnorm(4L))
  }))
  rbind(panel, data.frame(id = 101:108, time = 1, y1 = stats::rnorm(8L, 80),
                          y2 = stats::rnorm(8L, 80), one = 1, w = stats::rnorm(8L)))
})

# The value of `expr` and the classes of the warnings it raised, muffled.
with_warnings <- function(expr) {
  raised <- character()
  value <- withCallingHandlers(expr, warning = function(condition) {
    raised <<- c(raised, class(condition)[1L])
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = raised)
}

release_fit <- function(data, ...) {
  quietly(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
              n_starts = 2, seed = 1, ...))
}

test_that("an origin no move leaves is flagged and not counted by both engines (TRANS-1)", {
  skip_on_cran()
  fit_with <- function(...) {
    lta(release_empty_origin, vars = c("y1", "y2"), id = "id", n_profiles = 3,
        time = "time", n_starts = 2, seed = 1, ...)
  }
  homogeneous <- with_warnings(fit_with())
  expect_true("latents_empty_transition_row" %in% homogeneous$warnings)
  homogeneous <- homogeneous$value
  # 6 means + 6 variances + 2 initial logits + 4 informed transition logits.
  expect_identical(as.integer(homogeneous$n_parameters), 18L)
  general <- with_warnings(fit_with(weights = "one"))
  expect_true("latents_empty_transition_row" %in% general$warnings)
  general <- general$value
  expect_identical(as.integer(general$n_parameters), 18L)
  expect_equal(general$log_likelihood, homogeneous$log_likelihood, tolerance = 1e-6)
  expect_equal(general$aic, -2 * general$log_likelihood + 2 * 18)
  # The empty row is not its start values but equal probabilities, flagged.
  transitions <- get_results(general, "transitions")
  empty <- subset(transitions, from == "profile_3")
  expect_equal(empty$probability, rep(1 / 3, 3L))
  expect_false(any(empty$estimated))
  expect_true(all(subset(transitions, from != "profile_3")$estimated))
  # Wald inference withholds only the flat logits.
  inference <- parameter_inference(general)
  flat <- subset(inference, block == "transition" & from == "profile_3")
  expect_false(any(flat$estimated))
  expect_true(all(is.na(flat$standard_error)))
  informed <- subset(inference, block == "transition" & from != "profile_3")
  expect_true(all(informed$estimated))
  expect_true(all(is.finite(informed$standard_error)))
  # The same with a transition covariate: its slope out of the empty origin
  # is as flat as the intercept, so 2 x 2 fewer logits than 26.
  covariate <- with_warnings(fit_with(transition_covariates = "w"))
  expect_true("latents_empty_transition_row" %in% covariate$warnings)
  expect_identical(as.integer(covariate$value$n_parameters), 22L)
})

test_that("zero-weight units beyond the weighted occasions leave an occasion-measurement fit unchanged (TRANS-3)", {
  skip_on_cran()
  reference <- release_fit(release_panel, weights = "wt", measurement = "occasion")
  extended <- with_warnings(
    lta(release_zero_weight, vars = c("y1", "y2"), id = "id", n_profiles = 2,
        time = "time", n_starts = 2, seed = 1, weights = "wt", measurement = "occasion"))
  # An occasion only units of weight zero reach is absent, not thin.
  expect_false("latents_thin_profiles" %in% extended$warnings)
  extended <- extended$value
  expect_equal(extended$log_likelihood, reference$log_likelihood, tolerance = 1e-8)
  expect_identical(extended$n_parameters, reference$n_parameters)
})

test_that("the transitions table leaves out moves only zero-weight units make (TRANS-4)", {
  skip_on_cran()
  reference <- release_fit(release_panel, weights = "wt", transitions = "occasion")
  extended <- release_fit(release_zero_weight, weights = "wt", transitions = "occasion")
  expect_equal(extended$log_likelihood, reference$log_likelihood, tolerance = 1e-8)
  table <- get_results(extended, "transitions")
  expect_false(anyNA(table$probability))
  expect_identical(sort(unique(table$occasion)), 2:4)
  expect_equal(table$probability, get_results(reference, "transitions")$probability,
               tolerance = 1e-6)
  expect_true(all(table$estimated))
  # Homogeneous transitions averaged over occasions are not poisoned either.
  pooled <- release_fit(release_zero_weight, weights = "wt")
  expect_false(anyNA(get_results(pooled, "transitions")$probability))
})

test_that("get_results() on a general transition fit keeps the shared contract (TRANS-5, API-4, API-5)", {
  skip_on_cran()
  general <- release_fit(release_panel, weights = "wt")
  homogeneous <- release_fit(release_panel)
  expect_error(get_results(general, "nope"), class = "latents_bad_argument")
  expect_error(get_results(general, "classification"), class = "latents_bad_argument")
  # Not one name at all: missing, several or not text.
  expect_error(get_results(general, NA_character_), class = "latents_bad_argument")
  expect_error(get_results(general, c("transitions", "profiles")),
               class = "latents_bad_argument")
  expect_error(get_results(general, 1), class = "latents_bad_argument")
  expect_identical(get_results(general, NULL), get_results(general, "transitions"))
  expect_error(get_results(general, "transitions", foo = 1), class = "latents_bad_argument")
  expect_error(as.data.frame(general, foo = 1), class = "latents_bad_argument")
  expect_error(get_results(general, "transitions", level = 2), class = "latents_bad_argument")
  expect_error(parameter_inference(general, level = 2), class = "latents_bad_argument")
  # One assignments layout, whichever engine fitted the model.
  assigned <- get_results(general, "assignments")
  # The weighted fit's data columns include its weights, as for multilpa().
  expect_identical(setdiff(names(assigned), "wt"),
                   names(get_results(homogeneous, "assignments")))
  expect_identical(assigned$profile, get_results(homogeneous, "assignments")$profile)
  expect_type(assigned$profile, "integer")
  expect_identical(assigned$id, release_panel$id)
  joined <- get_results(general, "assignments", data = release_panel)
  expect_identical(joined$state, release_panel$state)
  recovery <- get_results(general, "assignments", data = release_panel, truth = "state")
  expect_equal(recovery, get_results(homogeneous, "assignments", data = release_panel,
                                     truth = "state"))
  expect_error(get_results(general, "assignments",
                           data = release_panel[rev(seq_len(nrow(release_panel))), ]),
               class = "latents_bad_inference_data")
  expect_true("iterations" %in% names(get_results(general, "starts")))
})

test_that("summary() of a general transition fit is a summary, not a list dump (API-6)", {
  skip_on_cran()
  general <- release_fit(release_panel, transitions = "occasion")
  # The occasion fit sits on a probability bound, which the coefficient
  # tables announce; that notice is not the subject here.
  summarised <- quietly(summary(general),
                        c(.multilpa_expected_warnings, "latents_no_standard_errors"))
  expect_s3_class(summarised, "summary_multilpa_lta")
  expect_identical(get_results(summarised, "transitions"),
                   get_results(general, "transitions"))
  expect_identical(as.data.frame(summarised), get_results(general, "transitions"))
  expect_error(get_results(summarised, "nope"), class = "latents_bad_argument")
  expect_error(get_results(summarised, NA_character_), class = "latents_bad_argument")
  expect_error(get_results(summarised, c("transitions", "profiles")),
               class = "latents_bad_argument")
  expect_error(summary(general, foo = 1), class = "latents_bad_argument")
  printed <- capture.output(print(summarised))
  expect_false(any(grepl("^\\$", printed)))
  expect_true(any(grepl("-- transitions", printed, fixed = TRUE)))
})

test_that("an unused factor level is dropped from transition covariates (TRANS-7)", {
  skip_on_cran()
  data <- transform(release_panel,
                    group = factor(ifelse(id %% 2 == 0, "a", "b"), levels = c("a", "b", "c")))
  unused <- release_fit(data, initial_covariates = "group")
  dropped <- release_fit(transform(data, group = droplevels(group)),
                         initial_covariates = "group")
  expect_equal(unused$log_likelihood, dropped$log_likelihood, tolerance = 1e-10)
  expect_identical(unused$n_parameters, dropped$n_parameters)
})

test_that("numbered text labels that sort out of order are refused (TRANS-8)", {
  labelled <- function(labels) {
    expect_null(.lta_check_time_order(labels, "time"))
  }
  expect_error(.lta_check_time_order(paste0("T", 1:10), "time"), class = "latents_bad_time")
  expect_error(.lta_check_time_order(paste("week", 1:10), "time"), class = "latents_bad_time")
  expect_error(.lta_check_time_order(paste0(1:10, "h"), "time"), class = "latents_bad_data")
  labelled(paste0("T", 1:9))
  labelled(sprintf("T%02d", 1:10))
  labelled(c("pre", "post", "follow-up"))
  labelled(c("2024-01-05", "2024-11-05"))
  labelled(c("a1", "b10", "c2"))
  long <- do.call(rbind, lapply(seq_len(4L), function(i) {
    data.frame(id = i, time = paste0("T", 1:10), y1 = seq_len(10) + i, y2 = 10:1 * i)
  }))
  expect_error(lta(long, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time"),
               class = "latents_bad_time")
})

test_that("the TRANS-8 refusal's two remedies both give the numeric-time fit", {
  skip_on_cran()
  # "T3", "T6", "T9", "T12": "T12" sorts first, so the labels are refused.
  labels <- paste0("T", release_panel$time * 3L)
  expect_error(release_fit(transform(release_panel, time = labels)),
               class = "latents_bad_time")
  numeric_time <- release_fit(release_panel)
  ordered <- release_fit(transform(release_panel, time = factor(
    labels, levels = paste0("T", c(3L, 6L, 9L, 12L)))))
  padded <- release_fit(transform(release_panel,
                                  time = sprintf("T%02d", release_panel$time * 3L)))
  expect_equal(ordered$log_likelihood, numeric_time$log_likelihood, tolerance = 1e-10)
  expect_equal(padded$log_likelihood, numeric_time$log_likelihood, tolerance = 1e-10)
})

test_that("bootstrap_lrt() refuses transition fits with different drop_empty (TRANS-9)", {
  skip_on_cran()
  null_model <- release_fit(release_panel)
  alternative <- release_fit(release_panel, transitions = "occasion", drop_empty = FALSE)
  expect_error(bootstrap_lrt(null_model, alternative, data = release_panel, iter = 1),
               class = "latents_incomparable_models")
})

test_that("confint() of a transition fit refuses a parm that is not an estimate (INF-5)", {
  skip_on_cran()
  fit <- quietly(lta(engagement_small, c("browse", "lectures"), "student",
                     n_profiles = 2, time = "sequence", n_starts = 2, seed = 1))
  bounds <- confint(fit)
  expect_identical(confint(fit, parm = 2), bounds[2L, , drop = FALSE])
  expect_identical(confint(fit, parm = rownames(bounds)[3L]), bounds[3L, , drop = FALSE])
  expect_error(confint(fit, parm = 99), class = "latents_bad_argument")
  expect_error(confint(fit, parm = 0), class = "latents_bad_argument")
  expect_error(confint(fit, parm = "no_such_parameter"), class = "latents_bad_argument")
})

test_that("a transition grid raises an input error every candidate shares (API-2)", {
  skip_on_cran()
  expect_error(enumerate_classes(engagement_small, c("browse", "lectures"), "student",
                                 n_profiles = 2:3, time = "nope", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  # A failure of one candidate stays in the grid, and its print says so.
  grid <- quietly(enumerate_classes(release_panel, c("y1", "y2"), "id",
                                    n_profiles = c(2, 1000), time = "time",
                                    n_starts = 1, seed = 1))
  expect_identical(sum(!is.na(as.data.frame(grid)$error)), 1L)
  expect_output(print(grid), "1 candidate(s) failed to fit", fixed = TRUE)
})

test_that("a transition or profile grid raises an argument error every candidate shares (API-2 handoff)", {
  skip_on_cran()
  vars <- c("browse", "lectures")
  expect_error(enumerate_classes(engagement_small, vars, "student", n_profiles = 2:3,
                                 time = "sequence", weights = "nope", n_starts = 1,
                                 seed = 1),
               class = "latents_bad_weights")
  expect_error(enumerate_classes(engagement_small, vars, "student", n_profiles = 2:3,
                                 time = "sequence", min_variance = -1, n_starts = 1,
                                 seed = 1),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(engagement_small, vars, "student", n_profiles = 1:2,
                                 min_variance = -1, n_starts = 1, seed = 1),
               class = "latents_bad_argument")
})

test_that("a general transition fit classes an unknown table or covariance type (INF-2, API-5)", {
  skip_on_cran()
  general <- release_fit(release_panel, weights = "wt")
  expect_error(get_results(general, "classification_errors"), class = "latents_bad_argument")
  expect_error(get_results(general, "bch_weights"), class = "latents_bad_argument")
  expect_error(get_results(general, "transitions", vcov_type = "sandwich"),
               class = "latents_bad_argument")
  expect_error(parameter_inference(general, vcov_type = "sandwich"),
               class = "latents_bad_argument")
  expect_error(parameter_inference(general, adjust = "nope"), class = "latents_bad_argument")
  homogeneous <- release_fit(release_panel)
  expect_error(bootstrap_lrt(homogeneous, homogeneous, data = release_panel, iter = 0),
               class = "latents_bad_argument")
})
