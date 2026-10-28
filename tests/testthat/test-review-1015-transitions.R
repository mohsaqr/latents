# Review 2026-10-15 of the transition-model files: regressions for the
# defects that review confirmed (each failed before its fix).

review_panel <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(8)
  n <- 60L
  occasions <- 4L
  w <- stats::rnorm(n)
  rows <- do.call(rbind, lapply(seq_len(n), function(j) {
    z <- stats::rnorm(occasions)
    state <- integer(occasions)
    state[1L] <- 1L + (stats::runif(1L) < stats::plogis(0.3 + 0.8 * w[j]))
    for (t in seq_len(occasions)[-1L]) {  # simulation only
      move <- stats::plogis(-1 + 0.9 * z[t])
      state[t] <- if (stats::runif(1L) < move) 3L - state[t - 1L] else state[t - 1L]
    }
    data.frame(id = j, time = seq_len(occasions), z = z, w = w[j],
               y = stats::rnorm(occasions, c(-1, 1)[state], 0.7))
  }))
  rows$my.z <- rows$z
  rows
})

review_general_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- lta(review_panel, "y", "id", n_profiles = 2, time = "time",
                     transition_covariates = "z", n_starts = 1, seed = 1)
    }
    cached
  }
})

test_that("a dotted covariate name gives the coefficient table without a warning", {
  fit <- lta(review_panel, "y", "id", n_profiles = 2, time = "time",
             transition_covariates = "my.z", n_starts = 1, seed = 1)
  expect_no_warning(table <- get_results(fit, "transition_coefficients"))
  expect_setequal(table$term, c("(Intercept)", "my.z"))
  expect_equal(nrow(table), 4L)
})

test_that("get_results() on a general transition fit refuses an argument it does not use", {
  # Refused, as every other family refuses it, rather than warned about and
  # dropped (release review API-4).
  fit <- review_general_fit()
  expect_error(get_results(fit, "transitions", estimated = FALSE),
               class = "latents_bad_argument")
  expect_no_warning(get_results(fit, "transitions", level = 0.9))
})

test_that("vcov() warns when given the other family's covariance argument", {
  general <- review_general_fit()
  expect_warning(ignored <- vcov(general, vcov_type = "robust"),
                 class = "latents_unused_argument")
  expect_identical(ignored, vcov(general, type = "observed"))
  homogeneous <- lta(review_panel, "y", "id", n_profiles = 2, time = "time",
                     n_starts = 1, seed = 1)
  expect_warning(ignored <- vcov(homogeneous, type = "robust"),
                 class = "latents_unused_argument")
  expect_identical(ignored, vcov(homogeneous, vcov_type = "observed"))
  expect_false(identical(ignored, vcov(homogeneous, vcov_type = "robust")))
})

test_that("a one-class general fit prints 'group class' in the singular", {
  expect_output(print(review_general_fit()), "1 group class \\(")
})

test_that("an emptied stayer class is reported as a small class", {
  fit <- NULL
  expect_warning(
    fit <- lta(review_panel, "y", "id", n_profiles = 2, time = "time",
               mover_stayer = TRUE, n_starts = 2, seed = 1),
    class = "latents_small_classes")
  # The fixture's stayer class is empty: the warning describes this fit.
  expect_lt(min(colSums(fit$group_posteriors)), 1)
})

test_that("a transition enumeration's get_results() names an argument it ignores", {
  grid <- quietly(enumerate_classes(
    review_panel, "y", "id", n_profiles = 2, time = "time", n_starts = 1, seed = 1))
  expect_warning(table <- get_results(grid, "candidates", sort = TRUE),
                 class = "latents_unused_argument")
  expect_identical(table, get_results(grid, "candidates"))
})

# Follow-up decisions (2026-10-15): refusals and tables the lead asked for.

test_that("a transition covariate at a skipped grid position is refused", {
  gappy <- review_panel[!(review_panel$id %in% c(3, 7) & review_panel$time == 2), ]
  expect_error(lta(gappy, "y", "id", n_profiles = 2, time = "time", occasions = "grid",
                   transition_covariates = "z", n_starts = 1, seed = 1),
               class = "latents_bad_data", regexp = "group 3 at 2, group 7 at 2")
  # Either remedy the message names fits; a gap without covariates still fits.
  expect_s3_class(lta(gappy, "y", "id", n_profiles = 2, time = "time",
                      transition_covariates = "z", n_starts = 1, seed = 1),
                  "multilpa_lta")
  expect_s3_class(lta(gappy, "y", "id", n_profiles = 2, time = "time", occasions = "grid",
                      transitions = "occasion", n_starts = 1, seed = 1),
                  "multilpa_lta")
  # The skipped rows added back with missing indicators: the covariate is then
  # observed, and shifting it by a constant only moves the intercepts, so the
  # likelihood is unchanged (the zero fill broke this invariance).
  filled <- review_panel
  filled$y[filled$id %in% c(3, 7) & filled$time == 2] <- NA
  shifted <- transform(filled, z = z + 5)
  fits <- lapply(list(filled, shifted), \(d) lta(d, "y", "id", n_profiles = 2,
    time = "time", occasions = "grid", missing = "fiml", transition_covariates = "z",
    n_starts = 1, seed = 1))
  expect_equal(fits[[1L]]$log_likelihood, fits[[2L]]$log_likelihood, tolerance = 1e-6)
})

test_that("numbers stored as text are refused when their alphabetical order is wrong", {
  long <- local({
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, globalenv()), add = TRUE)
    set.seed(4)
    do.call(rbind, lapply(1:30, function(j) {
      data.frame(id = j, t = 1:11, y = stats::rnorm(11, rep(c(-1.5, 1.5), c(5, 6)), 0.6))
    }))
  })
  numeric_fit <- lta(long, "y", "id", n_profiles = 2, time = "t", n_starts = 1, seed = 1)
  long$text <- as.character(long$t)
  expect_error(lta(long, "y", "id", n_profiles = 2, time = "text", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(lta(long, "y", "id", n_profiles = 2, time = "text", n_starts = 1, seed = 1,
                   transitions = "occasion"), class = "latents_bad_data")
  long$padded <- sprintf("%02d", long$t)
  padded_fit <- lta(long, "y", "id", n_profiles = 2, time = "padded", n_starts = 1, seed = 1)
  expect_equal(padded_fit$log_likelihood, numeric_fit$log_likelihood)
  long$date <- format(as.Date("2024-01-01") + 30 * (long$t - 1))
  date_fit <- lta(long, "y", "id", n_profiles = 2, time = "date", n_starts = 1, seed = 1)
  expect_equal(date_fit$log_likelihood, numeric_fit$log_likelihood)
})

test_that("a categorical general transition fit reports its response probabilities", {
  panel <- local({
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, globalenv()), add = TRUE)
    set.seed(5)
    within(review_panel, {
      grade <- factor(ifelse(y + stats::rnorm(length(y), 0, 0.8) > 0, "high", "low"))
    })
  })
  fit <- lta(panel, c("y", "grade"), "id", n_profiles = 2, time = "time",
             categorical = "grade", transition_covariates = "z", n_starts = 1, seed = 1)
  table <- get_results(fit, "responses")
  expect_named(table, c("occasion", "profile", "indicator", "category", "probability"))
  expect_equal(nrow(table), 4L)
  expect_setequal(table$category, c("high", "low"))
  expect_equal(as.vector(tapply(table$probability, table$profile, sum)), c(1, 1))
  expect_true("responses" %in% names(get_results(fit, "all")))
  expect_true("responses" %in% names(get_results(summary(fit), "all")))
  # Occasion-specific measurement: one block per occasion.
  by_occasion <- lta(panel, c("y", "grade"), "id", n_profiles = 2, time = "time",
                     categorical = "grade", measurement = "occasion", n_starts = 1, seed = 1)
  expect_equal(sort(unique(get_results(by_occasion, "responses")$occasion)), 1:4)
  # A fit without categorical indicators has no such table.
  expect_false("responses" %in% names(get_results(review_general_fit(), "all")))
})

test_that("a weighted covariate fit averages its model-implied tables with the weights", {
  weighted <- transform(review_panel, weight = 1 + (id %% 3 == 0) * 4 + (w > 0) * 2)
  fit <- lta(weighted, "y", "id", n_profiles = 2, time = "time", weights = "weight",
             transition_covariates = "z", initial_covariates = "w", n_starts = 1, seed = 1)
  coefficients <- coef(fit)
  unit <- weighted[weighted$time == 1, ]
  # Initial: logit of profile 1 against profile 2 (the last), on w.
  initial <- stats::plogis(coefficients[["initial.group_class_1.profile_1.(Intercept)"]] +
                             coefficients[["initial.group_class_1.profile_1.w"]] * unit$w)
  table <- get_results(fit, "initial")
  expect_equal(table$probability[table$profile == "profile_1"],
               sum(unit$weight * initial) / sum(unit$weight), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(table$probability[table$profile == "profile_1"],
                                mean(initial))))
  # Moves into occasion 3 from profile 1: log odds of leaving against staying, on z.
  arriving <- weighted[weighted$time == 3, ]
  leave <- stats::plogis(
    coefficients[["transition.group_class_1.profile_1->profile_2.(Intercept)"]] +
      coefficients[["transition.group_class_1.profile_1->profile_2.z"]] * arriving$z)
  moves <- subset(get_results(fit, "transitions"), occasion == 3L & from == "profile_1" &
                    to == "profile_2")
  expect_equal(moves$probability, sum(arriving$weight * leave) / sum(arriving$weight),
               tolerance = 1e-12)
})

test_that("the general transition plot uses the shared probability ramp and legible ink", {
  skip_if_not_installed("ggplot2")
  fit <- review_general_fit()
  plot <- plot(fit)
  built <- ggplot2::ggplot_build(plot)
  fill <- built$plot$scales$get_scales("fill")
  expect_identical(fill$map(c(0, 1)), c("#FFFFFF", "#0072B2"))
  text <- ggplot2::layer_data(plot, 2L)
  expect_identical(text$colour, .gg_ink(.gg_probability_fill(plot$data$probability)))
})

test_that("multilpa(time = ) refuses numbers stored as text in the wrong order", {
  long <- data.frame(id = rep(1:6, each = 11L), t = rep(1:11, 6L),
                     y = rep(c(rep(-2, 5L), rep(2, 6L)), 6L) + sin(seq_len(66L)) / 4)
  long$text <- as.character(long$t)
  expect_error(multilpa(long, "y", "id", n_profiles = 2, time = "text", n_starts = 1,
                        seed = 1), class = "latents_bad_time")
  long$padded <- sprintf("%02d", long$t)
  # Six sequences: a class with less than one member is expected.
  fit <- quietly(multilpa(long, "y", "id", n_profiles = 2, time = "padded", n_starts = 1,
                          seed = 1))
  expect_identical(unique(get_results(fit, "sequences")$time), sprintf("%02d", 1:11))
})
