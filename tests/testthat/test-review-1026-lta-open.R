# Regression tests for the three open latent transition findings of the
# 2026-10-26 review (validation/review-2026-10-26/area-reports/lta.md):
# LTA-F1, a rank-deficient covariate design was fitted and its unidentified
# logits counted; LTA-F2, the logits of occasions after the last observed one
# were counted; LTA-F3, an occasion with nothing observed entered the
# transition network with its prior-modal profile.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

open_lta_data <- function(n = 40L, n_occasions = 4L, seed = 11L) {
  withr::with_seed(seed, {
    do.call(rbind, lapply(seq_len(n), function(i) {
      states <- Reduce(function(previous, t) {
        if (stats::runif(1L) < 0.75) previous else 3L - previous
      }, seq_len(n_occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
      data.frame(id = i, time = seq_len(n_occasions), z = stats::rnorm(1L),
                 w = stats::rnorm(n_occasions),
                 y1 = stats::rnorm(n_occasions, c(0, 2)[states], 0.6),
                 y2 = stats::rnorm(n_occasions, c(0, 1.5)[states], 0.6))
    }))
  })
}

open_lta_fit <- function(data, ...) {
  quietly(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
              n_starts = 2, seed = 1, ...))
}

# Hand count for K = 2 profiles, two Gaussian indicators with varying
# variances (8 measurement parameters per block), one group class: initial
# logits p0 * (K - 1), transition logits p * (K - 1) * K.
open_lta_count <- function(p0, p, blocks = 1L) {
  blocks * 8L + p0 * 1L + p * 1L * 2L
}

# LTA-F1 ---------------------------------------------------------------------

test_that("a rank-deficient initial or transition design is refused", {
  data <- open_lta_data()
  data$constant <- 0.1
  data$occasion_number <- data$time
  expect_error(open_lta_fit(data, transition_covariates = "constant"),
               class = "latents_bad_covariate")
  expect_error(open_lta_fit(data, initial_covariates = "constant"),
               class = "latents_bad_covariate")
  # A covariate that takes one value per occasion is the occasion intercepts
  # over again (the report's 17 counted vs 15 identified parameters).
  expect_error(open_lta_fit(data, transitions = "occasion",
                            transition_covariates = "occasion_number"),
               class = "latents_bad_covariate")
  # Order 2: the first-order block governs only the move into occasion 2, so a
  # covariate constant at occasion 2 is unidentified there, though it varies
  # over the later moves.
  data$late <- ifelse(data$time == 2L, 1, data$w)
  expect_error(open_lta_fit(data, order = 2L, transition_covariates = "late"),
               class = "latents_bad_covariate")
})

test_that("the rank decision depends on the pattern of a design, not its units", {
  data <- open_lta_data()
  # Rank deficient in any units: refused at every scale and origin.
  invisible(lapply(c(1e-8, 1, 1e8), function(scale) {
    data$occasion_number <- data$time * scale + 3
    expect_error(open_lta_fit(data, transitions = "occasion",
                              transition_covariates = "occasion_number"),
                 class = "latents_bad_covariate")
    data$constant <- scale
    expect_error(open_lta_fit(data, initial_covariates = "constant"),
                 class = "latents_bad_covariate")
  }))
  # Full rank in any units and at any origin (a covariate held at 1e8 is not
  # the intercept): accepted, with the same count.
  fits <- Map(function(scale, origin) {
    data$scaled <- data$w * scale + origin
    open_lta_fit(data, transition_covariates = "scaled")
  }, c(1e-6, 1e6, 1), c(5, -3, 1e8))
  expect_equal(vapply(fits, `[[`, integer(1), "n_parameters"),
               rep(open_lta_count(p0 = 1L, p = 2L), 3L))
  # The likelihood is compared at moderate origins only: at 1e8 the fit itself
  # (not the rank check) reaches a lower maximum, reported separately.
  expect_equal(fits[[1L]]$log_likelihood, fits[[2L]]$log_likelihood, tolerance = 1e-6)
})

# LTA-F2 ---------------------------------------------------------------------

test_that("an lta() with no observed indicator after the first occasion is refused", {
  data <- open_lta_data()
  data[data$time > 1L, c("y1", "y2")] <- NA
  expect_error(open_lta_fit(data, transitions = "occasion", missing = "fiml"),
               class = "latents_bad_transition")
})

test_that("an all-missing trailing occasion leaves the count and the criteria unchanged", {
  data <- open_lta_data()
  # Ten groups gain a fifth occasion at which nothing is observed. The
  # likelihood is flat in the move into it, so it adds no parameter.
  trailing <- data[data$time == 4L & data$id <= 10L, ]
  trailing$time <- 5L
  trailing[c("y1", "y2")] <- NA
  padded <- rbind(data, trailing)
  base <- open_lta_fit(data, transitions = "occasion", missing = "fiml")
  extended <- open_lta_fit(padded, transitions = "occasion", missing = "fiml")
  expect_identical(base$n_parameters, open_lta_count(p0 = 1L, p = 3L))
  expect_identical(extended$n_parameters, base$n_parameters)
  expect_equal(extended$log_likelihood, base$log_likelihood, tolerance = 1e-6)
  expect_equal(c(extended$aic, extended$bic), c(base$aic, base$bic), tolerance = 1e-6)
  # Occasion-specific measurement: the empty occasion's measurement block is
  # flat as well (every row there has log density 0).
  base_occasion <- open_lta_fit(data, transitions = "occasion", missing = "fiml",
                                measurement = "occasion")
  extended_occasion <- open_lta_fit(padded, transitions = "occasion", missing = "fiml",
                                    measurement = "occasion")
  expect_identical(base_occasion$n_parameters,
                   open_lta_count(p0 = 1L, p = 3L, blocks = 4L))
  expect_identical(extended_occasion$n_parameters, base_occasion$n_parameters)
})

# LTA-F3 ---------------------------------------------------------------------

open_sequence_fit <- function(data) {
  quietly(multilpa(data, vars = c("y1", "y2"), id = "id", time = "time",
                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
                   missing = "fiml"))
}

test_that("a network from a fit without an ordering is still refused", {
  # Guards the rewritten sequence builder, which now reads the ordering itself.
  data <- open_lta_data()
  fit <- quietly(multilpa(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                          n_group_classes = 1, n_starts = 1, seed = 1))
  expect_error(.multilpa_modal_sequences(fit), class = "latents_no_time")
})

test_that("an occasion with nothing observed is treated as an absent one", {
  data <- open_lta_data(n = 40L, n_occasions = 4L, seed = 81L)
  emptied <- data$time == 2L & data$id %in% c(3L, 7L, 12L)
  with_empty <- data
  with_empty[emptied, c("y1", "y2")] <- NA
  without <- data[!emptied, ]
  fit_empty <- open_sequence_fit(with_empty)
  fit_absent <- open_sequence_fit(without)
  # The empty rows add nothing to the likelihood, so the two fits agree.
  expect_equal(fit_empty$log_likelihood, fit_absent$log_likelihood, tolerance = 1e-8)
  expect_identical(fit_empty$subject_profiles[!emptied], fit_absent$subject_profiles)
  sequences_empty <- .multilpa_modal_sequences(fit_empty)
  sequences_absent <- .multilpa_modal_sequences(fit_absent)
  # No state for an occasion with nothing observed: one state per classified row.
  expect_identical(sum(!is.na(as.matrix(sequences_empty$states))),
                   sum(!emptied))
  expect_identical(sequences_empty, sequences_absent)
  skip_if_not_installed("tna")
  expect_equal(get_tna(fit_empty)$weights, get_tna(fit_absent)$weights)
})

test_that("the sequence table leaves out an occasion with nothing observed", {
  data <- open_lta_data(n = 40L, n_occasions = 4L, seed = 81L)
  emptied <- data$time == 2L & data$id %in% c(3L, 7L, 12L)
  data[emptied, c("y1", "y2")] <- NA
  fit <- open_sequence_fit(data)
  long <- get_results(fit, "sequences")
  expect_identical(nrow(long), sum(!emptied))
  expect_false(any(long$group %in% c(3L, 7L, 12L) & long$time == 2L))
  wide <- get_results(fit, "sequences", format = "wide")
  expect_true(all(is.na(wide[wide$group %in% c(3L, 7L, 12L), "time_2"])))
  # The convention of 0.9.37 and earlier keeps every row.
  kept <- quietly(multilpa(data, vars = c("y1", "y2"), id = "id", time = "time",
                           n_profiles = 2, n_group_classes = 1, n_starts = 2,
                           seed = 1, missing = "fiml", drop_empty = FALSE))
  expect_identical(nrow(get_results(kept, "sequences")), nrow(data))
})

# LTA-F2, empty occasion followed by an observed one --------------------------

test_that("an empty occasion before an observed one warns under occasion transitions", {
  data <- open_lta_data()
  data[data$time == 2L, c("y1", "y2")] <- NA
  expect_warning(
    withCallingHandlers(
      lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
          transitions = "occasion", missing = "fiml", n_starts = 1, seed = 1),
      warning = function(condition) {
        if (!inherits(condition, "latents_unidentified_occasion") &&
            inherits(condition, .multilpa_expected_warnings)) {
          invokeRestart("muffleWarning")
        }
      }),
    class = "latents_unidentified_occasion")
  # Homogeneous transitions identify the one matrix from the other moves.
  expect_no_warning(
    quietly(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
                missing = "fiml", n_starts = 1, seed = 1)),
    class = "latents_unidentified_occasion")
})

# ML-2 on the transition path -------------------------------------------------

test_that("zero-weight rows cannot mask a constant indicator in lta()", {
  data <- open_lta_data()
  data$sampling <- ifelse(data$id <= 5L, 0, 1)
  # y2 varies only in the rows that carry no weight.
  data$y2[data$id > 5L] <- 1
  expect_error(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                   time = "time", weights = "sampling", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  # The same data without the zero-weight groups is refused the same way:
  # the decision does not depend on rows that do not enter the likelihood.
  expect_error(lta(subset(data, sampling > 0), vars = c("y1", "y2"), id = "id",
                   n_profiles = 2, time = "time", weights = "sampling",
                   n_starts = 1, seed = 1),
               class = "latents_bad_data")
})

# LTA-F5 ------------------------------------------------------------------------

# The logit coefficients under covariates c_new = unit * c + origin: the slope
# becomes slope / unit and the intercept intercept - origin * slope / unit.
open_chart_map <- function(names_all, covariates, unit, origin) {
  map <- diag(length(names_all))
  pattern <- sprintf("\\.(%s)$", paste(covariates, collapse = "|"))
  slopes <- grep(pattern, names_all)
  intercepts <- match(sub(pattern, ".(Intercept)", names_all[slopes]), names_all)
  map[cbind(slopes, slopes)] <- 1 / unit
  map[cbind(intercepts, slopes)] <- -origin / unit
  map
}

test_that("covariate logits do not depend on a covariate's origin or units", {
  data <- open_lta_data()
  fit_at <- function(unit, origin) {
    data$moving <- data$w * unit + origin
    data$starting <- data$z * unit + origin
    open_lta_fit(data, transition_covariates = "moving",
                 initial_covariates = "starting")
  }
  reference <- fit_at(1, 0)
  estimates <- coef(reference)
  covariance <- vcov(reference)
  invisible(Map(function(unit, origin) {
    fit <- fit_at(unit, origin)
    expect_true(fit$converged)
    expect_lt(abs(fit$log_likelihood - reference$log_likelihood),
              1e-8 * abs(reference$log_likelihood))
    map <- open_chart_map(names(estimates), c("moving", "starting"), unit, origin)
    expected <- drop(map %*% estimates)
    expect_lt(max(abs(coef(fit) - expected) / (1 + abs(expected))), 1e-6)
    expected_se <- sqrt(diag(map %*% covariance %*% t(map)))
    expect_lt(max(abs(sqrt(diag(vcov(fit))) - expected_se) / expected_se), 1e-4)
  }, c(1, 1e-4), c(1e8, 0)))
})

test_that("the sequence plot marks an occasion that observes nothing as its own legend entry", {
  skip_if_not_installed("ggplot2")
  withr::local_seed(3)
  panel <- data.frame(person = rep(1:12, each = 5), beep = rep(1:5, 12))
  state <- (panel$person + panel$beep) %% 2
  panel$a <- 2 * state + stats::rnorm(60, sd = 0.6)
  panel$b <- -1.5 * state + stats::rnorm(60, sd = 0.7)
  panel[panel$person == 3 & panel$beep == 3, c("a", "b")] <- NA
  fit <- quietly(lta(panel, c("a", "b"), "person", n_profiles = 2, time = "beep",
                     n_starts = 2, seed = 1, missing = "fiml"))
  built <- ggplot2::ggplot_build(plot(fit, what = "sequences"))
  tiles <- built$data[[1L]]
  # One grey cell for the empty occasion, drawn in the reserved grey, and the
  # 59 observed occasions in profile colours.
  expect_identical(sum(tiles$fill == "#999999"), 1L)
  expect_identical(nrow(tiles), 60L)
  labels <- unlist(lapply(built$data, function(layer) layer$label))
  expect_true("–" %in% labels)
  fills <- ggplot2::ggplot_build(plot(fit, what = "sequences"))$plot$scales$get_scales("fill")
  expect_true("Nothing observed" %in% names(fills$palette.cache %||% fills$get_breaks() |> stats::setNames(fills$get_breaks())))
})
