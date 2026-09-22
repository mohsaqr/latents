.surface_data <- function(seed = 5L) {
  set.seed(seed)
  school <- rep(seq_len(16), each = 8L)
  high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8L)
  x <- stats::rnorm(128)
  profile <- ifelse(stats::runif(128) <
                      stats::plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
  data.frame(school = school, wave = rep(seq_len(8), times = 16), x = x,
             y1 = stats::rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
             y2 = stats::rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7))
}

.surface_fit <- function(data, ...) {
  multilpa(data, c("y1", "y2"), "school", n_profiles = 2,
                 n_group_classes = 2, profile_covariates = "x",
                 n_starts = 4, seed = 1, ...)
}

test_that("coef reports every parameter rather than silently nothing", {
  data <- .surface_data()
  fit <- .surface_fit(data)
  estimates <- coef(fit)

  # Before this method existed, coef() fell through to coef.default() and
  # returned NULL, which looks like success.
  expect_type(estimates, "double")
  expect_length(estimates, fit$n_parameters)
  expect_false(is.null(estimates))
  expect_named(estimates)
  expect_equal(unname(estimates), parameter_inference(fit, data)$estimate)
  # means are in input units, variances in their own
  expect_equal(unname(estimates[grep("^measurement[.]mean[.]", names(estimates))]),
               as.vector(t(fit$means)))
  expect_true(all(estimates[grep("^measurement[.]variance[.]", names(estimates))] > 0))
  # A mean and a variance on the same indicator must not collide.
  expect_equal(anyDuplicated(names(estimates)), 0L)
  expect_equal(anyDuplicated(rownames(vcov(fit, data))), 0L)
})

test_that("confint matches the tidy interval and honours parm", {
  data <- .surface_data()
  fit <- .surface_fit(data)
  intervals <- confint(fit, data = data)
  inference <- parameter_inference(fit, data)

  expect_equal(dim(intervals), c(fit$n_parameters, 2L))
  expect_equal(unname(intervals[, 1L]), inference$conf_low)
  expect_equal(unname(intervals[, 2L]), inference$conf_high)
  expect_identical(rownames(intervals), names(coef(fit)))
  # a subset by name and by index agree
  chosen <- rownames(intervals)[c(1L, 9L)]
  expect_equal(confint(fit, parm = chosen, data = data),
               intervals[chosen, , drop = FALSE])
  expect_equal(confint(fit, parm = c(1L, 9L), data = data),
               intervals[c(1L, 9L), , drop = FALSE])
  expect_error(confint(fit, parm = "absent", data = data), "must identify")
  # `data` is optional: the fit stores the indicators, group index and designs
  # the information is rebuilt from, so the verb must not demand them back.
  expect_identical(confint(fit), confint(fit, data = data))
  # a wider level gives a wider interval
  wide <- confint(fit, data = data, level = 0.99)
  expect_true(all(wide[, 2L] - wide[, 1L] > intervals[, 2L] - intervals[, 1L]))
})

test_that("the shared diagnostics accept a covariate fit", {
  data <- .surface_data()
  fit <- .surface_fit(data)

  entropy <- get_results(fit, "entropy")
  expect_equal(nrow(entropy), 2L)
  expect_true(all(entropy$relative_entropy >= 0 & entropy$relative_entropy <= 1))

  criteria <- get_results(fit, "information_criteria", format = "long")
  expect_true(all(c("aic", "bic") %in% criteria$criterion))
  expect_equal(subset(criteria, criterion == "deviance")$value,
               -2 * fit$log_likelihood)

  classification <- get_results(fit, "classification", level = "both")
  expect_setequal(unique(classification$level), c("individuals", "groups"))
  expect_equal(sum(classification$n_modal[classification$level == "individuals"]),
               fit$n_observations)
  expect_equal(sum(classification$n_modal[classification$level == "groups"]),
               fit$n_groups)
})

test_that("a covariate fit can carry and report its ordering", {
  data <- .surface_data()
  bare <- .surface_fit(data)
  timed <- .surface_fit(data, time = "wave")

  expect_null(bare$time)
  expect_error(get_results(bare, "sequences"), class = "multilpa_no_time")
  expect_identical(timed$time, "wave")
  # the ordering is metadata and must not move an estimate
  expect_equal(bare$log_likelihood, timed$log_likelihood)
  expect_equal(bare$profile_coefficients, timed$profile_coefficients)

  long <- get_results(timed, "sequences")
  expect_named(long, c("group", "group_class", "time", "profile"))
  expect_equal(nrow(long), nrow(data))
  # The wide form carries the identifier as a column rather than as a row name,
  # so assert the contract by name instead of by width.
  expect_named(get_results(timed, "sequences", format = "wide"),
               c("group", "group_class", paste0("wave_", 1:8)))
  summary_table <- get_results(timed, "sequence_summary")
  expect_equal(sum(summary_table$groups), timed$n_groups)
  expect_equal(sum(summary_table$observations), nrow(data))
})

test_that("plot draws what a covariate model has and refuses what it has not", {
  data <- .surface_data()
  fit <- .surface_fit(data, time = "wave")
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)

  grDevices::png(file, width = 900, height = 650)
  expect_silent(plot(fit))
  expect_silent(plot(fit, what = "sequences"))
  grDevices::dev.off()
  expect_gt(file.size(file), 1000)

  # Prevalence is a function of each unit's covariates, so there is no single
  # vector to draw and the request is refused rather than averaged.
  expect_error(draw(plot(fit, what = "probabilities")),
               class = "multilpa_nothing_to_plot")
})

test_that("direct labels are measured in the weight they are drawn", {
  # strwidth() in the regular weight under-reserves the right margin, which
  # clipped the last character of every direct label. Text measurement needs an
  # open device: without one R opens Rplots.pdf and the widths come back from a
  # device nobody chose, so the measuring happens inside draw() like the drawing.
  draw({
    regular <- graphics::strwidth("Profile 1 (98%)", units = "inches", cex = 0.78)
    bold <- graphics::strwidth("Profile 1 (98%)", units = "inches", cex = 0.78,
                               font = 2L)
    expect_gt(bold, regular)
    expect_gt(.multilpa_label_margin("Profile 1 (98%)", 0.78),
              regular / graphics::par("csi") + 1.6)
  })
})
