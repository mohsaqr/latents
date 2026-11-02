# Review 2026-11-01, F5: confint() of a covariate fit holding measurement
# blocks by `fixed` named every natural parameter, held ones included, and
# assigned those names to the shorter interval matrix the inference table
# gives ("length of 'dimnames' [1] not equal to array extent"). The intervals
# are now named from the table's rows, default to the estimated coefficients,
# and a held coefficient is refused with `latents_held_parameter`, as the
# covariate-free confint.multilpa() refuses it. vcov() keeps the same rows.

skip_on_cran()

review_1101_case <- function() {
  set.seed(3)
  n <- 600L
  z <- stats::rnorm(n)
  state <- vapply(z, \(value) sample.int(3L, 1L, prob = exp(c(0.8, 0, -0.8) * value)),
                  integer(1))
  centres <- rbind(c(-2, 0), c(0, 2), c(2, -2))
  data <- data.frame(a = centres[state, 1L] + stats::rnorm(n),
                     b = centres[state, 2L] + stats::rnorm(n), z = z)
  plain <- lpa(data, c("a", "b"), n_profiles = 3, model = "VVI", n_starts = 3, seed = 1)
  held_fit <- function(fixed, model = "VVI") {
    lpa(data, c("a", "b"), n_profiles = 3, model = model, profile_covariates = "z",
        start = starting_values(plain), fixed = fixed, n_starts = 1, seed = 1)
  }
  list(measurement = held_fit("measurement"),
       means = held_fit("means"),
       variances = held_fit("variances"),
       free = lpa(data, c("a", "b"), n_profiles = 3, model = "VVI",
                  profile_covariates = "z", n_starts = 1, seed = 1))
}

review_1101_names <- function(table) {
  stem <- paste(table$level, table$parameter, table$outcome, sep = ".")
  ifelse(is.na(table$term), stem, paste(stem, table$term, sep = "."))
}

# The interval check a reader can redo by hand from the inference table.
expect_intervals_match_table <- function(fit) {
  table <- parameter_inference(fit)
  intervals <- confint(fit)
  expect_identical(rownames(intervals), review_1101_names(table))
  expect_identical(colnames(intervals), c("2.5%", "97.5%"))
  expect_equal(unname(intervals[, "2.5%"]), table$conf_low, tolerance = 1e-8)
  expect_equal(unname(intervals[, "97.5%"]), table$conf_high, tolerance = 1e-8)
  # Membership coefficients are unbounded: their interval is the plain Wald one.
  is_coefficient <- table$parameter == "coefficient"
  expect_gt(sum(is_coefficient), 0L)
  critical <- stats::qnorm(0.975)
  expect_equal(unname(intervals[is_coefficient, "2.5%"]),
               table$estimate[is_coefficient] - critical * table$standard_error[is_coefficient],
               tolerance = 1e-8)
  expect_equal(unname(intervals[is_coefficient, "97.5%"]),
               table$estimate[is_coefficient] + critical * table$standard_error[is_coefficient],
               tolerance = 1e-8)
  # vcov() covers the same rows on the natural scale, and its diagonal is the
  # table's squared standard errors.
  covariance <- vcov(fit)
  expect_identical(rownames(covariance), rownames(intervals))
  expect_identical(colnames(covariance), rownames(intervals))
  expect_equal(sqrt(unname(diag(covariance))), table$standard_error, tolerance = 1e-8)
  expect_identical(nrow(vcov(fit, scale = "unconstrained")), nrow(intervals))
  # coef() keeps every coefficient of the model, held ones included.
  expect_true(all(rownames(intervals) %in% names(coef(fit))))
  invisible(table)
}

cases <- review_1101_case()

test_that("default confint() of a fully held measurement model returns the membership rows", {
  fit <- cases$measurement
  expect_identical(fit$fixed, c("means", "variances"))
  table <- expect_intervals_match_table(fit)
  expect_identical(nrow(confint(fit)), 4L)
  expect_identical(unique(table$parameter), "coefficient")
  # The 16 model coefficients stay in coef(); 12 of them were held.
  expect_length(coef(fit), 16L)
  expect_identical(length(setdiff(names(coef(fit)), rownames(confint(fit)))), 12L)
})

test_that("partially held blocks keep the estimated measurement block's intervals", {
  means <- expect_intervals_match_table(cases$means)
  expect_identical(cases$means$fixed, "means")
  expect_setequal(unique(means$parameter), c("variance", "coefficient"))
  variances <- expect_intervals_match_table(cases$variances)
  expect_identical(cases$variances$fixed, "variances")
  expect_setequal(unique(variances$parameter), c("mean", "coefficient"))
})

test_that("named and numeric parm select estimated rows on every fit", {
  lapply(cases[c("measurement", "means", "variances")], \(fit) {
    intervals <- confint(fit)
    wanted <- grep("^profile\\.coefficient\\.", rownames(intervals), value = TRUE)
    expect_identical(confint(fit, parm = wanted), intervals[wanted, , drop = FALSE])
    positions <- match(wanted, names(coef(fit)))
    expect_identical(confint(fit, parm = positions), intervals[wanted, , drop = FALSE])
  })
  # A different level reaches the table through `level`, by name.
  wide <- confint(cases$measurement, parm = "profile.coefficient.profile_1.z", level = 0.9)
  table <- parameter_inference(cases$measurement, level = 0.9)
  expect_identical(colnames(wide), c("5%", "95%"))
  expect_equal(unname(wide[1L, ]),
               unlist(subset(table, outcome == "profile_1" & term == "z",
                             select = c(conf_low, conf_high)), use.names = FALSE),
               tolerance = 1e-8)
})

test_that("a held coefficient is refused by name or index with latents_held_parameter", {
  fit <- cases$measurement
  expect_error(confint(fit, parm = "measurement.mean.profile_1.a"),
               class = "latents_held_parameter")
  expect_error(confint(fit, parm = 1), class = "latents_held_parameter")
  expect_error(confint(fit, parm = c("profile.coefficient.profile_1.z",
                                     "measurement.variance.profile_2.b")),
               class = "latents_held_parameter")
  expect_error(confint(cases$means, parm = "measurement.mean.profile_3.b"),
               class = "latents_held_parameter")
  expect_error(confint(cases$variances, parm = "measurement.variance.profile_1.a"),
               class = "latents_held_parameter")
  # An estimated block of a partially held fit is not refused.
  expect_identical(nrow(confint(cases$means, parm = "measurement.variance.profile_1.a")), 1L)
  # Unknown names and indices outside coef() stay a bad argument.
  expect_error(confint(fit, parm = "measurement.mean.profile_9.a"),
               class = "latents_bad_argument")
  expect_error(confint(fit, parm = 17), class = "latents_bad_argument")
  expect_error(confint(fit, parm = Inf), class = "latents_bad_argument")
})

test_that("a fit holding nothing keeps every coefficient in confint(), vcov() and coef()", {
  fit <- cases$free
  expect_intervals_match_table(fit)
  expect_identical(rownames(confint(fit)), names(coef(fit)))
  expect_identical(rownames(vcov(fit)), names(coef(fit)))
  expect_identical(rownames(vcov(fit, scale = "unconstrained")),
                   names(coef(fit, scale = "unconstrained")))
})
