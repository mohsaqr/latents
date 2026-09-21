# Centring the indicators before the measurement model sees them. The point of
# testing it here is the boundary it creates: the fit holds one scale and the
# caller holds another, and every verb that crosses that boundary has to know
# which side it is on.

activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

centred_fit <- function(centering = "person", n_profiles = 3L, ...) {
  multilpa(course_engagement, activity, "student", n_profiles = n_profiles,
           n_group_classes = 1, n_starts = 5, seed = 1, centering = centering,
           ...)
}

test_that("person centring subtracts each unit's own mean", {
  fit <- centred_fit()
  expect_identical(fit$centering, "person")
  expect_equal(dim(fit$centering_offsets), c(106L, 5L))
  # The offsets are the group means of the data that was handed in.
  expected <- vapply(activity, function(v) {
    unname(tapply(course_engagement[[v]], course_engagement$student, mean))
  }, numeric(106L))
  expect_equal(unname(fit$centering_offsets), unname(expected))
  # A person-centred indicator sums to zero within every unit, which is the
  # defining property: all the between-unit variation is gone.
  stored <- fit$indicator_data
  within <- vapply(seq_len(ncol(stored)), function(column) {
    max(abs(tapply(stored[, column], course_engagement$student, sum)))
  }, numeric(1))
  expect_true(all(within < 1e-10))
})

test_that("grand centring moves the origin and nothing else", {
  fit <- centred_fit("grand")
  expect_equal(nrow(fit$centering_offsets), 1L)
  expect_equal(unname(fit$centering_offsets[1L, ]),
               unname(colMeans(course_engagement[activity])))
  plain <- centred_fit("none")
  # Shifting every observation by a constant shifts the means and leaves the
  # spread, the likelihood and the classification where they were.
  expect_equal(fit$log_likelihood, plain$log_likelihood)
  expect_equal(fit$variances, plain$variances)
  expect_equal(fit$subject_posteriors, plain$subject_posteriors)
})

test_that("the default leaves the fit exactly as it was", {
  explicit <- centred_fit("none")
  implicit <- multilpa(course_engagement, activity, "student", n_profiles = 3,
                       n_group_classes = 1, n_starts = 5, seed = 1)
  expect_identical(explicit$centering, "none")
  expect_null(explicit$centering_offsets)
  expect_equal(explicit$log_likelihood, implicit$log_likelihood)
  expect_equal(explicit$means, implicit$means)
})

test_that("centring reproduces the transform done by hand", {
  by_hand <- course_engagement
  by_hand[activity] <- lapply(activity, function(v) {
    by_hand[[v]] - stats::ave(by_hand[[v]], by_hand$student, FUN = mean)
  })
  manual <- multilpa(by_hand, activity, "student", n_profiles = 3,
                     n_group_classes = 1, n_starts = 5, seed = 1)
  fit <- centred_fit()
  expect_equal(fit$log_likelihood, manual$log_likelihood)
  expect_equal(fit$means, manual$means)
  expect_equal(fit$variances, manual$variances)
})

test_that("the caller's own scale comes back from the fit", {
  fit <- centred_fit()
  # `get_data("data")` is what the model was fitted to, expressed as it was
  # handed in: the offsets are on the fit, so the round trip is exact.
  returned <- get_data(fit, "data")
  expect_equal(returned[activity], course_engagement[activity],
               ignore_attr = TRUE)
  expect_equal(returned$student, course_engagement$student)
})

test_that("alignment is still checked against the frame the caller has", {
  fit <- centred_fit()
  # The fit holds centred values and the caller holds raw ones; the check has
  # to compare like with like or it would reject every correct frame.
  joined <- get_data(fit, "assignments", data = course_engagement)
  expect_identical(nrow(joined), nrow(course_engagement))
  expect_true(all(c("profile", "browse") %in% names(joined)))
  # A reordered frame must still be caught.
  shuffled <- course_engagement[rev(seq_len(nrow(course_engagement))), ]
  expect_error(get_data(fit, "assignments", data = shuffled),
               class = "multilpa_bad_inference_data")
})

test_that("the diagnostics work on the scale the model was fitted on", {
  fit <- centred_fit()
  # Bivariate residuals compare an observed association with a model-implied
  # one, so handing them raw columns for a centred fit would compare two
  # different scales. Supplying the frame must match supplying nothing.
  from_fit <- get_data(fit, "residuals")
  supplied <- get_data(fit, "residuals", data = course_engagement)
  expect_equal(from_fit, supplied)
  # And the residuals must be the ones the by-hand transform gives.
  by_hand <- course_engagement
  by_hand[activity] <- lapply(activity, function(v) {
    by_hand[[v]] - stats::ave(by_hand[[v]], by_hand$student, FUN = mean)
  })
  manual <- multilpa(by_hand, activity, "student", n_profiles = 3,
                     n_group_classes = 1, n_starts = 5, seed = 1)
  expect_equal(from_fit$residual, get_data(manual, "residuals")$residual)
})

test_that("inference works on a centred fit and reports the centred scale", {
  fit <- centred_fit(n_profiles = 2)
  inference <- parameter_inference(fit)
  expect_s3_class(inference, "data.frame")
  # The means are deviations, so they straddle zero rather than sitting at the
  # raw indicator level.
  measurement <- subset(inference, level == "measurement" & parameter == "mean")
  expect_lt(min(measurement$estimate), 0)
  expect_gt(max(measurement$estimate), 0)
})

test_that("centring that leaves nothing to model is refused", {
  # One observation per unit: person centring removes all of it.
  single <- course_engagement[!duplicated(course_engagement$student), ]
  expect_error(
    multilpa(single, activity, "student", n_profiles = 2, n_group_classes = 1,
             n_starts = 1, seed = 1, centering = "person"),
    class = "multilpa_bad_data")
})

test_that("the fit says which scale it is on", {
  fit <- centred_fit()
  expect_identical(get_data(fit, "model")$centering, "person")
  expect_output(print(fit), "person-centred", fixed = TRUE)
  plain <- centred_fit("none")
  expect_identical(get_data(plain, "model")$centering, "none")
  expect_failure(expect_output(print(plain), "centred", fixed = TRUE))
})

test_that("centring and a constrained structure compose", {
  fit <- centred_fit(volume = "equal", shape = "varying", max_iter = 5000)
  expect_identical(fit$centering, "person")
  expect_identical(fit$covariance_structure, "EVI")
  expect_equal(nrow(get_data(fit, "profiles")), 3L * length(activity))
})
