# predict() for multilpa fits: training reproduction, new-group invariance,
# densities, refusals, and agreement with mclust's E-step.

activity <- c("browse", "lectures", "forum_read")
quiet_single <- function(expression) {
  withCallingHandlers(expression, latents_single_level = function(w) {
    invokeRestart("muffleMessage")
  })
}

test_that("predicting the training rows reproduces the fit's posteriors", {
  gappy <- latents:::.mixture_with_seed(2, {
    d <- course_engagement
    d$browse[sample(nrow(d), 80)] <- NA
    d
  })
  activities <- c("time_with_friends", "on_social_media", "tv_video_games")
  esm <- subset(student_esm, day <= 3)
  fits <- list(
    two_level = multilpa(course_engagement, activity, "student", 2, 2,
                         n_starts = 2, seed = 1),
    full = multilpa(course_engagement, activity, "student", 2, 2,
                    covariance_model = "full", n_starts = 2, seed = 1),
    fiml = multilpa(gappy, activity, "student", 2, 2, missing = "fiml",
                    n_starts = 2, seed = 1),
    person = multilpa(course_engagement, activity, "student", 2, 2,
                      centering = "person", n_starts = 2, seed = 1),
    lca = multilca(esm, activities, "student", 2, 2, n_starts = 2, seed = 1),
    noise = quiet_single(multilpa(course_engagement, activity, NULL, 2,
                                  noise = TRUE, n_starts = 2, seed = 1)))
  data_for <- list(two_level = course_engagement, full = course_engagement,
                   fiml = gappy, person = course_engagement, lca = esm,
                   noise = course_engagement)
  invisible(lapply(names(fits), function(label) {
    fit <- fits[[label]]
    predicted <- predict(fit, data_for[[label]], type = "posterior")
    columns <- grep("^probability_", names(predicted), value = TRUE)
    expected <- cbind(fit$subject_posteriors,
                      if (isTRUE(fit$noise)) fit$noise_posteriors)
    expect_equal(unname(as.matrix(predicted[columns])), unname(expected),
                 tolerance = 1e-12, info = label)
    # Without `newdata` the rows are rebuilt from the fit; re-centring them
    # can differ from the original in the last bit.
    expect_equal(predict(fit), predict(fit, data_for[[label]]),
                 tolerance = 1e-12, info = label)
  }))
})

test_that("a single-level fit's densities sum to its log likelihood", {
  fit <- quiet_single(multilpa(course_engagement, activity, NULL, 3,
                               n_starts = 2, seed = 1))
  densities <- predict(fit, course_engagement, type = "density")
  expect_equal(sum(densities$log_density), fit$log_likelihood,
               tolerance = 1e-10)
})

test_that("a new group is classified from its own rows alone", {
  fit <- multilpa(subset(course_engagement, student <= 60), activity,
                  "student", 2, 2, n_starts = 2, seed = 1)
  newcomers <- subset(course_engagement, student > 100)
  together <- predict(fit, newcomers, type = "posterior")
  alone <- predict(fit, subset(newcomers, student == 101), type = "posterior")
  expect_equal(subset(together, student == 101, select = -row),
               subset(alone, select = -row), ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_true(all(c("row", "student", "profile", "posterior", "group_class",
                    "group_posterior") %in% names(together)))
})

test_that("predictions match mclust's E-step at the same parameters", {
  skip_if_not_installed("mclust")
  training <- stats::setNames(iris[1:4], c("sl", "sw", "pl", "pw"))
  fit <- quiet_single(multilpa(training, names(training), NULL, 3,
                               covariance_model = "full", n_starts = 2,
                               seed = 1))
  parameters <- list(
    pro = drop(fit$profile_probabilities), mean = t(unname(fit$means)),
    variance = list(modelName = "VVV", d = 4L, G = 3L,
                    sigma = unname(fit$covariances),
                    cholsigma = array(vapply(1:3, function(k) {
                      chol(fit$covariances[, , k])
                    }, numeric(16)), c(4L, 4L, 3L))))
  # `estep()` resolves its per-model worker from the calling frame.
  reference <- do.call(mclust::estep,
                       list(data = as.matrix(training), modelName = "VVV",
                            parameters = parameters),
                       envir = asNamespace("mclust"))
  predicted <- predict(fit, training, type = "posterior")
  expect_equal(unname(as.matrix(predicted[paste0("probability_profile_", 1:3)])),
               unname(reference$z), tolerance = 1e-10)
  expect_equal(sum(predict(fit, training, type = "density")$log_density),
               reference$loglik, tolerance = 1e-10)
})

test_that("new rows the fit cannot handle are refused by class", {
  fit <- multilpa(course_engagement, activity, "student", 2, 2, n_starts = 1,
                  seed = 1)
  expect_error(predict(fit, subset(course_engagement, select = -browse)),
               class = "latents_bad_data")
  expect_error(predict(fit, subset(course_engagement, select = -student)),
               class = "latents_bad_data")
  gappy <- course_engagement
  gappy$browse[1] <- NA
  expect_error(predict(fit, gappy), class = "latents_bad_data")
  activities <- c("time_with_friends", "on_social_media", "tv_video_games")
  esm <- subset(student_esm, day <= 3)
  lca <- multilca(esm, activities, "student", 2, 2, n_starts = 1, seed = 1)
  strange <- esm
  strange$time_with_friends <- as.character(strange$time_with_friends)
  strange$time_with_friends[1] <- "sometimes"
  expect_error(predict(lca, strange), class = "latents_bad_data")
})
