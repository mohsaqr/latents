# predict() for multilpa fits: training reproduction, new-group invariance,
# densities, refusals, and agreement with mclust's E-step.

activity <- c("browse", "lectures", "forum_read")
quiet_single <- function(expression) {
  withCallingHandlers(expression, latents_single_level = function(w) {
    invokeRestart("muffleMessage")
  })
}

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
