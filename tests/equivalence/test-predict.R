# predict() for multilpa fits: training reproduction, new-group invariance,
# densities, refusals, and agreement with mclust's E-step.

activity <- c("browse", "lectures", "forum_read")
quiet_single <- function(expression) {
  withCallingHandlers(expression, latents_single_level = function(w) {
    invokeRestart("muffleMessage")
  })
}

test_that("predictions match mclust's E-step at the same parameters", {
  skip_without_reference("mclust")
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
  columns <- paste0("probability_profile_", 1:3)
  posterior <- unname(as.matrix(predict(fit, training, type = "posterior")[columns]))
  colnames(posterior) <- columns
  expected <- unname(reference$z)
  colnames(expected) <- columns
  # Absolute: the iris species are well separated, so posteriors run down to
  # 1e-12 and below. 1e-10, the original tolerance.
  expect_reference(posterior, expected, reference = "mclust",
                   covers = c("multilpa()", "predict.multilpa()"),
                   quantity = "posterior profile probabilities, VVV", tolerance = 1e-10,
                   measure = "absolute",
                   note = "Reference: mclust::estep(modelName = 'VVV') at the parameters latents estimated.")
  expect_reference(sum(predict(fit, training, type = "density")$log_density),
                   reference$loglik, reference = "mclust",
                   covers = c("multilpa()", "predict.multilpa()"),
                   quantity = "log likelihood of the training rows, VVV", tolerance = 1e-10,
                   measure = "relative")
})
