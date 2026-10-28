# Model-labelled controls must keep every candidate when dots are forwarded.
test_that("comparison wrappers preserve restored model labels and BIC weights", {
  withr::local_seed(2812)
  data <- data.frame(y = rnorm(100L))
  fits <- list(
    lpa(data, "y", n_profiles = 1L, n_starts = 1L, seed = 1L),
    lpa(data, "y", n_profiles = 2L, n_starts = 2L, seed = 1L),
    lpa(data, "y", n_profiles = 2L, model = "EEI", n_starts = 2L, seed = 1L))
  forward <- function(...) compare_models(...)
  nested <- function(...) forward(...)
  append_model <- list(
    n_profiles = function(...) compare_models(..., n_profiles = fits[[3L]]),
    n_group_classes = function(...) compare_models(..., n_group_classes = fits[[3L]]),
    cutoff = function(...) compare_models(..., cutoff = fits[[3L]]))
  invisible(lapply(names(append_model), function(label) {
    # Put the restored model between ordinary arguments: retain caller order
    # as well as all three rows, their attached fits and normalized weights.
    arguments <- stats::setNames(fits, c("first", label, "third"))
    reference <- do.call(compare_models, arguments)
    expect_identical(reference$model, names(arguments))
    invisible(lapply(list(forward, nested), function(wrapper) {
      actual <- do.call(wrapper, arguments)
      expect_identical(actual, reference)
    }))
    appended <- append_model[[label]](first = fits[[1L]], second = fits[[2L]])
    expected <- do.call(compare_models,
                         stats::setNames(fits, c("first", "second", label)))
    expect_identical(appended, expected)
  }))
  several <- list(n_profiles = fits[[1L]], second = fits[[2L]], cutoff = fits[[3L]])
  expect_identical(do.call(nested, several), do.call(compare_models, several))
})
