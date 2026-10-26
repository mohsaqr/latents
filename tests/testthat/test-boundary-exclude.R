# `boundary = "exclude"`: starts with a variance at `min_variance` are set
# aside. Two clusters and one distant observation, fitted with three profiles:
# the highest likelihood isolates that observation in a profile of its own
# with every variance at the bound, which "exclude" sets aside. The tests use
# the kmeans search, whose starts these data were chosen for.

boundary_data <- function() {
  set.seed(3)
  data.frame(y1 = c(stats::rnorm(60, 0), stats::rnorm(60, 4), 15),
             y2 = c(stats::rnorm(60, 0), stats::rnorm(60, 4), 15),
             x = stats::rnorm(121))
}

boundary_fit <- function(data, ..., search = "kmeans") {
  withCallingHandlers(
    lpa(data, c("y1", "y2"), n_profiles = 3, model = "VVI", n_starts = 5, seed = 1,
        search = search, ...),
    latents_boundary = function(w) invokeRestart("muffleWarning"),
    latents_thin_profiles = function(w) invokeRestart("muffleWarning"),
    latents_small_classes = function(w) invokeRestart("muffleWarning"))
}

test_that("the default keeps the singleton spike; exclude takes the best start without one", {
  data <- boundary_data()
  default <- boundary_fit(data)
  excluded <- boundary_fit(data, boundary = "exclude")
  starts <- get_results(excluded, "starts")
  expect_true(default$boundary)
  expect_identical(min(tabulate(default$subject_profiles, 3L)), 1L)
  expect_false(excluded$boundary)
  expect_equal(excluded$log_likelihood, max(subset(starts, !boundary)$log_likelihood),
               tolerance = 1e-8)
  expect_gt(default$log_likelihood, excluded$log_likelihood)
  # The starts are the same runs either way; only the reported one differs.
  expect_identical(get_results(default, "starts"), starts)
})

test_that("the default is unchanged by the new argument", {
  data <- boundary_data()
  expect_identical(boundary_fit(data)$parameters,
                   boundary_fit(data, boundary = "allow")$parameters)
})

test_that("exclude combines with select_start = \"converged\"", {
  data <- boundary_data()
  both <- boundary_fit(data, boundary = "exclude", select_start = "converged")
  starts <- get_results(both, "starts")
  expect_false(both$boundary)
  expect_equal(both$log_likelihood,
               max(subset(starts, !boundary & converged)$log_likelihood), tolerance = 1e-8)
})

test_that("exclude reaches enumerate_lpa() and the covariate model", {
  data <- boundary_data()
  enumeration <- withCallingHandlers(
    enumerate_lpa(data, c("y1", "y2"), n_profiles = 2:3, model = "VVI", n_starts = 5,
                  seed = 1, boundary = "exclude", search = "kmeans"),
    latents_boundary = function(w) invokeRestart("muffleWarning"),
    latents_thin_profiles = function(w) invokeRestart("muffleWarning"),
    latents_small_classes = function(w) invokeRestart("muffleWarning"))
  three <- boundary_fit(data, boundary = "exclude")
  candidate <- subset(get_results(enumeration, "candidates"), n_profiles == 3)
  expect_false(candidate$boundary)
  expect_equal(candidate$log_likelihood, three$log_likelihood, tolerance = 1e-8)
  covariate <- boundary_fit(data, boundary = "exclude", profile_covariates = "x")
  expect_false(covariate$boundary)
  expect_true(boundary_fit(data, profile_covariates = "x")$boundary)
})

test_that("categorical indicators alone are unaffected", {
  set.seed(4)
  items <- as.data.frame(matrix(stats::rbinom(300 * 4, 1, 0.5), 300))
  names(items) <- paste0("u", 1:4)
  default <- lca(items, names(items), n_classes = 2, n_starts = 4, seed = 1)
  excluded <- lca(items, names(items), n_classes = 2, n_starts = 4, seed = 1,
                    boundary = "exclude")
  expect_identical(default$parameters, excluded$parameters)
})

test_that("when every start is at the bound the best is kept, with a classed warning", {
  expect_warning(
    index <- .multilpa_select_start(c(-10, -5, -7), rep(TRUE, 3), "likelihood",
                                    interior = rep(FALSE, 3)),
    class = "latents_all_starts_boundary")
  expect_identical(index, 2L)
  expect_identical(.multilpa_select_start(c(-10, -5, -7), rep(TRUE, 3), "likelihood",
                                          interior = c(TRUE, FALSE, TRUE)), 3L)
  expect_error(.multilpa_select_start(c(-10, -5), c(TRUE, TRUE), "likelihood",
                                      interior = TRUE))
})

test_that("an unknown boundary value is refused", {
  expect_error(boundary_fit(boundary_data(), boundary = "nonsense"),
               "should be one of")
})
