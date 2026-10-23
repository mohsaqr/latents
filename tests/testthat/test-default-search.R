# The default start search (`search = NULL`): emEM when every profile shares
# one covariance matrix, kmeans otherwise. `search = "kmeans"` reproducing
# the fits of 0.96.2 is checked against the recorded panel by
# equivalence/search/; here, the resolution and its propagation.

default_search_items <- c("cognitive_strategies", "intrinsic_value", "self_efficacy")
default_search_quietly <- function(expr) suppressWarnings(expr)

test_that("the default search follows the covariance structure", {
  expect_identical(.latents_default_search("EEI", "y"), "emEM")
  expect_identical(.latents_default_search("EEE", "y"), "emEM")
  expect_identical(.latents_default_search("EII", "y"), "emEM")
  invisible(lapply(c("VVI", "VVV", "VII", "VEI", "EVI", "VVE", "EEV"), function(structure) {
    expect_identical(.latents_default_search(structure, "y"), "kmeans")
  }))
  # Without a continuous indicator there is no covariance to share.
  expect_identical(.latents_default_search("EEI", character()), "kmeans")
})

test_that("fits record the search their structure resolves to", {
  fit_with <- function(model) default_search_quietly(
    lpa(srl, default_search_items, 2, model = model, n_starts = 3, seed = 1))
  expect_identical(fit_with("EEI")$search, "emEM")
  expect_identical(fit_with("EEE")$search, "emEM")
  expect_identical(fit_with("VVI")$search, "kmeans")
  expect_identical(fit_with("VVV")$search, "kmeans")
  categorical <- data.frame(a = rep(1:2, 60), b = rep(c(1, 2, 2, 3), 30))
  expect_identical(default_search_quietly(lca(categorical, c("a", "b"), 2, n_starts = 3,
                                              seed = 1))$search, "kmeans")
})

test_that("the default is exactly the named search it resolves to", {
  strip <- function(fit) { fit$call <- NULL; fit }
  shared <- default_search_quietly(lpa(srl, default_search_items, 2, model = "EEI",
                                       n_starts = 3, seed = 1))
  named <- default_search_quietly(lpa(srl, default_search_items, 2, model = "EEI",
                                      n_starts = 3, seed = 1, search = "emEM"))
  expect_identical(strip(shared), strip(named))
  separate <- default_search_quietly(lpa(srl, default_search_items, 2, model = "VVI",
                                         n_starts = 3, seed = 1))
  kmeans <- default_search_quietly(lpa(srl, default_search_items, 2, model = "VVI",
                                       n_starts = 3, seed = 1, search = "kmeans"))
  expect_identical(strip(separate), strip(kmeans))
})

test_that("a named search overrides the structure's default", {
  fit <- default_search_quietly(lpa(srl, default_search_items, 2, model = "EEI",
                                    n_starts = 3, seed = 1, search = "kmeans"))
  expect_identical(fit$search, "kmeans")
  expect_null(get_results(fit, "starts")$stage)
})

test_that("an invalid search is refused by class", {
  expect_error(lpa(srl, default_search_items, 2, search = "fast"),
               class = "latents_bad_argument")
  expect_error(lpa(srl, default_search_items, 2, search = c("kmeans", "emEM")),
               class = "latents_bad_argument")
  expect_error(lpa(srl, default_search_items, 2, search = NA_character_),
               class = "latents_bad_argument")
})

test_that("both stages of a staged fit resolve the default as direct fits do", {
  staged <- default_search_quietly(fit_staged(
    course_engagement, c("browse", "lectures", "forum_read"), id = "student",
    n_profiles = 2, n_group_classes = 2, variance_model = "equal", n_starts = 3, seed = 1))
  expect_identical(staged$stage_one$search, "emEM")
  expect_identical(staged$search, "emEM")
  named <- default_search_quietly(fit_staged(
    course_engagement, c("browse", "lectures", "forum_read"), id = "student",
    n_profiles = 2, n_group_classes = 2, variance_model = "equal", n_starts = 3, seed = 1,
    search = "kmeans"))
  expect_identical(named$stage_one$search, "kmeans")
  expect_identical(named$search, "kmeans")
})

test_that("refits replay the resolved search", {
  fit <- default_search_quietly(lpa(srl, default_search_items, 2, model = "EEI",
                                    n_starts = 3, seed = 1))
  expect_identical(.latents_refit_controls(fit)$search, "emEM")
})
