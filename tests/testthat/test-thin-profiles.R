# A profile with fewer effective members than its own means and variances is
# reported by a `latents_thin_profiles` warning; its estimate is not changed.

thin_profile_data <- function(scale = 1) {
  set.seed(3)
  truth <- sample.int(3L, 360L, replace = TRUE, prob = c(0.45, 0.35, 0.2))
  data.frame(a = scale * stats::rnorm(360L, c(-1.2, 1.2, 0)[truth]),
             b = scale * stats::rnorm(360L, c(-1.2, 1.2, 1.92)[truth]))
}

test_that("a profile's own parameters are its means and what its spread adds", {
  expect_identical(.latents_profile_parameters(5L, "varying", "full"), 20)
  expect_identical(.latents_profile_parameters(4L, "varying", "diagonal"), 8L)
  expect_identical(.latents_profile_parameters(5L, "equal", "full"), 0)
  expect_identical(.latents_profile_parameters(0L, "varying", "diagonal"), 0)
  # A constrained structure counts the parameters one more profile adds.
  expect_equal(.latents_profile_parameters(5L, "varying", "full", "VVV"), 20)
  expect_equal(.latents_profile_parameters(4L, "varying", "diagonal", "VVI"), 8)
  expect_equal(.latents_profile_parameters(3L, "varying", "full", "VEV"), 3 + 1 + 3)
  expect_equal(.latents_profile_parameters(3L, "varying", "full", "EEE"), 0)
})

test_that("the warning names the thin profiles and stays silent otherwise", {
  expect_warning(thin <- .latents_warn_thin_profiles(c(40, 3.9, 12), 4),
                 class = "latents_thin_profiles")
  expect_identical(thin, 2L)
  expect_no_warning(.latents_warn_thin_profiles(c(40, 3.9), 0))
  expect_no_warning(.latents_warn_thin_profiles(3.9, 4))
  expect_no_warning(.latents_warn_thin_profiles(c(40, 30), 4))
  message <- tryCatch(.latents_warn_thin_profiles(c(5.1, 40, 3.2), 8),
                      latents_thin_profiles = conditionMessage)
  expect_true(startsWith(message, "Profiles 1, 3: 5.1, 3.2 members for 8"))
})

# The emEM search reaches the spurious profile of these data (3.9 members);
# the kmeans search, from fewer starts, does not.
test_that("a spurious profile is reported, whatever the scale", {
  fits <- lapply(c(1, 10), function(scale) {
    data <- thin_profile_data(scale)
    expect_warning(fit <- lpa(data, c("a", "b"), 3, n_starts = 4, seed = 1, tol = 1e-12,
                              max_iter = 10000, search = "emEM"),
                   class = "latents_thin_profiles")
    fit
  })
  counts <- lapply(fits, function(fit) sort(unname(fit$effective_profile_counts)))
  expect_lt(counts[[1]][1], 4)
  expect_equal(counts[[2]], counts[[1]], tolerance = 1e-6)
})

test_that("a well-determined fit raises no thin-profile warning", {
  expect_no_warning(quietly(lpa(srl, names(srl), 3, seed = 1)),
                    class = "latents_thin_profiles")
})

test_that("a VVV fit with a spurious profile is reported", {
  expect_warning(quietly(lpa(srl, names(srl), 4, model = "VVV", seed = 1, search = "emEM"),
                         classes = "latents_boundary"),
                 class = "latents_thin_profiles")
})
