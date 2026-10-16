# Review 2026-10-15: items settled by the lead after the parallel review.

# Every warning class a call raises, with the warnings muffled.
warning_classes <- function(expr) {
  classes <- character()
  withCallingHandlers(expr, warning = function(w) {
    classes <<- c(classes, class(w))
    invokeRestart("muffleWarning")
  })
  unique(classes)
}

test_that("a transition grid passes its candidates' warnings on, as the profile grid does", {
  skip_on_cran()
  vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  classes <- warning_classes(
    grid <- enumerate_classes(course_engagement, vars, "student", time = "sequence",
                              n_profiles = 2, n_group_classes = 1, n_starts = 1,
                              max_iter = 3, seed = 1))
  expect_true("latents_unconverged" %in% classes)
  expect_false(get_results(grid, "candidates")$converged)
})

test_that("a family grid passes its candidates' warnings on", {
  skip_on_cran()
  truth <- list(means = rbind(c(-1, -1), c(1, 1)),
                between = rbind(c(0.25, 0.4), c(0.4, 0.25)),
                within = c(1, 1), weights = c(0.5, 0.5))
  data <- additive_draw(truth, rep(8L, 40L), seed = 21)
  classes <- warning_classes(
    enumerate_classes(data, c("y1", "y2"), "group", family = "additive",
                      n_group_classes = 2, n_starts = 1, max_iter = 2, seed = 1))
  expect_true("latents_unconverged" %in% classes)
})

test_that("lta() converges by default where the tighter tolerance needs over 1000 updates", {
  skip_on_cran()
  # A three-profile model of two indicators takes about 1,500 updates at tol = 1e-10.
  fit <- lta(course_engagement, c("browse", "lectures"), "student", time = "sequence",
             n_profiles = 3, n_starts = 1, seed = 1)
  expect_true(fit$converged)
  expect_gt(fit$iterations, 1000)
  expect_identical(formals(lta)$max_iter, 5000L)
})

test_that("a model comparison marks no row and circles no point as best", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  vars <- c("browse", "lectures")
  two <- lpa(course_engagement, vars, n_profiles = 2, n_starts = 1, seed = 1)
  three <- lpa(course_engagement, vars, n_profiles = 3, n_starts = 1, seed = 1)
  comparison <- compare_models(two = two, three = three)
  expect_false(any(grepl("lowest", utils::capture.output(print(comparison)))))
  drawn <- plot(comparison)
  points <- vapply(drawn$layers, \(layer) inherits(layer$geom, "GeomPoint"), logical(1))
  expect_identical(sum(points), 1L)
})

test_that("noise and shared categories take the Okabe-Ito grey, which no class reuses", {
  expect_false(.gg_reserved_grey %in% .gg_class_palette)
  expect_true(all(c(.gg_class_palette, .gg_reserved_grey) %in% .gg_okabe_ito))
  # Seven profiles reach the grey's slot in the plain palette.
  noisy <- list(n_profiles = 7L, effective_profile_counts = rep(10, 7L),
                noise = TRUE, noise_posteriors = rep(0.5, 10L))
  expect_false(.gg_reserved_grey %in% .gg_profile_key(noisy)$colour)
  plain <- modifyList(noisy, list(noise = FALSE))
  expect_true(.gg_reserved_grey %in% .gg_profile_key(plain)$colour)
})
