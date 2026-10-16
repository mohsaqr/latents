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

test_that("an occasion-varying lta() plot orders panels by occasion number and wraps labels", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  set.seed(4)
  panel_data <- data.frame(id = rep(seq_len(40), each = 11), time = rep(seq_len(11), 40))
  panel_data$y <- stats::rnorm(nrow(panel_data), ifelse(stats::runif(nrow(panel_data)) < 0.5, -1.5, 1.5))
  fit <- quietly(lta(panel_data, "y", "id", n_profiles = 2, time = "time",
                     transitions = "occasion", n_starts = 1, seed = 1, max_iter = 20),
                 c("latents_unconverged", "latents_boundary", "latents_small_classes"))
  panels <- levels(ggplot2::ggplot_build(plot(fit))$plot$data$panel)
  expect_identical(panels, paste0("group_class_1\noccasion ", 2:11))
})

test_that("a legend of more than four keys wraps into rows instead of running off the edge", {
  skip_if_not_installed("ggplot2")
  key <- data.frame(label = factor(paste("Class", 1:5), levels = paste("Class", 1:5)),
                    colour = .gg_class_palette[1:5], shape = 15:19)
  wrapped <- .gg_scales(key, legend = TRUE)
  expect_identical(wrapped$colour$guide$params$nrow, 2)
  four <- key[1:4, ]
  four$label <- droplevels(four$label)
  expect_identical(.gg_scales(four, legend = TRUE)$colour$guide, "legend")
  expect_identical(.gg_scales(key, legend = FALSE)$colour$guide, "none")
})

test_that("a grid over group classes at one profile count puts the classes on the x-axis", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  grid <- quietly(
    enumerate_classes(course_engagement, vars, "student", n_profiles = 3,
                      n_group_classes = 2:4, n_starts = 1, seed = 1),
    c("latents_boundary", "latents_small_classes", "latents_unconverged"))
  built <- ggplot2::ggplot_build(plot_enumeration(grid))
  points <- built$data[[2L]]
  expect_setequal(unique(points$x), 2:4)
  expect_match(built$plot$labels$x, "group classes")
  # A grid over profiles keeps profiles on the x-axis.
  profiles <- quietly(
    enumerate_classes(course_engagement, vars, "student", n_profiles = 2:3,
                      n_group_classes = 2, n_starts = 1, seed = 1),
    c("latents_boundary", "latents_small_classes", "latents_unconverged"))
  expect_identical(ggplot2::ggplot_build(plot_enumeration(profiles))$plot$labels$x,
                   "Number of profiles")
})
