named_plot_fits <- function() {
  set.seed(11)
  group <- rep(seq_len(24L), each = 8L)
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = if (j %% 2L == 0L) c(0.85, 0.15) else c(0.2, 0.8))
  }, integer(1))
  continuous <- data.frame(
    g = group,
    a = stats::rnorm(length(group), c(-2, 2)[profile], 1),
    b = stats::rnorm(length(group), c(40, 60)[profile], 6)
  )
  list(
    lpa = multilpa(continuous, c("a", "b"), "g", 2, 2, n_starts = 1, seed = 5),
    lca = lca(student_esm, c("time_with_friends", "on_social_media", "sports"),
              n_classes = 2, n_starts = 1, seed = 1, tol = 1e-10),
    lta = lta(subset(course_engagement, student <= 40),
              c("browse", "lectures"), "student", n_profiles = 2,
              time = "sequence", n_starts = 1, seed = 1),
    grid = enumerate_lpa(continuous, c("a", "b"), n_profiles = 1:3,
                         n_starts = 1, seed = 3)
  )
}

same_drawing <- function(named, generic) {
  expect_equal(named$data, generic$data)
  expect_identical(named$labels$title, generic$labels$title)
}

test_that("each named plot draws what plot(what = ) draws", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fits <- named_plot_fits()
  same_drawing(expect_plot(plot_profiles(fits$lpa)), plot(fits$lpa))
  same_drawing(plot_profiles(fits$lpa, scale = "standardized", labels = FALSE),
               plot(fits$lpa, scale = "standardized", labels = FALSE))
  same_drawing(expect_plot(plot_responses(fits$lca, category = "first")),
               plot(fits$lca, what = "responses", category = "first"))
  same_drawing(expect_plot(plot_prevalence(fits$lpa)),
               plot(fits$lpa, what = "probabilities"))
  same_drawing(expect_plot(plot_diagnostics(fits$lpa, what = "avepp")),
               plot(fits$lpa, what = "avepp"))
  same_drawing(expect_plot(plot_transitions(fits$lta)),
               plot(fits$lta, what = "transitions"))
  # A transition fit takes the measurement views too, without intervals.
  same_drawing(plot_profiles(fits$lta), plot(fits$lta, what = "profiles"))
  same_drawing(expect_plot(plot_tree(fits$grid)),
               plot(fits$grid, what = "tree"))
  # Titles are the caller's to set through the named function as well.
  titled <- plot_profiles(fits$lpa, main = "custom", subtitle = "sub")
  expect_identical(titled$labels$title, "custom")
  expect_identical(titled$labels$subtitle, "sub")
})

test_that("plot_diagnostics() returns every supported view, as diagnostics does", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fits <- named_plot_fits()
  every <- plot_diagnostics(fits$lpa)
  expect_s3_class(every, "latents_plots")
  expect_named(every, c("sizes", "posteriors", "entropy", "avepp"))
  # plot() on a diagnostics object delegates here, so the two cannot differ.
  through_tables <- plot(diagnostics(fits$lpa, plots = FALSE))
  expect_equal(lapply(every, `[[`, "data"), lapply(through_tables, `[[`, "data"))
  # One profile leaves nothing for posteriors or entropy to separate.
  set.seed(2)
  one <- lpa(data.frame(a = stats::rnorm(60)), "a", n_profiles = 1,
             n_starts = 1, seed = 1)
  expect_named(plot_diagnostics(one), c("sizes", "avepp"))
  expect_error(plot_diagnostics(one, what = "entropy"),
               class = "latents_nothing_to_plot")
})

test_that("named plots refuse what a fit cannot draw, by condition class", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fits <- named_plot_fits()
  # An argument the fit's own plot() lacks is named, not dropped.
  expect_error(plot_profiles(fits$lta, statistic = "median"),
               class = "latents_bad_argument")
  expect_error(plot_profiles(fits$lta, intervals = FALSE),
               class = "latents_bad_argument")
  expect_error(plot_responses(fits$lpa), class = "latents_no_categorical")
  expect_error(plot_prevalence(fits$lta), class = "latents_nothing_to_plot")
  expect_error(plot_transitions(fits$lpa), class = "latents_bad_argument")
  expect_error(plot_profiles(fits$lpa, statistic = "average"), "should be one of")
  # Their arguments are their own formals, so a misspelt one is refused by R.
  expect_error(plot_profiles(fits$lpa, palette = "x"), "unused argument")
  expect_error(plot_tree(fits$lpa), "multilpa_enumeration")
})

test_that("enumeration plots choose criteria with `what`, and 0.9.21's `criterion` still works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fits <- named_plot_fits()
  grid <- fits$grid
  same_drawing(plot(grid, what = "bic_groups"),
               plot_enumeration(grid, what = "bic_groups"))
  same_drawing(plot(grid, what = "all"), plot_enumeration(grid, what = "all"))
  # The view's former name still draws the default criteria.
  same_drawing(plot(grid, what = "enumeration"), plot(grid))
  same_drawing(plot(grid, what = "tree"), plot_tree(grid))
  expect_error(plot(grid, what = c("tree", "aic")), class = "latents_bad_argument")
  # 0.9.21 spellings, named and positional, draw what they drew then.
  same_drawing(plot(grid, criterion = "aic"), plot(grid, what = "aic"))
  same_drawing(plot(grid, what = "enumeration", criterion = "aic"),
               plot(grid, what = "aic"))
  same_drawing(plot(grid, "enumeration", c("aic", "bic_groups")),
               plot(grid, what = c("aic", "bic_groups")))
  expect_warning(drawn <- plot(grid, what = "tree", criterion = "aic"),
                 class = "latents_unused_argument")
  same_drawing(drawn, plot_tree(grid))
  expect_error(plot(grid, what = "bic_groups", criterion = "aic"),
               class = "latents_bad_argument")
  # plot_enumeration() is newer than 0.9.21 and has only `what`.
  expect_error(plot_enumeration(grid, criterion = "aic"), "unused argument")

  set.seed(1)
  ratings <- data.frame(
    team = rep(seq_len(40), each = 6),
    climate = rep(stats::rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
      stats::rnorm(240))
  families <- quietly(enumerate_classes(ratings, "climate", "team",
                                        family = c("additive", "dispersion"),
                                        n_group_classes = 1:2, n_starts = 1,
                                        seed = 1))
  transitions <- quietly(enumerate_classes(
    subset(course_engagement, student <= 40), c("browse", "lectures"),
    "student", time = "sequence", n_profiles = 2:3, n_group_classes = 1,
    n_starts = 1, seed = 1))
  for_each_grid <- list(families, transitions)
  invisible(lapply(for_each_grid, function(other) {
    drawn <- expect_plot(plot(other, what = "aic"))
    expect_identical(drawn$labels$y, "AIC")
    # 0.9.21 named this argument `criterion`.
    same_drawing(plot(other, criterion = "aic"), plot(other, what = "aic"))
    expect_error(plot(other, what = "bic", criterion = "aic"),
                 class = "latents_bad_argument")
    # 0.9.21 ignored an argument it did not use; it now says so and draws.
    expect_warning(plot(other, colour = "red"),
                   class = "latents_unused_argument")
  }))
})
