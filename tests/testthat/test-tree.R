# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

skip_if_not_installed("ggplot2")

# A two-level grid small enough to fit quickly: three numbers of profiles at
# one, two and three group classes.
tree_grid <- local({
  activity <- c("browse", "lectures", "forum_read")
  data <- subset(course_engagement, student <= 40)
  quietly(enumerate_classes(data, activity, "student", n_profiles = 1:3,
                            n_group_classes = 1:3, model = "VVI", n_starts = 2,
                            seed = 1))
})

test_that("the grid's trees are separate plots, one chosen by its count", {
  trees <- expect_plots(plot_tree(tree_grid))
  expect_named(trees, c("groups1", "groups2", "groups3"))
  one <- expect_plot(plot_tree(tree_grid, n_group_classes = 2))
  expect_match(one$labels$title, "2 group classes")
  expect_named(plot_tree(tree_grid, n_group_classes = c(1, 3)),
               c("groups1", "groups3"))
  # combine = TRUE is the one-page figure 0.9.21 drew, a panel per tree.
  combined <- expect_plot(plot_tree(tree_grid, combine = TRUE))
  expect_identical(nrow(ggplot2::ggplot_build(combined)$layout$layout), 3L)
})

test_that("plot(), plot_enumeration() and plot_tree() draw the same trees", {
  same_tree <- function(a, b) {
    expect_equal(plot_layer(a, 2L), plot_layer(b, 2L))
    expect_identical(a$labels$title, b$labels$title)
  }
  same_tree(plot(tree_grid, what = "tree", n_group_classes = 2),
            plot_tree(tree_grid, n_group_classes = 2))
  same_tree(plot_enumeration(tree_grid, what = "tree")[["groups3"]],
            plot_tree(tree_grid)[["groups3"]])
  same_tree(plot(tree_grid, what = "tree", combine = TRUE),
            plot_tree(tree_grid, combine = TRUE))
  expect_error(plot_enumeration(tree_grid, what = c("tree", "aic")),
               class = "latents_bad_argument")
  expect_error(plot(tree_grid, what = "tree", not_an_argument = 1),
               class = "latents_bad_argument")
})

test_that("a class keeps the colour of the class it matches in the first tree", {
  trees <- plot_tree(tree_grid)
  invisible(lapply(c("groups2", "groups3"), \(name) {
    count <- as.integer(sub("groups", "", name))
    invisible(lapply(2:3, \(k) {
      fit <- get_candidate(tree_grid, n_profiles = k, n_group_classes = count)
      reference <- get_candidate(tree_grid, n_profiles = k, n_group_classes = 1)
      matched <- .latents_match_classes(crossprod(fit$subject_posteriors,
                                                  reference$subject_posteriors))
      # The bars as drawn, with each class's number and its row.
      drawn <- trees[[name]]$layers[[2L]]$data
      reference_bars <- trees[["groups1"]]$layers[[2L]]$data
      # Rows sit at their count; the first tree starts at one profile, the
      # others at two, so its row k is the others' row k - 1.
      here <- drawn[drawn$row == k - 1L, ]
      there <- reference_bars[reference_bars$row == k, ]
      expect_identical(nrow(here), k)
      expect_identical(here$fill, there$fill[match(matched[here$class], there$class)])
    }))
  }))
})

test_that("anchoring gives a relabelled solution the first tree's colours", {
  # The same two classes in both trees, numbered the other way round in the
  # second: anchoring must swap its colours so each class keeps its own.
  posteriors <- cbind(c(0.9, 0.8, 0.1, 0.2), c(0.1, 0.2, 0.9, 0.8))
  tree <- function(p, colours) {
    list(rows = data.frame(key = c(1L, 2L)), weights = rep(1, 4),
         posteriors = list(matrix(1, 4, 1), p),
         colours = list("#000000", colours))
  }
  anchored <- .latents_tree_anchor(list(
    tree(posteriors, c("#E69F00", "#56B4E9")),
    tree(posteriors[, 2:1], c("#E69F00", "#56B4E9"))))
  expect_identical(anchored[[2L]]$colours[[2L]], c("#56B4E9", "#E69F00"))
  # A row whose count the first tree lacks keeps its own colours.
  unmatched <- tree(posteriors[, 2:1], c("#E69F00", "#56B4E9"))
  unmatched$rows$key <- c(1L, 5L)
  expect_identical(.latents_tree_anchor(list(
    tree(posteriors, c("#E69F00", "#56B4E9")), unmatched))[[2L]]$colours[[2L]],
    c("#E69F00", "#56B4E9"))
})

test_that("greedy matching takes the largest overlaps first", {
  overlap <- matrix(c(0.1, 0.6, 0.0,
                      0.5, 0.2, 0.1,
                      0.0, 0.1, 0.4), 3, byrow = TRUE)
  expect_identical(.latents_match_classes(overlap), c(2L, 1L, 3L))
  # More classes than reference classes: the extra one is left unmatched.
  expect_identical(.latents_match_classes(overlap[, 1:2]), c(2L, 1L, 0L))
})

test_that("the group-level tree divides groups into group classes", {
  trees <- expect_plots(plot_tree(tree_grid, level = "groups"))
  expect_named(trees, c("profiles2", "profiles3"))
  flows <- get_results(tree_grid, "tree", level = "groups", n_profiles = 3)
  expect_setequal(flows$from_candidate, c("profiles3_groups1", "profiles3_groups2"))
  fit <- get_candidate(tree_grid, n_profiles = 3, n_group_classes = 3)
  counts <- get_results(fit, "counts")
  expected <- counts$effective_proportion[counts$level == "groups"]
  into_three <- unique(flows[flows$to_candidate == "profiles3_groups3",
                             c("to_class", "to_share")])
  expect_equal(into_three$to_share[order(into_three$to_class)], expected,
               tolerance = 1e-10)
})

test_that("two chosen solutions are compared at either level", {
  chosen <- c("profiles3_groups1", "profiles3_groups3")
  expect_plot(plot_tree(tree_grid, candidates = chosen))
  expect_plot(plot_tree(tree_grid, candidates = chosen, level = "groups"))
  flows <- get_results(tree_grid, "tree", candidates = chosen)
  expect_identical(nrow(flows), 9L)
  # The shared mass is a joint distribution: it sums to one, and each class's
  # flows sum to its own share.
  expect_equal(sum(flows$mass), 1, tolerance = 1e-10)
  expect_equal(as.vector(tapply(flows$share_of_from, flows$from_class, sum)),
               rep(1, 3), tolerance = 1e-10)
  expect_equal(as.vector(tapply(flows$share_of_to, flows$to_class, sum)),
               rep(1, 3), tolerance = 1e-10)
  expect_error(plot_tree(tree_grid, candidates = chosen[c(1L, 1L)]),
               class = "latents_bad_argument")
})

test_that("one level is compared across numbers of group classes", {
  # A single number of profiles is one level of every tree: the rows are that
  # level's solutions, one per number of group classes.
  compared <- expect_plot(plot_tree(tree_grid, n_profiles = 3,
                                    n_group_classes = c(1, 3)))
  expect_match(compared$labels$title, "3 profiles compared")
  expect_identical(compared$labels$y, "Number of group classes")
  flows <- get_results(tree_grid, "tree", n_profiles = 3, n_group_classes = c(1, 3))
  expect_setequal(c(flows$from_candidate, flows$to_candidate),
                  c("profiles3_groups1", "profiles3_groups3"))
  expect_equal(sum(flows$mass), 1, tolerance = 1e-10)
  # Each class of the second row takes the colour of the first-row class it
  # shares the most cases with.
  bars <- compared$layers[[2L]]$data
  first <- get_candidate(tree_grid, n_profiles = 3, n_group_classes = 1)
  second <- get_candidate(tree_grid, n_profiles = 3, n_group_classes = 3)
  matched <- .latents_match_classes(crossprod(second$subject_posteriors,
                                              first$subject_posteriors))
  top <- bars[bars$row == 1L, ]
  bottom <- bars[bars$row == 2L, ]
  expect_identical(bottom$fill, top$fill[match(matched[bottom$class], top$class)])
  every <- expect_plot(plot_tree(tree_grid, n_profiles = 2))
  expect_identical(length(unique(every$layers[[2L]]$data$row)), 3L)
  groups <- expect_plot(plot_tree(tree_grid, level = "groups", n_group_classes = 2))
  expect_identical(groups$labels$y, "Number of profiles")
  # One solution is nothing to compare.
  expect_error(plot_tree(tree_grid, n_profiles = 3, n_group_classes = 2),
               class = "latents_no_tree")
})

test_that("from and to name two solutions by their class counts", {
  by_counts <- expect_plot(plot_tree(tree_grid, from = c(3, 1), to = c(3, 3)))
  by_names <- plot_tree(tree_grid, candidates = c("profiles3_groups1",
                                                  "profiles3_groups3"))
  expect_equal(by_counts$layers[[2L]]$data, by_names$layers[[2L]]$data)
  expect_identical(get_results(tree_grid, "tree", from = c(3, 1), to = c(3, 3)),
                   get_results(tree_grid, "tree", candidates = c("profiles3_groups1",
                                                                 "profiles3_groups3")))
  # One count is one group class; a name works as well.
  expect_plot(plot_tree(tree_grid, from = 2, to = "profiles3_groups2"))
  expect_plot(plot(tree_grid, what = "tree", from = c(2, 2), to = c(3, 2)))
  expect_error(plot_tree(tree_grid, from = c(3, 3)), class = "latents_bad_argument")
  expect_error(plot_tree(tree_grid, from = c(3, 3), to = c(3, 3)),
               class = "latents_bad_argument")
  expect_error(plot_tree(tree_grid, from = c(3, 9), to = c(3, 3)),
               class = "latents_unknown_candidate")
  # A repeated count is that count once.
  expect_equal(plot_tree(tree_grid, n_profiles = c(3, 3))$layers[[2L]]$data,
               plot_tree(tree_grid, n_profiles = 3)$layers[[2L]]$data)
})

test_that("tree arguments that name nothing are refused by class", {
  expect_error(plot_tree(tree_grid, n_group_classes = 9),
               class = "latents_unknown_candidate")
  expect_error(plot_tree(tree_grid, candidates = c("profiles2_groups1", "nope")),
               class = "latents_unknown_candidate")
  expect_error(plot_tree(tree_grid, candidates = "profiles2_groups1"),
               class = "latents_bad_argument")
  expect_error(plot_tree(tree_grid, candidates = c("profiles2_groups1",
                                                   "profiles3_groups1"),
                         n_profiles = 2),
               class = "latents_bad_argument")
  expect_error(plot_tree(tree_grid, candidates = c("profiles1_groups2",
                                                   "profiles2_groups2")),
               class = "latents_failed_candidate")
  single <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                  n_profiles = 1:2, n_starts = 1, seed = 1))
  expect_error(plot_tree(single, level = "groups"), class = "latents_no_tree")
  expect_error(plot_tree(single, n_profiles = 2), class = "latents_nothing_to_plot")
  # A grid with no tree leaves the table out of every table rather than failing.
  expect_false("tree" %in% names(get_results(
    quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                          n_profiles = 2, n_starts = 1, seed = 1)), "all")))
})

test_that("the shared mass is the weighted cross-product of the posteriors", {
  # Calibration against a hand computation: an implementation that ignored
  # which units the classes share (outer product of the shares) would keep
  # the margins and fail here.
  flows <- get_results(tree_grid, "tree", from = c(2, 1), to = c(3, 1))
  first <- get_candidate(tree_grid, n_profiles = 2, n_group_classes = 1)$subject_posteriors
  second <- get_candidate(tree_grid, n_profiles = 3, n_group_classes = 1)$subject_posteriors
  expected <- crossprod(first, second) / nrow(first)
  expect_equal(flows$mass, expected[cbind(flows$from_class, flows$to_class)],
               tolerance = 1e-12)
  independent <- outer(colSums(first), colSums(second)) / nrow(first)^2
  expect_gt(max(abs(expected - independent)), 0.05)
})

test_that("both levels use the counted units and sampling weights of the counts", {
  data <- subset(course_engagement, student <= 40)
  set.seed(4)
  weights <- data.frame(student = unique(data$student), sw = stats::runif(40, 0.2, 3))
  weighted <- quietly(enumerate_classes(merge(data, weights),
                                        c("browse", "lectures", "forum_read"),
                                        "student", n_profiles = 2,
                                        n_group_classes = 1:3, model = "VVI",
                                        weights = "sw", n_starts = 1, seed = 1))
  fit <- get_candidate(weighted, n_profiles = 2, n_group_classes = 3)
  counts <- get_results(fit, "counts")
  groups <- get_results(weighted, "tree", level = "groups", n_profiles = 2)
  drawn <- unique(groups[groups$to_candidate == "profiles2_groups3",
                         c("to_class", "to_share")])
  expect_equal(drawn$to_share[order(drawn$to_class)],
               counts$effective_proportion[counts$level == "groups"],
               tolerance = 1e-10)
  # Rows with every indicator missing are left out, as from the counts.
  empty <- data
  empty[1:30, c("browse", "lectures", "forum_read")] <- NA
  missing_grid <- quietly(enumerate_lpa(empty, c("browse", "lectures", "forum_read"),
                                        n_profiles = 1:2, model = "VVI",
                                        missing = "fiml", n_starts = 1, seed = 1))
  fit_two <- get_candidate(missing_grid, n_profiles = 2)
  counted <- fit_two$informative
  flows <- get_results(missing_grid, "tree")
  expected <- colSums(fit_two$subject_posteriors[counted, , drop = FALSE]) / sum(counted)
  into_two <- unique(flows[, c("to_class", "to_share")])
  expect_equal(into_two$to_share[order(into_two$to_class)], unname(expected),
               tolerance = 1e-10)
})

test_that("a noise component is refused, and every band with mass is drawn", {
  noisy <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                 n_profiles = 1:2, model = "VVI", noise = TRUE,
                                 n_starts = 1, seed = 1))
  expect_error(plot_tree(noisy), class = "latents_unsupported_noise")
  expect_false("tree" %in% names(get_results(noisy, "all")))
  tree <- plot_tree(tree_grid, n_group_classes = 1)
  flows <- get_results(tree_grid, "tree", n_group_classes = 1)
  expect_identical(length(unique(tree$layers[[1L]]$data$group)),
                   sum(flows$mass > 0))
})

test_that("classes in one row have distinct colours", {
  many <- quietly(enumerate_lpa(as.data.frame(scale(iris[1:4])), names(iris)[1:4],
                                n_profiles = 7:9, model = "VVI", n_starts = 1,
                                seed = 1))
  bars <- plot_tree(many)$layers[[2L]]$data
  expect_false(any(tapply(bars$fill, bars$row, anyDuplicated) > 0))
  expect_identical(.latents_distinct_colours(c("#E69F00", "#E69F00", "#56B4E9")),
                   c("#E69F00", "#009E73", "#56B4E9"))
})

test_that("matching maximises the shared mass and never pairs disjoint classes", {
  overlap <- matrix(c(0.4, 0.35, 0.25, 0), 2, byrow = TRUE)
  expect_identical(.latents_match_classes(overlap), c(2L, 1L))
  disjoint <- matrix(c(0.5, 0, 0, 0), 2, byrow = TRUE)
  expect_identical(.latents_match_classes(disjoint), c(1L, 0L))
})

test_that("the tree's arguments are checked, positionally and by class", {
  # Titles stay the second and third arguments, as before the tree had others.
  titled <- plot_tree(tree_grid, "My tree", "Sub", n_group_classes = 2)
  expect_identical(titled$labels$title, "My tree")
  expect_identical(titled$labels$subtitle, "Sub")
  expect_error(plot_tree(tree_grid, level = "wrong"), class = "latents_bad_argument")
  expect_error(plot_tree(tree_grid, n_profiles = "2"), class = "latents_bad_argument")
  expect_error(plot_tree(tree_grid, combine = NA), class = "latents_bad_argument")
  expect_error(plot(tree_grid, what = "tree", combine = NA),
               class = "latents_bad_argument")
  expect_plot(plot_enumeration(tree_grid, what = "tree", n_group_classes = 2))
  expect_plot(plot_enumeration(tree_grid, what = "tree", level = "groups",
                               n_profiles = 3))
  expect_error(plot_enumeration(tree_grid, what = "aic", n_group_classes = 2),
               class = "latents_bad_argument")
  unconverged <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                       n_profiles = 2:3, model = "VVI",
                                       max_iter = 1, n_starts = 1, seed = 1))
  expect_error(plot_tree(unconverged, candidates = c("profiles2", "profiles3")),
               class = "latents_unconverged")
  single <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                  n_profiles = 2:3, model = "VVI", n_starts = 1,
                                  seed = 1))
  expect_error(plot_tree(single, candidates = c("profiles2", "profiles3"),
                         level = "groups"), class = "latents_no_tree")
})

test_that("model resolves from and to on a grid crossed with structures", {
  crossed <- quietly(enumerate_lpa(engagement_small, c("browse", "lectures"),
                                   n_profiles = 2:3, model = c("EEI", "VVI"),
                                   n_starts = 1, seed = 1))
  expect_error(plot_tree(crossed, from = 2, to = 3),
               class = "latents_unknown_candidate")
  chosen <- plot_tree(crossed, from = 2, to = 3, model = "VVI")
  expect_equal(chosen$layers[[2L]]$data,
               plot_tree(crossed, candidates = c("profiles2_VVI",
                                                 "profiles3_VVI"))$layers[[2L]]$data)
})

test_that("a family grid's tree divides its groups, as the fits count them", {
  family <- quietly(enumerate_classes(subset(course_engagement, student <= 30),
                                      c("browse", "lectures"), "student",
                                      family = "additive", n_group_classes = 1:3,
                                      n_starts = 1, seed = 1))
  trees <- expect_plots(plot_tree(family))
  expect_named(trees, c("additive", "additive_equal"))
  flows <- get_results(family, "tree", model = "additive")
  expect_setequal(c(flows$from_candidate, flows$to_candidate),
                  c("additive_groups1", "additive_groups2", "additive_groups3"))
  fit <- .latents_tree_fit(family, match("additive_groups3",
                                         .latents_tree_table_of(family)$candidate))
  expected <- colSums(fit$group_posteriors) / nrow(fit$group_posteriors)
  into_three <- unique(flows[flows$to_candidate == "additive_groups3",
                             c("to_class", "to_share")])
  expect_equal(into_three$to_share[order(into_three$to_class)], unname(expected),
               tolerance = 1e-10)
  # One count names a number of group classes on a family grid.
  expect_identical(
    get_results(family, "tree", from = 2, to = 3, model = "additive"),
    get_results(family, "tree", candidates = c("additive_groups2", "additive_groups3")))
  expect_error(get_results(family, "tree", from = 2, to = 3),
               class = "latents_unknown_candidate")
})

test_that("a regression grid's tree divides the cases by their class posteriors", {
  regression <- quietly(enumerate_regressions(score ~ hours, study_hours,
                                              n_classes = 1:3, n_starts = 1,
                                              seed = 1))
  expect_plot(plot_tree(regression))
  flows <- get_results(regression, "tree")
  fit <- .latents_tree_fit(regression, 3L)
  posteriors <- .latents_regression_posteriors(fit, "individuals")
  into_three <- unique(flows[flows$to_candidate == "classes3", c("to_class", "to_share")])
  expect_equal(into_three$to_share[order(into_three$to_class)],
               unname(colMeans(posteriors)), tolerance = 1e-10)
  expect_equal(as.vector(tapply(flows$mass, flows$from_candidate, sum)), c(1, 1),
               tolerance = 1e-10)
  # The regression grid names its class counts n_classes; the tree takes it.
  flow_columns <- c("from_candidate", "from_class", "to_candidate", "to_class", "mass")
  expect_identical(get_results(regression, "tree", n_classes = 2:3)[flow_columns],
                   get_results(regression, "tree", from = 2, to = 3)[flow_columns])
  expect_identical(get_results(regression), get_results(regression, "fit"))
})
