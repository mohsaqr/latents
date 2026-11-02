# Wild-card verbs: one function per picture or table, any object, no `what`.
# Each view must be the one the existing method draws (equivalence), the
# views drawn from posteriors must carry the numbers of the diagnostic tables
# (calibration), and a view an object lacks must be refused by class.

skip_if_not_installed("ggplot2")

wild_lpa <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
                model = "EEI", n_starts = 2, seed = 1)
wild_two_level <- quietly(multilpa(
  subset(course_engagement, student <= 40), c("browse", "lectures", "forum_read"),
  "student", n_profiles = 2, n_group_classes = 2, n_starts = 2, seed = 1))
wild_regression <- mixture_regression(score ~ hours, data = study_hours,
                                      n_classes = 2, n_starts = 2, seed = 1)
wild_grid <- quietly(enumerate_lpa(faithful, c("eruptions", "waiting"),
                                   n_profiles = 1:3, model = "EEI",
                                   n_starts = 2, seed = 1))

layers_of <- function(plot) lapply(seq_along(plot$layers), \(i) plot_layer(plot, i))

test_that("a profile fit's views are the ones plot_diagnostics() and plot() draw", {
  draw({
    expect_equal(layers_of(plot_avepp(wild_lpa)),
                 layers_of(plot_diagnostics(wild_lpa, what = "avepp")))
    expect_equal(layers_of(plot_entropy(wild_lpa)),
                 layers_of(plot_diagnostics(wild_lpa, what = "entropy")))
    expect_equal(layers_of(plot_bars(wild_lpa)), layers_of(plot(wild_lpa, what = "bars")))
    expect_equal(layers_of(plot_pairs(wild_lpa)), layers_of(plot(wild_lpa, what = "pairs")))
    lapply(list(plot_sizes(wild_lpa), plot_posteriors(wild_lpa), plot_heatmap(wild_lpa),
                plot_raincloud(wild_lpa), plot_parallel(wild_lpa)), expect_plot)
  })
})

test_that("a fit with group classes draws both levels of a classification view", {
  draw({
    both <- plot_avepp(wild_two_level)
    expect_plots(both)
    expect_named(both, c("individuals", "groups"))
    expect_equal(layers_of(both$groups),
                 layers_of(plot_diagnostics(wild_two_level, what = "group_avepp")))
  })
})

test_that("a grid's views are the ones plot_enumeration() draws", {
  draw({
    expect_equal(layers_of(plot_criteria(wild_grid)),
                 layers_of(plot_enumeration(wild_grid, what = "all")))
    expect_equal(layers_of(plot_avepp(wild_grid)),
                 layers_of(plot_enumeration(wild_grid, what = "avepp")))
    expect_equal(layers_of(plot_sizes(wild_grid)),
                 layers_of(plot_enumeration(wild_grid, what = "smallest")))
    lapply(list(plot_entropy(wild_grid), plot_heatmap(wild_grid),
                plot_deltas(wild_grid), plot_trust(wild_grid)), expect_plot)
  })
})

test_that("views drawn from posteriors carry the diagnostic tables' numbers", {
  draw({
    # AvePP tiles: the "average_posteriors" table, cell by cell.
    tiles <- plot_layer(plot_avepp(wild_regression), 1L)
    table <- get_avepp(wild_regression)
    cells <- plot_avepp(wild_regression)$data
    expect_equal(sort(cells$average_posterior), sort(table$average_posterior),
                 tolerance = 1e-12)
    # Each assigned row of the matrix sums to one.
    expect_equal(as.vector(rowsum(cells$average_posterior, cells$row)), c(1, 1),
                 tolerance = 1e-12)
    expect_equal(nrow(tiles), 4L)

    # Sizes: the "classification" table's effective counts.
    sizes <- plot_sizes(wild_regression)$data
    classes <- get_classification(wild_regression)
    expect_equal(sizes$count[order(sizes$class)], unname(classes$estimated_n),
                 tolerance = 1e-10)
    expect_equal(sum(sizes$share), 1, tolerance = 1e-12)

    # Entropy: the subtitle states the "entropy" table's relative entropy.
    entropy <- plot_entropy(wild_regression)
    expect_match(entropy$labels$subtitle,
                 sprintf("%.3f", get_entropy(wild_regression)$relative_entropy),
                 fixed = TRUE)
    expect_true(all(entropy$data$value >= 0 & entropy$data$value <= 1))
  })
})

test_that("each fitted family gets every classification view", {
  growth <- mixture_regression(score ~ wave, subset(growth_scores, student <= 60),
                               n_classes = 2, id = "student", class_level = "group",
                               random = "intercept", n_starts = 2, seed = 1)
  draw({
    lapply(list(wild_regression, growth), \(fit) {
      lapply(list(plot_sizes(fit), plot_posteriors(fit), plot_entropy(fit),
                  plot_avepp(fit)), expect_plot)
    })
    expect_plot(plot_trajectories(growth))
    expect_plot(plot_individuals(growth))
    expect_plot(plot_classification(growth))
    expect_plot(plot_coefficients(wild_regression))
    expect_plot(plot_trajectories(wild_regression))
  })
})

test_that("plot_all() draws every view of fits, grids and other results", {
  draw({
    profile_views <- expect_plots(plot_all(wild_lpa))
    expect_false(anyDuplicated(names(profile_views)) > 0L)
    expect_true(all(c("profiles", "avepp", "entropy") %in% names(profile_views)))
    expect_true(all(c("group_avepp", "group_sizes") %in% names(plot_all(wild_two_level))))
    regression_views <- expect_plots(suppressMessages(plot_all(wild_regression)))
    expect_true(all(c("sizes", "entropy", "avepp", "posteriors") %in%
                      names(regression_views)))
    grid_views <- expect_plots(plot_all(wild_grid))
    expect_true(all(c("criteria", "entropy", "avepp", "trust") %in% names(grid_views)))
    compared <- compare_models(two = wild_lpa,
                               three = get_candidate(wild_grid, "profiles3"))
    expect_plots(plot_all(compared))
    expect_plot(plot_criteria(compared))
    # A diagnostics result is drawn from its fit.
    expect_equal(names(plot_all(diagnostics(wild_lpa, plots = FALSE))),
                 names(profile_views))
  })
})

test_that("tables are the diagnostics() tables, and a grid stacks its candidates", {
  quality <- diagnostics(wild_lpa, plots = FALSE)
  expect_identical(get_checks(wild_lpa), get_results(quality, "checks"))
  expect_identical(get_entropy(wild_lpa), get_results(quality, "entropy"))
  expect_identical(get_avepp(wild_lpa), get_results(quality, "average_posteriors"))
  expect_identical(get_classification(wild_lpa), get_results(quality, "classification"))
  expect_identical(get_residuals(wild_lpa), get_results(quality, "residuals"))
  expect_identical(get_replication(wild_lpa), get_results(quality, "replication"))
  expect_identical(get_entropy(quality), get_results(quality, "entropy"))

  stacked <- get_entropy(wild_grid)
  expect_identical(names(stacked)[1L], "candidate")
  expect_setequal(unique(stacked$candidate), c("profiles1", "profiles2", "profiles3"))
  one <- subset(stacked, candidate == "profiles2", select = -candidate)
  row.names(one) <- NULL
  expect_equal(one, get_entropy(get_candidate(wild_grid, "profiles2")))
  regression_grid <- quietly(enumerate_regressions(score ~ hours, study_hours,
                                                   n_classes = 1:2, n_starts = 2,
                                                   seed = 1))
  expect_setequal(unique(get_checks(regression_grid)$candidate),
                  c("classes1_groups1", "classes2_groups1"))
})

test_that("a view or table an object lacks is refused by class", {
  draw({
    expect_error(plot_bars(wild_regression), class = "latents_nothing_to_plot")
    expect_error(plot_trajectories(wild_lpa), class = "latents_nothing_to_plot")
    expect_error(plot_posteriors(wild_grid), class = "latents_nothing_to_plot")
    transitions <- quietly(lta(subset(course_engagement, student <= 30),
                               c("browse", "lectures"), "student", n_profiles = 2,
                               time = "sequence", n_starts = 1, seed = 1))
    # The transition plot() has no raincloud; refused by class, not match.arg().
    expect_error(plot_raincloud(transitions), class = "latents_nothing_to_plot")
    expect_error(plot_avepp(data.frame(a = 1)), class = "latents_not_a_fit")
    expect_error(plot_all(list()), class = "latents_not_a_fit")
  })
  expect_error(get_residuals(wild_regression), class = "latents_bad_argument")
  expect_error(get_checks(data.frame(a = 1)), class = "latents_not_a_fit")
  expect_error(overview(data.frame(a = 1)), class = "latents_not_a_fit")
  expect_error(overview(wild_lpa, plots = "yes"), class = "latents_bad_argument")
})

test_that("overview() prints any object and returns it invisibly", {
  draw({
    printed <- capture.output(returned <- withVisible(overview(wild_regression)))
    expect_false(returned$visible)
    expect_identical(returned$value, wild_regression)
    expect_true(any(grepl("Diagnostics", printed, fixed = TRUE)))
    expect_gt(length(capture.output(overview(wild_grid, plots = FALSE))), 5L)
    expect_gt(length(capture.output(overview(wild_lpa, plots = FALSE))), 5L)
  })
})

# Review 2026-11-02 (validation/review-2026-11-02-wildcards/REPORT.md).

test_that("a grid's failed candidates are left out of its tables (finding 1)", {
  d <- data.frame(y = rep(c(-1, 1), 15))
  grid <- suppressWarnings(enumerate_regressions(y ~ 1, d, n_classes = c(1, 40),
                                                 n_starts = 1, seed = 1))
  expect_s3_class(grid$fits[[2L]], "condition")
  checks <- get_checks(grid)
  expect_identical(unique(checks$candidate), "classes1_groups1")
  expect_equal(checks[names(checks) != "candidate"],
               get_checks(grid$fits[[1L]]), ignore_attr = "row.names")
})

test_that("posteriors of other families are each eligible unit's assigned class (findings 2, 3)", {
  set.seed(5)
  d <- data.frame(g = rep(1:20, each = 5), t = rep(0:4, 20))
  d$y <- c(-2, 2)[1 + d$g %% 2] + 0.3 * d$t + stats::rnorm(100)
  d$w <- rep(c(0, 1, 3, 1), length.out = 20)[d$g]
  weighted <- suppressWarnings(mixture_regression(
    y ~ t, d, n_classes = 2, id = "g", n_group_classes = 2, weights = "w",
    n_starts = 1, seed = 1, vcov_type = "none"))
  trajectories <- mixture_regression(score ~ wave, subset(growth_scores, student <= 30),
                                     n_classes = 2, id = "student",
                                     class_level = "group", n_starts = 1, seed = 1)
  draw({
    both <- plot_posteriors(weighted)
    expect_named(both, c("individuals", "groups"))
    units <- get_entropy(weighted)
    expect_equal(c(nrow(both$individuals$data), nrow(both$groups$data)),
                 units$n_units[match(c("individuals", "groups"), units$level)])
    # The weighted mean certainty drawn for each class (the segment layer) is
    # the table's average posterior of that class.
    classes <- subset(get_classification(weighted), level == "individuals")
    means <- plot_layer(both$individuals, 3L)
    expect_equal(sort(unique(means$x)), sort(classes$average_posterior),
                 tolerance = 1e-10)
    expect_true(all(both$individuals$data$value >= 0.5 - 1e-12))
    expect_true(all(c("posteriors_individuals", "posteriors_groups") %in%
                      names(suppressMessages(plot_all(weighted)))))
    certainty <- plot_posteriors(trajectories)$data
    expect_equal(nrow(certainty), 30L)
    expect_true(all(certainty$value >= 0.5 - 1e-12))
    expect_false(identical(certainty, plot_classification(trajectories)$data))
  })
})

test_that("one class still has a one-by-one AvePP (finding 4)", {
  one <- mixture_regression(score ~ hours, study_hours, n_classes = 1,
                            n_starts = 1, seed = 1, vcov_type = "none")
  draw({
    cells <- plot_avepp(one)$data
    expect_equal(nrow(cells), 1L)
    expect_equal(cells$average_posterior, get_avepp(one)$average_posterior)
    expect_true("avepp" %in% names(suppressMessages(plot_all(one))))
    expect_error(plot_entropy(one), class = "latents_nothing_to_plot")
  })
})

test_that("plot_all() keeps an object's own refusal (finding 5)", {
  boot <- quietly(parameter_inference(wild_lpa, method = "bootstrap", iter = 3,
                                      n_starts = 1, seed = 1))
  reordered <- boot[rev(seq_len(nrow(boot))), , drop = FALSE]
  direct <- tryCatch(plot(reordered), error = identity)
  draw(wildcard <- tryCatch(plot_all(reordered), error = identity))
  expect_s3_class(wildcard, "latents_bad_argument")
  expect_false(inherits(wildcard, "latents_nothing_to_plot"))
  expect_identical(conditionMessage(wildcard), conditionMessage(direct))
})

test_that("a refusal recommends only what this object can draw (finding 6)", {
  transitions <- quietly(lta(subset(course_engagement, student <= 10),
                             c("browse", "lectures"), "student", n_profiles = 2,
                             time = "sequence", n_starts = 1, seed = 1))
  refusal <- tryCatch(plot_raincloud(transitions), error = identity)
  expect_s3_class(refusal, "latents_nothing_to_plot")
  advice <- sub(".*For this object: ", "", conditionMessage(refusal))
  expect_false(grepl("plot_raincloud|plot_parallel|plot_pairs", advice))
  expect_match(advice, "plot_avepp()", fixed = TRUE)
})

test_that("no refusal of a wild-card verb points to a `what =` argument", {
  categorical <- quietly(multilca(subset(student_esm, day <= 1),
                                  c("time_with_friends", "on_social_media",
                                    "tv_video_games"), "student", n_profiles = 2,
                                  n_group_classes = 1, n_starts = 1, seed = 1))
  noisy <- quietly(lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
                       noise = TRUE, n_starts = 2, seed = 1))
  objects <- list(wild_lpa, categorical, noisy, wild_regression, wild_grid)
  verbs <- c("plot_sizes", "plot_entropy", "plot_posteriors", "plot_avepp",
             "plot_heatmap", "plot_bars", "plot_raincloud", "plot_parallel",
             "plot_pairs", "plot_trajectories", "plot_individuals",
             "plot_classification", "plot_coefficients", "plot_criteria",
             "plot_deltas", "plot_trust", "get_checks", "get_entropy",
             "get_classification", "get_avepp", "get_residuals", "get_replication")
  draw({
    messages <- unlist(lapply(objects, \(object) vapply(verbs, \(verb) tryCatch({
      get(verb)(object)
      ""
    }, error = conditionMessage), character(1))))
  })
  refused <- messages[nzchar(messages)]
  expect_gt(length(refused), 20L)
  expect_false(any(grepl("what *=", refused)))
  expect_true(all(grepl("For this object: ", refused, fixed = TRUE)))
  # The categorical fit is pointed to the views it has.
  raincloud <- tryCatch(plot_raincloud(categorical), error = conditionMessage)
  expect_match(raincloud, "plot_heatmap()", fixed = TRUE)
  expect_match(raincloud, "plot_responses()", fixed = TRUE)
  expect_error(plot_avepp(noisy), class = "latents_nothing_to_plot")
  expect_error(get_residuals(noisy), class = "latents_unsupported_noise")
})

test_that("grid tables have their own verbs, and the summary warns about edge minima", {
  edge <- quietly(enumerate_lpa(faithful, c("eruptions", "waiting"), n_profiles = 1:2,
                                model = "EEI", n_starts = 3, seed = 1))
  expect_identical(get_best(edge), get_results(edge, "best"))
  # get_candidates() is the candidates table without its all-empty columns.
  full <- get_results(edge, "candidates")
  kept <- get_candidates(edge)
  expect_identical(kept, full[names(kept)])
  expect_true(all(vapply(full[setdiff(names(full), names(kept))],
                         \(column) all(is.na(column) | column %in% ""), logical(1))))
  expect_false(any(vapply(kept, \(column) all(is.na(column)), logical(1))))
  expect_true(all(c("candidate", "n_profiles", "bic_groups") %in% names(kept)))
  expect_identical(get_knees(edge), get_results(edge, "knees"))
  printed <- capture.output(print(summary(edge)))
  expect_true(any(grepl("edge of the grid", printed, fixed = TRUE)))
  expect_false(any(grepl("what *=", printed)))
  expect_lt(length(printed), 40L)
  # Not every criterion at the edge: no edge warning.
  expect_false(any(grepl("edge of the grid", capture.output(print(summary(wild_grid))),
                         fixed = TRUE)))
  expect_error(get_best(wild_lpa), class = "latents_not_a_fit")
})

# Review 2026-11-02, latest (validation/review-2026-11-02-latest/REPORT.md).

test_that("grid edges are those of the requested grid, failures included (finding 1)", {
  grid <- data.frame(candidate = paste0("profiles", 1:3), n_profiles = 1:3,
                     n_group_classes = 1L, converged = c(TRUE, TRUE, FALSE),
                     n_best_replicated = c(2L, 2L, NA_integer_),
                     avepp_below = c(0L, 0L, NA_integer_), group_avepp_below = NA_integer_)
  inside <- capture.output(.multilpa_print_enumeration_warnings(
    grid, data.frame(candidate = "profiles2")))
  expect_false(any(grepl("edge of the grid", inside, fixed = TRUE)))
  grid$converged <- c(FALSE, FALSE, TRUE)
  grid$n_best_replicated <- c(NA, NA, 2L)
  at_edge <- capture.output(.multilpa_print_enumeration_warnings(
    grid, data.frame(candidate = "profiles3")))
  expect_true(any(grepl("edge of the grid", at_edge, fixed = TRUE)))
})

test_that("grid getters read tables without building the others (findings 2, 3)", {
  d <- data.frame(y = rep(c(-1, 1), 15))
  one <- suppressWarnings(enumerate_regressions(y ~ 1, d, n_classes = 1, n_starts = 1,
                                                seed = 1))
  expect_s3_class(get_candidates(one), "data.frame")
  families <- quietly(enumerate_classes(
    subset(course_engagement, student <= 30), c("browse", "lectures"), "student",
    family = "additive", n_group_classes = 1:3, n_starts = 1, seed = 1))
  expect_identical(get_tree(families), get_results(families, "tree"))
  expect_error(get_knees(families), class = "latents_bad_argument")
})

test_that("the grid summary respects an explicit `rows` (finding 4)", {
  summarised <- summary(wild_grid)
  all_rows <- capture.output(print(summarised))
  one_row <- capture.output(print(summarised, rows = 1))
  headings <- capture.output(print(summarised, rows = 0))
  criteria <- "^ +(aic|kic|bic|sabic|caic|awe|icl|clc) "
  expect_equal(sum(grepl(criteria, all_rows)), nrow(get_best(wild_grid)))
  expect_equal(sum(grepl(criteria, one_row)), 1L)
  expect_equal(sum(grepl(criteria, headings)), 0L)
})
