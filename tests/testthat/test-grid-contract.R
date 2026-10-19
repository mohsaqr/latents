# One contract for every enumeration grid: every grid verb either draws a
# view the grid offers or refuses it by class, and plot_views(grid), the
# verbs and the registry behind them agree.

skip_if_not_installed("ggplot2")

# Small grids of the four kinds, built once and only when a test needs them.
contract_cache <- new.env(parent = emptyenv())
contract_grids <- function() {
  if (!is.null(contract_cache$grids)) return(contract_cache$grids)
  engagement <- subset(course_engagement, student <= 40)
  activity <- c("browse", "lectures")
  families_data <- additive_draw(
    list(means = rbind(c(-1, -1), c(1, 1)),
         between = rbind(c(0.25, 0.4), c(0.4, 0.25)),
         within = c(1, 1), weights = c(0.5, 0.5)),
    rep(8L, 40L), seed = 21)
  scores <- subset(growth_scores, student %in% unique(growth_scores$student)[1:40])
  contract_cache$grids <- list(
    profiles = quietly(enumerate_classes(
      engagement, activity, "student", n_profiles = 1:3, n_group_classes = 1:2,
      model = "VVI", n_starts = 1, seed = 1)),
    transitions = quietly(enumerate_classes(
      engagement, activity, "student", time = "sequence", n_profiles = 2:3,
      n_group_classes = 1:2, n_starts = 1, seed = 1)),
    families = quietly(enumerate_classes(
      families_data, c("y1", "y2"), "group", family = c("additive", "dispersion"),
      n_group_classes = 1:3, n_starts = 1, seed = 1)),
    regressions = quietly(enumerate_regressions(
      score ~ wave, scores, n_classes = 1:2, n_group_classes = 1:2,
      id = "student", n_starts = 1, seed = 1)))
  contract_cache$grids
}

# The outcome of one call: "plot", "plots", the class of a package
# condition, or "unclassed" for any other error.
call_outcome <- function(expr) {
  result <- tryCatch(suppressMessages(expr), error = identity)
  if (inherits(result, "ggplot")) {
    ggplot2::ggplot_build(result)
    return("plot")
  }
  if (inherits(result, "latents_plots")) {
    invisible(lapply(result, ggplot2::ggplot_build))
    return("plots")
  }
  if (inherits(result, "error")) {
    package <- grep("^(latents|multilpa)_", class(result), value = TRUE)
    return(if (length(package) > 0L) package[1L] else
      paste("unclassed:", conditionMessage(result)))
  }
  paste("returned", class(result)[1L])
}

test_that("every grid verb draws or refuses by class, on every kind of grid", {
  skip_on_cran()
  grids <- contract_grids()
  expect_identical(vapply(grids, \(grid) class(grid)[1L], character(1)),
                   c(profiles = "multilpa_enumeration",
                     transitions = "latents_transition_enumeration",
                     families = "latents_family_enumeration",
                     regressions = "latents_regression_enumeration"))
  views <- c("all", "entropy", "avepp", "diagnostics", "deltas", "heatmap",
             "smallest", "trust", "selection", "tree")
  outcomes <- do.call(rbind, lapply(names(grids), function(kind) {
    x <- grids[[kind]]
    calls <- c(
      list(plot = quote(plot(x)), plot_all = quote(plot(x, what = "all")),
           plot_enumeration = quote(plot_enumeration(x)),
           plot_selection = quote(plot_selection(x)),
           plot_tree = quote(plot_tree(x))),
      stats::setNames(lapply(views, \(view) bquote(plot_enumeration(x, what = .(view)))),
                      paste0("plot_enumeration_", views)),
      stats::setNames(lapply(views, \(view) bquote(plot(x, what = .(view)))),
                      paste0("plot_", views)))
    data.frame(kind = kind, call = names(calls),
               outcome = vapply(calls, \(call) call_outcome(eval(call)), character(1)),
               stringsAsFactors = FALSE)
  }))
  drawn <- outcomes$outcome %in% c("plot", "plots")
  refused <- grepl("^(latents|multilpa)_", outcomes$outcome)
  expect_identical(outcomes[!(drawn | refused), ], outcomes[integer(), ])
  # The only refusals of these grids are views their kind does not offer.
  expect_setequal(unique(outcomes$outcome[refused]), "latents_unsupported_view")
  # Every kind draws its default criteria, with plot() and plot_enumeration().
  defaults <- outcomes[outcomes$call %in% c("plot", "plot_all", "plot_enumeration"), ]
  expect_true(all(defaults$outcome == "plot"))
})

test_that("the registry, plot_views() and the verbs agree on every grid", {
  skip_on_cran()
  grids <- contract_grids()
  invisible(lapply(grids, function(x) {
    listed <- plot_views(x)
    expect_s3_class(listed, "data.frame")
    expect_named(listed, c("type", "group", "supported", "description", "reason"))
    expect_identical(anyDuplicated(listed$type), 0L)
    registry <- .latents_grid_views(x)
    expect_setequal(listed$type[listed$supported],
                    c(unname(registry$criteria), "all", registry$views))
    expect_setequal(listed$type[!listed$supported],
                    names(registry$refused) %||% character())
    expect_true(all(is.na(listed$reason[listed$supported])))
    expect_true(all(nzchar(listed$reason[!listed$supported])))
    invisible(lapply(listed$type[listed$supported], function(view) {
      expect_true(call_outcome(plot_enumeration(x, what = view)) %in% c("plot", "plots"),
                  label = sprintf("%s draws %s", class(x)[1L], view))
    }))
    invisible(lapply(listed$type[!listed$supported], function(view) {
      expect_error(plot_enumeration(x, what = view), class = "latents_unsupported_view")
      # The refusal names what the grid does offer.
      expect_error(plot_enumeration(x, what = view), "This grid offers \"")
    }))
  }))
})

test_that("each kind offers the views its candidates carry", {
  skip_on_cran()
  grids <- contract_grids()
  offered <- lapply(grids, \(x) plot_views(x)$type[plot_views(x)$supported])
  views <- c("entropy", "avepp", "diagnostics", "deltas", "heatmap", "smallest",
             "trust", "selection", "tree")
  expect_true(all(views %in% offered$profiles))
  expect_true(all(views %in% offered$transitions))
  expect_setequal(intersect(views, offered$families), c("entropy", "avepp", "diagnostics"))
  expect_setequal(setdiff(views, offered$regressions), "tree")
  expect_identical(setdiff(offered$transitions, c(views, "all")),
                   c("aic", "bic", "bic_individual"))
  expect_identical(setdiff(offered$regressions, c(views, "all")),
                   c("aic", "bic", "bic_rows", "sabic", "icl"))
  # plot_tree() and plot_selection() refuse what the grid does not offer.
  expect_error(plot_tree(grids$families), class = "latents_unsupported_view")
  expect_error(plot_tree(grids$regressions), class = "latents_unsupported_view")
  expect_error(plot_selection(grids$families), class = "latents_unsupported_view")
  expect_error(plot_selection(grids$families, what = "trust"),
               class = "latents_unsupported_view")
  # The refusal keeps the class callers caught before the contract.
  expect_error(plot_tree(grids$families), class = "latents_bad_argument")
  expect_error(plot_tree(grids$families), "\"entropy\"")
  # A transition grid draws its trees, at both levels.
  expect_plots(plot_tree(grids$transitions))
  expect_plots(plot_tree(grids$transitions, level = "groups"))
})

test_that("an object that is not a grid, or a criterion it lacks, is refused by class", {
  skip_on_cran()
  grids <- contract_grids()
  fit <- get_candidate(grids$transitions, n_profiles = 2, n_group_classes = 1)
  expect_error(plot_enumeration(fit), class = "latents_bad_argument")
  expect_error(plot_selection(fit), class = "latents_bad_argument")
  expect_error(plot_tree(fit), class = "latents_bad_argument")
  expect_error(plot_views(fit), class = "latents_bad_argument")
  expect_error(plot_enumeration(grids$transitions, what = "bic_groups"),
               class = "latents_unknown_criterion")
  expect_error(plot_enumeration(grids$regressions, what = "bic_groups"),
               class = "latents_bad_argument")
  expect_error(plot(grids$regressions, what = "bic_individual"),
               class = "latents_unknown_criterion")
  expect_error(plot_enumeration(grids$families, what = "deltas", using = "bic"),
               class = "latents_unsupported_view")
  expect_error(plot_enumeration(grids$transitions, what = "deltas", using = "bic_groups"),
               class = "latents_bad_argument")
  expect_error(plot(grids$regressions, what = "tree", n_profiles = 2),
               class = "latents_bad_argument")
})

test_that("transition and family plots keep their 0.9.21 calls", {
  skip_on_cran()
  grids <- contract_grids()
  invisible(lapply(grids[c("transitions", "families")], function(x) {
    table <- as.data.frame(x)
    # plot(x) draws the grid's summary, as plot_enumeration(x) does: every
    # default criterion as a panel of one plot (0.9.21 drew BIC alone).
    default <- expect_plot(plot(x))
    expect_equal(plot_layer(default, 1L), plot_layer(plot_enumeration(x), 1L))
    expect_identical(nrow(ggplot2::ggplot_build(default)$layout$layout), 3L)
    invisible(lapply(c("bic", "aic", "bic_individual"), function(criterion) {
      drawn <- expect_plot(plot(x, what = criterion))
      expect_equal(sort(plot_layer(drawn, 1L)$y),
                   sort(table[[criterion]][table$converged]))
      expect_identical(plot_layer(plot(x, criterion = criterion), 1L),
                       plot_layer(drawn, 1L))
      # 0.9.21 took the criterion as the second argument, by position.
      expect_identical(plot_layer(plot(x, criterion), 1L), plot_layer(drawn, 1L))
    }))
    # Unique abbreviations, as match.arg() took them.
    expect_identical(plot(x, what = "bic_i")$labels$y,
                     plot(x, what = "bic_individual")$labels$y)
    expect_identical(plot(x, criterion = "bic_i")$labels$y,
                     plot(x, what = "bic_individual")$labels$y)
    expect_identical(plot(x, what = "a")$labels$y, "AIC")
    # The 0.9.21 default vector, passed explicitly, is BIC.
    expect_identical(plot(x, criterion = c("bic", "aic", "bic_individual"))$labels$y,
                     "BIC")
    expect_identical(plot(x, main = "Mine")$labels$title, "Mine")
    expect_error(plot(x, what = "aic", criterion = "bic"), class = "latents_bad_argument")
    expect_error(plot(x, what = 3), class = "latents_bad_argument")
    # 0.9.21 ignored any other argument; it is still drawn, with a warning.
    expect_warning(drawn <- plot(x, colour = "red"), class = "latents_unused_argument")
    expect_plot(drawn)
    expect_identical(plot_layer(plot(x, what = "all", combine = FALSE)[[1L]], 1L)$y,
                     plot_layer(plot_enumeration(x, what = "aic"), 1L)$y)
  }))
  expect_identical(plot(grids$families)$labels$x, "Number of group classes")
  expect_identical(plot(grids$transitions)$labels$x, "Number of profiles")
})

test_that("a regression grid has a criteria plot and a tidy summary", {
  skip_on_cran()
  x <- contract_grids()$regressions
  table <- as.data.frame(x)
  expect_identical(summary(x), table)
  expect_s3_class(summary(x), "data.frame")
  expect_warning(summary(x, digits = 2), class = "latents_unused_argument")
  drawn <- expect_plot(plot(x))
  expect_identical(levels(plot_layer(drawn, 1L)$PANEL), c("1", "2", "3", "4"))
  expect_identical(drawn$labels$x, "Number of classes")
  expect_equal(sort(plot_layer(plot(x, what = "icl"), 1L)$y),
               sort(table$icl[table$converged %in% TRUE]))
  expect_error(plot(x, what = c("aic", "deltas")), class = "latents_bad_argument")
  expect_error(plot(x, colour = "red"), class = "latents_bad_argument")
})

test_that("the grid table and quality numbers are the grids' own", {
  skip_on_cran()
  grids <- contract_grids()
  regression <- .latents_grid_table(grids$regressions)
  stored <- as.data.frame(grids$regressions)
  expect_identical(regression$n_profiles, as.integer(stored$n_classes))
  expect_identical(regression$bic, stored$bic)
  expect_identical(anyDuplicated(regression$candidate), 0L)
  # The smallest-class view reads the table's smallest_share, for two or more
  # classes.
  shares <- .latents_grid_smallest_shares(grids$regressions)$individuals
  expect_equal(shares, ifelse(stored$n_classes >= 2, stored$smallest_share, NA))
  # The quality of a transition grid is its fits' posteriors summarised: the
  # individual entropy is the one the table stores, and each AvePP is a mean
  # posterior probability.
  quality <- .latents_grid_quality_table(grids$transitions)
  expect_identical(quality$profile_entropy, as.data.frame(grids$transitions)$entropy)
  fit <- get_candidate(grids$transitions, n_profiles = 3, n_group_classes = 2)
  posteriors <- fit$subject_posteriors
  modal <- max.col(posteriors, ties.method = "first")
  by_hand <- min(vapply(seq_len(ncol(posteriors)), \(k) mean(posteriors[modal == k, k]),
                        numeric(1)))
  at <- which(quality$n_profiles == 3 & quality$n_group_classes == 2)
  expect_equal(quality$avepp_min[at], by_hand)
  in_range <- c(quality$avepp_min, quality$group_avepp_min, quality$group_entropy)
  in_range <- in_range[!is.na(in_range)]
  expect_true(all(in_range >= 0 & in_range <= 1))
  # A family grid's quality is its group classes': one group class has none.
  families <- .latents_grid_quality_table(grids$families)
  expect_true(all(is.na(families$profile_entropy)))
  expect_true(all(is.na(families$group_avepp_min[families$n_group_classes == 1])))
  expect_true(all(is.finite(families$group_avepp_min[families$n_group_classes > 1])))
  # A criterion equal to one before it is one panel: bic_rows is bic when the
  # rows are the units, and differs when persons hold several rows.
  all_criteria <- unname(.latents_grid_views(grids$regressions)$criteria)
  two_level <- .latents_grid_criteria_titles(
    grids$regressions, .latents_grid_table(grids$regressions), all_criteria)
  expect_identical(names(two_level), c("AIC", "BIC", "BIC (rows)", "SABIC", "ICL"))
  single <- quietly(enumerate_regressions(score ~ hours, study_hours, n_classes = 1:2,
                                          n_starts = 1, seed = 1))
  single_table <- .latents_grid_table(single)
  one_level <- .latents_grid_criteria_titles(single, single_table, all_criteria)
  expect_identical(unname(one_level), c("aic", "bic", "sabic", "icl"))
  # Reordering the candidates does not change which panels are drawn.
  reversed <- single_table[rev(seq_len(nrow(single_table))), ]
  expect_identical(.latents_grid_criteria_titles(single, reversed, all_criteria),
                   one_level)
  expect_identical(levels(plot_layer(plot(single, what = "all"), 1L)$PANEL),
                   c("1", "2", "3", "4"))
})
