# The replicates of a bootstrap table: as.data.frame(x, what = "replicates")
# and plot(), which draws the estimates with their percentile intervals. The
# table's standard errors and percentile intervals are read off the same
# replicates, so the tidy replicates must reproduce them.

skip_on_cran()

bootstrap_case <- local({
  made <- NULL
  function() {
    if (!is.null(made)) return(made)
    set.seed(31)
    n <- 240L
    truth <- sample.int(3L, n, replace = TRUE)
    data <- data.frame(a = stats::rnorm(n, c(0, 2.5, 5)[truth]),
                       b = stats::rnorm(n, c(0, 3, 1)[truth], c(1, 0.7, 1.2)[truth]))
    fit <- lpa(data, c("a", "b"), n_profiles = 3, model = "VVI", n_starts = 3, seed = 1)
    made <<- list(fit = fit, data = data,
                  boot = parameter_inference(fit, method = "bootstrap", iter = 25,
                                             n_starts = 2, seed = 1))
    made
  }
})

test_that("a bootstrap table keeps its rows and gains a class", {
  boot <- bootstrap_case()$boot
  expect_s3_class(boot, c("latents_bootstrap", "data.frame"), exact = TRUE)
  table <- as.data.frame(boot)
  expect_identical(class(table), "data.frame")
  expect_identical(names(table), names(boot))
  expect_null(attr(table, "replicates"))
  # vcov() still reads the replicates' covariance.
  expect_identical(dim(attr(boot, "covariance")), rep(nrow(boot), 2L))
})

test_that("the tidy replicates reproduce the standard errors and intervals", {
  boot <- bootstrap_case()$boot
  replicates <- as.data.frame(boot, what = "replicates")
  expect_identical(names(replicates),
                   c("replicate", "level", "outcome", "term", "parameter", "estimate"))
  # One row per kept replicate and parameter.
  expect_identical(nrow(replicates), attr(boot, "n_valid") * nrow(boot))
  key <- with(replicates, paste(level, outcome, term, parameter))
  expect_false(anyDuplicated(paste(replicates$replicate, key)) > 0L)
  # Each parameter's replicates give its standard error and its interval.
  rows <- with(boot, paste(level, outcome, term, parameter))
  spread <- vapply(rows, function(row) stats::sd(replicates$estimate[key == row]), numeric(1))
  expect_equal(unname(spread), boot$standard_error, tolerance = 1e-12)
  limits <- vapply(rows, function(row) {
    stats::quantile(replicates$estimate[key == row], c(0.025, 0.975), names = FALSE)
  }, numeric(2))
  expect_equal(limits[1L, ], boot$conf_low, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(limits[2L, ], boot$conf_high, tolerance = 1e-12, ignore_attr = TRUE)
  # Probabilities of the profiles sum to one in every replicate.
  shares <- subset(replicates, parameter == "probability" & level == "profile")
  expect_equal(as.vector(tapply(shares$estimate, shares$replicate, sum)),
               rep(1, attr(boot, "n_valid")), tolerance = 1e-10)
})

test_that("a subset of the table no longer has replicates and says so", {
  boot <- bootstrap_case()$boot
  expect_error(as.data.frame(subset(boot, parameter == "mean"), what = "replicates"),
               class = "latents_bad_argument")
  expect_error(as.data.frame(boot, what = "draws"), class = "latents_bad_argument")
  # Reordering or repeating every row keeps the attributes; the replicates
  # must not be paired with the moved rows.
  reversed <- boot[rev(seq_len(nrow(boot))), , drop = FALSE]
  expect_error(as.data.frame(reversed, what = "replicates"), class = "latents_bad_argument")
  expect_error(plot(reversed), class = "latents_bad_argument")
  repeated <- boot[c(1L, seq_len(nrow(boot) - 1L)), , drop = FALSE]
  expect_error(as.data.frame(repeated, what = "replicates"), class = "latents_bad_argument")
})

test_that("as.data.frame() keeps the row names, or takes the ones given", {
  boot <- bootstrap_case()$boot
  expect_identical(row.names(as.data.frame(boot)), row.names(boot))
  given <- sprintf("p%02d", seq_len(nrow(boot)))
  expect_identical(row.names(as.data.frame(boot, row.names = given)), given)
})

test_that("plot() draws each estimate with its percentile interval", {
  skip_if_not_installed("ggplot2")
  boot <- bootstrap_case()$boot
  table <- subset(as.data.frame(boot), parameter == "mean")
  means <- plot(boot, what = "mean")
  expect_s3_class(means, "ggplot")
  # One line per profile across the indicators; the points and whiskers are
  # the table's estimates and percentile bounds, in its order.
  frame <- means$data
  expect_identical(levels(frame$panel), "mean")
  expect_identical(frame$series, table$outcome)
  expect_identical(as.character(frame$at), table$term)
  expect_equal(frame$estimate, table$estimate, tolerance = 0)
  expect_equal(frame$conf_low, table$conf_low, tolerance = 0)
  expect_equal(frame$conf_high, table$conf_high, tolerance = 0)
  expect_match(means$labels$subtitle, "95% percentile intervals from 25 refitted", fixed = TRUE)
  # The single group class's probability is one in every replicate: it is
  # named in the caption, not drawn.
  shares <- plot(boot, what = "probability")
  expect_identical(nrow(shares$data), 3L)
  expect_match(shares$labels$caption, "group_class_1", fixed = TRUE)
  everything <- plot(boot)
  expect_identical(levels(everything$data$panel), c("mean", "variance", "probability"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(print(everything))
  expect_error(plot(boot, what = "slope"), class = "latents_bad_argument")
  expect_error(plot(subset(boot, parameter == "mean")), class = "latents_bad_argument")
})

test_that("the layout splits a kind by level and keeps one point per parameter", {
  table <- data.frame(level = c(rep("measurement", 4L), rep("profile", 4L), "group", "group"),
                      outcome = c("profile_1", "profile_1", "profile_2", "profile_2",
                                  "profile_1", "profile_2", "profile_1", "profile_2",
                                  "group_class_1", "group_class_2"),
                      term = c("a", "b", "a", "b", "group_class_1", "group_class_1",
                               "group_class_2", "group_class_2", NA, NA),
                      parameter = c(rep("mean", 4L), rep("probability", 6L)),
                      estimate = seq(0.1, 1, by = 0.1))
  table$conf_low <- table$estimate - 0.05
  table$conf_high <- table$estimate + 0.05
  columns <- list(labels = c("level", "outcome", "term", "parameter"), kind = "parameter")
  frame <- .latents_bootstrap_layout(table, columns, seq_len(nrow(table)))
  expect_identical(levels(frame$panel), c("mean", "profile probability",
                                          "group probability"))
  # A mean's indicator is the axis and its profile the series; the group
  # classes' probabilities are one series along the group classes.
  expect_identical(as.character(frame$at),
                   c("a", "b", "a", "b", "group_class_1", "group_class_1",
                     "group_class_2", "group_class_2", "group_class_1", "group_class_2"))
  expect_identical(frame$series, c("profile_1", "profile_1", "profile_2", "profile_2",
                                   "profile_1", "profile_2", "profile_1", "profile_2",
                                   NA, NA))
  expect_identical(frame$series_points, rep(2L, 10L))
  expect_false(anyDuplicated(paste(frame$panel, frame$series, frame$at)) > 0L)
  expect_equal(frame$conf_low, table$conf_low, tolerance = 0)
  # With one indicator only the profile varies: the profiles are the axis.
  one <- .latents_bootstrap_layout(table, columns, c(1L, 3L))
  expect_identical(as.character(one$at), c("profile_1", "profile_2"))
  expect_true(all(is.na(one$series)))
})

test_that("a transition bootstrap gets the same accessor", {
  skip_if_not_installed("ggplot2")
  set.seed(8)
  n <- 60L
  first <- sample.int(2L, n, replace = TRUE)
  second <- ifelse(stats::runif(n) < 0.8, first, 3L - first)
  long <- data.frame(id = rep(seq_len(n), each = 2L), time = rep(1:2, n),
                     y = stats::rnorm(2L * n, c(-1.5, 1.5)[c(rbind(first, second))]))
  transition <- lta(long, "y", id = "id", time = "time", n_profiles = 2, n_starts = 2,
                    seed = 1)
  boot <- parameter_inference(transition, data = long, method = "bootstrap", iter = 8,
                              n_starts = 1, seed = 1)
  expect_s3_class(boot, "latents_bootstrap")
  replicates <- as.data.frame(boot, what = "replicates")
  expect_identical(nrow(replicates), attr(boot, "n_valid") * nrow(boot))
  expect_s3_class(plot(boot), "ggplot")
})

test_that("a second-order transition bootstrap keeps the previous profile apart", {
  skip_if_not_installed("ggplot2")
  set.seed(41)
  n <- 120L
  states <- matrix(NA_integer_, n, 4L)
  states[, 1L] <- sample.int(2L, n, replace = TRUE)
  invisible(lapply(2:4, function(t) {
    states[, t] <<- ifelse(stats::runif(n) < 0.85, states[, t - 1L], 3L - states[, t - 1L])
  }))
  s <- as.vector(t(states))
  long <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(1:4, n),
                     y1 = c(0, 3)[s] + stats::rnorm(length(s)),
                     y2 = c(1, -1)[s] + stats::rnorm(length(s)))
  fit <- suppressWarnings(lta(long, c("y1", "y2"), "id", 2, time = "time", n_starts = 3,
                              seed = 1, transitions = "occasion", order = 2))
  boot <- parameter_inference(fit, data = long, method = "bootstrap", iter = 12,
                              n_starts = 1, seed = 3)
  expect_true("previous" %in% names(boot))
  replicates <- as.data.frame(boot, what = "replicates")
  expect_identical(names(replicates), c("replicate", "block", "group_class", "previous",
                                        "from", "to", "term", "estimate"))
  key <- with(replicates, paste(block, group_class, previous, from, to, term))
  rows <- with(boot, paste(block, group_class, previous, from, to, term))
  expect_false(anyDuplicated(rows) > 0L)
  spread <- vapply(rows, function(row) stats::sd(replicates$estimate[key == row]), numeric(1))
  expect_equal(unname(spread), boot$standard_error, tolerance = 1e-12)
  # The previous profile keeps a second-order transition's parameters apart
  # on the plot: no two share a panel, a series and a position.
  drawn <- plot(boot)$data
  expect_false(anyDuplicated(paste(drawn$panel, drawn$series, drawn$at)) > 0L)
  expect_true(all(grepl("profile_. \u00b7 profile_.", drawn$series[drawn$panel == "transition2"])))
})
