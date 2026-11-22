# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# descriptives(): counts, means, standard deviations, ranges and the one-way
# random-effects intraclass correlation of the model variables, for a data
# frame and for the data a fit was built from. References: base R's summaries
# (stats), and the ICC(1) from the mean squares of stats::aov() with the
# unbalanced-design group size n0 (Shrout and Fleiss 1979; Snijders and
# Bosker 2012, eq. 3.12): (MSB - MSW) / (MSB + (n0 - 1) MSW).

# The descriptive table a reader would write in base R, one row per variable.
reference_descriptives <- function(data, vars, id) {
  rows <- lapply(vars, \(variable) {
    values <- data[[variable]]
    kept <- !is.na(values)
    observed <- values[kept]
    groups <- factor(data[[id]][kept])
    squares <- summary(stats::aov(observed ~ groups))[[1L]][["Mean Sq"]]
    sizes <- as.vector(table(groups))
    k <- length(sizes)
    n0 <- (sum(sizes) - sum(sizes^2) / sum(sizes)) / (k - 1)
    data.frame(variable = variable, n = length(observed), n_missing = sum(!kept), mean = mean(observed),
               sd = stats::sd(observed), min = min(observed), max = max(observed),
               n_distinct = length(unique(observed)), n_groups = k,
               icc = (squares[[1L]] - squares[[2L]]) / (squares[[1L]] + (n0 - 1) * squares[[2L]]))
  })
  do.call(rbind, rows)
}

two_level_data <- function() {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(2)
  # Unequal group sizes, so that n0 differs from the mean group size.
  sizes <- rep(c(5L, 8L, 11L), 10L)
  g <- rep(seq_along(sizes), sizes)
  u <- stats::rnorm(length(sizes))[g]
  n <- length(g)
  data <- data.frame(g = g, y1 = stats::rnorm(n, ifelse(seq_len(n) %% 2 == 0, 0, 3)) + u,
                     y2 = stats::rnorm(n) + u, z = stats::rnorm(n))
  data$y2[c(4L, 50L)] <- NA
  data
}

test_that("a data frame's descriptives equal base R and the ANOVA intraclass correlation", {
  data <- two_level_data()
  # Sums, square roots and one ratio of mean squares: rounding only.
  expect_reference(descriptives(data, vars = c("y1", "y2", "z"), id = "g"),
                   reference_descriptives(data, c("y1", "y2", "z"), "g"), reference = "stats",
                   covers = c("descriptives()", "descriptives.data.frame()"),
                   quantity = "descriptives of three variables, unequal groups", tolerance = 1e-10,
                   measure = "relative", key = "variable",
                   note = "ICC(1) from the mean squares of stats::aov(); n0 = (N - sum n_j^2 / N) / (k - 1).")
})

test_that("a fit's descriptives are those of the data it was fitted to", {
  data <- two_level_data()
  complete <- subset(data, !is.na(y2))
  profile_fit <- multilpa(complete, c("y1", "y2"), id = "g", n_profiles = 2, n_group_classes = 1,
                          n_starts = 2, seed = 1)
  expected <- reference_descriptives(complete, c("y1", "y2"), "g")
  expect_reference(descriptives(profile_fit), expected, reference = "stats",
                   covers = c("multilpa()", "descriptives.multilpa()"), quantity = "descriptives of a multilpa() fit",
                   tolerance = 1e-10, measure = "relative", key = "variable")
  covariate_fit <- multilpa(complete, c("y1", "y2"), id = "g", n_profiles = 2, n_group_classes = 1,
                            profile_covariates = "z", n_starts = 2, seed = 1)
  expect_reference(descriptives(covariate_fit), expected, reference = "stats",
                   covers = c("multilpa()", "descriptives.multilpa_covariates()"),
                   quantity = "descriptives of a fit with covariates", tolerance = 1e-10, measure = "relative",
                   key = "variable")
  sequences <- subset(course_engagement, student <= 30)
  transition_fit <- lta(sequences, c("browse", "lectures"), id = "student", n_profiles = 2, time = "sequence",
                        n_starts = 2, seed = 1)
  expect_reference(descriptives(transition_fit), reference_descriptives(sequences, c("browse", "lectures"), "student"),
                   reference = "stats", covers = c("lta()", "descriptives.multilpa_transitions()"),
                   quantity = "descriptives of an lta() fit", tolerance = 1e-10, measure = "relative", key = "variable")
})
