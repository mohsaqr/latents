test_that("native group identifiers retain the independently calculated ICC", {
  data <- data.frame(id = rep(c(0.3, 0.1 + 0.2), each = 3L),
                     y = c(0, 1, 2, 10, 11, 12))
  expect_equal(length(unique(data$id)), 2L)
  # Between MS = 150, within MS = 1, n0 = 3: (150 - 1)/(150 + 2).
  expect_equal(.multilpa_icc(data$y, data$id), 149 / 152, tolerance = 1e-14)
  expect_equal(descriptives(data, "y", id = "id")$icc, 149 / 152,
               tolerance = 1e-14)
  expect_equal(.multilpa_icc(data$y, rep(c("a", "b"), each = 3)), 149 / 152)
})

test_that("shared Gaussian covariance models still diagnose a covariance floor", {
  expect_warning(fit <- lpa(data.frame(y = rep(c(-1, 1), each = 3)), "y", 2,
    model = "EEI", n_starts = 1, seed = 1), class = "latents_boundary")
  expect_true(fit$converged)
  expect_true(fit$boundary)
  # The default bound is a millionth of the indicator's variance, 1.2 here
  # (review 2026-10-26 round 2, K2-2); it was an absolute 1e-6.
  expect_equal(fit$min_variance, 1e-6 * stats::var(rep(c(-1, 1), each = 3)))
  expect_equal(as.numeric(fit$variances), rep(fit$min_variance, 2))
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
})

test_that("TNA sequences preserve numeric IDs with the same printed label", {
  data <- data.frame(id = rep(c(0.3, 0.1 + 0.2), each = 3L),
                     time = rep(1:3, 2), y = c(0, 1, 2, 10, 11, 12))
  fit <- multilpa(data, "y", "id", 1, 1, time = "time", n_starts = 1, seed = 1)
  sequences <- .multilpa_modal_sequences(fit)
  expect_equal(dim(sequences$states), c(2L, 3L))
  expect_equal(length(sequences$group_class), 2L)
  expect_true(all(as.matrix(sequences$states) == "profile_1"))
  if (requireNamespace("tna", quietly = TRUE)) {
    expect_s3_class(get_tna(fit), "tna")
  }
})

test_that("trust plots distinguish unavailable replication evidence", {
  set.seed(252)
  data <- data.frame(z = rnorm(200))
  state <- rbinom(200, 1, plogis(0.4 * data$z))
  data$y <- rnorm(200, ifelse(state == 1, 2, -2))
  grid <- enumerate_classes(data, "y", id = NULL, n_profiles = 2,
                             n_group_classes = 1, profile_covariates = "z",
                             n_starts = 1, seed = 1)
  # Covariate fits now record replication (RP-1); the unavailable state is
  # what an older saved grid carries, so it is set on the table directly.
  grid$table$n_best_replicated <- NA_integer_
  frame <- .gg_selection_trust(grid, NULL, NULL)$data
  unknown <- is.na(grid$table$n_best_replicated) & grid$table$converged &
    !grid$table$boundary & is.finite(grid$table$log_likelihood)
  expect_true(any(unknown))
  expect_true(all(as.character(frame$status[unknown]) == "replication unavailable"))
})

test_that("mover stayer trees count posterior classes including the stayer", {
  set.seed(141)
  data <- data.frame(id = rep(1:20, each = 4), time = rep(1:4, 20))
  data$y <- rnorm(80, sample(c(-3, 3), 80, replace = TRUE), 0.7)
  grid <- suppressWarnings(enumerate_classes(data, "y", "id", time = "time",
    n_profiles = 2:3, n_group_classes = 1, mover_stayer = TRUE,
    n_starts = 1, seed = 1))
  expect_true(all(grid$table$converged))
  expect_true(all(vapply(grid$fits, function(fit) ncol(fit$group_posteriors),
                         integer(1)) == 2L))
  table <- .latents_tree_table_of(grid)
  expect_equal(table$n_group_classes, c(2L, 2L))
  trees <- .latents_tree_specs(grid, "groups", NULL, NULL, NULL, names(grid$fits))
  expect_length(trees, 1L)
  expect_true(all(grepl("2 group classes", trees[[1]]$rows$label, fixed = TRUE)))
})

test_that("family grids pass supported emEM controls to their direct fits", {
  set.seed(251)
  data <- data.frame(id = rep(1:12, each = 4),
                     y = rep(rnorm(12), each = 4) + rnorm(48))
  grid <- enumerate_classes(data, "y", "id", family = "additive",
    n_group_classes = 1, search = "emEM", n_screen = 3, n_starts = 1, seed = 1)
  direct <- multilpa(data, "y", "id", family = "additive", n_group_classes = 1,
    search = "emEM", n_screen = 3, n_starts = 1, seed = 1)
  expect_equal(grid$fits[[1L]]$log_likelihood, direct$log_likelihood,
               tolerance = 1e-12)
})

test_that("additive BLRT checks exact fitted data and compatible variance bounds", {
  set.seed(151)
  data <- data.frame(id = rep(1:12, each = 4), y = 1e9 + rnorm(48))
  fit <- multilpa(data, "y", "id", family = "additive", n_group_classes = 1,
                   n_starts = 1, seed = 1)
  changed <- data
  changed$y[1L] <- changed$y[1L] + 1
  expect_error(bootstrap_lrt(fit, fit, data = changed, iter = 2),
               class = "latents_bad_inference_data")
  stricter <- fit
  stricter$min_variance <- fit$min_variance * 2
  expect_false(.additive_nested(fit, stricter))
  expect_true(.additive_nested(stricter, fit))
})
