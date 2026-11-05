skip_on_cran()

test_that("grand moments agree with weighted population moments and replication", {
  x <- cbind(a = c(0, 2, 5, 100), b = c(1, NA, 9, -100))
  weights <- c(1, 2, 4, 0)
  transform <- function(w) .multilpa_center_indicators(
    x, seq_len(nrow(x)), nrow(x), "grand", "grand", row_weights = w)
  actual <- transform(weights)
  center <- vapply(seq_len(ncol(x)), function(j) {
    kept <- !is.na(x[, j]) & weights > 0
    stats::weighted.mean(x[kept, j], weights[kept])
  }, numeric(1))
  scale <- vapply(seq_len(ncol(x)), function(j) {
    kept <- !is.na(x[, j]) & weights > 0
    sqrt(stats::weighted.mean((x[kept, j] - center[j])^2, weights[kept]))
  }, numeric(1))
  expect_equal(unname(actual$offsets[1, ]), center)
  expect_equal(unname(actual$scales[1, ]), scale)
  expect_equal(transform(weights * 17), actual)
  duplicated <- x[rep(seq_len(nrow(x)), weights), , drop = FALSE]
  repeated <- .multilpa_center_indicators(duplicated, seq_len(nrow(duplicated)),
    nrow(duplicated), "grand", "grand")
  expect_equal(actual$offsets, repeated$offsets)
  expect_equal(actual$scales, repeated$scales)
  expect_equal(actual$x, sweep(sweep(x, 2L, center), 2L, scale, "/"))
})

test_that("person centering refuses partial singletons and ignores absent groups", {
  x <- cbind(a = c(4, 1, 3), b = c(9, 2, 4))
  expect_error(.multilpa_center_indicators(x, c(1L, 2L, 2L), 2L, "person"),
               "at least two observed", class = "latents_bad_data")
  result <- .multilpa_center_indicators(x, c(1L, 2L, 2L), 2L, "person",
                                        carrying = c(FALSE, TRUE))
  expect_equal(result$x[2:3, ], matrix(c(-1, 1, -1, 1), 2L,
                                      dimnames = list(NULL, c("a", "b"))))
  partial <- rbind(x[2:3, ], c(2, 3), c(NA, 4))
  expect_error(.multilpa_center_indicators(partial, c(1L, 1L, 2L, 2L), 2L, "person"),
               "`a`", class = "latents_bad_data")
})

test_that("weighted grand fits reproduce integer-expanded data with missing cells", {
  set.seed(110511)
  n <- 180L
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = rnorm(n, c(-2, 2)[class], .6),
                     b = rnorm(n, c(2, -2)[class], 1), w = c(1, 3)[class])
  data$a[seq(3, n, 13)] <- NA_real_
  data$b[seq(7, n, 17)] <- NA_real_
  fit <- function(data, weights = NULL) quietly(lpa(data, c("a", "b"), 2,
    model = "EII", weights = weights, missing = "fiml", centering = "grand",
    scaling = "grand", n_starts = 3, seed = 12, tol = 1e-10))
  weighted <- fit(data, "w")
  repeated <- fit(data[rep(seq_len(n), data$w), , drop = FALSE])
  expect_equal(weighted$centering_offsets, repeated$centering_offsets)
  expect_equal(weighted$scaling_factors, repeated$scaling_factors)
  expect_equal(weighted$means, repeated$means, tolerance = 1e-5)
  expect_equal(weighted$variances, repeated$variances, tolerance = 1e-5)
  expect_equal(weighted$log_likelihood * sum(data$w) / n,
               repeated$log_likelihood, tolerance = 1e-7)
})

test_that("categorical bounds inform diagnostics and LRT without changing BLRT flags", {
  set.seed(1106)
  n <- 400L
  class <- sample(1:2, n, TRUE)
  data <- data.frame(a = ifelse(class == 1L, "x", "y"),
    b = ifelse(class == 1L, "u", "v"),
    c = ifelse(runif(n) < c(.85, .2)[class], "r", "s"),
    d = ifelse(runif(n) < c(.9, .3)[class], "p", "q"),
    e = ifelse(runif(n) < c(.75, .2)[class], "m", "n"))
  fit <- function(...) quietly(lca(data, names(data), 2, n_starts = 2,
                                  seed = 1105, tol = 1e-10, ...))
  null <- fit()
  alternative <- fit(local_dependence = "c:d")
  expect_false(null$boundary)
  expect_true(.latents_boundary_active(null))
  checks <- get_checks(null)
  expect_identical(checks$status[checks$check == "estimates at a bound"], "check")
  expect_warning(compare_models(null, alternative, test = "lrt"),
                 class = "latents_boundary")
  expect_error(parameter_inference(null), class = "latents_boundary_fit")
  grid <- quietly(enumerate_lca(data, names(data), n_classes = 2,
                                n_starts = 2, seed = 1105))
  expect_true(all(get_results(grid, "candidates")$boundary))
})

test_that("categorical simulators retain non-syntactic and joint item names", {
  set.seed(1116)
  id <- rep(seq_len(45L), each = 5L)
  group <- sample(1:2, 45L, TRUE)
  class <- 1L + rbinom(length(id), 1L, c(.15, .85)[group[id]])
  data <- data.frame(id = id,
    a = ifelse(runif(length(id)) < c(.2, .8)[class], "x", "y"),
    b = ifelse(runif(length(id)) < c(.7, .3)[class], "u", "v"),
    y = rnorm(length(id), c(-1.5, 1.5)[class]))
  fit <- function(groups) quietly(multilpa(data, c("a", "b", "y"), "id",
    n_profiles = 2, n_group_classes = groups, categorical = c("a", "b"),
    local_dependence = "a:b", n_starts = 3, seed = 1105, max_iter = 3000))
  null <- fit(1)
  alternative <- fit(2)
  expect_identical(names(.multilpa_simulate(null)), c(null$vars, null$id))
  result <- quietly(bootstrap_lrt(null, alternative, iter = 3, n_starts = 2, seed = 1105))
  expect_equal(as.data.frame(result)$n_valid, 3L)
  expect_true(is.finite(as.data.frame(result)$p_value))
  names(data)[names(data) == "a"] <- "item one"
  names(data)[names(data) == "b"] <- "item-two"
  unusual <- quietly(lca(data, c("item one", "item-two"), 1, n_starts = 1))
  expect_identical(names(.multilpa_simulate(unusual)), c(unusual$vars, unusual$id))
})

test_that("general transition measurement bounds also reach diagnostic checks", {
  set.seed(110518)
  n <- 240L
  class <- sample(1:2, n, TRUE)
  data <- data.frame(id = rep(seq_len(80L), each = 3L), time = rep(1:3, 80L), w = 1,
    a = ifelse(class == 1, "yes", "no"), b = ifelse(class == 1, "up", "down"),
    c = ifelse(runif(n) < c(.85, .2)[class], "yes", "no"),
    d = ifelse(runif(n) < c(.9, .3)[class], "yes", "no"))
  fit <- quietly(lta(data, c("a", "b", "c", "d"), "id", 2, "time",
    categorical = c("a", "b", "c", "d"), weights = "w", n_group_classes = 1,
    n_starts = 2, seed = 18))
  expect_s3_class(fit, "multilpa_lta")
  expect_false(fit$boundary)
  expect_true(.latents_boundary_active(fit))
  checks <- get_checks(fit)
  expect_identical(checks$status[checks$check == "estimates at a bound"], "check")
})

test_that("classification gaps compare weighted modal and posterior shares", {
  set.seed(110512)
  class <- rep(1:2, each = 60L)
  data <- data.frame(a = rnorm(120L, c(-6, 6)[class], .2),
                     b = rnorm(120L, c(6, -6)[class], .2), w = c(1, 5)[class])
  fit <- quietly(lpa(data, c("a", "b"), 2, weights = "w", n_starts = 1, seed = 4))
  table <- get_checks(fit)
  expect_lt(table$value[table$check == "modal vs estimated proportions"], 1e-5)
  basis <- .multilpa_standardization(fit, c("a", "b"))
  weights <- unname(fit$sampling_weights)
  expected <- vapply(c("a", "b"), function(name) stats::weighted.mean(data[[name]], weights),
                       numeric(1))
  expect_equal(basis$centre, unname(expected))
  expect_equal(basis$spread, unname(vapply(c("a", "b"), function(name) {
    sqrt(stats::weighted.mean((data[[name]] - expected[name])^2, weights))
  }, numeric(1))))
  expect_warning(described <- descriptives(fit), class = "latents_unweighted_descriptives")
  expect_equal(described$mean, vapply(data[c("a", "b")], mean, numeric(1)),
               ignore_attr = TRUE)
})

test_that("printed small slopes keep digits beside a large intercept", {
  table <- data.frame(estimate = c(200, .004, -.003))
  actual <- .latents_cells(table, .latents_column("Estimate", "estimate", "number"))
  expect_identical(actual, c("200.0000", "0.0040", "-0.0030"))
  names <- data.frame(term = c("all", "common"))
  expect_identical(.latents_cells(names, .latents_column("Term", "term", "label")),
                    names$term)
})

test_that("standardized references agree under unit weights and integer replication", {
  set.seed(110516)
  class <- rep(0:1, each = 30L)
  data <- data.frame(a = rnorm(60, 3 * class, .6),
                     b = rnorm(60, -2 * class, .8), w = ifelse(class == 0, 3L, 1L))
  fit <- function(data, weights = NULL) quietly(lpa(data, c("a", "b"), 2,
    weights = weights, n_starts = 2, seed = 17))
  weighted <- fit(data, "w")
  repeated <- fit(data[rep(seq_len(nrow(data)), data$w), , drop = FALSE])
  expect_equal(.multilpa_standardization(weighted, c("a", "b")),
               .multilpa_standardization(repeated, c("a", "b")))
  data$one <- 1
  expect_equal(.multilpa_standardization(fit(data, "one"), c("a", "b")),
               .multilpa_standardization(fit(data), c("a", "b")))
})

test_that("zero-weight missing Gaussian values do not enter the data contract", {
  set.seed(110513)
  data <- data.frame(a = rnorm(80), b = rnorm(80), w = c(0, rep(1, 79)))
  data$a[1] <- NA_real_
  fit <- quietly(lpa(data, c("a", "b"), 1, weights = "w", n_starts = 1))
  control <- quietly(lpa(data[-1, ], c("a", "b"), 1, n_starts = 1))
  expect_equal(fit$log_likelihood, control$log_likelihood)
  expect_equal(fit$means, control$means)
  expect_equal(nrow(predict(fit)), 80L)
})

test_that("unconverged sensitivity refits cannot be called the best optimum", {
  set.seed(110514)
  class <- rep(1:2, each = 60L)
  data <- data.frame(a = rnorm(120, c(-1, 1)[class]), b = rnorm(120, c(1, -1)[class]))
  fit <- quietly(lpa(data, c("a", "b"), 2, n_starts = 1, seed = 9))
  result <- quietly(sensitivity(fit, seeds = 1:3, n_starts = 1, max_iter = 1,
                                tol = 1e-12))
  expect_true(all(!result$converged))
  expect_true(all(is.na(result$optimum)))
  expect_false(any(result$best))
})

test_that("cross-level prevalence and sample counts refer to weighted individuals", {
  set.seed(1127)
  data <- do.call(rbind, lapply(seq_len(30L), function(group) {
    type <- 1L + group %% 2L
    n <- c(3L, 9L)[type]
    class <- sample(1:2, n, TRUE, prob = if (type == 1L) c(.9, .1) else c(.2, .8))
    data.frame(id = group, a = rnorm(n, c(-2, 2)[class], .6),
               b = rnorm(n, c(-1, 1)[class], .6), w = if (group <= 3) 0 else type)
  }))
  invisible(lapply(c("full_cross_level", "restricted_cross_level"), function(family) {
    fit <- quietly(multilpa(data, c("a", "b"), "id", 2, 2, family = family,
                            weights = "w", n_starts = 2, seed = 1127))
    weights <- fit$sampling_weights[fit$sufficient_statistics$index]
    expected <- colSums(fit$subject_posteriors * weights) / sum(weights)
    profile <- unique(subset(get_results(fit, "profiles"), level == "individual")
                      [c("profile", "prevalence")])
    expect_equal(profile$prevalence, unname(expected))
    table <- get_results(fit, "fit")
    expect_equal(table$n_groups, length(unique(data$id[data$w > 0])))
    expect_equal(table$n_obs, sum(data$w > 0))
    expect_output(print(fit), sprintf("%d groups, %d observations", table$n_groups, table$n_obs))
  }))
})
