# Regression tests for the 2026-10-20 review: preserve data identities and
# sampling-weight units across paired measurement and tree adapters.

skip_on_cran()

test_that("generated pair names cannot replace caller columns", {
  data <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
  data[["a:b"]] <- seq_len(nrow(data))
  pairs <- .latents_parse_pairs("a:b")
  expect_error(.latents_join_pairs(data, pairs), class = "latents_bad_data")
  expect_identical(data[["a:b"]], seq_len(nrow(data)))
  expect_error(.latents_join_pairs(data[c("a", "a:b")], pairs), class = "latents_bad_data")
  expect_error(lca(data[c("a", "c")], c("a", "b", "c"), 1,
                   local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(lca(data, c("a", "b", "c", "a:b"), 1,
                   local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(multilpa(data, c("a", "b", "c"), "a:b", 1, 1,
                        categorical = c("a", "b", "c"),
                        local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(multilpa(data, c("a", "b", "c"), NULL, 1, 1,
                        categorical = c("a", "b", "c"), profile_covariates = "a:b",
                        local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
})

test_that("paired inputs receive the original categorical validation", {
  data <- expand.grid(a = c(0.3, 0.1 + 0.2), b = 0:1, c = 0:1)
  fit <- function(data, paired) lca(data, c("a", "b", "c"), 1,
    local_dependence = if (paired) "a:b", n_starts = 1, seed = 1)
  invalid <- list(data, transform(data, a = Inf), transform(data, a = 1),
                  transform(data, a = as.Date("2020-01-01") + seq_len(nrow(data))))
  matrix_column <- data
  matrix_column$a <- matrix(seq_len(2 * nrow(data)), ncol = 2)
  list_column <- data
  list_column$a <- as.list(seq_len(nrow(data)))
  invalid <- c(invalid, list(matrix_column, list_column))
  invisible(lapply(invalid, function(frame) {
    expect_error(fit(frame, TRUE), class = "latents_bad_categorical")
    expect_error(fit(frame, FALSE), class = "latents_bad_categorical")
  }))
})

test_that("joined data remain idempotent, serializable, predictable and refittable", {
  data <- expand.grid(a = c(2, 10), b = c(FALSE, TRUE), c = 0:1)
  data$a <- ordered(data$a, levels = c(10, 2))
  pairs <- .latents_parse_pairs("a:b")
  joined <- .latents_join_pairs(data, pairs)
  expect_identical(.latents_join_pairs(joined, pairs), joined)
  expect_identical(.latents_join_pairs(unserialize(serialize(joined, NULL)), pairs), joined)
  fit <- lca(data, c("a", "b", "c"), 1, local_dependence = "a:b", n_starts = 1, seed = 1)
  expected <- predict(fit, type = "posterior")
  expect_equal(predict(fit, data, type = "posterior"), expected)
  expect_equal(predict(fit, joined, type = "posterior"), expected)
  expect_equal(predict(fit, data[1L, , drop = FALSE], type = "posterior"),
               expected[1L, , drop = FALSE])
  frame <- .multilpa_resolve_data(fit, NULL)
  expect_equal(predict(fit, frame, type = "posterior"), expected)
  refit <- do.call(multilpa, c(list(data = frame, n_starts = 1, seed = 1),
                              .multilpa_refit_arguments(fit)))
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-10)
  expect_equal(refit$n_parameters, fit$n_parameters)
  # An existing generated column must not hide changes to its source items.
  joined$a <- rev(joined$a)
  expect_identical(as.character(.latents_join_pairs(joined, pairs)[["a:b"]]),
                   paste(joined$a, joined$b, sep = ":"))
  bad <- data
  bad$a <- c(0.3, 0.1 + 0.2, rep(1, nrow(data) - 2L))
  expect_error(predict(fit, bad, type = "posterior"), class = "latents_bad_categorical")
})

test_that("tree adapters respect row, group and person posterior weights", {
  set.seed(318)
  data <- data.frame(y = c(rnorm(80, -4, .3), rnorm(20, 4, .3)),
                     w = c(rep(1, 80), rep(16, 20)))
  grid <- quietly(enumerate_regressions(y ~ 1, data, n_classes = 1:2,
    weights = "w", n_starts = 2, seed = 27))
  fit <- .latents_tree_fit(grid, 2L)
  draw <- .latents_tree_posteriors(fit, "individuals")
  expect_equal(draw$weights, fit$spec$row_weights)
  expected <- unname(colSums(fit$expectation$tau * fit$spec$row_weights) /
                       sum(fit$spec$row_weights))
  flows <- get_results(grid, "tree")
  expect_equal(flows$to_share[order(flows$to_class)], expected, tolerance = 1e-12)
  expect_equal(flows$mass[order(flows$to_class)], expected, tolerance = 1e-12)
  # Unequal rows per group ensure a recycled group weight cannot pass.
  grouped <- data[rep(seq_len(nrow(data)), rep(c(2, 3), 50)), ]
  grouped$id <- rep(seq_len(nrow(data)), rep(c(2, 3), 50))
  grouped$x <- sequence(rep(c(2, 3), 50))
  grouped$y <- grouped$y + rnorm(nrow(grouped), sd = 0.2)
  group_fit <- quietly(mixture_regression(y ~ 1, grouped, 2, id = "id",
    class_level = "group", weights = "w", n_starts = 1, seed = 1))
  expect_equal(.latents_tree_posteriors(group_fit, "individuals")$weights,
               group_fit$spec$sampling_weights)
  # This deliberately separated fixture checks posterior units; it need not
  # estimate uncertainty for nearly pure group-class memberships.
  two <- quietly(mixture_regression(y ~ 1, grouped, 2, id = "id",
    n_group_classes = 2, weights = "w", vcov_type = "none", n_starts = 1, seed = 1))
  expect_equal(.latents_tree_posteriors(two, "individuals")$weights, two$spec$row_weights)
  expect_equal(.latents_tree_posteriors(two, "groups")$weights, two$spec$sampling_weights)
  growth <- quietly(mixture_regression(y ~ x, grouped, 2, id = "id",
    class_level = "group", random = ~ 1, weights = "w", n_starts = 1, seed = 1))
  expect_s3_class(growth, "latents_growth_mixture")
  expect_equal(.latents_tree_posteriors(growth, "individuals")$weights,
               growth$spec$sampling_weights)
})

test_that("joining preserves original indicator identity validation", {
  data <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
  expect_error(lca(data, c("a", "a", "b", "c"), 1,
                   local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(lca(data, c("a", "b", "b", "c"), 1,
                   local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(multilpa(data, c("a", "b", "c"), NULL, 1, 1,
                        categorical = c("a", "b", "b", "c"),
                        local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_categorical")
  expect_error(multilpa(data, c("a", "b", "c"), "a", 1, 1,
                        categorical = c("a", "b", "c"),
                        local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
})

test_that("pair provenance survives base row combinations without trusting stale metadata", {
  data <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
  pairs <- .latents_parse_pairs("a:b")
  joined <- .latents_join_pairs(data, pairs)
  combined <- rbind(joined[1:3, ], joined[4:8, ])
  expect_equal(as.character(.latents_join_pairs(combined, pairs)[["a:b"]]),
               as.character(joined[["a:b"]]))
  # Single-column selection drops data-frame attributes, retaining the factor
  # provenance; serialization preserves either source of provenance.
  subset <- joined[c("a", "b", "a:b")]
  expect_equal(as.character(.latents_join_pairs(subset, pairs)[["a:b"]]),
               as.character(joined[["a:b"]]))
  expect_equal(as.character(.latents_join_pairs(unserialize(serialize(combined, NULL)),
                                                pairs)[["a:b"]]),
               as.character(joined[["a:b"]]))
  # Raw input with coincidentally correct joint labels is still a collision.
  raw <- data
  raw[["a:b"]] <- factor(paste(raw$a, raw$b, sep = ":"))
  expect_error(.latents_join_pairs(raw, pairs), class = "latents_bad_data")
  # Replacing a generated column leaves stale frame metadata, which cannot
  # authorize overwriting the replacement, even after column subsetting.
  replaced <- combined
  replaced[["a:b"]] <- factor(rep("unrelated", nrow(replaced)))
  expect_error(.latents_join_pairs(replaced, pairs), class = "latents_bad_data")
  expect_error(.latents_join_pairs(replaced[c("a", "b", "a:b")], pairs),
               class = "latents_bad_data")
})

test_that("known fits rebuild pairs after ordinary data operations lose provenance", {
  data <- expand.grid(a = 0:1, b = 0:1, c = 0:1)
  fit <- lca(data, c("a", "b", "c"), 1, local_dependence = "a:b", n_starts = 1, seed = 1)
  joined <- .latents_join_pairs(data, .latents_parse_pairs("a:b"))
  no_attributes <- joined
  attr(no_attributes, "latents_pairs") <- NULL
  attr(no_attributes[["a:b"]], "latents_pair") <- NULL
  variants <- list(subset(joined, a == 0), merge(joined, data.frame(c = 0:1, extra = 1:2)),
                   rbind(joined, joined)[c("a", "b", "a:b", "c")], no_attributes)
  invisible(lapply(variants, function(frame) {
    raw <- frame[c("a", "b", "c")]
    expect_equal(predict(fit, frame, type = "posterior"),
                 predict(fit, raw, type = "posterior"))
    resolved <- .multilpa_resolve_data(fit, frame)
    expect_equal(as.character(resolved[["a:b"]]), paste(frame$a, frame$b, sep = ":"))
  }))
  # Fitting has no established meaning for an existing generated-name column,
  # so matching its labels to raw sources does not authorize its replacement.
  expect_error(lca(no_attributes, c("a", "b", "c"), 1,
                   local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(lca(joined[c("a:b", "c")], c("a", "b", "c"), 1,
                   local_dependence = "a:b", n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_equal(predict(fit, joined[c("a:b", "c")], type = "posterior"),
               predict(fit, data, type = "posterior"))
})
