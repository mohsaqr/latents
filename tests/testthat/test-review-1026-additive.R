# Review 2026-10-26: calibration and invariance checks of the group-class
# (additive) families, the cross-level families and the covariance-structure
# charts. References are independent: dense multivariate normal densities,
# numerical derivatives, brute-force likelihoods and mclust's parameter counts.
# The comparisons with independent references are in tests/equivalence/test-review-1026-additive.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

review_additive_draw <- function(seed, n_groups = 60L) {
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(seed)
  sizes <- sample(2:7, n_groups, replace = TRUE)
  classes <- sample(1:2, n_groups, replace = TRUE, prob = c(0.6, 0.4))
  means <- rbind(c(0, 1), c(2.5, -1))
  between <- rbind(c(0.5, 0.8), c(1.2, 0.4))
  within <- rbind(c(1, 1.5), c(0.6, 2))
  do.call(rbind, lapply(seq_len(n_groups), function(j) {
    h <- classes[j]
    intercept <- stats::rnorm(2L, means[h, ], sqrt(between[h, ]))
    noise <- matrix(stats::rnorm(sizes[j] * 2L), sizes[j]) %*% diag(sqrt(within[h, ]))
    ratings <- sweep(noise, 2L, intercept, "+")
    colnames(ratings) <- c("y1", "y2")
    data.frame(g = sprintf("g%03d", j), ratings)
  }))
}

review_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("relabelling the group classes permutes the natural standard errors", {
  data <- review_additive_draw(1L)
  fit <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                               family = "additive_dispersion", n_starts = 3L,
                               seed = 1L, tol = 1e-10))
  swapped <- fit
  swapped$means <- fit$means[2:1, ]
  swapped$between_variances <- fit$between_variances[2:1, ]
  swapped$within_variances <- fit$within_variances[2:1, ]
  swapped$group_probabilities <- stats::setNames(fit$group_probabilities[2:1],
                                                 names(fit$group_probabilities))
  rownames(swapped$means) <- rownames(swapped$between_variances) <-
    rownames(swapped$within_variances) <- names(fit$group_probabilities)
  table <- function(object, labels) {
    result <- .additive_parameter_table(object, .additive_inference(object, "observed",
                                                                    warn_weak = FALSE))
    stats::setNames(result$standard_error,
                    paste(result$parameter, result$level, labels[result$group_class],
                          result$indicator))
  }
  original <- table(fit, c(group_class_1 = "A", group_class_2 = "B"))
  relabelled <- table(swapped, c(group_class_1 = "B", group_class_2 = "A"))
  expect_equal(relabelled[names(original)], original, tolerance = 1e-9)
})

test_that("an exactly constant within-group indicator is flagged as zero scatter", {
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(1)
  n_groups <- 20L
  # 0.1, 0.7, ... are not representable, so their group means round and the
  # scatter is of order 1e-32 rather than zero.
  data <- data.frame(g = rep(seq_len(n_groups), each = 3L),
                     y1 = stats::rnorm(3L * n_groups),
                     y2 = rep(rep(c(0.1, 0.7, 1.3, 2.9), length.out = n_groups),
                              each = 3L))
  stats <- .additive_prepare(data, c("y1", "y2"), "g")
  expect_identical(stats$zero_within_scatter, c(y1 = FALSE, y2 = TRUE))
  fit <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 1L,
                               family = "additive", n_starts = 1L, seed = 1L))
  expect_identical(fit$zero_within_scatter, c(y1 = FALSE, y2 = TRUE))
})

test_that("a higher alternative variance floor breaks nesting, by class", {
  data <- review_additive_draw(1L, n_groups = 30L)
  null_model <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 1L,
                                      family = "additive", n_starts = 1L, seed = 1L))
  alternative <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                                       family = "additive", n_starts = 2L, seed = 1L,
                                       min_variance = 1e-3))
  expect_true(.additive_nested(null_model, review_quiet(multilpa(
    data, c("y1", "y2"), "g", n_group_classes = 2L, family = "additive",
    n_starts = 2L, seed = 1L))))
  expect_false(.additive_nested(null_model, alternative))
  expect_error(bootstrap_lrt(null_model, alternative, iter = 1L, seed = 1L),
               class = "latents_bad_argument")
})

test_that("the same data with rows reordered within groups is accepted (ADD-A2)", {
  data <- review_additive_draw(31L, n_groups = 40L)
  fit <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2,
                               family = "additive", n_starts = 2, seed = 1))
  withr::local_seed(2)
  # Rows shuffled inside each group: the sufficient statistics are the same
  # numbers summed in another order.
  reordered <- do.call(rbind, lapply(split(data, data$g), function(group) {
    group[sample.int(nrow(group)), , drop = FALSE]
  }))
  row.names(reordered) <- NULL
  reordered <- reordered[match(unique(data$g), reordered$g) |> (\(first) {
    order(match(reordered$g, unique(data$g)))
  })(), ]
  same <- review_quiet(parameter_inference(fit, data = reordered))
  expect_equal(same$estimate, review_quiet(parameter_inference(fit))$estimate)
  # A changed rating is still refused.
  changed <- data
  changed$y1[[1L]] <- changed$y1[[1L]] + 0.01
  expect_error(parameter_inference(fit, data = changed),
               class = "latents_bad_inference_data")
  null_model <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 1,
                                      family = "additive", n_starts = 1, seed = 1))
  expect_s3_class(review_quiet(bootstrap_lrt(null_model, fit, iter = 2L, n_starts = 1L,
                                             seed = 1L, data = reordered)),
                  "multilpa_bootstrap_lrt")
})
