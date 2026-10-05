# Dispersion and additive-dispersion families: the same engine with
# class-specific within-group variances. References: dense covariance
# densities, numerical derivatives, closed-form nesting relations.

families_truth <- list(
  means = rbind(c(-0.8, 0.3), c(0.9, -0.4)),
  between = rbind(c(0.3, 0.3), c(0.3, 0.3)),
  within = rbind(c(0.4, 0.5), c(1.6, 2.0)),
  weights = c(0.5, 0.5))
families_data <- local({
  # Class-specific within variances: draw each group's ratings with its
  # class's residual SDs.
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(31)
  sizes <- rep(c(4L, 6L, 8L, 10L), 40L)
  classes <- sample.int(2L, length(sizes), replace = TRUE)
  do.call(rbind, lapply(seq_along(sizes), function(j) {
    h <- classes[j]
    intercept <- stats::rnorm(2L, families_truth$means[h, ],
                              sqrt(families_truth$between[h, ]))
    ratings <- sweep(matrix(stats::rnorm(sizes[j] * 2L), sizes[j], 2L) %*%
                       diag(sqrt(families_truth$within[h, ])), 2L, intercept, "+")
    data.frame(group = sprintf("g%03d", j), class = h,
               y1 = ratings[, 1L], y2 = ratings[, 2L])
  }))
})
families_vars <- c("y1", "y2")

test_that("class-specific within variances: compact density equals dense", {
  stats <- .additive_prepare(families_data, families_vars, "group")
  compact <- .additive_log_density(stats, families_truth)
  blocks <- additive_blocks(families_data, families_vars, "group")
  dense <- t(vapply(blocks, \(x) vapply(1:2, \(h) additive_dense_density(
    x, families_truth$means[h, ], families_truth$between[h, ],
    families_truth$within[h, ]), numeric(1)), numeric(2)))
  expect_equal(unname(compact), dense, tolerance = 1e-10)
})

test_that("blocks are shared or class-specific as each family says", {
  dispersion <- multilpa(families_data, families_vars, "group",
                         n_group_classes = 2, family = "dispersion", seed = 3)
  expect_equal(dispersion$means[1, ], dispersion$means[2, ])
  expect_equal(dispersion$between_variances[1, ], dispersion$between_variances[2, ])
  expect_false(isTRUE(all.equal(dispersion$within_variances[1, ],
                                dispersion$within_variances[2, ])))
  expect_identical(dispersion$n_parameters, 2L + 2L + 4L + 1L)
  both <- multilpa(families_data, families_vars, "group", n_group_classes = 2,
                   family = "additive_dispersion", seed = 3)
  expect_identical(both$n_parameters, 4L + 4L + 4L + 1L)
  both_equal <- multilpa(families_data, families_vars, "group",
                         n_group_classes = 2, family = "additive_dispersion",
                         between_variance = "equal", seed = 3)
  expect_identical(both_equal$n_parameters, 4L + 2L + 4L + 1L)
  table <- get_results(dispersion)
  expect_identical(subset(table, parameter == "mean")$group_class,
                   c("shared", "shared"))
  expect_identical(subset(table, level == "within")$group_class,
                   rep(c("group_class_1", "group_class_2"), each = 2L))
  expect_identical(nrow(table), length(coef(dispersion)))
  expect_identical(names(coef(dispersion)), rownames(vcov(dispersion)))
})

test_that("nested families order their likelihoods", {
  fits <- lapply(c("additive", "dispersion", "additive_dispersion"), \(family) {
    multilpa(families_data, families_vars, "group", n_group_classes = 2,
             family = family, between_variance = "equal", tol = 1e-10, seed = 3)
  })
  names(fits) <- c("additive", "dispersion", "additive_dispersion")
  # Additive-dispersion contains both, so it cannot fit worse than either.
  expect_gte(fits$additive_dispersion$log_likelihood,
             fits$additive$log_likelihood - 1e-6)
  expect_gte(fits$additive_dispersion$log_likelihood,
             fits$dispersion$log_likelihood - 1e-6)
  # With one class every family is the same model.
  one <- vapply(c("additive", "dispersion", "additive_dispersion"), \(family) {
    multilpa(families_data, families_vars, "group", n_group_classes = 1,
             family = family, between_variance = "equal", tol = 1e-12,
             seed = 1)$log_likelihood
  }, numeric(1))
  expect_equal(unname(one), rep(one[[1]], 3), tolerance = 1e-8)
})

test_that("the dispersion family recovers class-specific within variances", {
  fit <- multilpa(families_data, families_vars, "group", n_group_classes = 2,
                  family = "dispersion", seed = 3)
  recovery <- get_results(fit, "recovery", data = families_data, truth = "class")
  agreement <- sum(tapply(recovery$n, recovery$assigned, max)) / sum(recovery$n)
  expect_gt(agreement, 0.9)
  # Labels are arbitrary: the class with the smaller within variance should
  # match the truth's quiet class (0.4, 0.5), the other the noisy one.
  quiet <- which.min(fit$within_variances[, 1])
  expect_equal(unname(fit$within_variances[quiet, ]), c(0.4, 0.5), tolerance = 0.3)
  expect_equal(unname(fit$within_variances[3 - quiet, ]), c(1.6, 2.0),
               tolerance = 0.3)
})

test_that("dispersion refuses class-specific between variances", {
  expect_error(multilpa(families_data, families_vars, "group",
                        n_group_classes = 2, family = "dispersion",
                        between_variance = "varying"),
               class = "latents_bad_argument")
  expect_s3_class(multilpa(families_data, families_vars, "group",
                           n_group_classes = 2, family = "dispersion",
                           between_variance = "equal", seed = 1),
                  "multilpa_additive")
})

test_that("print, summary and plots work for the new families", {
  skip_if_not_installed("ggplot2")
  fits <- lapply(c("dispersion", "additive_dispersion"), \(family) {
    multilpa(families_data, families_vars, "group", n_group_classes = 2,
             family = family, seed = 3)
  })
  expect_output(print(fits[[1]]), "Dispersion group-class model")
  expect_output(print(fits[[2]]), "Additive-dispersion group-class model")
  lapply(fits, \(fit) {
    expect_s3_class(summary(fit), "summary_multilpa_additive")
    plots <- lapply(c("means", "variances", "intercepts"), \(view) plot(fit, what = view))
    expect_true(all(vapply(plots, ggplot2::is_ggplot, logical(1))))
  })
})
