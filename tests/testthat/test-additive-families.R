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

test_that("scores equal numerical derivatives for every family structure", {
  skip_if_not_installed("numDeriv")
  stats <- .additive_prepare(families_data, families_vars, "group")
  classes <- c("group_class_1", "group_class_2")
  point <- list(means = rbind(c(-0.5, 0.2), c(0.7, -0.1)),
                between = rbind(c(0.2, 0.4), c(0.35, 0.25)),
                within = rbind(c(0.6, 0.8), c(1.3, 1.7)), weights = c(0.4, 0.6))
  structures <- list(.additive_structure("dispersion"),
                     .additive_structure("additive_dispersion", "varying"),
                     .additive_structure("additive_dispersion", "equal"))
  lapply(structures, \(structure) {
    # Project the point onto the structure: shared blocks take one row.
    shared <- function(m, mode) if (mode == "equal") m[c(1, 1), ] else m
    projected <- list(means = shared(point$means, structure$means),
                      between = shared(point$between, structure$between),
                      within = shared(point$within, structure$within),
                      weights = point$weights)
    theta <- .additive_pack(projected, structure, classes, families_vars)
    analytic <- colSums(.additive_group_scores(stats, projected, structure,
                                               classes, families_vars))
    numeric_gradient <- numDeriv::grad(\(v) .additive_expectation(
      stats, .additive_unpack(v, 2L, 2L, structure))$log_likelihood, theta)
    expect_equal(unname(analytic), numeric_gradient, tolerance = 1e-7)
  })
})

test_that("fitted points are stationary and EM is monotone in the new families", {
  skip_if_not_installed("numDeriv")
  stats <- .additive_prepare(families_data, families_vars, "group")
  cases <- list(list("dispersion", NULL), list("additive_dispersion", "varying"),
                list("additive_dispersion", "equal"))
  lapply(cases, \(case) {
    fit <- multilpa(families_data, families_vars, "group", n_group_classes = 2,
                    family = case[[1]],
                    between_variance = case[[2]] %||% "equal",
                    tol = 1e-10, seed = 3)
    expect_true(fit$converged)
    expect_false(any(fit$between_zero))
    expect_true(all(diff(fit$log_likelihood_history) >=
                      -1e-10 * (1 + abs(fit$log_likelihood))))
    theta <- coef(fit, scale = "unconstrained")
    gradient <- numDeriv::grad(\(v) .additive_expectation(
      stats, .additive_unpack(v, 2L, 2L, fit$structure))$log_likelihood, theta)
    expect_lt(max(abs(gradient)), 1e-3)
  })
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

test_that("shared means are pooled by precision when classes are unequal", {
  skip_if_not_installed("numDeriv")
  # 75/25 classes: an unweighted average of class means would not be the
  # maximum, so stationarity pins down the pooling weights.
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(77)
  sizes <- rep(c(5L, 9L), 80L)
  classes <- ifelse(stats::runif(length(sizes)) < 0.75, 1L, 2L)
  sd_within <- c(0.6, 1.8)
  data <- do.call(rbind, lapply(seq_along(sizes), function(j) {
    intercept <- stats::rnorm(1L, 0.5, 0.5)
    data.frame(group = j, y = intercept +
                 stats::rnorm(sizes[j], 0, sd_within[classes[j]]))
  }))
  # EM alone, without the Newton finish that would mask a wrong M-step: its
  # fixed point must be stationary and every step must raise the likelihood.
  stats <- .additive_prepare(data, "y", "group")
  structure <- .additive_structure("dispersion")
  start <- .additive_initialize(stats, 2L, structure, 1e-6, 1L)
  em <- .additive_em(stats, start, structure, 1e-6, 20000L, 1e-14,
                     matrix(FALSE, 2L, 1L))
  expect_true(em$converged)
  expect_gt(max(em$parameters$weights), 0.65)
  expect_true(all(diff(em$history) >= -1e-10 * (1 + abs(em$history))[-1L]))
  theta <- .additive_pack(em$parameters, structure,
                          c("group_class_1", "group_class_2"), "y")
  gradient <- numDeriv::grad(\(v) .additive_expectation(
    stats, .additive_unpack(v, 2L, 1L, structure))$log_likelihood, theta)
  expect_lt(max(abs(gradient)), 1e-4)
})
