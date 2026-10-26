# Review 2026-10-26: calibration and invariance checks of the group-class
# (additive) families, the cross-level families and the covariance-structure
# charts. References are independent: dense multivariate normal densities,
# numerical derivatives, brute-force likelihoods and mclust's parameter counts.

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

review_structures <- function() {
  list(additive = .additive_structure("additive", "varying"),
       additive_equal = .additive_structure("additive", "equal"),
       dispersion = .additive_structure("dispersion"),
       additive_dispersion = .additive_structure("additive_dispersion", "varying"),
       additive_dispersion_equal = .additive_structure("additive_dispersion", "equal"))
}

review_parameters <- function() {
  list(means = rbind(c(0.1, 0.9), c(2.4, -0.8)),
       between = rbind(c(0.4, 0.7), c(1.1, 0.5)),
       within = rbind(c(1.1, 1.4), c(0.7, 1.9)),
       weights = c(0.55, 0.45))
}

review_transform_fit <- function(fit, shift, scale) {
  moved <- fit
  moved$means <- fit$means * scale + shift
  moved$between_variances <- fit$between_variances * scale^2
  moved$within_variances <- fit$within_variances * scale^2
  stats <- fit$sufficient_statistics
  stats$averages <- stats$averages * scale + shift
  stats$scatter <- stats$scatter * scale^2
  stats$pooled_scatter <- stats$pooled_scatter * scale^2
  moved$sufficient_statistics <- stats
  moved
}

test_that("analytic scores and information match numerical derivatives in every family", {
  skip_if_not_installed("numDeriv")
  data <- review_additive_draw(1L)
  vars <- c("y1", "y2")
  stats <- .additive_prepare(data, vars, "g")
  class_names <- paste0("group_class_", 1:2)
  invisible(lapply(names(review_structures()), function(name) {
    structure <- review_structures()[[name]]
    parameters <- .additive_check_parameters(review_parameters(), 2L)
    parameters$means <- .additive_pool(parameters$means, matrix(1, 2, 2), structure$means)
    parameters$between <- .additive_pool(parameters$between, matrix(1, 2, 2),
                                         structure$between)
    parameters$within <- .additive_pool(parameters$within, matrix(1, 2, 2),
                                        structure$within)
    theta <- .additive_pack(parameters, structure, class_names, vars)
    expect_length(theta, .additive_count_parameters(2L, 2L, structure))
    expect_equal(.additive_pack(.additive_unpack(theta, 2L, 2L, structure), structure,
                                class_names, vars), theta, tolerance = 1e-14,
                 info = name)
    log_likelihood <- function(point) {
      .additive_expectation(stats, .additive_unpack(point, 2L, 2L, structure))$log_likelihood
    }
    analytic <- colSums(.additive_group_scores(
      stats, .additive_unpack(theta, 2L, 2L, structure), structure, class_names, vars))
    expect_equal(unname(analytic), numDeriv::grad(log_likelihood, unname(theta)),
                 tolerance = 1e-6, info = name)
    expect_equal(unname(.additive_information(stats, theta, structure, class_names,
                                              vars, 1e-4)),
                 -numDeriv::hessian(log_likelihood, unname(theta)),
                 tolerance = 1e-6, info = name)
  }))
})

test_that("standard errors are translation invariant and scale equivariant", {
  skip_if_not_installed("numDeriv")
  data <- review_additive_draw(1L)
  fit <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                               family = "additive_dispersion", n_starts = 3L,
                               seed = 1L, tol = 1e-10))
  expect_true(fit$converged)
  standard_errors <- function(object, type) {
    .additive_parameter_table(object, .additive_inference(object, type,
                                                          warn_weak = FALSE))$standard_error
  }
  kind <- .additive_parameter_table(fit, NULL)$parameter
  factor <- ifelse(kind == "mean", 10, ifelse(kind == "variance", 100, 1))
  invisible(lapply(c("observed", "robust", "opg"), function(type) {
    base <- standard_errors(fit, type)
    # A raw-mean step once made the information indefinite at this origin.
    expect_equal(standard_errors(review_transform_fit(fit, 1e6, 1), type), base,
                 tolerance = 1e-8, info = type)
    expect_equal(standard_errors(review_transform_fit(fit, 0, 10), type),
                 factor * base, tolerance = 1e-8, info = type)
  }))
  # The observed covariance is the inverse of minus the numerical Hessian.
  structure <- fit$structure
  stats <- fit$sufficient_statistics
  theta <- coef(fit, scale = "unconstrained")
  log_likelihood <- function(point) {
    .additive_expectation(stats, .additive_unpack(point, 2L, 2L, structure))$log_likelihood
  }
  reference <- solve(-numDeriv::hessian(log_likelihood, unname(theta)))
  observed <- .additive_inference(fit, "observed", warn_weak = FALSE)$vcov
  expect_equal(unname(sqrt(diag(observed))),
               sqrt(diag(reference)), tolerance = 1e-6)
})

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

test_that("fitted maxima cannot be improved and integer weights equal duplicated groups", {
  data <- review_additive_draw(2L, n_groups = 50L)
  vars <- c("y1", "y2")
  fit <- review_quiet(multilpa(data, vars, "g", n_group_classes = 2L,
                               family = "dispersion", n_starts = 3L, seed = 1L,
                               tol = 1e-10))
  structure <- fit$structure
  stats <- fit$sufficient_statistics
  log_likelihood <- function(point) {
    .additive_expectation(stats, .additive_unpack(point, 2L, 2L, structure))$log_likelihood
  }
  polished <- stats::optim(coef(fit, scale = "unconstrained"), log_likelihood,
                           method = "BFGS",
                           control = list(fnscale = -1, reltol = 1e-14, maxit = 500L))
  expect_lt(polished$value - fit$log_likelihood, 1e-8)

  ids <- unique(data$g)
  copies <- rep(1:3, length.out = length(ids))
  data$w <- copies[match(data$g, ids)]
  duplicated_data <- do.call(rbind, lapply(seq_along(ids), function(j) {
    block <- subset(data, g == ids[j])
    do.call(rbind, lapply(seq_len(copies[j]), function(k) {
      transform(block, g = paste0(g, "_", k))
    }))
  }))
  weighted <- review_quiet(multilpa(data, vars, "g", n_group_classes = 2L,
                                    family = "additive_dispersion",
                                    between_variance = "equal", n_starts = 3L,
                                    seed = 1L, tol = 1e-12, weights = "w"))
  copied <- review_quiet(multilpa(duplicated_data, vars, "g", n_group_classes = 2L,
                                  family = "additive_dispersion",
                                  between_variance = "equal", n_starts = 3L,
                                  seed = 1L, tol = 1e-12))
  expect_equal(unname(weighted$means), unname(copied$means), tolerance = 1e-8)
  expect_equal(unname(weighted$within_variances), unname(copied$within_variances),
               tolerance = 1e-8)
  expect_equal(unname(weighted$between_variances), unname(copied$between_variances),
               tolerance = 1e-8)
  # Weights are normalized to the number of groups, so the copied likelihood
  # is the weighted one times the mean weight.
  expect_equal(copied$log_likelihood, mean(copies) * weighted$log_likelihood,
               tolerance = 1e-10)
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

test_that("cross-level likelihoods equal their brute-force definitions", {
  data <- review_additive_draw(3L, n_groups = 30L)
  vars <- c("y1", "y2")
  ids <- unique(data$g)
  brute <- function(fit) {
    n_profiles <- nrow(fit$profile_means)
    n_classes <- length(fit$group_probabilities)
    profile_density <- function(ratings) {
      matrix(vapply(seq_len(n_profiles), function(k) {
        exp(rowSums(stats::dnorm(ratings,
          matrix(fit$profile_means[k, ], nrow(ratings), 2L, byrow = TRUE),
          matrix(sqrt(fit$profile_variances[k, ]), nrow(ratings), 2L, byrow = TRUE),
          log = TRUE)))
      }, numeric(nrow(ratings))), nrow(ratings))
    }
    group_density <- function(average) {
      vapply(seq_len(n_classes), function(h) {
        exp(sum(stats::dnorm(average, fit$group_means[h, ],
                             sqrt(fit$group_mean_variances[h, ]), log = TRUE)))
      }, numeric(1))
    }
    per_group <- vapply(ids, function(id) {
      ratings <- as.matrix(subset(data, g == id, vars))
      members <- profile_density(ratings)
      classes <- group_density(colMeans(ratings))
      if (identical(fit$family, "full_cross_level")) {
        log(sum(fit$group_probabilities * classes *
                  apply(members %*% t(fit$composition), 2L, prod)))
      } else {
        sum(log(members %*% fit$profile_prevalences)) +
          log(sum(fit$group_probabilities * classes))
      }
    }, numeric(1))
    sum(per_group)
  }
  fits <- lapply(c("full_cross_level", "restricted_cross_level"), function(family) {
    review_quiet(multilpa(data, vars, "g", n_profiles = 2L, n_group_classes = 2L,
                          family = family, n_starts = 2L, seed = 1L, tol = 1e-10,
                          max_iter = 3000L))
  })
  invisible(lapply(fits, function(fit) {
    expect_equal(fit$log_likelihood, brute(fit), tolerance = 1e-10, info = fit$family)
    expect_equal(fit$bic, -2 * fit$log_likelihood + fit$n_parameters * log(length(ids)))
    expect_equal(fit$bic_individual,
                 -2 * fit$log_likelihood + fit$n_parameters * log(nrow(data)))
  }))
  # 2 x 2 profile means and variances, prevalences (per class when full),
  # 2 x 2 group means and variances, one class logit.
  expect_identical(fits[[1L]]$n_parameters, 19L)
  expect_identical(fits[[2L]]$n_parameters, 18L)
})

test_that("covariance-structure charts have analytic Jacobians and mclust's counts", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("mclust")
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(7)
  codes <- setdiff(.multilpa_structures(), .multilpa_separable_structures())
  invisible(lapply(codes, function(code) {
    chart <- .multilpa_chart(code, 3L, 3L)
    expect_identical(chart$width, as.integer(mclust::nVarParams(code, d = 3L, G = 3L)),
                     info = code)
    anchor <- if (identical(chart$form, "rotation")) {
      lapply(seq_len(chart$rotation_rows), function(row) {
        qr.Q(qr(matrix(stats::rnorm(9L), 3L)))
      })
    }
    phi <- stats::rnorm(chart$width, sd = 0.3)
    mapped <- .multilpa_chart_map(phi, chart, anchor, derivatives = TRUE)
    numerical <- numDeriv::jacobian(function(point) {
      unlist(.multilpa_chart_map(point, chart, anchor)$sigma)
    }, phi)
    expect_equal(do.call(rbind, mapped$jacobian), numerical, tolerance = 1e-7,
                 info = code)
    # Every profile's covariance is symmetric positive definite.
    expect_true(all(vapply(mapped$sigma, function(sigma) {
      isTRUE(all.equal(sigma, t(sigma))) &&
        min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values) > 0
    }, logical(1))), info = code)
  }))
})

test_that("the covariance gradient is the derivative of the weighted Gaussian log density", {
  skip_if_not_installed("numDeriv")
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(3)
  x <- matrix(stats::rnorm(40L), 20L, 2L)
  weights <- stats::runif(20L)
  covariance <- matrix(c(1.3, 0.4, 0.4, 0.8), 2L)
  parameters <- list(means = matrix(c(0.2, -0.1), 1L),
                     covariances = array(covariance, c(2L, 2L, 1L)))
  expectation <- list(subject_posteriors = matrix(weights), gaussian_moments = list(NULL))
  gradient <- .multilpa_profile_gradients(parameters, expectation, x)
  log_density <- function(entries) {
    sigma <- matrix(entries, 2L)
    residuals <- sweep(x, 2L, parameters$means[1L, ])
    -0.5 * sum(weights * (2 * log(2 * pi) +
                            as.numeric(determinant(sigma)$modulus) +
                            rowSums((residuals %*% t(solve(sigma))) * residuals)))
  }
  expect_equal(as.vector(gradient$covariance[[1L]]),
               numDeriv::grad(log_density, as.vector(covariance)), tolerance = 1e-8)
  mean_density <- function(mean) {
    residuals <- sweep(x, 2L, mean)
    -0.5 * sum(weights * rowSums((residuals %*% solve(covariance)) * residuals))
  }
  expect_equal(as.vector(gradient$means),
               numDeriv::grad(mean_density, c(0.2, -0.1)), tolerance = 1e-8)
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
