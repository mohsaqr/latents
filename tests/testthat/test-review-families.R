# Independent checks added in the 2026-10-01 model review.
review_family_data <- function() {
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(701)
    sizes <- rep(c(4L, 7L, 10L), 30L)
    class <- rep(c(1L, 1L, 2L), 30L)
    do.call(rbind, lapply(seq_along(sizes), function(j) {
      mean <- c(-1.6, -1.2) + 3.2 * (class[j] == 2L)
      intercept <- rnorm(2L, mean, sqrt(c(0.4, 0.6)))
      rating <- sweep(matrix(rnorm(sizes[j] * 2L), sizes[j], 2L),
                      2L, sqrt(c(0.5, 0.8) * ifelse(class[j] == 1L, 1, 3)), "*")
      rating <- sweep(rating, 2L, intercept, "+")
      data.frame(id = j, y1 = rating[, 1L], y2 = rating[, 2L],
                 w = c(0.4, 1.3, 2.8)[1L + (j %% 3L)])
    }))
}

test_that("reviewed family scores and information match dense Gaussian derivatives", {
  skip_if_not_installed("numDeriv")
  data <- review_family_data()
  vars <- c("y1", "y2")
  stats <- .additive_prepare(data, vars, "id", "w")
  blocks <- lapply(unique(data$id), function(j) as.matrix(data[data$id == j, vars]))
  cases <- list(c("additive", "varying"), c("additive", "equal"),
                c("dispersion", "equal"), c("additive_dispersion", "varying"),
                c("additive_dispersion", "equal"))
  invisible(lapply(cases, function(case) {
    structure <- .additive_structure(case[1L], case[2L])
    shared <- function(value, mode) if (mode == "equal") value[c(1L, 1L), ] else value
    point <- list(means = shared(rbind(c(-1, -0.5), c(1, 0.5)), structure$means),
                  within = shared(rbind(c(0.6, 0.8), c(1.5, 1.8)), structure$within),
                  between = shared(rbind(c(0.3, 0.4), c(0.7, 0.5)), structure$between),
                  weights = c(0.7, 0.3))
    classes <- c("group_class_1", "group_class_2")
    theta <- .additive_pack(point, structure, classes, vars)
    # Build a full covariance for every group instead of using the engine's
    # group sufficient-statistic density or likelihood routine.
    objective <- function(theta) {
      p <- .additive_unpack(theta, 2L, 2L, structure)
      sum(stats$sampling_weights * vapply(blocks, function(x) {
        joint <- vapply(seq_len(2L), function(h) {
          additive_dense_density(x, p$means[h, ], p$between[h, ], p$within[h, ]) +
            log(p$weights[h])
        }, numeric(1))
        max(joint) + log(sum(exp(joint - max(joint))))
      }, numeric(1)))
    }
    scores <- .additive_group_scores(stats, point, structure, classes, vars)
    expect_equal(unname(colSums(scores)), numDeriv::grad(objective, theta),
                 tolerance = 1e-6, info = paste(case, collapse = "/"))
    information <- .additive_information(stats, theta, structure, classes, vars, 1e-4)
    expect_equal(unname(information), -numDeriv::hessian(objective, theta),
                 tolerance = 1e-4, info = paste(case, collapse = "/"))
    expect_equal(.additive_expectation(stats, point)$log_likelihood, objective(theta),
                 tolerance = 1e-10)
  }))
})

test_that("group-class bootstrap refuses changed scatter and group sizes", {
  data <- review_family_data()
  null <- multilpa(data, c("y1", "y2"), "id", n_group_classes = 1,
                   family = "additive", seed = 1, n_starts = 1)
  alternative <- multilpa(data, c("y1", "y2"), "id", n_group_classes = 2,
                          family = "additive", seed = 1, n_starts = 2)
  # Same group averages, different raw-rating likelihood.
  changed <- data
  changed$y1[1:2] <- changed$y1[1:2] + c(-0.5, 0.5)
  expect_error(bootstrap_lrt(null, alternative, data = changed, iter = 1),
               class = "latents_bad_inference_data")
  enlarged <- rbind(data, transform(data[1L, ],
                                   y1 = mean(data$y1[data$id == 1L]),
                                   y2 = mean(data$y2[data$id == 1L])))
  other <- multilpa(enlarged, c("y1", "y2"), "id", n_group_classes = 2,
                    family = "additive", seed = 1, n_starts = 2)
  expect_error(bootstrap_lrt(null, other, iter = 1),
               class = "latents_bad_inference_data")
})

test_that("weighted family EM is stationary before Newton polishing", {
  data <- review_family_data()
  cases <- list(c("additive", "varying"), c("additive", "equal"),
                c("dispersion", "equal"), c("additive_dispersion", "varying"),
                c("additive_dispersion", "equal"))
  invisible(lapply(cases, function(case) {
    ratings <- data
    if (case[1L] == "dispersion") {
      # Shared true locations leave the classes distinguished by their spread.
      ratings[, c("y1", "y2")] <- ratings[, c("y1", "y2")] -
        3.2 * (ratings$id %% 3L == 0L)
    }
    stats <- .additive_prepare(ratings, c("y1", "y2"), "id", "w")
    structure <- .additive_structure(case[1L], case[2L])
    start <- .additive_initialize(stats, 2L, structure, 1e-6, 1L)
    fit <- .additive_em(stats, start, structure, 1e-6, 5000L, 1e-12,
                        matrix(FALSE, 2L, 2L))
    expect_true(fit$converged)
    expect_gt(min(fit$parameters$between), 0.1)
    expect_gte(min(diff(fit$history)), -1e-8)
    # These scores were independently checked against the dense likelihood
    # for every restriction above; Newton cannot repair this EM-only check.
    scores <- .additive_group_scores(stats, fit$parameters, structure,
                                     c("group_class_1", "group_class_2"), stats$vars)
    expect_lt(max(abs(colSums(scores))), 1e-3)
  }))
})

test_that("one-class equivalent group-class models cannot form a bootstrap LRT", {
  data <- review_family_data()
  additive <- multilpa(data, c("y1", "y2"), "id", n_group_classes = 1,
                       family = "additive", seed = 1, n_starts = 1)
  dispersion <- multilpa(data, c("y1", "y2"), "id", n_group_classes = 1,
                         family = "dispersion", seed = 1, n_starts = 1)
  expect_error(bootstrap_lrt(additive, dispersion, iter = 1),
               class = "latents_bad_argument")
})

test_that("family inference validates finite levels and positive difference steps", {
  fit <- multilpa(review_family_data(), c("y1", "y2"), "id", n_group_classes = 1,
                  family = "additive", seed = 1, n_starts = 1)
  invisible(lapply(list(0, 1, NA_real_, Inf, c(0.8, 0.9)), function(level) {
    expect_error(parameter_inference(fit, level = level), "level")
    expect_error(confint(fit, level = level), "level")
    expect_error(summary(fit, level = level), "level")
  }))
  invisible(lapply(list(0, -1, NA_real_, Inf, c(1e-4, 1e-5)), function(step) {
    expect_error(parameter_inference(fit, step = step), "step")
  }))
})

test_that("family bootstrap validates controls before fitting replicates", {
  data <- review_family_data()
  null <- multilpa(data, c("y1", "y2"), "id", n_group_classes = 1,
                   family = "additive", seed = 1, n_starts = 1)
  alternative <- multilpa(data, c("y1", "y2"), "id", n_group_classes = 2,
                          family = "additive", seed = 1, n_starts = 2)
  controls <- list(iter = c(Inf, NA, 0, 1.2), n_starts = c(Inf, NA, 0, 1.2),
                   max_iter = c(Inf, NA, 0, 1.2), tol = c(Inf, NA, 0, -1))
  invisible(lapply(names(controls), function(name) {
    invisible(lapply(controls[[name]], function(value) {
      expect_error(do.call(bootstrap_lrt, c(list(null_model = null,
                                                alternative_model = alternative),
                                           stats::setNames(list(value), name))), name)
    }))
  }))
})

# Direct evaluation of the manifest working likelihood, using only the
# returned parameters, raw ratings and ordinary univariate Gaussian densities.
review_cross_loglik <- function(fit, data, prevalence = NULL) {
  ids <- unique(data$id)
  weights <- if (is.null(fit$sampling_weights)) rep(1, length(ids)) else
    unname(fit$sampling_weights)
  if (is.null(prevalence)) {
    prevalence <- if (fit$family == "full_cross_level") fit$composition else {
      row_weights <- weights[match(data$id, ids)]
      matrix(colSums(fit$subject_posteriors * row_weights) / sum(row_weights),
             length(fit$group_probabilities), nrow(fit$profile_means), byrow = TRUE)
    }
  }
  sum(weights * vapply(ids, function(id) {
    x <- as.matrix(data[data$id == id, fit$vars, drop = FALSE])
    members <- matrix(vapply(seq_len(nrow(fit$profile_means)), function(k) {
      exp(rowSums(vapply(seq_along(fit$vars), function(r) {
        dnorm(x[, r], fit$profile_means[k, r], sqrt(fit$profile_variances[k, r]),
               log = TRUE)
      }, numeric(nrow(x)))))
    }, numeric(nrow(x))), nrow(x))
    joint <- vapply(seq_along(fit$group_probabilities), function(h) {
      log(fit$group_probabilities[h]) +
        sum(dnorm(colMeans(x), fit$group_means[h, ],
                   sqrt(fit$group_mean_variances[h, ]), log = TRUE)) +
        sum(log(drop(members %*% prevalence[h, ])))
    }, numeric(1))
    max(joint) + log(sum(exp(joint - max(joint))))
  }, numeric(1)))
}

test_that("all cross-level variance restrictions reproduce weighted raw formulas", {
  skip_if_not_installed("numDeriv")
  data <- review_family_data()
  grid <- expand.grid(family = c("restricted_cross_level", "full_cross_level"),
                       within = c("equal", "varying"),
                       between = c("equal", "varying"), stringsAsFactors = FALSE)
  invisible(lapply(seq_len(nrow(grid)), function(i) {
    fit <- multilpa(data, c("y1", "y2"), "id", n_profiles = 2,
                    n_group_classes = 2, family = grid$family[i], weights = "w",
                    variance_model = grid$within[i], between_variance = grid$between[i],
                    n_starts = 2, seed = 2, max_iter = 5000, tol = 1e-12)
    expect_true(fit$converged)
    expect_error(parameter_inference(fit), "ratings and group means computed",
                 class = "latents_unsupported_inference")
    expect_equal(fit$log_likelihood, review_cross_loglik(fit, data), tolerance = 1e-8)
    expect_equal(unname(rowSums(fit$group_posteriors)), rep(1, 90L), tolerance = 1e-12)
    expect_equal(unname(rowSums(fit$subject_posteriors)), rep(1, nrow(data)),
                 tolerance = 1e-12)
    if (grid$within[i] == "equal") {
      expect_equal(fit$profile_variances[1L, ], fit$profile_variances[2L, ])
    }
    if (grid$between[i] == "equal") {
      expect_equal(fit$group_mean_variances[1L, ], fit$group_mean_variances[2L, ])
    }
    row_weights <- unname(fit$sampling_weights)[match(data$id, unique(data$id))]
    prevalence <- if (fit$family == "full_cross_level") fit$composition else
      matrix(colSums(fit$subject_posteriors * row_weights) / sum(row_weights),
             2L, 2L, byrow = TRUE)
    block <- function(m, mode) if (mode == "equal") m[1L, ] else as.vector(m)
    pieces <- list(as.vector(fit$profile_means),
                   log(block(fit$profile_variances, grid$within[i])),
                   as.vector(fit$group_means),
                   log(block(fit$group_mean_variances, grid$between[i])),
                   log(prevalence[if (fit$family == "full_cross_level") 1:2 else 1L, 2L] /
                         prevalence[if (fit$family == "full_cross_level") 1:2 else 1L, 1L]),
                   log(fit$group_probabilities[2L] / fit$group_probabilities[1L]))
    ends <- cumsum(c(0L, lengths(pieces)))
    objective <- function(theta) {
      part <- function(k) theta[seq.int(ends[k] + 1L, ends[k + 1L])]
      expand <- function(values, mode) {
        if (mode == "equal") matrix(values, 2L, 2L, byrow = TRUE) else matrix(values, 2L)
      }
      moved <- fit
      moved$profile_means[] <- part(1L)
      moved$profile_variances <- exp(expand(part(2L), grid$within[i]))
      moved$group_means[] <- part(3L)
      moved$group_mean_variances <- exp(expand(part(4L), grid$between[i]))
      weight <- plogis(part(6L))
      moved$group_probabilities <- c(1 - weight, weight)
      share <- rep(plogis(part(5L)), length.out = 2L)
      review_cross_loglik(moved, data, cbind(1 - share, share))
    }
    expect_lt(max(abs(numDeriv::grad(objective, unlist(pieces)))), 2e-3)
  }))
})

test_that("weighted family covariance is the sandwich of weighted group scores", {
  fit <- multilpa(review_family_data(), c("y1", "y2"), "id", n_group_classes = 1,
                  family = "additive", weights = "w", seed = 1, n_starts = 1,
                  tol = 1e-12)
  inference <- parameter_inference(fit)
  bread <- solve(-attr(inference, "hessian"))
  expected <- bread %*% crossprod(attr(inference, "group_scores")) %*% bread
  expect_identical(attr(inference, "vcov_type"), "robust")
  expect_equal(attr(inference, "covariance_unconstrained"), expected, tolerance = 1e-10)
  expect_equal(vcov(fit, scale = "unconstrained"), expected, tolerance = 1e-10)
  expect_error(vcov(fit, type = "observed"), class = "latents_unsupported_weights")
  expect_error(parameter_inference(fit, method = "bootstrap"),
               class = "latents_unsupported_inference")
})

test_that("restricted cross-level weights preserve a similarly named indicator", {
  data <- review_family_data()
  ordinary <- multilpa(data, c("y1", "y2"), "id", n_profiles = 2,
                       n_group_classes = 2, family = "restricted_cross_level",
                       weights = "w", n_starts = 2, seed = 1, tol = 1e-10)
  names(data)[names(data) == "y1"] <- ".sampling_weight"
  named <- multilpa(data, c(".sampling_weight", "y2"), "id", n_profiles = 2,
                    n_group_classes = 2, family = "restricted_cross_level",
                    weights = "w", n_starts = 2, seed = 1, tol = 1e-10)
  expect_equal(named$log_likelihood, ordinary$log_likelihood, tolerance = 1e-10)
  expect_equal(unname(named$profile_means), unname(ordinary$profile_means),
               tolerance = 1e-10)
  expect_equal(named$log_likelihood, review_cross_loglik(named, data), tolerance = 1e-8)
})

test_that("cross-level counts refuse overflow before integer conversion", {
  data <- review_family_data()
  invisible(lapply(c("n_profiles", "n_group_classes", "n_starts", "max_iter"),
                    function(name) {
    args <- list(data = data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                 family = "full_cross_level")
    args[[name]] <- .Machine$integer.max + 1
    expect_error(do.call(multilpa, args), class = "latents_bad_argument")
  }))
})

test_that("family enumeration rejects missing variance restrictions explicitly", {
  expect_error(enumerate_classes(review_family_data(), c("y1", "y2"), "id",
                                 family = "additive", n_group_classes = 1,
                                 between_variance = c("equal", NA_character_)),
               class = "latents_bad_argument")
})
