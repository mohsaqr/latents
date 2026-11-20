# Comparisons with independent references moved from tests/testthat/test-review-1026-additive.R,
# so that they run in CI's equivalence job rather than under R CMD check.

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

# Muffle only the warnings and notices the package raises on purpose
# (quietly(), tests/testthat/helper-expected-warnings.R); any other warning
# still reaches testthat.
review_quiet <- function(expr) quietly(expr)

review_structures <- function() {
  list(additive = latents:::.additive_structure("additive", "varying"),
       additive_equal = latents:::.additive_structure("additive", "equal"),
       dispersion = latents:::.additive_structure("dispersion"),
       additive_dispersion = latents:::.additive_structure("additive_dispersion", "varying"),
       additive_dispersion_equal = latents:::.additive_structure("additive_dispersion", "equal"))
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
  skip_without_reference("numDeriv")
  data <- review_additive_draw(1L)
  vars <- c("y1", "y2")
  stats <- latents:::.additive_prepare(data, vars, "g")
  class_names <- paste0("group_class_", 1:2)
  sides <- lapply(names(review_structures()), function(name) {
    structure <- review_structures()[[name]]
    parameters <- latents:::.additive_check_parameters(review_parameters(), 2L)
    parameters$means <- latents:::.additive_pool(parameters$means, matrix(1, 2, 2), structure$means)
    parameters$between <- latents:::.additive_pool(parameters$between, matrix(1, 2, 2),
                                                   structure$between)
    parameters$within <- latents:::.additive_pool(parameters$within, matrix(1, 2, 2),
                                                  structure$within)
    theta <- latents:::.additive_pack(parameters, structure, class_names, vars)
    expect_length(theta, latents:::.additive_count_parameters(2L, 2L, structure))
    expect_equal(latents:::.additive_pack(latents:::.additive_unpack(theta, 2L, 2L, structure),
                                          structure, class_names, vars), theta,
                 tolerance = 1e-14, info = name)
    log_likelihood <- function(point) {
      latents:::.additive_expectation(stats, latents:::.additive_unpack(point, 2L, 2L, structure))$log_likelihood
    }
    analytic <- colSums(latents:::.additive_group_scores(
      stats, latents:::.additive_unpack(theta, 2L, 2L, structure), structure, class_names, vars))
    information <- latents:::.additive_information(stats, theta, structure, class_names, vars, 1e-4)
    hessian <- -numDeriv::hessian(log_likelihood, unname(theta))
    cells <- expand.grid(row = names(theta), column = names(theta), stringsAsFactors = FALSE)
    list(score = list(
           actual = data.frame(family = name, parameter = names(theta), score = unname(analytic)),
           expected = data.frame(family = name, parameter = names(theta),
                                 score = numDeriv::grad(log_likelihood, unname(theta)))),
         information = list(
           actual = data.frame(family = name, cells, information = as.vector(unname(information))),
           expected = data.frame(family = name, cells, information = as.vector(hessian))))
  })
  stack <- function(part, side) do.call(rbind, lapply(sides, \(s) s[[part]][[side]]))
  # numDeriv's error is absolute near zero; a relative measure would judge
  # the rounding noise of a near-zero score or information entry. The
  # original tolerance, 1e-6, for both.
  expect_reference(stack("score", "actual"), stack("score", "expected"), reference = "numDeriv",
                   covers = "multilpa()", quantity = "score off the optimum, every additive family",
                   tolerance = 1e-6, measure = "scaled", key = c("family", "parameter"),
                   note = "Package side: the internal analytic group scores (.additive_group_scores) behind the standard errors of multilpa(family = 'additive' / 'dispersion' / 'additive_dispersion'); reference: numDeriv::grad() of the package's likelihood.")
  expect_reference(stack("information", "actual"), stack("information", "expected"),
                   reference = "numDeriv", covers = c("multilpa()", "vcov.multilpa_additive()"),
                   quantity = "observed information off the optimum, every additive family",
                   tolerance = 1e-6, measure = "scaled", key = c("family", "row", "column"),
                   note = "Package side: the internal observed information (.additive_information) behind vcov() of the additive families; reference: minus numDeriv::hessian() of the package's likelihood.")
})

test_that("standard errors are translation invariant and scale equivariant", {
  skip_without_reference("numDeriv")
  data <- review_additive_draw(1L)
  fit <- review_quiet(multilpa(data, c("y1", "y2"), "g", n_group_classes = 2L,
                               family = "additive_dispersion", n_starts = 3L,
                               seed = 1L, tol = 1e-10))
  expect_true(fit$converged)
  standard_errors <- function(object, type) {
    latents:::.additive_parameter_table(object, latents:::.additive_inference(
      object, type, warn_weak = FALSE))$standard_error
  }
  kind <- latents:::.additive_parameter_table(fit, NULL)$parameter
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
    latents:::.additive_expectation(stats, latents:::.additive_unpack(point, 2L, 2L, structure))$log_likelihood
  }
  reference <- solve(-numDeriv::hessian(log_likelihood, unname(theta)))
  observed <- latents:::.additive_inference(fit, "observed", warn_weak = FALSE)$vcov
  # Standard errors are positive and of order 0.1 to 1: relative, the
  # original tolerance 1e-6 (the package's information is a central
  # difference of the analytic scores with step 1e-4).
  expect_reference(setNames(sqrt(diag(observed)), names(theta)),
                   setNames(sqrt(diag(reference)), names(theta)), reference = "numDeriv",
                   covers = c("multilpa()", "coef.multilpa_additive()", "vcov.multilpa_additive()"),
                   quantity = "observed standard errors, unconstrained scale, additive dispersion",
                   tolerance = 1e-6, measure = "relative",
                   note = "Package side: the internal observed covariance (.additive_inference) that vcov(fit, type = 'observed', scale = 'unconstrained') returns; reference: the inverse of minus numDeriv::hessian() of the package's likelihood at coef(fit, scale = 'unconstrained').")
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
    latents:::.additive_expectation(stats, latents:::.additive_unpack(point, 2L, 2L, structure))$log_likelihood
  }
  polished <- stats::optim(coef(fit, scale = "unconstrained"), log_likelihood,
                           method = "BFGS",
                           control = list(fnscale = -1, reltol = 1e-14, maxit = 500L))
  # BFGS started at the fit only climbs, so the original one-sided bound
  # (polished - fitted < 1e-8) is this absolute difference.
  expect_reference(as.numeric(logLik(fit)), polished$value, reference = "stats",
                   covers = c("multilpa()", "logLik.multilpa_additive()", "coef.multilpa_additive()"),
                   quantity = "maximum of the log likelihood, dispersion family",
                   tolerance = 1e-8, measure = "absolute",
                   note = "Reference: stats::optim(method = 'BFGS') maximising the package's likelihood from coef(fit, scale = 'unconstrained').")

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
  families <- c("full_cross_level", "restricted_cross_level")
  fits <- setNames(lapply(families, function(family) {
    review_quiet(multilpa(data, vars, "g", n_profiles = 2L, n_group_classes = 2L,
                          family = family, n_starts = 2L, seed = 1L, tol = 1e-10,
                          max_iter = 3000L))
  }), families)
  expect_reference(vapply(fits, \(fit) as.numeric(logLik(fit)), numeric(1)),
                   vapply(fits, brute, numeric(1)),
                   reference = "independent: brute-force cross-level likelihood, stats::dnorm()",
                   covers = c("multilpa()", "logLik.multilpa_cross_level()"),
                   quantity = "log likelihood, by cross-level family",
                   tolerance = 1e-10, measure = "relative")
  invisible(lapply(fits, function(fit) {
    expect_equal(fit$bic, -2 * fit$log_likelihood + fit$n_parameters * log(length(ids)))
    expect_equal(fit$bic_individual,
                 -2 * fit$log_likelihood + fit$n_parameters * log(nrow(data)))
  }))
  # 2 x 2 profile means and variances, prevalences (per class when full),
  # 2 x 2 group means and variances, one class logit.
  expect_reference(vapply(fits, \(fit) as.numeric(fit$n_parameters), numeric(1)),
                   c(full_cross_level = 19, restricted_cross_level = 18),
                   reference = "independent: parameter count from the model definition",
                   covers = "multilpa()", quantity = "number of parameters, by cross-level family",
                   tolerance = 0, measure = "absolute")
})

test_that("covariance-structure charts have mclust's parameter counts", {
  skip_without_reference("mclust")
  codes <- setdiff(latents:::.multilpa_structures(), latents:::.multilpa_separable_structures())
  widths <- vapply(codes, function(code) latents:::.multilpa_chart(code, 3L, 3L)$width, integer(1))
  counts <- vapply(codes, function(code) as.integer(mclust::nVarParams(code, d = 3L, G = 3L)),
                   integer(1))
  # Integer counts: exact.
  expect_reference(widths, counts, reference = "mclust",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "covariance parameters, 3 profiles, 3 indicators, by structure",
                   tolerance = 0, measure = "absolute",
                   note = "Package side: the width of the internal covariance chart (.multilpa_chart) on which parameter_inference() of a covariance-structure fit is computed; reference: mclust::nVarParams().")
})

test_that("covariance-structure charts have analytic Jacobians and positive definite maps", {
  skip_without_reference("numDeriv")
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old_seed, globalenv()), add = TRUE)
  set.seed(7)
  codes <- setdiff(latents:::.multilpa_structures(), latents:::.multilpa_separable_structures())
  sides <- lapply(codes, function(code) {
    chart <- latents:::.multilpa_chart(code, 3L, 3L)
    anchor <- if (identical(chart$form, "rotation")) {
      lapply(seq_len(chart$rotation_rows), function(row) {
        qr.Q(qr(matrix(stats::rnorm(9L), 3L)))
      })
    }
    phi <- stats::rnorm(chart$width, sd = 0.3)
    mapped <- latents:::.multilpa_chart_map(phi, chart, anchor, derivatives = TRUE)
    numerical <- numDeriv::jacobian(function(point) {
      unlist(latents:::.multilpa_chart_map(point, chart, anchor)$sigma)
    }, phi)
    # Every profile's covariance is symmetric positive definite.
    expect_true(all(vapply(mapped$sigma, function(sigma) {
      isTRUE(all.equal(sigma, t(sigma))) &&
        min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values) > 0
    }, logical(1))), info = code)
    cells <- expand.grid(entry = seq_len(nrow(numerical)), coordinate = seq_len(ncol(numerical)))
    list(actual = data.frame(structure = code, cells,
                             derivative = as.vector(do.call(rbind, mapped$jacobian))),
         expected = data.frame(structure = code, cells, derivative = as.vector(numerical)))
  })
  # Jacobian entries are often exactly zero, where numDeriv's error is
  # absolute: scaled. The original tolerance, 1e-7.
  expect_reference(do.call(rbind, lapply(sides, `[[`, "actual")),
                   do.call(rbind, lapply(sides, `[[`, "expected")), reference = "numDeriv",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "Jacobian of the covariance chart, by structure",
                   tolerance = 1e-7, measure = "scaled", key = c("structure", "entry", "coordinate"),
                   note = "Package side: the internal analytic chart Jacobian (.multilpa_chart_map) used by parameter_inference() of a covariance-structure fit; reference: numDeriv::jacobian() of the package's chart map. entry runs over the stacked profile covariances, coordinate over the chart.")
})

test_that("the covariance gradient is the derivative of the weighted Gaussian log density", {
  skip_without_reference("numDeriv")
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
  gradient <- latents:::.multilpa_profile_gradients(parameters, expectation, x)
  log_density <- function(entries) {
    sigma <- matrix(entries, 2L)
    residuals <- sweep(x, 2L, parameters$means[1L, ])
    -0.5 * sum(weights * (2 * log(2 * pi) +
                            as.numeric(determinant(sigma)$modulus) +
                            rowSums((residuals %*% t(solve(sigma))) * residuals)))
  }
  mean_density <- function(mean) {
    residuals <- sweep(x, 2L, mean)
    -0.5 * sum(weights * rowSums((residuals %*% solve(covariance)) * residuals))
  }
  note <- "Package side: the internal profile gradients (.multilpa_profile_gradients) behind the scores of full-covariance fits; reference: numDeriv::grad() of the weighted Gaussian log density written in this test."
  # numDeriv's error is absolute near zero: scaled. The original tolerance, 1e-8.
  expect_reference(as.vector(gradient$covariance[[1L]]),
                   numDeriv::grad(log_density, as.vector(covariance)), reference = "numDeriv",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "gradient of the weighted log density in the covariance entries",
                   tolerance = 1e-8, measure = "scaled", note = note)
  expect_reference(as.vector(gradient$means),
                   numDeriv::grad(mean_density, c(0.2, -0.1)), reference = "numDeriv",
                   covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "gradient of the weighted log density in the means",
                   tolerance = 1e-8, measure = "scaled", note = note)
})
