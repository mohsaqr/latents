# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

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
  skip_on_cran()
  skip_without_reference("numDeriv")
  data <- review_family_data()
  vars <- c("y1", "y2")
  stats <- latents:::.additive_prepare(data, vars, "id", "w")
  blocks <- lapply(unique(data$id), function(j) as.matrix(data[data$id == j, vars]))
  cases <- list(c("additive", "varying"), c("additive", "equal"),
                c("dispersion", "equal"), c("additive_dispersion", "varying"),
                c("additive_dispersion", "equal"))
  sides <- lapply(cases, function(case) {
    label <- paste(case, collapse = "/")
    structure <- latents:::.additive_structure(case[1L], case[2L])
    shared <- function(value, mode) if (mode == "equal") value[c(1L, 1L), ] else value
    point <- list(means = shared(rbind(c(-1, -0.5), c(1, 0.5)), structure$means),
                  within = shared(rbind(c(0.6, 0.8), c(1.5, 1.8)), structure$within),
                  between = shared(rbind(c(0.3, 0.4), c(0.7, 0.5)), structure$between),
                  weights = c(0.7, 0.3))
    classes <- c("group_class_1", "group_class_2")
    theta <- latents:::.additive_pack(point, structure, classes, vars)
    # Build a full covariance for every group instead of using the engine's
    # group sufficient-statistic density or likelihood routine. The vector is
    # read back with the package's parameter map (.additive_unpack) and the
    # normalised sampling weights (.additive_prepare); the density is the
    # dense one of helper-additive-dense.R.
    objective <- function(theta) {
      p <- latents:::.additive_unpack(theta, 2L, 2L, structure)
      sum(stats$sampling_weights * vapply(blocks, function(x) {
        joint <- vapply(seq_len(2L), function(h) {
          additive_dense_density(x, p$means[h, ], p$between[h, ], p$within[h, ]) +
            log(p$weights[h])
        }, numeric(1))
        max(joint) + log(sum(exp(joint - max(joint))))
      }, numeric(1)))
    }
    scores <- latents:::.additive_group_scores(stats, point, structure, classes, vars)
    information <- latents:::.additive_information(stats, theta, structure, classes, vars, 1e-4)
    cells <- expand.grid(row = names(theta), column = names(theta), stringsAsFactors = FALSE)
    list(log_likelihood = c(latents:::.additive_expectation(stats, point)$log_likelihood,
                            objective(theta)),
         score = list(
           actual = data.frame(family = label, parameter = names(theta),
                               score = unname(colSums(scores))),
           expected = data.frame(family = label, parameter = names(theta),
                                 score = numDeriv::grad(objective, theta))),
         information = list(
           actual = data.frame(family = label, cells, information = as.vector(unname(information))),
           expected = data.frame(family = label, cells,
                                 information = as.vector(-numDeriv::hessian(objective, theta)))))
  })
  labels <- vapply(cases, paste, character(1), collapse = "/")
  stack <- function(part, side) do.call(rbind, lapply(sides, \(s) s[[part]][[side]]))
  note <- "Package side: the internal weighted E-step, analytic scores and observed information (.additive_expectation, .additive_group_scores, .additive_information) behind multilpa(family = ..., weights = ) and its vcov()."
  expect_reference(setNames(vapply(sides, \(s) s$log_likelihood[[1L]], numeric(1)), labels),
                   setNames(vapply(sides, \(s) s$log_likelihood[[2L]], numeric(1)), labels),
                   reference = "independent: weighted dense Cholesky likelihood (helper-additive-dense.R)",
                   covers = "multilpa()", quantity = "weighted log likelihood, by family",
                   tolerance = 1e-10, measure = "relative", note = note)
  # numDeriv's error is absolute near zero; a relative measure would judge
  # the rounding noise of a near-zero score. The original tolerance, 1e-6.
  expect_reference(stack("score", "actual"), stack("score", "expected"), reference = "numDeriv",
                   covers = "multilpa()", quantity = "weighted score off the optimum, by family",
                   tolerance = 1e-6, measure = "scaled", key = c("family", "parameter"),
                   note = paste(note, "Reference: numDeriv::grad() of the weighted dense likelihood."))
  # 1e-4, the original tolerance: the package's information is a central
  # difference of the analytic scores with step 1e-4, the reference a
  # second numerical derivative of the dense likelihood. Scaled, because
  # entries run down to zero.
  expect_reference(stack("information", "actual"), stack("information", "expected"),
                   reference = "numDeriv", covers = c("multilpa()", "vcov.multilpa_additive()"),
                   quantity = "weighted observed information off the optimum, by family",
                   tolerance = 1e-4, measure = "scaled", key = c("family", "row", "column"),
                   note = paste(note, "Reference: minus numDeriv::hessian() of the weighted dense likelihood."))
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
  skip_on_cran()
  skip_without_reference("numDeriv")
  data <- review_family_data()
  grid <- expand.grid(family = c("restricted_cross_level", "full_cross_level"),
                       within = c("equal", "varying"),
                       between = c("equal", "varying"), stringsAsFactors = FALSE)
  sides <- lapply(seq_len(nrow(grid)), function(i) {
    label <- paste(grid$family[i], grid$within[i], grid$between[i], sep = "/")
    fit <- multilpa(data, c("y1", "y2"), "id", n_profiles = 2,
                    n_group_classes = 2, family = grid$family[i], weights = "w",
                    variance_model = grid$within[i], between_variance = grid$between[i],
                    n_starts = 2, seed = 2, max_iter = 5000, tol = 1e-12)
    expect_true(fit$converged)
    expect_error(parameter_inference(fit), "ratings and group means computed",
                 class = "latents_unsupported_inference")
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
    gradient <- numDeriv::grad(objective, unlist(pieces))
    list(log_likelihood = c(as.numeric(logLik(fit)), review_cross_loglik(fit, data)),
         gradient = data.frame(model = label, parameter = seq_along(gradient), gradient = gradient))
  })
  labels <- paste(grid$family, grid$within, grid$between, sep = "/")
  # The original tolerance, 1e-8 relative on a sum of 90 weighted group terms.
  expect_reference(setNames(vapply(sides, \(s) s$log_likelihood[[1L]], numeric(1)), labels),
                   setNames(vapply(sides, \(s) s$log_likelihood[[2L]], numeric(1)), labels),
                   reference = "independent: weighted manifest cross-level likelihood, stats::dnorm()",
                   covers = c("multilpa()", "logLik.multilpa_cross_level()"),
                   quantity = "weighted log likelihood, by family / within / between variance",
                   tolerance = 1e-8, measure = "relative")
  gradients <- do.call(rbind, lapply(sides, `[[`, "gradient"))
  # A maximum has a zero gradient: the package side is that zero, the
  # reference the numerical gradient of the raw-formula likelihood at the
  # fitted point. The original bound, max |gradient| < 2e-3, absolute.
  expect_reference(transform(gradients, gradient = 0), gradients, reference = "numDeriv",
                   covers = "multilpa()",
                   quantity = "gradient of the raw-formula likelihood at the fitted maximum",
                   tolerance = 2e-3, measure = "absolute", key = c("model", "parameter"),
                   note = "Reference: numDeriv::grad() of review_cross_loglik(), written in this file from univariate dnorm() densities. Parameters: profile means, log profile variances, group means, log group-mean variances, prevalence logits, class logit.")
})
