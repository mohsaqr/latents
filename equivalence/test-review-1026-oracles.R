# Review 2026-10-26: comparisons of latents against external packages that
# the shipped tests do not declare (depmixS4, mvtnorm, nlme, poLCA, Rmpfr,
# sandwich). Moved here from tests/testthat/test-review-1026-*.R so the
# package gains no Suggests; run with `Rscript equivalence/run.R review-1026`.

# helper from tests/testthat/test-review-1026-additive.R
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

# helper from tests/testthat/test-review-1026-growth.R
review_growth_data <- function(n = 70L, seed = 21L, classes = 2L) {
  withr::with_seed(seed, {
    times <- 0:4
    class <- sample(seq_len(classes), n, TRUE)
    effects <- t(t(chol(matrix(c(1, 0.2, 0.2, 0.15), 2L))) %*%
                   matrix(stats::rnorm(2L * n), 2L))
    out <- data.frame(id = rep(seq_len(n), each = length(times)),
                      time = rep(times, n), x = stats::rnorm(n * length(times)))
    person_class <- class[out$id]
    out$y <- c(2, 6)[person_class] + c(0.5, -0.4)[person_class] * out$time +
      0.3 * out$x + effects[out$id, 1L] + effects[out$id, 2L] * out$time +
      stats::rnorm(nrow(out), 0, 0.7)
    out[stats::runif(nrow(out)) > 0.2, ]
  })
}

# helper from tests/testthat/test-review-1026-growth.R
review_growth_fit <- function(data, n_classes = 2, ...) {
  mixture_regression(y ~ time + x, data, n_classes = n_classes, id = "id",
                     class_level = "group", common = ~ x, random = ~ 1 + time,
                     n_starts = 2, seed = 1, vcov_type = "none", ...)
}

# helper from tests/testthat/test-review-1026-lta.R
review_lta_data <- function(n = 30L, n_occasions = 3L, seed = 1L, stay = 0.75) {
  withr::with_seed(seed, {
    do.call(rbind, lapply(seq_len(n), function(i) {
      states <- Reduce(function(previous, t) {
        if (stats::runif(1L) < stay) previous else 3L - previous
      }, seq_len(n_occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
      data.frame(id = i, time = seq_len(n_occasions), z = stats::rnorm(1L),
                 w = stats::rnorm(n_occasions),
                 y1 = stats::rnorm(n_occasions, c(0, 2)[states], 0.7),
                 y2 = stats::rnorm(n_occasions, c(0, 1.5)[states], 0.7))
    }))
  })
}

# helper from tests/testthat/test-review-1026-additive.R
review_parameters <- function() {
  list(means = rbind(c(0.1, 0.9), c(2.4, -0.8)),
       between = rbind(c(0.4, 0.7), c(1.1, 0.5)),
       within = rbind(c(1.1, 1.4), c(0.7, 1.9)),
       weights = c(0.55, 0.45))
}

# helper from tests/testthat/test-review-1026-categorical.R
review_two_class_items <- function(n, seed) {
  withr::local_seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(.4, .6))
  probabilities <- list(rbind(c(.8, .2), c(.2, .8)),
                        rbind(c(.7, .2, .1), c(.1, .3, .6)),
                        rbind(c(.6, .4), c(.3, .7)),
                        rbind(c(.5, .3, .2), c(.2, .2, .6)))
  items <- lapply(probabilities, function(p) vapply(class, function(k)
    sample(seq_len(ncol(p)), 1L, prob = p[k, ]), integer(1)))
  stats::setNames(as.data.frame(items), paste0("i", 1:4))
}

# from tests/testthat/test-review-1026-lta.R
test_that("lta() agrees with depmixS4 at the same parameters", {
  skip_if_not_installed("depmixS4")
  data <- review_lta_data(n = 30L, n_occasions = 4L, seed = 23L)
  data <- subset(data, !(id %in% 1:5 & time == 4))
  data <- data[order(data$id, data$time), ]
  fit <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
             n_starts = 2, seed = 1)
  model <- depmixS4::depmix(list(y1 ~ 1, y2 ~ 1), data = data, nstates = 2,
                            family = list(stats::gaussian(), stats::gaussian()),
                            ntimes = as.vector(table(data$id)))
  response <- function(k) c(fit$means[k, 1L], sqrt(fit$variances[k, 1L]),
                            fit$means[k, 2L], sqrt(fit$variances[k, 2L]))
  model <- depmixS4::setpars(model, c(fit$initial_probabilities[1L, ],
                                      as.vector(t(fit$transition_probabilities[, , 1L])),
                                      response(1L), response(2L)))
  passes <- depmixS4::forwardbackward(model)
  expect_equal(passes$logLike, fit$log_likelihood, tolerance = 1e-10)
  expect_equal(unname(passes$gamma),
               unname(fit$subject_posteriors), tolerance = 1e-10)
})

# from tests/testthat/test-review-1026-growth.R
test_that("one class reproduces the linear mixed model of nlme::lme (ML)", {
  skip_on_cran()
  skip_if_not_installed("nlme")
  data <- review_growth_data(n = 80L, seed = 11L, classes = 1L)
  fit <- review_growth_fit(data, n_classes = 1)
  reference <- nlme::lme(y ~ time + x, random = ~ 1 + time | id, data = data,
                         method = "ML",
                         control = nlme::lmeControl(maxIter = 500L, msMaxIter = 500L,
                                                    tolerance = 1e-10, msTol = 1e-12))
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)), tolerance = 1e-7)
  expect_equal(unname(c(fit$params$beta[, 1L], fit$params$common)),
               unname(nlme::fixef(reference)), tolerance = 1e-4)
  expect_equal(unname(fit$params$random_covariance[[1L]]),
               unname(matrix(as.numeric(nlme::getVarCov(reference)), 2L)), tolerance = 1e-3)
  # Named by its class since 0.98.0.
  expect_equal(unname(fit$params$sigma2), reference$sigma^2, tolerance = 1e-4)
  # BIC counts persons (the independent units), not observations.
  fit_table <- get_results(fit, "fit")
  expect_equal(fit_table$bic, -2 * fit$log_likelihood +
                 log(length(unique(data$id))) * fit$n_parameters)
})

# from tests/testthat/test-review-1026-growth.R
test_that("an outcome in small units fits like the same outcome in ordinary units", {
  skip_on_cran()
  skip_if_not_installed("nlme")
  data <- review_growth_data(n = 80L, seed = 11L, classes = 1L)
  data$y <- data$y * 1e-5
  # Before the relative jitter, EM stopped with "likelihood decreased".
  fit <- review_growth_fit(data, n_classes = 1)
  reference <- nlme::lme(y ~ time + x, random = ~ 1 + time | id, data = data, method = "ML")
  # The quasi-Newton finish is not scale-free (reported), so only 0.05 here.
  expect_lt(abs(fit$log_likelihood - as.numeric(stats::logLik(reference))), 0.05)
})

# from tests/testthat/test-review-1026-additive.R
test_that("sufficient-statistic density equals the dense Gaussian density of the ratings", {
  skip_if_not_installed("mvtnorm")
  data <- review_additive_draw(1L, n_groups = 25L)
  stats <- .additive_prepare(data, c("y1", "y2"), "g")
  parameters <- review_parameters()
  dense <- vapply(1:2, function(h) {
    vapply(seq_len(stats$n_groups), function(j) {
      ratings <- as.matrix(subset(data, g == stats$ids[j], c("y1", "y2")))
      n <- nrow(ratings)
      sum(vapply(1:2, function(r) {
        covariance <- parameters$within[h, r] * diag(n) +
          parameters$between[h, r] * matrix(1, n, n)
        mvtnorm::dmvnorm(ratings[, r], rep(parameters$means[h, r], n), covariance,
                         log = TRUE)
      }, numeric(1)))
    }, numeric(1))
  }, numeric(stats$n_groups))
  expect_equal(unname(.additive_log_density(stats, parameters)), dense,
               tolerance = 1e-12)
  mixture <- sum(log(exp(dense) %*% parameters$weights))
  expect_equal(.additive_expectation(stats, parameters)$log_likelihood, mixture,
               tolerance = 1e-12)
})

# from tests/testthat/test-review-1026-kernel.R
test_that("FIML densities and conditional moments are the multivariate normal ones", {
  skip_if_not_installed("mvtnorm")
  set.seed(21)
  n <- 40L
  d <- 3L
  x <- matrix(stats::rnorm(n * d), n, d)
  x[sample(n * d, 30L)] <- NA
  x[5L, ] <- NA
  covariances <- array(c(2, .6, .3, .6, 1, .2, .3, .2, 1.5,
                         1, -.4, 0, -.4, 2, .5, 0, .5, 1), c(d, d, 2L))
  parameters <- list(means = matrix(c(0, 1, .5, -.5, -1, .3), 2L, d),
                     variances = t(apply(covariances, 3L, diag)),
                     covariances = covariances)
  moments <- .multilpa_gaussian_moments(x, parameters)
  invisible(lapply(seq_len(2L), function(profile) {
    sigma <- covariances[, , profile]
    mu <- parameters$means[profile, ]
    invisible(lapply(seq_len(n), function(i) {
      observed <- which(!is.na(x[i, ]))
      missing <- which(is.na(x[i, ]))
      reference <- if (length(observed)) {
        mvtnorm::dmvnorm(x[i, observed], mu[observed],
                         matrix(sigma[observed, observed], length(observed)), log = TRUE)
      } else 0
      expect_equal(moments$log_density[i, profile], reference, tolerance = 1e-12)
      if (length(missing)) {
        regression <- if (length(observed)) {
          sigma[missing, observed, drop = FALSE] %*%
            solve(sigma[observed, observed, drop = FALSE])
        } else matrix(0, length(missing), 0L)
        mean_given <- mu[missing] + if (length(observed)) {
          as.vector(regression %*% (x[i, observed] - mu[observed]))
        } else 0
        covariance_given <- sigma[missing, missing, drop = FALSE] -
          regression %*% sigma[observed, missing, drop = FALSE]
        adjustment <- Filter(function(block) i %in% block$rows,
                             moments$moments[[profile]]$adjustments)[[1L]]$covariance
        expect_equal(moments$moments[[profile]]$expected[i, missing], mean_given,
                     tolerance = 1e-12)
        expect_equal(adjustment[missing, missing, drop = FALSE], covariance_given,
                     tolerance = 1e-12, ignore_attr = TRUE)
      }
    }))
  }))
})

# from tests/testthat/test-review-1026-categorical.R
test_that("Poisson offset stays within 1e-7 of the exact log ratio at counts near 1e17", {
  skip_if_not_installed("Rmpfr")
  mpfr <- Rmpfr::mpfr
  counts <- c(1e17, 1e17 + 16000, 3)
  means <- c(1e17, 1e17 * (1 + 1e-9), 3)
  extra <- .latents_set_count_model(
    .latents_prepare_extra(data.frame(y = counts), character(), "y"))
  offset <- .latents_extra_log_density(extra, list(count_means = matrix(means, 3)),
                                       3, 3, offset = TRUE)
  relative <- offset$log_density - offset$log_density[, 1]
  mu <- mpfr(means, 200)
  exact <- t(vapply(counts, function(value) {
    y <- mpfr(value, 200)
    as.numeric(y * log(mu) - mu - (y * log(mu[1]) - mu[1]))
  }, numeric(3)))
  expect_lt(max(abs(relative[1:2, 2] - exact[1:2, 2])), 1e-7)
  expect_equal(relative[, 3] / exact[, 3], rep(1, 3), tolerance = 1e-12)
})

# from tests/testthat/test-review-1026-categorical.R
test_that("negative-binomial dispersion derivatives match exact Rmpfr values", {
  skip_if_not_installed("Rmpfr")
  mpfr <- Rmpfr::mpfr
  exact <- function(y, mu, alpha) {
    a <- mpfr(alpha, 300)
    m <- mpfr(mu, 300)
    r <- 1 / a
    j <- mpfr(seq_len(y) - 1, 300)
    d <- (if (y > 0) sum(1 / (r + j)) else 0) + log(r / (r + m)) + (m - y) / (r + m)
    d_prime <- (if (y > 0) -sum(1 / (r + j)^2) else 0) + 1 / r - 1 / (r + m) -
      (m - y) / (r + m)^2
    c(as.numeric(-r * d), as.numeric(r * d + r^2 * d_prime))
  }
  # Covers the regular, small-alpha series and large-size asymptotic branches.
  grid <- expand.grid(y = c(0, 1, 5, 30, 200), mu = c(.5, 3, 40, 300),
                      alpha = c(1e-9, 1e-6, 1e-4, 2e-3, .05, 1, 5))
  errors <- vapply(seq_len(nrow(grid)), function(i) {
    point <- grid[i, ]
    got <- .latents_negative_binomial_dispersion_derivatives(point$y, point$mu, point$alpha)
    want <- exact(point$y, point$mu, point$alpha)
    c(abs(got$score - want[1]) / max(abs(want[1]), 1e-300),
      abs(got$curvature - want[2]) / max(abs(want[2]), 1e-300))
  }, numeric(2))
  expect_lt(max(errors), 1e-8)
})

# from tests/testthat/test-review-1026-categorical.R
test_that("latent class fit agrees with poLCA and is an EM fixed point", {
  skip_if_not_installed("poLCA")
  data <- review_two_class_items(400, 3)
  fit <- lca(data, vars = names(data), n_classes = 2, n_starts = 5, seed = 1,
             tol = 1e-12, max_iter = 20000)
  polca <- poLCA::poLCA
  reference <- polca(cbind(i1, i2, i3, i4) ~ 1, data, nclass = 2, nrep = 5,
                     verbose = FALSE, maxiter = 20000, tol = 1e-14)
  expect_equal(fit$log_likelihood, reference$llik, tolerance = 1e-9)
  expect_identical(as.integer(fit$n_parameters), as.integer(reference$npar))
  # Hand E-step from the fitted parameters, then the closed-form M-step.
  responses <- fit$response_probabilities
  joint <- vapply(1:2, function(k) log(fit$profile_probabilities[k]) +
    Reduce(`+`, lapply(names(data), function(item) log(responses[[item]][k, data[[item]]]))),
    numeric(nrow(data)))
  posterior <- exp(joint - apply(joint, 1L, max))
  posterior <- posterior / rowSums(posterior)
  expect_equal(unname(fit$subject_posteriors), unname(posterior), tolerance = 1e-12)
  updated <- .multilpa_categorical_maximize(as.matrix(data), posterior,
                                            c(2L, 3L, 2L, 3L), 1e-10)
  # One further EM step moves a fit stopped at tol = 1e-12 by about 1e-8.
  expect_equal(lapply(updated, unname), unname(lapply(responses, unname)),
               tolerance = 1e-6)
  expect_equal(unname(rowSums(responses$i2)), c(1, 1), tolerance = 1e-14)
})

# from tests/testthat/test-review-1026-inference.R
test_that("the CR1 sandwich composition equals sandwich::vcovCL(type = 'HC0')", {
  skip_if_not_installed("sandwich")
  set.seed(2)
  n_groups <- 12L
  data <- data.frame(cluster = rep(seq_len(n_groups), each = 8), x = stats::rnorm(96))
  data$y <- stats::rbinom(96, 1, stats::plogis(0.3 + 0.8 * data$x +
                                                 rep(stats::rnorm(n_groups, 0, 0.5), each = 8)))
  model <- stats::glm(y ~ x, stats::binomial, data)
  meat <- crossprod(rowsum(sandwich::estfun(model), data$cluster))
  ours <- .inference_sandwich(unname(stats::vcov(model)), meat, n_groups)
  reference <- sandwich::vcovCL(model, cluster = ~cluster, type = "HC0", cadjust = TRUE)
  expect_equal(ours, unname(reference), tolerance = 1e-12)
  expect_identical(.inference_cluster_factor(n_groups), n_groups / (n_groups - 1))
  expect_error(.inference_cluster_factor(1))
})
