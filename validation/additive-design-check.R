# Numerical checks of the proposed additive model, not a package estimator.
# Run from the project root: Rscript validation/additive-design-check.R
# Uses only base R. No reports or reference outputs are overwritten.

#' Evaluate the proposed diagonal additive model
#' @param observations List of complete group-by-indicator matrices.
#' @param parameters Group-class means, between variances, within variances, weights.
#' @return Marginal likelihood and conditional group-intercept moments.
additive_design_evaluate <- function(observations, parameters) {
  stopifnot(is.list(observations), length(observations) > 0L,
            is.list(parameters), is.matrix(parameters$means),
            identical(dim(parameters$means), dim(parameters$between)),
            all(parameters$between >= 0), all(parameters$within > 0),
            all(parameters$weights > 0),
            abs(sum(parameters$weights) - 1) < 1e-12,
            length(parameters$within) == ncol(parameters$means),
            length(parameters$weights) == nrow(parameters$means),
            all(vapply(observations, function(x) {
              is.matrix(x) && is.numeric(x) && all(is.finite(x)) &&
                nrow(x) > 0L && ncol(x) == ncol(parameters$means)
            }, logical(1))))
  groups <- lapply(observations, function(x) {
    n <- nrow(x)
    average <- colMeans(x)
    scatter <- colSums(sweep(x, 2L, average, "-")^2)
    components <- lapply(seq_len(nrow(parameters$means)), function(h) {
      w <- parameters$within
      t <- parameters$between[h, ]
      mu <- parameters$means[h, ]
      denominator <- w + n * t
      list(log_density = -0.5 * sum(n * log(2 * pi) +
             (n - 1) * log(w) + log(denominator) + scatter / w +
             n * (average - mu)^2 / denominator),
           mean = mu + n * t / denominator * (average - mu),
           variance = t * w / denominator)
    })
    log_joint <- vapply(components, function(z) z$log_density, numeric(1)) +
      log(parameters$weights)
    offset <- max(log_joint)
    mass <- exp(log_joint - offset)
    list(log_likelihood = offset + log(sum(mass)),
         posterior = mass / sum(mass), components = components,
         n = n, average = average, scatter = scatter)
  })
  list(log_likelihood = sum(vapply(groups, function(z) z$log_likelihood,
                                 numeric(1))), groups = groups)
}

#' Evaluate a group density through its full covariance matrix
#' @param x Complete numeric group matrix.
#' @param mu Indicator means.
#' @param between Between-group indicator variances.
#' @param within Within-group indicator variances.
#' @return Component log density from direct Cholesky calculations.
additive_design_dense <- function(x, mu, between, within) {
  stopifnot(is.matrix(x), all(is.finite(x)), length(mu) == ncol(x),
            length(between) == ncol(x), length(within) == ncol(x),
            all(between >= 0), all(within > 0))
  sum(vapply(seq_len(ncol(x)), function(r) {
    covariance <- diag(within[r], nrow(x)) +
      matrix(between[r], nrow(x), nrow(x))
    residual <- x[, r] - mu[r]
    factor <- chol(covariance)
    -0.5 * (nrow(x) * log(2 * pi) + 2 * sum(log(diag(factor))) +
              sum(residual * solve(covariance, residual)))
  }, numeric(1)))
}

#' Take one interior EM step using conditional intercept moments
#' @param evaluated Result of additive_design_evaluate().
#' @return Proposed parameters with class-specific between variances.
additive_design_step <- function(evaluated) {
  stopifnot(is.list(evaluated), length(evaluated$groups) > 0L)
  groups <- evaluated$groups
  h_count <- length(groups[[1L]]$posterior)
  membership <- colSums(do.call(rbind, lapply(groups, function(z) z$posterior)))
  means <- do.call(rbind, lapply(seq_len(h_count), function(h) {
    Reduce(`+`, lapply(groups, function(z) {
      z$posterior[h] * z$components[[h]]$mean
    })) / membership[h]
  }))
  between <- do.call(rbind, lapply(seq_len(h_count), function(h) {
    Reduce(`+`, lapply(groups, function(z) {
      z$posterior[h] * (z$components[[h]]$variance +
                         (z$components[[h]]$mean - means[h, ])^2)
    })) / membership[h]
  }))
  within <- Reduce(`+`, lapply(groups, function(z) {
    Reduce(`+`, lapply(seq_len(h_count), function(h) {
      z$posterior[h] * (z$scatter + z$n *
        ((z$average - z$components[[h]]$mean)^2 +
           z$components[[h]]$variance))
    }))
  })) / sum(vapply(groups, function(z) z$n, numeric(1)))
  list(means = means, between = between, within = within,
       weights = membership / length(groups))
}

set.seed(20260930)
parameters <- list(means = rbind(c(-1, 0.5), c(1.4, -0.8)),
                   between = rbind(c(0.3, 0.7), c(0.8, 0.2)),
                   within = c(0.6, 1.1), weights = c(0.4, 0.6))
sizes <- c(1L, 2L, 3L, 5L, 7L, 4L, 6L, 2L, 8L, 3L, 5L, 4L)
truth <- sample.int(2L, length(sizes), replace = TRUE, prob = parameters$weights)
observations <- lapply(seq_along(sizes), function(j) {
  intercept <- rnorm(2L, parameters$means[truth[j], ],
                     sqrt(parameters$between[truth[j], ]))
  sweep(sweep(matrix(rnorm(sizes[j] * 2L), sizes[j], 2L), 2L,
              sqrt(parameters$within), "*"), 2L, intercept, "+")
})
str(observations[1:2])
print(head(observations[[5L]]))
print(list(groups = length(observations), rows = sum(sizes),
           parameters = names(parameters), class = class(observations)))

evaluated <- additive_design_evaluate(observations, parameters)
compact <- unlist(lapply(evaluated$groups, function(z) {
  vapply(z$components, function(component) component$log_density, numeric(1))
}))
dense <- unlist(lapply(observations, function(x) {
  vapply(seq_len(2L), function(h) additive_design_dense(
    x, parameters$means[h, ], parameters$between[h, ], parameters$within), numeric(1))
}))
stopifnot(max(abs(compact - dense)) < 1e-10,
          all(vapply(evaluated$groups, function(z) {
            abs(sum(z$posterior) - 1) < 1e-12
          }, logical(1))))

# Independent one-dimensional integration over the group's random intercept.
x <- observations[[3L]][, 1L]
integral <- integrate(function(b) {
  vapply(b, function(value) prod(dnorm(x, value, sqrt(parameters$within[1L]))) *
    dnorm(value, parameters$means[1L, 1L], sqrt(parameters$between[1L, 1L])),
    numeric(1))
}, -Inf, Inf, rel.tol = 1e-11, abs.tol = 1e-25)
integration_error <- abs(log(integral$value) - additive_design_dense(
  matrix(x, ncol = 1L), parameters$means[1L, 1L],
  parameters$between[1L, 1L], parameters$within[1L]))
stopifnot(integration_error < 1e-9)

# Compare conditional intercept moments with general Gaussian conditioning.
moment_error <- max(unlist(lapply(seq_along(observations), function(j) {
  x <- observations[[j]]
  unlist(lapply(seq_len(2L), function(h) {
    vapply(seq_len(2L), function(r) {
      t <- parameters$between[h, r]
      covariance <- diag(parameters$within[r], nrow(x)) + matrix(t, nrow(x), nrow(x))
      gain <- rep(t, nrow(x)) %*% solve(covariance)
      mean <- parameters$means[h, r] + drop(gain %*% (x[, r] - parameters$means[h, r]))
      variance <- t - sum(gain) * t
      component <- evaluated$groups[[j]]$components[[h]]
      max(abs(c(mean - component$mean[r], variance - component$variance[r])))
    }, numeric(1))
  }))
})))
stopifnot(moment_error < 1e-10)

updated <- additive_design_step(evaluated)
increment <- additive_design_evaluate(observations, updated)$log_likelihood -
  evaluated$log_likelihood
stopifnot(increment >= -1e-10, all(updated$within > 0), all(updated$between > 0))

permuted <- parameters
permuted$means <- parameters$means[2:1, , drop = FALSE]
permuted$between <- parameters$between[2:1, , drop = FALSE]
permuted$weights <- parameters$weights[2:1]
label_error <- abs(additive_design_evaluate(observations, permuted)$log_likelihood -
                     evaluated$log_likelihood)
row_error <- abs(additive_design_evaluate(lapply(observations, function(x) {
  x[rev(seq_len(nrow(x))), , drop = FALSE]
}), parameters)$log_likelihood - evaluated$log_likelihood)
stopifnot(label_error < 1e-10, row_error < 1e-10)

zero_between <- parameters
zero_between$between[,] <- 0
zero_eval <- additive_design_evaluate(observations, zero_between)
zero_dense <- unlist(lapply(observations, function(x) {
  vapply(seq_len(2L), function(h) sum(vapply(seq_len(2L), function(r) {
    sum(dnorm(x[, r], parameters$means[h, r], sqrt(parameters$within[r]), log = TRUE))
  }, numeric(1))), numeric(1))
}))
zero_error <- max(abs(unlist(lapply(zero_eval$groups, function(z) {
  vapply(z$components, function(component) component$log_density, numeric(1))
})) - zero_dense))
stopifnot(zero_error < 1e-10)

# Singleton groups cannot separate within from between variance.
singletons <- lapply(observations, function(x) x[1L, , drop = FALSE])
shifted <- parameters
shifted$within <- shifted$within + 0.05
shifted$between <- shifted$between - 0.05
singleton_error <- abs(additive_design_evaluate(singletons, parameters)$log_likelihood -
                         additive_design_evaluate(singletons, shifted)$log_likelihood)
stopifnot(singleton_error < 1e-10,
          abs(additive_design_evaluate(observations, shifted)$log_likelihood -
                evaluated$log_likelihood) > 1e-4)

# Check balanced one-class interior and boundary MLEs against a separate
# optimizer using the dense covariance likelihood, not the EM equations.
balanced_errors <- vapply(list(seq(-2, 3, length.out = 6L),
                              seq(-0.05, 0.05, length.out = 6L)), function(shifts) {
  balanced <- lapply(shifts, function(shift) {
    matrix(shift + c(-1, -0.3, 0.3, 1), ncol = 1L)
  })
  averages <- vapply(balanced, mean, numeric(1))
  mu <- mean(averages)
  scatter <- sum(vapply(balanced, function(z) sum((z - mean(z))^2), numeric(1)))
  spread <- sum((averages - mu)^2)
  within <- scatter / (length(balanced) * 3L)
  between <- spread / length(balanced) - within / 4L
  if (between < 0) {
    between <- 0
    within <- (scatter + 4L * spread) / (length(balanced) * 4L)
  }
  target <- c(mu, log(within), between)
  objective <- function(theta) {
    -sum(vapply(balanced, function(z) additive_design_dense(
      z, theta[1L], theta[3L], exp(theta[2L])), numeric(1)))
  }
  optimized <- optim(target + c(0.1, 0.1, 0.1), objective,
                     method = "L-BFGS-B", lower = c(-Inf, -Inf, 0),
                     control = list(factr = 1e4, pgtol = 1e-9, ndeps = rep(1e-5, 3L)))
  stopifnot(optimized$convergence == 0L,
            max(abs(optimized$par - target)) < 1e-4,
            abs(optimized$value - objective(target)) < 1e-8)
  abs(optimized$value - objective(target))
}, numeric(1))

print(data.frame(check = c("compact versus dense", "numerical integration",
                          "conditional moments", "class permutation",
                          "row permutation", "zero between variance", "singleton ridge",
                          "balanced interior MLE", "balanced boundary MLE"),
                 absolute_error = c(max(abs(compact - dense)), integration_error,
                                    moment_error, label_error, row_error, zero_error,
                                    singleton_error, balanced_errors)))
cat(sprintf("Interior EM log-likelihood increment: %.10f\n", increment))
cat("All additive design checks passed; recovery and coverage remain untested.\n")
cat(sprintf("Runtime: %s; synthetic seed: 20260930\n", R.version.string))
