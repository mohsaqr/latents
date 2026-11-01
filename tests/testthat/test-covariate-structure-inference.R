# Wald inference for membership-covariate models in the ten covariance
# structures that tie the profiles together. The package estimates them in
# its chart (log volumes, log shapes, Cayley rotations; structure-inference.R).
# The check here is a likelihood written independently for two indicators:
# Sigma_k = lambda_k R(angle_k) diag(a_k, 1 / a_k) R(angle_k)', with each of
# the volume, shape and angle shared or per profile as the structure says,
# and a logit on z for membership. Means and membership coefficients do not
# depend on how the spread is parameterized, so their standard errors must
# equal the inverse numDeriv Hessian of this likelihood; the variances and
# covariances follow by the delta method through its own Jacobian.

skip_on_cran()
skip_if_not_installed("numDeriv")

structure_case <- function() {
  set.seed(23)
  n <- 500L
  z <- stats::rnorm(n)
  state <- 1L + stats::rbinom(n, 1L, stats::plogis(-0.2 + 1.2 * z))
  first <- matrix(c(1, 0.35, 0.35, 0.8), 2L)
  second <- matrix(c(1.7, -0.4, -0.4, 1.1), 2L)
  noise <- matrix(stats::rnorm(2L * n), n, 2L)
  rows <- t(vapply(seq_len(n), function(i) {
    c(0, 0, 2.2, 1.6)[c(1L, 2L) + 2L * (state[i] - 1L)] +
      drop(noise[i, ] %*% chol(if (state[i] == 1L) first else second))
  }, numeric(2)))
  data.frame(y1 = rows[, 1L], y2 = rows[, 2L], z = z)
}

# Which of volume, shape and orientation each structure frees per profile;
# NA for a shape or orientation the structure fixes (spherical, axis-aligned).
structure_layout <- function(structure) {
  letters <- strsplit(structure, "")[[1L]]
  free <- function(letter) if (letter == "I") NA else letter == "V"
  list(volume = letters[[1L]] == "V", shape = free(letters[[2L]]),
       angle = free(letters[[3L]]))
}

# The independent parameter vector: means (profile by profile), log volume(s),
# log shape(s), angle(s), membership intercept and slope.
independent_model <- function(structure, data) {
  layout <- structure_layout(structure)
  widths <- c(means = 4L, volume = if (layout$volume) 2L else 1L,
              shape = if (is.na(layout$shape)) 0L else if (layout$shape) 2L else 1L,
              angle = if (is.na(layout$angle)) 0L else if (layout$angle) 2L else 1L,
              membership = 2L)
  ends <- cumsum(widths)
  block <- function(name, theta) theta[seq_len(widths[[name]]) + ends[[name]] - widths[[name]]]
  per_profile <- function(values, k) if (length(values) == 0L) 0 else values[min(k, length(values))]
  covariances <- function(theta) lapply(1:2, function(k) {
    angle <- per_profile(block("angle", theta), k)
    rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L)
    shape <- exp(per_profile(block("shape", theta), k))
    exp(per_profile(block("volume", theta), k)) *
      rotation %*% diag(c(shape, 1 / shape)) %*% t(rotation)
  })
  y <- as.matrix(data[c("y1", "y2")])
  log_likelihood <- function(theta) {
    means <- matrix(block("means", theta), 2L, byrow = TRUE)
    membership <- block("membership", theta)
    first_share <- stats::plogis(membership[[1L]] + membership[[2L]] * data$z)
    sigma <- covariances(theta)
    density <- vapply(1:2, function(k) {
      centred <- sweep(y, 2L, means[k, ])
      precision <- solve(sigma[[k]])
      exp(-0.5 * rowSums((centred %*% precision) * centred)) /
        (2 * pi * sqrt(det(sigma[[k]])))
    }, numeric(nrow(y)))
    sum(log(first_share * density[, 1L] + (1 - first_share) * density[, 2L]))
  }
  list(layout = layout, widths = widths, block = block, covariances = covariances,
       log_likelihood = log_likelihood)
}

# The independent coordinates at a fit's estimate: volumes and shapes from
# the eigenvalues, angles from the eigenvectors (a shared orientation from
# the first profile's), then polished to the independent likelihood's own
# maximum.
independent_estimate <- function(model, fit) {
  sigma <- lapply(1:2, function(k) fit$covariances[, , k] %||% diag(fit$variances[k, ]))
  axes <- lapply(sigma, function(block) eigen(block, symmetric = TRUE)$vectors)
  reference <- if (isFALSE(model$layout$angle)) axes[1L] else axes
  along <- lapply(1:2, function(k) {
    frame <- reference[[min(k, length(reference))]]
    diag(crossprod(frame, sigma[[k]] %*% frame))
  })
  volumes <- vapply(along, function(values) sqrt(prod(values)), numeric(1))
  shapes <- vapply(1:2, function(k) log(along[[k]][[1L]] / volumes[[k]]), numeric(1))
  angles <- vapply(reference, function(frame) atan2(frame[2L, 1L], frame[1L, 1L]), numeric(1))
  take <- function(values, free) if (is.na(free)) numeric(0) else if (free) values else values[1L]
  start <- c(as.vector(t(fit$means)), take(log(volumes), model$layout$volume),
             take(shapes, model$layout$shape), take(angles, model$layout$angle),
             as.vector(fit$profile_coefficients))
  stats::optim(start, function(theta) -model$log_likelihood(theta), method = "BFGS",
               control = list(reltol = 1e-15, maxit = 2000L))
}

constrained <- c("EII", "VII", "VEI", "EVI", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV")

test_that("covariate fits in the constrained structures have Wald errors of their likelihood", {
  data <- structure_case()
  invisible(lapply(constrained, function(structure) {
    fit <- quietly(lpa(data, c("y1", "y2"), n_profiles = 2, model = structure,
                       profile_covariates = "z", n_starts = 10, seed = 1, tol = 1e-12,
                       max_iter = 20000))
    model <- independent_model(structure, data)
    optimum <- independent_estimate(model, fit)
    # The two parameterizations reach the same maximum.
    expect_equal(-optimum$value, fit$log_likelihood, tolerance = 1e-8,
                 label = paste(structure, "log likelihood"))
    expect_identical(length(optimum$par), as.integer(fit$n_parameters),
                     label = paste(structure, "parameters"))
    covariance <- solve(-numDeriv::hessian(model$log_likelihood, optimum$par))
    table <- parameter_inference(fit)
    means <- subset(table, parameter == "mean")
    expect_equal(means$standard_error, sqrt(diag(covariance))[seq_len(4L)], tolerance = 1e-4,
                 label = paste(structure, "mean standard errors"))
    coefficients <- subset(table, parameter == "coefficient")
    membership <- length(optimum$par) - 1:0
    expect_equal(coefficients$standard_error, sqrt(diag(covariance))[membership],
                 tolerance = 1e-4, label = paste(structure, "coefficient standard errors"))
    # The spread entries the table reports, by the delta method through the
    # independent map, in the table's order.
    spread <- subset(table, parameter %in% c("variance", "covariance"))
    natural <- function(theta) {
      blocks <- model$covariances(theta)
      shared <- all(spread$outcome == "shared")
      entries <- if (all(spread$parameter == "covariance")) {
        function(block) block[lower.tri(block, diag = TRUE)]
      } else diag
      unlist(lapply(if (shared) blocks[1L] else blocks, entries))
    }
    expect_equal(spread$estimate, natural(optimum$par), tolerance = 1e-5,
                 label = paste(structure, "spread estimates"))
    jacobian <- numDeriv::jacobian(natural, optimum$par)
    expect_equal(spread$standard_error,
                 sqrt(diag(jacobian %*% covariance %*% t(jacobian))), tolerance = 1e-4,
                 label = paste(structure, "spread standard errors"))
  }))
})

test_that("the chart has one coordinate per parameter and maps back to the fit", {
  data <- structure_case()
  fit <- quietly(lpa(data, c("y1", "y2"), n_profiles = 2, model = "VEV",
                     profile_covariates = "z", n_starts = 5, seed = 1))
  unconstrained <- coef(fit, scale = "unconstrained")
  expect_identical(length(unconstrained), as.integer(fit$n_parameters))
  expect_true(any(startsWith(names(unconstrained), "measurement.rotation")))
  expect_identical(dim(vcov(fit, scale = "unconstrained")), rep(length(unconstrained), 2L))
  natural <- coef(fit)
  covariance_rows <- grepl("covariance", names(natural), fixed = TRUE)
  expect_equal(unname(natural[covariance_rows]),
               c(fit$covariances[, , 1L][lower.tri(diag(2), diag = TRUE)],
                 fit$covariances[, , 2L][lower.tri(diag(2), diag = TRUE)]), tolerance = 1e-10)
  # Robust errors are available as for the separable structures.
  robust <- parameter_inference(fit, vcov_type = "robust")
  expect_true(all(is.finite(robust$standard_error)))
})
