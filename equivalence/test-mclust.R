# Covariance-structure M-steps, parameter counts and likelihoods against mclust.
activity <- c("browse", "lectures", "forum_read")

test_that("every M-step reproduces mclust's own, or beats it", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  set.seed(3)
  n <- 400L; d <- 3L; k <- 3L
  X <- matrix(stats::rnorm(n * d), n, d)
  X[1:150, ] <- X[1:150, ] %*% diag(c(3, 1, 0.5)) + 2
  X[151:300, ] <- X[151:300, ] %*% diag(c(0.5, 2, 1)) - 1
  z <- matrix(stats::runif(n * k), n, k)
  z <- z / rowSums(z)
  weights <- colSums(z)
  means <- t(crossprod(z, X) / weights)
  scatter <- lapply(seq_len(k), function(profile) {
    residuals <- sweep(X, 2L, means[, profile], "-")
    crossprod(residuals, residuals * z[, profile])
  })
  # The M-step minimises this; whichever solution attains less is the better
  # maximiser, whatever the two disagree about elementwise.
  objective <- function(sigma) {
    sum(vapply(seq_len(k), function(profile) {
      block <- matrix(sigma[, , profile], d)
      weights[profile] *
        as.numeric(determinant(block, logarithm = TRUE)$modulus) +
        sum(chol2inv(chol(block)) * scatter[[profile]])
    }, numeric(1)))
  }
  invisible(lapply(multilpa:::.multilpa_structures(), function(code) {
    # `mstep()` resolves its per-model worker from the calling frame, so the
    # call is evaluated in mclust's own namespace.
    reference <- do.call(mclust::mstep, list(modelName = code, data = X, z = z),
                         envir = asNamespace("mclust"))
    theirs <- reference$parameters$variance$sigma
    ours <- multilpa:::.multilpa_structure_covariances(scatter, weights, code,
                                                       1e-12)
    # The twelve with a closed form or a scalar fixed point agree exactly. The
    # two whose orientation has no closed form need only be no worse.
    if (code %in% c("EVE", "VVE")) {
      expect_lte(objective(ours), objective(theirs) + 1e-8, label = code)
      expect_equal(ours, theirs, tolerance = 1e-2, ignore_attr = TRUE,
                   info = code)
    } else {
      expect_equal(ours, theirs, tolerance = 1e-8, ignore_attr = TRUE,
                   info = code)
    }
  }))
})

test_that("the parameter counts are the ones mclust reports", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  # `Mclust()` resolves `mclustBIC` from the calling frame, so a namespaced
  # call needs it bound here.
  mclustBIC <- mclust::mclustBIC
  count <- multilpa:::.multilpa_structure_parameters
  X <- as.matrix(course_engagement[activity])
  codes <- multilpa:::.multilpa_structures()
  invisible(lapply(codes, function(code) {
    reference <- mclust::Mclust(X, G = 3L, modelNames = code, verbose = FALSE)
    skip_if(is.null(reference))
    # mclust's `df` counts the means, the mixing proportions and the spread;
    # this package's counter is the spread alone.
    spread <- reference$df - 3L * length(activity) - 2L
    expect_identical(count(code, 3L, length(activity)), as.integer(spread),
                     info = code)
  }))
})

test_that("the likelihood is the mixture likelihood, at mclust's own estimate", {
  skip_if_not_installed("mclust")
  skip_on_cran()
  mclustBIC <- mclust::mclustBIC
  frame <- course_engagement[c("student", activity)]
  X <- as.matrix(frame[activity])
  reference <- mclust::Mclust(X, G = 3L, modelNames = "EVI", verbose = FALSE)
  skip_if(is.null(reference))
  means <- t(reference$parameters$mean)
  variances <- t(vapply(seq_len(3L), function(k) {
    diag(matrix(reference$parameters$variance$sigma[, , k], length(activity)))
  }, numeric(length(activity))))
  proportions <- reference$parameters$pro

  # The mixture log likelihood written from its definition.
  densities <- vapply(seq_len(3L), function(k) {
    proportions[k] * exp(rowSums(stats::dnorm(
      X, rep(means[k, ], each = nrow(X)),
      rep(sqrt(variances[k, ]), each = nrow(X)), log = TRUE)))
  }, numeric(nrow(X)))
  by_hand <- sum(log(rowSums(densities)))

  # `max_iter = 0` evaluates at the supplied values without moving them, so
  # this compares two likelihood functions rather than two optimisers.
  held <- multilpa(frame, activity, "student", n_profiles = 3L,
                   n_group_classes = 1, n_starts = 1, max_iter = 0, seed = 1,
                   volume = "equal", shape = "varying",
                   start = list(means = means, variances = variances,
                                profile_probabilities = matrix(proportions, 1L, 3L),
                                group_probabilities = 1))
  expect_equal(held$log_likelihood, by_hand, tolerance = 1e-10)
})

test_that("full covariance mixture limits reproduce mclust VVV and EEE", {
  skip_if_not_installed("mclust")
  mclustBIC <- mclust::mclustBIC
  set.seed(65)
  data <- data.frame(group = rep(seq_len(40), each = 4),
                     x = c(rnorm(80, -3), rnorm(80, 3)), y = rnorm(160))
  data$y <- data$y + .6 * data$x
  x <- as.matrix(data[c("x", "y")])
  invisible(lapply(c(varying = "VVV", equal = "EEE"), function(model) {
    reference <- mclust::Mclust(x, G = 2, modelNames = model, verbose = FALSE,
                                 control = mclust::emControl(tol = c(1e-12, 1e-12)))
    start <- list(means = t(reference$parameters$mean),
                   covariances = reference$parameters$variance$sigma,
                   profile_probabilities = matrix(reference$parameters$pro, 1L),
                   group_probabilities = 1)
    fit <- multilpa(data, c("x", "y"), "group", 2, 1, covariance_model = "full",
                      variance_model = if (model == "VVV") "varying" else "equal",
                      n_starts = 1, start = start, tol = 1e-13)
    expect_equal(fit$log_likelihood, as.numeric(reference$loglik), tolerance = 1e-7)
    expect_equal(unname(fit$means), unname(start$means), tolerance = 2e-5)
    expect_equal(unname(fit$covariances), unname(start$covariances), tolerance = 3e-5)
    expect_equal(unname(fit$subject_posteriors), unname(reference$z), tolerance = 2e-5)
    expect_equal(fit$n_parameters, if (model == "VVV") 11 else 8)
  }))
})

