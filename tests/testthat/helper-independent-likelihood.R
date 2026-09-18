# Independent numerical references for tiny synthetic examples. These deliberately
# use direct probability arithmetic and exhaustive assignments, rather than the
# fitted model's log-domain likelihood or E-step. Do not use on large datasets.

#' Enumerate every latent assignment in each small group
#' @param x Numeric matrix of indicators.
#' @param group Cluster identifiers in row order.
#' @param parameters List of means, variances, profile and group probabilities.
#' @return Independently calculated likelihood, posteriors and joint posteriors.
enumerate_latent_assignments <- function(x, group, parameters) {
  stopifnot(is.matrix(x), is.numeric(x), !anyNA(x),
            length(group) == nrow(x), !anyNA(group), is.list(parameters),
            all(c("means", "variances", "profile_probabilities",
                  "group_probabilities") %in% names(parameters)))
  means <- parameters$means
  variances <- parameters$variances
  profile_probabilities <- parameters$profile_probabilities
  group_probabilities <- parameters$group_probabilities
  n_profiles <- nrow(means)
  n_group_classes <- length(group_probabilities)
  stopifnot(identical(dim(means), dim(variances)), ncol(means) == ncol(x),
            all(variances > 0), all(is.finite(means)),
            all(is.finite(variances)),
            identical(dim(profile_probabilities), c(n_group_classes, n_profiles)),
            all(profile_probabilities >= 0), all(group_probabilities >= 0),
            max(abs(rowSums(profile_probabilities) - 1)) < 1e-8,
            abs(sum(group_probabilities) - 1) < 1e-8)
  group_indices <- lapply(unique(group), function(value) {
    stopifnot(length(value) == 1L)
    which(group == value)
  })
  stopifnot(max(lengths(group_indices)) <= 12L)
  group_results <- lapply(group_indices, function(indices) {
    stopifnot(is.integer(indices), length(indices) > 0L)
    assignments <- as.matrix(expand.grid(rep(list(seq_len(n_profiles)),
                                             length(indices))))
    densities <- vapply(seq_len(n_profiles), function(profile) {
      stopifnot(length(profile) == 1L)
      apply(x[indices, , drop = FALSE], 1L, function(observation) {
        stopifnot(is.numeric(observation), length(observation) == ncol(x))
        prod(stats::dnorm(observation, means[profile, ],
                          sqrt(variances[profile, ])))
      })
    }, numeric(length(indices)))
    dim(densities) <- c(length(indices), n_profiles)
    assignment_weights <- vapply(seq_len(n_group_classes), function(group_class) {
      stopifnot(length(group_class) == 1L)
      apply(assignments, 1L, function(assignment) {
        stopifnot(length(assignment) == length(indices))
        group_probabilities[group_class] *
          prod(profile_probabilities[group_class, assignment]) *
          prod(densities[cbind(seq_along(indices), assignment)])
      })
    }, numeric(nrow(assignments)))
    dim(assignment_weights) <- c(nrow(assignments), n_group_classes)
    likelihood <- sum(assignment_weights)
    stopifnot(is.finite(likelihood), likelihood > 0)
    normalized_weights <- assignment_weights / likelihood
    joint <- array(0, c(length(indices), n_group_classes, n_profiles))
    invisible(lapply(seq_len(n_profiles), function(profile) {
      stopifnot(length(profile) == 1L)
      joint[, , profile] <<- t(assignments == profile) %*% normalized_weights
      NULL
    }))
    list(log_likelihood = log(likelihood), group = colSums(normalized_weights),
         subject = apply(joint, c(1L, 3L), sum), joint = joint)
  })
  subject_posteriors <- matrix(0, nrow(x), n_profiles)
  joint_posteriors <- array(0, c(nrow(x), n_group_classes, n_profiles))
  invisible(lapply(seq_along(group_indices), function(group) {
    stopifnot(length(group) == 1L)
    subject_posteriors[group_indices[[group]], ] <<- group_results[[group]]$subject
    joint_posteriors[group_indices[[group]], , ] <<- group_results[[group]]$joint
    NULL
  }))
  list(log_likelihood = sum(vapply(group_results, `[[`, numeric(1L),
                                   "log_likelihood")),
       subject_posteriors = subject_posteriors,
       group_posteriors = do.call(rbind, lapply(group_results, `[[`, "group")),
       joint_posteriors = joint_posteriors)
}

#' Calculate one EM update using exhaustive-assignment posteriors
#' @param x Numeric indicator matrix.
#' @param group Cluster identifiers.
#' @param parameters Starting parameters.
#' @param variance_model Whether variances vary by profile or are equal.
#' @param min_variance Positive variance floor.
#' @return Updated parameter list.
enumerated_em_update <- function(x, group, parameters,
                                 variance_model = "varying", min_variance = 1e-6) {
  stopifnot(is.matrix(x), is.numeric(x), length(group) == nrow(x),
            is.list(parameters), variance_model %in% c("varying", "equal"),
            length(min_variance) == 1L, min_variance > 0)
  posterior <- enumerate_latent_assignments(x, group, parameters)
  weights <- posterior$subject_posteriors
  means <- crossprod(weights, x) / colSums(weights)
  squared_residual_sums <- t(matrix(vapply(seq_len(ncol(weights)), function(profile) {
    stopifnot(length(profile) == 1L)
    colSums(weights[, profile] * sweep(x, 2L, means[profile, ])^2)
  }, numeric(ncol(x))), nrow = ncol(x), ncol = ncol(weights)))
  variances <- squared_residual_sums / colSums(weights)
  if (variance_model == "equal") {
    variances <- matrix(rep(colSums(squared_residual_sums) / nrow(x),
                            each = nrow(means)), nrow = nrow(means))
  }
  variances <- pmax(variances, min_variance)
  profile_counts <- apply(posterior$joint_posteriors, c(2L, 3L), sum)
  list(means = means, variances = variances,
       profile_probabilities = profile_counts / rowSums(profile_counts),
       group_probabilities = colMeans(posterior$group_posteriors))
}

#' Independently optimize a two-profile, two-group, one-indicator likelihood
#' @param x Numeric single-column indicator matrix.
#' @param group Cluster identifiers.
#' @param start Starting model parameters.
#' @return An optim result and unpacked parameter estimates.
optim_multilpa_reference <- function(x, group, start) {
  stopifnot(is.matrix(x), is.numeric(x), ncol(x) == 1L, !anyNA(x),
            length(group) == nrow(x), is.list(start),
            identical(dim(start$means), c(2L, 1L)),
            identical(dim(start$variances), c(2L, 1L)),
            identical(dim(start$profile_probabilities), c(2L, 2L)),
            length(start$group_probabilities) == 2L)
  initial <- c(start$means, log(start$variances),
               stats::qlogis(start$profile_probabilities[, 1L]),
               stats::qlogis(start$group_probabilities[1L]))
  unpack <- function(theta) {
    stopifnot(is.numeric(theta), length(theta) == 7L)
    proportions <- stats::plogis(theta[5:6])
    group_probability <- stats::plogis(theta[7L])
    list(means = matrix(theta[1:2], 2L, 1L),
         variances = matrix(exp(theta[3:4]), 2L, 1L),
         profile_probabilities = cbind(proportions, 1 - proportions),
         group_probabilities = c(group_probability, 1 - group_probability))
  }
  indices <- lapply(unique(group), function(value) {
    stopifnot(length(value) == 1L)
    which(group == value)
  })
  objective <- function(theta) {
    stopifnot(is.numeric(theta), length(theta) == 7L)
    parameters <- unpack(theta)
    densities <- vapply(seq_len(2L), function(profile) {
      stopifnot(length(profile) == 1L)
      stats::dnorm(x[, 1L], parameters$means[profile, 1L],
                   sqrt(parameters$variances[profile, 1L]))
    }, numeric(nrow(x)))
    within_group_type <- densities %*% t(parameters$profile_probabilities)
    likelihoods <- vapply(indices, function(rows) {
      stopifnot(is.integer(rows), length(rows) > 0L)
      sum(parameters$group_probabilities *
            apply(within_group_type[rows, , drop = FALSE], 2L, prod))
    }, numeric(1L))
    if (any(!is.finite(likelihoods)) || any(likelihoods <= 0)) return(1e100)
    -sum(log(likelihoods))
  }
  result <- stats::optim(initial, objective, method = "BFGS",
                          control = list(maxit = 2000L, reltol = 1e-12))
  list(optim = result, parameters = unpack(result$par))
}
