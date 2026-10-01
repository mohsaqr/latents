# SQUAREM acceleration of the multilpa EM map.
#
# Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods
# for accelerating the convergence of any EM algorithm. Scandinavian Journal
# of Statistics, 35, 335-353.
#
# One cycle takes two EM steps from theta0 (theta1, theta2), forms the
# first and second differences r = theta1 - theta0 and v = theta2 - theta1 - r,
# extrapolates to theta' = theta0 - 2 alpha r + alpha^2 v with the steplength
# alpha = -|r| / |v| (scheme S3), and finishes with one EM step from theta'.
# That last step is what keeps a constrained covariance structure on its
# constraint set: the extrapolated point may leave it, the M-step output never
# does. alpha = -1 gives theta' = theta2, plain EM, which is the fallback.

#' The parameter blocks EM updates and SQUAREM may extrapolate
#' @noRd
.multilpa_accelerated_blocks <- function() {
  c("means", "variances", "covariances", "profile_probabilities",
    "group_probabilities", "response_probabilities", "ordinal_intercepts",
    "ordinal_locations", "count_means", "count_dispersion")
}

#' Flatten the updatable blocks of a parameter list
#' @param parameters A parameter list.
#' @return A numeric vector, blocks in `.multilpa_accelerated_blocks()` order.
#' @noRd
.multilpa_flatten_parameters <- function(parameters) {
  blocks <- parameters[intersect(.multilpa_accelerated_blocks(),
                                 names(parameters))]
  unlist(blocks, use.names = FALSE)
}

#' Refill a parameter list's updatable blocks from a flat vector
#'
#' The inverse of `.multilpa_flatten_parameters()`: every array, matrix and
#' nested list keeps its shape and names; only the numbers change. Blocks
#' outside the updatable set (a noise density) are carried over unchanged.
#'
#' @param template A parameter list with the target shapes.
#' @param values A numeric vector of the flattened length.
#' @return A parameter list.
#' @noRd
.multilpa_restore_parameters <- function(template, values) {
  position <- 0L
  fill <- function(block) {
    if (is.list(block)) return(lapply(block, fill))
    size <- length(block)
    block[] <- values[position + seq_len(size)]
    position <<- position + size
    block
  }
  names_present <- intersect(.multilpa_accelerated_blocks(), names(template))
  template[names_present] <- lapply(template[names_present], fill)
  stopifnot("the flat vector must match the template's size" =
              position == length(values))
  template
}

#' Can the E-step be evaluated at this extrapolated point?
#'
#' Probabilities must be positive and variances positive, and every covariance
#' matrix positive definite. The extrapolation preserves the probabilities'
#' row sums, because the differences it combines sum to zero.
#'
#' @param parameters A parameter list.
#' @return A single logical.
#' @noRd
.multilpa_valid_parameters <- function(parameters) {
  probabilities <- unlist(parameters[c("profile_probabilities",
                                       "group_probabilities",
                                       "response_probabilities")],
                          use.names = FALSE)
  if (any(!is.finite(probabilities)) || any(probabilities <= 0) ||
      any(probabilities > 1 + 1e-12)) {
    return(FALSE)
  }
  if (length(c(parameters$count_means, parameters$count_dispersion)) > 0L &&
      (any(!is.finite(c(parameters$count_means, parameters$count_dispersion))) ||
       any(c(parameters$count_means, parameters$count_dispersion) <= 0))) {
    return(FALSE)
  }
  if (length(parameters$variances) > 0L &&
      (any(!is.finite(parameters$variances)) || any(parameters$variances <= 0))) {
    return(FALSE)
  }
  covariances <- parameters$covariances
  if (is.null(covariances)) return(TRUE)
  all(vapply(seq_len(dim(covariances)[3L]), function(profile) {
    block <- matrix(covariances[, , profile], dim(covariances)[1L])
    all(is.finite(block)) &&
      !inherits(tryCatch(chol(block), error = function(error) error), "error")
  }, logical(1)))
}

#' One SQUAREM cycle
#'
#' @param step A function mapping a parameter list and its E-step to the next
#'   parameter list and E-step: `list(parameters, expectation)`.
#' @param parameters Current parameters.
#' @param expectation Their E-step.
#' @param evaluate A function giving the E-step of a parameter list.
#' @return A list with `parameters`, `expectation`, `history` (the three
#'   likelihoods along the accepted path) and `accelerated` (whether the
#'   extrapolated point was kept).
#' @noRd
.multilpa_squarem_cycle <- function(step, parameters, expectation, evaluate) {
  first <- step(parameters, expectation)
  second <- step(first$parameters, first$expectation)
  theta0 <- .multilpa_flatten_parameters(parameters)
  r <- .multilpa_flatten_parameters(first$parameters) - theta0
  v <- .multilpa_flatten_parameters(second$parameters) - theta0 - 2 * r
  plain <- list(parameters = second$parameters,
                expectation = second$expectation,
                history = c(first$expectation$log_likelihood,
                            second$expectation$log_likelihood,
                            second$expectation$log_likelihood),
                accelerated = FALSE)
  norm_v <- sqrt(sum(v^2))
  if (!is.finite(norm_v) || norm_v <= 0) return(plain)
  alpha <- min(-1, -sqrt(sum(r^2)) / norm_v)
  # Backtrack towards alpha = -1 (plain EM, always valid) until the E-step is
  # defined at the extrapolated point. Halving the distance to -1 reaches a
  # valid point in a bounded number of tries.
  candidate <- NULL
  tries <- 0L
  while (is.null(candidate) && alpha < -1 && tries < 30L) {
    theta <- theta0 - 2 * alpha * r + alpha^2 * v
    proposal <- .multilpa_restore_parameters(parameters, theta)
    if (.multilpa_valid_parameters(proposal)) candidate <- proposal else
      alpha <- (alpha - 1) / 2
    tries <- tries + 1L
  }
  if (is.null(candidate)) return(plain)
  # The E-step at the extrapolated point can still fail numerically (a
  # density underflow on every component); that point is simply not taken.
  landed <- tryCatch(step(candidate, evaluate(candidate)),
                     error = function(error) NULL)
  if (is.null(landed) ||
      !is.finite(landed$expectation$log_likelihood) ||
      landed$expectation$log_likelihood < second$expectation$log_likelihood) {
    return(plain)
  }
  list(parameters = landed$parameters, expectation = landed$expectation,
       history = c(first$expectation$log_likelihood,
                   second$expectation$log_likelihood,
                   landed$expectation$log_likelihood),
       accelerated = TRUE)
}
