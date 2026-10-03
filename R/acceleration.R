# What SQUAREM needs from the profile models' parameter lists: which blocks
# it may extrapolate, how to flatten and refill them, and whether an
# extrapolated point is valid. The cycle itself is the kernel's
# (`.latents_squarem_cycle()` in kernel-em.R), shared by every engine that
# accelerates.
#
# Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods
# for accelerating the convergence of any EM algorithm. Scandinavian Journal
# of Statistics, 35, 335-353.

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
