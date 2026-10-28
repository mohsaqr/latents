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
    "group_probabilities", "initial_probabilities", "transition_probabilities",
    "response_probabilities", "ordinal_intercepts",
    "ordinal_locations", "count_means", "count_dispersion", "direct_effects")
}

#' Flatten the updatable blocks of a parameter list
#' @param parameters A parameter list.
#' @param blocks Names of the blocks to flatten.
#' @return A numeric vector, blocks in `.multilpa_accelerated_blocks()` order.
#' @noRd
.multilpa_flatten_parameters <- function(parameters,
                                         blocks = .multilpa_accelerated_blocks()) {
  blocks <- parameters[intersect(blocks, names(parameters))]
  unlist(blocks, use.names = FALSE)
}

#' SQUAREM coordinates in reference standard-deviation units
#'
#' Only the acceleration chart changes: the E-step, M-step and covariance
#' constraints still operate in the supplied units. Non-Gaussian blocks keep
#' their existing coordinates. The reference spread is fixed across starts
#' and resumed cycles, using observed rows of positive sampling weight.
#' @param x Centred continuous indicators.
#' @param parameters A parameter list defining the block shapes.
#' @param row_weights Optional sampling weights in row order.
#' @param extra Optional extra-indicator and direct-effect carrier.
#' @return The flatten, restore and valid functions for the EM kernel.
#' @noRd
.multilpa_accelerator <- function(x, parameters, row_weights = NULL, extra = NULL) {
  rows <- if (is.null(row_weights)) rep(TRUE, nrow(x)) else row_weights > 0
  scales <- vapply(seq_len(ncol(x)), function(column) {
    stats::sd(x[rows, column], na.rm = TRUE)
  }, numeric(1))
  # A constant or singly observed column has no empirical unit. Its bound
  # supplies a positive reference in variance units instead.
  unusable <- !is.finite(scales) | scales <= 0
  if (any(unusable)) {
    scales[unusable] <- sqrt(colMeans(parameters$variances)[unusable])
  }
  chart <- .multilpa_restore_parameters(
    parameters, rep(1, length(.multilpa_flatten_parameters(parameters))))
  chart$means <- sweep(chart$means, 2L, scales, "*")
  chart$variances <- sweep(chart$variances, 2L, scales^2, "*")
  if (!is.null(chart$covariances)) {
    chart$covariances <- sweep(chart$covariances, c(1L, 2L), outer(scales, scales), "*")
  }
  if (!is.null(chart$direct_effects)) {
    chart$direct_effects[] <- scales[extra$direct$effects$column] /
      extra$direct$scales
  }
  units <- .multilpa_flatten_parameters(chart)
  list(flatten = function(point) .multilpa_flatten_parameters(point) / units,
       restore = function(point, values) {
         .multilpa_restore_parameters(point, values * units)
       },
       valid = .multilpa_valid_parameters)
}

#' Refill a parameter list's updatable blocks from a flat vector
#'
#' The inverse of `.multilpa_flatten_parameters()`: every array, matrix and
#' nested list keeps its shape and names; only the numbers change. Blocks
#' outside the updatable set (a noise density) are carried over unchanged.
#'
#' @param template A parameter list with the target shapes.
#' @param values A numeric vector of the flattened length.
#' @param blocks Names of the blocks to restore.
#' @return A parameter list.
#' @noRd
.multilpa_restore_parameters <- function(template, values,
                                         blocks = .multilpa_accelerated_blocks()) {
  position <- 0L
  fill <- function(block) {
    if (is.list(block)) return(lapply(block, fill))
    size <- length(block)
    block[] <- values[position + seq_len(size)]
    position <<- position + size
    block
  }
  names_present <- intersect(blocks, names(template))
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
                                       "initial_probabilities",
                                       "transition_probabilities",
                                       "response_probabilities")],
                          use.names = FALSE)
  if (any(!is.finite(probabilities)) || any(probabilities <= 0) ||
      any(probabilities > 1 + 1e-12)) {
    return(FALSE)
  }
  normalized <- function(totals) all(is.finite(totals) & abs(totals - 1) <= 1e-8)
  if (!is.null(parameters$profile_probabilities) &&
      !normalized(rowSums(parameters$profile_probabilities))) return(FALSE)
  if (!is.null(parameters$response_probabilities) &&
      !all(vapply(parameters$response_probabilities, function(block) {
        normalized(rowSums(block))
      }, logical(1)))) return(FALSE)
  if (!is.null(parameters$initial_probabilities) &&
      !normalized(rowSums(parameters$initial_probabilities))) return(FALSE)
  if (!is.null(parameters$transition_probabilities) &&
      !normalized(apply(parameters$transition_probabilities, c(1L, 3L), sum))) {
    return(FALSE)
  }
  if (!is.null(parameters$group_probabilities) &&
      !normalized(sum(parameters$group_probabilities))) return(FALSE)
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


#' SQUAREM maps for homogeneous probabilities or general transition logits
#'
#' General states carry one measurement block per measurement occasion and
#' unconstrained initial/transition logits. The fixed stayer mask is retained
#' by restore and never extrapolated.
#' @param acceleration SQUAREM or plain EM.
#' @param general Whether the state uses logits and nested measurement blocks.
#' @return NULL or the flatten, restore and valid functions for the EM kernel.
#' @noRd
.lta_accelerator <- function(acceleration, general = FALSE) {
  if (identical(acceleration, "none")) return(NULL)
  if (!general) return(list(flatten = .multilpa_flatten_parameters,
                           restore = .multilpa_restore_parameters,
                           valid = .multilpa_valid_parameters))
  blocks <- c("measurement", "initial", "transition", "transition2",
              "group_probabilities")
  list(flatten = function(point) .multilpa_flatten_parameters(point, blocks),
       restore = function(point, values) .multilpa_restore_parameters(point, values, blocks),
       valid = function(point) {
         all(is.finite(.multilpa_flatten_parameters(point, blocks))) &&
           .multilpa_valid_parameters(list(group_probabilities = point$group_probabilities)) &&
           all(vapply(point$measurement, .multilpa_valid_parameters, logical(1)))
       })
}
