# Wald inference for the general transition model (multilpa_lta). The
# estimation scale: per measurement block, profile means and log variances;
# group-class logits (baseline: class 1); initial-distribution logit
# coefficients (reference: last profile); transition logit coefficients
# (reference: staying). By Fisher's identity each group's score is the
# expected complete-data score under its posteriors.

#' Free design columns of a transition coefficient block
#' @param terms Design column names.
#' @param second Whether second-order transitions are fitted.
#' @param block First- or second-order transition block.
#' @param n_occasion_intercepts Number of generated occasion intercepts.
#' @return Integer indices excluding occasion intercepts that never apply.
#' @noRd
.lta_transition_columns <- function(terms, second, block = 1L, n_occasion_intercepts = NULL) {
  if (!second || !identical(terms[1L], "occasion_2")) return(seq_along(terms))
  if (is.null(n_occasion_intercepts)) {
    generated <- terms == paste0("occasion_", seq_along(terms) + 1L)
    n_occasion_intercepts <- match(FALSE, generated, nomatch = length(terms) + 1L) - 1L
  }
  occasion <- seq_along(terms) <= n_occasion_intercepts
  which(!occasion | if (block == 1L) terms == "occasion_2" else terms != "occasion_2")
}

#' Does LTA measurement need a covariance chart rather than log variances?
#' @param fit A fit or parameter specification.
#' @return A single logical.
#' @noRd
.lta_uses_covariance_chart <- function(fit) {
  length(fit$continuous) > 0L &&
    (!identical(fit$covariance_model %||% "diagonal", "diagonal") ||
       !is.null(fit$covariance_structure) && !fit$covariance_structure %in% c("EEI", "VVI"))
}

#' Core measurement view for covariance coordinates of one LTA occasion
#' @param block Measurement parameters of one occasion.
#' @param fit Fit or parameter specification carrying its covariance structure.
#' @return A core model view with no categorical or mixing free blocks.
#' @noRd
.lta_covariance_view <- function(block, fit) {
  k <- nrow(block$means)
  structure(list(means = block$means, variances = block$variances,
    covariances = block$covariances, continuous = fit$continuous, vars = fit$continuous,
    covariance_model = fit$covariance_model, covariance_structure = fit$covariance_structure,
    variance_model = fit$variance_model, n_profiles = k, n_group_classes = 1L,
    profile_probabilities = matrix(1 / k, 1L, k), group_probabilities = 1),
    class = "multilpa")
}

#' Covariance coordinates of a single LTA measurement block
#' @param block Measurement parameters.
#' @param fit Fit or parameter specification.
#' @return Named unconstrained covariance coordinates from the core charts.
#' @noRd
.lta_covariance_coordinates <- function(block, fit) {
  view <- .lta_covariance_view(block, fit)
  if (.multilpa_uses_chart(view)) return(.multilpa_chart_encode(view)$values)
  .multilpa_covariance_coordinates(view, "unconstrained")
}

#' Names of a single LTA covariance coordinate block
#' @param fit Fit carrying the measurement blocks.
#' @param block Index of its measurement block.
#' @param label Occasion label.
#' @return Names in the LTA coefficient grammar.
#' @noRd
.lta_covariance_names <- function(fit, block, label) {
  view <- .lta_covariance_view(fit$measurement[[block]], fit)
  coordinates <- .lta_covariance_coordinates(fit$measurement[[block]], fit)
  labels <- .multilpa_tidy_labels(names(coordinates), view)
  stem <- paste(labels$parameter, label, labels$outcome, sep = ".")
  ifelse(is.na(labels$term), stem, paste(stem, labels$term, sep = "."))
}

#' Names of the estimation-scale coefficients
#' @noRd
.lta_names <- function(fit) {
  profiles <- paste0("profile_", seq_len(fit$n_profiles))
  classes <- paste0("group_class_", seq_len(fit$n_group_classes))
  classes[fit$stayer %||% rep(FALSE, length(classes))] <- "stayers"
  blocks <- if (length(fit$measurement) == 1L) "all" else
    paste0("occasion_", seq_along(fit$measurement))
  shared <- identical(fit$variance_model, "equal")
  categorical <- fit$categorical
  ## One rule for every label (`.latents_term_component()`): a term naming one
  ## label -- the indicator of a mean or a log variance -- is that label's own
  ## text, and a term pairing two -- a covariance entry, from
  ## `.lta_covariance_names()` -- is spelled by `.latents_term()`. So `y:2`
  ## has the mean term `y:2` and the diagonal term `y%3A2:y%3A2`, exactly as
  ## `coef()` of the profile and transition models spells them.
  measurement <- unlist(lapply(seq_along(blocks), function(index) {
    b <- blocks[index]
    c(if (length(fit$continuous) > 0L) c(
        paste("mean", b, rep(profiles, each = length(fit$continuous)),
              fit$continuous, sep = "."),
        if (.lta_uses_covariance_chart(fit)) .lta_covariance_names(fit, index, b) else
          if (shared) paste("log_variance", b, fit$continuous, sep = ".") else
          paste("log_variance", b, rep(profiles, each = length(fit$continuous)),
                fit$continuous, sep = ".")),
      unlist(lapply(seq_along(categorical), function(i) {
        unlist(lapply(profiles, function(k) {
          paste("logit", b, k, categorical[i],
                paste0("category_", seq_len(fit$n_categories[[i]] - 1L)), sep = ".")
        }))
      })),
      .lta_extra_names(fit, b, profiles))
  }))
  initial_terms <- colnames(fit$designs$initial)
  transition_terms <- colnames(fit$designs$transition[[1L]])
  n_intercepts <- fit$layout$n_occasions - 1L
  first_terms <- transition_terms[.lta_transition_columns(
    transition_terms, (fit$order %||% 1L) >= 2L, n_occasion_intercepts = n_intercepts)]
  second_terms <- transition_terms[.lta_transition_columns(
    transition_terms, (fit$order %||% 1L) >= 2L, 2L, n_intercepts)]
  initial <- unlist(lapply(classes, function(h) {
    unlist(lapply(profiles[-fit$n_profiles], function(k) {
      paste("initial", h, k, initial_terms, sep = ".")
    }))
  }))
  movers <- classes[!(fit$stayer %||% rep(FALSE, length(classes)))]
  transition <- unlist(lapply(movers, function(h) {
    unlist(lapply(seq_len(fit$n_profiles), function(k) {
      unlist(lapply(setdiff(seq_len(fit$n_profiles), k), function(l) {
        paste("transition", h, paste0(profiles[k], "->", profiles[l]),
              first_terms, sep = ".")
      }))
    }))
  }))
  transition2 <- if ((fit$order %||% 1L) < 2L) character() else {
    unlist(lapply(movers, function(h) {
      unlist(lapply(seq_len(fit$n_profiles^2), function(pair) {
        i <- (pair - 1L) %/% fit$n_profiles + 1L
        j <- (pair - 1L) %% fit$n_profiles + 1L
        unlist(lapply(setdiff(seq_len(fit$n_profiles), j), function(l) {
          paste("transition2", h,
                paste0(profiles[i], "->", profiles[j], "->", profiles[l]),
                second_terms, sep = ".")
        }))
      }))
    }))
  }
  c(measurement, sprintf("logit.%s", classes[-1L]), initial, transition, transition2)
}

#' Pack a parameter list onto the estimation scale
#' @noRd
.lta_pack <- function(parameters, variance_model, extra = NULL, fit = parameters) {
  measurement <- unlist(lapply(parameters$measurement, function(block) {
    c(if (!is.null(block$means) && ncol(block$means) > 0L) c(
        as.vector(t(block$means)),
        if (.lta_uses_covariance_chart(fit)) unname(.lta_covariance_coordinates(block, fit)) else
          if (identical(variance_model, "equal")) log(block$variances[1L, ]) else
          log(as.vector(t(block$variances)))),
      # Response probabilities as logits against the last category, per
      # indicator, profile-major.
      unlist(lapply(block$response_probabilities, function(probabilities) {
        free <- seq_len(ncol(probabilities) - 1L)
        as.vector(t(log(probabilities[, free, drop = FALSE] /
                          probabilities[, ncol(probabilities)])))
      })),
      if (!is.null(extra)) unname(.latents_extra_coordinates(
        block, extra, nrow(block$means), "unconstrained")))
  }))
  weights <- parameters$group_probabilities
  # Initial: per class, per free profile, the design terms.
  initial <- unlist(lapply(seq_len(dim(parameters$initial)[3L]), function(h) {
    as.vector(matrix(parameters$initial[, , h], dim(parameters$initial)[1L]))
  }))
  movers <- which(!(parameters$stayer %||% rep(FALSE, dim(parameters$transition)[4L])))
  terms <- dimnames(parameters$transition)[[1L]]
  n_intercepts <- attr(parameters$transition, "n_occasion_intercepts")
  first_columns <- if (is.null(terms)) seq_len(dim(parameters$transition)[1L]) else
    .lta_transition_columns(terms, !is.null(parameters$transition2),
                            n_occasion_intercepts = n_intercepts)
  second_columns <- if (is.null(terms)) seq_len(dim(parameters$transition)[1L]) else
    .lta_transition_columns(terms, !is.null(parameters$transition2), 2L, n_intercepts)
  transition <- unlist(lapply(movers, function(h) {
    unlist(lapply(seq_len(dim(parameters$transition)[3L]), function(k) {
      as.vector(matrix(parameters$transition[first_columns, , k, h], length(first_columns)))
    }))
  }))
  transition2 <- if (is.null(parameters$transition2)) NULL else
    unlist(lapply(movers, function(h) {
      unlist(lapply(seq_len(dim(parameters$transition2)[3L]), function(pair) {
        as.vector(matrix(parameters$transition2[second_columns, , pair, h],
                         length(second_columns)))
      }))
    }))
  unname(c(measurement, log(weights[-1L] / weights[1L]), initial, transition,
           transition2))
}

#' Unpack the estimation-scale vector into a parameter list
#' @noRd
.lta_unpack <- function(theta, fit) {
  n_states <- fit$n_profiles
  d <- length(fit$continuous)
  n_classes <- fit$n_group_classes
  p0 <- ncol(fit$designs$initial)
  p <- ncol(fit$designs$transition[[1L]])
  shared <- identical(fit$variance_model, "equal")
  position <- 0L
  take <- function(n) {
    values <- theta[position + seq_len(n)]
    position <<- position + n
    values
  }
  # The chart coordinates of each block are kept, so its scores can be taken
  # with respect to exactly these coordinates (.lta_chart_scores()).
  chart_coordinates <- vector("list", length(fit$measurement))
  measurement <- lapply(seq_along(fit$measurement), function(b) {
    block <- list()
    if (d > 0L) {
      block$means <- matrix(take(n_states * d), n_states, d, byrow = TRUE)
      if (.lta_uses_covariance_chart(fit)) {
        view <- .lta_covariance_view(fit$measurement[[b]], fit)
        width <- .multilpa_coordinate_widths(view, "unconstrained")[["variances"]]
        coordinates <- take(width)
        chart_coordinates[[b]] <<- coordinates
        decoded <- .multilpa_decode(c(as.vector(t(block$means)), coordinates,
                                      numeric(n_states - 1L)), view)
        block$variances <- decoded$variances
        if (!is.null(decoded$covariances)) block$covariances <- decoded$covariances
      } else {
        block$variances <- if (shared) matrix(exp(take(d)), n_states, d, byrow = TRUE) else
          matrix(exp(take(n_states * d)), n_states, d, byrow = TRUE)
      }
    } else {
      block$means <- matrix(0, n_states, 0L)
      block$variances <- matrix(1, n_states, 0L)
    }
    if (length(fit$n_categories) > 0L) {
      block$response_probabilities <- lapply(fit$n_categories, function(n_cat) {
        logits <- cbind(matrix(take(n_states * (n_cat - 1L)), n_states, n_cat - 1L, byrow = TRUE), 0)
        probabilities <- exp(logits - .multilpa_row_max(logits))
        probabilities / rowSums(probabilities)
      })
    }
    if (!is.null(fit$extra_data)) {
      decoded <- .latents_extra_decode(
        take(.latents_extra_width(fit$extra_data, n_states)), fit$extra_data, n_states)
      block[names(decoded)] <- decoded
    }
    block
  })
  logits <- c(0, take(n_classes - 1L))
  initial <- array(take(p0 * (n_states - 1L) * n_classes), c(p0, n_states - 1L, n_classes))
  stayer <- fit$stayer %||% rep(FALSE, n_classes)
  movers <- which(!stayer)
  # Stayer classes carry no free transition coefficients (fixed identity).
  transition <- array(0, c(p, n_states - 1L, n_states, n_classes))
  terms <- colnames(fit$designs$transition[[1L]])
  n_intercepts <- fit$layout$n_occasions - 1L
  first_columns <- .lta_transition_columns(terms, (fit$order %||% 1L) >= 2L,
                                           n_occasion_intercepts = n_intercepts)
  second_columns <- .lta_transition_columns(terms, (fit$order %||% 1L) >= 2L, 2L, n_intercepts)
  transition[first_columns, , , movers] <- array(
    take(length(first_columns) * (n_states - 1L) * n_states * length(movers)),
    c(length(first_columns), n_states - 1L, n_states, length(movers)))
  dimnames(transition) <- list(terms, NULL, NULL, NULL)
  attr(transition, "n_occasion_intercepts") <- n_intercepts
  transition2 <- if ((fit$order %||% 1L) < 2L) NULL else {
    out <- array(0, c(p, n_states - 1L, n_states^2, n_classes))
    out[second_columns, , , movers] <- array(
      take(length(second_columns) * (n_states - 1L) * n_states^2 * length(movers)),
      c(length(second_columns), n_states - 1L, n_states^2, length(movers)))
    out
  }
  c(list(measurement = measurement, initial = initial, transition = transition),
    if (!is.null(transition2)) list(transition2 = transition2),
    list(group_probabilities = exp(logits - max(logits)) /
           sum(exp(logits - max(logits))), stayer = unname(stayer),
         continuous = fit$continuous, variance_model = fit$variance_model,
         covariance_model = fit$covariance_model,
         covariance_structure = fit$covariance_structure,
         chart_coordinates = chart_coordinates))
}

#' Per-group scores of the general transition model on the estimation scale
#' @noRd
.lta_group_scores <- function(theta, fit) {
  parameters <- .lta_unpack(theta, fit)
  x <- fit$x
  e <- .lta_expectation(x, fit$codes, fit$layout, fit$designs, parameters,
                        fit$occasion_of_row, unname(fit$sampling_weights),
                        fit$extra_data)
  n_states <- fit$n_profiles
  n_classes <- fit$n_group_classes
  n_units <- fit$n_groups
  shared <- identical(fit$variance_model, "equal")
  group_index <- fit$group_index
  measurement <- do.call(cbind, lapply(seq_along(parameters$measurement), function(b) {
    rows <- if (length(parameters$measurement) == 1L) seq_len(nrow(x)) else
      which(fit$occasion_of_row == b)
    block <- parameters$measurement[[b]]
    observed <- !is.na(x[rows, , drop = FALSE])
    response <- if (is.null(fit$codes)) NULL else {
      by_group <- .multilpa_response_scores(fit$codes[rows, , drop = FALSE],
                                            e$subject_posteriors[rows, , drop = FALSE],
                                            block$response_probabilities,
                                            group_index[rows])
      out <- matrix(0, n_units, ncol(by_group))
      out[unique(group_index[rows]), ] <- by_group
      out
    }
    extra <- if (is.null(fit$extra_data)) NULL else {
      by_group <- .latents_extra_scores(.latents_extra_rows(fit$extra_data, rows),
                                        e$subject_posteriors[rows, , drop = FALSE],
                                        block, group_index[rows])
      out <- matrix(0, n_units, ncol(by_group))
      out[unique(group_index[rows]), ] <- by_group
      out
    }
    if (ncol(x) == 0L) return(cbind(response, extra))
    if (.lta_uses_covariance_chart(fit)) {
      return(cbind(.lta_chart_scores(fit, b, block, parameters$chart_coordinates[[b]],
                                     x[rows, , drop = FALSE], e, rows, group_index,
                                     n_units),
                   response, extra))
    }
    pieces <- lapply(seq_len(n_states), function(k) {
      residuals <- sweep(x[rows, , drop = FALSE], 2L, block$means[k, ], "-")
      # A missing indicator contributes no score (it is integrated out).
      residuals[!observed] <- 0
      weight <- e$subject_posteriors[rows, k]
      list(mean = .lta_rowsum(sweep(residuals, 2L, block$variances[k, ], "/") * weight,
                              group_index[rows], n_units),
           log_variance = .lta_rowsum(0.5 * weight *
                                        (sweep(residuals^2, 2L, block$variances[k, ], "/") -
                                           observed),
                                      group_index[rows], n_units))
    })
    cbind(do.call(cbind, lapply(pieces, `[[`, "mean")),
          if (shared) Reduce(`+`, lapply(pieces, `[[`, "log_variance")) else
            do.call(cbind, lapply(pieces, `[[`, "log_variance")),
          response, extra)
  }))
  # A weighted posterior row sums to the person's weight, so the expected
  # share is that weight times the class probability; unweighted, it is one.
  logits <- e$group_posteriors[, -1L, drop = FALSE] -
    outer(rowSums(e$group_posteriors), parameters$group_probabilities[-1L])
  initial <- do.call(cbind, lapply(seq_len(n_classes), function(h) {
    counts <- e$moments[[h]]$initial
    fitted <- exp(e$per_class[[h]]$log_initial) * rowSums(counts)
    residual <- (counts - fitted)[, seq_len(n_states - 1L), drop = FALSE]
    do.call(cbind, lapply(seq_len(n_states - 1L), function(l) {
      residual[, l] * fit$designs$initial
    }))
  }))
  first_columns <- .lta_transition_columns(
    colnames(fit$designs$transition[[1L]]), (fit$order %||% 1L) >= 2L,
    n_occasion_intercepts = fit$layout$n_occasions - 1L)
  second_columns <- .lta_transition_columns(
    colnames(fit$designs$transition[[1L]]), (fit$order %||% 1L) >= 2L, 2L,
    fit$layout$n_occasions - 1L)
  movers <- which(!(fit$stayer %||% rep(FALSE, n_classes)))
  transition <- do.call(cbind, lapply(movers, function(h) {
    pairs <- e$moments[[h]]$pairs
    log_transition <- e$per_class[[h]]$log_transition
    do.call(cbind, lapply(seq_len(n_states), function(k) {
      others <- setdiff(seq_len(n_states), k)
      do.call(cbind, lapply(others, function(l) {
        Reduce(`+`, lapply(seq_along(pairs), function(t) {
          moves <- matrix(pairs[[t]][, k, ], n_units)
          residual <- moves[, l] - rowSums(moves) * exp(log_transition[[t]][, k, l])
          residual * fit$designs$transition[[t]][, first_columns, drop = FALSE]
        }))
      }))
    }))
  }))
  transition2 <- if (is.null(parameters$transition2)) NULL else {
    do.call(cbind, lapply(movers, function(h) {
      moves <- e$moments[[h]]$pairs2
      log_p <- e$per_class[[h]]$log_transition2
      do.call(cbind, lapply(seq_len(n_states^2), function(pair) {
        j <- (pair - 1L) %% n_states + 1L
        do.call(cbind, lapply(setdiff(seq_len(n_states), j), function(l) {
          Reduce(`+`, lapply(seq_along(moves), function(t) {
            counts <- matrix(moves[[t]][, pair, ], n_units)
            residual <- counts[, l] - rowSums(counts) * exp(log_p[[t]][, pair, l])
            residual * fit$designs$transition[[t + 1L]][, second_columns, drop = FALSE]
          }))
        }))
      }))
    }))
  }
  scores <- cbind(measurement, logits, initial, transition, transition2)
  dimnames(scores) <- list(fit$group_ids, .lta_names(fit))
  scores
}

#' Per-group scores of one measurement block in its covariance chart
#'
#' The gradient of the log likelihood with respect to each profile's mean and
#' unconstrained covariance (`.multilpa_profile_gradients()`, which also
#' handles indicators missing at random) is carried to the block's chart
#' coordinates by the chain rule, with the Jacobian of the coordinates-to-
#' covariance map. That map involves no data, so its Jacobian is taken by
#' central differences of the same decoding the parameters came from; every
#' structure, separable or constrained, is scored the same way.
#'
#' @param fit The transition fit (its blocks anchor the chart).
#' @param b The block's index.
#' @param block The block's decoded parameters.
#' @param phi Its chart coordinates.
#' @param x The block's rows of the indicator matrix.
#' @param e The expectation at these parameters.
#' @param rows The block's rows.
#' @param group_index,n_units Each row's group and the number of groups.
#' @return A groups-by-coordinates matrix: means, then the chart.
#' @noRd
.lta_chart_scores <- function(fit, b, block, phi, x, e, rows, group_index, n_units) {
  view <- .lta_covariance_view(fit$measurement[[b]], fit)
  k <- nrow(block$means)
  d <- ncol(x)
  covariance_of <- function(coordinates) {
    decoded <- .multilpa_decode(c(as.vector(t(block$means)), coordinates,
                                  numeric(k - 1L)), view)
    vapply(seq_len(k), function(profile) {
      as.vector(if (is.null(decoded$covariances)) diag(decoded$variances[profile, ], d) else
        decoded$covariances[, , profile])
    }, numeric(d * d))
  }
  columns <- lapply(seq_along(phi), function(j) {
    h <- 1e-6 * max(1, abs(phi[j]))
    (covariance_of(replace(phi, j, phi[j] + h)) -
       covariance_of(replace(phi, j, phi[j] - h))) / (2 * h)
  })
  jacobian <- lapply(seq_len(k), function(profile) {
    matrix(vapply(columns, function(column) column[, profile], numeric(d * d)),
           d * d, length(phi))
  })
  expectation <- list(subject_posteriors = e$subject_posteriors[rows, , drop = FALSE],
                      gaussian_moments = e$gaussian_moments[[b]])
  gradients <- .multilpa_profile_gradients(block, expectation, x, group_index[rows])
  by_group <- cbind(gradients$means,
                    Reduce(`+`, Map(`%*%`, gradients$covariance, jacobian)))
  out <- matrix(0, n_units, ncol(by_group))
  out[unique(group_index[rows]), ] <- by_group
  out
}

#' rowsum() into a fixed J rows (groups with no rows get zeros)
#' @noRd
.lta_rowsum <- function(values, groups, n_groups) {
  out <- matrix(0, n_groups, ncol(values))
  summed <- rowsum(values, groups, reorder = TRUE)
  out[as.integer(rownames(summed)), ] <- summed
  out
}

#' Local parameter scales for transition-model score differences
#'
#' Means and unlogged Cholesky entries use their measurement SDs; log,
#' probability and transition coordinates use unit scale. Constrained
#' covariance coordinates use the same chart scales as profile inference.
#' These cap a score-derived step when a nearly flat coordinate has a huge
#' inverse score norm, without introducing an absolute outcome-unit floor.
#' @param fit A general transition fit.
#' @return Named positive scales in `.lta_names()` order.
#' @noRd
.lta_inference_scales <- function(fit) {
  scales <- stats::setNames(rep(1, length(.lta_names(fit))), .lta_names(fit))
  if (!length(fit$continuous)) return(scales)
  invisible(lapply(seq_along(fit$measurement), function(b) {
    block <- fit$measurement[[b]]
    label <- if (length(fit$measurement) == 1L) "all" else paste0("occasion_", b)
    means <- startsWith(names(scales), paste0("mean.", label, "."))
    scales[means] <<- as.vector(t(sqrt(block$variances)))
    if (.lta_uses_covariance_chart(fit)) {
      view <- .lta_covariance_view(block, fit)
      covariance_scale <- if (.multilpa_uses_chart(view)) {
        .multilpa_chart_encode(view)$scale
      } else .multilpa_covariance_coordinate_scale(
        block$variances, length(fit$continuous),
        if (identical(fit$variance_model, "equal")) 1L else seq_len(fit$n_profiles))
      scales[.lta_covariance_names(fit, b, label)] <<- covariance_scale
    }
  }))
  scales
}

#' Estimation-scale covariance of a general transition fit
#' @noRd
.lta_inference <- function(fit, vcov_type = c("observed", "robust", "opg"),
                           step = 1e-4) {
  vcov_type <- match.arg(vcov_type)
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) || step <= 0) {
    stop(errorCondition("`step` must be a positive finite number.",
                        class = "latents_bad_argument", call = NULL))
  }
  if (any(vapply(fit$measurement, function(block) {
    any(block$count_dispersion <= .latents_min_dispersion * (1 + 1e-6)) ||
      any(block$count_means <= 1e-8)
  }, logical(1)))) {
    stop(errorCondition("A count parameter is at its boundary; Wald inference does not apply.",
                        class = "latents_boundary_fit", call = NULL))
  }
  if (any(vapply(fit$measurement, function(block) {
    any(unlist(block$response_probabilities) <= 1e-6)
  }, logical(1)))) {
    stop(errorCondition(paste(
      "A response probability sits at its bound; Wald inference does not apply."),
      class = "latents_boundary_fit", call = NULL))
  }
  if (!isTRUE(fit$converged)) {
    stop(errorCondition("Inference requires a converged fit.",
                        class = "latents_no_converge", call = NULL))
  }
  if (isTRUE(fit$boundary)) {
    stop(errorCondition(paste(
      "A variance or covariance eigenvalue sits at min_variance, or a transition",
      "or initial probability at zero; Wald inference does not apply."),
      class = "latents_boundary_fit", call = NULL))
  }
  theta <- .lta_theta(fit)
  names(theta) <- .lta_names(fit)
  # Only the identified coordinates are differenced and inverted: a flat or
  # product-only logit, or the measurement block of an empty occasion, would
  # make the information singular for every coefficient (R2-LTA-2). The
  # others get NA rows and columns.
  identified <- .lta_identified_coordinates(fit)
  n_units <- sum(.multilpa_inference_groups(fit))
  if (vcov_type != "observed") {
    .multilpa_require_clusters(sum(.multilpa_inference_groups(fit, observed_only = TRUE)),
                               sum(identified),
                               "A robust or OPG transition covariance")
  }
  # The information is differenced, inverted and judged for singularity in the
  # chart the fit was searched in, with covariates centred and scaled
  # (`.lta_design_chart()`): in the supplied units a covariate held at 1e8 makes
  # the information nearly singular and the differences meaningless. The
  # results are mapped back exactly: theta = J phi, so the covariance is
  # J V J', the information A' I A and the scores S A, with A = J^-1 (LTA-F5).
  # J moves an intercept by its covariates' slopes and nothing else, and a
  # covariate column is uncounted only with its whole block, so the identified
  # block of J is closed: J and J^-1 restrict to it exactly.
  used <- .lta_used_rows(fit$layout, fit$x, fit$codes, fit$extra_data,
                         unname(fit$sampling_weights))
  chart <- .lta_design_chart(fit$designs, fit$layout, used)
  internal <- fit
  internal$designs <- .lta_chart_designs(fit$designs, chart)
  phi <- .lta_pack(.lta_map_coefficients(.lta_parameters(fit), chart, "to_chart"),
                   fit$variance_model, fit$extra_data)
  to_raw <- .lta_chart_jacobian(fit, chart, length(theta), "to_raw")
  to_chart <- .lta_chart_jacobian(fit, chart, length(theta), "to_chart")
  all_scores <- .lta_group_scores(phi, internal)
  chart_scores <- all_scores[, identified, drop = FALSE]
  # Differentiate displacements measured in each coordinate's score-based
  # standard-error scale. Covariate standardization alone leaves measurement
  # means in their supplied units: a fixed 1e-4 step can exceed a whole SD,
  # and comparing raw information eigenvalues can falsely declare singularity.
  # With phi = phi_hat + D u, scores are S D and information is D I D.
  # An all-zero score column still receives a finite scale so genuine missing
  # information is judged by the usual positive-definiteness/rank checks.
  parameter_scale <- 1 / sqrt(colSums(chart_scores^2))
  parameter_scale[!is.finite(parameter_scale) | parameter_scale <= 0] <- 1
  # Conditioning and derivative widths have different jobs: a nearly flat
  # logit can have a huge score-based scale, but differencing it must remain
  # local. In displacement coordinates this ratio caps its step; mean and
  # Cholesky caps retain their measurement units.
  derivative_scale <- pmin(parameter_scale, .lta_inference_scales(fit)[identified])
  scaled_scores <- sweep(chart_scores, 2L, parameter_scale, "*")
  total <- function(displacement) {
    point <- phi[identified] + displacement * parameter_scale
    colSums(.lta_group_scores(replace(phi, identified, point), internal))[identified] *
      parameter_scale
  }
  scaled_information <- .inference_score_information(
    total, numeric(sum(identified)), step, scheme = "five_point", h_rule = "scaled",
    scale = derivative_scale / parameter_scale)
  inverse <- .inference_invert(scaled_information, "min_le_eps_max")
  meat <- crossprod(scaled_scores)
  if (vcov_type != "observed") {
    .inference_require_rank(scaled_scores[
      .multilpa_inference_groups(fit, observed_only = TRUE), , drop = FALSE],
      "Too few groups for a robust or OPG covariance.")
  }
  scaled_covariance <- switch(vcov_type, observed = inverse,
                              robust = .inference_sandwich(inverse, meat, n_units),
                              opg = chol2inv(chol(meat)))
  chart_covariance <- scaled_covariance * tcrossprod(parameter_scale)
  chart_information <- scaled_information / tcrossprod(parameter_scale)
  kept_raw <- to_raw[identified, identified, drop = FALSE]
  kept_chart <- to_chart[identified, identified, drop = FALSE]
  covariance <- information <- matrix(NA_real_, length(theta), length(theta))
  covariance[identified, identified] <- kept_raw %*% chart_covariance %*% t(kept_raw)
  information[identified, identified] <- crossprod(kept_chart,
                                                   chart_information %*% kept_chart)
  # A flat coordinate's score is zero in every group.
  scores <- all_scores %*% to_chart
  dimnames(information) <- dimnames(covariance) <- list(names(theta), names(theta))
  colnames(scores) <- colnames(all_scores)
  list(theta = theta, vcov = covariance, information = information,
       group_scores = scores, gradient = colSums(scores), vcov_type = vcov_type,
       identified = stats::setNames(identified, names(theta)))
}

#' Which estimation-scale coordinates of a general transition fit are identified
#'
#' The logits of `.lta_identified_columns()` and the measurement blocks of
#' occasions at which something is observed; the group-class and initial
#' logits always are. The rest are counted in neither `n_parameters` nor the
#' information, and their standard errors are `NA`.
#'
#' @param fit A `multilpa_lta` fit.
#' @return A logical vector in `.lta_names()` order.
#' @noRd
.lta_identified_coordinates <- function(fit) {
  names_all <- .lta_names(fit)
  occupied <- fit$occupied %||% rep(TRUE, fit$layout$n_occasions)
  terms <- colnames(fit$designs$transition[[1L]])
  counted <- .lta_identified_columns(terms, fit$layout$n_occasions,
                                     fit$transitions %||% "homogeneous",
                                     fit$order %||% 1L, occupied)
  # Names are kind.owner.outcome.term; the owner of a measurement coordinate
  # is its block ("all" or "occasion_t"), of a logit its group class.
  pieces <- strsplit(names_all, ".", fixed = TRUE)
  kind <- vapply(pieces, `[`, character(1), 1L)
  owner <- vapply(pieces, function(p) if (length(p) >= 2L) p[[2L]] else "", character(1))
  term <- vapply(pieces, function(p) paste(p[-(1:3)], collapse = "."), character(1))
  identified <- rep(TRUE, length(names_all))
  by_occasion <- !kind %in% c("initial", "transition", "transition2") &
    grepl("^occasion_[0-9]+$", owner)
  identified[by_occasion] <- occupied[as.integer(sub("^occasion_", "", owner[by_occasion]))]
  first <- kind == "transition"
  second <- kind == "transition2"
  identified[first] <- term[first] %in% terms[counted$first]
  identified[second] <- term[second] %in% terms[counted$second]
  # The logits out of an origin no expected move leaves are flat
  # (.lta_empty_origins()); a fit saved before 0.97.5 records none.
  identified[(first | second) & names_all %in% .lta_empty_origin_names(fit)] <- FALSE
  identified
}

#' Estimation-scale names of the logits out of empty transition origins
#' @param fit A `multilpa_lta` fit.
#' @return A character vector of `.lta_names()` entries, possibly empty.
#' @noRd
.lta_empty_origin_names <- function(fit) {
  profiles <- paste0("profile_", seq_len(fit$n_profiles))
  name_block <- function(empty, kind, origin_label, origins_of) {
    if (is.null(empty) || !any(empty)) return(character())
    cells <- which(empty, arr.ind = TRUE)
    classes <- colnames(empty)
    unlist(lapply(seq_len(nrow(cells)), function(i) {
      origin <- cells[i, 1L]
      destinations <- setdiff(seq_len(fit$n_profiles), origins_of(origin))
      sprintf("%s.%s.%s->%s.", kind, classes[cells[i, 2L]], origin_label(origin),
              profiles[destinations])
    }))
  }
  stems <- c(
    name_block(fit$empty_transition_rows, "transition", function(k) profiles[k],
               function(k) k),
    name_block(fit$empty_second_order_rows, "transition2", function(pair) {
      paste0(profiles[(pair - 1L) %/% fit$n_profiles + 1L], "->",
             profiles[(pair - 1L) %% fit$n_profiles + 1L])
    }, function(pair) (pair - 1L) %% fit$n_profiles + 1L))
  names_all <- .lta_names(fit)
  names_all[Reduce(`|`, lapply(stems, startsWith, x = names_all), logical(length(names_all)))]
}

#' A fit's estimates in the engine's coordinates
#'
#' The packed vector the likelihood, scores and EM work with: continuous means
#' centred on `fit$centers`. [coef()] reports the same vector with the means
#' in the data's units.
#' @param fit A general transition fit (`multilpa_lta`).
#' @return A numeric vector in [.lta_names()] order.
#' @noRd
.lta_theta <- function(fit) {
  .lta_pack(.lta_parameters(fit), fit$variance_model, fit$extra_data)
}

#' The fit's parameters on the centred scale the engine works on
#' @noRd
.lta_parameters <- function(fit) {
  list(measurement = lapply(fit$measurement, function(block) {
         c(list(means = if (length(fit$continuous) > 0L)
                  unname(sweep(block$means, 2L, fit$centers, "-")) else
                    matrix(0, fit$n_profiles, 0L),
                variances = if (length(fit$continuous) > 0L) unname(block$variances) else
                  matrix(1, fit$n_profiles, 0L)),
           if (!is.null(block$covariances)) list(covariances = unname(block$covariances)),
           if (!is.null(block$response_probabilities))
             list(response_probabilities = lapply(block$response_probabilities, unname)),
           if (!is.null(block$count_means)) list(count_means = unname(block$count_means)),
           if (!is.null(block$count_dispersion))
             list(count_dispersion = unname(block$count_dispersion)),
           if (!is.null(block$ordinal_intercepts))
             list(ordinal_intercepts = lapply(block$ordinal_intercepts, unname),
                  ordinal_locations = unname(block$ordinal_locations)))
       }),
       initial = unname(fit$initial_coefficients),
       transition = fit$transition_coefficients,
       transition2 = if (is.null(fit$second_order_coefficients)) NULL else
         unname(fit$second_order_coefficients),
       group_probabilities = unname(fit$group_probabilities),
       stayer = unname(fit$stayer %||% rep(FALSE, length(fit$group_probabilities))),
       continuous = fit$continuous, variance_model = fit$variance_model,
       covariance_model = fit$covariance_model, covariance_structure = fit$covariance_structure)
}

#' Finish EM with a bounded quasi-Newton search on the exact likelihood
#'
#' EM converges slowly when profiles are weakly separated. From the EM point,
#' L-BFGS-B maximizes the log likelihood with the analytic gradient, keeping
#' log variances at or above log(min_variance); the point is kept only if the
#' likelihood rises.
#' @param spec The fit-shaped specification the pack/unpack and score helpers
#'   read.
#' @return A list with `parameters`, `log_likelihood`, `gradient` (max abs)
#'   and `converged`.
#' @noRd
.lta_quasi_newton <- function(spec, parameters, min_variance, tol) {
  theta <- .lta_pack(parameters, spec$variance_model, spec$extra_data)
  names_all <- .lta_names(spec)
  lower <- ifelse(startsWith(names_all, "log_variance."), log(min_variance), -Inf)
  lower[startsWith(names_all, "log_count_mean.")] <- log(1e-10)
  lower[startsWith(names_all, "log_count_dispersion.")] <- log(.latents_min_dispersion)
  value <- function(v) {
    -.lta_expectation(spec$x, spec$codes, spec$layout, spec$designs,
                      .lta_unpack(v, spec), spec$occasion_of_row,
                      spec$sampling_weights, spec$extra_data)$log_likelihood
  }
  gradient <- function(v) -colSums(.lta_group_scores(v, spec))
  # A transition never observed has its maximum at probability zero, where the
  # logit runs off and the gradient decays without vanishing; that boundary is
  # flagged separately (see .lta_probability_boundary()). Near it each
  # iteration gains little, so a stopping rule tied to `tol` (factr = tol /
  # eps) stops short: at tol = 1e-10 an LMest occasion model was left 2.8e-6
  # below its maximum with transitions 1e-4 off. The finish is therefore run
  # to L-BFGS-B's own tight rule, as before 0.9.13; `tol` still governs EM and
  # the code-52 check.
  # Means are searched in reference standard deviations, so the path does
  # not depend on the units of the indicators.
  row_weights <- .lta_row_weights(spec$layout, spec$sampling_weights, nrow(spec$x))
  parscale <- rep(1, length(theta))
  means <- startsWith(names_all, "mean.")
  if (any(means)) {
    parscale[means] <- rep(.latents_reference_scales(
      spec$x, row_weights, parameters$measurement[[1L]]$variances),
      length.out = sum(means))
  }
  finish <- .latents_quasi_newton(
    theta, lower, value, gradient, tol,
    list(maxit = 5000L, factr = 10, pgtol = 0),
    scale = .latents_observed_values(spec$x, spec$codes,
                                     .latents_extra_matrix(spec$extra_data),
                                     row_weights = row_weights),
    parscale = parscale)
  if (!finish$improved) {
    return(list(parameters = parameters, log_likelihood = -finish$start_value,
                gradient = max(abs(gradient(finish$par))), converged = FALSE))
  }
  final_gradient <- gradient(finish$par)
  free <- !(finish$par <= lower + 1e-12 & final_gradient > 0)
  list(parameters = .lta_unpack(finish$par, spec), log_likelihood = -finish$value,
       gradient = max(abs(final_gradient[free]), 0),
       converged = finish$converged)
}

#' Does any fitted initial or transition probability sit at zero?
#'
#' Read only where the weighted observed likelihood uses the probability; a probability
#' below 1e-6 marks a move (or start) that never occurs, whose logit has no
#' finite maximum.
#' @noRd
.lta_probability_boundary <- function(expectation, layout, stayer = NULL, used = NULL) {
  used <- used %||% list(initial = rep(TRUE, nrow(layout$slot)), moves = layout$within)
  initial_min <- function(log_initial) {
    if (!any(used$initial)) return(1)
    min(exp(log_initial[used$initial, , drop = FALSE]))
  }
  # A stayer class's zero moves are structural, not estimated: skip them.
  stayer <- stayer %||% rep(FALSE, length(expectation$per_class))
  any(vapply(seq_along(expectation$per_class), function(h) {
    class_terms <- expectation$per_class[[h]]
    if (isTRUE(stayer[h])) {
      return(initial_min(class_terms$log_initial) < 1e-6)
    }
    initial <- initial_min(class_terms$log_initial)
    first_occasions <- if (length(class_terms$log_transition2) > 0L) 1L else
      seq_along(class_terms$log_transition)
    moves <- vapply(first_occasions, function(t) {
      active <- used$moves[, t + 1L]
      if (!any(active)) return(1)
      min(exp(class_terms$log_transition[[t]][active, , , drop = FALSE]))
    }, numeric(1))
    second <- vapply(seq_along(class_terms$log_transition2), function(t) {
      active <- used$moves[, t + 2L]
      if (!any(active)) return(1)
      min(exp(class_terms$log_transition2[[t]][active, , , drop = FALSE]))
    }, numeric(1))
    min(initial, moves, second) < 1e-6
  }, logical(1)))
}

#' Estimation-scale names of a transition fit's ordinal and count coordinates
#' @param fit The fit or its specification.
#' @param block The measurement block label (`all` or `occasion_t`).
#' @param profiles Profile labels.
#' @return A character vector in `.latents_extra_coordinates()` order.
#' @noRd
.lta_extra_names <- function(fit, block, profiles) {
  extra <- fit$extra_data
  if (is.null(extra)) return(character())
  ordinal <- unlist(lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    indicator <- colnames(extra$ordinal)[[j]]
    c(paste("ordinal_intercept", block, indicator, extra$ordinal_levels[[j]][-1L],
            sep = "."),
      paste("ordinal_location", block, profiles[-length(profiles)], indicator,
            sep = "."))
  }))
  dispersion <- if (identical(extra$count_dispersion, "equal")) "shared" else profiles
  count <- unlist(lapply(colnames(extra$count) %||% character(), function(indicator) {
    c(paste("log_count_mean", block, profiles, indicator, sep = "."),
      if (.latents_negative_binomial(extra))
        paste("log_count_dispersion", block, dispersion, indicator, sep = "."))
  }))
  c(ordinal, count)
}
