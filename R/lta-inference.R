# Wald inference for the general transition model (multilpa_lta). The
# estimation scale: per measurement block, profile means and log variances;
# group-class logits (baseline: class 1); initial-distribution logit
# coefficients (reference: last profile); transition logit coefficients
# (reference: staying). By Fisher's identity each group's score is the
# expected complete-data score under its posteriors.

#' Names of the estimation-scale coefficients
#' @noRd
.lta_names <- function(fit) {
  profiles <- paste0("profile_", seq_len(fit$n_profiles))
  classes <- paste0("group_class_", seq_len(fit$n_group_classes))
  blocks <- if (length(fit$measurement) == 1L) "all" else
    paste0("occasion_", seq_along(fit$measurement))
  shared <- identical(fit$variance_model, "equal")
  categorical <- fit$categorical
  measurement <- unlist(lapply(blocks, function(b) {
    c(if (length(fit$continuous) > 0L) c(
        paste("mean", b, rep(profiles, each = length(fit$continuous)),
              fit$continuous, sep = "."),
        if (shared) paste("log_variance", b, fit$continuous, sep = ".") else
          paste("log_variance", b, rep(profiles, each = length(fit$continuous)),
                fit$continuous, sep = ".")),
      unlist(lapply(seq_along(categorical), function(i) {
        unlist(lapply(profiles, function(k) {
          paste("logit", b, k, categorical[i],
                paste0("category_", seq_len(fit$n_categories[[i]] - 1L)), sep = ".")
        }))
      })))
  }))
  initial_terms <- colnames(fit$designs$initial)
  transition_terms <- colnames(fit$designs$transition[[1L]])
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
              transition_terms, sep = ".")
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
                transition_terms, sep = ".")
        }))
      }))
    }))
  }
  c(measurement, sprintf("logit.%s", classes[-1L]), initial, transition, transition2)
}

#' Pack a parameter list onto the estimation scale
#' @noRd
.lta_pack <- function(parameters, variance_model) {
  measurement <- unlist(lapply(parameters$measurement, function(block) {
    c(if (!is.null(block$means) && ncol(block$means) > 0L) c(
        as.vector(t(block$means)),
        if (identical(variance_model, "equal")) log(block$variances[1L, ]) else
          log(as.vector(t(block$variances)))),
      # Response probabilities as logits against the last category, per
      # indicator, profile-major.
      unlist(lapply(block$response_probabilities, function(probabilities) {
        free <- seq_len(ncol(probabilities) - 1L)
        as.vector(t(log(probabilities[, free, drop = FALSE] /
                          probabilities[, ncol(probabilities)])))
      })))
  }))
  weights <- parameters$group_probabilities
  # Initial: per class, per free profile, the design terms.
  initial <- unlist(lapply(seq_len(dim(parameters$initial)[3L]), function(h) {
    as.vector(matrix(parameters$initial[, , h], dim(parameters$initial)[1L]))
  }))
  movers <- which(!(parameters$stayer %||% rep(FALSE, dim(parameters$transition)[4L])))
  transition <- unlist(lapply(movers, function(h) {
    unlist(lapply(seq_len(dim(parameters$transition)[3L]), function(k) {
      as.vector(matrix(parameters$transition[, , k, h], dim(parameters$transition)[1L]))
    }))
  }))
  transition2 <- if (is.null(parameters$transition2)) NULL else
    unlist(lapply(movers, function(h) {
      unlist(lapply(seq_len(dim(parameters$transition2)[3L]), function(pair) {
        as.vector(matrix(parameters$transition2[, , pair, h],
                         dim(parameters$transition2)[1L]))
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
  measurement <- lapply(seq_along(fit$measurement), function(b) {
    block <- list()
    if (d > 0L) {
      block$means <- matrix(take(n_states * d), n_states, d, byrow = TRUE)
      block$variances <- if (shared) matrix(exp(take(d)), n_states, d, byrow = TRUE) else
        matrix(exp(take(n_states * d)), n_states, d, byrow = TRUE)
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
    block
  })
  logits <- c(0, take(n_classes - 1L))
  initial <- array(take(p0 * (n_states - 1L) * n_classes), c(p0, n_states - 1L, n_classes))
  stayer <- fit$stayer %||% rep(FALSE, n_classes)
  movers <- which(!stayer)
  # Stayer classes carry no free transition coefficients (fixed identity).
  transition <- array(0, c(p, n_states - 1L, n_states, n_classes))
  transition[, , , movers] <- array(take(p * (n_states - 1L) * n_states * length(movers)),
                                    c(p, n_states - 1L, n_states, length(movers)))
  transition2 <- if ((fit$order %||% 1L) < 2L) NULL else {
    out <- array(0, c(p, n_states - 1L, n_states^2, n_classes))
    out[, , , movers] <- array(take(p * (n_states - 1L) * n_states^2 * length(movers)),
                               c(p, n_states - 1L, n_states^2, length(movers)))
    out
  }
  c(list(measurement = measurement, initial = initial, transition = transition),
    if (!is.null(transition2)) list(transition2 = transition2),
    list(group_probabilities = exp(logits - max(logits)) /
           sum(exp(logits - max(logits))), stayer = unname(stayer)))
}

#' Per-group scores of the general transition model on the estimation scale
#' @noRd
.lta_group_scores <- function(theta, fit) {
  parameters <- .lta_unpack(theta, fit)
  x <- fit$x
  e <- .lta_expectation(x, fit$codes, fit$layout, fit$designs, parameters,
                        fit$occasion_of_row, unname(fit$sampling_weights))
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
    if (ncol(x) == 0L) return(response)
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
          response)
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
          residual * fit$designs$transition[[t]]
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
            residual * fit$designs$transition[[t + 1L]]
          }))
        }))
      }))
    }))
  }
  scores <- cbind(measurement, logits, initial, transition, transition2)
  dimnames(scores) <- list(fit$group_ids, .lta_names(fit))
  scores
}

#' rowsum() into a fixed J rows (groups with no rows get zeros)
#' @noRd
.lta_rowsum <- function(values, groups, n_groups) {
  out <- matrix(0, n_groups, ncol(values))
  summed <- rowsum(values, groups, reorder = TRUE)
  out[as.integer(rownames(summed)), ] <- summed
  out
}

#' Estimation-scale covariance of a general transition fit
#' @noRd
.lta_inference <- function(fit, vcov_type = c("observed", "robust", "opg"),
                           step = 1e-4) {
  vcov_type <- match.arg(vcov_type)
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
  if (!identical(fit$covariance_model %||% "diagonal", "diagonal") ||
      !is.null(fit$covariance_structure) &&
      !fit$covariance_structure %in% c("EEI", "VVI")) {
    stop(errorCondition(paste(
      "Standard errors for the extended transition model are available for",
      "diagonal covariances (EEI, VVI)."),
      class = "latents_unsupported_inference", call = NULL))
  }
  if (isTRUE(fit$boundary)) {
    stop(errorCondition(paste(
      "A variance sits at min_variance or a transition or initial probability",
      "at zero; Wald inference does not apply."),
      class = "latents_boundary_fit", call = NULL))
  }
  theta <- .lta_pack(.lta_parameters(fit), fit$variance_model)
  names(theta) <- .lta_names(fit)
  scores <- .lta_group_scores(theta, fit)
  total <- function(v) colSums(.lta_group_scores(v, fit))
  columns <- vapply(seq_along(theta), function(k) {
    h <- step * max(1, abs(theta[k]))
    shift <- function(m) {
      point <- theta
      point[k] <- point[k] + m * h
      total(point)
    }
    (-shift(2) + 8 * shift(1) - 8 * shift(-1) + shift(-2)) / (12 * h)
  }, numeric(length(theta)))
  information <- -(columns + t(columns)) / 2
  dimnames(information) <- list(names(theta), names(theta))
  eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(eigenvalues)) || min(eigenvalues) <= 1e-10 * max(eigenvalues)) {
    stop(errorCondition("The observed information is not positive definite.",
                        class = "latents_singular_information", call = NULL))
  }
  inverse <- chol2inv(chol(information))
  meat <- crossprod(scores)
  if (vcov_type != "observed" && qr(scores)$rank < length(theta)) {
    stop(errorCondition("Too few groups for a robust or OPG covariance.",
                        class = "latents_too_few_groups", call = NULL))
  }
  covariance <- switch(vcov_type, observed = inverse,
                       robust = inverse %*% meat %*% inverse,
                       opg = chol2inv(chol(meat)))
  dimnames(covariance) <- dimnames(information)
  list(theta = theta, vcov = covariance, information = information,
       group_scores = scores, gradient = colSums(scores), vcov_type = vcov_type)
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
           if (!is.null(block$response_probabilities))
             list(response_probabilities = lapply(block$response_probabilities, unname)))
       }),
       initial = unname(fit$initial_coefficients),
       transition = unname(fit$transition_coefficients),
       transition2 = if (is.null(fit$second_order_coefficients)) NULL else
         unname(fit$second_order_coefficients),
       group_probabilities = unname(fit$group_probabilities),
       stayer = unname(fit$stayer %||% rep(FALSE, length(fit$group_probabilities))))
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
  theta <- .lta_pack(parameters, spec$variance_model)
  names_all <- .lta_names(spec)
  lower <- ifelse(startsWith(names_all, "log_variance."), log(min_variance), -Inf)
  theta <- pmax(theta, lower)
  value <- function(v) {
    -.lta_expectation(spec$x, spec$codes, spec$layout, spec$designs,
                      .lta_unpack(v, spec), spec$occasion_of_row,
                      spec$sampling_weights)$log_likelihood
  }
  gradient <- function(v) -colSums(.lta_group_scores(v, spec))
  start_value <- value(theta)
  # Convergence is the relative change of the log likelihood, as for EM: a
  # transition never observed has its maximum at probability zero, where the
  # logit runs off and the gradient decays without vanishing. That boundary is
  # flagged separately (see .lta_probability_boundary()).
  found <- stats::optim(theta, value, gradient, method = "L-BFGS-B", lower = lower,
                        control = list(maxit = 5000L,
                                       factr = max(tol / .Machine$double.eps, 1),
                                       pgtol = 0))
  if (!is.finite(found$value) || found$value > start_value) {
    return(list(parameters = parameters, log_likelihood = -start_value,
                gradient = max(abs(gradient(theta))), converged = FALSE))
  }
  final_gradient <- gradient(found$par)
  free <- !(found$par <= lower + 1e-12 & final_gradient > 0)
  list(parameters = .lta_unpack(found$par, spec), log_likelihood = -found$value,
       gradient = max(abs(final_gradient[free]), 0),
       converged = found$convergence == 0L)
}

#' Does any fitted initial or transition probability sit at zero?
#'
#' Read over every group and every occasion it is observed at; a probability
#' below 1e-6 marks a move (or start) that never occurs, whose logit has no
#' finite maximum.
#' @noRd
.lta_probability_boundary <- function(expectation, layout, stayer = NULL) {
  # A stayer class's zero moves are structural, not estimated: skip them.
  stayer <- stayer %||% rep(FALSE, length(expectation$per_class))
  any(vapply(seq_along(expectation$per_class), function(h) {
    class_terms <- expectation$per_class[[h]]
    if (isTRUE(stayer[h])) {
      return(min(exp(class_terms$log_initial)) < 1e-6)
    }
    initial <- min(exp(class_terms$log_initial))
    moves <- vapply(seq_along(class_terms$log_transition), function(t) {
      active <- layout$within[, t + 1L]
      if (!any(active)) return(1)
      min(exp(class_terms$log_transition[[t]][active, , , drop = FALSE]))
    }, numeric(1))
    second <- vapply(seq_along(class_terms$log_transition2), function(t) {
      active <- layout$within[, t + 2L]
      if (!any(active)) return(1)
      min(exp(class_terms$log_transition2[[t]][active, , , drop = FALSE]))
    }, numeric(1))
    min(initial, moves, second) < 1e-6
  }, logical(1)))
}
