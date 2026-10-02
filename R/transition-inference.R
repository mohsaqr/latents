# Wald inference for latent transition models.
#
# A transition model is the cross-sectional measurement model with a different
# dependence structure over a group's occasions: group-class shares, an initial
# profile distribution per group class and a transition matrix per group
# class. By the Fisher identity the observed-data score is the expected
# complete-data score, so every block's score is an expected count from the
# forward-backward pass minus the count its fitted probabilities imply. The
# measurement blocks score exactly as in multilpa(), weighted by the sequence
# posteriors, and are computed by the same helpers through a measurement view
# of the fit. The observed information is the numerical Jacobian of the
# analytic score, as for every other fit in the package.

#' A `multilpa`-shaped view of a transition fit's measurement model
#'
#' The coefficient encoders, the delta-method Jacobian and the full-covariance
#' group scores of the cross-sectional model read only measurement fields, the
#' group-class shares and the grouping. This view carries exactly those, so
#' the transition path uses the same code rather than a copy of it. Its
#' `profile_probabilities` are the fitted initial distribution, present only
#' because the cross-sectional encoder lays out a block of that shape; the
#' transition path discards that block.
#'
#' @param x A `multilpa_transitions` fit.
#' @param centers Column centres to take off the means, or `NULL`.
#' @return A list of class `multilpa`.
#' @noRd
.multilpa_transition_view <- function(x, centers = NULL) {
  stopifnot(inherits(x, "multilpa_transitions"))
  means <- x$means
  if (!is.null(centers) && ncol(means) > 0L) means <- sweep(means, 2L, centers, "-")
  structure(list(
    means = means, variances = x$variances, covariances = x$covariances,
    response_probabilities = x$response_probabilities,
    profile_probabilities = x$initial_probabilities,
    group_probabilities = x$group_probabilities,
    n_profiles = x$n_profiles, n_group_classes = x$n_group_classes,
    n_groups = x$n_groups, group_index = x$group_index, group_ids = x$group_ids,
    vars = x$vars, continuous = x$continuous, categorical = x$categorical,
    variance_model = x$variance_model, covariance_model = x$covariance_model,
    min_probability = x$min_probability), class = "multilpa")
}

#' Which cross-sectional coordinates the transition model shares
#'
#' Everything but the profile block, which the transition model replaces with
#' its initial and transition logits.
#' @param view A view from `.multilpa_transition_view()`.
#' @param scale `"natural"` or `"unconstrained"`.
#' @return A list of integer positions: `measurement`, `profile`, `group`.
#' @noRd
.multilpa_transition_shared_blocks <- function(view, scale) {
  widths <- .multilpa_coordinate_widths(view, scale)
  n_measurement <- sum(widths[c("means", "variances", "response_probabilities")])
  list(measurement = seq_len(n_measurement),
       profile = n_measurement + seq_len(widths[["profile"]]),
       group = n_measurement + widths[["profile"]] + seq_len(widths[["group"]]))
}

#' Encode a transition fit as one parameter vector
#'
#' Measurement, then group-class shares, then each group class's initial
#' distribution, then each group class's transition matrix row by row (from
#' profile, then to profile). On the unconstrained scale every distribution is
#' a baseline-category logit against its last profile.
#'
#' @param x A `multilpa_transitions` fit.
#' @param scale `"natural"` or `"unconstrained"`.
#' @param centers Column centres taken off the means, or `NULL`.
#' @return A named numeric vector; names in the internal `kind[index,...]`
#'   grammar.
#' @noRd
.multilpa_transition_encode <- function(x, scale, centers = NULL) {
  view <- .multilpa_transition_view(x, centers)
  raw <- .multilpa_coefficients_raw(view, scale)
  blocks <- .multilpa_transition_shared_blocks(view, scale)
  n_profiles <- x$n_profiles
  n_types <- x$n_group_classes
  natural <- identical(scale, "natural")
  kept <- if (natural) seq_len(n_profiles) else seq_len(n_profiles - 1L)
  contrast <- function(probabilities) {
    if (natural) probabilities[kept] else
      log(probabilities[kept]) - log(probabilities[n_profiles])
  }
  initial <- unlist(lapply(seq_len(n_types), function(type) {
    stats::setNames(contrast(x$initial_probabilities[type, ]),
                    sprintf("%s[%d,%d]", if (natural) "initial_probability" else
                      "initial_logit", type, kept))
  }))
  transition <- unlist(lapply(seq_len(n_types), function(type) {
    unlist(lapply(seq_len(n_profiles), function(from) {
      stats::setNames(
        contrast(x$transition_probabilities[from, , type]),
        sprintf("%s[%d,%d,%d]", if (natural) "transition_probability" else
          "transition_logit", type, from, kept))
    }))
  }))
  c(raw[c(blocks$measurement, blocks$group)], initial, transition)
}

#' Decode a transition parameter vector into the parameters the E-step reads
#' @param theta Unconstrained vector in `.multilpa_transition_encode()` order.
#' @param view The measurement view the vector was encoded from.
#' @param n_measurement,n_group Widths of the measurement and group blocks.
#' @return A parameter list for `.multilpa_transition_expectation()`.
#' @noRd
.multilpa_transition_decode <- function(theta, view, n_measurement, n_group) {
  n_profiles <- view$n_profiles
  n_types <- view$n_group_classes
  n_free <- n_profiles - 1L
  measurement <- theta[seq_len(n_measurement)]
  group <- theta[n_measurement + seq_len(n_group)]
  initial <- theta[n_measurement + n_group + seq_len(n_types * n_free)]
  transition <- theta[n_measurement + n_group + n_types * n_free +
                        seq_len(n_types * n_profiles * n_free)]
  ## The cross-sectional decoder rebuilds the measurement and the group shares;
  ## the profile block it expects is filled with zeros and discarded.
  parameters <- .multilpa_decode(c(measurement, numeric(n_types * n_free), group),
                                 view)
  parameters$profile_probabilities <- NULL
  softmax <- function(logits) {
    full <- c(logits, 0)
    exp(full - .multilpa_log_sum_exp(matrix(full, 1L)))
  }
  parameters$initial_probabilities <- t(matrix(vapply(seq_len(n_types), function(type) {
    softmax(initial[(type - 1L) * n_free + seq_len(n_free)])
  }, numeric(n_profiles)), n_profiles, n_types))
  parameters$transition_probabilities <- array(vapply(seq_len(n_types), function(type) {
    rows <- vapply(seq_len(n_profiles), function(from) {
      at <- ((type - 1L) * n_profiles + from - 1L) * n_free
      softmax(transition[at + seq_len(n_free)])
    }, numeric(n_profiles))
    t(rows)
  }, numeric(n_profiles * n_profiles)), c(n_profiles, n_profiles, n_types))
  parameters
}

#' Per-group scores of a transition model's observed-data log likelihood
#'
#' @param theta Unconstrained vector in `.multilpa_transition_encode()` order.
#' @param x Centred continuous indicator matrix, `NA` where unobserved.
#' @param layout The sequence layout.
#' @param view The measurement view.
#' @param widths Named widths `measurement` and `group`.
#' @param codes Categorical codes, or `NULL`.
#' @return A groups-by-parameters matrix whose column sums are the gradient of
#'   the log likelihood.
#' @noRd
.multilpa_transition_group_scores <- function(theta, x, layout, view, widths,
                                              codes = NULL) {
  parameters <- .multilpa_transition_decode(theta, view, widths[["measurement"]],
                                            widths[["group"]])
  expectation <- .multilpa_transition_expectation(x, layout, parameters, codes)
  n_profiles <- view$n_profiles
  n_types <- view$n_group_classes
  group_index <- view$group_index
  ## Measurement: the cross-sectional scores, weighted by the sequence
  ## posteriors of each observation.
  measurement <- if (identical(view$covariance_model, "full")) {
    .multilpa_full_measurement_group_score(parameters, expectation, view)
  } else {
    observed <- !is.na(x)
    blocks <- lapply(seq_len(n_profiles), function(profile) {
      residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
      residuals[!observed] <- 0
      weights <- expectation$subject_posteriors[, profile]
      list(means = rowsum(sweep(residuals, 2L, parameters$variances[profile, ], "/") *
                            weights, group_index, reorder = FALSE),
           variances = rowsum(0.5 * (sweep(residuals^2, 2L,
             parameters$variances[profile, ], "/") - observed) * weights,
             group_index, reorder = FALSE))
    })
    list(means = do.call(cbind, lapply(blocks, `[[`, "means")),
         covariances = if (identical(view$variance_model, "equal")) {
           Reduce(`+`, lapply(blocks, `[[`, "variances"))
         } else do.call(cbind, lapply(blocks, `[[`, "variances")))
  }
  response <- .multilpa_response_scores(codes, expectation$subject_posteriors,
                                        parameters$response_probabilities,
                                        group_index)
  group <- sweep(expectation$group_posteriors, 2L, parameters$group_probabilities,
                 "-")[, seq_len(n_types - 1L), drop = FALSE]
  counts <- .multilpa_transition_group_counts(x, layout, parameters, codes,
                                              expectation$group_posteriors)
  free <- seq_len(n_profiles - 1L)
  initial <- do.call(cbind, lapply(seq_len(n_types), function(type) {
    first <- counts$initial[[type]]
    (first - outer(rowSums(first), parameters$initial_probabilities[type, ]))[
      , free, drop = FALSE]
  }))
  transition <- do.call(cbind, lapply(seq_len(n_types), function(type) {
    pairs <- counts$transition[[type]]
    do.call(cbind, lapply(seq_len(n_profiles), function(from) {
      ## Columns of `pairs` are laid out as the transition matrix is,
      ## column-major: from profile `from` to profile `to` is column
      ## (to - 1) * K + from.
      moves <- pairs[, (seq_len(n_profiles) - 1L) * n_profiles + from, drop = FALSE]
      (moves - outer(rowSums(moves),
                     parameters$transition_probabilities[from, , type]))[
        , free, drop = FALSE]
    }))
  }))
  cbind(measurement$means, measurement$covariances, response, group, initial,
        transition)
}

#' Expected first-occasion and move counts, one row per group
#'
#' The per-group terms that `.multilpa_sequence_moments()` sums over groups,
#' kept apart because the sandwich estimators need each group's own score.
#'
#' @param x Centred continuous indicator matrix.
#' @param layout The sequence layout.
#' @param parameters Decoded parameters.
#' @param codes Categorical codes, or `NULL`.
#' @param group_posteriors Groups-by-group-classes posteriors.
#' @return A list with `initial` (per group class, groups-by-profiles) and
#'   `transition` (per group class, groups-by-profile-pairs, column-major as
#'   the transition matrix).
#' @noRd
.multilpa_transition_group_counts <- function(x, layout, parameters, codes,
                                              group_posteriors) {
  n_profiles <- nrow(parameters$means)
  n_groups <- nrow(layout$slot)
  log_density <- if (ncol(x) > 0L) {
    .multilpa_gaussian_moments(x, parameters)$log_density
  } else matrix(0, nrow(x), n_profiles)
  if (!is.null(codes)) {
    log_density <- log_density +
      .multilpa_categorical_log_density(codes, parameters$response_probabilities)
  }
  emission <- .multilpa_sequence_emission(log_density, layout)$log_density
  origin <- rep(seq_len(n_profiles), times = n_profiles)
  destination <- rep(seq_len(n_profiles), each = n_profiles)
  per_type <- lapply(seq_len(ncol(group_posteriors)), function(type) {
    transition <- matrix(parameters$transition_probabilities[, , type],
                         n_profiles, n_profiles)
    pass <- .multilpa_forward_backward(emission, layout,
                                       parameters$initial_probabilities[type, ],
                                       transition)
    weights <- group_posteriors[, type]
    initial <- exp(pass$alpha[[1L]] + pass$beta[[1L]] - pass$log_scaled) * weights
    log_transition <- rep(log(transition), each = n_groups)
    moves <- if (layout$n_occasions < 2L) {
      matrix(0, n_groups, n_profiles * n_profiles)
    } else Reduce(`+`, lapply(seq_len(layout$n_occasions)[-1L], function(occasion) {
      active <- weights * layout$within[, occasion]
      log_from <- pass$alpha[[occasion - 1L]] - pass$log_scaled
      log_to <- pass$beta[[occasion]] + emission[[occasion]]
      joint <- log_from[, origin, drop = FALSE] + log_to[, destination, drop = FALSE] +
        log_transition
      active * exp(joint)
    }))
    list(initial = initial, transition = moves)
  })
  list(initial = lapply(per_type, `[[`, "initial"),
       transition = lapply(per_type, `[[`, "transition"))
}

#' Delta-method Jacobian from a transition fit's logits to its probabilities
#' @param x A `multilpa_transitions` fit.
#' @param view Its measurement view.
#' @return A natural-by-unconstrained matrix.
#' @noRd
.multilpa_transition_jacobian <- function(x, view) {
  natural_blocks <- .multilpa_transition_shared_blocks(view, "natural")
  free_blocks <- .multilpa_transition_shared_blocks(view, "unconstrained")
  shared <- .multilpa_inference_jacobian(view)[
    c(natural_blocks$measurement, natural_blocks$group),
    c(free_blocks$measurement, free_blocks$group), drop = FALSE]
  n_profiles <- x$n_profiles
  n_types <- x$n_group_classes
  ## Each distribution is its own simplex, with the derivative diag(p) - p p'
  ## restricted to the free (non-reference) logits.
  simplex <- function(probabilities) {
    derivative <- diag(probabilities, nrow = length(probabilities)) -
      tcrossprod(probabilities)
    derivative[, seq_len(n_profiles - 1L), drop = FALSE]
  }
  distributions <- c(
    lapply(seq_len(n_types), function(type) x$initial_probabilities[type, ]),
    unlist(lapply(seq_len(n_types), function(type) {
      lapply(seq_len(n_profiles), function(from) {
        x$transition_probabilities[from, , type]
      })
    }), recursive = FALSE))
  blocks <- lapply(distributions, simplex)
  sequence_block <- matrix(0, n_profiles * length(blocks),
                           (n_profiles - 1L) * length(blocks))
  invisible(lapply(seq_along(blocks), function(index) {
    rows <- (index - 1L) * n_profiles + seq_len(n_profiles)
    columns <- (index - 1L) * (n_profiles - 1L) + seq_len(n_profiles - 1L)
    sequence_block[rows, columns] <<- blocks[[index]]
  }))
  top <- cbind(shared, matrix(0, nrow(shared), ncol(sequence_block)))
  bottom <- cbind(matrix(0, nrow(sequence_block), ncol(shared)), sequence_block)
  rbind(top, bottom)
}

#' Tidy labels for a transition fit's coefficients
#' @param names Internal names from `.multilpa_transition_encode()`.
#' @param view The measurement view.
#' @return A data frame with `level`, `outcome`, `term` and `parameter`.
#' @noRd
.multilpa_transition_labels <- function(names, view) {
  prefix <- sub("\\[.*$", "", names)
  sequence <- prefix %in% c("initial_probability", "initial_logit",
                            "transition_probability", "transition_logit")
  labels <- data.frame(level = character(length(names)),
                       outcome = character(length(names)),
                       term = character(length(names)),
                       parameter = character(length(names)),
                       stringsAsFactors = FALSE)
  if (any(!sequence)) {
    labels[!sequence, ] <- .multilpa_tidy_labels(names[!sequence], view)[
      c("level", "outcome", "term", "parameter")]
  }
  if (any(sequence)) {
    inside <- strsplit(sub("^[^\\[]*\\[", "", sub("\\]$", "", names[sequence])),
                       ",", fixed = TRUE)
    type <- vapply(inside, `[[`, character(1), 1L)
    initial <- prefix[sequence] %in% c("initial_probability", "initial_logit")
    to <- vapply(inside, function(part) part[[length(part)]], character(1))
    from <- vapply(inside, function(part) if (length(part) == 3L) part[[2L]] else
      NA_character_, character(1))
    ## The decomposition coef() already uses for these parameters: the profile
    ## moved to (or started in) is the outcome, the group class and, for a
    ## move, the profile moved from are the term.
    labels$level[sequence] <- "profile"
    labels$outcome[sequence] <- paste0("profile_", to)
    labels$term[sequence] <- ifelse(initial, paste0("group_class_", type),
                                    paste0("group_class_", type, ":profile_", from))
    labels$parameter[sequence] <- prefix[sequence]
  }
  labels
}

#' Wald inference for a latent transition model
#'
#' Standard errors, tests and intervals for every parameter of an [lta()]
#' fit: the measurement model, the group-class shares, each group class's
#' initial profile distribution and its transition matrix. The scores follow
#' the Fisher identity (each is an expected count from the forward-backward
#' pass minus the count the fitted probabilities imply), and the observed
#' information is the numerical Jacobian of that analytic score.
#'
#' @param x A fitted `multilpa_transitions` model from [lta()].
#' @param data Optional. The data frame the model was fitted to. The fit stores
#'   its indicators, codes, groups and occasions, so it is not needed; when
#'   given it must reproduce them.
#' @param level Confidence level of the intervals.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` (observed information), `"robust"` (the
#'   sandwich clustered on groups, which are the sequences) or `"opg"` (the
#'   outer product of the group scores).
#' @param adjust Multiplicity correction for `p_adjusted`, one of the methods
#'   [stats::p.adjust()] accepts.
#' @param method `"wald"` (the default) or `"bootstrap"`. The bootstrap
#'   resamples persons with replacement, refits the same specification, and
#'   matches each replicate's profiles and group classes to the original's;
#'   `standard_error` is the replicates' standard deviation, `conf_low` and
#'   `conf_high` their percentile interval, and `statistic`, `p_value` and
#'   `p_adjusted` are `NA`. It needs `data`.
#' @param iter,n_starts,max_iter,tol,seed Bootstrap resamples, the controls of
#'   each refit, and an optional seed (the caller's random state is restored).
#' @param boundary Accepted for the generic. Transition inference holds
#'   nothing at a bound: with no bound active `"fix"` gives the same table as
#'   `"error"`, and a fit with an active bound raises `latents_boundary_fit`
#'   under either.
#' @param ... Unused.
#' @return A base `data.frame` with one row per natural parameter and the
#'   columns of [parameter_inference()]: `level`, `outcome`, `term`,
#'   `parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
#'   `p_adjusted`, `conf_low` and `conf_high`, named as [coef()] names them.
#'   An initial probability (`parameter = "initial_probability"`, `level =
#'   "profile"`) has `outcome` the profile and `term` the group class; a
#'   transition probability (`"transition_probability"`) has `outcome` the
#'   profile moved to and `term` `group_class_h:profile_j`, the group class
#'   and the profile moved from. Probabilities and variances carry no Wald
#'   test. Attributes `covariance` (natural scale),
#'   `covariance_unconstrained`, `vcov_type` and `level`.
#' @section Conditions:
#'   `latents_no_converge` for an unconverged fit; `latents_boundary_fit` when
#'   a variance, a response, initial or transition probability, or a
#'   group-class share sits at its bound, or a transition row was never
#'   informed (a profile no group occupies before its last occasion);
#'   `latents_singular_information` when the information cannot be inverted;
#'   `latents_bad_inference_data` when a supplied `data` does not reproduce
#'   the fit; `latents_too_few_groups` for robust or OPG errors with no more
#'   groups than parameters; `latents_bad_argument` for `method = "bootstrap"`
#'   without `data`.
#' @examples
#' fit <- lta(subset(course_engagement, student <= 40),
#'            c("browse", "lectures", "forum_read"), "student",
#'            n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
#' parameter_inference(fit)
#' @rdname parameter_inference.multilpa_transitions
#' @export
parameter_inference.multilpa_transitions <- function(x, data = NULL, level = 0.95,
                                                     step = 1e-4,
                                                     vcov_type = c("observed", "robust", "opg"),
                                                     adjust = .multilpa_p_adjust_methods,
                                                     method = c("wald", "bootstrap"),
                                                     iter = 199L, n_starts = 10L,
                                                     max_iter = 1000L, tol = 1e-8,
                                                     seed = NULL,
                                                     boundary = c("error", "fix")) {
  stopifnot(inherits(x, "multilpa_transitions"),
            "`level` must be a single number in (0, 1)" =
              is.numeric(level) && length(level) == 1L && is.finite(level) &&
              level > 0 && level < 1,
            "`step` must be a single positive number" =
              is.numeric(step) && length(step) == 1L && is.finite(step) && step > 0)
  # Nothing is held at a bound here and a bound-active fit is refused, so
  # "fix" and "error" agree wherever a table is returned.
  match.arg(boundary)
  adjust <- match.arg(adjust)
  if (identical(match.arg(method), "bootstrap")) {
    return(.lta_bootstrap_dispatch(x, data, level, iter, n_starts, max_iter, tol,
                                   seed, adjust))
  }
  vcov_type <- match.arg(vcov_type)
  inference <- .multilpa_transition_covariance(x, data, step, vcov_type)
  estimates <- .multilpa_transition_encode(x, "natural")
  covariance <- inference$covariance
  standard_errors <- sqrt(pmax(diag(covariance), 0))
  critical <- stats::qnorm((1 + level) / 2)
  statistic <- unname(estimates / standard_errors)
  labels <- .multilpa_transition_labels(names(estimates), inference$view)
  result <- data.frame(
    labels, estimate = unname(estimates), standard_error = unname(standard_errors),
    statistic = statistic, p_value = 2 * stats::pnorm(-abs(statistic)),
    conf_low = unname(estimates - critical * standard_errors),
    conf_high = unname(estimates + critical * standard_errors),
    row.names = NULL, stringsAsFactors = FALSE)
  continuous <- .multilpa_continuous_names(x)
  bounds <- .multilpa_wald_bounds(
    result$estimate, result$standard_error,
    .multilpa_interval_kind(result$parameter, result$term, continuous), critical)
  result$conf_low <- bounds$low
  result$conf_high <- bounds$high
  bounded <- result$parameter %in% c("variance", "probability", "response",
                                     "initial_probability",
                                     "transition_probability") |
    (result$parameter == "covariance" &
       result$term %in% paste(continuous, continuous, sep = ":"))
  result$statistic[bounded] <- NA_real_
  result$p_value[bounded] <- NA_real_
  result <- .multilpa_adjust_p(result, adjust)
  names_natural <- .multilpa_parameter_names(labels)
  dimnames(covariance) <- list(names_natural, names_natural)
  attributes(result) <- c(attributes(result), list(
    covariance = covariance,
    covariance_unconstrained = inference$covariance_unconstrained,
    level = level, step = step, method = "wald", vcov_type = vcov_type))
  result
}

#' Covariance of a transition fit's parameters
#' @param x A `multilpa_transitions` fit.
#' @param data Optional fitting data.
#' @param step Finite-difference step.
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`.
#' @return A list with `covariance` (natural), `covariance_unconstrained` and
#'   the measurement `view`.
#' @noRd
.multilpa_transition_covariance <- function(x, data, step, vcov_type) {
  .multilpa_check_transition_regularity(x, vcov_type)
  x_matrix <- x$indicator_data
  if (is.null(x_matrix)) x_matrix <- matrix(numeric(0), x$n_observations, 0L)
  x_matrix <- matrix(as.numeric(x_matrix), nrow(x_matrix), ncol(x_matrix))
  if (!is.null(data)) .multilpa_check_transition_data(x, data)
  centers <- if (ncol(x_matrix) > 0L) colMeans(x_matrix, na.rm = TRUE) else numeric(0)
  centred <- if (ncol(x_matrix) > 0L) sweep(x_matrix, 2L, centers, "-") else x_matrix
  codes <- x$categorical_data
  layout <- .multilpa_sequence_layout(x$group_index, x$time_values, x$n_groups,
                                      x$occasions)
  view <- .multilpa_transition_view(x, centers)
  blocks <- .multilpa_transition_shared_blocks(view, "unconstrained")
  widths <- c(measurement = length(blocks$measurement), group = length(blocks$group))
  theta <- .multilpa_transition_encode(x, "unconstrained", centers)
  stopifnot("the coordinates must match the parameters the fit counts" =
              length(theta) == x$n_parameters)
  negative <- function(value) {
    -.multilpa_transition_expectation(
      centred, layout,
      .multilpa_transition_decode(value, view, widths[["measurement"]],
                                  widths[["group"]]), codes)$log_likelihood
  }
  if (abs(-negative(unname(theta)) - x$log_likelihood) >
      1e-8 * (1 + abs(x$log_likelihood))) {
    stop(errorCondition(
      "The stored data do not reproduce the fitted log likelihood.",
      class = "latents_bad_inference_data", call = NULL))
  }
  gradient <- function(value) {
    -colSums(.multilpa_transition_group_scores(value, centred, layout, view,
                                               widths, codes))
  }
  ## Each coordinate on its own unit, as in the cross-sectional path: means in
  ## their indicator's within-profile standard deviation, logits as they are.
  scale <- rep(1, length(theta))
  n_means <- length(x$means)
  if (n_means > 0L) scale[seq_len(n_means)] <- as.vector(t(sqrt(x$variances)))
  dimension <- length(.multilpa_continuous_names(x))
  if (identical(x$covariance_model, "full") && dimension > 0L) {
    covariance_scale <- .multilpa_covariance_coordinate_scale(
      x$variances, dimension,
      if (identical(x$variance_model, "equal")) 1L else seq_len(x$n_profiles))
    scale[n_means + seq_along(covariance_scale)] <- covariance_scale
  }
  start <- unname(theta)
  information <- .multilpa_observed_hessian(
    function(displacement) negative(start + displacement * scale),
    function(displacement) gradient(start + displacement * scale) * scale,
    scale, step)
  scaled_inverse <- information$inverse
  if (vcov_type %in% c("robust", "opg")) {
    scores <- sweep(.multilpa_transition_group_scores(start, centred, layout, view,
                                                      widths, codes),
                    2L, scale, "*")
    cross <- .multilpa_cross_product(scores)
    scaled_inverse <- if (identical(vcov_type, "opg")) {
      .multilpa_opg_inverse(cross)
    } else scaled_inverse %*% cross %*% scaled_inverse
  }
  covariance_unconstrained <- scaled_inverse * tcrossprod(scale)
  dimnames(covariance_unconstrained) <- list(names(theta), names(theta))
  jacobian <- .multilpa_transition_jacobian(x, view)
  list(covariance = jacobian %*% covariance_unconstrained %*% t(jacobian),
       covariance_unconstrained = covariance_unconstrained, view = view)
}

#' Refuse Wald inference where it is not defined for a transition fit
#' @param x A `multilpa_transitions` fit.
#' @param vcov_type The requested covariance.
#' @return `NULL`, invisibly; raises on the first broken condition.
#' @noRd
.multilpa_check_transition_regularity <- function(x, vcov_type) {
  refuse <- function(message, class) {
    stop(errorCondition(message, class = class, call = NULL))
  }
  if (!isTRUE(x$converged)) {
    refuse("Inference requires a converged fit.", "latents_no_converge")
  }
  if (!is.null(x$covariance_structure) &&
      !x$covariance_structure %in% c("EEI", "VVI", "EEE", "VVV")) {
    refuse(paste(
      "Standard errors for transition models are available for the EEI, VVI,",
      "EEE and VVV structures; this fit uses", x$covariance_structure, "."),
      "latents_unsupported_inference")
  }
  if (vcov_type %in% c("robust", "opg") && x$n_groups <= x$n_parameters) {
    refuse(sprintf(paste(
      "Robust and OPG inference need more groups than parameters; this fit has",
      "%d groups and %d parameters."), x$n_groups, x$n_parameters),
      "latents_too_few_groups")
  }
  if (isTRUE(x$boundary)) {
    refuse("Wald inference is unavailable for a bound-active fit.",
           "latents_boundary_fit")
  }
  if (any(x$empty_transition_rows)) {
    refuse(paste(
      "A transition row was never informed: some profile is occupied by no",
      "group before its last occasion, so its outgoing probabilities are not",
      "estimated and have no standard error."), "latents_boundary_fit")
  }
  floor <- .multilpa_probability_floor(x)
  probabilities <- c(unlist(x$response_probabilities, use.names = FALSE),
                     x$initial_probabilities, x$transition_probabilities,
                     x$group_probabilities)
  if (any(probabilities <= floor)) {
    refuse(paste(
      "Wald inference is unavailable when a probability sits at its lower",
      "bound: its logit is at minus infinity, where the likelihood has no",
      "curvature to invert."), "latents_boundary_fit")
  }
  invisible(NULL)
}

#' Check that supplied data reproduce a transition fit's stored data
#' @param x A `multilpa_transitions` fit.
#' @param data The supplied data frame.
#' @return `NULL`, invisibly; raises `latents_bad_inference_data`.
#' @noRd
.multilpa_check_transition_data <- function(x, data) {
  continuous <- .multilpa_continuous_names(x)
  same <- function(left, right) identical(unname(left), unname(right))
  ok <- is.data.frame(data) && nrow(data) == x$n_observations &&
    all(c(x$vars, x$id, x$time) %in% names(data)) &&
    identical(match(data[[x$id]], x$group_values), x$group_index) &&
    (length(continuous) == 0L ||
       same(as.matrix(data[continuous]), x$indicator_data))
  if (ok && length(x$categorical %||% character()) > 0L) {
    encoded <- .multilpa_encode_categorical(data[x$categorical])
    ok <- same(encoded$codes, x$categorical_data) &&
      identical(encoded$levels, x$categorical_levels)
  }
  if (ok) ok <- identical(.multilpa_time_values(data, x$time, x$id, x$vars), x$time_values)
  if (!ok) {
    stop(errorCondition(
      "`data` must reproduce the indicators, categories, groups, occasions and row order of the fit.",
      class = "latents_bad_inference_data", call = NULL))
  }
  invisible(NULL)
}

#' @rdname parameter_inference.multilpa_transitions
#' @param object A fitted `multilpa_transitions` model.
#' @param scale `"natural"` (probabilities) or `"unconstrained"` (logits).
#' @export
vcov.multilpa_transitions <- function(object, data = NULL, step = 1e-4,
                                      vcov_type = c("observed", "robust", "opg"),
                                      scale = c("natural", "unconstrained"), ...) {
  vcov_type <- match.arg(vcov_type)
  scale <- match.arg(scale)
  inference <- .multilpa_transition_covariance(object, data, step, vcov_type)
  if (identical(scale, "unconstrained")) {
    covariance <- inference$covariance_unconstrained
    names_used <- .multilpa_parameter_names(.multilpa_transition_labels(
      rownames(covariance), inference$view))
  } else {
    covariance <- inference$covariance
    names_used <- .multilpa_parameter_names(.multilpa_transition_labels(
      names(.multilpa_transition_encode(object, "natural")), inference$view))
  }
  dimnames(covariance) <- list(names_used, names_used)
  covariance
}

#' @rdname parameter_inference.multilpa_transitions
#' @param parm Parameters to report, as names or positions of the natural
#'   coefficients; all by default.
#' @export
confint.multilpa_transitions <- function(object, parm, level = 0.95, data = NULL,
                                         ...) {
  inference <- parameter_inference(object, data = data, level = level, ...)
  names_used <- .multilpa_parameter_names(inference)
  bounds <- cbind(inference$conf_low, inference$conf_high)
  dimnames(bounds) <- list(names_used, paste0(
    format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), " %"))
  if (missing(parm)) return(bounds)
  bounds[parm, , drop = FALSE]
}
