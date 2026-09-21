#' Pack a covariate fit's free parameters into one vector
#'
#' Variances are carried on the log scale so the optimizer and the numerical
#' Hessian both work on an unbounded space; the delta method returns them.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param means,variances,beta,gamma Optional replacements, for decoding.
#' @return A numeric vector: means, log variances, profile logits, group logits.
#' @noRd
.multilpa_cov_encode <- function(object) {
  variances <- if (identical(object$variance_model, "equal")) {
    object$variances[1L, , drop = FALSE]
  } else object$variances
  ## The fit reports means in input units; the likelihood is evaluated on
  ## centred indicators, so the centre comes off here and goes back on in the
  ## reported estimates. A location shift leaves the standard errors alone.
  centred_means <- sweep(object$means, 2L, object$center, "-")
  spread <- if (.multilpa_cov_is_full(object)) {
    .multilpa_cov_cholesky_coordinates(object)
  } else as.vector(t(log(variances)))
  c(as.vector(t(centred_means)), spread,
    as.vector(object$profile_coefficients), as.vector(object$group_coefficients))
}

#' Does this covariate fit carry unrestricted residual covariances?
#' @param object A fitted `multilpa_covariates` model.
#' @return A single logical.
#' @noRd
.multilpa_cov_is_full <- function(object) {
  identical(object$covariance_model %||% "diagonal", "full")
}

#' Which rows of the spread block are estimated separately
#' @param object A fitted `multilpa_covariates` model.
#' @return Integer profile indices, or `1L` when the spread is shared.
#' @noRd
.multilpa_cov_spread_rows <- function(object) {
  if (identical(object$variance_model, "equal")) 1L else seq_len(object$n_profiles)
}

#' Residual covariances as log-Cholesky coordinates
#'
#' The same parameterization the covariate-free model uses: the lower triangle
#' of the Cholesky factor, with its diagonal logged so the parameter is
#' unconstrained and the covariance it implies is positive definite by
#' construction.
#'
#' @param object A fitted full-covariance `multilpa_covariates` model.
#' @return A numeric vector, one block per estimated profile.
#' @noRd
.multilpa_cov_cholesky_coordinates <- function(object) {
  dimension <- ncol(object$means)
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  unlist(lapply(.multilpa_cov_spread_rows(object), function(profile) {
    factor <- t(chol(matrix(object$covariances[, , profile], dimension, dimension)))
    diag(factor) <- log(diag(factor))
    factor[lower]
  }), use.names = FALSE)
}

#' Rebuild a covariance array from log-Cholesky coordinates
#' @param values Numeric vector of coordinates, one block per estimated profile.
#' @param dimension Number of continuous indicators.
#' @param n_profiles Number of profiles the array must carry.
#' @param shared Whether one block is shared across profiles.
#' @return An indicators-by-indicators-by-profiles array.
#' @noRd
.multilpa_cov_cholesky_decode <- function(values, dimension, n_profiles, shared) {
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  width <- sum(lower)
  blocks <- lapply(seq_len(if (shared) 1L else n_profiles), function(index) {
    factor <- matrix(0, dimension, dimension)
    factor[lower] <- values[(index - 1L) * width + seq_len(width)]
    diag(factor) <- exp(diag(factor))
    tcrossprod(factor)
  })
  if (shared) blocks <- rep(blocks, n_profiles)
  array(unlist(blocks, use.names = FALSE), c(dimension, dimension, n_profiles))
}

#' Unpack a parameter vector for a covariate fit
#' @param theta Numeric vector as built by [.multilpa_cov_encode()].
#' @param object The fit supplying the shapes.
#' @return A list with `parameters` (means, variances), `beta` and `gamma`.
#' @noRd
.multilpa_cov_decode <- function(theta, object) {
  n_profiles <- object$n_profiles
  n_indicators <- length(.multilpa_continuous_names(object))
  n_variance_rows <- if (identical(object$variance_model, "equal")) 1L else n_profiles
  at <- 0L
  take <- function(n) {
    values <- theta[at + seq_len(n)]
    at <<- at + n
    values
  }
  means <- matrix(take(n_profiles * n_indicators), n_profiles, n_indicators,
                  byrow = TRUE)
  covariances <- NULL
  if (.multilpa_cov_is_full(object)) {
    width <- n_indicators * (n_indicators + 1L) / 2L
    covariances <- .multilpa_cov_cholesky_decode(
      take(length(.multilpa_cov_spread_rows(object)) * width), n_indicators,
      n_profiles, n_variance_rows == 1L)
    variances <- t(matrix(vapply(seq_len(n_profiles), function(profile) {
      diag(matrix(covariances[, , profile], n_indicators, n_indicators))
    }, numeric(n_indicators)), n_indicators, n_profiles))
  } else {
    variances <- exp(matrix(take(n_variance_rows * n_indicators), n_variance_rows,
                            n_indicators, byrow = TRUE))
    if (n_variance_rows == 1L) {
      variances <- matrix(variances, n_profiles, n_indicators, byrow = TRUE)
    }
  }
  beta <- matrix(take(length(object$profile_coefficients)),
                 nrow(object$profile_coefficients), ncol(object$profile_coefficients))
  gamma <- matrix(take(length(object$group_coefficients)),
                  nrow(object$group_coefficients), ncol(object$group_coefficients))
  parameters <- list(means = means, variances = variances)
  if (!is.null(covariances)) parameters$covariances <- covariances
  list(parameters = parameters, beta = beta, gamma = gamma)
}

#' Per-group score contributions for a covariate fit
#'
#' The observed-data score equals the expected complete-data score at the
#' maximum, so each block is a posterior-weighted residual: observed
#' responsibility minus model-implied probability for the two multinomial
#' logits, and the usual Gaussian residuals for the measurement model. Rows are
#' summed within group because groups are the independent units.
#'
#' @param theta Parameter vector.
#' @param x Centred indicator matrix.
#' @param object The fit, supplying designs and shapes.
#' @return A matrix with one row per group and one column per parameter.
#' @noRd
.multilpa_cov_group_scores <- function(theta, x, object) {
  pieces <- .multilpa_cov_decode(theta, object)
  expectation <- .multilpa_cov_expectation(
    x, object$group_index, pieces$parameters, object$profile_design,
    object$group_design, pieces$beta, pieces$gamma)
  n_groups <- nrow(object$group_design)
  n_profiles <- object$n_profiles
  n_group_classes <- object$n_group_classes
  posteriors <- expectation$subject_posteriors
  equal <- identical(object$variance_model, "equal")

  ## Measurement: means, then the spread on the scale it is estimated on.
  if (.multilpa_cov_is_full(object)) {
    measurement <- .multilpa_cov_full_scores(x, pieces$parameters, posteriors, object)
    mean_block <- measurement$means
    variance_block <- measurement$covariances
  } else {
  mean_block <- do.call(cbind, lapply(seq_len(n_profiles), function(k) {
    residual <- sweep(x, 2L, pieces$parameters$means[k, ], "-")
    weighted <- residual * posteriors[, k]
    rowsum(sweep(weighted, 2L, pieces$parameters$variances[k, ], "/"),
           object$group_index, reorder = FALSE)
  }))
  variance_pieces <- lapply(seq_len(n_profiles), function(k) {
    residual <- sweep(x, 2L, pieces$parameters$means[k, ], "-")
    standardized <- sweep(residual^2, 2L, pieces$parameters$variances[k, ], "/")
    rowsum(0.5 * (standardized - 1) * posteriors[, k], object$group_index,
           reorder = FALSE)
  })
  variance_block <- if (equal) Reduce(`+`, variance_pieces) else
    do.call(cbind, variance_pieces)
  }

  ## Profile logits: the design row depends on the group class, so the
  ## contribution is summed over classes weighted by the group's posterior.
  group_weight <- expectation$group_posteriors[object$group_index, , drop = FALSE]
  beta_block <- do.call(cbind, lapply(seq_len(ncol(pieces$beta)), function(k) {
    contribution <- Reduce(`+`, lapply(seq_len(n_group_classes), function(h) {
      residual <- expectation$joint[[h]][, k] -
        group_weight[, h] * expectation$profile_priors[[h]][, k]
      object$profile_design[[h]] * residual
    }))
    rowsum(contribution, object$group_index, reorder = FALSE)
  }))

  ## Group logits: one row per group already.
  gamma_block <- do.call(cbind, lapply(seq_len(ncol(pieces$gamma)), function(h) {
    object$group_design * (expectation$group_posteriors[, h] -
                             expectation$group_priors[, h])
  }))

  scores <- cbind(mean_block, variance_block, beta_block, gamma_block)
  stopifnot("scores must be one row per group" = nrow(scores) == n_groups,
            "scores must be one column per parameter" = ncol(scores) == length(theta))
  scores
}

#' Per-group measurement scores in log-Cholesky coordinates
#'
#' The same quantities the covariate-free model forms, but kept one row per
#' group rather than summed, because groups are the independent units the
#' sandwich resamples over. Covariate fits are complete-data only, so the
#' conditional-moment corrections the missing-data model needs are all zero and
#' the scatter is formed directly from the residuals.
#'
#' @param x Centred continuous indicator matrix.
#' @param parameters Decoded means and covariances.
#' @param posteriors Observation-by-profile responsibilities.
#' @param object The fit, supplying the group index and shapes.
#' @return A list with `means` and `covariances`, each a groups-by-parameters
#'   matrix in encode order.
#' @noRd
.multilpa_cov_full_scores <- function(x, parameters, posteriors, object) {
  dimension <- ncol(x)
  n_groups <- nrow(object$group_design)
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  index <- object$group_index
  per_profile <- lapply(seq_len(object$n_profiles), function(k) {
    covariance <- matrix(parameters$covariances[, , k], dimension, dimension)
    precision <- chol2inv(chol(covariance))
    weights <- posteriors[, k]
    residual <- sweep(x, 2L, parameters$means[k, ], "-")
    mean_score <- rowsum(residual * weights, index, reorder = FALSE) %*% precision
    # Every entry of the per-group scatter is a weighted sum of a product of
    # two residual columns, so all of them come from one rowsum.
    pairs <- residual[, rep(seq_len(dimension), each = dimension), drop = FALSE] *
      residual[, rep(seq_len(dimension), times = dimension), drop = FALSE]
    scatter <- rowsum(pairs * weights, index, reorder = FALSE)
    totals <- as.vector(rowsum(weights, index, reorder = FALSE))
    list(means = mean_score, scatter = scatter, totals = totals,
         precision = precision, covariance = covariance)
  })
  spread_rows <- .multilpa_cov_spread_rows(object)
  covariance_block <- do.call(cbind, lapply(spread_rows, function(profile) {
    factor <- t(chol(matrix(parameters$covariances[, , profile], dimension, dimension)))
    contributors <- if (identical(object$variance_model, "equal"))
      seq_len(object$n_profiles) else profile
    rows <- lapply(seq_len(n_groups), function(g) {
      gradient <- Reduce(`+`, lapply(contributors, function(k) {
        piece <- per_profile[[k]]
        scatter <- matrix(piece$scatter[g, ], dimension, dimension)
        0.5 * (piece$precision %*% scatter %*% piece$precision -
                 piece$totals[g] * piece$precision)
      }))
      # Chain from the covariance to the factor, then to its logged diagonal.
      factor_score <- 2 * gradient %*% factor
      diag(factor_score) <- diag(factor_score) * diag(factor)
      factor_score[lower]
    })
    matrix(unlist(rows, use.names = FALSE), n_groups, sum(lower), byrow = TRUE)
  }))
  list(means = do.call(cbind, lapply(per_profile, `[[`, "means")),
       covariances = covariance_block)
}

#' Tidy inference for a covariate fit
#'
#' Standard errors, tests and intervals for every free parameter of
#' `multilpa(profile_covariates = )`, at all three levels at once: the
#' measurement model, the
#' profile logits and the group-class logits. Without them a membership
#' coefficient cannot be reported, because nothing distinguishes a real effect
#' from separation.
#'
#' @param x A fitted `multilpa_covariates` model.
#' @param data Optional. The data frame the model was fitted to; when omitted
#'   it is rebuilt from the indicators, group index and designs the fit stores.
#'   Supplying it is the stronger check that the caller still holds that frame.
#' @param level Confidence level for the intervals, between 0 and 1.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` uses the observed information;
#'   `"robust"` uses the group-clustered sandwich, which is the appropriate
#'   choice when the measurement model may be misspecified.
#' @param adjust Multiplicity correction applied to `p_value` to produce
#'   `p_adjusted`, defaulting to `"none"`; see [parameter_inference()].
#' @return A base `data.frame` with one row per free parameter and the same
#'   columns, in the same order, that [parameter_inference()] returns for every
#'   other fitted class: `level` (`"measurement"`, `"profile"` or `"group"`),
#'   `outcome` (which profile or group class the coefficient predicts, or which
#'   profile a measurement parameter belongs to), `term`, `parameter`
#'   (`"mean"`, `"variance"`, `"covariance"` or `"coefficient"`), `estimate`,
#'   `standard_error`, `statistic`, `p_value`, `p_adjusted`, `conf_low` and
#'   `conf_high`. Variances and covariances are reported in their natural units,
#'   with standard errors carried through the delta method from the log and
#'   log-Cholesky coordinates on which they are estimated.
#' @details Class-membership coefficients are on the multinomial logit scale
#'   with the final profile and the final group class as references, so a
#'   coefficient is a log odds against that reference. The tests are Wald tests
#'   and inherit the usual caveat: they are unreliable for a coefficient driven
#'   to the boundary by separation, which the size of the estimate and its
#'   standard error together will reveal.
#' @examples
#' # Two separated profiles, membership driven by `x` and by the group's class.
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(
#'   runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#' )
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                 n_group_classes = 2, profile_covariates = "x",
#'                 n_starts = 2, seed = 1)
#' parameter_inference(fit, example_data)
#' @rdname parameter_inference
#' @export
parameter_inference.multilpa_covariates <- function(x, data = NULL, level = 0.95,
                                                    step = 1e-4,
                                                    vcov_type = c("observed", "robust"),
                                                    adjust = .multilpa_p_adjust_methods) {
  stopifnot(
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1
  )
  vcov_type <- match.arg(vcov_type)
  adjust <- match.arg(adjust)
  covariance <- .multilpa_cov_covariance(x, data, step, vcov_type)
  .multilpa_cov_inference_frame(x, .multilpa_cov_encode(x), covariance,
                                level, vcov_type, adjust)
}

#' Covariance of a covariate fit's free parameters
#' @param object A fitted `multilpa_covariates` model.
#' @param data Optional. The data frame the model was fitted to; when omitted
#'   it is rebuilt from the indicators, group index and designs the fit stores.
#'   Supplying it is the stronger check that the caller still holds that frame.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` or `"robust"`.
#' @return A square matrix, ordered as [.multilpa_cov_encode()].
#' @noRd
.multilpa_cov_covariance <- function(object, data = NULL, step, vcov_type) {
  stopifnot(
    "`object` must be a fitted `multilpa_covariates` model" =
      inherits(object, "multilpa_covariates"))
  # A verb must not demand what the object already owns. The fit stores the
  # indicator matrix, the group index and both designs, which is everything the
  # observed information is rebuilt from; passing `data` remains the stronger
  # check that the caller still holds the frame that produced the estimates.
  if (is.null(data)) data <- .multilpa_cov_stored_data(object)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`step` must be a single positive number" =
      is.numeric(step) && length(step) == 1L && is.finite(step) && step > 0
  )
  if (is.null(object$profile_design) || is.null(object$group_design)) {
    stop(errorCondition(
      "This fit predates covariate inference; refit with the current version.",
      class = "multilpa_unsupported_inference", call = NULL))
  }
  if (!is.null(object$categorical) && length(object$categorical) > 0L) {
    stop(errorCondition(paste(
      "Standard errors are not available for a covariate fit with categorical",
      "indicators. The measurement block has no score implemented for its",
      "response probabilities, and reporting the other blocks alone would",
      "understate the parameter count."),
      class = "multilpa_unsupported_inference", call = NULL))
  }
  .multilpa_check_regularity(object, vcov_type)
  columns <- unique(c(object$vars, object$id,
                      object$profile_covariates, object$group_covariates))
  if (!all(columns %in% names(data)) || anyDuplicated(names(data)) ||
      nrow(data) != object$n_observations ||
      !all(vapply(data[setdiff(columns, object$id)], function(column) {
        is.numeric(column) && is.null(dim(column)) && all(is.finite(column))
      }, logical(1))) ||
      !identical(match(data[[object$id]], object$group_values), object$group_index)) {
    stop(errorCondition(
      "`data` must contain the original finite indicators, covariates and group identifiers in fitting order.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  first_rows <- match(seq_len(object$n_groups), object$group_index)
  designs <- .multilpa_cov_designs(data, object$vars,
    object$profile_covariates, object$group_covariates, first_rows,
    object$n_group_classes)
  same <- function(left, right) identical(unname(left), unname(right))
  if ((!is.null(object$indicator_data) &&
       !same(as.matrix(data[object$vars]), object$indicator_data)) ||
      !all(vapply(seq_along(designs$profile_design), function(index) {
        same(designs$profile_design[[index]], object$profile_design[[index]])
      }, logical(1))) || !same(designs$w, object$group_design) ||
      any(vapply(object$group_covariates, function(name) {
        any(data[[name]] != data[[name]][first_rows][object$group_index])
      }, logical(1)))) {
    stop(errorCondition("`data` must reproduce the original indicators and covariates.",
                        class = "multilpa_bad_inference_data", call = NULL))
  }
  x <- sweep(as.matrix(data[object$vars]), 2L, object$center, "-")
  theta <- .multilpa_cov_encode(object)
  objective <- function(parameters) {
    pieces <- .multilpa_cov_decode(parameters, object)
    -.multilpa_cov_expectation(x, object$group_index, pieces$parameters,
                               object$profile_design, object$group_design,
                               pieces$beta, pieces$gamma)$log_likelihood
  }
  reproduced <- -objective(theta)
  if (abs(reproduced - object$log_likelihood) >
      1e-6 * (1 + abs(object$log_likelihood))) {
    stop(errorCondition(
      "`data` do not reproduce the fitted log likelihood; supply the fitting data.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  gradient <- function(parameters) {
    -colSums(.multilpa_cov_group_scores(parameters, x, object))
  }
  ## Every block of the parameter vector is put on its own unit before the
  ## finite-difference Hessian sees it, or the information matrix's condition
  ## number tracks the units the indicators and covariates happened to arrive
  ## in and is then declared singular for no other reason.
  scale <- c(as.vector(t(sqrt(object$variances))),
             rep(1, length(theta) - length(object$means)))
  ## The off-diagonal Cholesky coordinates carry indicator units; the logged
  ## diagonal does not. This is the same scaling the covariate-free path in
  ## `parameter_inference()` applies, shared through one helper.
  dimension <- length(.multilpa_continuous_names(object))
  if (.multilpa_cov_is_full(object) && dimension > 0L) {
    covariance_scale <- .multilpa_covariance_coordinate_scale(
      object$variances, dimension, .multilpa_cov_spread_rows(object))
    scale[length(object$means) + seq_along(covariance_scale)] <- covariance_scale
  }
  ## A coefficient's unit is the reciprocal of its predictor's.
  profile_scale <- 1 / .multilpa_design_scale(do.call(rbind, object$profile_design))
  group_scale <- 1 / .multilpa_design_scale(object$group_design)
  membership <- length(theta) - length(object$profile_coefficients) -
    length(object$group_coefficients)
  scale[membership + seq_len(length(theta) - membership)] <-
    c(rep(profile_scale, ncol(object$profile_coefficients)),
      rep(group_scale, ncol(object$group_coefficients)))
  information <- .multilpa_observed_hessian(
    function(displacement) objective(theta + displacement * scale),
    function(displacement) gradient(theta + displacement * scale) * scale, scale, step)
  covariance <- information$inverse
  if (identical(vcov_type, "robust")) {
    scores <- sweep(.multilpa_cov_group_scores(theta, x, object), 2L, scale, "*")
    covariance <- covariance %*% .multilpa_cross_product(scores) %*% covariance
  }
  covariance <- covariance * tcrossprod(scale)
  dimnames(covariance) <- list(.multilpa_cov_parameter_names(object),
                               .multilpa_cov_parameter_names(object))
  covariance
}

#' Names for the free parameters of a covariate fit
#' @param object A fitted `multilpa_covariates` model.
#' @return A character vector, ordered as [.multilpa_cov_encode()].
#' @noRd
.multilpa_cov_parameter_names <- function(object) {
  ## One spelling for every class in the package: the same four tidy columns,
  ## in the same order, serialised by the same helper. The kind belongs in the
  ## name because a mean and a variance share a level, an outcome and a term,
  ## so leaving it out makes them indistinguishable and gives vcov() duplicate
  ## dimnames.
  .multilpa_parameter_names(.multilpa_cov_labels(object))
}

#' Covariance matrix of a covariate fit
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param data Optional. The data frame the model was fitted to; when omitted
#'   it is rebuilt from the indicators, group index and designs the fit stores.
#'   Supplying it is the stronger check that the caller still holds that frame.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` or `"robust"`.
#' @param scale Which parameter scale the covariance is on, matching
#'   [vcov.multilpa()]. `"natural"`, the default, is the covariance of the
#'   estimates [coef()] and [parameter_inference()] report: variances and
#'   residual covariances in their own units, carried from the estimation scale
#'   by the delta method. `"unconstrained"` is the covariance on the scale the
#'   model is estimated on, with log variances and log-Cholesky coordinates.
#'   Means and membership coefficients are the same on both scales.
#' @param ... Ignored, present for generic compatibility.
#' @return A square numeric matrix with one row and column per free parameter,
#'   named as [coef()] names them and ordered measurement means, measurement
#'   variances or covariances, profile logits, group logits.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(
#'   runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#' )
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                 n_group_classes = 2, profile_covariates = "x",
#'                 n_starts = 2, seed = 1)
#' vcov(fit, example_data, scale = "unconstrained")
#' @export
vcov.multilpa_covariates <- function(object, data = NULL, step = 1e-4,
                                     vcov_type = c("observed", "robust"),
                                     scale = c("natural", "unconstrained"), ...) {
  scale <- match.arg(scale)
  covariance <- .multilpa_cov_covariance(object, data, step, match.arg(vcov_type))
  if (identical(scale, "unconstrained")) {
    ## The kind in the name has to name the scale, or a log variance is served
    ## under a name that says `variance`.
    names <- .multilpa_parameter_names(.multilpa_cov_estimation_labels(object))
    dimnames(covariance) <- list(names, names)
    return(covariance)
  }
  theta <- .multilpa_cov_encode(object)
  jacobian <- .multilpa_cov_natural_jacobian(object, theta,
                                             .multilpa_cov_labels(object))
  natural <- jacobian %*% covariance %*% t(jacobian)
  dimnames(natural) <- dimnames(covariance)
  natural
}

#' Tidy labels of a covariate fit on the scale it is estimated on
#'
#' The same rows as [.multilpa_cov_labels()], with the spread block renamed to
#' the coordinate it really is: a log variance, or a Cholesky entry whose
#' diagonal is logged. This is the vocabulary [coef.multilpa()] already uses for
#' the covariate-free model, so both classes spell an estimation-scale
#' coordinate the same way.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @return A data frame with `level`, `outcome`, `term` and `parameter`.
#' @noRd
.multilpa_cov_estimation_labels <- function(object) {
  labels <- .multilpa_cov_labels(object)
  on_diagonal <- vapply(strsplit(labels$term, ":", fixed = TRUE), function(pair) {
    length(pair) == 2L && identical(pair[[1L]], pair[[2L]])
  }, logical(1))
  labels$parameter <- ifelse(
    labels$parameter == "variance", "log_variance",
    ifelse(labels$parameter == "covariance",
           ifelse(on_diagonal, "log_cholesky", "cholesky"), labels$parameter))
  labels
}

#' Assemble the tidy inference table
#' @param object The fit; `theta` the estimates; `covariance` their covariance.
#' @param level Confidence level. `vcov_type` is recorded as an attribute.
#' @return One row per free parameter, at every level.
#' @noRd
.multilpa_cov_inference_frame <- function(object, theta, covariance, level,
                                          vcov_type, adjust) {
  labels <- .multilpa_cov_labels(object)
  stopifnot("every parameter must be labelled" = nrow(labels) == length(theta))

  ## Variances and covariances are estimated as logs and log-Cholesky
  ## coordinates; one Jacobian returns every block, and the covariance it
  ## implies, to the natural units a reader can compare with the data.
  estimate <- .multilpa_cov_natural_estimate(object, theta, labels)
  jacobian <- .multilpa_cov_natural_jacobian(object, theta, labels)
  errors <- sqrt(pmax(diag(jacobian %*% covariance %*% t(jacobian)), 0))
  is_variance <- labels$parameter == "variance"
  is_covariance <- labels$parameter == "covariance"

  quantile <- stats::qnorm(1 - (1 - level) / 2)
  statistic <- estimate / errors
  result <- data.frame(
    labels, estimate = estimate, standard_error = errors,
    statistic = statistic,
    p_value = 2 * stats::pnorm(-abs(statistic)),
    conf_low = estimate - quantile * errors,
    conf_high = estimate + quantile * errors,
    row.names = NULL, stringsAsFactors = FALSE
  )
  ## A Wald test of a variance against zero is meaningless; the interval is not.
  ## The same applies to a covariance on the diagonal, which is a variance; an
  ## off-diagonal covariance is a real hypothesis and keeps its test.
  on_diagonal <- is_covariance &
    vapply(strsplit(labels$term, ":", fixed = TRUE), function(pair) {
      length(pair) == 2L && identical(pair[1L], pair[2L])
    }, logical(1))
  result$statistic[is_variance | on_diagonal] <- NA_real_
  result$p_value[is_variance | on_diagonal] <- NA_real_
  result <- .multilpa_adjust_p(result, adjust)
  attr(result, "vcov_type") <- vcov_type
  attr(result, "level") <- level
  result
}

#' Natural-unit estimates of a covariate fit's free parameters
#'
#' Means come back to input units, variances out of logs and covariances out of
#' log-Cholesky coordinates; the membership coefficients are already the logits
#' they are reported as.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param theta The estimation-scale parameter vector.
#' @param labels The tidy labels of `theta`, in the same order.
#' @return A numeric vector the same length as `theta`.
#' @noRd
.multilpa_cov_natural_estimate <- function(object, theta, labels) {
  estimate <- theta
  is_mean <- labels$parameter == "mean"
  is_variance <- labels$parameter == "variance"
  is_covariance <- labels$parameter == "covariance"
  estimate[is_mean] <- estimate[is_mean] +
    rep(object$center, times = object$n_profiles)
  estimate[is_variance] <- exp(theta[is_variance])
  if (any(is_covariance)) {
    dimension <- length(.multilpa_continuous_names(object))
    width <- dimension * (dimension + 1L) / 2L
    starts <- which(is_covariance)[seq(1L, sum(is_covariance), by = width)]
    estimate[is_covariance] <- unlist(lapply(starts, function(first) {
      .multilpa_cov_cholesky_natural(theta[first + seq_len(width) - 1L], dimension)
    }), use.names = FALSE)
  }
  estimate
}

#' Delta-method Jacobian from estimation to natural coordinates
#'
#' A mean is a location shift and a membership coefficient is reported as
#' estimated, so both differentiate to one. A variance is `exp()` of its
#' coordinate. A covariance block needs a whole matrix rather than a scalar
#' derivative, because every entry depends on several log-Cholesky coordinates.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param theta The estimation-scale parameter vector.
#' @param labels The tidy labels of `theta`, in the same order.
#' @return A square matrix, rows natural coordinates and columns estimation
#'   coordinates, in the order [.multilpa_cov_encode()] packs them.
#' @noRd
.multilpa_cov_natural_jacobian <- function(object, theta, labels) {
  jacobian <- diag(1, length(theta))
  is_variance <- labels$parameter == "variance"
  is_covariance <- labels$parameter == "covariance"
  diag(jacobian)[is_variance] <- exp(theta[is_variance])
  if (any(is_covariance)) {
    dimension <- length(.multilpa_continuous_names(object))
    width <- dimension * (dimension + 1L) / 2L
    starts <- which(is_covariance)[seq(1L, sum(is_covariance), by = width)]
    invisible(lapply(starts, function(first) {
      at <- first + seq_len(width) - 1L
      jacobian[at, at] <<- .multilpa_cov_cholesky_jacobian(theta[at], dimension)
    }))
  }
  jacobian
}

#' The covariance a log-Cholesky block implies, in vech order
#' @param values Log-Cholesky coordinates for one profile.
#' @param dimension Number of continuous indicators.
#' @return The lower triangle of the implied covariance matrix.
#' @noRd
.multilpa_cov_cholesky_natural <- function(values, dimension) {
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  factor <- matrix(0, dimension, dimension)
  factor[lower] <- values
  diag(factor) <- exp(diag(factor))
  tcrossprod(factor)[lower]
}

#' Jacobian of the log-Cholesky to covariance map
#'
#' Computed by central differences on an explicit, cheap and smooth map, rather
#' than derived by hand: the derivative of a matrix product through a logged
#' diagonal is easy to get subtly wrong, and the map costs one small
#' multiplication to evaluate.
#'
#' @param values Log-Cholesky coordinates for one profile.
#' @param dimension Number of continuous indicators.
#' @return A square matrix of derivatives of the covariance entries with
#'   respect to the coordinates, both in vech order.
#' @noRd
.multilpa_cov_cholesky_jacobian <- function(values, dimension) {
  step <- 1e-6
  columns <- lapply(seq_along(values), function(index) {
    up <- values; down <- values
    up[index] <- up[index] + step
    down[index] <- down[index] - step
    (.multilpa_cov_cholesky_natural(up, dimension) -
       .multilpa_cov_cholesky_natural(down, dimension)) / (2 * step)
  })
  matrix(unlist(columns, use.names = FALSE), length(values), length(values))
}

#' Level, outcome, term and kind for every free parameter
#' @param object A fitted `multilpa_covariates` model.
#' @return A data frame in the order [.multilpa_cov_encode()] packs them.
#' @noRd
.multilpa_cov_labels <- function(object) {
  vars <- object$vars
  profiles <- paste0("profile_", seq_len(object$n_profiles))
  variance_rows <- if (identical(object$variance_model, "equal")) "shared" else profiles
  block <- function(term, outcome, level, parameter) {
    grid <- expand.grid(term = term, outcome = outcome, stringsAsFactors = FALSE)
    data.frame(level = rep(level, nrow(grid)), outcome = grid$outcome, term = grid$term,
               parameter = rep(parameter, nrow(grid)), stringsAsFactors = FALSE)
  }
  continuous <- .multilpa_continuous_names(object)
  spread <- if (.multilpa_cov_is_full(object)) {
    positions <- which(lower.tri(matrix(0, length(continuous), length(continuous)),
                                 diag = TRUE), arr.ind = TRUE)
    pairs <- sprintf("%s:%s", continuous[positions[, 1L]], continuous[positions[, 2L]])
    block(pairs, variance_rows, "measurement", "covariance")
  } else block(continuous, variance_rows, "measurement", "variance")
  rbind(
    block(continuous, profiles, "measurement", "mean"),
    spread,
    block(rownames(object$profile_coefficients),
          colnames(object$profile_coefficients), "profile", "coefficient"),
    block(rownames(object$group_coefficients),
          colnames(object$group_coefficients), "group", "coefficient")
  )
}

#' Estimated parameters of a covariate fit
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param scale Natural coefficients, or the unconstrained coordinates the model
#'   is estimated on: log variances for a diagonal fit and log-Cholesky
#'   coordinates for a full-covariance fit. Matches [coef.multilpa()].
#' @param ... Ignored, present for generic compatibility.
#' @return A named numeric vector of every free parameter, in the order
#'   [parameter_inference()] reports them: measurement means, measurement
#'   variances or residual covariances, profile logits, then group-class logits.
#'   Names are `level.parameter.outcome.term`, the same four-part decomposition
#'   [parameter_inference()] reports as columns and the same spelling every
#'   fitted class in this package uses; it is what keeps a mean and a variance
#'   on the same indicator distinguishable, and the `parameter` part names the
#'   scale, so a log variance is never served under a name that says
#'   `variance`. On the natural scale variances and residual covariances are in
#'   their own units and agree with [parameter_inference()]'s `estimate` column
#'   exactly.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(
#'   runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#' )
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                 n_group_classes = 2, profile_covariates = "x",
#'                 n_starts = 2, seed = 1)
#' coef(fit)
#' @export
coef.multilpa_covariates <- function(object, scale = c("natural", "unconstrained"),
                                     ...) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"))
  scale <- match.arg(scale)
  theta <- .multilpa_cov_encode(object)
  if (identical(scale, "unconstrained")) {
    return(stats::setNames(
      theta, .multilpa_parameter_names(.multilpa_cov_estimation_labels(object))))
  }
  ## Means are stored centred for the likelihood and reported in input units;
  ## the spread block is estimated as logs or log-Cholesky coordinates and comes
  ## back through the same map [parameter_inference()] reports it through, so
  ## the two can never disagree.
  stats::setNames(
    .multilpa_cov_natural_estimate(object, theta, .multilpa_cov_labels(object)),
    .multilpa_cov_parameter_names(object))
}

#' Wald confidence intervals for a covariate fit
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param parm Optional parameter names or indices; defaults to all of them.
#' @param level Confidence level strictly between zero and one.
#' @param data Optional. The data frame the model was fitted to; when omitted
#'   it is rebuilt from the indicators, group index and designs the fit stores.
#'   Supplying it is the stronger check that the caller still holds that frame.
#' @param ... Passed to [parameter_inference()], so `vcov_type = "robust"` and
#'   `step` reach it.
#' @return A two-column matrix of Wald intervals, one row per requested
#'   parameter, named as [coef()] names them. Bounds are on the natural scale
#'   and are not constrained to respect a variance's positivity or a
#'   probability's range.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(
#'   runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#' )
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                 n_group_classes = 2, profile_covariates = "x",
#'                 n_starts = 2, seed = 1)
#' confint(fit, data = example_data)
#' @export
confint.multilpa_covariates <- function(object, parm, level = 0.95, data = NULL,
                                        ...) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"),
            "`data` must be a data frame when supplied" =
              is.null(data) || is.data.frame(data))
  inference <- parameter_inference(object, data, level = level, ...)
  names_all <- .multilpa_cov_parameter_names(object)
  intervals <- cbind(inference$conf_low, inference$conf_high)
  dimnames(intervals) <- list(names_all, paste0(
    format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%"))
  if (missing(parm)) return(intervals)
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(parm < 1) || any(parm > length(names_all)) ||
        any(parm != floor(parm))) {
      stop("`parm` must identify existing parameters by name or index.")
    }
    parm <- names_all[parm]
  }
  if (!is.character(parm) || anyNA(parm) || !all(parm %in% names_all)) {
    stop("`parm` must identify existing parameters by name or index.")
  }
  intervals[parm, , drop = FALSE]
}
