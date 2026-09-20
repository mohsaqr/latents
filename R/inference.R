#' Serialise a tidy parameter decomposition into one name per parameter
#'
#' Every class in this package decomposes a parameter the same way --
#' `level`, `parameter`, `outcome`, `term` -- and those four columns are what
#' [parameter_inference()] reports. A `coef()` name is that decomposition
#' written on one line, dot separated and in that order, so a name can always be
#' read back into the columns and two fits of different families spell the same
#' concept the same way. `level`, `parameter` and `outcome` never contain a dot,
#' so the first three separators are unambiguous and everything after the third
#' is the term. A term that does not apply -- a group-class probability has no
#' term -- is left off rather than spelled `NA`.
#'
#' @param labels A data frame with `level`, `parameter`, `outcome` and `term`.
#' @return A character vector, one name per row of `labels`.
#' @noRd
.multilpa_parameter_names <- function(labels) {
  stopifnot(
    "`labels` must carry level, parameter, outcome and term" =
      is.data.frame(labels) &&
      all(c("level", "parameter", "outcome", "term") %in% names(labels))
  )
  stem <- paste(labels$level, labels$parameter, labels$outcome, sep = ".")
  ifelse(is.na(labels$term), stem, paste(stem, labels$term, sep = "."))
}

#' The tidy decomposition of every coefficient of a fitted mixture
#' @param object A fitted `multilpa` model.
#' @param scale Natural estimates or unconstrained fitting coordinates.
#' @return A data frame with `level`, `outcome`, `term` and `parameter`, one row
#'   per coefficient, in the order [.multilpa_coefficients()] returns them.
#' @noRd
.multilpa_coefficient_labels <- function(object, scale) {
  .multilpa_tidy_labels(names(.multilpa_coefficients_raw(object, scale)), object)
}

#' Parameterize a fitted mixture
#' @param object A fitted `multilpa` model.
#' @param scale Natural estimates or unconstrained fitting coordinates.
#' @return A named numeric parameter vector, named as
#'   [.multilpa_parameter_names()] spells them.
#' @noRd
.multilpa_coefficients <- function(object, scale) {
  values <- .multilpa_coefficients_raw(object, scale)
  stats::setNames(unname(values),
                  .multilpa_parameter_names(.multilpa_tidy_labels(names(values), object)))
}

#' Parameterize a fitted mixture, in the internal generation grammar
#'
#' The generated `kind[index,...]` labels are an implementation detail: they are
#' compact to build alongside the values and they are parsed straight back into
#' tidy columns by [.multilpa_tidy_labels()]. Nothing user facing sees them;
#' [.multilpa_coefficients()] renames them before they leave the package.
#'
#' @param object A fitted `multilpa` model.
#' @param scale Natural estimates or unconstrained fitting coordinates.
#' @return A named numeric parameter vector.
#' @noRd
.multilpa_coefficients_raw <- function(object, scale) {
  stopifnot(inherits(object, "multilpa"), scale %in% c("natural", "unconstrained"))
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  ## Only the Gaussian indicators carry means and variances; a mixed or
  ## all-categorical fit has fewer of these than it has indicators.
  continuous <- .multilpa_continuous_names(object)
  n_indicators <- length(continuous)
  mean_names <- as.vector(t(outer(seq_len(n_profiles), continuous,
    function(profile, indicator) sprintf("mean[%s,%s]", profile, indicator))))
  variance_names <- if (object$variance_model == "equal") {
    sprintf("variance[shared,%s]", continuous)
  } else sub("^mean", "variance", mean_names)
  variance_values <- if (object$variance_model == "equal") object$variances[1L, ] else
    as.vector(t(object$variances))
  profile_indices <- if (scale == "natural") seq_len(n_profiles) else seq_len(n_profiles - 1L)
  group_indices <- if (scale == "natural") seq_len(n_types) else seq_len(n_types - 1L)
  profile_names <- as.vector(t(outer(seq_len(n_types), profile_indices,
    function(id, profile) sprintf("profile_probability[%s,%s]", id, profile))))
  group_names <- sprintf("group_probability[%s]", group_indices)
  profile_values <- as.vector(t(object$profile_probabilities[, profile_indices, drop = FALSE]))
  group_values <- object$group_probabilities[group_indices]
  if (scale == "unconstrained") {
    variance_values <- log(variance_values)
    variance_names <- sub("^variance", "log_variance", variance_names)
    profile_values <- as.vector(t(log(object$profile_probabilities[, profile_indices, drop = FALSE]) -
      log(object$profile_probabilities[, n_profiles])))
    group_values <- log(group_values) - log(object$group_probabilities[n_types])
    profile_names <- sub("^profile_probability", "profile_logit", profile_names)
    group_names <- sub("^group_probability", "group_logit", group_names)
  }
  if (identical(object$covariance_model, "full")) {
    covariance_coordinates <- .multilpa_covariance_coordinates(object, scale)
    variance_values <- unname(covariance_coordinates)
    variance_names <- names(covariance_coordinates)
  }
  response <- .multilpa_response_coordinates(object, scale)
  stats::setNames(c(as.vector(t(object$means)), variance_values,
                   unname(response), profile_values, group_values),
                 c(mean_names, variance_names, names(response),
                   profile_names, group_names))
}

#' Encode categorical response probabilities
#'
#' Natural coordinates are the probabilities themselves; unconstrained
#' coordinates are multinomial logits against each indicator's last category,
#' which is the same simplex device the mixing probabilities already use.
#'
#' @param object A fitted `multilpa` model.
#' @param scale `"natural"` or `"unconstrained"`.
#' @return A named numeric vector, empty when no indicator is categorical.
#' @noRd
.multilpa_response_coordinates <- function(object, scale) {
  blocks <- object$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) return(stats::setNames(numeric(0), character(0)))
  prefix <- if (scale == "natural") "response" else "response_logit"
  unlist(lapply(names(blocks), function(indicator) {
    block <- blocks[[indicator]]
    categories <- colnames(block) %||% as.character(seq_len(ncol(block)))
    keep <- if (scale == "natural") seq_len(ncol(block)) else seq_len(ncol(block) - 1L)
    values <- if (scale == "natural") block[, keep, drop = FALSE] else
      log(block[, keep, drop = FALSE]) - log(block[, ncol(block)])
    labels <- as.vector(t(outer(seq_len(nrow(block)), categories[keep],
      function(profile, category) sprintf("%s[%s,%s,%s]", prefix, profile,
                                          indicator, category))))
    stats::setNames(as.vector(t(values)), labels)
  }), use.names = TRUE)
}

#' Rebuild response probabilities from multinomial logits
#' @param theta Unconstrained values for the categorical block, in encode order.
#' @param object The fit supplying the shapes and names.
#' @return A named list of profiles-by-categories probability matrices.
#' @noRd
.multilpa_response_decode <- function(theta, object) {
  blocks <- object$response_probabilities
  at <- 0L
  result <- lapply(blocks, function(block) {
    n_free <- nrow(block) * (ncol(block) - 1L)
    logits <- cbind(matrix(theta[at + seq_len(n_free)], nrow(block),
                           ncol(block) - 1L, byrow = TRUE), 0)
    at <<- at + n_free
    probabilities <- exp(sweep(logits, 1L, .multilpa_log_sum_exp(logits), "-"))
    dimnames(probabilities) <- dimnames(block)
    probabilities
  })
  names(result) <- names(blocks)
  result
}

#' How many free parameters the categorical block contributes
#' @param object A fitted `multilpa` model.
#' @return A single integer, zero when no indicator is categorical.
#' @noRd
.multilpa_response_width <- function(object) {
  blocks <- object$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) return(0L)
  sum(vapply(blocks, function(block) nrow(block) * (ncol(block) - 1L), numeric(1)))
}

#' Encode unrestricted Gaussian covariance matrices
#' @param object A fitted full-covariance model.
#' @param scale Natural covariance or log-Cholesky coordinates.
#' @return A named vector of nonredundant covariance parameters.
#' @noRd
.multilpa_covariance_coordinates <- function(object, scale) {
  stopifnot(inherits(object, "multilpa"), identical(object$covariance_model, "full"),
            scale %in% c("natural", "unconstrained"))
  dimension <- length(.multilpa_continuous_names(object))
  if (dimension == 0L) return(stats::setNames(numeric(0), character(0)))
  profiles <- if (object$variance_model == "equal") 1L else seq_len(object$n_profiles)
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  positions <- which(lower, arr.ind = TRUE)
  unlist(lapply(profiles, function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance <- matrix(object$covariances[, , profile], dimension, dimension)
    values <- if (scale == "natural") covariance else t(chol(covariance))
    if (scale == "unconstrained") diag(values) <- log(diag(values))
    prefix <- if (scale == "natural") "covariance" else "cholesky"
    labels <- sprintf("%s[%s,%s,%s]", prefix,
      if (object$variance_model == "equal") "shared" else profile,
      .multilpa_continuous_names(object)[positions[, 1L]],
      .multilpa_continuous_names(object)[positions[, 2L]])
    if (scale == "unconstrained") {
      labels[positions[, 1L] == positions[, 2L]] <- sub("^cholesky", "log_cholesky",
        labels[positions[, 1L] == positions[, 2L]])
    }
    stats::setNames(values[lower], labels)
  }), use.names = TRUE)
}

#' Widths of the coefficient blocks, in the order they are encoded
#'
#' [.multilpa_coefficients_raw()] lays a fit out as means, spread, categorical
#' responses, profile mixing and group mixing, in that order, on either scale.
#' Counting those widths arithmetically -- rather than encoding the fit and
#' measuring the result -- is what lets the free-coordinate map below be built
#' without a Cholesky factorisation on every likelihood evaluation.
#'
#' @param object A fitted `multilpa` model.
#' @param scale `"natural"` or `"unconstrained"`.
#' @return A named numeric vector with `means`, `variances`,
#'   `response_probabilities`, `profile` and `group`, whose sum is the length of
#'   [.multilpa_coefficients()] on that scale.
#' @noRd
.multilpa_coordinate_widths <- function(object, scale) {
  stopifnot("`object` must be a fitted multilpa model" = inherits(object, "multilpa"),
            "`scale` must be natural or unconstrained" =
              length(scale) == 1L && scale %in% c("natural", "unconstrained"))
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  n_continuous <- length(.multilpa_continuous_names(object))
  spread <- if (identical(object$covariance_model, "full")) {
    n_continuous * (n_continuous + 1L) / 2L
  } else n_continuous
  blocks <- object$response_probabilities %||% list()
  c(means = n_profiles * n_continuous,
    variances = spread * if (object$variance_model == "equal") 1L else n_profiles,
    response_probabilities = if (scale == "natural") {
      sum(vapply(blocks, length, numeric(1)))
    } else .multilpa_response_width(object),
    profile = n_types * if (scale == "natural") n_profiles else n_profiles - 1L,
    group = if (scale == "natural") n_types else n_types - 1L)
}

#' The coordinates a fit actually estimates
#'
#' A `fixed` fit holds whole measurement blocks at the values `start` supplied,
#' so those coordinates are constants of the likelihood, not parameters. They
#' must still enter the likelihood -- the model is conditional on them, not
#' without them -- but they carry no score, no row of the information matrix and
#' no row of the inference table. Because `fixed` names whole blocks and the
#' blocks are contiguous in the encoding, the free set is simply the complement
#' of the held ranges; there is no partially-held block whose free remainder
#' would need a separate identification argument.
#'
#' @param object A fitted `multilpa` model.
#' @param scale `"natural"` or `"unconstrained"`.
#' @return An integer vector of positions into [.multilpa_coefficients()] on that
#'   scale, ascending, naming the coordinates this fit estimated.
#' @noRd
.multilpa_free_index <- function(object, scale) {
  widths <- .multilpa_coordinate_widths(object, scale)
  total <- as.integer(sum(widths))
  held <- object$fixed %||% character()
  if (length(held) == 0L) return(seq_len(total))
  measurement <- c("means", "variances", "response_probabilities")
  unknown <- setdiff(held, measurement)
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "This fit holds %s, which is not a measurement block inference can restore.",
      paste(sprintf("`%s`", unknown), collapse = " or ")),
      class = "multilpa_bad_fixed", call = NULL))
  }
  offsets <- cumsum(c(0, widths[measurement]))
  frozen <- unlist(lapply(seq_along(measurement), function(position) {
    if (measurement[[position]] %in% held) {
      offsets[[position]] + seq_len(widths[[measurement[[position]]]])
    } else integer(0)
  }), use.names = FALSE)
  setdiff(seq_len(total), as.integer(frozen))
}

#' Restore the held coordinates around a free parameter vector
#'
#' The returned closure scatters free values into a template taken from the fit
#' itself, so a held block enters every likelihood, score and Hessian evaluation
#' at exactly the value it was held at.
#'
#' @param template Full-length unconstrained coordinates of the fit whose held
#'   values are to be restored; the centred fit, when the likelihood is
#'   evaluated on centred indicators.
#' @param free Positions of the free coordinates within `template`.
#' @return A function of the free values returning the full coordinate vector.
#' @noRd
.multilpa_restore_held <- function(template, free) {
  stopifnot("`template` must be numeric" = is.numeric(template),
            "`free` must index `template`" = is.numeric(free) &&
              length(free) <= length(template) &&
              (length(free) == 0L || (min(free) >= 1L && max(free) <= length(template))))
  template <- unname(template)
  function(values) {
    stopifnot("free values must match the free coordinate count" =
                is.numeric(values) && length(values) == length(free))
    full <- template
    full[free] <- values
    full
  }
}

#' Decode unconstrained mixture coordinates
#' @param theta Parameter vector with log variances and baseline-category
#'   logits, carrying every coordinate of the model including any held block.
#' @param object A fitted model defining the parameter dimensions.
#' @return The four parameter blocks used by the expectation step.
#' @noRd
.multilpa_decode <- function(theta, object) {
  stopifnot(is.numeric(theta), inherits(object, "multilpa"),
            "`theta` must carry one value per model coordinate, held blocks included" =
              length(theta) ==
              sum(.multilpa_coordinate_widths(object, "unconstrained")))
  theta <- unname(theta)
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  n_indicators <- length(.multilpa_continuous_names(object))
  n_means <- n_profiles * n_indicators
  n_covariance <- if (identical(object$covariance_model, "full"))
    n_indicators * (n_indicators + 1L) / 2L else n_indicators
  n_variances <- n_covariance * if (object$variance_model == "equal") 1L else n_profiles
  n_logits <- n_types * (n_profiles - 1L)
  means <- matrix(theta[seq_len(n_means)], n_profiles, byrow = TRUE)
  covariances <- NULL
  if (identical(object$covariance_model, "full") && n_indicators > 0L) {
    lower <- lower.tri(matrix(0, n_indicators, n_indicators), diag = TRUE)
    covariance_list <- lapply(seq_len(n_profiles), function(profile) {
      stopifnot(is.numeric(profile), length(profile) == 1L)
      block <- if (object$variance_model == "equal") 1L else profile
      factor <- matrix(0, n_indicators, n_indicators)
      factor[lower] <- theta[n_means + (block - 1L) * n_covariance + seq_len(n_covariance)]
      diag(factor) <- exp(diag(factor))
      tcrossprod(factor)
    })
    covariances <- array(unlist(covariance_list, use.names = FALSE),
      c(n_indicators, n_indicators, n_profiles))
    variances <- t(matrix(vapply(covariance_list, diag, numeric(n_indicators)),
      n_indicators, n_profiles))
  } else {
    variances <- matrix(exp(theta[n_means + seq_len(n_variances)]),
      n_profiles, n_indicators, byrow = TRUE)
  }
  n_response <- .multilpa_response_width(object)
  response_probabilities <- if (n_response == 0L) NULL else
    .multilpa_response_decode(theta[n_means + n_variances + seq_len(n_response)],
                              object)
  offset <- n_means + n_variances + n_response
  profile_logits <- cbind(matrix(theta[offset + seq_len(n_logits)],
    n_types, n_profiles - 1L, byrow = TRUE), 0)
  profile_probabilities <- exp(sweep(profile_logits, 1L,
    .multilpa_log_sum_exp(profile_logits), "-"))
  group_logits <- c(theta[offset + n_logits + seq_len(n_types - 1L)], 0)
  group_probabilities <- exp(group_logits - max(group_logits))
  parameters <- list(means = means, variances = variances,
    profile_probabilities = profile_probabilities,
    group_probabilities = group_probabilities / sum(group_probabilities))
  if (!is.null(covariances)) parameters$covariances <- covariances
  if (!is.null(response_probabilities)) {
    parameters$response_probabilities <- response_probabilities
  }
  parameters
}

#' Score contributions of the categorical response probabilities
#'
#' For one indicator the complete-data score of a multinomial logit is the
#' posterior-weighted difference between the category actually observed and the
#' probability the model gives it, exactly as for the mixing weights. Missing
#' codes contribute nothing.
#'
#' @param codes Integer code matrix, one column per categorical indicator.
#' @param posteriors Individual profile posteriors.
#' @param blocks The fitted response probabilities, one matrix per indicator.
#' @param group_index Group indices, or `NULL` for a pooled total.
#' @return A matrix with one row per group (or one row overall) and one column
#'   per free response parameter, in [.multilpa_coefficients()] order.
#' @noRd
.multilpa_response_scores <- function(codes, posteriors, blocks, group_index = NULL) {
  if (is.null(blocks) || length(blocks) == 0L) {
    rows <- if (is.null(group_index)) 1L else max(group_index)
    return(matrix(numeric(0), rows, 0L))
  }
  pieces <- lapply(seq_along(blocks), function(indicator) {
    block <- blocks[[indicator]]
    code <- codes[, indicator]
    observed <- !is.na(code)
    free <- seq_len(ncol(block) - 1L)
    columns <- lapply(seq_len(nrow(block)), function(profile) {
      weight <- posteriors[, profile]
      vapply(free, function(category) {
        indicator_value <- as.numeric(observed & code == category)
        weight * (indicator_value - observed * block[profile, category])
      }, numeric(length(code)))
    })
    do.call(cbind, columns)
  })
  contribution <- do.call(cbind, pieces)
  if (is.null(group_index)) {
    matrix(colSums(contribution), 1L, ncol(contribution))
  } else rowsum(contribution, group_index, reorder = FALSE)
}

#' Differentiate the observed mixture likelihood
#' @param theta Unconstrained parameter vector.
#' @param x Centered numeric indicator matrix, possibly containing NA.
#' @param object Fitted model defining dimensions and groups.
#' @return The gradient of the negative observed-data log likelihood.
#' @noRd
.multilpa_score <- function(theta, x, object, codes = NULL) {
  stopifnot(is.numeric(theta), is.matrix(x), inherits(object, "multilpa"))
  parameters <- .multilpa_decode(theta, object)
  expectation <- .multilpa_expectation(x, object$group_index, parameters, codes)
  if (identical(object$covariance_model, "full")) {
    measurement <- .multilpa_full_measurement_score(parameters, expectation, object)
    mean_score <- measurement$means
    variance_score <- measurement$covariances
  } else {
  observed <- !is.na(x)
  measurement_scores <- lapply(seq_len(object$n_profiles), function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
    residuals[!observed] <- 0
    weights <- expectation$subject_posteriors[, profile]
    list(means = colSums(sweep(residuals, 2L, parameters$variances[profile, ], "/") * weights),
      variances = 0.5 * colSums((sweep(residuals^2, 2L,
        parameters$variances[profile, ], "/") - observed) * weights))
  })
  mean_score <- unlist(lapply(measurement_scores, `[[`, "means"), use.names = FALSE)
  variance_score <- if (object$variance_model == "equal") {
    Reduce(`+`, lapply(measurement_scores, `[[`, "variances"))
  } else unlist(lapply(measurement_scores, `[[`, "variances"), use.names = FALSE)
  }
  profile_score <- unlist(lapply(seq_len(object$n_group_classes), function(id) {
    stopifnot(is.numeric(id), length(id) == 1L)
    counts <- colSums(expectation$joint[[id]])
    (counts - sum(counts) * parameters$profile_probabilities[id, ])[
      seq_len(object$n_profiles - 1L)]
  }), use.names = FALSE)
  group_score <- (colSums(expectation$group_posteriors) -
    object$n_groups * parameters$group_probabilities)[seq_len(object$n_group_classes - 1L)]
  response_score <- as.vector(.multilpa_response_scores(
    codes, expectation$subject_posteriors, parameters$response_probabilities))
  -c(mean_score, variance_score, response_score, profile_score, group_score)
}

#' Gaussian scores in log-Cholesky coordinates
#' @param parameters Decoded Gaussian and mixing parameters.
#' @param expectation Nested expectation step including Gaussian moments.
#' @param object A fitted full-covariance model.
#' @return Mean and log-Cholesky score vectors.
#' @noRd
.multilpa_full_measurement_score <- function(parameters, expectation, object) {
  stopifnot(is.list(parameters), is.list(expectation), inherits(object, "multilpa"),
            identical(object$covariance_model, "full"))
  dimension <- length(.multilpa_continuous_names(object))
  if (dimension == 0L) return(list(means = numeric(0), covariances = numeric(0)))
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  component_scores <- lapply(seq_len(object$n_profiles), function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance <- matrix(parameters$covariances[, , profile], dimension, dimension)
    precision <- chol2inv(chol(covariance))
    moments <- expectation$gaussian_moments[[profile]]
    weights <- expectation$subject_posteriors[, profile]
    residuals <- sweep(moments$expected, 2L, parameters$means[profile, ], "-")
    scatter <- crossprod(residuals, residuals * weights) +
      Reduce(`+`, lapply(moments$adjustments, function(adjustment) {
        stopifnot(is.list(adjustment))
        adjustment$covariance * sum(weights[adjustment$rows])
      }))
    list(means = as.vector(precision %*% colSums(residuals * weights)),
      covariance = 0.5 * (precision %*% scatter %*% precision - sum(weights) * precision))
  })
  profiles <- if (object$variance_model == "equal") 1L else seq_len(object$n_profiles)
  covariance_scores <- lapply(profiles, function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance_score <- if (object$variance_model == "equal") {
      Reduce(`+`, lapply(component_scores, `[[`, "covariance"))
    } else component_scores[[profile]]$covariance
    factor <- t(chol(matrix(parameters$covariances[, , profile], dimension, dimension)))
    factor_score <- 2 * covariance_score %*% factor
    diag(factor_score) <- diag(factor_score) * diag(factor)
    factor_score[lower]
  })
  list(means = unlist(lapply(component_scores, `[[`, "means"), use.names = FALSE),
       covariances = unlist(covariance_scores, use.names = FALSE))
}

#' Delta-method Jacobian for mixture parameters
#' @param object A fitted model.
#' @return Jacobian mapping unconstrained to natural coefficients.
#' @noRd
.multilpa_inference_jacobian <- function(object) {
  stopifnot(inherits(object, "multilpa"))
  natural <- .multilpa_coefficients(object, "natural")
  theta <- .multilpa_coefficients(object, "unconstrained")
  n_means <- length(object$means)
  full_covariance <- identical(object$covariance_model, "full")
  dimension <- length(.multilpa_continuous_names(object))
  n_covariance <- if (full_covariance) dimension * (dimension + 1L) / 2L else dimension
  n_variances <- if (dimension == 0L) 0L else
    n_covariance * if (object$variance_model == "equal") 1L else object$n_profiles
  n_measurement <- n_means + n_variances
  blocks <- object$response_probabilities %||% list()
  natural_response <- sum(vapply(blocks, length, numeric(1)))
  free_response <- .multilpa_response_width(object)
  natural_offset <- n_measurement + natural_response
  free_offset <- n_measurement + free_response
  jacobian <- matrix(0, length(natural), length(theta), dimnames = list(names(natural), names(theta)))
  diag(jacobian)[seq_len(n_means)] <- 1
  if (full_covariance && dimension > 0L) {
    lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
    positions <- which(lower, arr.ind = TRUE)
    invisible(lapply(seq_len(n_variances / n_covariance), function(profile) {
      stopifnot(is.numeric(profile), length(profile) == 1L)
      factor <- t(chol(matrix(object$covariances[, , profile], dimension, dimension)))
      block <- vapply(seq_len(n_covariance), function(index) {
        stopifnot(is.numeric(index), length(index) == 1L)
        derivative <- matrix(0, dimension, dimension)
        row <- positions[index, 1L]
        column <- positions[index, 2L]
        derivative[row, column] <- if (row == column) factor[row, column] else 1
        covariance_derivative <- derivative %*% t(factor) + factor %*% t(derivative)
        covariance_derivative[lower]
      }, numeric(n_covariance))
      indices <- n_means + (profile - 1L) * n_covariance + seq_len(n_covariance)
      jacobian[indices, indices] <<- block
    }))
  } else {
    diag(jacobian)[n_means + seq_len(n_variances)] <- natural[n_means + seq_len(n_variances)]
  }
  ## Each profile-by-indicator row is its own simplex, with the same
  ## diag(p) - p p' derivative the mixing weights use.
  natural_at <- n_measurement
  free_at <- n_measurement
  invisible(lapply(blocks, function(block) {
    invisible(lapply(seq_len(nrow(block)), function(profile) {
      probabilities <- block[profile, ]
      derivative <- diag(probabilities, nrow = length(probabilities)) -
        tcrossprod(probabilities)
      rows <- natural_at + seq_len(ncol(block))
      columns <- free_at + seq_len(ncol(block) - 1L)
      jacobian[rows, columns] <<- derivative[, seq_len(ncol(block) - 1L), drop = FALSE]
      natural_at <<- natural_at + ncol(block)
      free_at <<- free_at + ncol(block) - 1L
    }))
  }))
  invisible(lapply(seq_len(object$n_group_classes), function(id) {
    stopifnot(is.numeric(id), length(id) == 1L)
    probabilities <- object$profile_probabilities[id, ]
    block <- diag(probabilities, nrow = length(probabilities)) - tcrossprod(probabilities)
    natural_rows <- natural_offset + (id - 1L) * object$n_profiles + seq_len(object$n_profiles)
    free_columns <- free_offset + (id - 1L) * (object$n_profiles - 1L) + seq_len(object$n_profiles - 1L)
    jacobian[natural_rows, free_columns] <<- block[, seq_len(object$n_profiles - 1L), drop = FALSE]
  }))
  probabilities <- object$group_probabilities
  block <- diag(probabilities, nrow = length(probabilities)) - tcrossprod(probabilities)
  natural_rows <- natural_offset + object$n_group_classes * object$n_profiles + seq_len(object$n_group_classes)
  free_columns <- free_offset + object$n_group_classes * (object$n_profiles - 1L) + seq_len(object$n_group_classes - 1L)
  jacobian[natural_rows, free_columns] <- block[, seq_len(object$n_group_classes - 1L), drop = FALSE]
  jacobian
}

#' Observed-information inference for multilevel LPA
#'
#' Computes ordinary maximum-likelihood standard errors from the inverse
#' observed Hessian and delta-method Wald confidence intervals. These are local,
#' asymptotic estimates conditional on the selected numbers of classes. They are
#' not robust sandwich errors and do not account for model selection. Mixture
#' boundary solutions and unidentified Hessians do not admit this calculation.
#'
#' @param x A converged `multilpa` fit with inactive variance bounds.
#' @param data Optional. The data frame the model was fitted to; when omitted
#'   it is rebuilt from the indicators, identifiers and occasions the fit
#'   stores, which round-trip exactly. Supplying it is the stronger check that
#'   the caller still holds that frame.
#'   A supplied frame is checked exactly against the stored indicators and
#'   identifiers; for an older fit without stored indicators, the likelihood
#'   provides a weaker consistency check.
#' @param level Confidence level strictly between zero and one.
#' @param step Positive finite-difference step for the analytic score derivative.
#' @param vcov_type `"observed"` inverts the observed information.
#'   `"robust"` returns the Huber-White sandwich `A^-1 B A^-1`, where `B`
#'   accumulates the outer product of the per-group scores. Groups are the
#'   independent units, so robust errors relax the assumption that the Gaussian
#'   within-group model is correctly specified. They do not relax the assumption
#'   that groups are independent.
#' @param adjust Multiplicity correction applied to `p_value` to produce
#'   `p_adjusted`, one of the methods [stats::p.adjust()] accepts. The default
#'   `"none"` leaves the two columns equal: a correction changes what a p-value
#'   means, so it is applied only when it is asked for, and the method that was
#'   applied is recorded in the `adjust` attribute. The correction is taken
#'   over the tests the table actually reports; bounded parameters carry no test
#'   and do not count towards the family.
#' @return A base `data.frame` with one row per reported parameter and the
#'   columns `level` (`"measurement"`, `"profile"` or `"group"`), `outcome`,
#'   `term`, `parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
#'   `p_adjusted`, `conf_low` and `conf_high`. Estimates are on the natural
#'   scale: variances in their own units, probabilities as probabilities. A
#'   parameter whose null value is on the boundary of its space -- a variance, a
#'   probability, a residual variance on the covariance diagonal -- carries `NA`
#'   for `statistic`, `p_value` and `p_adjusted`, because a Wald test against
#'   that boundary is not a question worth asking; the interval still is. All
#'   class probabilities are reported; their sum constraints make the natural
#'   covariance singular. Wald intervals are not clipped to parameter bounds.
#'
#'   A fit made with `fixed` reports only the parameters it estimated: a held
#'   measurement block is a constant of this likelihood, so it contributes no
#'   row here and no row or column to the information matrix. The table is on
#'   the natural scale, where every class probability is reported, so it has one
#'   row per estimated natural coefficient: `n_parameters` rows plus one for
#'   each set of probabilities whose reference category the estimation scale
#'   drops. The held values still enter the likelihood at
#'   the values they were held at, so every standard error, interval and p-value
#'   is conditional on that measurement solution and does not propagate its
#'   uncertainty. See [fit_staged()] for what that conditioning means.
#'
#'   Diagnostics of the fit as a whole travel as attributes rather than as
#'   columns repeated down every row: `covariance` and `covariance_unconstrained`
#'   (natural and estimation-scale covariance of the estimates), `hessian`,
#'   `gradient`, `scaled_score`, `condition_ratio`, `level`, `step`, `vcov_type`,
#'   `fixed` (the measurement blocks this fit held, `character()` when none) and
#'   `adjust`. When `vcov_type` is `"robust"`, `group_scores` holds the
#'   groups-by-parameters score matrix and `scaling_correction` the MLR scaling
#'   correction factor `tr(A^-1 B) / q`, used for scaled likelihood-ratio
#'   difference tests.
#' @section Conditions:
#'   `multilpa_no_converge` for an unconverged fit, `multilpa_boundary_fit` when
#'   a variance, a response probability or a mixing probability sits at its
#'   bound, `multilpa_singular_information` when the observed information cannot
#'   be inverted, `multilpa_no_free_parameters` when `fixed` held every
#'   parameter, `multilpa_bad_inference_data` when a supplied `data` does not
#'   reproduce the fit, and `multilpa_too_few_groups` for `vcov_type = "robust"`
#'   with fewer independent groups than reported parameters. A fit whose score
#'   is still far from zero is reported with a `multilpa_unconverged` warning
#'   rather than refused.
#' @examples
#' set.seed(42)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 10),
#'   score_a = rnorm(100), score_b = rnorm(100)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' parameter_inference(fit)
#'
#' # Many tests in one table: name the correction, do not apply one by stealth.
#' parameter_inference(fit, adjust = "BH")
#' @export
parameter_inference <- function(x, data = NULL, level = 0.95, step = 1e-4,
                                vcov_type = c("observed", "robust"),
                                adjust = .multilpa_p_adjust_methods) {
  UseMethod("parameter_inference")
}

#' @rdname parameter_inference
#' @export
parameter_inference.multilpa <- function(x, data = NULL, level = 0.95, step = 1e-4,
                             vcov_type = c("observed", "robust"),
                             adjust = .multilpa_p_adjust_methods) {
  stopifnot(inherits(x, "multilpa"))
  data <- .multilpa_resolve_data(x, data)
  stopifnot(is.data.frame(data),
            is.numeric(level), length(level) == 1L, is.finite(level), level > 0, level < 1,
            is.numeric(step), length(step) == 1L, is.finite(step), step > 0)
  vcov_type <- match.arg(vcov_type)
  adjust <- match.arg(adjust)
  .multilpa_check_regularity(x, vcov_type)
  free <- .multilpa_free_index(x, "unconstrained")
  free_natural <- .multilpa_free_index(x, "natural")
  if (length(free) == 0L) {
    stop(errorCondition(
      "This fit holds every parameter it has, so there is nothing to report a standard error for.",
      class = "multilpa_no_free_parameters", call = NULL))
  }
  stopifnot("the free coordinates must match the parameters the fit counts" =
              length(free) == x$n_parameters)
  full_theta <- .multilpa_coefficients(x, "unconstrained")
  theta <- full_theta[free]
  prepared <- .multilpa_inference_matrix(x, data)
  observed <- prepared$x
  centered_object <- x
  centered_object$means <- sweep(x$means, 2L, prepared$centers, "-")
  centered_full <- .multilpa_coefficients(centered_object, "unconstrained")
  ## A held block is a constant of this likelihood, not a parameter: it is put
  ## back at its held value on every evaluation, so the model is conditional on
  ## the measurement rather than fitted without it.
  restore <- .multilpa_restore_held(centered_full, free)
  centered_theta <- unname(centered_full)[free]
  codes <- prepared$codes
  objective <- function(parameters) {
    stopifnot(is.numeric(parameters))
    -.multilpa_expectation(observed, x$group_index,
      .multilpa_decode(restore(parameters), x), codes)$log_likelihood
  }
  score <- function(parameters) {
    stopifnot(is.numeric(parameters))
    .multilpa_score(restore(parameters), observed, x, codes)[free]
  }
  fitted_likelihood <- -objective(centered_theta)
  if (abs(fitted_likelihood - x$log_likelihood) > 1e-8 * (1 + abs(x$log_likelihood))) {
    stop(errorCondition(
      "data do not reproduce the fitted log likelihood; supply the original fitting data.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  parameter_scale <- c(as.vector(t(sqrt(x$variances))),
                       rep(1, length(full_theta) - length(x$means)))
  dimension <- length(.multilpa_continuous_names(x))
  if (identical(x$covariance_model, "full") && dimension > 0L) {
    # Shared with the covariate path, which lacked this scaling and reported a
    # singular information for an indicator in large units. One rule, one place.
    covariance_scale <- .multilpa_covariance_coordinate_scale(
      x$variances, dimension,
      if (x$variance_model == "equal") 1L else seq_len(x$n_profiles))
    parameter_scale[length(x$means) + seq_along(covariance_scale)] <- covariance_scale
  }
  parameter_scale <- parameter_scale[free]
  scaled_objective <- function(displacement) {
    stopifnot(is.numeric(displacement))
    objective(centered_theta + displacement * parameter_scale)
  }
  scaled_gradient <- function(displacement) {
    stopifnot(is.numeric(displacement))
    score(centered_theta + displacement * parameter_scale) * parameter_scale
  }
  information <- .multilpa_observed_hessian(scaled_objective, scaled_gradient,
                                          parameter_scale, step)
  scaled_hessian <- information$scaled
  hessian <- information$natural
  condition_ratio <- information$condition_ratio
  scaled_inverse <- information$inverse
  group_scores <- NULL
  scaling_correction <- NA_real_
  if (identical(vcov_type, "robust")) {
    group_scores <- .multilpa_group_scores(restore(centered_theta), observed,
                                           x, codes)[, free, drop = FALSE]
    scaled_cross <- .multilpa_cross_product(sweep(group_scores, 2L, parameter_scale, "*"))
    scaling_correction <- sum(diag(scaled_inverse %*% scaled_cross)) / length(theta)
    scaled_inverse <- scaled_inverse %*% scaled_cross %*% scaled_inverse
  }
  covariance_unconstrained <- scaled_inverse * tcrossprod(parameter_scale)
  dimnames(hessian) <- dimnames(covariance_unconstrained) <- list(names(theta), names(theta))
  ## The natural coordinates of a held block depend only on the held
  ## unconstrained ones, so restricting the Jacobian to the free rows and
  ## columns is the delta method for the free parameters, not an approximation
  ## of it.
  jacobian <- .multilpa_inference_jacobian(x)[free_natural, free, drop = FALSE]
  covariance <- jacobian %*% covariance_unconstrained %*% t(jacobian)
  standard_errors <- sqrt(pmax(diag(covariance), 0))
  estimates <- .multilpa_coefficients(x, "natural")[free_natural]
  critical <- stats::qnorm((1 + level) / 2)
  intervals <- cbind(estimates - critical * standard_errors, estimates + critical * standard_errors)
  colnames(intervals) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%")
  gradient <- stats::setNames(score(centered_theta), names(theta))
  scaled_score <- max(abs(gradient * parameter_scale))
  ## How far the reported estimate sits from the stationary point, measured in
  ## its own standard errors: a score `g` against information `I` displaces the
  ## estimate by `g / I`, and the standard error is `sqrt(1 / I)`, so `g * SE`
  ## is that displacement in standard-error units. It is dimensionless and it
  ## accounts for the sample size, which a fixed cut on the scaled score does
  ## not -- at the default `tol = 1e-8` on 720 observations the scaled score is
  ## 1.19e-02 and used to warn, while the displacement it implies is 1.7e-03,
  ## under a fifth of one percent of a standard error and unable to move any
  ## reported figure. A genuinely loose fit still trips it: `tol = 1e-4` on the
  ## same data gives 5.0e-02, five percent of a standard error.
  standard_errors_unconstrained <- sqrt(pmax(diag(covariance_unconstrained), 0))
  score_displacement <- max(abs(gradient) * standard_errors_unconstrained)
  if (is.finite(score_displacement) && score_displacement > 0.01) {
    warning(warningCondition(sprintf(
      paste("The fitted likelihood still carries a score worth %.1f%% of a",
            "standard error; refit with a tighter `tol` before using Wald",
            "inference."), 100 * score_displacement),
      class = "multilpa_unconverged", call = NULL))
  }
  statistic <- unname(estimates / standard_errors)
  result <- data.frame(
    .multilpa_coefficient_labels(x, "natural")[free_natural, , drop = FALSE],
    estimate = unname(estimates), standard_error = unname(standard_errors),
    statistic = statistic, p_value = 2 * stats::pnorm(-abs(statistic)),
    conf_low = unname(intervals[, 1L]), conf_high = unname(intervals[, 2L]),
    row.names = NULL, stringsAsFactors = FALSE)
  ## A Wald test of a variance against zero is not a question worth asking; the
  ## interval still is.
  continuous <- .multilpa_continuous_names(x)
  bounded <- result$parameter %in% c("variance", "probability", "response") |
    (result$parameter == "covariance" &
       result$term %in% paste(continuous, continuous, sep = ":"))
  result$statistic[bounded] <- NA_real_
  result$p_value[bounded] <- NA_real_
  result <- .multilpa_adjust_p(result, adjust)
  ## Diagnostics of the fit, not of any one parameter, so they travel as
  ## attributes rather than as columns repeated down every row.
  attributes(result) <- c(attributes(result), list(
    covariance = covariance, covariance_unconstrained = covariance_unconstrained,
    hessian = hessian, gradient = gradient, scaled_score = scaled_score,
    score_displacement = score_displacement,
    condition_ratio = condition_ratio, level = level, step = step,
    vcov_type = vcov_type, group_scores = group_scores,
    scaling_correction = scaling_correction,
    fixed = x$fixed %||% character()))
  result
}

#' Extract multilevel LPA coefficients
#' @param object A fitted `multilpa` model.
#' @param scale Natural coefficients or unconstrained log variances (diagonal),
#'   log-Cholesky coordinates (full covariance), and baseline-category logits.
#' @param ... Reserved for generic compatibility.
#' @return A named numeric vector of every coefficient the model has, held ones
#'   included: a block `fixed` held is part of the model and is reported here,
#'   even though it has no standard error and no interval. Natural coefficients
#'   include all mixing probabilities; unconstrained coefficients exclude their
#'   reference categories. Names are `level.parameter.outcome.term`, the same four-part
#'   decomposition [parameter_inference()] reports as columns and the same
#'   spelling every fitted class in this package uses, so a name written for one
#'   fit means the same thing for another. A parameter with no term -- a
#'   group-class probability -- carries the first three parts only.
#'   [parameter_inference()] is the tidy form and the one to prefer.
#' @examples
#' set.seed(3)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 8),
#'   score_a = stats::rnorm(80), score_b = stats::rnorm(80)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' coef(fit)
#' @export
#' @importFrom stats coef
coef.multilpa <- function(object, scale = c("natural", "unconstrained"), ...) {
  stopifnot(inherits(object, "multilpa"))
  .multilpa_coefficients(object, match.arg(scale))
}

#' Extract multilevel LPA covariance estimates
#' @param object A fitted `multilpa` model.
#' @param data Optional. The data frame the model was fitted to; when omitted
#'   it is rebuilt from the indicators, identifiers and occasions the fit
#'   stores, which round-trip exactly. Supplying it is the stronger check that
#'   the caller still holds that frame.
#' @param scale Which parameter scale the covariance is on. `"natural"`, the
#'   default, is the covariance of the estimates [coef()] reports -- variances
#'   in their own units and probabilities as probabilities -- carried from the
#'   estimation scale by the delta method with
#'   `.multilpa_inference_jacobian()`. `"unconstrained"` is the covariance on
#'   the scale the model is actually estimated on: log variances for a diagonal
#'   fit, log-Cholesky coordinates for a full-covariance fit, and
#'   baseline-category logits for the mixing and response probabilities.
#' @param ... Additional arguments passed to [parameter_inference()].
#' @return A square numeric matrix with one row and column per *estimated*
#'   coefficient, named as [coef()] names them, on the scale `scale` asks for. A
#'   fit made with `fixed` held part of its measurement model at supplied
#'   values; those coefficients were not estimated here, so they carry no row or
#'   column. With `scale = "unconstrained"` the matrix is `n_parameters` square,
#'   for any fit; on the natural scale it is larger by one row and column for
#'   each set of probabilities whose reference category the estimation scale
#'   drops, and singular by construction, because each set of probabilities sums
#'   to one. [confint()] and the
#'   `conf_low`/`conf_high` columns of [parameter_inference()] are built from
#'   the natural-scale matrix.
#' @examples
#' set.seed(3)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 8),
#'   score_a = stats::rnorm(80), score_b = stats::rnorm(80)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' vcov(fit, data = example_data, scale = "unconstrained")
#' @export
#' @importFrom stats vcov
vcov.multilpa <- function(object, data = NULL, scale = c("natural", "unconstrained"), ...) {
  stopifnot(inherits(object, "multilpa"))
  scale <- match.arg(scale)
  # `data` is no longer required: a fit carries the columns it was built from,
  # and `parameter_inference()` falls back to them.
  information <- parameter_inference(object, data, ...)
  if (scale == "natural") attr(information, "covariance") else
    attr(information, "covariance_unconstrained")
}

#' Wald confidence intervals for multilevel LPA coefficients
#' @param object A fitted `multilpa` model.
#' @param parm Optional coefficient names or indices; defaults to every
#'   coefficient the fit estimated. Naming a coefficient that `fixed` held
#'   raises `multilpa_held_parameter`, because a held value has no interval.
#' @param level Confidence level strictly between zero and one.
#' @param data Optional, exactly as for [vcov()].
#' @param ... Additional arguments passed to [parameter_inference()].
#' @return A two-column matrix of Wald intervals on the natural scale, one row
#'   per requested coefficient and named as [coef()] names them. Bounds are not
#'   clipped to the probability or variance parameter space. For a fit made with
#'   `fixed`, the default rows are the estimated coefficients only and the
#'   intervals are conditional on the held measurement solution.
#' @examples
#' set.seed(3)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 8),
#'   score_a = stats::rnorm(80), score_b = stats::rnorm(80)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' confint(fit, data = example_data)
#' @export
#' @importFrom stats confint
confint.multilpa <- function(object, parm, level = 0.95, data = NULL, ...) {
  stopifnot(inherits(object, "multilpa"), is.numeric(level), length(level) == 1L,
            is.finite(level), level > 0, level < 1)
  estimates <- coef.multilpa(object)
  covariance <- vcov.multilpa(object, data = data, ...)
  ## A held measurement coefficient is reported by coef() -- it is part of the
  ## model -- but it was not estimated here and has no interval. The default is
  ## therefore what the covariance actually covers.
  estimated <- rownames(covariance)
  if (missing(parm)) parm <- estimated
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm)) || any(parm != floor(parm)) ||
        any(parm < 1) || any(parm > length(estimates))) stop("Invalid coefficient indices in parm.")
    parm <- names(estimates)[parm]
  }
  if (!is.character(parm) || anyNA(parm) || !all(parm %in% names(estimates))) {
    stop("parm must identify existing coefficient names or indices.")
  }
  held <- setdiff(parm, estimated)
  if (length(held) > 0L) {
    stop(errorCondition(sprintf(
      "%s %s held fixed by this fit, so %s carr%s no confidence interval.",
      paste(sprintf("`%s`", held), collapse = ", "),
      if (length(held) == 1L) "was" else "were",
      if (length(held) == 1L) "it" else "they",
      if (length(held) == 1L) "ies" else "y"),
      class = "multilpa_held_parameter", call = NULL))
  }
  standard_errors <- sqrt(pmax(diag(covariance)[parm], 0))
  critical <- stats::qnorm((1 + level) / 2)
  intervals <- cbind(estimates[parm] - critical * standard_errors,
                     estimates[parm] + critical * standard_errors)
  colnames(intervals) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%")
  intervals
}

#' Check that a fit admits Wald inference at all
#'
#' Boundary solutions, unconverged fits and zero mixing probabilities each put
#' the estimate where the asymptotic normal approximation does not hold, so
#' they are refused rather than reported with a number.
#'
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_check_regularity <- function(object, vcov_type) {
  if (identical(vcov_type, "robust") && object$n_groups <= object$n_parameters) {
    stop(errorCondition(sprintf(
      "Robust inference needs more groups than parameters; this fit has %d groups and %d parameters.",
      object$n_groups, object$n_parameters),
      class = "multilpa_too_few_groups", call = NULL))
  }
  if (!isTRUE(object$converged)) {
    stop(errorCondition("Inference requires a converged fit.",
                        class = "multilpa_no_converge", call = NULL))
  }
  ## A bound-active value that this fit held fixed is a constant of the
  ## conditional likelihood, not an estimate sitting on its boundary, so it does
  ## not disqualify the parameters that were actually estimated.
  held <- object$fixed %||% character()
  if (isTRUE(object$boundary) && !("variances" %in% held)) {
    stop(errorCondition("Wald inference is unavailable for a bound-active fit.",
                        class = "multilpa_boundary_fit", call = NULL))
  }
  response <- if ("response_probabilities" %in% held) NULL else
    unlist(object$response_probabilities, use.names = FALSE)
  if (length(response) > 0L &&
      any(response <= (object$min_probability %||% 0) * (1 + 1e-7))) {
    stop(errorCondition(
      "Wald inference is unavailable for bound-active categorical response probabilities.",
      class = "multilpa_boundary_fit", call = NULL))
  }
  if (any(object$profile_probabilities <= 0) || any(object$group_probabilities <= 0)) {
    stop(errorCondition(
      "Wald inference requires strictly positive mixing probabilities.",
      class = "multilpa_boundary_fit", call = NULL))
  }
  invisible(NULL)
}

#' Recover and centre the indicator matrix the fit was built on
#'
#' The Hessian is evaluated at the fitted estimates, so it is only meaningful
#' if `data` really is the data that produced them. Stored indicators are
#' compared exactly; older fits without them fall back to the likelihood check
#' the caller performs next.
#'
#' @return A list with the centred matrix `x` and the `centers` removed from it.
#' @noRd
.multilpa_inference_matrix <- function(object, data) {
  if (nrow(data) != object$n_observations ||
      !all(c(object$vars, object$id) %in% names(data)) ||
      anyDuplicated(names(data)) ||
      !all(vapply(data[, .multilpa_continuous_names(object), drop = FALSE],
                  is.numeric, logical(1)))) {
    stop(errorCondition(
      "data must contain the original numeric indicators and group column.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  group_index <- match(data[[object$id]], object$group_values)
  if (!identical(group_index, object$group_index)) {
    stop(errorCondition(
      "data must retain the original group identifiers and row order.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  ## Only the Gaussian indicators form the numeric matrix; the categorical ones
  ## are re-encoded from `data` rather than taken from the fit, so that supplying
  ## the wrong categorical columns is caught rather than silently accepted.
  continuous <- .multilpa_continuous_names(object)
  ## as.matrix() of a zero-column frame is logical, and the expectation demands
  ## a numeric matrix, so an all-categorical fit needs the mode forced.
  x <- if (length(continuous) == 0L) matrix(numeric(0), nrow(data), 0L) else
    as.matrix(data[, continuous, drop = FALSE])
  if (any(is.infinite(x)) || any(is.nan(x)) ||
      (anyNA(x) && !identical(object$missing, "fiml"))) {
    stop(errorCondition(
      "data contain unsupported missing or non-finite indicators.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  if (!is.null(object$indicator_data) && ncol(x) > 0L &&
      !identical(x, object$indicator_data)) {
    stop(errorCondition(
      "data must reproduce the original indicator data, including row order and names.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  codes <- NULL
  if (length(object$categorical %||% character()) > 0L) {
    encoded <- .multilpa_encode_categorical(data[, object$categorical, drop = FALSE])
    codes <- encoded$codes
    if (!identical(unname(codes), unname(object$categorical_data)) ||
        !identical(encoded$levels, object$categorical_levels)) {
      stop(errorCondition(
        "data must reproduce the original categorical indicators, including their categories and row order.",
        class = "multilpa_bad_inference_data", call = NULL))
    }
  }
  centers <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  list(x = if (ncol(x) > 0L) sweep(x, 2L, centers, "-") else x,
       centers = centers, codes = codes)
}

#' Differentiate the observed information and check that it is usable
#'
#' Curvature is tested after scaling, so that indicators on different units do
#' not by themselves make an identified model look singular.
#'
#' @return A list with the `scaled` and `natural` Hessians, the scaled
#'   `condition_ratio`, and the `inverse` of the scaled Hessian.
#' @noRd
.multilpa_observed_hessian <- function(scaled_objective, scaled_gradient,
                                     parameter_scale, step) {
  n_parameters <- length(parameter_scale)
  scaled <- tryCatch(
    stats::optimHess(rep(0, n_parameters), scaled_objective, gr = scaled_gradient,
                     control = list(ndeps = rep(step, n_parameters))),
    error = function(error) {
      stop(sprintf("Observed Hessian failed: %s", conditionMessage(error)))
    })
  natural <- scaled / tcrossprod(parameter_scale)
  eigenvalues <- eigen(natural, symmetric = TRUE, only.values = TRUE)$values
  scaled_eigenvalues <- eigen(scaled, symmetric = TRUE, only.values = TRUE)$values
  condition_ratio <- min(scaled_eigenvalues) / max(scaled_eigenvalues)
  if (any(!is.finite(eigenvalues)) || !is.finite(condition_ratio) ||
      min(scaled_eigenvalues) <= 0 || condition_ratio <= 1e-10) {
    stop(errorCondition(
      "Observed information is not positive definite or is numerically singular; Wald inference is unavailable.",
      class = "multilpa_singular_information", call = NULL))
  }
  inverse <- tryCatch(solve(scaled), error = function(error) {
    stop(sprintf("Observed information could not be inverted: %s",
                 conditionMessage(error)))
  })
  list(scaled = scaled, natural = natural, condition_ratio = condition_ratio,
       inverse = inverse)
}

#' Split a generated parameter name into tidy columns
#'
#' The names are machine-generated with a strict grammar -- `mean[profile,
#' indicator]`, `variance[profile,indicator]`, `covariance[profile,a,b]`,
#' `profile_probability[group_class,profile]`, `group_probability[group_class]`,
#' and their unconstrained counterparts `log_variance`, `cholesky`,
#' `log_cholesky`, `response_logit`, `profile_logit` and `group_logit` -- so
#' they can be taken apart reliably. Reporting the pieces as columns means a
#' caller never has to parse a label to learn which level a parameter belongs
#' to, which is the whole point of a tidy result.
#'
#' @param names Character vector of generated parameter names.
#' @param object The fit, supplying the profile and group-class counts.
#' @return A data frame with `level`, `outcome`, `term` and `parameter`, one row
#'   per name, in the order given.
#' @noRd
.multilpa_tidy_labels <- function(names, object) {
  stopifnot("`names` must be character" = is.character(names))
  prefix <- sub("\\[.*$", "", names)
  unknown <- setdiff(prefix, c(.multilpa_measurement_kinds,
                               "profile_probability", "profile_logit",
                               "group_probability", "group_logit"))
  stopifnot("every parameter name must use a known kind" = length(unknown) == 0L)
  inside <- sub("^[^\\[]*\\[", "", sub("\\]$", "", names))
  parts <- strsplit(inside, ",", fixed = TRUE)
  first <- vapply(parts, function(p) p[[1L]], character(1))
  second <- vapply(parts, function(p) if (length(p) >= 2L) p[[2L]] else NA_character_,
                   character(1))
  third <- vapply(parts, function(p) if (length(p) >= 3L) p[[3L]] else NA_character_,
                  character(1))
  profile_label <- function(value) ifelse(value == "shared", "shared",
                                          paste0("profile_", value))
  class_label <- function(value) paste0("group_class_", value)
  measurement <- prefix %in% .multilpa_measurement_kinds
  within_profile <- prefix %in% c("profile_probability", "profile_logit")
  data.frame(
    level = ifelse(measurement, "measurement",
                   ifelse(within_profile, "profile", "group")),
    outcome = ifelse(measurement, profile_label(first),
                     ifelse(within_profile,
                            profile_label(second), class_label(first))),
    term = ifelse(prefix %in% .multilpa_paired_kinds,
                  paste(second, third, sep = ":"),
                  ifelse(measurement, second,
                         ifelse(within_profile,
                                class_label(first), NA_character_))),
    parameter = ifelse(prefix %in% c("profile_probability", "group_probability"),
                       "probability",
                       ifelse(prefix %in% c("profile_logit", "group_logit"),
                              "logit", prefix)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' Parameter kinds that belong to the measurement model
#' @noRd
.multilpa_measurement_kinds <- c("mean", "variance", "log_variance",
                                 "covariance", "cholesky", "log_cholesky",
                                 "response", "response_logit")

#' Parameter kinds whose term names a pair, written `a:b`
#' @noRd
.multilpa_paired_kinds <- c("covariance", "cholesky", "log_cholesky",
                            "response", "response_logit")

#' Multiplicity corrections an inference verb accepts
#'
#' The same methods [stats::p.adjust()] knows, with `"none"` first so that
#' `match.arg()` makes the uncorrected column the default. Applying a correction
#' silently would change what every p-value in the table means.
#'
#' @noRd
.multilpa_p_adjust_methods <- c("none",
                                setdiff(stats::p.adjust.methods, "none"))

#' Add the adjusted p-value column an inference table reports
#'
#' The family is the tests the table actually carries: a bounded parameter has
#' no test, is `NA` in `p_value`, and must not inflate the correction for the
#' parameters that do.
#'
#' @param result An inference table carrying `statistic` and `p_value`.
#' @param adjust One of [.multilpa_p_adjust_methods].
#' @return `result` with `p_adjusted` inserted after `p_value` and the
#'   `adjust` attribute recorded.
#' @noRd
.multilpa_adjust_p <- function(result, adjust) {
  stopifnot("`result` must carry a `p_value` column" = "p_value" %in% names(result))
  result$p_adjusted <- stats::p.adjust(result$p_value, method = adjust)
  ordered <- c("level", "outcome", "term", "parameter", "estimate",
               "standard_error", "statistic", "p_value", "p_adjusted",
               "conf_low", "conf_high")
  stopifnot("the inference table must carry exactly the documented columns" =
              setequal(names(result), ordered))
  result <- result[, ordered, drop = FALSE]
  attr(result, "adjust") <- adjust
  result
}
