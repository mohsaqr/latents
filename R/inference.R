#' Escape delimiters within one component of a paired parameter term
#' @param x Character indicator or category labels.
#' @return Labels with reserved delimiters escaped; ordinary labels unchanged.
#' @noRd
.latents_parameter_component <- function(x) {
  stopifnot(is.character(x))
  substitutions <- c("%" = "%25", ":" = "%3A", "," = "%2C",
                     "[" = "%5B", "]" = "%5D", "\r" = "%0D")
  Reduce(function(value, delimiter) {
    gsub(delimiter, substitutions[[delimiter]], value, fixed = TRUE)
  }, names(substitutions), init = x)
}

#' Undo `.latents_parameter_component()`
#'
#' Every `%` in an escaped component starts one of its three-character codes,
#' so each code is found where it was written; `%25` is decoded last, so a
#' label that itself reads like a code (`"a%3Ab"`, escaped `"a%253Ab"`) comes
#' back unchanged.
#' @param x Escaped components.
#' @return The original labels.
#' @noRd
.latents_parameter_decode <- function(x) {
  stopifnot("`x` must be character" = is.character(x))
  codes <- c("%3A" = ":", "%2C" = ",", "%5B" = "[", "%5D" = "]",
             "%0D" = "\r", "%25" = "%")
  Reduce(function(value, code) gsub(code, codes[[code]], value, fixed = TRUE),
         names(codes), init = x)
}

#' One label as it is written inside a paired term
#'
#' The rule for every label a user reads in a `term` column or a coefficient
#' name: it is the label's own text. The single exception is a term that pairs
#' two labels -- a covariance `a:b`, a category `indicator:category`, a direct
#' effect `indicator:covariate`. There the pair is split at its one bare `:`, so
#' a `:` inside either label is written `%3A`, and `%` itself `%25`; otherwise
#' indicator `"a:b"` with category `"c"` and indicator `"a"` with category
#' `"b:c"` would share the term `a:b:c`. Commas and brackets are only delimiters
#' of the internal `kind[...]` grammar (`.latents_parameter_component()`), never
#' of a term, so they keep their own text.
#' @param x Character labels.
#' @return `x` with `%` and `:` escaped; other labels unchanged.
#' @noRd
.latents_term_component <- function(x) {
  stopifnot("`x` must be character" = is.character(x))
  gsub(":", "%3A", gsub("%", "%25", x, fixed = TRUE), fixed = TRUE)
}

#' The term that pairs two labels
#'
#' The one spelling of a paired term, used by every family that labels a
#' parameter and by every accessor and plot that looks one up, so a producer
#' and its consumers cannot drift apart.
#' @param first,second Character labels, recycled against each other.
#' @return A character vector of terms `first:second`.
#' @noRd
.latents_term <- function(first, second) {
  paste(.latents_term_component(as.character(first)),
        .latents_term_component(as.character(second)), sep = ":")
}

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
#'   per coefficient, in the order `.multilpa_coefficients()` returns them.
#' @noRd
.multilpa_coefficient_labels <- function(object, scale) {
  .multilpa_tidy_labels(names(.multilpa_coefficients_raw(object, scale)), object)
}

#' Parameterize a fitted mixture
#' @param object A fitted `multilpa` model.
#' @param scale Natural estimates or unconstrained fitting coordinates.
#' @return A named numeric parameter vector, named as
#'   `.multilpa_parameter_names()` spells them.
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
#' tidy columns by `.multilpa_tidy_labels()`. Nothing user facing sees them;
#' `.multilpa_coefficients()` renames them before they leave the package.
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
  ## Indicator names are escaped like every other label in the `kind[...]`
  ## grammar, so a comma or colon in a name cannot split it (review 2026-10-26,
  ## INF-4); `.multilpa_tidy_labels()` gives the text back.
  continuous <- .latents_parameter_component(.multilpa_continuous_names(object))
  n_indicators <- length(continuous)
  mean_names <- as.vector(t(outer(seq_len(n_profiles), continuous,
    function(profile, indicator) sprintf("mean[%s,%s]", profile, indicator))))
  shared <- .multilpa_spread_shared(object)
  variance_names <- if (shared) {
    sprintf("variance[shared,%s]", continuous)
  } else sub("^mean", "variance", mean_names)
  variance_values <- if (shared) object$variances[1L, ] else
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
  ## A constrained structure is estimated in its own chart: log volumes, log
  ## shapes and orientations, one coordinate per parameter it has.
  if (scale == "unconstrained" && .multilpa_uses_chart(object)) {
    chart <- .multilpa_chart_encode(object)$values
    variance_values <- unname(chart)
    variance_names <- names(chart)
  }
  response <- .multilpa_response_coordinates(object, scale)
  extra <- .latents_extra_coordinates(object, object$extra_data, n_profiles, scale)
  stats::setNames(c(as.vector(t(object$means)), variance_values,
                   unname(response), unname(extra), profile_values, group_values),
                 c(mean_names, variance_names, names(response), names(extra),
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
                                          .latents_parameter_component(indicator),
                                          .latents_parameter_component(category)))))
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
  shared <- .multilpa_spread_shared(object)
  profiles <- if (shared) 1L else seq_len(object$n_profiles)
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  positions <- which(lower, arr.ind = TRUE)
  unlist(lapply(profiles, function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    covariance <- matrix(object$covariances[, , profile], dimension, dimension)
    values <- if (scale == "natural") covariance else t(chol(covariance))
    if (scale == "unconstrained") diag(values) <- log(diag(values))
    prefix <- if (scale == "natural") "covariance" else "cholesky"
    names <- .latents_parameter_component(.multilpa_continuous_names(object))
    labels <- sprintf("%s[%s,%s,%s]", prefix,
      if (shared) "shared" else profile,
      names[positions[, 1L]], names[positions[, 2L]])
    if (scale == "unconstrained") {
      labels[positions[, 1L] == positions[, 2L]] <- sub("^cholesky", "log_cholesky",
        labels[positions[, 1L] == positions[, 2L]])
    }
    stats::setNames(values[lower], labels)
  }), use.names = TRUE)
}

#' Widths of the coefficient blocks, in the order they are encoded
#'
#' `.multilpa_coefficients_raw()` lays a fit out as means, spread, categorical
#' responses, profile mixing and group mixing, in that order, on either scale.
#' Counting those widths arithmetically -- rather than encoding the fit and
#' measuring the result -- is what lets the free-coordinate map below be built
#' without a Cholesky factorisation on every likelihood evaluation.
#'
#' @param object A fitted `multilpa` model.
#' @param scale `"natural"` or `"unconstrained"`.
#' @return A named numeric vector with `means`, `variances`,
#'   `response_probabilities`, `profile` and `group`, whose sum is the length of
#'   `.multilpa_coefficients()` on that scale.
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
  spread_width <- if (scale == "unconstrained" && .multilpa_uses_chart(object)) {
    .multilpa_structure_parameters(object$covariance_structure, n_profiles,
                                   n_continuous)
  } else spread * if (.multilpa_spread_shared(object)) 1L else n_profiles
  c(means = n_profiles * n_continuous,
    variances = spread_width,
    response_probabilities = if (scale == "natural") {
      sum(vapply(blocks, length, numeric(1)))
    } else .multilpa_response_width(object),
    extra = .latents_extra_width(object$extra_data, n_profiles),
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
#' @return An integer vector of positions into `.multilpa_coefficients()` on that
#'   scale, ascending, naming the coordinates this fit estimated.
#' @noRd
.multilpa_free_index <- function(object, scale) {
  widths <- .multilpa_coordinate_widths(object, scale)
  total <- as.integer(sum(widths))
  # A covariance pattern (local dependence) holds the other covariances at
  # zero; their coordinates, in the same lower-triangle order in both scales,
  # are constants of the likelihood like a held block.
  pattern <- object$covariance_pattern
  patterned <- if (is.null(pattern)) integer() else {
    lower <- lower.tri(pattern, diag = TRUE)
    zero <- which(!pattern[lower])
    width <- sum(lower)
    blocks <- widths[["variances"]] %/% width
    as.integer(widths[["means"]] + as.vector(outer(zero, (seq_len(blocks) - 1L) * width, "+")))
  }
  held <- object$fixed %||% character()
  if (length(held) == 0L) return(setdiff(seq_len(total), patterned))
  measurement <- c("means", "variances", "response_probabilities")
  unknown <- setdiff(held, measurement)
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "This fit holds %s, which is not a measurement block inference can restore.",
      paste(sprintf("`%s`", unknown), collapse = " or ")),
      class = "latents_bad_fixed", call = NULL))
  }
  offsets <- cumsum(c(0, widths[measurement]))
  frozen <- unlist(lapply(seq_along(measurement), function(position) {
    if (measurement[[position]] %in% held) {
      offsets[[position]] + seq_len(widths[[measurement[[position]]]])
    } else integer(0)
  }), use.names = FALSE)
  setdiff(seq_len(total), c(as.integer(frozen), patterned))
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
  n_variances <- .multilpa_coordinate_widths(object, "unconstrained")[["variances"]]
  n_logits <- n_types * (n_profiles - 1L)
  means <- matrix(theta[seq_len(n_means)], n_profiles, byrow = TRUE)
  covariances <- NULL
  if (.multilpa_uses_chart(object)) {
    spread <- .multilpa_chart_decode(theta[n_means + seq_len(n_variances)], object)
    variances <- spread$variances
    covariances <- spread$covariances
  } else if (identical(object$covariance_model, "full") && n_indicators > 0L) {
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
  n_extra <- .latents_extra_width(object$extra_data, n_profiles)
  extra_blocks <- if (n_extra == 0L) NULL else
    .latents_extra_decode(theta[n_means + n_variances + n_response + seq_len(n_extra)],
                          object$extra_data, n_profiles)
  offset <- n_means + n_variances + n_response + n_extra
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
  if (!is.null(extra_blocks)) parameters[names(extra_blocks)] <- extra_blocks
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
#'   per free response parameter, in `.multilpa_coefficients()` order.
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
    # A half-observed pair row observes a set of joint categories; its
    # complete-data indicator is replaced by its expectation given the set,
    # the joint probability renormalized within the set (review 2026-10-26
    # round 2, C2-2).
    partial <- .multilpa_partial_of(codes, indicator)
    half <- if (is.null(partial)) NULL else which(!is.na(partial$index))
    sets <- if (is.null(partial)) NULL else .multilpa_partial_sets(partial, ncol(block))
    columns <- lapply(seq_len(nrow(block)), function(profile) {
      weight <- posteriors[, profile]
      within_set <- if (is.null(partial)) NULL else {
        conditional <- sweep(sets, 2L, block[profile, ], "*")
        (conditional / rowSums(conditional))[partial$index[half], , drop = FALSE]
      }
      vapply(free, function(category) {
        indicator_value <- as.numeric(observed & code == category)
        if (is.null(partial)) {
          return(weight * (indicator_value - observed * block[profile, category]))
        }
        indicator_value[half] <- within_set[, category]
        seen <- observed
        seen[half] <- TRUE
        weight * (indicator_value - seen * block[profile, category])
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
  expectation <- .multilpa_expectation(x, object$group_index, parameters, codes,
                                       unname(object$sampling_weights),
                                       object$extra_data)
  if (.multilpa_uses_chart(object)) {
    measurement <- .multilpa_chart_scores(
      parameters, expectation, x, object,
      theta[length(parameters$means) +
              seq_len(.multilpa_coordinate_widths(object, "unconstrained")[["variances"]])])
    mean_score <- as.vector(measurement$means)
    variance_score <- as.vector(measurement$covariances)
  } else if (identical(object$covariance_model, "full")) {
    measurement <- .multilpa_full_measurement_score(parameters, expectation, object)
    mean_score <- measurement$means
    variance_score <- measurement$covariances
  } else {
  # The expectation removes the direct effects itself; the residuals here
  # must be of the indicators with them removed too.
  shifted <- .latents_direct_shift(x, object$extra_data, parameters)
  observed <- !is.na(x)
  measurement_scores <- lapply(seq_len(object$n_profiles), function(profile) {
    stopifnot(is.numeric(profile), length(profile) == 1L)
    residuals <- sweep(shifted, 2L, parameters$means[profile, ], "-")
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
    sum(object$sampling_weights %||% object$n_groups) *
      parameters$group_probabilities)[seq_len(object$n_group_classes - 1L)]
  response_score <- as.vector(.multilpa_response_scores(
    codes, expectation$subject_posteriors, parameters$response_probabilities))
  extra_score <- as.vector(.latents_extra_scores(
    object$extra_data, expectation$subject_posteriors, parameters, x = x))
  -c(mean_score, variance_score, response_score, extra_score, profile_score,
     group_score)
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
  n_variances <- .multilpa_coordinate_widths(object, "natural")[["variances"]]
  n_chart <- .multilpa_coordinate_widths(object, "unconstrained")[["variances"]]
  n_measurement <- n_means + n_variances
  blocks <- object$response_probabilities %||% list()
  natural_response <- sum(vapply(blocks, length, numeric(1)))
  free_response <- .multilpa_response_width(object)
  n_extra <- .latents_extra_width(object$extra_data, object$n_profiles)
  natural_offset <- n_measurement + natural_response + n_extra
  free_offset <- n_means + n_chart + free_response + n_extra
  jacobian <- matrix(0, length(natural), length(theta), dimnames = list(names(natural), names(theta)))
  if (n_extra > 0L) {
    jacobian[n_measurement + natural_response + seq_len(n_extra),
             n_means + n_chart + free_response + seq_len(n_extra)] <-
      .latents_extra_jacobian(object, object$extra_data, object$n_profiles)
  }
  diag(jacobian)[seq_len(n_means)] <- 1
  if (.multilpa_uses_chart(object)) {
    jacobian[n_means + seq_len(n_variances), n_means + seq_len(n_chart)] <-
      .multilpa_chart_natural_jacobian(object)
  } else if (full_covariance && dimension > 0L) {
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
  free_at <- n_means + n_chart
  invisible(lapply(blocks, function(block) {
    invisible(lapply(seq_len(nrow(block)), function(profile) {
      rows <- natural_at + seq_len(ncol(block))
      columns <- free_at + seq_len(ncol(block) - 1L)
      jacobian[rows, columns] <<- .inference_simplex_jacobian(block[profile, ], "last")
      natural_at <<- natural_at + ncol(block)
      free_at <<- free_at + ncol(block) - 1L
    }))
  }))
  invisible(lapply(seq_len(object$n_group_classes), function(id) {
    stopifnot(is.numeric(id), length(id) == 1L)
    natural_rows <- natural_offset + (id - 1L) * object$n_profiles + seq_len(object$n_profiles)
    free_columns <- free_offset + (id - 1L) * (object$n_profiles - 1L) + seq_len(object$n_profiles - 1L)
    jacobian[natural_rows, free_columns] <<-
      .inference_simplex_jacobian(object$profile_probabilities[id, ], "last")
  }))
  natural_rows <- natural_offset + object$n_group_classes * object$n_profiles + seq_len(object$n_group_classes)
  free_columns <- free_offset + object$n_group_classes * (object$n_profiles - 1L) + seq_len(object$n_group_classes - 1L)
  jacobian[natural_rows, free_columns] <-
    .inference_simplex_jacobian(object$group_probabilities, "last")
  jacobian
}

#' Standard errors, confidence intervals and tests for model parameters
#'
#' Returns a table of the estimated parameters of a fitted model with their
#' standard errors, Wald tests and confidence intervals. By default the
#' standard errors are maximum-likelihood errors from the inverse observed
#' information, with delta-method Wald intervals; `method = "bootstrap"`
#' resamples groups instead. The estimates are local and asymptotic,
#' conditional on the selected numbers of classes, and do not account for
#' model selection. Boundary solutions and unidentified Hessians are refused
#' by the Wald path.
#'
#' @param x A converged fit from [multilpa()], [lpa()], [lca()] or
#'   [multilca()] (with or without covariates), [lta()], the additive family,
#'   or a growth mixture. A fit with an active bound raises
#'   `latents_boundary_fit` and an unconverged fit `latents_no_converge`.
#'   Mixture regressions are refused with `latents_bad_argument` (their
#'   inference is in [get_results()], [confint()] and [vcov()]), and the
#'   cross-level families with `latents_unsupported_inference`.
#' @param data Optional data frame the model was fitted to. When omitted it is
#'   rebuilt from the indicators, identifiers and occasions the fit stores. A
#'   supplied frame is checked exactly against the stored indicators and
#'   identifiers; for an older fit without stored indicators, the likelihood
#'   gives a weaker consistency check.
#' @param level Confidence level, strictly between zero and one. Default
#'   `0.95`.
#' @param step Positive finite-difference step for the analytic score
#'   derivative. Default `1e-4`.
#' @param vcov_type Covariance estimator for the Wald path.
#'   `"observed"` (default) inverts the observed information.
#'   `"robust"` returns the Huber-White sandwich `A^-1 B A^-1`, where `B` is
#'   the outer product of the per-group scores, multiplied by the CR1
#'   small-sample factor `G / (G - 1)` for `G` groups (the default of
#'   `sandwich::vcovCL()`, applied to every robust standard error in the
#'   package). Mplus's MLR applies no such factor, so these errors are
#'   `sqrt(G / (G - 1))` times Mplus's. Groups are the independent units:
#'   robust errors relax the assumption that the within-group model is
#'   correctly specified, not the assumption that groups are independent.
#'   `"opg"` inverts the outer product of the per-group scores (the BHHH
#'   estimator), as `glca` reports; it agrees with `"observed"` in large
#'   samples and can differ in small ones, where `"observed"` is usually more
#'   reliable. `"robust"` and `"opg"` need more groups than parameters; `"opg"`
#'   also needs group scores of full rank, since it inverts their outer
#'   product, while the sandwich does not. For a
#'   fit with sampling weights the default is `"robust"`, and any other value
#'   raises `latents_unsupported_weights`. For a fit with `centering =
#'   "person"` or `scaling = "person"` the default is also `"robust"`, and
#'   `"observed"` or `"opg"` raise `latents_unsupported_inference`: a unit's
#'   transformed rows are not independent observations, and the inverse
#'   information understated their errors (by up to 22% in a simulation),
#'   while the sandwich over units is valid. Under `centering = "grand"` or
#'   `scaling = "grand"` the sample means and standard deviations the
#'   indicators were transformed with are estimated jointly with the model,
#'   as estimating equations, and their variability is added to every
#'   covariance type: a centred mean `mu_k - xbar` has the standard error of
#'   that difference, not that of `mu_k`. `vcov_type` and `step` are not used
#'   by `method = "bootstrap"`.
#' @param adjust Multiplicity correction applied to `p_value` to give
#'   `p_adjusted`: `"none"` (default) or one of the [stats::p.adjust()]
#'   methods. The method applied is stored in the `adjust` attribute. The
#'   correction is taken over the tests the table reports; bounded parameters
#'   carry no test and are not counted.
#' @param method `"wald"` (default) inverts the observed information in the
#'   estimation coordinates; it is available for all fourteen covariance
#'   structures. EEI, VVI, EEE and VVV are estimated in log variances or
#'   log-Cholesky coordinates; the ten structures that constrain volume,
#'   shape or orientation across profiles are estimated in log volumes, log
#'   shapes with the determinant-one constraint built into the coordinates,
#'   and Cayley-transform orientations (or, for VEE and EVV, determinant-one
#'   log-Cholesky factors), so the information has one row per free
#'   parameter. The reported variances and covariances are obtained from
#'   these coordinates by the delta method, and their natural-scale
#'   covariance is singular, as for the class probabilities. `"bootstrap"`
#'   resamples groups with replacement, refits within the same covariance
#'   structure, realigns each refit's labels, and reports the percentile
#'   interval and the standard deviation of the replicates; it is also
#'   available for all fourteen structures and assumes independent groups, as
#'   `vcov_type = "robust"` does.
#' @param iter Number of resamples for `method = "bootstrap"`, at least two.
#'   Default `199`; a 95% percentile interval is more stable with 999 or more.
#' @param n_starts,max_iter,tol Passed to each bootstrap refit, as in
#'   [multilpa()].
#' @param seed Optional seed for the resampling: a whole number accepted by
#'   `set.seed()`; any other value raises `latents_bad_argument` on either
#'   path. The Wald path uses no random numbers. A supplied integer seed restores
#'   the caller's `.Random.seed` on exit; `NULL` uses its current stream.
#'   R does not expose the Box-Muller cached normal
#'   variate, so that cache cannot be restored; see [base::Random].
#' @param boundary Handling of categorical response probabilities at their
#'   lower bound, `min_probability`, where the information is singular.
#'   `"error"` (default) raises `latents_boundary_fit`. `"fix"` holds those
#'   probabilities at the bound and reports every other parameter conditional
#'   on them, as Mplus does; the held probabilities are omitted from the table
#'   (or, for a covariate fit, reported with `NA` standard errors) and named in
#'   the `fixed_at_bound` attribute.
#' @return A base `data.frame` with one row per reported parameter and the
#'   columns `level` (`"measurement"`, `"profile"` or `"group"`), `outcome`,
#'   `term`, `parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
#'   `p_adjusted`, `conf_low` and `conf_high`. Estimates are on the natural
#'   scale (variances in their own units, probabilities as probabilities). A
#'   parameter whose null value lies on the boundary of its space (a variance,
#'   a probability, a covariance diagonal, a count mean) has `NA` `statistic`,
#'   `p_value` and `p_adjusted`, but still has an interval. All class
#'   probabilities are reported; their sum constraints make the natural
#'   covariance singular. Wald intervals for probabilities are formed on the
#'   logit scale and for variances (and covariance diagonals) and count means
#'   on the log scale, then
#'   back-transformed, so they stay within the parameter's range; every other
#'   interval is `estimate +/- z * standard_error`. `confint()` returns the
#'   same intervals. A probability at or below 1e-6 is treated as on its bound
#'   (refused, or held by `boundary = "fix"`).
#'
#'   A fit made with `fixed` reports only the parameters it estimated: a held
#'   measurement block contributes no row to the table or to the information
#'   matrix. The table has one row per estimated natural coefficient:
#'   `n_parameters` rows plus one for each set of probabilities whose
#'   reference category the estimation scale drops. Standard errors,
#'   intervals and p-values are conditional on the held measurement solution
#'   and do not include its uncertainty (see [fit_staged()]). For a
#'   [fit_staged()] result, `method = "bootstrap"` refits both stages in each
#'   resample, so the table includes the measurement parameters and that
#'   uncertainty is propagated (attribute `stages_resampled` is 2).
#'
#'   With `method = "bootstrap"` the table has the same columns and rows;
#'   `estimate` is the fitted value, `standard_error` the standard deviation
#'   of the replicates and `conf_low`/`conf_high` the percentile interval.
#'   `statistic`, `p_value` and `p_adjusted` are `NA`. Its attributes are
#'   `covariance` (the covariance of the replicates, returned by `vcov()`),
#'   `level`, `method`, `iter`, `n_valid`, `replicates` (the kept replicates,
#'   one row each), `messages` (why a resample was dropped, `NA` where it was
#'   not), `structure` and `fixed`; `covariance_unconstrained` is `NULL`.
#'
#'   For the Wald path, diagnostics of the whole fit are attributes:
#'   `covariance` and `covariance_unconstrained` (natural and
#'   estimation-scale covariance of the estimates), `hessian`, `gradient`,
#'   `scaled_score`, `score_displacement` (the largest score times its
#'   estimation-scale standard error, which the `latents_unconverged` warning
#'   uses), `condition_ratio`, `level`, `step`, `method`, `vcov_type`, `fixed`
#'   (the measurement blocks the fit held, `character()` when none),
#'   `fixed_at_bound` (the probabilities held by `boundary = "fix"`) and
#'   `adjust`. When `vcov_type` is `"robust"` or `"opg"`, `group_scores` holds
#'   the groups-by-parameters score matrix and `scaling_correction` the MLR
#'   scaling correction factor `tr(A^-1 B) / q` (with `B` without the CR1
#'   factor), used for scaled likelihood-ratio difference tests; it is `NA`
#'   otherwise.
#' @details For a covariate fit (`multilpa(profile_covariates = )`), the
#'   table covers every free parameter at all three levels: the measurement
#'   model, the profile logits and the group-class logits. `outcome` is the
#'   profile or group class a coefficient predicts, or the profile a
#'   measurement parameter belongs to, and `parameter` is `"mean"`,
#'   `"variance"`, `"covariance"`, `"response"` or `"coefficient"`. Each
#'   categorical indicator contributes one `"response"` row per profile,
#'   indicator and category, named `indicator:category`: the category
#'   probability, with its standard error from the multinomial logits and no
#'   Wald test. Class-membership coefficients are on the
#'   multinomial logit scale, with the final profile and the final group class
#'   as references, so each coefficient is a log odds against that reference.
#'   The tests are Wald tests, which are invalid for a coefficient driven to
#'   the boundary by separation or an empty class. A fit with a membership
#'   logit above 20 in coordinates independent of the covariates' units and
#'   origin (roughly, a slope per standard deviation of its covariate and an
#'   intercept at the covariates' means) raises `latents_boundary_fit`; the
#'   fit itself warns with `latents_extreme_coefficients`. Only
#'   `method = "wald"` is implemented for a covariate fit; `"bootstrap"`
#'   raises `latents_unsupported_inference`.
#'
#'   For a latent transition model from [lta()], the table covers the
#'   measurement model, the group-class shares, each group class's initial
#'   profile distribution and its transition matrix. An initial probability
#'   (`parameter = "initial_probability"`, `level = "profile"`) has `outcome`
#'   the profile and `term` the group class; a transition probability
#'   (`"transition_probability"`) has `outcome` the profile moved to and
#'   `term` `group_class_h:profile_j`, the group class and the profile moved
#'   from. The scores follow the Fisher identity and the observed information
#'   is the numerical Jacobian of that analytic score; `vcov_type = "robust"`
#'   clusters on the groups, which are the sequences. `method = "bootstrap"`
#'   resamples persons, refits and matches each replicate's profiles and group
#'   classes to the original's; it needs `data` (`latents_bad_argument`
#'   otherwise). A transition row with no information (a profile no group
#'   occupies before its last occasion) raises `latents_boundary_fit`.
#'
#'   For a growth mixture ([mixture_regression()] with `random`), the
#'   table has one row per natural parameter: `block` (`"coefficient"`,
#'   `"membership"` and, for a multilevel growth mixture,
#'   `"group_membership"`, `"residual_sd"`, `"random_variance"`,
#'   `"random_sd"`, `"random_correlation"`, `"share"`), `class` (`"all"` for
#'   a quantity shared by every class), `term`, `estimate`, `standard_error`,
#'   `statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`. Wald
#'   rows use the fit's covariance (`vcov_type`) with the delta method;
#'   intervals of standard deviations and variances are formed on the log
#'   scale, of correlations on the Fisher scale and of shares on the logit
#'   scale, and only coefficients and logits are tested against zero.
#'   `method = "bootstrap"` resamples persons, refits each resample from
#'   `n_starts` starts and matches its classes to the original by their
#'   posterior memberships of the original persons; it is refused for a
#'   multilevel growth mixture (`latents_unsupported_inference`). In a
#'   simulation of 1000 datasets per design, Wald intervals covered
#'   0.936-0.965 with well-separated classes and 0.885-0.943 with weakly
#'   separated ones, whose standard errors were 6-30% too small
#'   (validation/GROWTH_SIMULATION.md in the source repository); the
#'   bootstrap is preferable for weakly separated classes. A class
#'   whose random-effect covariance is degenerate, or (with correlated
#'   residuals) a residual standard deviation or residual process held on its
#'   bound, raises `latents_boundary_fit` under `boundary = "error"`, naming
#'   which; under `"fix"` those rows have no standard error or interval. `data`, when
#'   given, must be the fitted data; `step` is not used (the information is
#'   the derivative of the analytic score).
#'
#'   For a general transition fit (an [lta()] fit with covariate,
#'   occasion-varying or second-order transitions), the table has one row per
#'   transition, initial and second-order logit: `block`, `group_class`,
#'   `from`, `to`, `term`, `estimate`, `standard_error`, `statistic`,
#'   `p_value`, `p_adjusted`, `conf_low`, `conf_high`, `estimated` (and
#'   `previous` for a second-order fit). `estimated` is `FALSE` for a logit
#'   the likelihood does not identify: a move into an occasion at which no
#'   row observes an indicator, under `transitions = "occasion"` (or every
#'   second-order logit when nothing is observed from occasion 3 on). Such a
#'   logit is not counted in `n_parameters`, AIC or BIC, and its standard
#'   error and test are `NA`. `method = "bootstrap"` resamples persons, refits, and
#'   matches every replicate's profiles (on the measurement) and group classes
#'   to the original's, rebasing the initial logits when the reference profile
#'   moves; it needs `data`, and works where Wald inference is refused (for
#'   example full covariance structures).
#'
#'   Transition inference (both kinds of [lta()] fit) fixes no
#'   parameter at a bound: without an active bound `boundary = "fix"` gives
#'   the same table as `"error"`, and a fit with an active bound raises
#'   `latents_boundary_fit` under either.
#' @section Computation: The observed information takes two evaluations of
#'   the score over all the data for each free parameter, which makes it the
#'   slowest part of Wald inference. For profile, class, covariate and
#'   transition models it is computed once per R session for a fit, its data
#'   and `step`, and reused by later calls, including those made by the plots
#'   that draw intervals and by `get_results()`. Any change to the fit or its
#'   data computes it anew, and the result is the same either way.
#' @section Conditions:
#'   Errors: `latents_no_converge` for an unconverged fit;
#'   `latents_boundary_fit` when a variance, response probability or mixing
#'   probability is at its bound; `latents_singular_information` when the
#'   observed information cannot be inverted; `latents_no_free_parameters`
#'   when `fixed` held every parameter; `latents_bad_inference_data` when a
#'   supplied `data` does not reproduce the fit; `latents_too_few_groups` for
#'   `vcov_type = "robust"` or `"opg"` when the fit has no more independent
#'   groups than free parameters. The Wald path raises
#'   `latents_unsupported_inference` for a fit made with `multilpa(prior = )`
#'   (a posterior mode, where the score is not zero), and either path raises
#'   `latents_unsupported_noise` for a fit with a noise component. A fit whose
#'   score is still far from zero gives a `latents_unconverged` warning.
#'   An invalid `level`, `step`, `iter`, `n_starts`, `max_iter` or `tol`
#'   raises `latents_bad_argument`, and a `data` that is not a data frame
#'   `latents_bad_data`. A fitted mixture regression is refused with
#'   `latents_bad_argument` naming `get_results(x, "coefficients")`, which
#'   reports its standard errors and tests; an object that is not a fitted
#'   model raises `latents_not_a_fit`.
#'
#'   `method = "bootstrap"` also raises `latents_unsupported_inference` for an
#'   additive or covariate fit, and for a
#'   fit holding a measurement block through `multilpa(fixed = )`, since
#'   resampling the data does not resample the held values (a [fit_staged()]
#'   result is bootstrapped by refitting both stages instead), and
#'   `latents_bootstrap_failed` when fewer than two resamples gave a usable
#'   fit; its message gives the first reason. Resamples that fail or do not
#'   converge are dropped with a `latents_bootstrap_dropped` warning stating
#'   how many.
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
#'
#' # A structure that constrains the shape across profiles (VEI) is charted in
#' # its own free coordinates, so it has Wald standard errors too.
#' shaped <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                    n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
#'                    volume = "varying", shape = "equal", orientation = "axis")
#' parameter_inference(shaped)
#'
#' # The bootstrap resamples schools and refits inside the same family.
#' parameter_inference(shaped, method = "bootstrap", iter = 10, n_starts = 1,
#'                     seed = 1)
#'
#' # A covariate fit: membership driven by `x` and by the group's class.
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(
#'   runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#' )
#' covariate_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' covariate_fit <- multilpa(covariate_data, c("y1", "y2"), "school",
#'                           n_profiles = 2, n_group_classes = 2,
#'                           profile_covariates = "x", n_starts = 2, seed = 1)
#' parameter_inference(covariate_fit, covariate_data)
#'
#' # A latent transition model: initial and transition probabilities too.
#' moves <- lta(subset(course_engagement, student <= 40),
#'              c("browse", "lectures", "forum_read"), "student",
#'              n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
#' parameter_inference(moves)
#' @export
parameter_inference <- function(x, data = NULL, level = 0.95, step = 1e-4,
                                vcov_type = c("observed", "robust", "opg"),
                                adjust = .multilpa_p_adjust_methods,
                                method = c("wald", "bootstrap"), iter = 199L,
                                n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                                seed = NULL, boundary = c("error", "fix")) {
  UseMethod("parameter_inference")
}

#' @noRd
#' @export
parameter_inference.default <-function(x, ...) {
  # A fit of this package without a method is told so, not that it is not a
  # fit of this package (review 2026-10-27, INF-2).
  if (inherits(x, "latents_mixture_regression")) {
    stop(errorCondition(paste(
      "parameter_inference() does not take a fitted mixture regression (from",
      "mixture_regression()): its coefficients, standard errors and tests are",
      "get_results(x, \"coefficients\"), and confint() and vcov() give its",
      "intervals and covariance."), class = "latents_bad_argument", call = NULL))
  }
  if (inherits(x, .latents_model_classes)) {
    stop(errorCondition(sprintf(
      "parameter_inference() does not take %s.", .latents_describe(x)),
      class = "latents_bad_argument", call = NULL))
  }
  # Every family of fit has a method or is refused above, so anything else is
  # not a fit at all.
  .latents_require_fit(x, substitute(x), expected = "a fitted model of this package")
}

#' @noRd
#' @export
parameter_inference.multilpa <-function(x, data = NULL, level = 0.95, step = 1e-4,
                             vcov_type = c("observed", "robust", "opg"),
                             adjust = .multilpa_p_adjust_methods,
                             method = c("wald", "bootstrap"), iter = 199L,
                             n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                             seed = NULL, boundary = c("error", "fix")) {
  stopifnot(inherits(x, "multilpa"))
  .multilpa_refuse_noise(x, "parameter_inference()")
  boundary <- .latents_match_formal(boundary)
  # Classed, as three_step() and sensitivity() class the same mistakes
  # (review 2026-10-27, INF-6).
  if (!is.null(data)) .latents_check_data_frame(data, substitute(data))
  data <- .multilpa_resolve_data(x, data)
  .latents_require(
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1,
    "`step` must be a single positive number" =
      is.numeric(step) && length(step) == 1L && is.finite(step) && step > 0)
  vcov_defaulted <- missing(vcov_type)
  vcov_type <- .latents_weighted_vcov(!is.null(x$sampling_weights),
                                      .latents_match_formal(vcov_type), vcov_defaulted)
  adjust <- .latents_match_formal(adjust)
  method <- .latents_match_formal(method)
  # A seed set.seed() would not take is refused on either path, so the
  # argument means the same thing whichever method is asked for.
  .multilpa_check_seed(seed)
  if (identical(method, "bootstrap")) {
    .latents_require(
      "`iter` must be a single integer of at least two" =
        is.numeric(iter) && length(iter) == 1L && is.finite(iter) &&
        iter >= 2 && iter == floor(iter),
      "`n_starts` must be a single positive integer" =
        is.numeric(n_starts) && length(n_starts) == 1L && is.finite(n_starts) &&
        n_starts >= 1 && n_starts == floor(n_starts),
      "`max_iter` must be a single positive integer" =
        is.numeric(max_iter) && length(max_iter) == 1L && is.finite(max_iter) &&
        max_iter >= 1 && max_iter == floor(max_iter),
      "`tol` must be a single positive number" =
        is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0)
    data <- .multilpa_validate_bootstrap_data(x, data)
    .latents_local_seed(seed, after = FALSE)
    return(.multilpa_bootstrap_inference(x, data, level, as.integer(iter),
                                         as.integer(n_starts),
                                         as.integer(max_iter), tol, adjust))
  }
  vcov_type <- .multilpa_transformed_vcov(x, vcov_type, vcov_defaulted)
  .multilpa_check_regularity(x, vcov_type, boundary)
  .multilpa_check_structure_inference(x)
  if (!is.null(x$prior)) {
    stop(errorCondition(paste(
      "Wald standard errors invert the curvature of the likelihood at its",
      "maximum, and a fit made with `prior` sits at the posterior mode instead,",
      "where the likelihood's score is not zero. Use `method = \"bootstrap\"`,",
      "which refits every resample with the same prior."),
      class = "latents_unsupported_inference", call = NULL))
  }
  free <- .multilpa_free_index(x, "unconstrained")
  free_natural <- .multilpa_free_index(x, "natural")
  counted_free <- length(free)
  if (length(free) == 0L) {
    stop(errorCondition(
      "This fit holds every parameter it has, so there is nothing to report a standard error for.",
      class = "latents_no_free_parameters", call = NULL))
  }
  stopifnot("the free coordinates must match the parameters the fit counts" =
              counted_free == x$n_parameters)
  ## Response probabilities on their bound are held there, like a block the
  ## fit held, and the rest is inferred conditionally on them.
  at_bound <- if (identical(boundary, "fix")) .multilpa_response_bound(x) else
    list(natural = character(), unconstrained = character())
  free <- setdiff(free, match(at_bound$unconstrained,
                              names(.multilpa_coefficients(x, "unconstrained"))))
  free_natural <- setdiff(free_natural, match(
    at_bound$natural, names(.multilpa_coefficients(x, "natural"))))
  full_theta <- .multilpa_coefficients(x, "unconstrained")
  theta <- full_theta[free]
  prepared <- .multilpa_inference_matrix(x, data)
  observed <- prepared$x
  centered_object <- x
  centered_object$means <- sweep(x$means, 2L, prepared$centers, "-")
  # Evaluate intercepts at the predictor centres. Scaling slope coordinates
  # alone leaves a large predictor origin nearly collinear with the means.
  # The shear below maps centred coordinates back to the reported origin,
  # including every intercept/effect covariance.
  back_transform <- diag(length(full_theta))
  direct <- x$extra_data$direct
  if (!is.null(direct)) {
    centered_object$extra_data <- .latents_center_direct(
      x$extra_data, observed, x$sampling_weights[x$group_index])
    centers <- centered_object$extra_data$direct$centers
    centered_direct <- centered_object$extra_data$direct
    effects <- centered_direct$effects
    shift <- as.vector(.latents_direct_matrix(centered_direct, x$direct_effects,
                                               ncol(x$means)) %*% centers)
    centered_object$means <- sweep(centered_object$means, 2L, shift, "+")
    widths <- .multilpa_coordinate_widths(x, "unconstrained")
    last <- sum(widths[c("means", "variances", "response_probabilities", "extra")])
    positions <- last - nrow(effects) + seq_len(nrow(effects))
    mean_positions <- outer((seq_len(x$n_profiles) - 1L) * ncol(x$means),
                             effects$column, "+")
    back_transform[cbind(as.vector(mean_positions),
                           rep(positions, each = x$n_profiles))] <-
      -rep(centers[effects$covariate_column], each = x$n_profiles)
  }
  back_transform <- back_transform[free, free, drop = FALSE]
  # This is an identity plus a mean/effect shear, whose square is zero.
  forward_transform <- 2 * diag(length(free)) - back_transform
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
      .multilpa_decode(restore(parameters), x), codes,
      unname(x$sampling_weights), centered_object$extra_data)$log_likelihood
  }
  score <- function(parameters) {
    stopifnot(is.numeric(parameters))
    .multilpa_score(restore(parameters), observed, centered_object, codes)[free]
  }
  fitted_likelihood <- -objective(centered_theta)
  if (abs(fitted_likelihood - x$log_likelihood) > 1e-8 * (1 + abs(x$log_likelihood))) {
    stop(errorCondition(
      "data do not reproduce the fitted log likelihood; supply the original fitting data.",
      class = "latents_bad_inference_data", call = NULL))
  }
  parameter_scale <- c(as.vector(t(sqrt(x$variances))),
                       rep(1, length(full_theta) - length(x$means)))
  dimension <- length(.multilpa_continuous_names(x))
  if (.multilpa_uses_chart(x)) {
    chart_scale <- .multilpa_chart_encode(x)$scale
    parameter_scale[length(x$means) + seq_along(chart_scale)] <- chart_scale
  } else if (identical(x$covariance_model, "full") && dimension > 0L) {
    # Shared with the covariate path, which lacked this scaling and reported a
    # singular information for an indicator in large units. One rule, one place.
    covariance_scale <- .multilpa_covariance_coordinate_scale(
      x$variances, dimension,
      if (x$variance_model == "equal") 1L else seq_len(x$n_profiles))
    parameter_scale[length(x$means) + seq_along(covariance_scale)] <- covariance_scale
  }
  # A direct effect is in indicator units per covariate unit, so its finite
  # differences are taken on the indicator's spread over the covariate's; at
  # unit scale a covariate in units of 1e4 gave a standard error 10% too large.
  if (!is.null(direct)) {
    spread <- sqrt(colMeans(x$variances))[effects$column]
    parameter_scale[positions] <- spread / centered_direct$scales
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
  information <- .multilpa_observed_hessian(
    scaled_objective, scaled_gradient, parameter_scale, step,
    cache_key = list(path = "profile", fit = x, observed = observed,
                     codes = codes, free = free))
  scaled_hessian <- information$scaled
  hessian <- information$natural
  condition_ratio <- information$condition_ratio
  scaled_inverse <- information$inverse
  group_scores <- NULL
  group_scores_conditioned <- NULL
  scaling_correction <- NA_real_
  if (vcov_type %in% c("robust", "opg")) {
    group_scores <- .multilpa_group_scores(restore(centered_theta), observed,
                                           centered_object, codes)[, free, drop = FALSE]
    group_scores_conditioned <- sweep(group_scores, 2L, parameter_scale, "*")
    # Only the OPG estimator inverts the outer product of the scores. The
    # sandwich never does: a deficient rank there is an exact constraint the
    # estimates obey (person centring makes the profile-weighted means sum to
    # zero within every unit), whose zero variance is the right answer, and
    # the rule of more groups than parameters (.multilpa_check_regularity())
    # already guards against too few units (review 2026-10-26 round 2, INF2-2).
    if (identical(vcov_type, "opg")) {
      .inference_require_rank(group_scores_conditioned[
        .multilpa_inference_groups(x, observed_only = TRUE), , drop = FALSE],
        "The contributing group scores are rank deficient; OPG inference is unavailable.")
    }
    scaled_cross <- .multilpa_cross_product(group_scores_conditioned)
    scaling_correction <- sum(diag(scaled_inverse %*% scaled_cross)) / length(theta)
    scaled_inverse <- if (identical(vcov_type, "opg")) {
      .multilpa_opg_inverse(scaled_cross)
    } else .inference_sandwich(scaled_inverse, scaled_cross,
                                sum(.multilpa_inference_groups(x)))
    group_scores <- group_scores %*% forward_transform
  }
  centered_covariance <- scaled_inverse * tcrossprod(parameter_scale)
  covariance_unconstrained <- back_transform %*% centered_covariance %*% t(back_transform)
  # A grand mean or standard deviation the indicators were transformed with
  # is an estimate too; its variability is added through the stacked
  # estimating equations (review 2026-10-26 round 2, INF2-1). The bread is the
  # inverse the model's own covariance was formed with.
  bread <- if (identical(vcov_type, "opg")) scaled_inverse else information$inverse
  covariance_unconstrained <- .multilpa_transform_covariance(
    x, covariance_unconstrained,
    back_transform %*% (bread * tcrossprod(parameter_scale)) %*% t(back_transform),
    free, codes, robust = identical(vcov_type, "robust"))
  hessian <- t(forward_transform) %*% hessian %*% forward_transform
  dimnames(hessian) <- dimnames(covariance_unconstrained) <- list(names(theta), names(theta))
  ## The natural coordinates of a held block depend only on the held
  ## unconstrained ones, so restricting the Jacobian to the free rows and
  ## columns is the delta method for the free parameters, not an approximation
  ## of it.
  jacobian <- .multilpa_inference_jacobian(x)[free_natural, free, drop = FALSE]
  covariance <- .inference_delta_covariance(jacobian, covariance_unconstrained)
  standard_errors <- sqrt(pmax(diag(covariance), 0))
  estimates <- .multilpa_coefficients(x, "natural")[free_natural]
  critical <- .inference_critical(level, "upper")
  intervals <- cbind(estimates - critical * standard_errors, estimates + critical * standard_errors)
  colnames(intervals) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%")
  centered_gradient <- score(centered_theta)
  gradient <- stats::setNames(as.vector(centered_gradient %*% forward_transform), names(theta))
  scaled_score <- max(abs(centered_gradient * parameter_scale))
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
  standard_errors_unconstrained <- sqrt(pmax(diag(centered_covariance), 0))
  score_displacement <- max(abs(centered_gradient) * standard_errors_unconstrained)
  if (is.finite(score_displacement) && score_displacement > 0.01) {
    warning(warningCondition(sprintf(
      paste("The fitted likelihood still carries a score worth %.1f%% of a",
            "standard error; refit with a tighter `tol` before using Wald",
            "inference."), 100 * score_displacement),
      class = "latents_unconverged", call = NULL))
  }
  statistic <- unname(estimates / standard_errors)
  result <- data.frame(
    .multilpa_coefficient_labels(x, "natural")[free_natural, , drop = FALSE],
    estimate = unname(estimates), standard_error = unname(standard_errors),
    statistic = statistic, p_value = 2 * stats::pnorm(-abs(statistic)),
    conf_low = unname(intervals[, 1L]), conf_high = unname(intervals[, 2L]),
    row.names = NULL, stringsAsFactors = FALSE)
  continuous <- .multilpa_continuous_names(x)
  bounds <- .multilpa_wald_bounds(
    result$estimate, result$standard_error,
    .multilpa_interval_kind(result$parameter, result$term, continuous), critical)
  result$conf_low <- bounds$low
  result$conf_high <- bounds$high
  ## A Wald test of a variance against zero is not a question worth asking; the
  ## interval still is.
  bounded <- result$parameter %in% c("variance", "probability", "response",
                                      "count_mean", "count_dispersion") |
    (result$parameter == "covariance" &
       result$term %in% .latents_term(continuous, continuous))
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
    method = "wald", vcov_type = vcov_type, group_scores = group_scores,
    hessian_conditioned = scaled_hessian,
    group_scores_conditioned = group_scores_conditioned,
    scaling_correction = scaling_correction,
    fixed = x$fixed %||% character(),
    fixed_at_bound = at_bound$natural))
  result
}

#' Response probabilities sitting on their lower bound
#'
#' A probability held at `min_probability` has a logit on the edge of the
#' space, where the likelihood is flat and the information singular. A row's
#' logits are measured against its last category, `log(p_c / p_last)`.
#'
#' - A category on the floor is held, with its logit.
#' - When the last category is on the floor, every logit is near `log(1 / floor)`
#'   and the likelihood is flat in one direction only: the common shift that
#'   moves `p_last` alone. The differences between the logits, which are the
#'   ratios of the interior probabilities, stay identified. So one logit is
#'   held, that of the row's largest category, which fixes `log(p_last / p_max)`
#'   and leaves the others free. This is the bound category measured against an
#'   off-bound one. Holding the whole row instead withheld the errors of
#'   interior probabilities (0.943 and 0.057 on gss82) and made every other
#'   error conditional on them and too small (class shares 0.0396 against
#'   0.0521; review 2026-10-26, INF-6).
#' - A category whose every other category is on the floor sits on the upper
#'   bound, one minus the floors, and is held too. In a binary item that is the
#'   complement of a floored probability, which otherwise had an error of
#'   exactly zero when it was the last category (INF-7).
#'
#' @param object A fitted model with `response_probabilities`.
#' @return A list of two character vectors of parameter names, `natural`
#'   (probabilities) and `unconstrained` (logits), in the
#'   `level.parameter.outcome.term` spelling every fitted class uses.
#' @noRd
.multilpa_response_bound <- function(object) {
  blocks <- object$response_probabilities %||% list()
  floor <- .multilpa_probability_floor(object)
  pieces <- lapply(names(blocks), function(indicator) {
    block <- blocks[[indicator]]
    categories <- colnames(block)
    at <- block <= floor
    others_at <- rowSums(at) - at
    held <- at | others_at == ncol(block) - 1L
    largest <- col(block) == max.col(block, ties.method = "first")
    held_logit <- (at | (at[, ncol(block)] & largest))[, -ncol(block), drop = FALSE]
    ## Spelled as coef() spells the same parameter, by the same helper.
    spell <- function(positions, kind) {
      sprintf("measurement.%s.profile_%d.%s", kind, positions[, 1L],
              .latents_term(rep(indicator, nrow(positions)),
                            categories[positions[, 2L]]))
    }
    list(natural = spell(which(held, arr.ind = TRUE), "response"),
         unconstrained = spell(which(held_logit, arr.ind = TRUE), "response_logit"))
  })
  list(natural = unlist(lapply(pieces, `[[`, "natural")) %||% character(),
       unconstrained = unlist(lapply(pieces, `[[`, "unconstrained")) %||% character())
}

#' The probability below which an estimate is on its boundary
#'
#' EM approaches a zero probability slowly, so a fit rarely lands exactly on
#' `min_probability`; it stops a little above it, at a logit of -15 or -20
#' where the likelihood is flat. Treating only the exact floor as the bound
#' left such an estimate free, and the observed information was then singular.
#' At or below 1e-6 a probability's expected count is far below one
#' observation at any realistic sample size, so the data cannot tell it from
#' zero and it is treated as on its bound.
#' @param object A fitted model.
#' @return A single probability.
#' @noRd
.multilpa_probability_floor <- function(object) {
  max((object$min_probability %||% 1e-10) * (1 + 1e-7), 1e-6)
}

#' @rdname latents-model-methods
#' @export
#' @importFrom stats coef
coef.multilpa <- function(object, scale = c("natural", "unconstrained"), ...) {
  stopifnot(inherits(object, "multilpa"))
  .multilpa_coefficients(object, .latents_match_formal(scale))
}

#' @rdname latents-model-methods
#' @export
#' @importFrom stats vcov
vcov.multilpa <- function(object, data = NULL, scale = c("natural", "unconstrained"), ...) {
  stopifnot(inherits(object, "multilpa"))
  .multilpa_reject_inference_arguments(list(...), "vcov()")
  scale <- .latents_match_formal(scale)
  # `data` is no longer required: a fit carries the columns it was built from,
  # and `parameter_inference()` falls back to them.
  information <- parameter_inference(object, data, ...)
  if (scale == "natural") return(attr(information, "covariance"))
  covariance <- attr(information, "covariance_unconstrained")
  # The bootstrap refits each resample in its own coordinates and reads the
  # natural parameters only, so it has no estimation-scale covariance. Saying
  # so is the answer; a NULL would flow silently into the caller's next step.
  if (is.null(covariance)) {
    stop(errorCondition(paste(
      "The bootstrap has no covariance on the estimation scale: it reads the",
      "natural parameters of each refit. Use `scale = \"natural\"`, or the",
      "Wald path (`method = \"wald\"`) for the estimation scale."),
      class = "latents_unsupported_inference", call = NULL))
  }
  covariance
}

#' @rdname latents-model-methods
#' @export
#' @importFrom stats confint
confint.multilpa <- function(object, parm, level = 0.95, data = NULL, ...) {
  stopifnot(inherits(object, "multilpa"))
  .latents_require(
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1)
  .multilpa_reject_inference_arguments(list(...), "confint()")
  estimates <- coef.multilpa(object)
  ## One call, not one through `vcov()`: the bootstrap interval is the
  ## percentile interval this table already carries, and reaching it through
  ## the covariance would rebuild it as `estimate +/- z * se`, a different
  ## interval from the one `parameter_inference()` reports for the same fit.
  information <- parameter_inference(object, data = data, level = level, ...)
  covariance <- attr(information, "covariance")
  ## A held measurement coefficient is reported by coef() -- it is part of the
  ## model -- but it was not estimated here and has no interval. The default is
  ## therefore what the covariance actually covers.
  estimated <- rownames(covariance)
  if (missing(parm)) parm <- estimated
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm)) || any(parm != floor(parm)) ||
        any(parm < 1) || any(parm > length(estimates))) {
      stop(errorCondition("Invalid coefficient indices in parm.",
                          class = "latents_bad_argument", call = NULL))
    }
    parm <- names(estimates)[parm]
  }
  if (!is.character(parm) || anyNA(parm) || !all(parm %in% names(estimates))) {
    stop(errorCondition("parm must identify existing coefficient names or indices.",
                        class = "latents_bad_argument", call = NULL))
  }
  held <- setdiff(parm, estimated)
  if (length(held) > 0L) {
    stop(errorCondition(sprintf(
      "%s %s held fixed by this fit, so %s carr%s no confidence interval.",
      paste(sprintf("`%s`", held), collapse = ", "),
      if (length(held) == 1L) "was" else "were",
      if (length(held) == 1L) "it" else "they",
      if (length(held) == 1L) "ies" else "y"),
      class = "latents_held_parameter", call = NULL))
  }
  # The intervals parameter_inference() reports, whichever path made them:
  # percentile for the bootstrap, and for Wald ones formed where the parameter
  # is unbounded, so the two verbs can never disagree.
  intervals <- cbind(information$conf_low, information$conf_high)
  rownames(intervals) <- rownames(covariance)
  intervals <- intervals[parm, , drop = FALSE]
  colnames(intervals) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%")
  intervals
}

#' Refuse an argument that `parameter_inference()` does not take
#'
#' `vcov()` and `confint()` pass `...` on to [parameter_inference()], where a
#' misnamed argument (`type = "robust"`, the name the transition and
#' regression fits use) stopped with R's unclassed "unused argument" error.
#' It is refused here with `latents_bad_argument` instead, naming it.
#' @param extra `list(...)` as the method received it.
#' @param verb The method, for the message.
#' @return `NULL`, invisibly, or `latents_bad_argument`.
#' @noRd
.multilpa_reject_inference_arguments <- function(extra, verb) {
  taken <- setdiff(names(formals(parameter_inference)), c("x", "data"))
  supplied <- names(extra) %||% rep("", length(extra))
  # Partial names are matched as parameter_inference() matches them.
  known <- nzchar(supplied) & !is.na(pmatch(supplied, taken, duplicates.ok = TRUE))
  .multilpa_reject_extra_arguments(
    extra[!known], verb,
    paste("Its other arguments go to parameter_inference(); choose the",
          "covariance with `vcov_type =`."))
}

#' Check that a fit admits Wald inference at all
#'
#' Boundary solutions, unconverged fits and zero mixing probabilities each put
#' the estimate where the asymptotic normal approximation does not hold, so
#' they are refused rather than reported with a number.
#'
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_check_regularity <- function(object, vcov_type, boundary = "error") {
  n_groups <- sum(.multilpa_inference_groups(object, observed_only = TRUE))
  if (vcov_type %in% c("robust", "opg") && n_groups <= object$n_parameters) {
    stop(errorCondition(sprintf(
      "Robust and OPG inference need more groups than parameters; this fit has %d groups and %d parameters.",
      n_groups, object$n_parameters),
      class = "latents_too_few_groups", call = NULL))
  }
  if (!isTRUE(object$converged)) {
    stop(errorCondition("Inference requires a converged fit.",
                        class = "latents_no_converge", call = NULL))
  }
  ## A bound-active value that this fit held fixed is a constant of the
  ## conditional likelihood, not an estimate sitting on its boundary, so it does
  ## not disqualify the parameters that were actually estimated.
  held <- object$fixed %||% character()
  if (isTRUE(object$boundary) && !("variances" %in% held)) {
    stop(errorCondition("Wald inference is unavailable for a bound-active fit.",
                        class = "latents_boundary_fit", call = NULL))
  }
  ## `boundary = "fix"` asks for the bound-active response probabilities to be
  ## held at their bound and the rest estimated conditionally on them.
  response <- if ("response_probabilities" %in% held ||
                  identical(boundary, "fix")) NULL else
    unlist(object$response_probabilities, use.names = FALSE)
  if (length(response) > 0L &&
      any(response <= .multilpa_probability_floor(object))) {
    stop(errorCondition(
      paste("Wald inference is unavailable for bound-active categorical response",
            "probabilities. `boundary = \"fix\"` holds them at their bound and",
            "reports the other parameters conditionally on them."),
      class = "latents_boundary_fit", call = NULL))
  }
  if (length(object$count_means) > 0L && any(object$count_means <= 1e-8)) {
    stop(errorCondition(
      "Wald inference is unavailable when a count mean is estimated at zero.",
      class = "latents_boundary_fit", call = NULL))
  }
  if (length(object$count_dispersion) > 0L &&
      any(object$count_dispersion <= .latents_min_dispersion * (1 + 1e-6))) {
    stop(errorCondition(paste(
      "A negative-binomial dispersion is at its floor (the Poisson limit), where",
      "Wald inference does not apply; refit with `count_model = \"poisson\"`."),
      class = "latents_boundary_fit", call = NULL))
  }
  if (any(object$profile_probabilities <= 0) || any(object$group_probabilities <= 0)) {
    stop(errorCondition(
      "Wald inference requires strictly positive mixing probabilities.",
      class = "latents_boundary_fit", call = NULL))
  }
  invisible(NULL)
}

#' Give the caller's frame a single-level fit's internal row numbering
#'
#' A single-level fit numbered its rows itself; the caller's frame has no such
#' column, and its row order is what the numbering was. Frames that already
#' carry the column, or have a different number of rows, are returned as given
#' for the caller's own checks to judge.
#'
#' @param object A fitted model of this package.
#' @param data The caller's frame.
#' @return `data`, with the fit's identifier column added when it was absent.
#' @noRd
.multilpa_number_single_level <- function(object, data) {
  if (isTRUE(object$single_level) && !object$id %in% names(data) &&
      nrow(data) == object$n_observations) {
    data[[object$id]] <- object$group_values[seq_len(nrow(data))]
  }
  data
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
  data <- .multilpa_number_single_level(object, data)
  if (nrow(data) != object$n_observations ||
      !all(c(object$vars, object$id) %in% names(data)) ||
      anyDuplicated(names(data)) ||
      !all(vapply(data[, .multilpa_continuous_names(object), drop = FALSE],
                  is.numeric, logical(1)))) {
    stop(errorCondition(
      "data must contain the original numeric indicators and group column.",
      class = "latents_bad_inference_data", call = NULL))
  }
  group_index <- match(data[[object$id]], object$group_values)
  if (!identical(group_index, object$group_index)) {
    stop(errorCondition(
      "data must retain the original group identifiers and row order.",
      class = "latents_bad_inference_data", call = NULL))
  }
  ## Only the Gaussian indicators form the numeric matrix; the categorical ones
  ## are re-encoded from `data` rather than taken from the fit, so that supplying
  ## the wrong categorical columns is caught rather than silently accepted.
  continuous <- .multilpa_continuous_names(object)
  ## as.matrix() of a zero-column frame is logical, and the expectation demands
  ## a numeric matrix, so an all-categorical fit needs the mode forced.
  x <- if (length(continuous) == 0L) matrix(numeric(0), nrow(data), 0L) else
    .multilpa_center_like(object, as.matrix(data[, continuous, drop = FALSE]))
  # A supplied column is the fitted indicator when it passes the identity
  # rule every verb uses (`.multilpa_indicator_reproduced()`): exactly, as
  # the fit's own transform of it, or as a write.csv() round trip of the
  # original or of its reconstruction. The stored fitting values are then
  # used, so the numbers do not depend on which copy was supplied; a changed
  # value is refused (review 2026-10-26 round 2, INF-3).
  if (!is.null(object$indicator_data) && length(continuous) > 0L) {
    reconstructed <- .multilpa_uncentered_indicators(object)
    reproduced <- vapply(continuous, function(name) {
      .multilpa_indicator_reproduced(object, name, data[[name]], reconstructed)
    }, logical(1))
    if (!all(reproduced)) {
      stop(errorCondition(sprintf(paste(
        "data must reproduce the original indicator data, including row order;",
        "%s differ%s from the fit."),
        paste(sprintf("`%s`", continuous[!reproduced]), collapse = ", "),
        if (sum(!reproduced) == 1L) "s" else ""),
        class = "latents_bad_inference_data", call = NULL))
    }
    x <- object$indicator_data
  }
  if (any(is.infinite(x)) || any(is.nan(x)) ||
      (anyNA(x) && !identical(object$missing, "fiml"))) {
    stop(errorCondition(
      "data contain unsupported missing or non-finite indicators.",
      class = "latents_bad_inference_data", call = NULL))
  }
  if (!is.null(object$indicator_data) && ncol(x) > 0L &&
      !isTRUE(all.equal(x, object$indicator_data, tolerance = 0))) {
    stop(errorCondition(
      "data must reproduce the original indicator data, including row order and names.",
      class = "latents_bad_inference_data", call = NULL))
  }
  codes <- NULL
  # Rows of weight zero do not enter the likelihood, so their categorical and
  # ordinal values are re-encoded as the fit encoded them (CORE-2).
  contributing <- if (is.null(object$sampling_weights)) NULL else
    unname(object$sampling_weights)[object$group_index] > 0
  if (length(object$categorical %||% character()) > 0L) {
    encoded <- .multilpa_encode_categorical(data[, object$categorical, drop = FALSE],
                                            contributing = contributing)
    codes <- encoded$codes
    if (!identical(unname(codes), unname(object$categorical_data)) ||
        !identical(encoded$levels, object$categorical_levels)) {
      stop(errorCondition(
        "data must reproduce the original categorical indicators, including their categories and row order.",
        class = "latents_bad_inference_data", call = NULL))
    }
  }
  if (!is.null(object$extra_data)) {
    extra <- .latents_prepare_extra(data, object$ordinal, object$count,
                                    object$missing %||% "error",
                                    contributing = contributing)
    # A carrier can hold direct effects alone, with no indicator to compare.
    if (length(c(object$ordinal, object$count)) > 0L &&
        (!identical(unname(extra$ordinal), unname(object$extra_data$ordinal)) ||
           !identical(extra$ordinal_levels, object$extra_data$ordinal_levels) ||
           !identical(unname(extra$count), unname(object$extra_data$count)))) {
      stop(errorCondition(
        "data must reproduce the original ordinal and count indicators.",
        class = "latents_bad_inference_data", call = NULL))
    }
    direct <- object$extra_data$direct
    stored_covariates <- direct$covariates
    # The covariates are compared as the row-order check compares them: exactly
    # or as their write.csv() serialization (INF-3).
    serialized <- !is.null(direct) &&
      all(colnames(stored_covariates) %in% names(data)) &&
      all(vapply(colnames(stored_covariates), function(name) {
        .multilpa_alignment_values(data[[name]], stored_covariates[, name])
      }, logical(1)))
    if (!is.null(direct) && !serialized && !identical(
      unname(.latents_direct_like(direct, data)$covariates), unname(direct$covariates))) {
      stop(errorCondition(
        "data must reproduce the original direct-effect covariates and row order.",
        class = "latents_bad_inference_data", call = NULL))
    }
  }
  centers <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  list(x = if (ncol(x) > 0L) sweep(x, 2L, centers, "-") else x,
       centers = centers, codes = codes)
}

#' Validate the original observations before any bootstrap resampling
#'
#' A bootstrap reports uncertainty for the fit's estimates, so its data must
#' reproduce the same observations, grouping, ordering and normalized sampling
#' weights. A likelihood check also protects older fits without stored rows.
#' @param object A fitted `multilpa` model.
#' @param data The supplied or reconstructed original data frame.
#' @return The checked data, including a single-level fit's internal row ID.
#' @noRd
.multilpa_validate_bootstrap_data <- function(object, data) {
  stopifnot(inherits(object, "multilpa"), is.data.frame(data))
  data <- .multilpa_number_single_level(object, data)
  prepared <- .multilpa_inference_matrix(object, data)
  .multilpa_check_alignment(object, data)
  if (!is.null(object$weights)) {
    supplied_weights <- tryCatch(
      .latents_sampling_weights(data, object$weights, object$group_index,
                               object$n_groups),
      error = function(condition) {
        stop(errorCondition(paste("data must reproduce the fitted sampling weights:",
                                  conditionMessage(condition)),
                            class = "latents_bad_inference_data", call = NULL))
      })
    if (!isTRUE(all.equal(unname(supplied_weights),
                          unname(object$sampling_weights), tolerance = 1e-12))) {
      stop(errorCondition("data must reproduce the fitted sampling weights.",
                          class = "latents_bad_inference_data", call = NULL))
    }
  }
  parameters <- .multilpa_decode(.multilpa_coefficients(object, "unconstrained"),
                                 object)
  if (length(prepared$centers)) {
    parameters$means <- sweep(parameters$means, 2L, prepared$centers, "-")
  }
  likelihood <- .multilpa_expectation(
    prepared$x, object$group_index, parameters, prepared$codes,
    unname(object$sampling_weights), object$extra_data)$log_likelihood
  if (!is.finite(likelihood) ||
      abs(likelihood - object$log_likelihood) >
        1e-8 * (1 + abs(object$log_likelihood))) {
    stop(errorCondition("data must reproduce the fitted log likelihood.",
                        class = "latents_bad_inference_data", call = NULL))
  }
  data
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
#' Labels inside the brackets are escaped (`.latents_parameter_component()`);
#' the `term` column gives them back as their own text, except that a term
#' pairing two labels escapes `:` and `%` within each (`.latents_term()`).
#' Coefficient names are serialised from these columns, so `coef()`, `vcov()`,
#' `confint()` and the inference table all spell a label one way.
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
  ## The components are escaped for this grammar; a term shows the labels'
  ## own text, escaping only what a paired term needs (`.latents_term()`).
  data.frame(
    level = ifelse(measurement, "measurement",
                   ifelse(within_profile, "profile", "group")),
    outcome = ifelse(measurement, profile_label(first),
                     ifelse(within_profile,
                            profile_label(second), class_label(first))),
    term = ifelse(prefix %in% .multilpa_paired_kinds,
                  .latents_term(.latents_parameter_decode(second),
                                .latents_parameter_decode(third)),
                  ifelse(measurement, .latents_parameter_decode(second),
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
                                 "log_volume", "log_shape", "rotation",
                                 "shape_cholesky", "log_shape_cholesky",
                                 "response", "response_logit",
                                 "ordinal_intercept", "ordinal_location",
                                 "count_mean", "log_count_mean",
                                 "count_dispersion", "log_count_dispersion",
                                 "direct_effect")

#' Parameter kinds whose term names a pair, written `a:b`
#' @noRd
.multilpa_paired_kinds <- c("covariance", "cholesky", "log_cholesky",
                            "rotation", "shape_cholesky", "log_shape_cholesky",
                            "response", "response_logit", "ordinal_intercept",
                            "direct_effect")

#' Multiplicity corrections an inference verb accepts
#'
#' The same methods [stats::p.adjust()] knows, with `"none"` first so that
#' `match.arg()` makes the uncorrected column the default. Applying a correction
#' silently would change what every p-value in the table means.
#'
#' @noRd
.multilpa_p_adjust_methods <- c("none",
                                setdiff(stats::p.adjust.methods, "none"))

