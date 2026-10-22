# Direct effects of covariates on continuous indicators.
#
# A direct effect moves an indicator's mean with a covariate within every
# class, by the same amount in each: y_ij = mu_kj + b_j' z_i + e_ij, e_ij with
# the class's variance. Latent GOLD fits it with `y <- 1 + Cluster + z;`. The
# covariate does not change who is in which class (that is
# `profile_covariates`); it removes from the indicator a shift that would
# otherwise be read as a difference between classes (differential item
# functioning; Bandeen-Roche et al. 1997).
#
# The effects ride with the ordinal and count indicators in the fit's `extra`
# carrier, so every evaluation of the likelihood sees the covariates: the
# Gaussian block measures x - Z B' instead of x (.latents_measurement_log_density),
# the effects are the last coordinates of the extra block, their natural and
# unconstrained values coincide, and the M-step updates them once per
# iteration given the new means and variances, an ECM step (Meng and Rubin
# 1993) that keeps the likelihood increasing. The covariates enter with their
# values as given, so `means` are the class means at a covariate value of 0.

#' Parse direct-effect specifications
#'
#' @param direct_effects A character vector of `"indicator ~ covariate"`
#'   terms; several covariates are joined with `+`.
#' @return A data frame with `indicator` and `covariate`, one row per effect,
#'   in the order given.
#' @noRd
.latents_parse_direct_effects <- function(direct_effects) {
  if (!is.character(direct_effects) || anyNA(direct_effects)) {
    stop(errorCondition("`direct_effects` must be a character vector without missing values.",
                        class = "latents_bad_argument", call = NULL))
  }
  # strsplit() drops a trailing empty field ("y ~ z +" would read as "y ~ z"),
  # so a space is appended to keep it and have it refused as empty.
  split <- function(text, separator) strsplit(paste0(text, " "), separator, fixed = TRUE)
  pieces <- split(direct_effects, "~")
  malformed <- lengths(pieces) != 2L
  sides <- lapply(pieces[!malformed], trimws)
  indicators <- vapply(sides, `[[`, character(1), 1L)
  covariates <- lapply(sides, function(side) trimws(split(side[[2L]], "+")[[1L]]))
  malformed[!malformed] <- !nzchar(indicators) |
    vapply(covariates, function(names) length(names) == 0L || !all(nzchar(names)), logical(1))
  if (any(malformed)) {
    stop(errorCondition(sprintf(paste(
      "`direct_effects` terms are written \"indicator ~ covariate\", several",
      "covariates joined with `+`; %s is not."),
      paste(sprintf("\"%s\"", direct_effects[malformed]), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  effects <- data.frame(indicator = rep(indicators, lengths(covariates)),
                        covariate = unlist(covariates), stringsAsFactors = FALSE)
  repeated <- duplicated(effects)
  if (any(repeated)) {
    stop(errorCondition(sprintf(
      "`direct_effects` names %s more than once.",
      paste(sprintf("\"%s ~ %s\"", effects$indicator[repeated], effects$covariate[repeated]),
            collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  effects
}

#' Refuse the settings direct effects are not fitted with
#' @param conflicts A named logical vector, `TRUE` where a setting is in use.
#' @return `NULL`, invisibly.
#' @noRd
.latents_check_direct_conflicts <- function(conflicts) {
  if (any(conflicts)) {
    stop(errorCondition(sprintf(paste(
      "Direct effects are fitted in a profile model with a diagonal covariance",
      "(EEI or VVI); they are not combined with %s."),
      paste(names(conflicts)[conflicts], collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  invisible(NULL)
}

#' Put prepared direct effects into the extra carrier
#' @param extra The carrier of the ordinal and count indicators, or `NULL`.
#' @param direct From `.latents_prepare_direct()`.
#' @return The carrier with `direct`; a carrier without ordinal or count
#'   indicators when there were none.
#' @noRd
.latents_add_direct <- function(extra, direct) {
  extra <- extra %||% list(ordinal = NULL, ordinal_levels = list(),
                           ordinal_categories = integer(), count = NULL,
                           count_model = "poisson", count_dispersion = "varying")
  extra$direct <- direct
  extra
}

#' Centre the covariates of direct effects for fitting
#'
#' The model y = m + B z is the same model as y = (m + B c) + B (z - c), so
#' EM can run on centred covariates, where the intercepts and the effects are
#' nearly uncorrelated, and `.latents_uncenter_direct()` moves the means back.
#'
#' @param extra The extra carrier, or `NULL`.
#' @param x Continuous indicators, used to select the observations of each
#'   effect's target; `NULL` uses every row.
#' @param row_weights Optional sampling weights in observation order.
#' @return The carrier with one centred internal covariate column per effect,
#'   its target-specific `centers` and `scales`. The original carrier is kept
#'   on the reported fit, retaining the caller's covariate names and columns.
#' @noRd
.latents_center_direct <- function(extra, x = NULL, row_weights = NULL) {
  stopifnot(is.null(extra) || is.list(extra), is.null(x) || is.matrix(x),
            is.null(row_weights) || is.numeric(row_weights))
  if (is.null(extra$direct)) return(extra)
  effects <- extra$direct$effects
  values <- extra$direct$covariates[, effects$covariate_column, drop = FALSE]
  weights <- row_weights %||% rep(1, nrow(values))
  stopifnot(length(weights) == nrow(values),
            is.null(x) || nrow(x) == nrow(values))
  # A predictor value on a missing target or a zero-weight row supplies no
  # information about this effect. Including it in the centre can make a
  # regular intercept/slope pair arbitrarily collinear. Each target gets its
  # own chart, even when its predictor is shared with another target.
  observed <- if (is.null(x)) matrix(TRUE, nrow(values), ncol(values)) else
    !is.na(x[, effects$column, drop = FALSE])
  observed <- observed & weights > 0
  centers <- vapply(seq_len(ncol(values)), function(column) {
    rows <- observed[, column]
    stats::weighted.mean(values[rows, column], weights[rows])
  }, numeric(1))
  scales <- vapply(seq_len(ncol(values)), function(column) {
    stats::sd(values[observed[, column], column])
  }, numeric(1))
  extra$direct$covariates <- sweep(values, 2L, centers, "-")
  extra$direct$effects$covariate_column <- seq_len(nrow(effects))
  extra$direct$centers <- centers
  extra$direct$scales <- scales
  extra
}

#' Express means fitted on centred covariates at a covariate value of zero
#' @param parameters The fitted parameters, means at the covariates' centres.
#' @param extra The carrier the fit was made with, from `.latents_center_direct()`.
#' @return `parameters` with `means` at zero covariate values.
#' @noRd
.latents_uncenter_direct <- function(parameters, extra) {
  centers <- extra$direct$centers
  if (is.null(centers)) return(parameters)
  shift <- .latents_direct_matrix(extra$direct, parameters$direct_effects,
                                  ncol(parameters$means)) %*% centers
  parameters$means <- sweep(parameters$means, 2L, as.vector(shift), "-")
  parameters
}

#' A fit's direct effects as a set, for comparing two fits
#' @param fit A fitted model.
#' @return The sorted `"indicator:covariate"` names; empty without effects.
#' @noRd
.latents_direct_key <- function(fit) {
  effects <- fit$extra_data$direct$effects
  if (is.null(effects)) return(character())
  sort(paste(effects$indicator, effects$covariate, sep = ":"))
}

#' The `direct_effects` terms that give a table of effects back
#' @param effects From `.latents_parse_direct_effects()`.
#' @return One `"indicator ~ covariate + ..."` term per indicator.
#' @noRd
.latents_direct_terms <- function(effects) {
  indicators <- unique(effects$indicator)
  vapply(indicators, function(indicator) {
    sprintf("%s ~ %s", indicator,
            paste(effects$covariate[effects$indicator == indicator], collapse = " + "))
  }, character(1), USE.NAMES = FALSE)
}

#' Validate direct effects against the data and prepare the covariates
#'
#' @param data The fitting data.
#' @param effects From `.latents_parse_direct_effects()`.
#' @param vars,continuous The fit's indicators and its continuous ones, in
#'   the order of the continuous matrix.
#' @return A list with `covariates` (rows by covariates, numeric) and
#'   `effects`, which adds `column` (the indicator's column of the continuous
#'   matrix) and `covariate_column`.
#' @noRd
.latents_prepare_direct <- function(data, effects, vars, continuous) {
  not_continuous <- setdiff(effects$indicator, continuous)
  if (length(not_continuous) > 0L) {
    stop(errorCondition(sprintf(paste(
      "A direct effect moves the mean of a continuous indicator; %s %s not a",
      "continuous indicator among `vars`."),
      paste(sprintf("`%s`", not_continuous), collapse = ", "),
      if (length(not_continuous) > 1L) "are" else "is"),
      class = "latents_bad_argument", call = NULL))
  }
  covariates <- unique(effects$covariate)
  indicators_as_covariates <- intersect(covariates, vars)
  if (length(indicators_as_covariates) > 0L) {
    stop(errorCondition(sprintf(
      "%s %s an indicator, so it cannot also be a covariate of a direct effect.",
      paste(sprintf("`%s`", indicators_as_covariates), collapse = ", "),
      if (length(indicators_as_covariates) > 1L) "are" else "is"),
      class = "latents_bad_argument", call = NULL))
  }
  absent <- setdiff(covariates, names(data))
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf("`data` has no column %s for a direct effect.",
                                paste(sprintf("`%s`", absent), collapse = ", ")),
                        class = "latents_bad_data", call = NULL))
  }
  values <- .latents_direct_covariates(data, covariates)
  constant <- covariates[vapply(covariates, function(name) {
    length(unique(values[, name])) < 2L
  }, logical(1))]
  if (length(constant) > 0L) {
    stop(errorCondition(sprintf(paste(
      "Covariate %s takes one value only, so its effect cannot be told apart",
      "from the class means."), paste(sprintf("`%s`", constant), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  effects$column <- match(effects$indicator, continuous)
  effects$covariate_column <- match(effects$covariate, covariates)
  list(covariates = values, effects = effects)
}

#' The covariate columns of a direct effect, checked
#' @param data A data frame holding the covariates.
#' @param covariates Their names.
#' @return A rows-by-covariates numeric matrix.
#' @noRd
.latents_direct_covariates <- function(data, covariates) {
  absent <- setdiff(covariates, names(data))
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf(
      "`data` lacks the direct-effect covariate%s %s.",
      if (length(absent) > 1L) "s" else "",
      paste(sprintf("`%s`", absent), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  unusable <- covariates[!vapply(data[covariates], function(value) {
    is.numeric(value) && is.null(dim(value)) && all(is.finite(value))
  }, logical(1))]
  if (length(unusable) > 0L) {
    stop(errorCondition(sprintf(paste(
      "A direct-effect covariate must be numeric and observed in every row;",
      "%s %s not. Code a categorical covariate as indicator columns."),
      paste(sprintf("`%s`", unusable), collapse = ", "),
      if (length(unusable) > 1L) "are" else "is"),
      class = "latents_bad_data", call = NULL))
  }
  matrix(as.numeric(as.matrix(data[covariates])), nrow(data), length(covariates),
         dimnames = list(NULL, covariates))
}

#' The direct effects as an indicators-by-covariates coefficient matrix
#' @param direct The `direct` entry of the extra carrier.
#' @param values The effects, one per row of `direct$effects`.
#' @param n_indicators Columns of the continuous matrix.
#' @return An `n_indicators` by covariates matrix, zero where no effect is fitted.
#' @noRd
.latents_direct_matrix <- function(direct, values, n_indicators) {
  coefficients <- matrix(0, n_indicators, ncol(direct$covariates))
  coefficients[cbind(direct$effects$column, direct$effects$covariate_column)] <- values
  coefficients
}

#' Remove the direct effects from the continuous indicators
#'
#' @param x The continuous matrix, `NA` where unobserved.
#' @param extra The extra carrier, or `NULL`.
#' @param parameters A parameter list, possibly with `direct_effects`.
#' @return `x - Z B'`; `x` itself when there are no effects.
#' @noRd
.latents_direct_shift <- function(x, extra, parameters) {
  direct <- extra$direct
  if (is.null(direct) || is.null(parameters$direct_effects)) return(x)
  x - tcrossprod(direct$covariates,
                 .latents_direct_matrix(direct, parameters$direct_effects, ncol(x)))
}

#' Update the direct effects given the class means and variances
#'
#' With the means and the diagonal variances held, the expected complete-data
#' log likelihood of indicator j is a weighted least squares in its effects:
#' sum_i sum_k r_ik (x_ij - mu_kj - b_j' z_i)^2 / s2_kj, whose minimum solves
#' (sum_i w_ij z_i z_i') b_j = sum_i z_i sum_k r_ik (x_ij - mu_kj) / s2_kj with
#' w_ij = sum_k r_ik / s2_kj. Rows where the indicator is unobserved carry no
#' information about it under a diagonal covariance.
#'
#' @param x The continuous matrix, unshifted.
#' @param extra The extra carrier.
#' @param posteriors Rows by profiles.
#' @param means,variances Profiles by indicators.
#' @return The effects, one per row of `extra$direct$effects`.
#' @noRd
.latents_direct_maximize <- function(x, extra, posteriors, means, variances) {
  direct <- extra$direct
  effects <- direct$effects
  updated <- numeric(nrow(effects))
  indicators <- unique(effects$column)
  solved <- lapply(indicators, function(column) {
    rows <- which(effects$column == column)
    observed <- !is.na(x[, column])
    precision <- sweep(posteriors[observed, , drop = FALSE], 2L, variances[, column], "/")
    weight <- rowSums(precision)
    residual <- rowSums(precision * outer(x[observed, column], means[, column], "-"))
    design <- direct$covariates[observed, effects$covariate_column[rows], drop = FALSE]
    as.vector(qr.solve(crossprod(design, design * weight), crossprod(design, residual)))
  })
  updated[unlist(lapply(indicators, function(column) which(effects$column == column)))] <-
    unlist(solved)
  updated
}

#' Scores of the direct effects
#'
#' The score of b_jc is sum_i z_ic sum_k r_ik (x~_ij - mu_kj) / s2_kj, where
#' x~ is the indicator with the effects removed: the mean score of indicator j
#' weighted by the covariate.
#'
#' @param x The continuous matrix, unshifted.
#' @param extra The extra carrier.
#' @param posteriors Rows by profiles (weighted, under sampling weights).
#' @param parameters The point, with `means`, `variances` and `direct_effects`.
#' @param group_index Group of each row, or `NULL` for the total.
#' @return Groups (or one row) by effects.
#' @noRd
.latents_direct_scores <- function(x, extra, posteriors, parameters, group_index = NULL) {
  direct <- extra$direct
  effects <- direct$effects
  shifted <- .latents_direct_shift(x, extra, parameters)
  contribution <- vapply(seq_len(nrow(effects)), function(row) {
    column <- effects$column[[row]]
    observed <- !is.na(shifted[, column])
    residual <- ifelse(observed, shifted[, column], 0)
    precision <- sweep(posteriors, 2L, parameters$variances[, column], "/")
    rowSums(precision * outer(residual, parameters$means[, column], "-")) * observed *
      direct$covariates[, effects$covariate_column[[row]]]
  }, numeric(nrow(x)))
  contribution <- matrix(contribution, nrow(x), nrow(effects))
  if (is.null(group_index)) {
    matrix(colSums(contribution), 1L, ncol(contribution))
  } else rowsum(contribution, group_index, reorder = FALSE)
}

#' Rebuild a fit's direct-effect covariates from a data frame
#' @param direct The fit's `direct` entry.
#' @param data A data frame with the covariate columns.
#' @return The `direct` entry with the covariates read from `data`.
#' @noRd
.latents_direct_like <- function(direct, data) {
  direct$covariates <- .latents_direct_covariates(data, colnames(direct$covariates))
  direct
}

#' The direct effects of a fit, with their standard errors
#'
#' @param x A fit with direct effects.
#' @param data The fitting data, or `NULL` for the columns the fit carries.
#' @return One row per effect: `indicator`, `covariate`, `estimate`,
#'   `standard_error`, `statistic`, `p_value`, `conf_low` and `conf_high`.
#' @noRd
.latents_direct_frame <- function(x, data = NULL) {
  effects <- x$extra_data$direct$effects
  if (is.null(effects)) {
    stop(errorCondition("This fit has no direct effects; fit them with `direct_effects`.",
                        class = "latents_no_direct_effects", call = NULL))
  }
  errors <- .multilpa_table_errors(x, data)
  term <- paste(effects$indicator, effects$covariate, sep = ":")
  wanted <- if (is.null(errors)) NULL else
    errors[errors$parameter == "direct_effect", , drop = FALSE]
  matched <- if (is.null(wanted)) rep(NA_integer_, length(term)) else match(term, wanted$term)
  column <- function(name) if (is.null(wanted)) rep(NA_real_, length(term)) else wanted[[name]][matched]
  data.frame(indicator = effects$indicator, covariate = effects$covariate,
             estimate = unname(x$direct_effects),
             standard_error = column("standard_error"), statistic = column("statistic"),
             p_value = column("p_value"), conf_low = column("conf_low"),
             conf_high = column("conf_high"), stringsAsFactors = FALSE)
}
