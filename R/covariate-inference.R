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
  c(as.vector(t(centred_means)), as.vector(t(log(variances))),
    as.vector(object$profile_coefficients), as.vector(object$group_coefficients))
}

#' Unpack a parameter vector for a covariate fit
#' @param theta Numeric vector as built by [.multilpa_cov_encode()].
#' @param object The fit supplying the shapes.
#' @return A list with `parameters` (means, variances), `beta` and `gamma`.
#' @noRd
.multilpa_cov_decode <- function(theta, object) {
  n_profiles <- object$n_profiles
  n_indicators <- length(object$indicators)
  n_variance_rows <- if (identical(object$variance_model, "equal")) 1L else n_profiles
  at <- 0L
  take <- function(n) {
    values <- theta[at + seq_len(n)]
    at <<- at + n
    values
  }
  means <- matrix(take(n_profiles * n_indicators), n_profiles, n_indicators,
                  byrow = TRUE)
  variances <- exp(matrix(take(n_variance_rows * n_indicators), n_variance_rows,
                          n_indicators, byrow = TRUE))
  if (n_variance_rows == 1L) {
    variances <- matrix(variances, n_profiles, n_indicators, byrow = TRUE)
  }
  beta <- matrix(take(length(object$profile_coefficients)),
                 nrow(object$profile_coefficients), ncol(object$profile_coefficients))
  gamma <- matrix(take(length(object$group_coefficients)),
                  nrow(object$group_coefficients), ncol(object$group_coefficients))
  list(parameters = list(means = means, variances = variances),
       beta = beta, gamma = gamma)
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

  ## Measurement: means, then log variances.
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

#' Tidy inference for a covariate fit
#'
#' Standard errors, tests and intervals for every free parameter of
#' [fit_covariates()], at all three levels at once: the measurement model, the
#' profile logits and the group-class logits. Without them a membership
#' coefficient cannot be reported, because nothing distinguishes a real effect
#' from separation.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param data The data frame the model was fitted to.
#' @param level Confidence level for the intervals, between 0 and 1.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` uses the observed information;
#'   `"robust"` uses the group-clustered sandwich, which is the appropriate
#'   choice when the measurement model may be misspecified.
#' @return A base `data.frame` with one row per free parameter and the columns
#'   `level` (`"measurement"`, `"profile"` or `"group"`), `outcome` (which
#'   profile or group class the coefficient predicts, or which profile a
#'   measurement parameter belongs to), `term`, `parameter` (`"mean"`,
#'   `"variance"` or `"coefficient"`), `estimate`, `standard_error`, `statistic`,
#'   `p_value`, `conf_low` and `conf_high`. Variances are reported in their
#'   natural units, with standard errors carried through the delta method from
#'   the log scale on which they are estimated.
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
#' profile <- ifelse(runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- fit_covariates(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                       n_group_classes = 2, profile_covariates = "x",
#'                       n_starts = 2, seed = 1)
#' parameter_inference(fit, example_data)
#' @rdname parameter_inference
#' @export
parameter_inference.multilpa_covariates <- function(object, data, level = 0.95,
                                                    step = 1e-4,
                                                    vcov_type = c("observed", "robust")) {
  stopifnot(
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1
  )
  vcov_type <- match.arg(vcov_type)
  covariance <- .multilpa_cov_covariance(object, data, step, vcov_type)
  .multilpa_cov_inference_frame(object, .multilpa_cov_encode(object), covariance,
                                level, vcov_type)
}

#' Covariance of a covariate fit's free parameters
#' @param object A fitted `multilpa_covariates` model.
#' @param data The data frame the model was fitted to.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` or `"robust"`.
#' @return A square matrix, ordered as [.multilpa_cov_encode()].
#' @noRd
.multilpa_cov_covariance <- function(object, data, step, vcov_type) {
  stopifnot(
    "`object` must be a fitted `multilpa_covariates` model" =
      inherits(object, "multilpa_covariates"),
    "`data` must be a data frame" = is.data.frame(data),
    "`step` must be a single positive number" =
      is.numeric(step) && length(step) == 1L && is.finite(step) && step > 0
  )
  if (is.null(object$profile_design) || is.null(object$group_design)) {
    stop(errorCondition(
      "This fit predates covariate inference; refit with the current version.",
      class = "multilpa_unsupported_inference", call = NULL))
  }
  if (!all(object$indicators %in% names(data))) {
    stop(errorCondition(
      "`data` must contain every indicator the model was fitted to.",
      class = "multilpa_bad_inference_data", call = NULL))
  }
  x <- sweep(as.matrix(data[object$indicators]), 2L, object$center, "-")
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
  scale <- rep(1, length(theta))
  information <- .multilpa_observed_hessian(
    function(displacement) objective(theta + displacement),
    function(displacement) gradient(theta + displacement), scale, step)
  covariance <- information$inverse
  if (identical(vcov_type, "robust")) {
    scores <- .multilpa_cov_group_scores(theta, x, object)
    covariance <- covariance %*% crossprod(scores) %*% covariance
  }
  dimnames(covariance) <- list(.multilpa_cov_parameter_names(object),
                               .multilpa_cov_parameter_names(object))
  covariance
}

#' Names for the free parameters of a covariate fit
#' @param object A fitted `multilpa_covariates` model.
#' @return A character vector, ordered as [.multilpa_cov_encode()].
#' @noRd
.multilpa_cov_parameter_names <- function(object) {
  labels <- .multilpa_cov_labels(object)
  ## The kind belongs in the name: a mean and a variance share a level, an
  ## outcome and a term, so leaving it out makes them indistinguishable and
  ## gives vcov() duplicate dimnames.
  paste(labels$level, labels$parameter, labels$outcome, labels$term, sep = ".")
}

#' Covariance matrix of a covariate fit
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param data The data frame the model was fitted to.
#' @param step Finite-difference step for the observed information.
#' @param vcov_type `"observed"` or `"robust"`.
#' @param ... Ignored, present for generic compatibility.
#' @return A square numeric matrix with one row and column per free parameter,
#'   named `level.parameter.outcome.term` and ordered measurement means, measurement
#'   variances, profile logits, group logits. Variances are on the log scale,
#'   which is where they are estimated; [covariate_inference()] returns them and
#'   their standard errors in natural units.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- fit_covariates(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                       n_group_classes = 2, profile_covariates = "x",
#'                       n_starts = 2, seed = 1)
#' dim(vcov(fit, example_data))
#' @export
vcov.multilpa_covariates <- function(object, data, step = 1e-4,
                                     vcov_type = c("observed", "robust"), ...) {
  .multilpa_cov_covariance(object, data, step, match.arg(vcov_type))
}

#' Assemble the tidy inference table
#' @param object The fit; `theta` the estimates; `covariance` their covariance.
#' @param level Confidence level. `vcov_type` is recorded as an attribute.
#' @return One row per free parameter, at every level.
#' @noRd
.multilpa_cov_inference_frame <- function(object, theta, covariance, level,
                                          vcov_type) {
  errors <- sqrt(pmax(diag(covariance), 0))
  indicators <- object$indicators
  profiles <- paste0("profile_", seq_len(object$n_profiles))
  equal <- identical(object$variance_model, "equal")
  variance_rows <- if (equal) "shared" else profiles

  labels <- .multilpa_cov_labels(object)
  stopifnot("every parameter must be labelled" = nrow(labels) == length(theta))

  ## Variances are estimated as logs; the delta method returns them and their
  ## errors to natural units, where a reader can compare them with the data.
  is_variance <- labels$parameter == "variance"
  is_mean <- labels$parameter == "mean"
  estimate <- theta
  estimate[is_mean] <- estimate[is_mean] + rep(object$center, times = object$n_profiles)
  estimate[is_variance] <- exp(theta[is_variance])
  errors[is_variance] <- errors[is_variance] * estimate[is_variance]

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
  result$statistic[is_variance] <- NA_real_
  result$p_value[is_variance] <- NA_real_
  attr(result, "vcov_type") <- vcov_type
  attr(result, "level") <- level
  result
}

#' Level, outcome, term and kind for every free parameter
#' @param object A fitted `multilpa_covariates` model.
#' @return A data frame in the order [.multilpa_cov_encode()] packs them.
#' @noRd
.multilpa_cov_labels <- function(object) {
  indicators <- object$indicators
  profiles <- paste0("profile_", seq_len(object$n_profiles))
  variance_rows <- if (identical(object$variance_model, "equal")) "shared" else profiles
  block <- function(term, outcome, level, parameter) {
    grid <- expand.grid(term = term, outcome = outcome, stringsAsFactors = FALSE)
    data.frame(level = level, outcome = grid$outcome, term = grid$term,
               parameter = parameter, stringsAsFactors = FALSE)
  }
  rbind(
    block(indicators, profiles, "measurement", "mean"),
    block(indicators, variance_rows, "measurement", "variance"),
    block(rownames(object$profile_coefficients),
          colnames(object$profile_coefficients), "profile", "coefficient"),
    block(rownames(object$group_coefficients),
          colnames(object$group_coefficients), "group", "coefficient")
  )
}

#' Estimated parameters of a covariate fit
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param ... Ignored, present for generic compatibility.
#' @return A named numeric vector of every free parameter, in the order
#'   [parameter_inference()] reports them: measurement means, measurement
#'   variances, profile logits, then group-class logits. Names are
#'   `level.parameter.outcome.term`, which is what keeps a mean and a variance
#'   on the same indicator distinguishable. Variances are in natural units.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- fit_covariates(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                       n_group_classes = 2, profile_covariates = "x",
#'                       n_starts = 2, seed = 1)
#' coef(fit)
#' @export
coef.multilpa_covariates <- function(object, ...) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"))
  theta <- .multilpa_cov_encode(object)
  labels <- .multilpa_cov_labels(object)
  ## Means are stored centred for the likelihood and reported in input units;
  ## variances are estimated as logs and reported in their own units.
  theta[labels$parameter == "mean"] <- as.vector(t(object$means))
  is_variance <- labels$parameter == "variance"
  theta[is_variance] <- exp(theta[is_variance])
  stats::setNames(theta, .multilpa_cov_parameter_names(object))
}

#' Wald confidence intervals for a covariate fit
#'
#' @param object A fitted `multilpa_covariates` model.
#' @param parm Optional parameter names or indices; defaults to all of them.
#' @param level Confidence level strictly between zero and one.
#' @param data The data frame the model was fitted to.
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
#' profile <- ifelse(runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- fit_covariates(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                       n_group_classes = 2, profile_covariates = "x",
#'                       n_starts = 2, seed = 1)
#' confint(fit, data = example_data)
#' @export
confint.multilpa_covariates <- function(object, parm, level = 0.95, data, ...) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"),
            "`data` must be supplied; it is needed to rebuild the information" =
              !missing(data) && is.data.frame(data))
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
