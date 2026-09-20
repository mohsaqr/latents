#' Baseline-category multinomial probabilities
#' @param design Numeric design matrix.
#' @param coefficients Coefficients for all except the final category.
#' @return Row probabilities.
#' @noRd
.multilpa_softmax <- function(design, coefficients) {
  exp(.multilpa_log_softmax(design, coefficients))
}

#' Log multinomial probabilities without losing finite rare-class logits
#' @param design Numeric design matrix.
#' @param coefficients Coefficients for all except the final category.
#' @return Matrix of log probabilities, one row per design row.
#' @noRd
.multilpa_log_softmax <- function(design, coefficients) {
  stopifnot(is.matrix(design), is.matrix(coefficients),
            ncol(design) == nrow(coefficients))
  scores <- cbind(design %*% coefficients, 0)
  scores <- sweep(scores, 1L, apply(scores, 1L, max), "-")
  sweep(scores, 1L, log(rowSums(exp(scores))), "-")
}

#' Fit weighted multinomial logits
#' @param design Design matrix.
#' @param counts Expected category counts in each row.
#' @param initial Initial coefficient matrix.
#' @return Updated coefficients and optimization diagnostics.
#' @noRd
.multilpa_weighted_logits <- function(design, counts, initial) {
  stopifnot(is.matrix(design), is.matrix(counts), is.matrix(initial),
            nrow(design) == nrow(counts), all(counts >= 0))
  if (ncol(counts) == 1L) return(list(coefficients = initial, converged = TRUE))
  unpack <- function(values) {
    stopifnot(is.numeric(values))
    matrix(values, ncol(design), ncol(counts) - 1L)
  }
  objective <- function(values) {
    stopifnot(is.numeric(values))
    -sum(counts * .multilpa_log_softmax(design, unpack(values)))
  }
  gradient <- function(values) {
    stopifnot(is.numeric(values))
    residual <- .multilpa_softmax(design, unpack(values)) * rowSums(counts) - counts
    as.vector(crossprod(design, residual[, seq_len(ncol(counts) - 1L), drop = FALSE]))
  }
  fit <- stats::optim(as.vector(initial), objective, gradient, method = "BFGS",
                      control = list(maxit = 500L, reltol = 1e-11))
  if (fit$value > objective(as.vector(initial)) + 1e-8) {
    stop("Multinomial M-step decreased the expected log likelihood.")
  }
  list(coefficients = unpack(fit$par), converged = fit$convergence == 0L)
}

#' Covariate-dependent nested expectation
#' @param x Complete indicator matrix.
#' @param group_index Group indices.
#' @param parameters Measurement parameters.
#' @param profile_design List of profile design matrices, one per group class.
#' @param group_design Group-level design matrix.
#' @param beta Profile logit coefficients.
#' @param gamma Group logit coefficients.
#' @return Log likelihood and posterior responsibilities.
#' @noRd
.multilpa_cov_expectation <- function(x, group_index, parameters, profile_design,
                                    group_design, beta, gamma, codes = NULL) {
  stopifnot(is.matrix(x), !anyNA(x), is.list(parameters), is.list(profile_design),
            is.matrix(group_design), is.matrix(beta), is.matrix(gamma))
  n_profiles <- nrow(parameters$means)
  # The measurement model is the covariate-free one; only the mixing weights
  # differ here. Residual covariances need the conditional moments, which the
  # maximization step then reuses, so they are computed once and carried.
  gaussian <- if (ncol(x) > 0L && !is.null(parameters$covariances)) {
    .multilpa_gaussian_moments(x, parameters)
  } else NULL
  log_density <- if (!is.null(gaussian)) gaussian$log_density else {
    matrix(vapply(seq_len(n_profiles), function(k) {
      residual <- sweep(x, 2L, parameters$means[k, ], "-")
      -0.5 * rowSums(sweep(residual^2, 2L, parameters$variances[k, ], "/") +
                      matrix(log(2 * pi * parameters$variances[k, ]),
                             nrow(x), ncol(x), byrow = TRUE))
    }, numeric(nrow(x))), nrow(x), n_profiles)
  }
  density_offset <- apply(log_density, 1L, max)
  log_density <- sweep(log_density, 1L, density_offset, "-")
  # Categorical indicators are conditionally independent of the continuous ones
  # given the profile, and are added after the Gaussian offset comes off so
  # their contribution keeps its precision.
  if (!is.null(codes)) {
    log_density <- log_density +
      .multilpa_categorical_log_density(codes, parameters$response_probabilities)
  }
  conditional <- lapply(profile_design, function(design) {
    log_prior <- .multilpa_log_softmax(design, beta)
    prior <- exp(log_prior)
    scores <- log_density + log_prior
    offset <- apply(scores, 1L, max)
    weights <- exp(sweep(scores, 1L, offset, "-"))
    total <- rowSums(weights)
    marginal <- offset + log(total)
    list(prior = prior, marginal = marginal,
         posterior = weights / total)
  })
  group_log_prior <- .multilpa_log_softmax(group_design, gamma)
  group_prior <- exp(group_log_prior)
  evidence <- matrix(vapply(seq_along(conditional), function(h) {
    as.vector(rowsum(conditional[[h]]$marginal, group_index, reorder = FALSE))
  }, numeric(nrow(group_design))), nrow(group_design))
  evidence_offset <- apply(evidence, 1L, max)
  scores <- sweep(evidence, 1L, evidence_offset, "-") + group_log_prior
  score_offset <- apply(scores, 1L, max)
  weights <- exp(sweep(scores, 1L, score_offset, "-"))
  total <- rowSums(weights)
  group_posteriors <- weights / total
  group_log_likelihood <- evidence_offset + score_offset + log(total) +
    as.vector(rowsum(density_offset, group_index, reorder = FALSE))
  joint <- lapply(seq_along(conditional), function(h) {
    conditional[[h]]$posterior * group_posteriors[group_index, h]
  })
  list(log_likelihood = sum(group_log_likelihood),
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors, subject_posteriors = Reduce(`+`, joint),
       joint = joint, group_priors = group_prior,
       profile_priors = lapply(conditional, `[[`, "prior"),
       gaussian_moments = gaussian$moments)
}

#' Fit multilevel LPA with class-membership covariates
#'
#' Numeric covariates predict individual-profile and group-class membership via
#' multinomial logits, using the final class as reference. Individual-profile
#' slopes are shared across group classes; profile intercepts differ by group
#' class. The measurement model is the one [multilpa()] fits and stays
#' invariant across group classes: Gaussian, categorical or mixed indicators,
#' with diagonal or unrestricted residual covariance. This is one-step maximum
#' likelihood, not regression on assigned classes. Covariates are used in their
#' supplied units; center or scale them beforehand if desired. Only complete
#' data are supported here; missing indicators are not integrated out.
#' @param data Data frame.
#' @param vars Names of continuous indicator columns.
#' @param id Name of group identifier column.
#' @param n_profiles Number of individual profiles.
#' @param n_group_classes Number of discrete group classes.
#' @param profile_covariates Names of numeric predictors of profile membership.
#' @param group_covariates Names of numeric predictors constant within groups.
#' @param variance_model Varying or equal diagonal variances across profiles.
#' @param n_starts Number of independent initializations.
#' @param max_iter Maximum EM iterations per start.
#' @param tol Relative log-likelihood convergence tolerance.
#' @param min_variance Explicit variance lower bound.
#' @param time Optional name of a column giving each observation's position
#'   within its group, stored so that [sequences()] and `plot(what =
#'   "sequences")` can read the assignments back in order. The model never uses
#'   it.
#' @param seed Optional seed; the caller's random state is restored.
#' @param covariance_model `"diagonal"` assumes indicators are independent
#'   given the profile; `"full"` estimates unrestricted within-profile residual
#'   covariances, as in [multilpa()].
#' @param categorical Character vector naming indicators to treat as
#'   categorical, modelled by unrestricted profile-specific response
#'   probabilities. Indicators not named here stay Gaussian, so naming a subset
#'   fits a mixed-mode measurement model.
#' @param min_probability Positive lower bound on every categorical response
#'   probability, defining a constrained maximum-likelihood problem in the same
#'   way `min_variance` does for Gaussian indicators.
#' @return An object of class `multilpa_covariates`. Read it with the verbs
#'   that describe it rather than by reaching into it: [as.data.frame()] gives
#'   the measurement model, the membership coefficients, and either set of
#'   posteriors, one row per profile-indicator, coefficient, individual or
#'   group; [summary()] gives the model-level fit summary and the restart
#'   diagnostics; [plot()] draws the measurement model; and
#'   [parameter_inference()] gives one row per free parameter with a standard
#'   error and an interval.
#' @details [parameter_inference()], [vcov()] and [confint()] cover Gaussian
#'   and full-covariance fits. They refuse a fit with categorical indicators
#'   with a `multilpa_unsupported_inference` condition, because no score is
#'   implemented for the response probabilities and reporting the other blocks
#'   alone would understate the parameter count.
#' @examples
#' set.seed(1)
#' d <- data.frame(group = rep(1:20, each = 10), z = rnorm(200))
#' d$y <- rnorm(200, ifelse(runif(200) < plogis(d$z), -3, 3))
#' fit <- fit_covariates(d, "y", "group", 2, 1,
#'                             profile_covariates = "z", n_starts = 2, seed = 1)
#' as.data.frame(fit, what = "coefficients")
#' @export
fit_covariates <- function(data, vars, id, n_profiles,
                                  n_group_classes = 2L,
                                  profile_covariates = character(),
                                  group_covariates = character(),
                                  variance_model = c("varying", "equal"),
                                  n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                                  min_variance = 1e-6, seed = NULL,
                                  time = NULL,
                                  covariance_model = c("diagonal", "full"),
                                  categorical = character(),
                                  min_probability = 1e-10) {
  stopifnot(is.data.frame(data), is.character(vars), is.character(id),
            is.character(profile_covariates), is.character(group_covariates),
            !anyDuplicated(profile_covariates), !anyDuplicated(group_covariates),
            all(c(profile_covariates, group_covariates) %in% names(data)),
            is.numeric(n_starts), length(n_starts) == 1L,
            is.finite(n_starts), n_starts >= 1, n_starts == as.integer(n_starts))
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  stopifnot(
    "`categorical` must be a character vector of indicator names" =
      is.character(categorical) && !anyNA(categorical),
    "`min_probability` must be a single number in (0, 1)" =
      is.numeric(min_probability) && length(min_probability) == 1L &&
      is.finite(min_probability) && min_probability > 0 && min_probability < 1)
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed),
              seed >= 0, seed <= .Machine$integer.max, seed == as.integer(seed))
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
              rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  .multilpa_cov_check_covariates(data, profile_covariates, group_covariates,
                               vars, id)
  # The base fit validates indicators/model sizes and supplies an initial mode.
  base <- multilpa(data, vars, id, n_profiles, n_group_classes,
                     variance_model, n_starts = 1L, max_iter = max_iter,
                     tol = tol, min_variance = min_variance,
                     covariance_model = covariance_model,
                     categorical = categorical,
                     min_probability = min_probability)
  group_index <- base$group_index
  first_rows <- match(seq_len(base$n_groups), group_index)
  if (length(group_covariates) && any(vapply(group_covariates, function(name) {
    any(data[[name]] != data[[name]][first_rows][group_index])
  }, logical(1)))) stop("Every group_covariate must be constant within each group.")
  if (n_profiles == 1L && length(profile_covariates)) {
    stop("Profile covariates require at least two profiles.")
  }
  if (n_group_classes == 1L && length(group_covariates)) {
    stop("Group covariates require at least two group classes.")
  }
  designs <- .multilpa_cov_designs(data, vars, profile_covariates,
                                 group_covariates, first_rows, n_group_classes,
                                 categorical, min_probability)
  x <- designs$x
  center <- designs$center
  control <- list(variance_model = variance_model, min_variance = min_variance,
                  max_iter = max_iter, tol = tol,
                  covariance_model = covariance_model,
                  min_probability = min_probability,
                  n_profile_covariates = length(profile_covariates),
                  n_group_covariates = length(group_covariates))
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch(.multilpa_cov_start(start_index, x, group_index, n_profiles,
                               n_group_classes, designs, base, control),
             error = function(error) list(error = conditionMessage(error)))
  })
  starts <- do.call(rbind, lapply(seq_along(attempts), function(i) {
    result <- attempts[[i]]
    data.frame(start = i, log_likelihood = if (is.null(result$expectation)) -Inf else
                 result$expectation$log_likelihood,
               converged = isTRUE(result$converged),
               iterations = result$iterations %||% 0L,
               error = result$error)
  }))
  if (!any(is.finite(starts$log_likelihood))) {
    stop(sprintf("All covariate starts failed: %s", paste(unique(starts$error), collapse = "; ")))
  }
  best_index <- which.max(starts$log_likelihood)
  result <- .multilpa_cov_assemble(
    best = attempts[[best_index]], best_index = best_index, starts = starts,
    designs = designs, base = base, data = data, vars = vars,
    id = id, profile_covariates = profile_covariates,
    group_covariates = group_covariates,
    # Stored as integers so the field has the same type it has on every other fit
    # class; a double here made `identical()` on any count derived from it fail.
    n_profiles = as.integer(n_profiles),
    n_group_classes = as.integer(n_group_classes), group_index = group_index,
    variance_model = variance_model, min_variance = min_variance,
    covariance_model = covariance_model, categorical = categorical,
    call = match.call())
  ## Set here rather than inside the assembler, where the name `time` would
  ## resolve to stats::time instead of this argument.
  result$time <- time
  result$time_values <- .multilpa_time_values(data, time, id, vars)
  ## The same effective counts the Gaussian fit carries, so the shared
  ## diagnostics and plot panels need no special case for this class.
  result$effective_profile_counts <- colSums(result$subject_posteriors)
  result$effective_group_counts <- colSums(result$group_posteriors)
  if (any(!is.finite(starts$log_likelihood))) {
    warning(warningCondition(paste(
      "Some covariate starts failed; summary() reports every start."),
      class = "multilpa_failed_starts"))
  }
  if (!result$converged) warning("Best covariate fit did not converge.")
  if (result$boundary) warning("A residual variance reached min_variance.")
  if (result$extreme_logits) warning("Extreme logit coefficients: inspect scaling, sparse classes and separation.")
  result
}

#' Print a covariate LPA fit
#' @param x A covariate LPA fit.
#' @param ... Reserved.
#' @return The model invisibly.
#' @examples
#' # print(fit)
#' @export
print.multilpa_covariates <- function(x, ...) {
  stopifnot(inherits(x, "multilpa_covariates"))
  cat(sprintf("Multilevel LPA with covariates: %d profiles, %d group classes\n",
              x$n_profiles, x$n_group_classes))
  cat(sprintf("Log likelihood %.6f; AIC %.3f; BIC (groups) %.3f; converged %s\n",
              x$log_likelihood, x$aic, x$bic, x$converged))
  invisible(x)
}

#' Summarize a covariate LPA fit
#'
#' Collects the model-level fit, the measurement model, the membership
#' regressions and the restart diagnostics into one object, so that none of
#' them has to be read out of the fit by hand.
#'
#' @param object A covariate LPA fit from [fit_covariates()].
#' @param ... Reserved for compatibility with `summary()`.
#' @return An object of class `summary_multilpa_covariates`, with a `print`
#'   method and an [as.data.frame()] accessor. `as.data.frame()` returns the
#'   one-row model summary by default; `what = "profiles"`, `"coefficients"`
#'   and `"starts"` return the measurement model, the membership coefficients
#'   and the restart diagnostics. The membership coefficients carry no standard
#'   errors here; [parameter_inference()] reports those.
#' @examples
#' set.seed(1)
#' d <- data.frame(group = rep(1:20, each = 10), z = rnorm(200))
#' d$y <- rnorm(200, ifelse(runif(200) < plogis(d$z), -3, 3))
#' fit <- fit_covariates(d, "y", "group", 2, 1,
#'                       profile_covariates = "z", n_starts = 2, seed = 1)
#' summary(fit)
#' as.data.frame(summary(fit), what = "coefficients")
#' @export
summary.multilpa_covariates <- function(object, ...) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"))
  ## Built field by field rather than by subsetting the fit by a name vector,
  ## which yields a silent NULL keyed `NA` for any name the fit does not carry.
  model <- data.frame(
    n_observations = object$n_observations,
    n_groups = object$n_groups,
    n_profiles = object$n_profiles,
    n_group_classes = object$n_group_classes,
    n_profile_covariates = length(object$profile_covariates),
    n_group_covariates = length(object$group_covariates),
    variance_model = object$variance_model,
    covariance_model = object$covariance_model %||% "diagonal",
    n_parameters = object$n_parameters,
    log_likelihood = object$log_likelihood,
    aic = object$aic,
    bic_groups = object$bic,
    bic_individual = object$bic_individual,
    converged = object$converged,
    boundary = object$boundary,
    extreme_logits = object$extreme_logits,
    n_starts = nrow(object$starts),
    row.names = NULL, stringsAsFactors = FALSE)
  result <- list(
    model = model,
    profiles = as.data.frame(object, what = "profiles"),
    coefficients = as.data.frame(object, what = "coefficients"),
    starts = object$starts,
    effective_profile_counts = object$effective_profile_counts,
    effective_group_counts = object$effective_group_counts,
    call = object$call)
  class(result) <- "summary_multilpa_covariates"
  result
}

#' Print a covariate LPA summary
#' @param x A `summary_multilpa_covariates` object.
#' @param digits Number of printed significant digits.
#' @param ... Passed to the underlying `data.frame` printing.
#' @return The summary, invisibly.
#' @examples
#' # After fitting: print(summary(fit), digits = 3)
#' @export
print.summary_multilpa_covariates <- function(x, digits = 4L, ...) {
  stopifnot("`x` must be a `summary_multilpa_covariates` object" =
              inherits(x, "summary_multilpa_covariates"),
            "`digits` must be a single number between 1 and 22" =
              is.numeric(digits) && length(digits) == 1L && is.finite(digits) &&
              digits >= 1 && digits <= 22)
  model <- x$model
  cat(sprintf("Multilevel LPA with covariates: %d profiles, %d group classes\n",
              model$n_profiles, model$n_group_classes))
  cat(sprintf("Individuals: %d; groups: %d; parameters: %d; converged: %s\n",
              model$n_observations, model$n_groups, model$n_parameters,
              model$converged))
  cat(sprintf("%d profile covariate(s); %d group covariate(s)\n",
              model$n_profile_covariates, model$n_group_covariates))
  cat("\nMeasurement model:\n")
  print(x$profiles, digits = digits, row.names = FALSE, ...)
  cat("\nMembership regressions (no standard errors; see parameter_inference()):\n")
  print(x$coefficients, digits = digits, row.names = FALSE, ...)
  cat("\nEffective individual memberships:\n")
  print(x$effective_profile_counts, digits = digits, ...)
  cat("\nEffective group memberships:\n")
  print(x$effective_group_counts, digits = digits, ...)
  cat(sprintf("\nLog likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              model$log_likelihood, model$aic, model$bic_groups,
              model$bic_individual))
  if (!model$converged) cat("WARNING: the best start did not converge.\n")
  if (model$boundary) cat("WARNING: a residual variance is at min_variance.\n")
  if (model$extreme_logits) {
    cat("WARNING: extreme logit coefficients; check scaling, sparse classes and separation.\n")
  }
  cat("\nStart diagnostics:\n")
  print(x$starts, digits = digits, row.names = FALSE, ...)
  invisible(x)
}

#' Tidy a covariate LPA summary
#' @param x A `summary_multilpa_covariates` object.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return: `"model"`, `"profiles"`,
#'   `"coefficients"` or `"starts"`.
#' @param ... Ignored.
#' @return A base `data.frame`. `"model"` has exactly one row describing the
#'   fit, with columns `n_observations`, `n_groups`, `n_profiles`,
#'   `n_group_classes`, `n_profile_covariates`, `n_group_covariates`,
#'   `variance_model`, `covariance_model`, `n_parameters`, `log_likelihood`,
#'   `aic`, `bic_groups`, `bic_individual`, `converged`, `boundary`,
#'   `extreme_logits` and `n_starts`. `"profiles"`, `"coefficients"` and
#'   `"starts"` have one row per profile-indicator, per membership coefficient,
#'   and per EM start respectively.
#' @examples
#' # After fitting: as.data.frame(summary(fit), what = "starts")
#' @export
as.data.frame.summary_multilpa_covariates <- function(x, row.names = NULL,
                                                      optional = FALSE,
                                                      what = c("model",
                                                               "profiles",
                                                               "coefficients",
                                                               "starts"), ...) {
  stopifnot("`x` must be a `summary_multilpa_covariates` object" =
              inherits(x, "summary_multilpa_covariates"))
  what <- match.arg(what)
  result <- switch(what, model = x$model, profiles = x$profiles,
                   coefficients = x$coefficients, starts = x$starts)
  row.names(result) <- row.names
  result
}

#' Rebuild the fitting data a covariate fit was estimated from
#'
#' A covariate fit stores its centred indicator matrix, both membership design
#' matrices and the group index, which between them carry every column the
#' likelihood reads. This reassembles that data frame so the inference verbs
#' need not demand data the object already owns. It is exact, not approximate:
#' the indicators come back as stored, and the covariates are read off the
#' designs they were built into.
#'
#' @param object A fitted `multilpa_covariates` model.
#' @return A `data.frame` carrying the indicators, the group identifier and
#'   every profile and group covariate, in fitting row order.
#' @noRd
.multilpa_cov_stored_data <- function(object) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"))
  if (is.null(object$indicator_data) || is.null(object$profile_design) ||
      is.null(object$group_design)) {
    stop(errorCondition(paste(
      "This fit does not store its indicators and membership designs, so the",
      "fitting data cannot be rebuilt; supply `data`."),
      class = "multilpa_no_indicator_data", call = NULL))
  }
  result <- as.data.frame(object$indicator_data)
  result[[object$id]] <- object$group_values[object$group_index]
  ## The profile design is the group-class indicator block followed by the
  ## profile covariates, in the order they were named.
  profile_block <- object$profile_design[[1L]][,
    -seq_len(object$n_group_classes), drop = FALSE]
  if (ncol(profile_block) > 0L) {
    colnames(profile_block) <- object$profile_covariates
    result <- cbind(result, as.data.frame(profile_block))
  }
  ## The group design is an intercept followed by the group covariates, held at
  ## one row per group, so it expands back through the group index.
  group_block <- object$group_design[, -1L, drop = FALSE]
  if (ncol(group_block) > 0L) {
    colnames(group_block) <- object$group_covariates
    result <- cbind(result,
                    as.data.frame(group_block[object$group_index, , drop = FALSE]))
  }
  result
}

#' Extract a covariate LPA log likelihood
#' @param object A covariate LPA fit.
#' @param ... Reserved.
#' @return A logLik object using independent groups for nobs.
#' @examples
#' # logLik(fit)
#' @export
logLik.multilpa_covariates <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_covariates"))
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' Count independent groups in a covariate LPA fit
#' @param object A covariate LPA fit.
#' @param ... Reserved.
#' @return Number of observed groups.
#' @examples
#' # nobs(fit)
#' @export
nobs.multilpa_covariates <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_covariates"))
  object$n_groups
}

#' Validate covariate columns
#'
#' Covariates enter the membership regressions directly, so a non-finite or
#' non-numeric column would propagate silently into the logits.
#'
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_cov_check_covariates <- function(data, profile_covariates,
                                         group_covariates, vars, id) {
  predictors <- unique(c(profile_covariates, group_covariates))
  if (length(predictors) && !all(vapply(data[predictors], function(column) {
    is.numeric(column) && is.null(dim(column)) && all(is.finite(column))
  }, logical(1)))) {
    stop("Covariates must be finite numeric columns without missing values.")
  }
  if (any(predictors %in% c(vars, id))) {
    stop("Covariates must be distinct from indicators and the group identifier.")
  }
  invisible(NULL)
}

#' Build the centered indicator matrix and both membership designs
#'
#' The profile design repeats one group-class indicator block per group class,
#' so that stacking it gives the weighted multinomial regression its rows. The
#' group design takes one row per group, at its first occurrence.
#'
#' @return A list with `x`, `center`, `w`, `profile_design` and `stacked_design`.
#' @noRd
.multilpa_cov_designs <- function(data, vars, profile_covariates,
                                group_covariates, first_rows, n_group_classes,
                                categorical = character(),
                                min_probability = 1e-10) {
  measurement <- .multilpa_prepare_indicators(data, vars, categorical,
                                              "error", min_probability)
  x <- measurement$x
  center <- if (ncol(x) > 0L) colMeans(x) else numeric(0)
  z <- as.matrix(data[profile_covariates])
  w <- cbind(`(Intercept)` = 1,
             as.matrix(data[first_rows, group_covariates, drop = FALSE]))
  profile_design <- lapply(seq_len(n_group_classes), function(group_class) {
    cbind(matrix(as.numeric(seq_len(n_group_classes) == group_class), nrow(data),
                 n_group_classes, byrow = TRUE), z)
  })
  stacked_design <- do.call(rbind, profile_design)
  if (qr(stacked_design)$rank != ncol(stacked_design) || qr(w)$rank != ncol(w)) {
    stop("Covariate design is rank deficient; remove constant or collinear predictors.")
  }
  list(x = sweep(x, 2L, center, "-"), center = center, w = w,
       profile_design = profile_design, stacked_design = stacked_design,
       continuous = measurement$continuous, codes = measurement$codes,
       n_categories = measurement$n_categories,
       categorical_levels = measurement$encoded$levels)
}

#' Run one covariate EM start to convergence
#'
#' The first start resumes from the covariate-free base fit; later starts use a
#' fresh random initialization. Membership coefficients begin at the implied
#' logits with zero slopes, so the first iteration reproduces the base model.
#'
#' @return A list with the fitted parameters, expectation, both coefficient
#'   matrices, and the convergence diagnostics for this start.
#' @noRd
.multilpa_cov_start <- function(start_index, x, group_index, n_profiles,
                              n_group_classes, designs, base, control) {
  parameters <- if (start_index == 1L) {
    resumed <- list(means = sweep(unname(as.matrix(base$means)), 2L,
                                  designs$center, "-"),
                    variances = unname(as.matrix(base$variances)),
                    profile_probabilities = unname(as.matrix(base$profile_probabilities)),
                    group_probabilities = unname(base$group_probabilities))
    # The covariate-free fit already carries whichever measurement blocks the
    # model has; resuming from it must not silently drop them.
    if (!is.null(base$covariances)) resumed$covariances <- unname(base$covariances)
    if (!is.null(base$response_probabilities)) {
      resumed$response_probabilities <- unname(lapply(base$response_probabilities,
                                                      function(block) unname(as.matrix(block))))
    }
    resumed
  } else {
    .multilpa_initialize(x, group_index, n_profiles, n_group_classes,
                       control$variance_model, control$min_variance, start_index,
                       control$covariance_model, designs$codes,
                       designs$n_categories, control$min_probability)
  }
  beta <- rbind(
    log(parameters$profile_probabilities[, seq_len(n_profiles - 1L), drop = FALSE] /
          parameters$profile_probabilities[, n_profiles]),
    matrix(0, control$n_profile_covariates, n_profiles - 1L))
  gamma <- rbind(
    matrix(log(parameters$group_probabilities[seq_len(n_group_classes - 1L)] /
                 parameters$group_probabilities[n_group_classes]), 1L),
    matrix(0, control$n_group_covariates, n_group_classes - 1L))
  expectation <- .multilpa_cov_expectation(x, group_index, parameters,
                                         designs$profile_design, designs$w,
                                         beta, gamma, designs$codes)
  history <- expectation$log_likelihood
  iteration <- 0L
  converged <- FALSE
  while (iteration < control$max_iter && !converged) {
    parameters <- .multilpa_maximization(x, expectation, control$variance_model,
                                       control$min_variance,
                                       control$covariance_model, designs$codes,
                                       designs$n_categories,
                                       control$min_probability)
    profile_update <- .multilpa_weighted_logits(designs$stacked_design,
                                              do.call(rbind, expectation$joint), beta)
    group_update <- .multilpa_weighted_logits(designs$w,
                                            expectation$group_posteriors, gamma)
    beta <- profile_update$coefficients
    gamma <- group_update$coefficients
    updated <- .multilpa_cov_expectation(x, group_index, parameters,
                                       designs$profile_design, designs$w,
                                       beta, gamma, designs$codes)
    change <- updated$log_likelihood - expectation$log_likelihood
    if (change < -1e-9 * (1 + abs(expectation$log_likelihood))) {
      stop("Covariate EM decreased the observed log likelihood.")
    }
    iteration <- iteration + 1L
    # Both the likelihood and the inner logit steps must have settled, so a
    # stalled regression cannot be reported as a converged fit.
    converged <- abs(change) <= control$tol * (1 + abs(expectation$log_likelihood)) &&
      profile_update$converged && group_update$converged
    expectation <- updated
    history <- c(history, expectation$log_likelihood)
  }
  list(parameters = parameters, expectation = expectation, beta = beta,
       gamma = gamma, iterations = iteration, converged = converged,
       history = history, error = NA_character_)
}

#' Assemble the covariate fit result
#' @return An `multilpa_covariates` object.
#' @noRd
.multilpa_cov_assemble <- function(best, best_index, starts, designs, base, data,
                                 vars, id, profile_covariates,
                                 group_covariates, n_profiles, n_group_classes,
                                 group_index, variance_model, min_variance,
                                 covariance_model = "diagonal",
                                 categorical = character(), call) {
  result <- c(best$parameters,
              best$expectation[c("log_likelihood", "group_log_likelihood",
                                 "group_posteriors", "subject_posteriors",
                                 "group_priors", "profile_priors")])
  n_indicators <- ncol(designs$x)
  result$means <- sweep(result$means, 2L, designs$center, "+")
  result$profile_coefficients <- best$beta
  dimnames(result$profile_coefficients) <- list(
    c(paste0("group_class_", seq_len(n_group_classes)), profile_covariates),
    if (n_profiles > 1L) paste0("profile_", seq_len(n_profiles - 1L)))
  result$group_coefficients <- best$gamma
  dimnames(result$group_coefficients) <- list(
    colnames(designs$w),
    if (n_group_classes > 1L) paste0("group_class_", seq_len(n_group_classes - 1L)))
  # Mixing probabilities from the Gaussian M-step are not constant priors here.
  result$profile_probabilities <- NULL
  result$group_probabilities <- NULL
  # The membership logits replace the mixing probabilities the covariate-free
  # count includes, so the measurement half is taken from the shared counter and
  # the two coefficient matrices are added in their place.
  measurement_parameters <- .multilpa_count_parameters(
    n_profiles, n_group_classes, n_indicators, designs$n_categories,
    variance_model, covariance_model) -
    ((n_group_classes - 1L) + n_group_classes * (n_profiles - 1L))
  result$n_parameters <- length(best$beta) + length(best$gamma) +
    measurement_parameters
  result$aic <- -2 * result$log_likelihood + 2 * result$n_parameters
  result$bic <- -2 * result$log_likelihood + log(base$n_groups) * result$n_parameters
  result$bic_individual <- -2 * result$log_likelihood +
    log(nrow(data)) * result$n_parameters
  result$n_observations <- nrow(data)
  result$n_informative <- base$n_informative
  result$indicator_data <- as.matrix(data[vars])
  result$n_groups <- base$n_groups
  result$n_profiles <- as.integer(n_profiles)
  result$n_group_classes <- as.integer(n_group_classes)
  result$group_index <- group_index
  result$group_values <- base$group_values
  result$vars <- vars
  result$id <- id
  result$profile_covariates <- profile_covariates
  result$group_covariates <- group_covariates
  result$variance_model <- variance_model
  result$covariance_model <- covariance_model
  result$categorical <- categorical
  result$continuous <- designs$continuous
  result$categorical_levels <- designs$categorical_levels
  result$categorical_data <- designs$codes
  result$n_categories <- designs$n_categories
  result$measurement_model <- if (is.null(designs$codes)) "gaussian" else
    if (n_indicators == 0L) "categorical" else "mixed"
  result$standard_deviations <- sqrt(result$variances)
  if (!is.null(result$response_probabilities)) {
    names(result$response_probabilities) <- categorical
    result$response_probabilities <- stats::setNames(
      lapply(seq_along(categorical), function(i) {
        block <- result$response_probabilities[[i]]
        dimnames(block) <- list(paste0("profile_", seq_len(n_profiles)),
                                designs$categorical_levels[[i]])
        block
      }), categorical)
  }
  ## Carried so that covariate_inference() can rebuild the likelihood and
  ## its scores without re-deriving the designs from the data.
  result$profile_design <- designs$profile_design
  result$group_design <- designs$w
  result$center <- designs$center
  result$min_variance <- min_variance
  result$converged <- best$converged
  result$iterations <- best$iterations
  result$starts <- starts
  result$best_start <- best_index
  result$log_likelihood_history <- best$history
  result$subject_profiles <- max.col(result$subject_posteriors, ties.method = "first")
  result$group_classes <- max.col(result$group_posteriors, ties.method = "first")
  result$boundary <- .multilpa_covariance_boundary(
    list(means = result$means, variances = result$variances,
         covariances = result$covariances), min_variance)
  result$extreme_logits <- any(abs(c(best$beta, best$gamma)) > 20)
  result$call <- call
  structure(result, class = "multilpa_covariates")
}
