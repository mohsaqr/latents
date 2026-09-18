#' Baseline-category multinomial probabilities
#' @param design Numeric design matrix.
#' @param coefficients Coefficients for all except the final category.
#' @return Row probabilities.
#' @noRd
.ml_lpa_softmax <- function(design, coefficients) {
  stopifnot(is.matrix(design), is.matrix(coefficients),
            ncol(design) == nrow(coefficients))
  scores <- cbind(design %*% coefficients, 0)
  exp(sweep(scores, 1L, .ml_lpa_log_sum_exp(scores), "-"))
}

#' Fit weighted multinomial logits
#' @param design Design matrix.
#' @param counts Expected category counts in each row.
#' @param initial Initial coefficient matrix.
#' @return Updated coefficients and optimization diagnostics.
#' @noRd
.ml_lpa_weighted_logits <- function(design, counts, initial) {
  stopifnot(is.matrix(design), is.matrix(counts), is.matrix(initial),
            nrow(design) == nrow(counts), all(counts >= 0))
  if (ncol(counts) == 1L) return(list(coefficients = initial, converged = TRUE))
  unpack <- function(values) {
    stopifnot(is.numeric(values))
    matrix(values, ncol(design), ncol(counts) - 1L)
  }
  objective <- function(values) {
    stopifnot(is.numeric(values))
    scores <- cbind(design %*% unpack(values), 0)
    sum(rowSums(counts) * .ml_lpa_log_sum_exp(scores) - rowSums(counts * scores))
  }
  gradient <- function(values) {
    stopifnot(is.numeric(values))
    residual <- .ml_lpa_softmax(design, unpack(values)) * rowSums(counts) - counts
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
.ml_lpa_cov_expectation <- function(x, group_index, parameters, profile_design,
                                    group_design, beta, gamma) {
  stopifnot(is.matrix(x), !anyNA(x), is.list(parameters), is.list(profile_design),
            is.matrix(group_design), is.matrix(beta), is.matrix(gamma))
  log_density <- matrix(vapply(seq_len(nrow(parameters$means)), function(k) {
    residual <- sweep(x, 2L, parameters$means[k, ], "-")
    -0.5 * rowSums(sweep(residual^2, 2L, parameters$variances[k, ], "/") +
                    matrix(log(2 * pi * parameters$variances[k, ]),
                           nrow(x), ncol(x), byrow = TRUE))
  }, numeric(nrow(x))), nrow(x))
  conditional <- lapply(profile_design, function(design) {
    prior <- .ml_lpa_softmax(design, beta)
    scores <- log_density + log(prior)
    marginal <- .ml_lpa_log_sum_exp(scores)
    list(prior = prior, marginal = marginal,
         posterior = exp(sweep(scores, 1L, marginal, "-")))
  })
  group_prior <- .ml_lpa_softmax(group_design, gamma)
  scores <- matrix(vapply(seq_along(conditional), function(h) {
    as.vector(rowsum(conditional[[h]]$marginal, group_index, reorder = FALSE)) +
      log(group_prior[, h])
  }, numeric(nrow(group_design))), nrow(group_design))
  group_log_likelihood <- .ml_lpa_log_sum_exp(scores)
  group_posteriors <- exp(sweep(scores, 1L, group_log_likelihood, "-"))
  joint <- lapply(seq_along(conditional), function(h) {
    conditional[[h]]$posterior * group_posteriors[group_index, h]
  })
  list(log_likelihood = sum(group_log_likelihood),
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors, subject_posteriors = Reduce(`+`, joint),
       joint = joint, group_priors = group_prior,
       profile_priors = lapply(conditional, `[[`, "prior"))
}

#' Fit multilevel LPA with class-membership covariates
#'
#' Numeric covariates predict individual-profile and group-class membership via
#' multinomial logits, using the final class as reference. Individual-profile
#' slopes are shared across group classes; profile intercepts differ by group
#' class. Measurement means and diagonal variances remain invariant. This is
#' one-step maximum likelihood, not regression on assigned classes. Covariates
#' are used in their supplied units; center/scale them beforehand if desired.
#' Only complete data and diagonal residual variances are supported here.
#' @param data Data frame.
#' @param indicators Names of continuous indicator columns.
#' @param cluster Name of group identifier column.
#' @param n_profiles Number of individual profiles.
#' @param n_group_classes Number of discrete group classes.
#' @param profile_covariates Names of numeric predictors of profile membership.
#' @param group_covariates Names of numeric predictors constant within groups.
#' @param variance_model Varying or equal diagonal variances across profiles.
#' @param n_starts Number of independent initializations.
#' @param max_iter Maximum EM iterations per start.
#' @param tol Relative log-likelihood convergence tolerance.
#' @param min_variance Explicit variance lower bound.
#' @param seed Optional seed; the caller's random state is restored.
#' @return An `ml_lpa_covariates` fit with coefficient matrices, priors,
#'   posteriors, likelihood, information criteria and start diagnostics.
#' @examples
#' set.seed(1)
#' d <- data.frame(group = rep(1:20, each = 10), z = rnorm(200))
#' d$y <- rnorm(200, ifelse(runif(200) < plogis(d$z), -3, 3))
#' fit <- fit_ml_lpa_covariates(d, "y", "group", 2, 1,
#'                             profile_covariates = "z", n_starts = 2, seed = 1)
#' fit$profile_coefficients
#' @export
fit_ml_lpa_covariates <- function(data, indicators, cluster, n_profiles,
                                  n_group_classes = 2L,
                                  profile_covariates = character(),
                                  group_covariates = character(),
                                  variance_model = c("varying", "equal"),
                                  n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                                  min_variance = 1e-6, seed = NULL) {
  stopifnot(is.data.frame(data), is.character(indicators), is.character(cluster),
            is.character(profile_covariates), is.character(group_covariates),
            !anyDuplicated(profile_covariates), !anyDuplicated(group_covariates),
            all(c(profile_covariates, group_covariates) %in% names(data)),
            is.numeric(n_starts), length(n_starts) == 1L,
            is.finite(n_starts), n_starts >= 1, n_starts == as.integer(n_starts))
  variance_model <- match.arg(variance_model)
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
  .ml_lpa_cov_check_covariates(data, profile_covariates, group_covariates,
                               indicators, cluster)
  # The base fit validates indicators/model sizes and supplies an initial mode.
  base <- fit_ml_lpa(data, indicators, cluster, n_profiles, n_group_classes,
                     variance_model, n_starts = 1L, max_iter = max_iter,
                     tol = tol, min_variance = min_variance)
  group_index <- base$group_index
  first_rows <- match(seq_len(base$n_groups), group_index)
  if (length(group_covariates) && any(vapply(group_covariates, function(name) {
    any(data[[name]] != data[[name]][first_rows][group_index])
  }, logical(1)))) stop("Every group_covariate must be constant within each cluster.")
  if (n_profiles == 1L && length(profile_covariates)) {
    stop("Profile covariates require at least two profiles.")
  }
  if (n_group_classes == 1L && length(group_covariates)) {
    stop("Group covariates require at least two group classes.")
  }
  designs <- .ml_lpa_cov_designs(data, indicators, profile_covariates,
                                 group_covariates, first_rows, n_group_classes)
  x <- designs$x
  center <- designs$center
  control <- list(variance_model = variance_model, min_variance = min_variance,
                  max_iter = max_iter, tol = tol,
                  n_profile_covariates = length(profile_covariates),
                  n_group_covariates = length(group_covariates))
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch(.ml_lpa_cov_start(start_index, x, group_index, n_profiles,
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
  result <- .ml_lpa_cov_assemble(
    best = attempts[[best_index]], best_index = best_index, starts = starts,
    designs = designs, base = base, data = data, indicators = indicators,
    cluster = cluster, profile_covariates = profile_covariates,
    group_covariates = group_covariates, n_profiles = n_profiles,
    n_group_classes = n_group_classes, group_index = group_index,
    variance_model = variance_model, min_variance = min_variance,
    call = match.call())
  if (any(!is.finite(starts$log_likelihood))) warning("Some covariate starts failed; inspect $starts.")
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
print.ml_lpa_covariates <- function(x, ...) {
  stopifnot(inherits(x, "ml_lpa_covariates"))
  cat(sprintf("Multilevel LPA with covariates: %d profiles, %d group classes\n",
              x$n_profiles, x$n_group_classes))
  cat(sprintf("Log likelihood %.6f; AIC %.3f; BIC (groups) %.3f; converged %s\n",
              x$log_likelihood, x$aic, x$bic, x$converged))
  invisible(x)
}

#' Extract a covariate LPA log likelihood
#' @param object A covariate LPA fit.
#' @param ... Reserved.
#' @return A logLik object using independent groups for nobs.
#' @examples
#' # logLik(fit)
#' @export
logLik.ml_lpa_covariates <- function(object, ...) {
  stopifnot(inherits(object, "ml_lpa_covariates"))
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
nobs.ml_lpa_covariates <- function(object, ...) {
  stopifnot(inherits(object, "ml_lpa_covariates"))
  object$n_groups
}

#' Validate covariate columns
#'
#' Covariates enter the membership regressions directly, so a non-finite or
#' non-numeric column would propagate silently into the logits.
#'
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.ml_lpa_cov_check_covariates <- function(data, profile_covariates,
                                         group_covariates, indicators, cluster) {
  predictors <- unique(c(profile_covariates, group_covariates))
  if (length(predictors) && !all(vapply(data[predictors], function(column) {
    is.numeric(column) && is.null(dim(column)) && all(is.finite(column))
  }, logical(1)))) {
    stop("Covariates must be finite numeric columns without missing values.")
  }
  if (any(predictors %in% c(indicators, cluster))) {
    stop("Covariates must be distinct from indicators and the cluster identifier.")
  }
  invisible(NULL)
}

#' Build the centered indicator matrix and both membership designs
#'
#' The profile design repeats one group-class indicator block per group class,
#' so that stacking it gives the weighted multinomial regression its rows. The
#' group design takes one row per cluster, at its first occurrence.
#'
#' @return A list with `x`, `center`, `w`, `profile_design` and `stacked_design`.
#' @noRd
.ml_lpa_cov_designs <- function(data, indicators, profile_covariates,
                                group_covariates, first_rows, n_group_classes) {
  x <- as.matrix(data[indicators])
  center <- colMeans(x)
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
       profile_design = profile_design, stacked_design = stacked_design)
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
.ml_lpa_cov_start <- function(start_index, x, group_index, n_profiles,
                              n_group_classes, designs, base, control) {
  parameters <- if (start_index == 1L) {
    list(means = sweep(base$means, 2L, designs$center, "-"),
         variances = base$variances,
         profile_probabilities = base$profile_probabilities,
         group_probabilities = base$group_probabilities)
  } else {
    .ml_lpa_initialize(x, group_index, n_profiles, n_group_classes,
                       control$variance_model, control$min_variance, start_index)
  }
  beta <- rbind(
    log(parameters$profile_probabilities[, seq_len(n_profiles - 1L), drop = FALSE] /
          parameters$profile_probabilities[, n_profiles]),
    matrix(0, control$n_profile_covariates, n_profiles - 1L))
  gamma <- rbind(
    matrix(log(parameters$group_probabilities[seq_len(n_group_classes - 1L)] /
                 parameters$group_probabilities[n_group_classes]), 1L),
    matrix(0, control$n_group_covariates, n_group_classes - 1L))
  expectation <- .ml_lpa_cov_expectation(x, group_index, parameters,
                                         designs$profile_design, designs$w,
                                         beta, gamma)
  history <- expectation$log_likelihood
  iteration <- 0L
  converged <- FALSE
  while (iteration < control$max_iter && !converged) {
    parameters <- .ml_lpa_maximization(x, expectation, control$variance_model,
                                       control$min_variance)
    profile_update <- .ml_lpa_weighted_logits(designs$stacked_design,
                                              do.call(rbind, expectation$joint), beta)
    group_update <- .ml_lpa_weighted_logits(designs$w,
                                            expectation$group_posteriors, gamma)
    beta <- profile_update$coefficients
    gamma <- group_update$coefficients
    updated <- .ml_lpa_cov_expectation(x, group_index, parameters,
                                       designs$profile_design, designs$w,
                                       beta, gamma)
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
#' @return An `ml_lpa_covariates` object.
#' @noRd
.ml_lpa_cov_assemble <- function(best, best_index, starts, designs, base, data,
                                 indicators, cluster, profile_covariates,
                                 group_covariates, n_profiles, n_group_classes,
                                 group_index, variance_model, min_variance, call) {
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
  result$n_parameters <- length(best$beta) + length(best$gamma) +
    n_profiles * n_indicators +
    if (variance_model == "varying") n_profiles * n_indicators else n_indicators
  result$aic <- -2 * result$log_likelihood + 2 * result$n_parameters
  result$bic <- -2 * result$log_likelihood + log(base$n_groups) * result$n_parameters
  result$bic_individual <- -2 * result$log_likelihood +
    log(nrow(data)) * result$n_parameters
  result$n_observations <- nrow(data)
  result$n_groups <- base$n_groups
  result$n_profiles <- n_profiles
  result$n_group_classes <- n_group_classes
  result$group_index <- group_index
  result$group_values <- base$group_values
  result$indicators <- indicators
  result$cluster <- cluster
  result$profile_covariates <- profile_covariates
  result$group_covariates <- group_covariates
  result$variance_model <- variance_model
  result$min_variance <- min_variance
  result$converged <- best$converged
  result$iterations <- best$iterations
  result$starts <- starts
  result$best_start <- best_index
  result$log_likelihood_history <- best$history
  result$subject_profiles <- max.col(result$subject_posteriors, ties.method = "first")
  result$group_classes <- max.col(result$group_posteriors, ties.method = "first")
  result$boundary <- any(result$variances <= min_variance * (1 + 1e-7))
  result$extreme_logits <- any(abs(c(best$beta, best$gamma)) > 20)
  result$call <- call
  structure(result, class = "ml_lpa_covariates")
}
