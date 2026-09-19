#' Standard normal Gaussian quadrature
#' @param n Number of nodes.
#' @return Nodes and probability weights for a standard normal distribution.
#' @noRd
.ri_quadrature <- function(n) {
  stopifnot(length(n) == 1L, is.numeric(n), is.finite(n), n >= 3L, n == as.integer(n))
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- sqrt(seq_len(n - 1L))
  decomposition <- eigen(jacobi + t(jacobi), symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1L, ]^2)
}

#' Unpack random-intercept model parameters
#' @param theta Unconstrained parameter vector.
#' @param k Number of profiles.
#' @param d Number of indicators.
#' @param variance_model Residual variance specification.
#' @return Model parameter list.
#' @noRd
.ri_unpack <- function(theta, k, d, variance_model) {
  stopifnot(is.numeric(theta), k >= 1L, d >= 1L,
            variance_model %in% c("equal", "varying"))
  n_var <- if (variance_model == "equal") d else k * d
  means <- matrix(theta[seq_len(k * d)], k, d)
  variances <- exp(theta[k * d + seq_len(n_var)])
  variances <- if (variance_model == "equal") matrix(rep(variances, each = k), k, d) else matrix(variances, k, d)
  random_sd <- exp(theta[k * d + n_var + 1L])
  logits <- c(if (k > 1L) theta[k * d + n_var + 1L + seq_len(k - 1L)] else numeric(), 0)
  probabilities <- exp(logits - max(logits))
  list(means = means, variances = variances, random_sd = random_sd,
       profile_probabilities = probabilities / sum(probabilities))
}

#' Evaluate a scalar random-intercept mixture likelihood
#' @param x Complete numeric indicator matrix.
#' @param group_index Contiguous integer group indices.
#' @param parameters Means, variances, random_sd, and profile_probabilities.
#' @param quadrature Standard normal quadrature rule.
#' @param posterior Whether to return conditional profile and intercept summaries.
#' @return Log likelihood and optional posterior summaries.
#' @noRd
.ri_evaluate <- function(x, group_index, parameters, quadrature, posterior = FALSE) {
  stopifnot(is.matrix(x), is.numeric(x), length(group_index) == nrow(x),
            is.list(parameters), is.list(quadrature), is.logical(posterior))
  k <- nrow(parameters$means)
  group_rows <- split(seq_len(nrow(x)), group_index)
  if (k == 1L) {
    residual <- sweep(x, 2L, parameters$means[1L, ], "-")
    inverse_variance <- 1 / parameters$variances[1L, ]
    tau2 <- parameters$random_sd^2
    group_results <- lapply(group_rows, function(rows) {
      precision <- length(rows) * sum(inverse_variance)
      weighted_sum <- sum(sweep(residual[rows, , drop = FALSE], 2L, inverse_variance, "*"))
      quadratic <- sum(sweep(residual[rows, , drop = FALSE]^2, 2L, inverse_variance, "*"))
      log_det <- length(rows) * sum(log(parameters$variances[1L, ])) + log1p(tau2 * precision)
      c(log_likelihood = -0.5 * (length(rows) * ncol(x) * log(2 * pi) + log_det +
          quadratic - tau2 * weighted_sum^2 / (1 + tau2 * precision)),
        random_mean = tau2 * weighted_sum / (1 + tau2 * precision),
        random_variance = tau2 / (1 + tau2 * precision))
    })
    values <- do.call(rbind, group_results)
    result <- list(log_likelihood = sum(values[, "log_likelihood"]),
                   group_log_likelihood = values[, "log_likelihood"])
    if (posterior) {
      result$subject_posteriors <- matrix(1, nrow(x), 1L)
      result$random_intercept_mean <- values[, "random_mean"]
      result$random_intercept_sd <- sqrt(values[, "random_variance"])
    }
    return(result)
  }
  nodes <- parameters$random_sd * quadrature$nodes
  log_profile <- lapply(seq_len(k), function(profile) {
    Reduce(`+`, lapply(seq_len(ncol(x)), function(indicator) {
      centered <- outer(x[, indicator] - parameters$means[profile, indicator], nodes, "-")
      -0.5 * (log(2 * pi * parameters$variances[profile, indicator]) +
                centered^2 / parameters$variances[profile, indicator])
    })) + log(parameters$profile_probabilities[profile])
  })
  maximum <- Reduce(pmax, log_profile)
  density <- maximum + log(Reduce(`+`, lapply(log_profile, function(value) exp(value - maximum))))
  group_density <- rowsum(density, group_index, reorder = FALSE)
  # rowsum's encounter order matches the contiguous indices created by the fitter.
  log_weight <- sweep(group_density, 2L, log(quadrature$weights), "+")
  log_likelihood <- .multilpa_log_sum_exp(log_weight)
  result <- list(log_likelihood = sum(log_likelihood), group_log_likelihood = log_likelihood)
  if (posterior) {
    node_posterior <- exp(log_weight - log_likelihood)
    subject_node <- node_posterior[group_index, , drop = FALSE]
    result$subject_posteriors <- matrix(vapply(log_profile, function(value) {
      rowSums(exp(value - density) * subject_node)
    }, numeric(nrow(x))), nrow(x), k)
    result$random_intercept_mean <- drop(node_posterior %*% nodes)
    result$random_intercept_sd <- sqrt(pmax(0, drop(node_posterior %*% nodes^2) - result$random_intercept_mean^2))
  }
  result
}

#' Fit LPA with one continuous group random intercept
#'
#' Fits a Gaussian profile mixture with a normally distributed group intercept
#' that adds to every indicator with loading fixed at one. Profile probabilities
#' are constant across groups. Residuals are conditionally independent. This is
#' a specific random-intercept model, not a general random-effects interface.
#' A single profile is integrated analytically; mixtures use fixed Gaussian
#' quadrature. A higher-order likelihood check diagnoses integration error but
#' does not refit the model. Increase both node counts if the check fails.
#'
#' @param data Data frame with complete numeric indicators.
#' @param indicators Character vector of indicator names.
#' @param group Name of the group identifier column.
#' @param n_profiles Positive integer number of profiles.
#' @param variance_model Profile-specific (varying) or shared (equal) residual variances.
#' @param n_starts Number of independently initialized optimizations.
#' @param max_iter Maximum optimizer iterations per start.
#' @param tol Relative optimizer convergence tolerance.
#' @param min_variance Lower bound for residual variances.
#' @param quadrature_nodes Number of Gaussian quadrature nodes used for fitting.
#' @param quadrature_check_nodes Larger number of nodes for diagnostic evaluation.
#' @param quadrature_tolerance Maximum acceptable absolute log-likelihood discrepancy.
#' @param seed Optional integer random seed, with caller RNG state restored.
#' @return An object of class `multilpa_random_intercept` containing parameter
#' estimates, posterior profile probabilities and group intercept moments,
#' likelihood, information criteria, optimizer starts, and integration diagnostics.
#' The default BIC uses independent groups; `bic_individual` uses people.
#' Standard errors are not provided. A near-zero random variance is a boundary.
#' For numerical search, intercept SD is restricted to data RMS times
#' `exp(c(-16, 5))`, and profile logits relative to the final profile to
#' `[-25, 25]`. Fits near those search bounds are flagged as boundaries.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(10), each = 4),
#'                            score = rnorm(40))
#' fit <- fit_random_intercept(example_data, "score", "group", 1,
#'                                  n_starts = 1)
#' print(fit)
#' @export
fit_random_intercept <- function(data, indicators, group, n_profiles,
    variance_model = c("varying", "equal"), n_starts = 5L, max_iter = 1000L,
    tol = 1e-8, min_variance = 1e-6, quadrature_nodes = 61L,
    quadrature_check_nodes = 121L, quadrature_tolerance = 1e-3, seed = NULL) {
  stopifnot(is.data.frame(data), is.character(indicators), length(indicators) >= 1L,
    !anyDuplicated(indicators), is.character(group), length(group) == 1L,
    all(c(indicators, group) %in% names(data)), !group %in% indicators,
    is.numeric(n_profiles), length(n_profiles) == 1L, is.finite(n_profiles),
    n_profiles >= 1L, n_profiles == as.integer(n_profiles),
    is.numeric(n_starts), length(n_starts) == 1L, is.finite(n_starts),
    n_starts >= 1L, n_starts == as.integer(n_starts),
    is.numeric(max_iter), length(max_iter) == 1L, is.finite(max_iter),
    max_iter >= 1L, max_iter == as.integer(max_iter),
    is.numeric(tol), length(tol) == 1L, is.finite(tol), tol > 0,
    is.numeric(min_variance), length(min_variance) == 1L, is.finite(min_variance), min_variance > 0,
    is.numeric(quadrature_tolerance), length(quadrature_tolerance) == 1L,
    is.finite(quadrature_tolerance), quadrature_tolerance > 0)
  variance_model <- match.arg(variance_model)
  quadrature <- .ri_quadrature(quadrature_nodes)
  check_quadrature <- .ri_quadrature(quadrature_check_nodes)
  stopifnot(quadrature_check_nodes > quadrature_nodes)
  stopifnot(all(vapply(data[indicators], is.numeric, logical(1))))
  x <- as.matrix(data[indicators])
  stopifnot(nrow(x) >= 2L, all(is.finite(x)),
            !anyNA(data[[group]]),
            is.numeric(data[[group]]) || is.character(data[[group]]) || is.factor(data[[group]]))
  if (is.numeric(data[[group]])) stopifnot(all(is.finite(data[[group]])))
  group_values <- unique(data[[group]])
  group_index <- match(data[[group]], group_values)
  group_sizes <- tabulate(group_index)
  stopifnot(length(group_values) >= 2L, any(group_sizes > 1L),
            n_profiles <= nrow(unique(x)))
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed), seed >= 0,
              seed <= .Machine$integer.max, seed == as.integer(seed))
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  center <- colMeans(x)
  x <- sweep(x, 2L, center, "-")
  k <- as.integer(n_profiles)
  d <- ncol(x)
  n_var <- if (variance_model == "equal") d else k * d
  scale <- sqrt(mean(x^2))
  stopifnot(is.finite(scale), scale > 0)
  objective <- function(theta) {
    stopifnot(is.numeric(theta))
    -.ri_evaluate(x, group_index, .ri_unpack(theta, k, d, variance_model), quadrature)$log_likelihood
  }
  starts <- lapply(seq_len(n_starts), function(index) {
    tryCatch({
      classification <- if (k == 1L) rep(1L, nrow(x)) else
        stats::kmeans(x, centers = k, algorithm = "Lloyd", iter.max = 100L)$cluster
      means <- t(matrix(vapply(seq_len(k), function(profile)
        colMeans(x[classification == profile, , drop = FALSE]), numeric(d)), d, k))
      residual <- x - means[classification, , drop = FALSE]
      variances <- t(matrix(vapply(seq_len(k), function(profile) {
        pmax(colMeans(residual[classification == profile, , drop = FALSE]^2), min_variance * 10)
      }, numeric(d)), d, k))
      if (variance_model == "equal") variances <- pmax(colMeans(residual^2), min_variance * 10)
      group_average <- rowsum(rowMeans(residual), group_index, reorder = FALSE)[, 1L] / group_sizes
      tau <- max(sqrt(max(stats::var(group_average) - mean(variances) * mean(1 / group_sizes) / d, 0)), scale * 0.05)
      if (index > 1L) tau <- tau * exp(stats::rnorm(1L, sd = 0.5))
      probabilities <- tabulate(classification, nbins = k) / nrow(x)
      theta <- c(as.vector(means), log(as.vector(variances)), log(tau),
                 if (k > 1L) log(probabilities[seq_len(k - 1L)] / probabilities[k]) else numeric())
      result <- stats::optim(theta, objective, method = "L-BFGS-B",
        lower = c(rep(-Inf, k * d), rep(log(min_variance), n_var), log(scale) - 16, rep(-25, k - 1L)),
        upper = c(rep(Inf, k * d + n_var), log(scale) + 5, rep(25, k - 1L)),
        control = list(maxit = max_iter, factr = tol / .Machine$double.eps,
                       ndeps = rep(1e-5, length(theta))))
      list(result = result, error = NA_character_)
    }, error = function(error) list(result = NULL, error = conditionMessage(error)))
  })
  values <- vapply(starts, function(start) if (is.null(start$result)) -Inf else -start$result$value, numeric(1))
  if (!any(is.finite(values))) stop("All random-intercept starts failed: ", paste(vapply(starts, `[[`, character(1), "error"), collapse = "; "))
  best_start <- which.max(values)
  best <- starts[[best_start]]$result
  parameters <- .ri_unpack(best$par, k, d, variance_model)
  estimates <- .ri_evaluate(x, group_index, parameters, quadrature, posterior = TRUE)
  check <- .ri_evaluate(x, group_index, parameters, check_quadrature)
  discrepancy <- abs(check$log_likelihood - estimates$log_likelihood)
  converged <- best$convergence == 0L
  if (!converged) warning("Best random-intercept start did not converge: ", best$message, call. = FALSE)
  if (any(!is.finite(values))) warning("Some random-intercept starts failed; inspect $starts.", call. = FALSE)
  if (discrepancy > quadrature_tolerance) warning("Quadrature likelihood check failed; increase quadrature_nodes and quadrature_check_nodes before interpreting this fit.", call. = FALSE)
  boundary_flags <- c(residual_variance = any(parameters$variances <= min_variance * (1 + 1e-5)),
    random_variance_near_zero = parameters$random_sd <= scale * 1e-5,
    random_sd_upper = parameters$random_sd >= scale * exp(5) * (1 - 1e-5),
    profile_logit = k > 1L && any(abs(best$par[k * d + n_var + 1L + seq_len(k - 1L)]) >= 25 - 1e-5))
  boundary <- any(boundary_flags)
  if (boundary) warning("Random-intercept fit is near a variance or numerical search boundary.", call. = FALSE)
  parameters$means <- sweep(parameters$means, 2L, center, "+")
  dimnames(parameters$means) <- dimnames(parameters$variances) <- list(paste0("profile_", seq_len(k)), indicators)
  q <- length(best$par)
  result <- c(parameters, estimates, list(call = match.call(), n_profiles = k,
    n_observations = nrow(x), n_groups = length(group_values), group_values = group_values,
    group_index = group_index, group_sizes = group_sizes, indicators = indicators,
    group = group, variance_model = variance_model, n_parameters = q,
    aic = -2 * estimates$log_likelihood + 2 * q,
    bic = -2 * estimates$log_likelihood + log(length(group_values)) * q,
    bic_individual = -2 * estimates$log_likelihood + log(nrow(x)) * q,
    converged = converged, boundary = boundary, boundary_flags = boundary_flags,
    optimization_bounds = list(residual_variance_min = min_variance,
      random_sd = scale * exp(c(-16, 5)), profile_logit = c(-25, 25)),
    optimizer_message = best$message, best_start = best_start,
    starts = data.frame(start = seq_len(n_starts), log_likelihood = values,
      converged = vapply(starts, function(start) !is.null(start$result) && start$result$convergence == 0L, logical(1)),
      error = vapply(starts, `[[`, character(1), "error")),
    integration = if (k == 1L) "analytic" else "fixed Gaussian quadrature",
    quadrature_nodes = quadrature_nodes, quadrature_check_nodes = quadrature_check_nodes,
    quadrature_log_likelihood_difference = discrepancy,
    quadrature_check_passed = discrepancy <= quadrature_tolerance))
  class(result) <- "multilpa_random_intercept"
  result
}

#' @export
print.multilpa_random_intercept <- function(x, ...) {
  stopifnot(inherits(x, "multilpa_random_intercept"))
  cat(sprintf("Random-intercept LPA: %d profiles, %d people, %d groups\n", x$n_profiles, x$n_observations, x$n_groups))
  cat(sprintf("Log likelihood: %.6f; random-intercept SD: %.6f\n", x$log_likelihood, x$random_sd))
  cat(sprintf("Integration: %s; likelihood check difference: %.3g\n", x$integration, x$quadrature_log_likelihood_difference))
  invisible(x)
}

#' @export
logLik.multilpa_random_intercept <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_random_intercept"))
  structure(object$log_likelihood, df = object$n_parameters, nobs = object$n_groups, class = "logLik")
}

#' @export
nobs.multilpa_random_intercept <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_random_intercept"))
  object$n_groups
}
