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
#' @param vars Character vector of indicator names.
#' @param id Name of the group identifier column.
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
#' @return An object of class `multilpa_random_intercept`. Read it with the
#' verbs that describe it rather than by reaching into it: [as.data.frame()]
#' gives the measurement model, the individual posteriors or the posterior
#' group intercepts; [summary()] gives the model-level fit summary, the
#' integration diagnostics and the restart diagnostics; [plot()] draws the
#' measurement model, the group intercepts, or the case-level entropy and
#' posterior panels; and [coef()] returns the free parameters as a named
#' vector.
#' The default BIC uses independent groups; `bic_individual` uses people.
#' Standard errors are not provided: [parameter_inference()], [vcov()] and
#' [confint()] raise a `multilpa_no_inference` condition on this class rather
#' than return an interval the package cannot compute.
#' A near-zero random variance is a boundary.
#' For numerical search, intercept SD is restricted to data RMS times
#' `exp(c(-16, 5))`, and profile logits relative to the final profile to
#' `[-25, 25]`. Fits near those search bounds are flagged as boundaries.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(10), each = 4),
#'                            score = rnorm(40))
#' fit <- fit_random_intercept(example_data, "score", "group",
#'                             n_profiles = 1, n_starts = 1, seed = 1)
#' print(fit)
#' @export
fit_random_intercept <- function(data, vars, id, n_profiles,
    variance_model = c("varying", "equal"), n_starts = 5L, max_iter = 1000L,
    tol = 1e-8, min_variance = 1e-6, quadrature_nodes = 61L,
    quadrature_check_nodes = 121L, quadrature_tolerance = 1e-3, seed = NULL) {
  stopifnot(is.data.frame(data), is.character(vars), length(vars) >= 1L,
    !anyDuplicated(vars), is.character(id), length(id) == 1L,
    all(c(vars, id) %in% names(data)), !id %in% vars,
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
  stopifnot(all(vapply(data[vars], is.numeric, logical(1))))
  x <- as.matrix(data[vars])
  stopifnot(nrow(x) >= 2L, all(is.finite(x)),
            !anyNA(data[[id]]),
            is.numeric(data[[id]]) || is.character(data[[id]]) || is.factor(data[[id]]))
  if (is.numeric(data[[id]])) stopifnot(all(is.finite(data[[id]])))
  group_values <- unique(data[[id]])
  group_index <- match(data[[id]], group_values)
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
  ## Kept in input units so the shared measurement panel and the tidy
  ## accessors can standardize against the observed indicators.
  indicator_data <- x
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
  if (!converged) {
    warning(warningCondition(
      paste0("Best random-intercept start did not converge: ", best$message),
      class = "multilpa_unconverged", call = NULL))
  }
  if (any(!is.finite(values))) {
    warning(warningCondition(
      "Some random-intercept starts failed; summary() reports every start.",
      class = "multilpa_failed_starts"))
  }
  if (discrepancy > quadrature_tolerance) {
    warning(warningCondition(
      "Quadrature likelihood check failed; increase quadrature_nodes and quadrature_check_nodes before interpreting this fit.",
      class = "multilpa_quadrature_check", call = NULL))
  }
  boundary_flags <- c(residual_variance = any(parameters$variances <= min_variance * (1 + 1e-5)),
    random_variance_near_zero = parameters$random_sd <= scale * 1e-5,
    random_sd_upper = parameters$random_sd >= scale * exp(5) * (1 - 1e-5),
    profile_logit = k > 1L && any(abs(best$par[k * d + n_var + 1L + seq_len(k - 1L)]) >= 25 - 1e-5))
  boundary <- any(boundary_flags)
  if (boundary) {
    warning(warningCondition(
      "Random-intercept fit is near a variance or numerical search boundary.",
      class = "multilpa_boundary", call = NULL))
  }
  parameters$means <- sweep(parameters$means, 2L, center, "+")
  dimnames(parameters$means) <- dimnames(parameters$variances) <- list(paste0("profile_", seq_len(k)), vars)
  q <- length(best$par)
  result <- c(parameters, estimates, list(call = match.call(), n_profiles = k,
    n_observations = nrow(x), n_groups = length(group_values), group_values = group_values,
    group_index = group_index, group_sizes = group_sizes, vars = vars,
    continuous = vars, indicator_data = indicator_data, center = center,
    effective_profile_counts = colSums(estimates$subject_posteriors),
    ## Every other family stores the modal assignment beside the posteriors, and
    ## The assignments table and `descriptives(by = "profile")` both document
    ## that they
    ## work here too. Without it the first died with an unclassed `cbind()` error.
    subject_profiles = max.col(estimates$subject_posteriors,
                               ties.method = "first"),
    id = id, variance_model = variance_model, n_parameters = q,
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

#' Print a continuous group random-intercept fit
#' @param x A fitted `multilpa_random_intercept` model.
#' @param rows How many rows of the printed table to show before truncating.
#' @param ... Reserved for compatibility with `print()`.
#' @return The input model, invisibly. Called for the side effect of printing
#'   the profile count and sample sizes, the log likelihood and the fitted
#'   random-intercept standard deviation, and the integration diagnostics.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' random_intercept <- fit_random_intercept(example_data,
#'                                          c("score_a", "score_b"), "group",
#'                                          n_profiles = 2, n_starts = 2,
#'                                          seed = 1)
#' print(random_intercept)
#' @export
print.multilpa_random_intercept <- function(x, rows = 20L, ...) {
  stopifnot(inherits(x, "multilpa_random_intercept"))
  cat(sprintf("Random-intercept LPA: %d profiles, %d people, %d groups\n", x$n_profiles, x$n_observations, x$n_groups))
  cat(sprintf("Log likelihood: %.6f; random-intercept SD: %.6f\n", x$log_likelihood, x$random_sd))
  cat(sprintf("Integration: %s; likelihood check difference: %.3g\n", x$integration, x$quadrature_log_likelihood_difference))
  .multilpa_print_primary(x, rows = rows)
  invisible(x)
}

#' Summarize a continuous group random-intercept fit
#'
#' Collects the model-level fit, the measurement model, the posterior group
#' intercepts, the integration diagnostics and the restart diagnostics into one
#' object, so that none of them has to be read out of the fit by hand.
#'
#' @param object A fitted `multilpa_random_intercept` model.
#' @param ... Reserved for compatibility with `summary()`.
#' @return An object of class `summary_multilpa_random_intercept`, with a
#'   `print` method and an [as.data.frame()] accessor. `as.data.frame()`
#'   returns the one-row model summary by default; `what = "profiles"`,
#'   `"random_intercepts"` and `"starts"` return the measurement model, the
#'   posterior group intercepts and the optimizer starts.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' fit <- fit_random_intercept(example_data, c("score_a", "score_b"), "group",
#'                             n_profiles = 2, n_starts = 2, seed = 1)
#' summary(fit)
#' get_results(fit, what = "random_intercepts")
#' @export
summary.multilpa_random_intercept <- function(object, ...) {
  stopifnot("`object` must be a fitted `multilpa_random_intercept` model" =
              inherits(object, "multilpa_random_intercept"))
  # One builder for the fit and for its summary, so `get_results(x, "model")`
  # and the printed header cannot describe the same fit differently.
  model <- .multilpa_intercept_fit_frame(object)
  result <- list(
    model = model,
    profiles = get_results(object, "profiles"),
    random_intercepts = get_results(object, "random_intercepts"),
    starts = object$starts,
    profile_probabilities = object$profile_probabilities,
    effective_profile_counts = object$effective_profile_counts,
    boundary_flags = object$boundary_flags,
    call = object$call)
  # Every table the fit can produce, built once here, so `get_results()` on the
  # summary serves the same tables the fit would and `print()` can show them
  # all without recomputing anything.
  result$tables <- get_results(object, "all")
  class(result) <- "summary_multilpa_random_intercept"
  result
}

#' Print a random-intercept LPA summary
#' @param x A `summary_multilpa_random_intercept` object.
#' @param digits Number of printed significant digits.
#' @param rows How many rows of each table to print. A longer table is shown
#'   to that depth, with its remaining row count and the `get_results()` call that
#'   returns it whole.
#' @param ... Passed to the underlying `data.frame` printing.
#' @return The summary, invisibly. Called for the side effect of printing the
#'   model line, the measurement model, the profile probabilities and effective
#'   memberships, the likelihood and information criteria, the integration
#'   diagnostics, and the restart diagnostics.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' random_intercept <- fit_random_intercept(example_data,
#'                                          c("score_a", "score_b"), "group",
#'                                          n_profiles = 2, n_starts = 2,
#'                                          seed = 1)
#' print(summary(random_intercept), digits = 3)
#' @export
print.summary_multilpa_random_intercept <- function(x, digits = 4L, rows = 10L,
                                                    ...) {
  stopifnot("`x` must be a `summary_multilpa_random_intercept` object" =
              inherits(x, "summary_multilpa_random_intercept"))
  .multilpa_check_print_arguments(digits, rows)
  model <- x$model
  cat(sprintf("Random-intercept LPA: %d profiles, %d people, %d groups\n",
              model$n_profiles, model$n_observations, model$n_groups))
  cat(sprintf("Parameters: %d; converged: %s; group random-intercept SD: %.6f\n",
              model$n_parameters, model$converged, model$random_intercept_sd))
  cat(sprintf("Log likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              model$log_likelihood, model$aic, model$bic_groups,
              model$bic_individual))
  cat(sprintf("Integration: %s with %d nodes; likelihood check difference %.3g (%s).\n",
              model$integration, model$quadrature_nodes,
              model$quadrature_log_likelihood_difference,
              if (isTRUE(model$quadrature_check_passed)) "passed" else "FAILED"))
  if (!model$converged) cat("WARNING: the best start did not converge.\n")
  if (model$boundary) {
    cat(sprintf("WARNING: at a boundary (%s).\n",
                paste(names(x$boundary_flags)[x$boundary_flags], collapse = ", ")))
  }
  cat("Standard errors are not available for this model family.\n")
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' Coerce a random-intercept model summary to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' A summary carries every table the object it describes can produce, and
#' [get_results()] names them.
#'
#' @param x An object of class `summary_multilpa_random_intercept`.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `multilpa_bad_argument`
#'   naming it, rather than being dropped.
#' @return A base `data.frame`: one row per profile and continuous indicator.
#' @seealso [get_results()] for every other table this summary holds.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- fit_random_intercept(
#'   example_data, c("score_a", "score_b"), "school", n_profiles = 2,
#'   n_starts = 1, seed = 1
#' )
#' as.data.frame(summary(fit))
#' @export
as.data.frame.summary_multilpa_random_intercept <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_random_intercept`" = inherits(x, "summary_multilpa_random_intercept"))
  .multilpa_coerce(x, row.names, list(...))
}

#' Plot a continuous group random-intercept fit
#'
#' Draws the measurement model exactly as [plot.multilpa()] draws it, the
#' posterior group intercepts with their posterior standard deviations, or the
#' case-level classification diagnostics. The intercept panel is the one thing
#' this model family has that a discrete group-class model does not, so it is
#' available here and nowhere else.
#'
#' @param x A fitted `multilpa_random_intercept` model.
#' @param what `"profiles"` (the default) draws the Gaussian measurement model,
#'   one line per profile across the indicators. `"random_intercepts"` draws
#'   one interval per group: its posterior mean intercept plus and minus one
#'   posterior standard deviation, with the groups ordered by that mean.
#'   `"entropy"` and `"posteriors"` draw the case-level classification
#'   diagnostics, as in [plot.multilpa()]; both need more than one profile and
#'   raise `multilpa_nothing_to_plot` on a single-profile fit.
#' @param scale For `what = "profiles"`, `"raw"` or `"standardized"`, as in
#'   [plot.multilpa()].
#' @param labels `TRUE` prints a direct label at the right end of each series.
#' @param main,subtitle Panel title and secondary line, or `NULL` for defaults.
#' @param palette,symbols,linetypes Series aesthetics, or `NULL` for the
#'   package defaults.
#' @param style A list of visual constants, as built by `.multilpa_style()`.
#' @param ... Further named visual constants, merged into `style`.
#' @return The fitted model, invisibly. Called for the side effect of drawing.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' fit <- fit_random_intercept(example_data, c("score_a", "score_b"), "group",
#'                             n_profiles = 2, n_starts = 2, seed = 1)
#' plot(fit)
#' plot(fit, what = "random_intercepts")
#' plot(fit, what = "entropy")
#' @export
plot.multilpa_random_intercept <- function(x, what = c("profiles",
                                                       "random_intercepts",
                                                       "entropy", "posteriors",
                                                       "all"),
                                           scale = c("raw", "standardized"),
                                           labels = TRUE, main = NULL,
                                           subtitle = NULL, palette = NULL,
                                           symbols = NULL, linetypes = NULL,
                                           style = .multilpa_style(), ...) {
  stopifnot("`x` must be a fitted `multilpa_random_intercept` model" =
              inherits(x, "multilpa_random_intercept"),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels))
  what <- match.arg(what)
  if (identical(what, "all")) {
    return(.multilpa_plot_every_view(x, match.call(), parent.frame()))
  }
  scale <- match.arg(scale)
  style <- utils::modifyList(style, list(...))
  previous <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA)
  if (identical(what, "profiles")) {
    ## The shared measurement panel reads `n_group_classes` only to build a
    ## default subtitle. This family has a continuous group effect instead of
    ## discrete classes, so the adapter supplies that field and this method
    ## always replaces the subtitle it would produce.
    panel_input <- x
    panel_input$n_group_classes <- 1L
    .multilpa_plot_profiles(panel_input, scale, labels, main,
      if (is.null(subtitle)) sprintf(
        "%d profiles; continuous group random intercept (SD %.3g); %s scale",
        x$n_profiles, x$random_sd,
        if (identical(scale, "standardized")) "standardized" else "input")
      else subtitle,
      palette, symbols, linetypes, style)
  } else if (identical(what, "random_intercepts")) {
    .multilpa_plot_random_intercepts(x, main, subtitle, palette, style)
  } else {
    ## The case-level diagnostics read only the posteriors and the effective
    ## profile counts, which this family carries like every other, so they are
    ## drawn by the shared helper rather than reimplemented here. It refuses a
    ## single-profile fit, where there is nothing to separate.
    .multilpa_plot_case_diagnostic(x, what, main, subtitle, palette, style)
  }
  invisible(x)
}

#' Draw the posterior group intercepts
#' @param x A fitted `multilpa_random_intercept` model.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Series aesthetics, or `NULL` for defaults.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_random_intercepts <- function(x, main, subtitle, palette, style) {
  frame <- get_results(x, "random_intercepts")
  order_by_mean <- order(frame$mean, frame$group)
  centre <- frame$mean[order_by_mean]
  spread <- frame$sd[order_by_mean]
  colours <- if (is.null(palette)) .multilpa_palette(2L) else
    rep(palette, length.out = 2L)
  positions <- seq_along(centre)
  span <- range(c(centre - spread, centre + spread, 0))
  padding <- 0.12 * max(diff(span), .Machine$double.eps)
  graphics::par(mar = style$margins)
  .multilpa_panel(xlim = c(0.4, length(centre) + 0.6),
    ylim = c(span[1L] - padding, span[2L] + padding),
    xlab = "Group, ordered by posterior mean intercept",
    ylab = "Posterior group intercept",
    main = if (is.null(main)) "Posterior group random intercepts" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d groups; estimated intercept SD %.3g; bars are +/- one posterior SD",
      length(centre), x$random_sd) else subtitle,
    x_at = positions, x_labels = as.character(frame$group[order_by_mean]),
    style = style)
  graphics::abline(h = 0, col = style$muted_colour, lwd = 1, lty = 3L)
  graphics::segments(positions, centre - spread, positions, centre + spread,
                     col = colours[1L], lwd = style$line_width)
  graphics::points(positions, centre, pch = 21L, bg = colours[1L],
                   col = style$panel_fill, cex = style$point_size, lwd = 1.4)
  invisible(NULL)
}

#' Free parameters of a random-intercept fit
#'
#' @param object A fitted `multilpa_random_intercept` model.
#' @param ... Reserved for compatibility with `coef()`.
#' @return A named numeric vector with one element per free parameter. Names
#'   follow the package-wide `level.parameter.outcome.term` grammar, which is
#'   the tidy decomposition [parameter_inference()] reports written on one
#'   line, so a name reads back into those columns:
#'   `measurement.mean.profile_k.indicator` for every profile mean,
#'   `measurement.variance.profile_k.indicator` (or
#'   `measurement.variance.shared.indicator` when `variance_model = "equal"`,
#'   where one spread is shared across profiles) for every residual variance,
#'   `group.standard_deviation.random_intercept` for the group intercept
#'   standard deviation, and `profile.logit.profile_k` for each profile logit
#'   relative to the final profile. Means, variances and the intercept standard
#'   deviation are in input units; a logit is named as a logit because that is
#'   the scale it is estimated and reported on. It has one element per free
#'   parameter the fit reports, so its length is the parameter count
#'   `summary()` prints. This is the base generic's contract, a plain named
#'   vector; for a tidy table of the measurement model use
#'   `as.data.frame(object)`.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' fit <- fit_random_intercept(example_data, c("score_a", "score_b"), "group",
#'                             n_profiles = 2, n_starts = 2, seed = 1)
#' coef(fit)
#' @export
coef.multilpa_random_intercept <- function(object, ...) {
  stopifnot("`object` must be a fitted `multilpa_random_intercept` model" =
              inherits(object, "multilpa_random_intercept"))
  vars <- object$vars
  n_profiles <- object$n_profiles
  profiles <- paste0("profile_", seq_len(n_profiles))
  probabilities <- object$profile_probabilities
  ## An equal-variance fit estimates one variance per indicator, shared across
  ## profiles; naming every profile's copy would report it several times and
  ## overstate the free parameter count.
  shared_variance <- identical(object$variance_model, "equal")
  values <- c(
    as.vector(t(object$means)),
    if (shared_variance) object$variances[1L, ] else as.vector(t(object$variances)),
    object$random_sd,
    if (n_profiles > 1L)
      log(probabilities[seq_len(n_profiles - 1L)] / probabilities[n_profiles])
    else numeric())
  labels <- rbind(
    data.frame(level = "measurement", parameter = "mean",
               outcome = rep(profiles, each = length(vars)),
               term = rep(vars, times = n_profiles)),
    data.frame(level = "measurement", parameter = "variance",
               outcome = if (shared_variance) "shared" else
                 rep(profiles, each = length(vars)),
               term = if (shared_variance) vars else
                 rep(vars, times = n_profiles)),
    data.frame(level = "group", parameter = "standard_deviation",
               outcome = "random_intercept", term = NA_character_),
    ## A single-profile fit estimates no logit at all, so it contributes no
    ## row rather than a zero-length column that `data.frame()` cannot build.
    if (n_profiles > 1L) data.frame(level = "profile", parameter = "logit",
               outcome = profiles[seq_len(n_profiles - 1L)],
               term = NA_character_) else NULL)
  stopifnot("every free parameter must be labelled" =
              nrow(labels) == length(values))
  ## The shared serialiser is the single point of truth for this grammar, so
  ## the classes of this package cannot drift into separate naming schemes.
  stats::setNames(unname(values), .multilpa_parameter_names(labels))
}

#' Standard errors are not available for a random-intercept fit
#'
#' @param object A fitted `multilpa_random_intercept` model. `vcov()` and
#'   `confint()` are base generics, so their first formal is `object`.
#' @param x The same fitted model, under the name this package's own verbs
#'   use; `parameter_inference()` is documented on this page too.
#' @param ... Ignored.
#' @return Nothing; always raises a `multilpa_no_inference` condition. The
#'   observed-information and sandwich machinery this package uses differentiates
#'   the discrete two-level likelihood analytically. That score does not cover
#'   the quadrature integral over a continuous group intercept, so no standard
#'   error is reported rather than one that is not the model's.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' random_intercept <- fit_random_intercept(example_data,
#'                                          c("score_a", "score_b"), "group",
#'                                          n_profiles = 2, n_starts = 2,
#'                                          seed = 1)
#' # The refusal is catchable by class, not by message text.
#' tryCatch(vcov(random_intercept),
#'          multilpa_no_inference = function(condition) {
#'            "no standard errors for this model family"
#'          })
#' @export
#' @importFrom stats vcov
vcov.multilpa_random_intercept <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_random_intercept"))
  stop(errorCondition(paste(
    "Standard errors are not available for a continuous group random-intercept",
    "model. Read the estimates with coef() and as.data.frame()."),
    class = "multilpa_no_inference", call = NULL))
}

#' @rdname vcov.multilpa_random_intercept
#' @param data Ignored; present for compatibility with the generic.
#' @export
parameter_inference.multilpa_random_intercept <- function(x, data = NULL, ...) {
  vcov(x)
}

#' @rdname vcov.multilpa_random_intercept
#' @param parm Ignored; present for compatibility with the generic.
#' @param level Ignored; present for compatibility with the generic.
#' @export
#' @importFrom stats confint
confint.multilpa_random_intercept <- function(object, parm, level = 0.95, ...) {
  vcov(object)
}

#' Extract a random-intercept LPA log likelihood
#' @param object A fitted `multilpa_random_intercept` model.
#' @param ... Reserved for compatibility with `logLik()`.
#' @return A `logLik` object with parameter count `df` and the number of
#'   observed groups as `nobs`, so `stats::BIC()` uses the group-count BIC.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' random_intercept <- fit_random_intercept(example_data,
#'                                          c("score_a", "score_b"), "group",
#'                                          n_profiles = 2, n_starts = 2,
#'                                          seed = 1)
#' logLik(random_intercept)
#' @export
logLik.multilpa_random_intercept <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_random_intercept"))
  structure(object$log_likelihood, df = object$n_parameters, nobs = object$n_groups, class = "logLik")
}

#' Count independent groups in a random-intercept LPA fit
#' @param object A fitted `multilpa_random_intercept` model.
#' @param ... Reserved for compatibility with `nobs()`.
#' @return A single integer: the number of observed groups, which are the
#'   independent units of this likelihood.
#' @examples
#' set.seed(8)
#' example_data <- data.frame(group = rep(seq_len(12), each = 5),
#'                            score_a = rnorm(60), score_b = rnorm(60))
#' random_intercept <- fit_random_intercept(example_data,
#'                                          c("score_a", "score_b"), "group",
#'                                          n_profiles = 2, n_starts = 2,
#'                                          seed = 1)
#' nobs(random_intercept)
#' @export
nobs.multilpa_random_intercept <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_random_intercept"))
  object$n_groups
}
