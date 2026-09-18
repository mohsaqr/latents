#' Row-wise log-sum-exp
#' @param log_values Numeric matrix of log weights.
#' @return Vector of log sums.
#' @noRd
.multilpa_log_sum_exp <- function(log_values) {
  stopifnot(is.matrix(log_values), is.numeric(log_values),
            nrow(log_values) > 0L, ncol(log_values) > 0L,
            !anyNA(log_values))
  row_max <- apply(log_values, 1L, max)
  if (any(!is.finite(row_max))) {
    stop("All component densities vanished, or a density overflowed.")
  }
  row_max + log(rowSums(exp(sweep(log_values, 1L, row_max, "-"))))
}

#' Compute nested posterior probabilities
#' @param x Observation-by-indicator matrix of continuous indicators, optionally
#'   containing NA. It may have zero columns when every indicator is categorical.
#' @param group_index Integer group indices in first-occurrence order.
#' @param parameters List of model parameters.
#' @param codes Optional integer matrix of categorical indicator codes.
#' @return Log likelihood, posteriors, and conditional class probabilities.
#' @noRd
.multilpa_expectation <- function(x, group_index, parameters, codes = NULL) {
  stopifnot(is.matrix(x), is.numeric(x), !any(is.infinite(x)),
            length(group_index) == nrow(x), is.list(parameters))
  n_profiles <- nrow(parameters$means)
  n_types <- length(parameters$group_probabilities)
  gaussian <- if (ncol(x) > 0L && (anyNA(x) || !is.null(parameters$covariances))) {
    .multilpa_gaussian_moments(x, parameters)
  } else NULL
  log_density <- if (!is.null(gaussian)) gaussian$log_density else vapply(seq_len(n_profiles), function(profile) {
    residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
    -0.5 * rowSums(sweep(residuals^2, 2L,
                        parameters$variances[profile, ], "/") +
                    matrix(log(2 * pi) + log(parameters$variances[profile, ]),
                           nrow(x), ncol(x), byrow = TRUE))
  }, numeric(nrow(x)))
  # Categorical indicators are conditionally independent of the continuous ones
  # given the profile, so their log densities simply add.
  if (!is.null(codes)) {
    log_density <- log_density +
      .multilpa_categorical_log_density(codes, parameters$response_probabilities)
  }
  # Each group type supplies a different prior over the same Gaussian profiles.
  conditional <- lapply(seq_len(n_types), function(group_type) {
    log_joint <- sweep(log_density, 2L,
                       log(parameters$profile_probabilities[group_type, ]), "+")
    log_marginal <- .multilpa_log_sum_exp(log_joint)
    list(log_marginal = log_marginal,
         posterior = exp(sweep(log_joint, 1L, log_marginal, "-")))
  })
  group_scores <- matrix(vapply(seq_len(n_types), function(group_type) {
    as.numeric(rowsum(conditional[[group_type]]$log_marginal,
                      group_index, reorder = FALSE)) +
      log(parameters$group_probabilities[group_type])
  }, numeric(max(group_index))), nrow = max(group_index), ncol = n_types)
  group_log_likelihood <- .multilpa_log_sum_exp(group_scores)
  group_posteriors <- exp(sweep(group_scores, 1L, group_log_likelihood, "-"))
  joint <- lapply(seq_len(n_types), function(group_type) {
    conditional[[group_type]]$posterior * group_posteriors[group_index, group_type]
  })
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum(group_log_likelihood)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors))) {
    stop("Non-finite likelihood or posterior probabilities.")
  }
  list(log_likelihood = log_likelihood,
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors,
       subject_posteriors = subject_posteriors, joint = joint,
       gaussian_moments = gaussian$moments)
}

#' Maximize the expected complete-data log likelihood
#' @param x Observation-by-indicator matrix.
#' @param expectation Nested posterior calculation.
#' @param variance_model Either varying or equal across profiles.
#' @param min_variance Lower bound on every indicator variance.
#' @param covariance_model Diagonal or full residual covariance.
#' @return Updated model parameters.
#' @noRd
.multilpa_maximization <- function(x, expectation, variance_model, min_variance,
                                  covariance_model = "diagonal", codes = NULL,
                                  n_categories = NULL, min_probability = 1e-10) {
  stopifnot(is.matrix(x), is.list(expectation),
            variance_model %in% c("varying", "equal"), min_variance > 0)
  weights <- colSums(expectation$subject_posteriors)
  group_weights <- colSums(expectation$group_posteriors)
  if (any(weights <= 0) || any(group_weights <= 0)) {
    stop("A profile or group class has zero effective membership.")
  }
  gaussian <- if (!is.null(expectation$gaussian_moments)) {
    .multilpa_maximize_moments(x, expectation, variance_model, min_variance, covariance_model)
  } else NULL
  if (ncol(x) == 0L) {
    means <- matrix(numeric(0), length(weights), 0L)
    variances <- matrix(numeric(0), length(weights), 0L)
  } else if (is.null(gaussian)) {
    means <- sweep(crossprod(expectation$subject_posteriors, x), 1L, weights, "/")
    variance_sums <- t(matrix(vapply(seq_along(weights), function(profile) {
      residuals <- sweep(x, 2L, means[profile, ], "-")
      colSums(residuals^2 * expectation$subject_posteriors[, profile])
    }, numeric(ncol(x))), nrow = ncol(x), ncol = length(weights)))
    variances <- if (variance_model == "varying") {
      sweep(variance_sums, 1L, weights, "/")
    } else {
      matrix(colSums(variance_sums) / nrow(x),
             length(weights), ncol(x), byrow = TRUE)
    }
    # This is the exact M-step for the stated variance >= min_variance constraint.
    # It is not a hidden ridge, and bound-active estimates are reported to callers.
    variances <- pmax(variances, min_variance)
  } else {
    means <- gaussian$means
    variances <- gaussian$variances
  }
  profile_counts <- t(matrix(vapply(expectation$joint, colSums, numeric(length(weights))),
                             nrow = length(weights), ncol = length(group_weights)))
  profile_probabilities <- profile_counts / rowSums(profile_counts)
  parameters <- list(means = means, variances = variances,
                     profile_probabilities = profile_probabilities,
                     group_probabilities = group_weights / sum(group_weights))
  if (!is.null(gaussian$covariances)) parameters$covariances <- gaussian$covariances
  if (!is.null(codes)) {
    parameters$response_probabilities <- .multilpa_categorical_maximize(
      codes, expectation$subject_posteriors, n_categories, min_probability)
  }
  if (any(!is.finite(unlist(parameters, use.names = FALSE)))) {
    stop("Non-finite parameters in the M-step; check indicator scales.")
  }
  parameters
}

#' Initialize a two-level Gaussian mixture
#' @param x Observation-by-indicator matrix.
#' @param group_index Integer group indices.
#' @param n_profiles Number of individual profiles.
#' @param n_types Number of group classes.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param start_index Restart number.
#' @param covariance_model Diagonal or full residual covariance.
#' @return Strictly positive initial mixture probabilities and Gaussian parameters.
#' @noRd
.multilpa_initialize <- function(x, group_index, n_profiles, n_types,
                               variance_model, min_variance, start_index,
                               covariance_model = "diagonal", codes = NULL,
                               n_categories = NULL, min_probability = 1e-10) {
  stopifnot(is.matrix(x), length(group_index) == nrow(x),
            n_profiles >= 1L, n_types >= 1L, start_index >= 1L,
            variance_model %in% c("varying", "equal"), min_variance > 0)
  if (ncol(x) > 0L && anyNA(x)) {
    initial_means <- colMeans(x, na.rm = TRUE)
    x <- matrix(vapply(seq_len(ncol(x)), function(indicator) {
      value <- x[, indicator]
      value[is.na(value)] <- initial_means[indicator]
      value
    }, numeric(nrow(x))), nrow(x), ncol(x))
  }
  overall_means <- colMeans(x)
  overall_variances <- colMeans(sweep(x, 2L, overall_means, "-")^2)
  if (n_profiles == 1L) {
    means <- matrix(overall_means, 1L)
    assignments <- rep(1L, nrow(x))
  } else {
    # Standardizing only for initialization stops a large-scale indicator from
    # determining all initial clusters; fitting remains in the original units.
    standardized <- sweep(sweep(x, 2L, overall_means, "-"),
                          2L, sqrt(overall_variances), "/")
    # Categorical indicators join the clustering as centered dummy columns, so
    # a mixed-mode model is initialized from all of its indicators at once.
    clustering_input <- if (is.null(codes)) standardized else
      cbind(standardized, .multilpa_categorical_design(codes, n_categories))
    initial_fit <- stats::kmeans(clustering_input, centers = n_profiles,
                                iter.max = 100L, algorithm = "Lloyd")
    if (!is.null(initial_fit$ifault) && initial_fit$ifault != 0L) {
      stop("Initial k-means did not converge.")
    }
    assignments <- initial_fit$cluster
    centers <- initial_fit$centers[, seq_len(ncol(x)), drop = FALSE]
    means <- sweep(sweep(centers, 2L, sqrt(overall_variances), "*"),
                   2L, overall_means, "+")
  }
  variances <- matrix(pmax(overall_variances, min_variance),
                      n_profiles, ncol(x), byrow = TRUE)
  indicator <- vapply(seq_len(n_profiles), function(profile) {
    as.numeric(assignments == profile)
  }, numeric(nrow(x)))
  group_counts <- rowsum(indicator, group_index, reorder = FALSE)
  # Pseudocounts are used only to initialize, never to estimate the final model.
  group_composition <- (group_counts + 0.5) /
    (rowSums(group_counts) + 0.5 * n_profiles)
  if (n_types == 1L) {
    profile_probabilities <- matrix(colMeans(indicator), 1L)
  } else if (start_index %% 2L == 1L &&
             nrow(unique(group_composition)) >= n_types) {
    group_fit <- stats::kmeans(group_composition, centers = n_types,
                              iter.max = 100L, algorithm = "Lloyd")
    if (!is.null(group_fit$ifault) && group_fit$ifault != 0L) {
      stop("Initial group k-means did not converge.")
    }
    profile_probabilities <- group_fit$centers
  } else {
    profile_probabilities <- matrix(stats::runif(n_types * n_profiles, 0.1, 1),
                                    n_types, n_profiles)
  }
  profile_probabilities <- profile_probabilities / rowSums(profile_probabilities)
  result <- list(means = means, variances = variances,
       profile_probabilities = profile_probabilities,
       group_probabilities = rep(1 / n_types, n_types))
  if (!is.null(codes)) {
    result$response_probabilities <- .multilpa_categorical_initialize(
      codes, assignments, n_profiles, n_categories, min_probability)
  }
  if (ncol(x) > 0L && covariance_model == "full") {
    residuals <- sweep(x, 2L, overall_means, "-")
    covariance <- .multilpa_bound_covariance(crossprod(residuals) / nrow(x), min_variance)
    result$covariances <- array(rep(covariance, n_profiles), c(ncol(x), ncol(x), n_profiles))
    result$variances <- matrix(diag(covariance), n_profiles, ncol(x), byrow = TRUE)
  }
  result
}

#' Validate user-supplied initial parameters
#' @param start Initial parameter list.
#' @param n_profiles Number of profiles.
#' @param n_types Number of group classes.
#' @param n_indicators Number of indicators.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param covariance_model Diagonal or full residual covariance.
#' @return A validated copy without dimension names.
#' @noRd
.multilpa_validate_start <- function(start, n_profiles, n_types, n_indicators,
                                   variance_model, min_variance, covariance_model = "diagonal") {
  stopifnot(is.list(start), n_profiles >= 1L, n_types >= 1L,
            n_indicators >= 1L, min_variance > 0)
  covariances <- NULL
  if (covariance_model == "full") {
    covariances <- start$covariances
    if (!is.numeric(covariances) ||
        !identical(as.integer(dim(covariances)), c(as.integer(n_indicators),
          as.integer(n_indicators), as.integer(n_profiles))) || any(!is.finite(covariances))) {
      stop("start$covariances must be a finite indicators x indicators x profiles array.")
    }
    invisible(lapply(seq_len(n_profiles), function(profile) {
      covariance <- matrix(covariances[, , profile], n_indicators, n_indicators)
      if (max(abs(covariance - t(covariance))) > 1e-12 ||
          min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) < min_variance * (1 - 1e-8)) {
        stop("Starting covariances must be symmetric with eigenvalues at least min_variance.")
      }
      if (variance_model == "equal" && max(abs(covariance - covariances[, , 1L])) > 1e-12) {
        stop("Starting covariances must be equal across profiles.")
      }
    }))
    variances <- t(matrix(vapply(seq_len(n_profiles), function(profile) {
      diag(matrix(covariances[, , profile], n_indicators, n_indicators))
    }, numeric(n_indicators)), n_indicators, n_profiles))
    if (!is.null(start$variances) && !isTRUE(all.equal(unname(start$variances), variances, tolerance = 1e-12))) {
      stop("Starting variances must agree with covariance diagonals.")
    }
    start$variances <- variances
    start$covariances <- NULL
  }
  required <- c("means", "variances", "profile_probabilities", "group_probabilities")
  if (!setequal(names(start), required) || anyDuplicated(names(start))) {
    stop("start must contain exactly means, variances, profile_probabilities, and group_probabilities.")
  }
  dimensions <- list(means = c(n_profiles, n_indicators),
                     variances = c(n_profiles, n_indicators),
                     profile_probabilities = c(n_types, n_profiles))
  invisible(lapply(names(dimensions), function(field) {
    value <- start[[field]]
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(as.integer(dim(value)), as.integer(dimensions[[field]])) ||
        any(!is.finite(value))) {
      stop(sprintf("start$%s must be a finite numeric matrix with dimensions %s.",
                   field, paste(dimensions[[field]], collapse = " x ")))
    }
  }))
  if (any(start$variances < min_variance)) {
    stop("start$variances must be at least min_variance.")
  }
  if (variance_model == "equal" &&
      any(abs(sweep(start$variances, 2L, start$variances[1L, ], "-")) > 1e-12)) {
    stop("start$variances must be equal across profiles for variance_model = 'equal'.")
  }
  group_probabilities <- start$group_probabilities
  if (!is.numeric(group_probabilities) || !is.null(dim(group_probabilities)) ||
      length(group_probabilities) != n_types || any(!is.finite(group_probabilities)) ||
      any(group_probabilities <= 0) || abs(sum(group_probabilities) - 1) > 1e-8 ||
      any(start$profile_probabilities <= 0) ||
      any(abs(rowSums(start$profile_probabilities) - 1) > 1e-8)) {
    stop("Starting probabilities must be strictly positive and sum to one (within each profile-probability row).")
  }
  # Normalize only roundoff admitted by the validation tolerance.
  result <- list(means = unname(start$means), variances = unname(start$variances),
       profile_probabilities = unname(start$profile_probabilities /
                                       rowSums(start$profile_probabilities)),
       group_probabilities = unname(group_probabilities / sum(group_probabilities)))
  if (!is.null(covariances)) result$covariances <- unname(covariances)
  result
}

#' Run one nested EM optimization
#' @param x Observation-by-indicator matrix.
#' @param group_index Integer group indices.
#' @param parameters Starting model parameters.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param max_iter Maximum M-steps.
#' @param tol Relative likelihood convergence tolerance.
#' @param covariance_model Diagonal or full residual covariance.
#' @return Final parameters, posterior calculation, and convergence history.
#' @noRd
.multilpa_em <- function(x, group_index, parameters, variance_model,
                        min_variance, max_iter, tol, covariance_model = "diagonal",
                        codes = NULL, n_categories = NULL, min_probability = 1e-10) {
  stopifnot(is.matrix(x), is.list(parameters), max_iter >= 1L, tol > 0,
            length(group_index) == nrow(x), min_variance > 0,
            variance_model %in% c("varying", "equal"))
  expectation <- .multilpa_expectation(x, group_index, parameters, codes)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  while (iteration < max_iter && !converged) {
    parameters <- .multilpa_maximization(x, expectation, variance_model, min_variance,
                                       covariance_model, codes, n_categories,
                                       min_probability)
    updated <- .multilpa_expectation(x, group_index, parameters, codes)
    improvement <- updated$log_likelihood - expectation$log_likelihood
    if (improvement < -1e-10 * (1 + abs(expectation$log_likelihood))) {
      stop("EM likelihood decreased beyond numerical roundoff.")
    }
    converged <- abs(improvement) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated$log_likelihood)
    expectation <- updated
  }
  list(parameters = parameters, expectation = expectation,
       converged = converged, iterations = iteration, history = history)
}

#' Fit a two-level latent profile model
#'
#' Fits individual Gaussian profiles nested within observed groups. A discrete
#' latent group class determines the profile proportions. Profile means and
#' variances are shared across group classes (measurement invariance), and
#' indicators are independent conditional on individual profile membership by
#' default. Full residual covariance and observed-data maximum likelihood for
#' missing indicators are available.
#'
#' @param data A data frame containing indicators and a group identifier.
#' @param indicators Unique character vector of continuous indicator column names.
#' @param group Name of the observed group identifier column. Character,
#'   factor, or numeric identifiers are supported; missing identifiers are not.
#' @param n_profiles Positive integer number of individual profiles.
#' @param n_group_classes Positive integer number of latent group classes.
#'   Use one for an ordinary, pooled LPA.
#' @param variance_model Either `"varying"` (profile-specific variances or
#'   covariance matrices) or `"equal"` (shared across profiles).
#' @param n_starts Positive integer number of EM starts. When `start` is supplied,
#'   it supplies the first start; remaining starts are random initializations.
#' @param max_iter Positive integer maximum number of EM updates per start.
#' @param tol Positive relative log-likelihood tolerance. Convergence requires
#'   absolute change no greater than `tol * (1 + abs(previous log likelihood))`.
#' @param min_variance Positive lower bound on each variance, or each covariance
#'   eigenvalue for full covariance, in squared input units. This defines a constrained maximum-likelihood problem. Bound-active
#'   estimates are explicitly reported and generate a warning.
#' @param seed Optional nonnegative integer random seed. With a supplied seed,
#'   the caller's random-number state is restored on exit.
#' @param start Optional list of `means` and `variances` (profiles by indicators),
#'   `profile_probabilities` (group classes by profiles), and
#'   `group_probabilities` (vector). Starting probabilities must be positive.
#'   For full covariance, supply `covariances` (indicators by indicators by
#'   profiles); `variances` may be omitted or must match their diagonals.
#' @param categorical Character vector naming indicators to treat as
#'   categorical. Each is modelled by unrestricted, profile-specific response
#'   probabilities over its observed categories, which is the latent class
#'   measurement model. Binary, ordinal and unordered indicators are all handled
#'   by the same unrestricted parameterization; numeric, integer, logical,
#'   character and factor columns are accepted. Indicators not named here stay
#'   Gaussian, so naming a subset fits a mixed-mode model.
#' @param min_probability Positive lower bound on every categorical response
#'   probability, defining a constrained maximum-likelihood problem in the same
#'   way `min_variance` does for Gaussian indicators.
#' @param missing `"error"` rejects missing indicators; `"fiml"` maximizes the
#'   observed-data likelihood under an ignorable missingness mechanism (MAR).
#'   Missing indicators are integrated out, not filled in for likelihood fitting.
#' @param covariance_model `"diagonal"` assumes conditional independence;
#'   `"full"` estimates within-profile residual covariances.
#' @return An `multilpa` object containing `means`, `variances`, optional
#'   `covariances` (indicators by indicators by profiles),
#'   `profile_probabilities`, `group_probabilities`, posterior matrices,
#'   classifications, log likelihood, information criteria, restart diagnostics,
#'   and convergence history. Individual rows retain their input order; groups
#'   retain first-occurrence order. `bic` and `bic_groups` use the observed group
#'   count; `bic_individual` uses the individual count. These are alternative
#'   conventions, not interchangeable criteria. `boundary` identifies variance
#'   bounds; `small_classes` flags effective memberships below one. No standard
#'   errors, likelihood-ratio tests, or guarantees of global optimality are given.
#' @details Infinite and constant observed indicators are rejected. Each
#'   indicator must have at least two distinct observed values. No rows are
#'   silently dropped. In FIML mode, fully missing individuals contribute no
#'   direct measurement likelihood but receive posterior probabilities from
#'   their group's information. Counts and information-criterion penalties
#'   retain these individuals and groups. Missing patterns can prevent parameter
#'   identification; no general identification guarantee is made. Initialization
#'   alone uses indicator-mean filling. EM uses conditional Gaussian sufficient
#'   statistics and observed marginal densities. More than one group class requires more than one profile
#'   and at least one group with multiple individuals; these checks are necessary
#'   but do not establish identification. The highest finite likelihood across
#'   starts is returned, even if that start did not converge; inspect `converged`
#'   and `starts`. Profile and group-class labels are arbitrary.
#' @references Vermunt, J. K. (2003). Multilevel latent class models.
#'   Sociological Methodology, 33, 213--239.
#'   doi:10.1111/j.0081-1750.2003.t01-1-00131.x.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(10), each = 12),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                  n_profiles = 1, n_group_classes = 1, n_starts = 1,
#'                  seed = 42)
#' summary(fit)
#' @export
#' @importFrom stats setNames
multilpa <- function(data, indicators, group, n_profiles,
                       n_group_classes = 2L, variance_model = c("varying", "equal"),
                       n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                       min_variance = 1e-6, seed = NULL, start = NULL,
                       missing = c("error", "fiml"),
                       covariance_model = c("diagonal", "full"),
                       categorical = character(), min_probability = 1e-10) {
  stopifnot(is.data.frame(data), is.character(indicators), is.character(group),
            "`categorical` must be a character vector of indicator names" =
              is.character(categorical) && !anyNA(categorical),
            "`min_probability` must be a single number in (0, 1)" =
              is.numeric(min_probability) && length(min_probability) == 1L &&
              is.finite(min_probability) && min_probability > 0 &&
              min_probability < 1)
  call <- match.call()
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  missing <- match.arg(missing)
  .multilpa_check_arguments(data, indicators, group, n_profiles, n_group_classes,
                          n_starts, max_iter, tol, min_variance, min_probability,
                          seed, categorical)
  measurement <- .multilpa_prepare_indicators(data, indicators, categorical,
                                            missing, min_probability)
  continuous <- measurement$continuous
  indicator_frame <- measurement$frame
  encoded <- measurement$encoded
  codes <- measurement$codes
  n_categories <- measurement$n_categories
  x <- measurement$x
  groups <- .multilpa_prepare_groups(data[[group]])
  group_values <- groups$values
  group_index <- groups$index
  group_ids <- groups$ids
  n_groups <- groups$n
  group_sizes <- groups$sizes
  distinct_rows <- if (is.null(codes)) nrow(unique(x)) else
    nrow(unique(cbind(x, codes)))
  if (n_profiles > distinct_rows) {
    stop("n_profiles cannot exceed the number of distinct observed indicator rows.")
  }
  if (n_group_classes > n_groups) stop("n_group_classes cannot exceed the number of groups.")
  if (n_group_classes > 1L && (n_profiles == 1L || all(group_sizes == 1L))) {
    stop("Multiple group classes are not identifiable with one profile or only singleton groups.")
  }
  if (!is.null(start)) {
    if (!is.null(codes)) {
      stop(errorCondition("Supplying `start` is not yet supported for categorical indicators.",
                          class = "multilpa_unsupported_start", call = NULL))
    }
    start <- .multilpa_validate_start(start, n_profiles, n_group_classes,
                                    ncol(x), variance_model, min_variance, covariance_model)
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  # Center to reduce cancellation in weighted means when indicators have large
  # offsets. This is a translation, with no density Jacobian or scale change.
  centers <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  x <- sweep(x, 2L, centers, "-")
  if (any(!is.finite(x[!is.na(x)]^2))) stop("Indicator scales overflow squared residuals; rescale the data.")
  if (!is.null(start)) start$means <- sweep(start$means, 2L, centers, "-")
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch({
      initial <- if (start_index == 1L && !is.null(start)) start else {
        .multilpa_initialize(x, group_index, n_profiles, n_group_classes,
                           variance_model, min_variance, start_index,
                           covariance_model, codes, n_categories, min_probability)
      }
      .multilpa_em(x, group_index, initial, variance_model, min_variance, max_iter,
                 tol, covariance_model, codes, n_categories, min_probability)
    }, error = function(error) list(error = conditionMessage(error)))
  })
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  if (!any(valid)) {
    stop(sprintf("All %d starts failed: %s", n_starts,
                 paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")))
  }
  scores <- vapply(attempts, function(attempt) {
    if (is.null(attempt$error)) attempt$expectation$log_likelihood else -Inf
  }, numeric(1))
  best_start <- which.max(scores)
  best <- attempts[[best_start]]
  parameters <- best$parameters
  parameters$means <- sweep(parameters$means, 2L, centers, "+")
  profile_names <- paste0("profile_", seq_len(n_profiles))
  type_names <- paste0("group_class_", seq_len(n_group_classes))
  parameters <- .multilpa_label_parameters(parameters, profile_names, continuous,
                                         categorical, encoded$levels)
  dimnames(parameters$profile_probabilities) <- list(type_names, profile_names)
  names(parameters$group_probabilities) <- type_names
  subject_posteriors <- best$expectation$subject_posteriors
  group_posteriors <- best$expectation$group_posteriors
  dimnames(subject_posteriors) <- list(rownames(data), profile_names)
  dimnames(group_posteriors) <- list(group_ids, type_names)
  n_parameters <- .multilpa_count_parameters(n_profiles, n_group_classes, ncol(x),
                                           n_categories, variance_model,
                                           covariance_model)
  log_likelihood <- best$expectation$log_likelihood
  starts <- do.call(rbind, lapply(seq_along(attempts), function(start_index) {
    attempt <- attempts[[start_index]]
    data.frame(start = start_index, log_likelihood = scores[start_index],
               converged = if (valid[start_index]) attempt$converged else FALSE,
               iterations = if (valid[start_index]) attempt$iterations else NA_integer_,
               error = if (valid[start_index]) NA_character_ else attempt$error,
               boundary = if (valid[start_index])
                 .multilpa_covariance_boundary(attempt$parameters, min_variance) else NA)
  }))
  boundary <- .multilpa_covariance_boundary(parameters, min_variance)
  small_classes <- any(colSums(subject_posteriors) < 1) || any(colSums(group_posteriors) < 1)
  result <- c(parameters, list(
    call = call, indicators = indicators, continuous = continuous,
    categorical = categorical,
    categorical_levels = encoded$levels, min_probability = min_probability,
    indicator_data = as.matrix(indicator_frame),
    categorical_data = codes,
    group = group, group_ids = group_ids,
    group_values = group_values, group_index = group_index,
    group_sizes = setNames(group_sizes, group_ids),
    n_observations = nrow(x), n_groups = n_groups, n_profiles = as.integer(n_profiles),
    n_group_classes = as.integer(n_group_classes), variance_model = variance_model,
    covariance_model = covariance_model, missing = missing,
    n_observed_by_indicator = setNames(
      c(colSums(!is.na(x)), if (is.null(codes)) NULL else colSums(!is.na(codes))),
      c(continuous, categorical)),
    min_variance = min_variance, standard_deviations = sqrt(parameters$variances),
    measurement_model = if (is.null(codes)) "gaussian" else
      if (ncol(x) == 0L) "categorical" else "mixed",
    subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
    subject_profiles = max.col(subject_posteriors, ties.method = "first"),
    group_classes = max.col(group_posteriors, ties.method = "first"),
    log_likelihood = log_likelihood,
    group_log_likelihood = setNames(best$expectation$group_log_likelihood, group_ids),
    n_parameters = n_parameters, aic = -2 * log_likelihood + 2 * n_parameters,
    bic = -2 * log_likelihood + log(n_groups) * n_parameters,
    bic_groups = -2 * log_likelihood + log(n_groups) * n_parameters,
    bic_individual = -2 * log_likelihood + log(nrow(x)) * n_parameters,
    converged = best$converged, iterations = best$iterations,
    log_likelihood_history = best$history, starts = starts, best_start = best_start,
    n_failed_starts = sum(!valid), boundary = boundary, small_classes = small_classes,
    effective_profile_counts = colSums(subject_posteriors),
    effective_group_counts = colSums(group_posteriors),
    n_best_replicated = sum(valid & abs(scores - log_likelihood) <=
                             1e-6 * (1 + abs(log_likelihood))),
    replication_tolerance = 1e-6 * (1 + abs(log_likelihood))))
  class(result) <- "multilpa"
  if (any(!valid)) warning(sprintf("%d of %d starts failed; inspect $starts$error.", sum(!valid), n_starts), call. = FALSE)
  if (!best$converged) warning("The best start did not converge; increase max_iter and inspect starts.", call. = FALSE)
  if (boundary) warning(if (covariance_model == "full")
    "A covariance eigenvalue reached min_variance; this is a bound-active constrained fit." else
    "A variance reached min_variance; this is a bound-active constrained fit.", call. = FALSE)
  if (small_classes) warning("A profile or group class has effective membership below one.", call. = FALSE)
  result
}

#' Validate the scalar arguments of a fit
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_check_arguments <- function(data, indicators, group, n_profiles,
                                    n_group_classes, n_starts, max_iter, tol,
                                    min_variance, min_probability, seed,
                                    categorical) {
  if (nrow(data) < 2L || length(indicators) < 1L || anyNA(indicators) ||
      anyDuplicated(indicators) || anyDuplicated(names(data)) ||
      length(group) != 1L || is.na(group) || group %in% indicators ||
      !all(c(indicators, group) %in% names(data))) {
    stop("Supply at least two rows, unique existing indicators, and one distinct group column.")
  }
  counts <- list(n_profiles = n_profiles, n_group_classes = n_group_classes,
                 n_starts = n_starts, max_iter = max_iter)
  invisible(lapply(names(counts), function(field) {
    value <- counts[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(sprintf("%s must be a positive integer.", field))
    }
  }))
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      !is.numeric(min_variance) || length(min_variance) != 1L ||
      !is.finite(min_variance) || min_variance <= 0) {
    stop("tol and min_variance must be finite positive numbers.")
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) ||
      seed > .Machine$integer.max)) stop("seed must be a nonnegative integer or NULL.")
  if (anyDuplicated(categorical) || !all(categorical %in% indicators)) {
    stop(errorCondition("`categorical` must name distinct indicators listed in `indicators`.",
                        class = "multilpa_bad_categorical", call = NULL))
  }
  invisible(NULL)
}

#' Split indicators into Gaussian and categorical blocks
#'
#' Returns the continuous indicator matrix, which may have zero columns, and the
#' encoded categorical codes, which may be `NULL`. Every measurement-side
#' contract is checked here so the fitting path can assume clean input.
#'
#' @return A list with `continuous`, `frame`, `x`, `encoded`, `codes` and
#'   `n_categories`.
#' @noRd
.multilpa_prepare_indicators <- function(data, indicators, categorical, missing,
                                       min_probability) {
  continuous <- setdiff(indicators, categorical)
  encoded <- if (length(categorical) > 0L) {
    .multilpa_encode_categorical(data[, categorical, drop = FALSE])
  } else NULL
  codes <- encoded$codes
  n_categories <- encoded$n_categories
  if (!is.null(codes) && min_probability * max(n_categories) >= 1) {
    stop(errorCondition(sprintf(
      "`min_probability` of %g leaves no room for an indicator with %d categories.",
      min_probability, max(n_categories)),
      class = "multilpa_bad_categorical", call = NULL))
  }
  if (!is.null(codes) && missing == "error" && anyNA(codes)) {
    stop("Indicators contain missing or non-finite values.")
  }
  frame <- data[, continuous, drop = FALSE]
  if (!all(vapply(frame, is.numeric, logical(1))) ||
      any(vapply(frame, function(value) !is.null(dim(value)), logical(1)))) {
    stop("Every indicator must be a numeric vector; factors are not continuous indicators.")
  }
  x <- matrix(as.matrix(frame), nrow = nrow(data), ncol = length(continuous),
              dimnames = list(NULL, continuous))
  if (any(is.infinite(x)) || any(is.nan(x)) || (missing == "error" && anyNA(x))) {
    stop("Indicators contain missing or non-finite values.")
  }
  if (any(colSums(!is.na(x)) == 0L)) stop("Every indicator must have observed values.")
  if (any(vapply(frame, function(value) length(unique(value[!is.na(value)])) < 2L,
                 logical(1)))) stop("Constant indicators cannot identify Gaussian profiles.")
  list(continuous = continuous, frame = frame, x = x, encoded = encoded,
       codes = codes, n_categories = n_categories)
}

#' Index the observed groups in first-occurrence order
#' @return A list with `values`, `index`, `ids`, `n` and `sizes`.
#' @noRd
.multilpa_prepare_groups <- function(raw_groups) {
  if (!(is.character(raw_groups) || is.factor(raw_groups) || is.numeric(raw_groups)) ||
      !is.null(dim(raw_groups)) || anyNA(raw_groups) ||
      (is.numeric(raw_groups) && any(!is.finite(raw_groups)))) {
    stop("group must contain nonmissing, finite numeric, factor, or character identifiers.")
  }
  values <- unique(raw_groups)
  index <- match(raw_groups, values)
  ids <- make.unique(as.character(values))
  list(values = values, index = index, ids = ids, n = length(ids),
       sizes = tabulate(index, nbins = length(ids)))
}

#' Attach readable dimension names to the fitted parameters
#' @return The parameter list, with profiles, indicators and categories named.
#' @noRd
.multilpa_label_parameters <- function(parameters, profile_names, continuous,
                                     categorical, levels) {
  dimnames(parameters$means) <- dimnames(parameters$variances) <-
    list(profile_names, continuous)
  if (!is.null(parameters$covariances)) {
    dimnames(parameters$covariances) <- list(continuous, continuous, profile_names)
  }
  if (!is.null(parameters$response_probabilities)) {
    parameters$response_probabilities <- stats::setNames(
      lapply(seq_along(categorical), function(index) {
        block <- parameters$response_probabilities[[index]]
        dimnames(block) <- list(profile_names, levels[[categorical[index]]])
        block
      }), categorical)
  }
  parameters
}

#' Count the free parameters of a fitted model
#'
#' Mixing proportions contribute `(H - 1) + H (K - 1)`. Gaussian indicators add
#' means plus variances or covariances, shared or profile-specific. Categorical
#' measurement is always profile-specific, so `variance_model` has no
#' categorical analogue and never reduces that block.
#'
#' @return A single number of free parameters.
#' @noRd
.multilpa_count_parameters <- function(n_profiles, n_group_classes, n_continuous,
                                     n_categories, variance_model,
                                     covariance_model) {
  covariance_parameters <- if (covariance_model == "full") {
    n_continuous * (n_continuous + 1) / 2
  } else n_continuous
  gaussian <- if (n_continuous == 0L) 0L else {
    n_profiles * n_continuous +
      if (variance_model == "varying") n_profiles * covariance_parameters else
        covariance_parameters
  }
  categorical <- if (is.null(n_categories)) 0L else
    .multilpa_categorical_parameters(n_profiles, n_categories)
  (n_group_classes - 1L) + n_group_classes * (n_profiles - 1L) +
    gaussian + categorical
}
