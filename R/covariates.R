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
  scores <- sweep(scores, 1L, .multilpa_row_max(scores), "-")
  sweep(scores, 1L, log(rowSums(exp(scores))), "-")
}

#' Fit weighted multinomial logits
#'
#' Newton-Raphson with the exact Hessian and step halving: the weighted
#' multinomial log likelihood is concave, so each accepted step cannot lower it
#' and the score falls quadratically to the tolerance. (BFGS, used before,
#' stalled above the tolerance on about half of the steps of an ordinary fit,
#' at small coefficients and differently on each platform.) The search runs on
#' unit-free coordinates, not on the covariates as supplied: a covariate held
#' in units a few orders of magnitude away from the rest would otherwise make
#' the Hessian badly scaled. Because the logit is the same model under any linear
#' rescaling of a predictor, the entire search is moved onto columns divided by
#' their root-mean-square (see `.multilpa_design_scale()`) and the coefficients
#' are divided back afterwards, so the caller receives estimates in the units
#' the covariates arrived in and the search itself no longer depends on those
#' units. Convergence is then judged by the score in those unit-free
#' coordinates rather than by the optimizer's own return code.
#'
#' @param design Design matrix.
#' @param counts Expected category counts in each row.
#' @param initial Initial coefficient matrix, in the units of `design`.
#' @param score_tol Relative tolerance for the unit-free score: the step counts
#'   as solved when the largest absolute scaled score is at most
#'   `score_tol * (1 + sum(counts))`, which is the total weight the score is
#'   accumulated over. The default is tight enough that the inner solve is not
#'   what stops the outer EM iteration; a step that cannot reach it is one whose
#'   coefficient is drifting towards infinity, and that is reported rather than
#'   passed off as a maximum.
#' @param max_newton Most Newton steps before the step is declared unconverged.
#' @return A list with `coefficients` (in the units of `design`), `converged`,
#'   and `scaled_score`, the largest absolute unit-free score at the solution.
#' @noRd
.multilpa_weighted_logits <- function(design, counts, initial,
                                      score_tol = 1e-8, max_newton = 100L) {
  stopifnot(
    "`design` must be a numeric matrix" =
      is.matrix(design) && is.numeric(design),
    "`counts` must be a numeric matrix" =
      is.matrix(counts) && is.numeric(counts),
    "`initial` must be a numeric matrix" =
      is.matrix(initial) && is.numeric(initial),
    "`design` and `counts` must have the same number of rows" =
      nrow(design) == nrow(counts),
    "`counts` must be non-negative" = all(counts >= 0),
    "`initial` must have one row per design column" =
      nrow(initial) == ncol(design),
    "`initial` must have one column per non-reference category" =
      ncol(initial) == ncol(counts) - 1L,
    "`score_tol` must be a single positive number" =
      is.numeric(score_tol) && length(score_tol) == 1L &&
      is.finite(score_tol) && score_tol > 0,
    "`max_newton` must be a single positive whole number" =
      is.numeric(max_newton) && length(max_newton) == 1L &&
      is.finite(max_newton) && max_newton >= 1 &&
      max_newton == as.integer(max_newton))
  ## A non-finite design or count would reach the optimizer as a silent `NaN`
  ## objective, which BFGS reports as a completed search.
  if (!all(is.finite(design)) || !all(is.finite(counts))) {
    stop(errorCondition(paste(
      "The membership design and the expected class counts must be finite;",
      "a non-finite entry reaches the multinomial optimizer as an undefined",
      "objective."), class = "latents_bad_data", call = NULL))
  }
  if (ncol(counts) == 1L) {
    return(list(coefficients = initial, converged = TRUE, scaled_score = 0))
  }
  scale <- .multilpa_design_scale(design)
  scaled_design <- sweep(design, 2L, scale, "/")
  n_free <- ncol(counts) - 1L
  unpack <- function(values) {
    stopifnot(is.numeric(values))
    matrix(values, ncol(design), n_free)
  }
  objective <- function(values) {
    stopifnot(is.numeric(values))
    -sum(counts * .multilpa_log_softmax(scaled_design, unpack(values)))
  }
  gradient <- function(values) {
    stopifnot(is.numeric(values))
    residual <- .multilpa_softmax(scaled_design, unpack(values)) *
      rowSums(counts) - counts
    as.vector(crossprod(scaled_design, residual[, seq_len(n_free), drop = FALSE]))
  }
  start <- as.vector(sweep(initial, 1L, scale, "*"))
  baseline <- objective(start)
  ## The score is a sum over the expected counts, so the tolerance it is judged
  ## against scales with them; the `1 +` keeps an empty step well defined.
  threshold <- score_tol * (1 + sum(counts))
  # Minus the Hessian of the log likelihood, block (a, b) over the free
  # categories: sum_i t_i p_ia (1[a = b] - p_ib) x_i x_i'.
  information <- function(values) {
    probabilities <- .multilpa_softmax(scaled_design, unpack(values))
    totals <- rowSums(counts)
    cells <- expand.grid(a = seq_len(n_free), b = seq_len(n_free))
    blocks <- lapply(seq_len(nrow(cells)), function(cell) {
      a <- cells$a[cell]
      b <- cells$b[cell]
      crossprod(scaled_design, scaled_design *
                  (totals * probabilities[, a] * ((a == b) - probabilities[, b])))
    })
    do.call(rbind, lapply(seq_len(n_free), function(a) {
      do.call(cbind, blocks[cells$a == a])
    }))
  }
  current <- start
  value <- baseline
  score <- max(abs(gradient(current)))
  iteration <- 0L
  # Newton steps refine one estimate in sequence, so there is nothing here to
  # vectorize. A ridge far below the data scale keeps an empty category
  # solvable without moving a determined solution.
  while (iteration < max_newton && score > threshold) {
    hessian <- information(current)
    ridge <- 1e-10 * (1 + max(abs(diag(hessian))))
    factor <- chol(hessian + diag(ridge, nrow(hessian)))
    direction <- -backsolve(factor, forwardsolve(t(factor), gradient(current)))
    ## Never move to a worse point: the M-step has to be monotone for EM's
    ## ascent guarantee to hold.
    step <- 1
    candidate_value <- Inf
    while (step > 1e-10) {
      candidate <- current + step * direction
      candidate_value <- objective(candidate)
      # Accept a step that does not raise the objective beyond rounding.
      if (is.finite(candidate_value) &&
          candidate_value <= value + 1e-12 * (1 + abs(value))) break
      step <- step / 2
    }
    if (!is.finite(candidate_value) ||
        candidate_value > value + 1e-12 * (1 + abs(value))) break
    current <- candidate
    value <- candidate_value
    score <- max(abs(gradient(current)))
    iteration <- iteration + 1L
  }
  if (value > baseline + 1e-8) {
    stop(errorCondition(
      "Multinomial M-step decreased the expected log likelihood.",
      class = "latents_no_converge", call = NULL))
  }
  list(coefficients = sweep(unpack(current), 1L, scale, "/"),
       converged = score <= threshold, scaled_score = score)
}

#' Covariate-dependent nested expectation
#' @param x Centred continuous indicator matrix, `NA` where unobserved.
#' @param group_index Group indices.
#' @param parameters Measurement parameters.
#' @param profile_design List of profile design matrices, one per group class.
#' @param group_design Group-level design matrix.
#' @param beta Profile logit coefficients.
#' @param gamma Group logit coefficients.
#' @return Log likelihood and posterior responsibilities.
#' @noRd
.multilpa_cov_expectation <- function(x, group_index, parameters, profile_design,
                                    group_design, beta, gamma, codes = NULL,
                                    sampling_weights = NULL, extra = NULL) {
  stopifnot(is.matrix(x), !any(is.infinite(x)), is.list(parameters),
            is.list(profile_design), is.matrix(group_design), is.matrix(beta),
            is.matrix(gamma))
  # The measurement model is the covariate-free one; only the mixing weights
  # differ here. Residual covariances and missing indicators need the
  # conditional moments, which the maximization step and the scores then
  # reuse, so they are computed once and carried.
  measurement <- .latents_measurement_log_density(x, parameters, codes, extra,
                                                  constant = "joint")
  result <- c(.multilpa_logit_structure(measurement$log_density, measurement$offset,
                                        group_index, profile_design, group_design,
                                        beta, gamma),
              list(gaussian_moments = measurement$moments))
  # Sampling weights scale the posteriors into weighted counts, exactly as in
  # the covariate-free model, so the logit steps and scores need no change.
  if (is.null(sampling_weights)) result else
    .multilpa_weigh(result, sampling_weights, group_index)
}

#' Fit multilevel LPA with class-membership covariates
#'
#' The estimator behind `multilpa(profile_covariates =, group_covariates =)`.
#' Numeric covariates predict individual-profile and group-class membership via
#' multinomial logits, using the final class as reference. Individual-profile
#' slopes are shared across group classes; profile intercepts differ by group
#' class. The measurement model is the one [multilpa()] fits and stays
#' invariant across group classes. This is one-step maximum likelihood, not
#' regression on assigned classes. Covariates are used in their supplied units.
#' Under `missing = "fiml"` missing indicators are integrated out of the
#' measurement density, exactly as in the covariate-free model; covariates must
#' still be complete.
#'
#' @param data,vars,id,n_profiles,n_group_classes As in [multilpa()].
#' @param profile_covariates,group_covariates Names of numeric predictors of
#'   profile membership, and of group-class membership.
#' @param variance_model,n_starts,max_iter,tol,min_variance As in [multilpa()].
#' @param seed,time,covariance_model,categorical,min_probability As in
#'   [multilpa()].
#' @param call The user-facing call to record on the result, so the object
#'   reports the `multilpa()` the caller wrote rather than this delegation.
#' @return An object of class `multilpa_covariates`.
#' @noRd
.multilpa_fit_covariates <- function(data, vars, id, n_profiles,
                                  n_group_classes = 2L,
                                  profile_covariates = character(),
                                  group_covariates = character(),
                                  variance_model = c("varying", "equal"),
                                  n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                                  min_variance = 1e-6, seed = NULL,
                                  time = NULL,
                                  covariance_model = c("diagonal", "full"),
                                  categorical = character(),
                                  min_probability = 1e-10,
                                  missing = c("error", "fiml"),
                                  profile_slopes = c("shared", "group_class"),
                                  select_start = "likelihood", call = NULL,
                                  weights = NULL, ordinal = character(),
                                  count = character(), count_model = "poisson",
                                  count_dispersion = "varying") {
  stopifnot(is.data.frame(data), is.character(vars), is.character(id),
            is.character(profile_covariates), is.character(group_covariates),
            !anyDuplicated(profile_covariates), !anyDuplicated(group_covariates),
            all(c(profile_covariates, group_covariates) %in% names(data)),
            is.numeric(n_starts), length(n_starts) == 1L,
            is.finite(n_starts), n_starts >= 1, n_starts == as.integer(n_starts))
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  missing <- match.arg(missing)
  profile_slopes <- match.arg(profile_slopes)
  stopifnot(
    "`categorical` must be a character vector of indicator names" =
      is.character(categorical) && !anyNA(categorical),
    "`min_probability` must be a single number in (0, 1)" =
      is.numeric(min_probability) && length(min_probability) == 1L &&
      is.finite(min_probability) && min_probability > 0 && min_probability < 1)
  .multilpa_check_seed(seed)
  .latents_local_seed(seed)
  .multilpa_cov_check_covariates(data, profile_covariates, group_covariates,
                               vars, id)
  # The base fit validates indicators/model sizes and supplies an initial mode.
  base <- multilpa(data, vars, id, n_profiles, n_group_classes,
                     variance_model = variance_model,
                     n_starts = 1L, max_iter = max_iter,
                     tol = tol, min_variance = min_variance,
                     covariance_model = covariance_model,
                     categorical = categorical,
                     min_probability = min_probability, missing = missing,
                     weights = weights, ordinal = ordinal, count = count,
                     count_model = count_model, count_dispersion = count_dispersion)
  group_index <- base$group_index
  sampling_weights <- unname(base$sampling_weights)
  first_rows <- match(seq_len(base$n_groups), group_index)
  if (length(group_covariates) && any(vapply(group_covariates, function(name) {
    any(data[[name]] != data[[name]][first_rows][group_index])
  }, logical(1)))) {
    stop(errorCondition("Every group_covariate must be constant within each group.",
                        class = "latents_bad_data", call = NULL))
  }
  if (n_profiles == 1L && length(profile_covariates)) {
    stop("Profile covariates require at least two profiles.")
  }
  if (n_group_classes == 1L && length(group_covariates)) {
    stop("Group covariates require at least two group classes.")
  }
  designs <- .multilpa_cov_designs(data, vars, profile_covariates,
                                 group_covariates, first_rows, n_group_classes,
                                 categorical, min_probability, missing,
                                 profile_slopes, ordinal, count, count_model,
                                 count_dispersion)
  x <- designs$x
  center <- designs$center
  control <- list(variance_model = variance_model, min_variance = min_variance,
                  max_iter = max_iter, tol = tol,
                  covariance_model = covariance_model,
                  min_probability = min_probability,
                  n_profile_covariates = ncol(designs$stacked_design) -
                    n_group_classes,
                  n_group_covariates = length(group_covariates),
                  sampling_weights = sampling_weights)
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
    stop(errorCondition(sprintf("All covariate starts failed: %s",
                                paste(unique(starts$error), collapse = "; ")),
                        class = "latents_all_starts_failed", call = NULL))
  }
  ## Restarts of a mixture routinely reach the same optimum under different
  ## label permutations, and then differ only in the last bits of the
  ## likelihood. Picking by `which.max()` alone makes the reported labelling
  ## turn on floating-point noise, so the maximum is taken up to a relative
  ## tolerance far below any difference that could mean anything, and the
  ## earliest start attaining it wins. The first start is the one resumed from
  ## the covariate-free fit, so ties resolve towards the reproducible mode.
  ## A converged start is preferred over an unconverged one that reached the
  ## same likelihood, because only the converged one is a maximum;
  ## `select_start = "converged"` goes further and refuses an unconverged start
  ## whenever any start converged.
  best_index <- .multilpa_select_start(starts$log_likelihood, starts$converged,
                                       select_start)
  best <- attempts[[best_index]]
  # Weighted posteriors are counts for the M-step; each unit is classified by
  # its own, unweighted posterior, and the likelihood stays the weighted one.
  weighted_expectation <- best$expectation
  if (!is.null(sampling_weights)) {
    best$expectation <- .multilpa_cov_expectation(
      x, group_index, best$parameters, designs$profile_design, designs$w,
      best$beta, best$gamma, designs$codes, extra = designs$extra)
    best$expectation$log_likelihood <- weighted_expectation$log_likelihood
  }
  result <- .multilpa_cov_assemble(
    best = best, best_index = best_index, starts = starts,
    designs = designs, base = base, data = data, vars = vars,
    id = id, profile_covariates = profile_covariates,
    group_covariates = group_covariates,
    # Stored as integers so the field has the same type it has on every other fit
    # class; a double here made `identical()` on any count derived from it fail.
    n_profiles = as.integer(n_profiles),
    n_group_classes = as.integer(n_group_classes), group_index = group_index,
    variance_model = variance_model, min_variance = min_variance,
    covariance_model = covariance_model, categorical = categorical,
    call = call %||% match.call())
  ## Set here rather than inside the assembler, where the name `time` would
  ## resolve to stats::time instead of this argument.
  result$time <- time
  result$time_values <- .multilpa_time_values(data, time, id, vars)
  ## The floor the response probabilities were held at; inference needs it to
  ## tell an estimate on that boundary from one merely close to it.
  result$min_probability <- min_probability
  ## Read by inference to accept, and integrate over, the same missing values.
  result$missing <- missing
  ## Read by inference to rebuild the same design.
  result$profile_slopes <- profile_slopes
  ## The same effective counts the Gaussian fit carries, so the shared
  ## diagnostics and plot panels need no special case for this class.
  result$effective_profile_counts <- colSums(result$subject_posteriors)
  result$weights <- weights
  if (!is.null(sampling_weights)) {
    result$sampling_weights <- stats::setNames(sampling_weights, base$group_ids)
    result$effective_profile_counts <- stats::setNames(
      colSums(weighted_expectation$subject_posteriors),
      names(result$effective_profile_counts))
  }
  result$effective_group_counts <- colSums(result$group_posteriors)
  if (!is.null(sampling_weights)) {
    result$effective_group_counts <- stats::setNames(
      colSums(weighted_expectation$group_posteriors),
      names(result$effective_group_counts))
  }
  if (any(!is.finite(starts$log_likelihood))) {
    warning(warningCondition(paste(
      "Some covariate starts failed; summary() reports every start."),
      class = "latents_failed_starts"))
  }
  if (!result$converged) {
    warning(warningCondition(paste(
      "The best covariate start had not converged when `max_iter` was reached,",
      "or its membership logits still carry a non-negligible score, so the",
      "returned estimate is not a maximum."),
      class = "latents_unconverged"))
  }
  if (result$boundary) {
    warning(warningCondition("A residual variance reached min_variance.",
                             class = "latents_boundary", call = NULL))
  }
  if (result$extreme_logits) {
    warning(warningCondition(
      "Extreme logit coefficients: inspect scaling, sparse classes and separation.",
      class = "latents_extreme_coefficients", call = NULL))
  }
  .latents_warn_poisson_limit(result$count_dispersion)
  result
}

#' Print a covariate LPA fit
#' @param x A covariate LPA fit.
#' @param rows How many rows of the printed table to show before truncating.
#' @param ... Reserved.
#' @return The model, invisibly. Called for the side effect of printing the
#'   class counts, the covariate counts, the log likelihood with the
#'   information criteria, and the convergence diagnostics.
#' @examples
#' fit <- multilpa(subset(course_engagement, student <= 40),
#'                 c("browse", "lectures", "forum_read"), "student",
#'                 n_profiles = 2, n_group_classes = 1,
#'                 profile_covariates = "previous_grade", n_starts = 2, seed = 1)
#' print(fit)
#' @export
print.multilpa_covariates <- function(x, rows = 20L, ...) {
  stopifnot(inherits(x, "multilpa_covariates"))
  cat(sprintf("Multilevel LPA with covariates: %d profiles, %d group classes\n",
              x$n_profiles, x$n_group_classes))
  .latents_print_extra(x)
  .latents_print_weights(x)
  cat(sprintf("Log likelihood %.6f; AIC %.3f; BIC (groups) %.3f; converged %s\n",
              x$log_likelihood, x$aic, x$bic, x$converged))
  .multilpa_print_primary(x, rows = rows)
  invisible(x)
}

#' Summarize a covariate LPA fit
#'
#' Collects the model-level fit, the measurement model, the membership
#' regressions and the restart diagnostics into one object, so that none of
#' them has to be read out of the fit by hand.
#'
#' @param object A covariate LPA fit from `multilpa(profile_covariates = )`.
#' @param ... Reserved for compatibility with `summary()`.
#' @return An object of class `summary_multilpa_covariates`, with a `print`
#'   method and an [as.data.frame()] accessor. `as.data.frame()` returns the
#'   one-row model summary by default; `what = "profiles"`, `"coefficients"`
#'   and `"starts"` return the measurement model, the membership coefficients
#'   and the restart diagnostics. The membership coefficients carry no standard
#'   errors here; [parameter_inference()] reports those.
#' @examples
#' fit <- multilpa(subset(course_engagement, student <= 40),
#'                 c("browse", "lectures", "forum_read"), "student",
#'                 n_profiles = 2, n_group_classes = 1,
#'                 profile_covariates = "previous_grade", n_starts = 2, seed = 1)
#' summary(fit)
#' get_results(fit, what = "coefficients")
#' @export
summary.multilpa_covariates <- function(object, ...) {
  stopifnot("`object` must be a fitted `multilpa_covariates` model" =
              inherits(object, "multilpa_covariates"))
  # One builder for the fit and for its summary, so `get_results(x, "model")`
  # and the printed header cannot describe the same fit differently.
  model <- .multilpa_covariate_fit_frame(object)
  result <- list(
    model = model,
    profiles = get_results(object, "profiles"),
    coefficients = get_results(object, "coefficients"),
    starts = object$starts,
    effective_profile_counts = object$effective_profile_counts,
    effective_group_counts = object$effective_group_counts,
    call = object$call)
  # Every table the fit can produce, built once here, so `get_results()` on the
  # summary serves the same tables the fit would and `print()` can show them
  # all without recomputing anything.
  result$tables <- get_results(object, "all")
  class(result) <- "summary_multilpa_covariates"
  result
}

#' Print a covariate LPA summary
#' @param x A `summary_multilpa_covariates` object.
#' @param digits Number of printed significant digits.
#' @param rows How many rows of each table to print. A longer table is shown
#'   to that depth, with its remaining row count and the `get_results()` call that
#'   returns it whole.
#' @param ... Passed to the underlying `data.frame` printing.
#' @return The summary, invisibly. Called for the side effect of printing the
#'   model line, the measurement model, the membership regressions, the
#'   effective memberships at both levels, the likelihood and information
#'   criteria, any warnings, and the restart diagnostics.
#' @examples
#' fit <- multilpa(subset(course_engagement, student <= 40),
#'                 c("browse", "lectures", "forum_read"), "student",
#'                 n_profiles = 2, n_group_classes = 1,
#'                 profile_covariates = "previous_grade", n_starts = 2, seed = 1)
#' print(summary(fit), digits = 3)
#' @export
print.summary_multilpa_covariates <- function(x, digits = 4L, rows = 10L, ...) {
  stopifnot("`x` must be a `summary_multilpa_covariates` object" =
              inherits(x, "summary_multilpa_covariates"))
  .multilpa_check_print_arguments(digits, rows)
  model <- x$model
  cat(sprintf("Multilevel LPA with covariates: %d profiles, %d group classes\n",
              model$n_profiles, model$n_group_classes))
  cat(sprintf("Individuals: %d; groups: %d; parameters: %d; converged: %s\n",
              model$n_observations, model$n_groups, model$n_parameters,
              model$converged))
  cat(sprintf("%d profile covariate(s); %d group covariate(s)\n",
              model$n_profile_covariates, model$n_group_covariates))
  cat(sprintf("Log likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              model$log_likelihood, model$aic, model$bic_groups,
              model$bic_individual))
  cat("Membership coefficients carry no standard errors here; parameter_inference() has them.\n")
  if (!model$converged) cat("WARNING: the best start did not converge.\n")
  if (model$boundary) cat("WARNING: a residual variance is at min_variance.\n")
  if (model$extreme_logits) {
    cat("WARNING: extreme logit coefficients; check scaling, sparse classes and separation.\n")
  }
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' Coerce a covariate-model summary to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' A summary carries every table the object it describes can produce, and
#' [get_results()] names them.
#'
#' @param x An object of class `summary_multilpa_covariates`.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `latents_bad_argument`
#'   naming it, rather than being dropped.
#' @return A base `data.frame`: one row per profile and continuous indicator.
#' @seealso [get_results()] for every other table this summary holds.
#' @examples
#' fit <- multilpa(
#'   subset(course_engagement, student <= 30),
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2,
#'   profile_covariates = "sequence", n_starts = 1, seed = 1
#' )
#' as.data.frame(summary(fit))
#' @export
as.data.frame.summary_multilpa_covariates <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa_covariates`" = inherits(x, "summary_multilpa_covariates"))
  .multilpa_coerce(x, row.names, list(...))
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
      class = "latents_no_indicator_data", call = NULL))
  }
  result <- data.frame(row.names = seq_len(object$n_observations))
  result <- cbind(result, as.data.frame(object$indicator_data))
  ## Categorical indicators come back from their codes through the one typed
  ## value stored for each category, so a factor returns as that factor and a
  ## number as that number.
  invisible(lapply(object$categorical %||% character(), function(name) {
    result[[name]] <<- object$categorical_values[[name]][
      object$categorical_data[, name]]
  }))
  extra <- .latents_draw_extra_columns(object$extra_data)
  if (length(extra) > 0L) result[names(extra)] <- extra
  result[[object$id]] <- object$group_values[object$group_index]
  ## The profile design is the group-class indicator block followed by the
  ## profile covariates, in the order they were named.
  ## With slopes by group class the covariates repeat once per class; group
  ## class 1's design holds them unmasked in the first of those blocks.
  profile_block <- object$profile_design[[1L]][,
    object$n_group_classes + seq_along(object$profile_covariates), drop = FALSE]
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
#' @return A `logLik` object carrying the maximized log likelihood, the free
#'   parameter count as `df`, and the number of observed groups as `nobs`, so
#'   `stats::BIC()` uses the group-count BIC.
#' @examples
#' fit <- multilpa(subset(course_engagement, student <= 40),
#'                 c("browse", "lectures", "forum_read"), "student",
#'                 n_profiles = 2, n_group_classes = 1,
#'                 profile_covariates = "previous_grade", n_starts = 2, seed = 1)
#' logLik(fit)
#' @export
logLik.multilpa_covariates <- function(object, ...) {
  stopifnot(inherits(object, "multilpa_covariates"))
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' Count independent groups in a covariate LPA fit
#' @param object A covariate LPA fit.
#' @param ... Reserved.
#' @return A single integer: the number of observed groups, which are the
#'   independent units of this likelihood.
#' @examples
#' fit <- multilpa(subset(course_engagement, student <= 40),
#'                 c("browse", "lectures", "forum_read"), "student",
#'                 n_profiles = 2, n_group_classes = 1,
#'                 profile_covariates = "previous_grade", n_starts = 2, seed = 1)
#' nobs(fit)
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
    stop(errorCondition(paste(
      "Covariates must be finite numeric columns without missing values;",
      "`missing = \"fiml\"` integrates out missing indicators, not missing",
      "covariates."), class = "latents_bad_data", call = NULL))
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
#' With `profile_slopes = "shared"` the covariates follow the indicators once,
#' so every group class uses the same slopes. With `"group_class"` they follow
#' once per group class, and group class `h`'s design carries them only in its
#' own block: the slopes are then the covariates' interactions with the group
#' classes. Either way the first block after the indicators holds the
#' covariates themselves in group class 1's design.
#'
#' @return A list with `x`, `center`, `w`, `profile_design` and `stacked_design`.
#' @noRd
.multilpa_cov_designs <- function(data, vars, profile_covariates,
                                group_covariates, first_rows, n_group_classes,
                                categorical = character(),
                                min_probability = 1e-10, missing = "error",
                                profile_slopes = "shared", ordinal = character(),
                                count = character(), count_model = "poisson",
                                count_dispersion = "varying") {
  measurement <- .multilpa_prepare_indicators(data, vars, categorical,
                                              missing, min_probability,
                                              other = c(ordinal, count))
  x <- measurement$x
  # The centre is only a location shift that is added back to the reported
  # means, so the observed values' mean serves; every indicator has some,
  # which the indicator check above guarantees.
  center <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  z <- as.matrix(data[profile_covariates])
  w <- cbind(`(Intercept)` = 1,
             as.matrix(data[first_rows, group_covariates, drop = FALSE]))
  profile_design <- lapply(seq_len(n_group_classes), function(group_class) {
    indicators <- matrix(as.numeric(seq_len(n_group_classes) == group_class),
                         nrow(data), n_group_classes, byrow = TRUE)
    slopes <- if (identical(profile_slopes, "group_class")) {
      do.call(cbind, lapply(seq_len(n_group_classes), function(block) {
        z * as.numeric(block == group_class)
      }))
    } else z
    cbind(indicators, slopes)
  })
  stacked_design <- do.call(rbind, profile_design)
  if (qr(stacked_design)$rank != ncol(stacked_design) || qr(w)$rank != ncol(w)) {
    stop("Covariate design is rank deficient; remove constant or collinear predictors.")
  }
  list(x = sweep(x, 2L, center, "-"), center = center, w = w,
       profile_design = profile_design, stacked_design = stacked_design,
       profile_slopes = profile_slopes,
       continuous = measurement$continuous, codes = measurement$codes,
       n_categories = measurement$n_categories,
       categorical_levels = measurement$encoded$levels,
       categorical_values = .multilpa_categorical_values(
         data, categorical, measurement$encoded),
       extra = .latents_set_count_model(
         .latents_prepare_extra(data, ordinal, count, missing), count_model,
         count_dispersion))
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
    if (!is.null(base$extra_data)) {
      resumed$ordinal_intercepts <- if (is.null(base$ordinal_intercepts)) NULL else
        unname(lapply(base$ordinal_intercepts, unname))
      resumed$ordinal_locations <- unname(base$ordinal_locations)
      resumed$count_means <- unname(base$count_means)
      resumed$count_dispersion <- unname(base$count_dispersion)
    }
    resumed
  } else {
    # The first start is Ward's hierarchical clustering, as in multilpa().
    .multilpa_initialize(x, group_index, n_profiles, n_group_classes,
                       control$variance_model, control$min_variance, start_index,
                       control$covariance_model, designs$codes,
                       designs$n_categories, control$min_probability,
                       hierarchical = start_index == 1L, extra = designs$extra)
  }
  beta <- rbind(
    log(parameters$profile_probabilities[, seq_len(n_profiles - 1L), drop = FALSE] /
          parameters$profile_probabilities[, n_profiles]),
    matrix(0, control$n_profile_covariates, n_profiles - 1L))
  gamma <- rbind(
    matrix(log(parameters$group_probabilities[seq_len(n_group_classes - 1L)] /
                 parameters$group_probabilities[n_group_classes]), 1L),
    matrix(0, control$n_group_covariates, n_group_classes - 1L))
  evaluate <- function(state) {
    .multilpa_cov_expectation(x, group_index, state$parameters,
                              designs$profile_design, designs$w, state$beta,
                              state$gamma, designs$codes,
                              control$sampling_weights, designs$extra)
  }
  maximize <- function(state, expectation) {
    parameters <- .multilpa_maximization(x, expectation, control$variance_model,
                                         control$min_variance,
                                         control$covariance_model, designs$codes,
                                         designs$n_categories,
                                         control$min_probability,
                                         previous = state$parameters,
                                         extra = designs$extra)
    profile_update <- .multilpa_weighted_logits(designs$stacked_design,
                                                do.call(rbind, expectation$joint),
                                                state$beta)
    group_update <- .multilpa_weighted_logits(designs$w,
                                              expectation$group_posteriors,
                                              state$gamma)
    # Both the likelihood and the inner logit steps must settle, so a stalled
    # regression cannot be reported as a converged fit.
    structure(list(parameters = parameters, beta = profile_update$coefficients,
                   gamma = group_update$coefficients),
              solved = profile_update$converged && group_update$converged)
  }
  ## Sweeps spent waiting only on the membership logits are bounded separately
  ## from `max_iter`: a coefficient drifting to infinity leaves that step
  ## unsolved for ever while the observed likelihood has stopped moving, and
  ## sweeping `max_iter` times against it costs a great deal and buys nothing.
  fit <- .latents_em(
    list(parameters = parameters, beta = beta, gamma = gamma), evaluate, maximize,
    max_iter = control$max_iter, tol = control$tol,
    settings = .latents_em_settings(
      decrease_tolerance = 1e-9,
      decrease_condition = "Covariate EM decreased the observed log likelihood.",
      max_stalled = 20L))
  list(parameters = fit$state$parameters, expectation = fit$expectation,
       beta = fit$state$beta, gamma = fit$state$gamma, iterations = fit$iterations,
       converged = fit$converged, likelihood_settled = fit$settled,
       history = fit$history, error = NA_character_)
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
  slope_terms <- if (identical(designs$profile_slopes, "group_class")) {
    paste0(rep(profile_covariates, times = n_group_classes), ":group_class_",
           rep(seq_len(n_group_classes), each = length(profile_covariates)))
  } else profile_covariates
  dimnames(result$profile_coefficients) <- list(
    # With one group class the class intercept is simply the intercept.
    c(if (n_group_classes == 1L) "(Intercept)" else
        paste0("group_class_", seq_len(n_group_classes)), slope_terms),
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
    measurement_parameters + .latents_extra_n_parameters(designs$extra, n_profiles)
  result$aic <- -2 * result$log_likelihood + 2 * result$n_parameters
  result$bic <- -2 * result$log_likelihood + log(base$n_groups) * result$n_parameters
  result$bic_individual <- -2 * result$log_likelihood +
    log(nrow(data)) * result$n_parameters
  result$n_observations <- nrow(data)
  result$n_informative <- base$n_informative
  result$indicator_data <- as.matrix(data[designs$continuous])
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
  result$continuous_types <- vapply(data[designs$continuous], typeof,
                                    character(1))
  result$categorical_levels <- designs$categorical_levels
  result$categorical_values <- designs$categorical_values
  result$categorical_data <- designs$codes
  result$n_categories <- designs$n_categories
  result$measurement_model <- .latents_measurement_model(
    n_indicators > 0L, !is.null(designs$codes),
    colnames(designs$extra$ordinal) %||% character(),
    colnames(designs$extra$count) %||% character())
  result$extra_data <- designs$extra
  result$ordinal <- colnames(designs$extra$ordinal) %||% character()
  result$count <- colnames(designs$extra$count) %||% character()
  result <- .latents_label_extra(result, designs$extra,
                                 paste0("profile_", seq_len(n_profiles)))
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
  ## Whether the likelihood itself stopped moving, apart from the membership
  ## logits. A logit drifting to infinity on an empty or separated class leaves
  ## the likelihood at its supremum, which is all a likelihood-ratio statistic
  ## reads; bootstrap_lrt() accepts such a replicate, inference does not.
  result$likelihood_converged <- isTRUE(best$likelihood_settled)
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
