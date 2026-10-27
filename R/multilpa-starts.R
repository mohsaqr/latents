#' Run one prepared profile-model start
#' @param start_index Original start index.
#' @param initial Starting values, or NULL when resuming.
#' @param pause,resume Screening budget and checkpoint.
#' @param context Forced EM inputs from setup.
#' @return One EM attempt.
#' @noRd
.multilpa_run_start <- function(start_index, initial, pause = Inf, resume = NULL,
                                context) {
  with(context, {
    # A resumed start continues from its checkpoint, not from new values.
    if (!is.null(resume)) initial <- resume$state else
      initial <- .latents_mask_covariances(initial, covariance_pattern)
    .multilpa_em(x, group_index, initial,
               variance_model, min_variance, max_iter,
               tol, covariance_model, codes, n_categories, min_probability,
               held, structure, prior_parameters,
               accelerate_em = identical(acceleration, "squarem"),
               sampling_weights = sampling_weights, extra = fitting_extra,
               covariance_pattern = covariance_pattern, pause = pause,
               resume = resume)
  })
}

#' Identify a start clear of the covariance lower bound
#' @param attempt One completed EM attempt.
#' @param context Forced variance-bound setting.
#' @return TRUE when the start is interior.
#' @noRd
.multilpa_start_interior <- function(attempt, context) {
  !.multilpa_covariance_boundary(attempt$parameters, context$min_variance)
}

#' Run the starts of a prepared profile fit with the shared search
#' @param plan A prepared fitting plan.
#' @param cores Core budget for a direct fit.
#' @return A fitted profile or latent class model.
#' @noRd
.multilpa_run_plan <- function(plan, cores) {
  stopifnot("plan is a prepared profile fit" =
              inherits(plan, "latents_multilpa_plan"))
  starts <- tryCatch(.latents_run_starts(
    plan$n_starts, plan$run, plan$select_start, initial = plan$initial,
    cores = cores, screen = plan$screen, interior = plan$interior),
    error = function(error) {
      .multilpa_report_failed_plan(plan, error)
      stop(error)
    })
  .multilpa_finish_plan(plan, starts)
}

#' Report covariance qualifications of an unsuccessful profile search
#' @param plan A prepared fitting plan.
#' @param error The selection error carrying its attempts' count records.
#' @return NULL, invisibly.
#' @noRd
.multilpa_report_failed_plan <- function(plan, error) {
  .lta_report_failed_plan(plan, error)
}

#' Finalise a profile fit in the parent and report covariance qualifications
#' @param plan A prepared fitting plan.
#' @param starts Shared start-selection result.
#' @return A fitted profile or latent class model.
#' @noRd
.multilpa_finish_plan <- function(plan, starts) {
  .lta_finish_plan(plan, starts)
}

#' Finalise a profile fit from its selected starts
#' @param starts Shared start-selection result.
#' @param context Forced setup inputs.
#' @return A fitted multilpa object.
#' @noRd
.multilpa_finalise <- function(starts, context) {
  with(context, {
    found <- starts
    attempts <- starts$attempts
    valid <- starts$valid
    scores <- starts$scores
    best_start <- starts$best_start
    final_stage <- .latents_final_stage(found)
    best <- attempts[[best_start]]
    parameters <- best$parameters
    # The weighted expectation carries weighted posteriors, which are counts for
    # the M-step; each unit's own classification is the unweighted posterior.
    weighted_expectation <- best$expectation
    if (!is.null(sampling_weights)) {
      best$expectation <- .multilpa_expectation(x, group_index, parameters, codes,
                                                extra = fitting_extra)
    }
    subject_posteriors <- best$expectation$subject_posteriors
    subject_profiles <- max.col(subject_posteriors, ties.method = "first")
    noise_fields <- list(noise = isTRUE(noise))
    if (isTRUE(noise)) {
      split <- .multilpa_split_noise(parameters, subject_posteriors)
      parameters <- split$parameters
      subject_posteriors <- split$posteriors
      subject_profiles <- split$classes
      noise_fields <- list(
        noise = TRUE, noise_probability = split$noise_probability,
        noise_posteriors = stats::setNames(split$noise_posteriors, rownames(data)),
        hypervolume = exp(-noise_log_density),
        n_noise = sum(subject_profiles == 0L))
    }
    parameters$means <- sweep(parameters$means, 2L, centers, "+")
    parameters <- .latents_uncenter_direct(parameters, fitting_extra)
    profile_names <- paste0("profile_", seq_len(n_profiles))
    type_names <- paste0("group_class_", seq_len(n_group_classes))
    parameters <- .multilpa_label_parameters(parameters, profile_names, continuous,
                                           categorical, encoded$levels)
    parameters <- .latents_label_extra(parameters, extra, profile_names)
    dimnames(parameters$profile_probabilities) <- list(type_names, profile_names)
    names(parameters$group_probabilities) <- type_names
    group_posteriors <- best$expectation$group_posteriors
    dimnames(subject_posteriors) <- list(rownames(data), profile_names)
    dimnames(group_posteriors) <- list(group_ids, type_names)
    # A noise component adds its mixing proportion and, as mclust counts it
    # (`nMclustParams(noise = TRUE)`), the hypervolume estimated from the data:
    # two, so that every criterion equals mclust's.
    n_parameters_unconstrained <- .multilpa_count_parameters(structure = structure,
      n_profiles, n_group_classes, ncol(x), n_categories, variance_model,
      covariance_model)
    if (isTRUE(noise)) n_parameters_unconstrained <- n_parameters_unconstrained + 2L
    # The covariances a pattern holds at zero are not parameters.
    if (!is.null(covariance_pattern)) {
      n_parameters_unconstrained <- n_parameters_unconstrained -
        .latents_pattern_zeros(covariance_pattern) *
        if (identical(variance_model, "equal")) 1L else n_profiles
    }
    n_parameters_unconstrained <- n_parameters_unconstrained +
      .latents_extra_n_parameters(extra, n_profiles)
    n_parameters <- n_parameters_unconstrained -
      .multilpa_fixed_parameters(fixed, n_profiles, ncol(x), n_categories,
                                 variance_model, covariance_model)
    ## A held fit is comparable with other held fits on `n_parameters` and with a
    ## joint fit on the count that includes the measurement it was handed. Only
    ## `fit_staged()` used to record the second one, so `print()` fell through
    ## `%||%` on a `multilpa(fixed = )` fit and reported the same number twice.
    n_parameters_with_measurement <- if (length(fixed) == 0L) NULL else
      n_parameters_unconstrained
    log_likelihood <- weighted_expectation$log_likelihood
    starts <- do.call(rbind, lapply(seq_along(attempts), function(start_index) {
      attempt <- attempts[[start_index]]
      data.frame(start = start_index, log_likelihood = scores[start_index],
                 converged = if (valid[start_index]) attempt$converged else FALSE,
                 iterations = if (valid[start_index]) attempt$iterations else NA_integer_,
                 error = if (valid[start_index]) NA_character_ else attempt$error,
                 boundary = if (valid[start_index])
                   .multilpa_covariance_boundary(attempt$parameters %||% attempt$state,
                                                 min_variance) else NA)
    }))
    starts <- .latents_search_columns(starts, found, n_starts, supplied = !is.null(start))
    # A row with no observed indicator contributes nothing to the likelihood and
    # nothing to any estimate, so it must not enlarge the sample size the
    # individual-level BIC penalizes against.
    observed_per_row <- rowSums(!is.na(x)) +
      (if (is.null(codes)) 0L else rowSums(!is.na(codes))) +
      if (is.null(.latents_extra_matrix(extra))) 0L else
        rowSums(!is.na(.latents_extra_matrix(extra)))
    # A unit of sampling weight zero contributes nothing either, so its rows
    # are not classified and do not count in any sample size (ML-1).
    if (!is.null(sampling_weights)) {
      observed_per_row[sampling_weights[group_index] <= 0] <- 0L
    }
    # Reported by print() and summary() rather than warned about: resampling verbs
    # refit the same data many times, and a warning would fire once per replicate.
    n_informative <- sum(observed_per_row > 0L)
    # A group whose every row observes nothing is left out of the group sample
    # size the same way, as Mplus leaves it out of its clusters.
    n_informative_groups <- length(unique(group_index[observed_per_row > 0L]))
    boundary <- .multilpa_covariance_boundary(parameters, min_variance)
    # Effective memberships are the counts the fit reports: weighted for a
    # weighted fit, so the flag, its warning and the printed counts agree.
    # Under `drop_empty`, the counts are of the rows and groups that were
    # classified, as Mplus counts them; otherwise of every row and group.
    counted_rows <- if (isTRUE(drop_empty)) observed_per_row > 0L else
      rep(TRUE, nrow(x))
    counted_groups <- if (isTRUE(drop_empty))
      seq_len(n_groups) %in% group_index[observed_per_row > 0L] else rep(TRUE, n_groups)
    effective_profile_counts <- if (is.null(sampling_weights))
      colSums(subject_posteriors[counted_rows, , drop = FALSE]) else
        setNames(colSums(weighted_expectation$subject_posteriors[counted_rows, , drop = FALSE]),
                 profile_names)
    effective_group_counts <- if (is.null(sampling_weights))
      colSums(group_posteriors[counted_groups, , drop = FALSE]) else
        setNames(colSums(weighted_expectation$group_posteriors[counted_groups, , drop = FALSE]),
                 type_names)
    # `drop_empty = FALSE` counts every group but those of weight zero, which
    # are in no sample size under either convention (ZW-5).
    group_sample_size <- if (isTRUE(drop_empty)) n_informative_groups else
      if (is.null(sampling_weights)) n_groups else sum(sampling_weights > 0)
    small_classes <- any(effective_profile_counts < 1) ||
      any(effective_group_counts < 1)
    result <- c(parameters, list(
      call = call, vars = vars, continuous = continuous,
      continuous_types = vapply(data[continuous], typeof, character(1)),
      ordinal = ordinal, count = count, extra_data = extra,
      categorical = categorical,
      categorical_levels = encoded$levels, min_probability = min_probability,
      categorical_values = .multilpa_categorical_values(data, categorical, encoded),
      indicator_data = as.matrix(indicator_frame),
      centering = centering, centering_offsets = centred$offsets,
      scaling = scaling, scaling_factors = centred$scales,
      categorical_data = codes,
      id = id, group_ids = group_ids, single_level = single_level,
      weights = weights,
      sampling_weights = if (is.null(sampling_weights)) NULL else
        setNames(sampling_weights, group_ids),
      time = time, time_values = time_values,
      group_values = group_values, group_index = group_index,
      group_sizes = setNames(group_sizes, group_ids),
      n_observations = nrow(x), n_informative = n_informative,
      informative = observed_per_row > 0L,
      n_informative_groups = n_informative_groups, drop_empty = drop_empty,
      n_groups = n_groups, n_profiles = as.integer(n_profiles),
      n_group_classes = as.integer(n_group_classes), variance_model = variance_model,
      covariance_model = covariance_model, covariance_structure = structure,
      prior = prior, acceleration = acceleration,
      prior_parameters = if (is.null(prior_parameters)) NULL else
        utils::modifyList(prior_parameters,
                          list(mean = unname(prior_parameters$mean + centers))),
      missing = missing, max_iter = max_iter, tol = tol,
      n_observed_by_indicator = setNames(
        c(colSums(!is.na(x) & carried_rows), if (is.null(codes)) NULL else
          colSums(!is.na(codes) & carried_rows)),
        c(continuous, categorical)),
      min_variance = min_variance, standard_deviations = sqrt(parameters$variances),
      fixed = fixed,
      measurement_model = .latents_measurement_model(ncol(x) > 0L, !is.null(codes),
                                                     ordinal, count),
      subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
      subject_profiles = subject_profiles,
      group_classes = max.col(group_posteriors, ties.method = "first"),
      log_likelihood = log_likelihood,
      group_log_likelihood = setNames(best$expectation$group_log_likelihood, group_ids),
      n_parameters = n_parameters,
      n_parameters_with_measurement = n_parameters_with_measurement,
      aic = -2 * log_likelihood + 2 * n_parameters,
      bic = -2 * log_likelihood + log(group_sample_size) * n_parameters,
      bic_groups = -2 * log_likelihood + log(group_sample_size) * n_parameters,
      bic_individual = -2 * log_likelihood + log(n_informative) * n_parameters,
      converged = best$converged, iterations = best$iterations,
      log_likelihood_history = best$history, starts = starts, best_start = best_start,
      n_failed_starts = sum(!valid), boundary = boundary, small_classes = small_classes,
      effective_profile_counts = effective_profile_counts,
      effective_group_counts = effective_group_counts,
      n_best_replicated = sum(valid & final_stage &
                                abs(scores - log_likelihood) <=
                               1e-6 * (1 + abs(log_likelihood))),
      replication_tolerance = 1e-6 * (1 + abs(log_likelihood))), noise_fields,
      .latents_search_record(search, n_starts, screen, select_start,
                             if (exclude_boundary) "exclude" else "allow"))
    class(result) <- "multilpa"
    # Each row's classification entropy within its group class, from which the
    # entropy ICL, AWE and CLC charge is summed (review 2026-10-26 round 2,
    # S-1). Taken from the unweighted expectation, as the posteriors are.
    result$classification_entropy_rows <- if (is.null(best$expectation$joint)) NULL else
      .latents_row_conditional_entropy(best$expectation$joint,
                                       best$expectation$group_posteriors, group_index)
    # Recorded only when used, so a fit without local dependence has exactly the
    # fields it always had.
    if (!is.null(covariance_pattern)) {
      result$covariance_pattern <- covariance_pattern
      result$local_dependence <- .latents_classify_pairs(
        .latents_parse_pairs(local_dependence), categorical)
    }
    if (any(!valid)) {
      warning(warningCondition(sprintf(
        "%d of %d starts failed; see get_results(fit, \"starts\").",
        sum(!valid), length(attempts)), class = "latents_failed_starts", call = NULL))
    }
    # max_iter = 0 is a deliberate evaluate-only call, so non-convergence is
    # expected rather than an anomaly worth reporting.
    if (!best$converged && max_iter > 0L) {
      warning(warningCondition(
        "The best start did not converge; increase max_iter and see get_results(fit, \"starts\").",
        class = "latents_unconverged", call = NULL))
    }
    if (boundary) {
      warning(warningCondition(if (covariance_model == "full")
        "A covariance eigenvalue reached min_variance; this is a bound-active constrained fit." else
        "A variance reached min_variance; this is a bound-active constrained fit.",
        class = "latents_boundary", call = NULL))
    }
    if (small_classes) {
      warning(warningCondition(
        "A profile or group class has effective membership below one.",
        class = "latents_small_classes", call = NULL))
    }
    # Classification can include categorical-only or retained empty rows;
    # those rows do not support a Gaussian measurement density.
    gaussian_rows <- rowSums(!is.na(x)) > 0L
    gaussian_counts <- colSums(weighted_expectation$subject_posteriors[
      gaussian_rows, seq_len(n_profiles), drop = FALSE])
    .latents_warn_thin_profiles(
      gaussian_counts,
      .latents_profile_parameters(ncol(x), variance_model, covariance_model, structure,
                                  covariance_pattern = covariance_pattern, fixed = fixed))
    result$gaussian_diagnostics <- .latents_gaussian_diagnostics(
      parameters, x, subject_posteriors,
      weights = if (is.null(sampling_weights)) NULL else sampling_weights[group_index],
      min_variance = min_variance)
    .latents_warn_poisson_limit(result$count_dispersion)
    result
  })
}
