# Tables and S3 methods for the general transition model (multilpa_lta).

utils::globalVariables(c("from", "to", "probability", "occasion"))

.lta_tables <- function(fit = NULL) {
  c("transitions", "transition_coefficients",
    if (!is.null(fit) && isTRUE(fit$order >= 2L))
      c("second_order_transitions", "second_order_coefficients"),
    "initial", "initial_coefficients", "profiles",
    if (length(fit$categorical) > 0L) "responses",
    if (!is.null(fit$extra_data$ordinal)) "ordinal",
    if (!is.null(fit$extra_data$count)) "count_means",
    "group_classes", "assignments", "fit", "starts",
    if (!is.null(fit$gaussian_diagnostics) && nrow(fit$gaussian_diagnostics) > 0L)
      "gaussian_diagnostics")
}

#' Coefficient table (log odds) with Wald columns when inference applies
#' @noRd
.lta_coefficient_table <- function(fit, block, level, vcov_type, inference = NULL) {
  # `parameter_inference()` passes the inference it computed with the caller's
  # `step`; the tables compute it at the default step.
  inference <- inference %||% tryCatch(.lta_inference(fit, vcov_type), error = function(condition) {
    refusals <- c("latents_unsupported_inference", "latents_no_converge",
                  "latents_boundary_fit", "latents_singular_information",
                  "latents_too_few_groups")
    if (!inherits(condition, refusals)) stop(condition)
    .multilpa_notice(sprintf("No standard errors in this table: %s",
                             conditionMessage(condition)),
                     class = "latents_no_standard_errors")
    NULL
  })
  names_all <- .lta_names(fit)
  theta <- .lta_theta(fit)
  selected <- startsWith(names_all, paste0(block, "."))
  # name layout: block.group_class.profile(s).term ; a term may contain dots,
  # so the pieces differ in length and are read by position, never bound.
  pieces <- strsplit(names_all[selected], ".", fixed = TRUE)
  table <- data.frame(
    group_class = vapply(pieces, `[`, character(1), 2L),
    outcome = vapply(pieces, `[`, character(1), 3L),
    term = vapply(pieces, function(p) paste(p[-(1:3)], collapse = "."), character(1)),
    estimate = theta[selected], stringsAsFactors = FALSE)
  if (identical(block, "transition")) {
    moves <- strsplit(table$outcome, "->", fixed = TRUE)
    table$from <- vapply(moves, `[`, character(1), 1L)
    table$to <- vapply(moves, `[`, character(1), 2L)
    table$outcome <- NULL
    table <- table[, c("group_class", "from", "to", "term", "estimate")]
  } else if (identical(block, "transition2")) {
    moves <- strsplit(table$outcome, "->", fixed = TRUE)
    table$previous <- vapply(moves, `[`, character(1), 1L)
    table$from <- vapply(moves, `[`, character(1), 2L)
    table$to <- vapply(moves, `[`, character(1), 3L)
    table$outcome <- NULL
    table <- table[, c("group_class", "previous", "from", "to", "term", "estimate")]
  } else {
    names(table)[names(table) == "outcome"] <- "profile"
  }
  se <- if (is.null(inference)) rep(NA_real_, nrow(table)) else
    sqrt(pmax(diag(inference$vcov)[selected], 0))
  z <- .inference_critical(level, "one_minus")
  table$standard_error <- unname(se)
  table$statistic <- table$estimate / table$standard_error
  table$p_value <- 2 * stats::pnorm(-abs(table$statistic))
  table$conf_low <- table$estimate - z * table$standard_error
  table$conf_high <- table$estimate + z * table$standard_error
  table$odds_ratio <- exp(table$estimate)
  # A logit the likelihood does not identify (a move into an occasion with
  # nothing observed, see .lta_identified_columns()) is reported, with no
  # standard error, but flagged rather than passed off as an estimate.
  table$estimated <- unname(.lta_identified_coordinates(fit)[selected])
  rownames(table) <- NULL
  table
}

#' Model-implied probabilities averaged over the groups observed at each
#' occasion (transitions) or over groups (initial distribution)
#'
#' A move no unit of positive weight makes (an occasion reached only by units
#' of weight zero) is absent from the likelihood, so it is left out of the
#' table, as it would be had those units not been supplied (TRANS-4).
#' `estimated` is `FALSE` for a row that rests on a logit the likelihood does
#' not identify (`.lta_identified_coordinates()`): the move into an empty
#' occasion, or out of an empty origin.
#' @noRd
.lta_probability_table <- function(fit, what) {
  parameters <- .lta_parameters(fit)
  # A weighted fit describes the population, so its groups are averaged with
  # their sampling weights; unweighted, every group counts once.
  unit_weights <- if (is.null(fit$sampling_weights)) NULL else unname(fit$sampling_weights)
  n_states <- fit$n_profiles
  profiles <- paste0("profile_", seq_len(n_states))
  classes <- names(fit$group_probabilities)
  if (identical(what, "initial")) {
    return(do.call(rbind, lapply(seq_len(fit$n_group_classes), function(h) {
      p <- exp(.lta_log_initial(fit$designs$initial,
                                matrix(parameters$initial[, , h], ncol(fit$designs$initial))))
      data.frame(group_class = classes[h], profile = profiles,
                 probability = .lta_group_mean(p, unit_weights), stringsAsFactors = FALSE)
    })))
  }
  by_occasion <- !identical(fit$transitions, "homogeneous") ||
    length(fit$transition_covariates) > 0L
  transition_occasions <- if ((fit$order %||% 1L) >= 2L) 1L else
    seq_along(fit$designs$transition)
  weights <- .lta_occasion_weights(fit, transition_occasions, unit_weights)
  transition_occasions <- transition_occasions[weights > 0]
  weights <- weights[weights > 0]
  identified <- stats::setNames(.lta_identified_coordinates(fit), .lta_names(fit))
  rows <- lapply(seq_len(fit$n_group_classes), function(h) {
    per_occasion <- lapply(transition_occasions, function(t) {
      active <- fit$layout$within[, t + 1L]
      log_p <- if (isTRUE(parameters$stayer[h])) {
        .lta_identity_transition(sum(active), n_states)
      } else .lta_log_transition(fit$designs$transition[[t]][active, , drop = FALSE],
                                 array(parameters$transition[, , , h],
                                       dim(parameters$transition)[1:3]))
      .lta_group_mean(exp(log_p), unit_weights[active])
    })
    terms <- lapply(transition_occasions, function(t) .lta_terms_used(fit, t, unit_weights))
    if (!by_occasion) {
      per_occasion <- list(Reduce(`+`, Map(`*`, per_occasion, weights)) / sum(weights))
      terms <- list(unique(unlist(terms)))
    }
    do.call(rbind, lapply(seq_along(per_occasion), function(t) {
      cells <- expand.grid(from = seq_len(n_states), to = seq_len(n_states))
      estimated <- vapply(seq_len(n_states), function(k) {
        .lta_logits_identified(identified, "transition", classes[h], profiles[k],
                               profiles[-k], terms[[t]])
      }, logical(1))
      data.frame(group_class = classes[h],
                 occasion = if (by_occasion) transition_occasions[t] + 1L else NA_integer_,
                 from = profiles[cells$from], to = profiles[cells$to],
                 probability = per_occasion[[t]][cbind(cells$from, cells$to)],
                 estimated = estimated[cells$from],
                 stringsAsFactors = FALSE)
    }))
  })
  table <- do.call(rbind, rows)
  table[order(table$group_class, table$occasion, table$from, table$to), , drop = FALSE] |>
    `rownames<-`(NULL)
}

#' Second-order transition probabilities, averaged over the groups observed
#' at each occasion from the third on
#' @noRd
.lta_second_order_table <- function(fit) {
  parameters <- .lta_parameters(fit)
  unit_weights <- if (is.null(fit$sampling_weights)) NULL else unname(fit$sampling_weights)
  n_states <- fit$n_profiles
  profiles <- paste0("profile_", seq_len(n_states))
  classes <- names(fit$group_probabilities)
  occasions <- seq_along(fit$designs$transition)[-1L]
  weights <- .lta_occasion_weights(fit, occasions, unit_weights)
  occasions <- occasions[weights > 0]
  weights <- weights[weights > 0]
  identified <- stats::setNames(.lta_identified_coordinates(fit), .lta_names(fit))
  terms <- unique(unlist(lapply(occasions, .lta_terms_used, fit = fit,
                                unit_weights = unit_weights)))
  do.call(rbind, lapply(seq_len(fit$n_group_classes), function(h) {
    averaged <- Reduce(`+`, Map(function(t, w) {
      active <- fit$layout$within[, t + 1L]
      log_p <- if (isTRUE(parameters$stayer[h])) {
        .lta_identity_transition2(sum(active), n_states)
      } else .lta_log_transition2(fit$designs$transition[[t]][active, , drop = FALSE],
                                    array(parameters$transition2[, , , h],
                                          dim(parameters$transition2)[1:3]))
      .lta_group_mean(exp(log_p), unit_weights[active]) * w
    }, occasions, weights)) / sum(weights)
    cells <- expand.grid(to = seq_len(n_states), pair = seq_len(n_states^2))
    previous <- (cells$pair - 1L) %/% n_states + 1L
    from <- (cells$pair - 1L) %% n_states + 1L
    estimated <- vapply(seq_len(n_states^2), function(pair) {
      i <- (pair - 1L) %/% n_states + 1L
      j <- (pair - 1L) %% n_states + 1L
      .lta_logits_identified(identified, "transition2", classes[h],
                             paste0(profiles[i], "->", profiles[j]), profiles[-j], terms)
    }, logical(1))
    data.frame(group_class = classes[h],
               previous = profiles[previous],
               from = profiles[from],
               to = profiles[cells$to],
               probability = averaged[cbind(cells$pair, cells$to)],
               estimated = estimated[cells$pair],
               stringsAsFactors = FALSE)
  }))
}

#' Weight of the moves into each occasion
#'
#' The groups whose span covers the move, summed with their sampling weights
#' (counted once each when unweighted). A move only units of weight zero make
#' weighs nothing and is absent from the likelihood.
#' @param fit A `multilpa_lta` fit.
#' @param occasions Indices into `fit$designs$transition` (the move into
#'   occasion `t + 1`).
#' @param unit_weights Sampling weight of each group, or `NULL`.
#' @return One number per occasion.
#' @noRd
.lta_occasion_weights <- function(fit, occasions, unit_weights) {
  vapply(occasions, function(t) {
    sum((unit_weights %||% 1) * fit$layout$within[, t + 1L])
  }, numeric(1))
}

#' Transition design terms that a move into an occasion reads
#' @param fit A `multilpa_lta` fit.
#' @param t Index into `fit$designs$transition`.
#' @param unit_weights Sampling weight of each group, or `NULL`.
#' @return The names of the design columns not zero on every group whose span
#'   covers the move (and that has positive weight).
#' @noRd
.lta_terms_used <- function(fit, t, unit_weights) {
  positive <- if (is.null(unit_weights)) TRUE else unit_weights > 0
  active <- fit$layout$within[, t + 1L] & positive
  design <- fit$designs$transition[[t]][active, , drop = FALSE]
  colnames(design)[colSums(abs(design)) > 0]
}

#' Are the logits behind one origin's row of probabilities identified?
#' @param identified Named logical, from `.lta_identified_coordinates()`.
#' @param kind `"transition"` or `"transition2"`.
#' @param class_label,origin_label The group class and origin as named.
#' @param destinations The destination profiles other than staying.
#' @param terms The design terms the row reads.
#' @return A single logical: `TRUE` when every such logit is identified (and
#'   for a stayer class, which has none).
#' @noRd
.lta_logits_identified <- function(identified, kind, class_label, origin_label,
                                   destinations, terms) {
  wanted <- as.vector(outer(sprintf("%s.%s.%s->%s.", kind, class_label, origin_label,
                                    destinations), terms, paste0))
  all(identified[intersect(wanted, names(identified))])
}

#' Average a groups-first array over its groups
#'
#' @param values A groups x ... array (or groups x K matrix).
#' @param weights Sampling weights of those groups, or `NULL` for a plain mean.
#' @return The remaining dimensions, averaged.
#' @noRd
.lta_group_mean <- function(values, weights = NULL) {
  dimensions <- seq_along(dim(values))[-1L]
  if (is.null(weights)) {
    return(if (length(dimensions) == 1L) colMeans(values) else
      apply(values, dimensions, mean))
  }
  apply(values, dimensions, function(v) sum(weights * v) / sum(weights))
}

#' Response probabilities of a general transition fit's categorical indicators
#'
#' The layout of the profile model's `"responses"` table, with an `occasion`
#' column (`NA` for invariant measurement) as in `"profiles"`.
#' @param fit A `multilpa_lta` fit.
#' @return One row per occasion block, indicator, category and profile.
#' @noRd
.lta_responses_table <- function(fit) {
  profiles <- paste0("profile_", seq_len(fit$n_profiles))
  do.call(rbind, lapply(seq_along(fit$measurement), function(b) {
    block <- fit$measurement[[b]]
    do.call(rbind, lapply(seq_along(fit$categorical), function(i) {
      indicator <- fit$categorical[[i]]
      probabilities <- block$response_probabilities[[i]]
      levels <- fit$categorical_levels[[indicator]] %||%
        as.character(seq_len(ncol(probabilities)))
      cells <- expand.grid(profile = seq_len(fit$n_profiles),
                           category = seq_len(ncol(probabilities)))
      data.frame(occasion = if (length(fit$measurement) == 1L) NA_integer_ else b,
                 profile = profiles[cells$profile], indicator = indicator,
                 category = levels[cells$category],
                 probability = probabilities[cbind(cells$profile, cells$category)],
                 stringsAsFactors = FALSE)
    }))
  }))
}

#' The fit as the shared assignment and alignment helpers read it
#'
#' `.multilpa_model_frame()`, `.multilpa_check_alignment()` and
#' `.multilpa_recovery_frame()` read a fit's identifiers, occasions, stored
#' indicators (with their centring) and group classes under the names the
#' profile and homogeneous transition fits use; this gives them those fields,
#' so the general fit's assignments follow the one documented layout.
#' @param fit A `multilpa_lta` fit.
#' @return A list of those fields.
#' @noRd
.lta_assignment_view <- function(fit) {
  continuous <- length(fit$continuous) > 0L
  list(id = fit$id, group_values = fit$group_values %||% fit$group_ids,
       group_index = fit$group_index, time = fit$time, time_values = fit$time_values,
       weights = fit$weights, sampling_weights = fit$sampling_weights,
       vars = fit$vars, n_observations = fit$n_observations,
       indicator_data = if (continuous) fit$x else NULL,
       centering_offsets = if (continuous)
         matrix(fit$centers, 1L, dimnames = list(NULL, colnames(fit$x))) else NULL,
       categorical_data = fit$codes, categorical_levels = fit$categorical_levels,
       extra_data = fit$extra_data, ordinal = fit$ordinal, count = fit$count,
       group_classes = max.col(fit$group_posteriors, ties.method = "first"))
}

#' One row per observation with its assigned profile and group class
#'
#' The layout of `get_results(fit, "assignments")` for every fit: the data
#' columns, `profile`, `group_class`, `uncertainty` and one
#' `posterior_profile_*` column per profile; with `truth`, the recovery table.
#' @param fit A `multilpa_lta` fit.
#' @param data The fitting data (with any further columns), or `NULL` for the
#'   columns the fit carries.
#' @param truth Column names of `data` holding known labels, or `NULL`.
#' @return A base `data.frame`.
#' @noRd
.lta_assignments <- function(fit, data = NULL, truth = NULL) {
  view <- .lta_assignment_view(fit)
  frame <- if (is.null(data)) .multilpa_model_frame(view) else {
    .latents_require("`data` must be a data frame" = is.data.frame(data))
    data
  }
  .multilpa_check_row_count(view, frame)
  .multilpa_check_alignment(view, frame)
  posteriors <- as.data.frame(unname(.multilpa_reported_posteriors(fit)))
  names(posteriors) <- sprintf("posterior_profile_%d", seq_len(fit$n_profiles))
  added <- data.frame(profile = max.col(fit$subject_posteriors, ties.method = "first"),
                      group_class = view$group_classes[fit$group_index],
                      uncertainty = pmax(1 - apply(fit$subject_posteriors, 1L, max), 0))
  # Keep the full input layout, but report classifications only for units
  # included by the same policy as the classification summaries.
  excluded <- !.multilpa_informative_rows(fit)
  added[excluded, c("profile", "uncertainty")] <- NA
  added$group_class[!.multilpa_informative_groups(fit)[fit$group_index]] <- NA_integer_
  clashes <- intersect(names(frame), c(names(added), names(posteriors)))
  if (length(clashes) > 0L) {
    stop(errorCondition(sprintf(
      "`data` already has a column named %s, which the assignments would overwrite.",
      paste(sprintf("`%s`", clashes), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  result <- cbind(frame, added, posteriors)
  row.names(result) <- NULL
  if (is.null(truth)) return(result)
  .multilpa_recovery_frame(result, truth, view)
}

#' Build one table of a general transition fit
#' @noRd
.lta_table <- function(fit, what, level = 0.95, vcov_type = "observed") {
  classes <- names(fit$group_probabilities)
  profiles <- paste0("profile_", seq_len(fit$n_profiles))
  switch(what,
    gaussian_diagnostics = fit$gaussian_diagnostics,
    transitions = .lta_probability_table(fit, "transitions"),
    initial = .lta_probability_table(fit, "initial"),
    transition_coefficients = .lta_coefficient_table(fit, "transition", level, vcov_type),
    second_order_transitions = .lta_second_order_table(fit),
    second_order_coefficients = .lta_coefficient_table(fit, "transition2", level, vcov_type),
    initial_coefficients = .lta_coefficient_table(fit, "initial", level, vcov_type),
    profiles = do.call(rbind, lapply(seq_along(fit$measurement), function(b) {
      block <- fit$measurement[[b]]
      cells <- expand.grid(profile = seq_len(fit$n_profiles),
                           indicator = seq_along(fit$continuous))
      data.frame(occasion = rep(if (length(fit$measurement) == 1L) NA_integer_ else b,
                                nrow(cells)),
                 profile = profiles[cells$profile],
                 indicator = fit$continuous[cells$indicator],
                 mean = block$means[cbind(cells$profile, cells$indicator)],
                 variance = block$variances[cbind(cells$profile, cells$indicator)],
                 stringsAsFactors = FALSE)
    })),
    responses = .lta_responses_table(fit),
    ordinal = .lta_extra_table(fit, "ordinal"),
    count_means = .lta_extra_table(fit, "count_means"),
    group_classes = data.frame(group_class = classes,
                               weight = unname(fit$group_probabilities),
                                count = unname(colSums(
                                  .multilpa_classified_group_posteriors(fit) *
                                    ((fit$sampling_weights %||% rep(1, nrow(fit$group_posteriors)))[
                                      .multilpa_informative_groups(fit)]))),
                               stringsAsFactors = FALSE),
    assignments = .lta_assignments(fit),
    fit = data.frame(n_profiles = fit$n_profiles, n_group_classes = fit$n_group_classes,
                     transitions = fit$transitions,
                     transition_covariates = paste(fit$transition_covariates, collapse = ", "),
                     initial_covariates = paste(fit$initial_covariates, collapse = ", "),
                     measurement = fit$measurement_model,
                     log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters,
                     aic = fit$aic, bic = fit$bic, bic_individual = fit$bic_individual,
                     n_groups = fit$n_groups, n_observations = fit$n_observations,
                     converged = fit$converged, boundary = fit$boundary,
                     n_best_replicated = fit$n_best_replicated, stringsAsFactors = FALSE),
    starts = fit$starts)
}

#' @noRd
#' @export
get_results.multilpa_lta <- function(x, what = "transitions", level = 0.95,
                                     vcov_type = c("observed", "robust", "opg"), ...) {
  .latents_require("`level` must be a single finite number in (0, 1)" =
    is.numeric(level) && length(level) == 1L && is.finite(level) && level > 0 && level < 1)
  # The one exact matcher of every family: an unknown name is a classed
  # refusal listing the tables, not match.arg()'s unclassed error (TRANS-5).
  what <- .multilpa_match_table(what, .lta_tables(x))
  # `assignments` alone takes further arguments; any other is refused, as the
  # profile families refuse it, rather than dropped (API-4).
  extra <- list(...)
  .multilpa_check_table_arguments(what, if (identical(what, "assignments"))
    c("data", "truth") else character(), extra)
  vcov_type <- .latents_weighted_vcov(
    .latents_is_weighted(x),
    .latents_match_arg(vcov_type, c("observed", "robust", "opg"), "vcov_type"),
    missing(vcov_type))
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.lta_tables(x), function(name) {
      .lta_table(x, name, level, vcov_type)
    }), .lta_tables(x)))
  }
  if (identical(what, "assignments")) return(do.call(.lta_assignments, c(list(x), extra)))
  .lta_table(x, what, level, vcov_type)
}

#' @noRd
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.multilpa_lta <- function(x, row.names = NULL, optional = FALSE,
                                       what = "transitions", ...) {
  get_results.multilpa_lta(x, what = what, ...)
}

#' @noRd
#' @export
print.multilpa_lta <- function(x, ...) {
  pieces <- c(if (isTRUE(x$mover_stayer)) "with a stayer class",
              if (isTRUE(x$order >= 2L)) "second-order transitions",
              if (!identical(x$transitions, "homogeneous")) "occasion-varying transitions",
              if (length(x$transition_covariates) > 0L)
                sprintf("transitions on %s", paste(x$transition_covariates, collapse = ", ")),
              if (length(x$initial_covariates) > 0L)
                sprintf("initial profiles on %s", paste(x$initial_covariates, collapse = ", ")),
              if (!identical(x$measurement_model, "invariant")) "occasion-specific measurement")
  if (length(pieces) == 0L) pieces <- "homogeneous transitions"
  cat(sprintf(paste0("Latent transition model: %d profiles, %d group %s (%s)\n",
                     "%d groups, %d observations, %d occasions; log-likelihood %.3f, ",
                     "%d parameters, BIC %.3f\nconverged: %s\n"),
              x$n_profiles, x$n_group_classes,
              if (x$n_group_classes == 1L) "class" else "classes",
              paste(pieces, collapse = "; "),
              x$n_groups, x$n_observations, x$n_occasions, x$log_likelihood,
              x$n_parameters, x$bic, if (x$converged) "yes" else "no"))
  .latents_print_extra(x)
  .latents_print_weights(x)
  .latents_print_gaussian_diagnostics(x$gaussian_diagnostics)
  cat(paste("Tables: get_results(x, what = ), e.g. \"transitions\",",
            "\"transition_coefficients\", \"profiles\"; summary(x).\n"))
  invisible(x)
}

#' Summary of a general transition fit
#'
#' Every table the fit offers, built once, under a class whose print shows
#' them section by section and whose `get_results()` serves them back, as the
#' other fits' summaries do (API-6).
#' @param object A `multilpa_lta` fit.
#' @param ... `level` and `vcov_type`, for the coefficient tables; anything
#'   else raises `latents_bad_argument`.
#' @return A `summary_multilpa_lta` object.
#' @noRd
#' @export
summary.multilpa_lta <- function(object, ...) {
  extra <- list(...)
  supplied <- names(extra) %||% rep("", length(extra))
  .multilpa_reject_extra_arguments(
    extra[is.na(match(supplied, c("level", "vcov_type")))], "summary()",
    "A transition fit's summary takes `level` and `vcov_type`.")
  structure(list(
    n_profiles = object$n_profiles, n_group_classes = object$n_group_classes,
    n_groups = object$n_groups, n_observations = object$n_observations,
    n_occasions = object$n_occasions, n_parameters = object$n_parameters,
    log_likelihood = object$log_likelihood, aic = object$aic, bic = object$bic,
    bic_individual = object$bic_individual, converged = object$converged,
    boundary = object$boundary,
    empty_transition_rows = any(object$empty_transition_rows) ||
      any(object$empty_second_order_rows),
    tables = do.call(get_results.multilpa_lta, c(list(object, "all"), extra))),
    class = "summary_multilpa_lta")
}

#' @noRd
#' @export
print.summary_multilpa_lta <- function(x, digits = 4L, rows = 10L, ...) {
  .latents_require("`x` must be a `summary_multilpa_lta` object" =
                     inherits(x, "summary_multilpa_lta"))
  .multilpa_check_print_arguments(digits, rows)
  cat(sprintf("Latent transition model: %d profiles and %d group %s\n",
              x$n_profiles, x$n_group_classes,
              if (x$n_group_classes == 1L) "class" else "classes"))
  cat(sprintf("Observations: %d; groups: %d; up to %d occasions\n",
              x$n_observations, x$n_groups, x$n_occasions))
  cat(sprintf("Parameters: %d; converged: %s\n", x$n_parameters, x$converged))
  cat(sprintf("Log likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              x$log_likelihood, x$aic, x$bic, x$bic_individual))
  if (!isTRUE(x$converged)) cat("WARNING: best start did not converge.\n")
  if (isTRUE(x$boundary)) cat("WARNING: the fit is at a boundary.\n")
  if (isTRUE(x$empty_transition_rows)) {
    cat("WARNING: a transition origin is never left; its logits are not estimated.\n")
  }
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' @noRd
#' @export
get_results.summary_multilpa_lta <- function(x, what = NULL, ...) {
  # The summary serves the tables it was built with, so get_results() on it
  # and on the fit are the same table.
  what <- .multilpa_match_table(what, names(x$tables))
  .multilpa_check_table_arguments(what, character(), list(...))
  if (identical(what, "all")) return(x$tables)
  x$tables[[what]]
}

#' @noRd
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.summary_multilpa_lta <- function(x, row.names = NULL, optional = FALSE, ...) {
  .multilpa_coerce(x, row.names, list(...))
}

#' @noRd
#' @export
coef.multilpa_lta <- function(object, ...) {
  # .lta_parameters() gives means in the engine's centred coordinates, which
  # inference needs; the user's are the data's own units.
  parameters <- .lta_parameters(object)
  if (length(object$continuous) > 0L) {
    parameters$measurement <- lapply(parameters$measurement, function(block) {
      block$means <- sweep(block$means, 2L, object$centers, "+")
      block
    })
  }
  stats::setNames(.lta_pack(parameters, object$variance_model, object$extra_data),
                  .lta_names(object))
}

#' @noRd
#' @param type Covariance type, as `vcov_type`.
#' @export
vcov.multilpa_lta <- function(object, type = c("observed", "robust", "opg"), ...) {
  # `vcov_type = "robust"`, the name parameter_inference() uses, used to fall
  # into `...` and return the observed-information covariance without a word.
  .latents_warn_unused_arguments(list(...), "vcov()",
                                 "Choose the covariance with `type =`.")
  .lta_inference(object, .latents_weighted_vcov(.latents_is_weighted(object),
                                                .latents_match_formal(type), missing(type)))$vcov
}

#' @noRd
#' @export
logLik.multilpa_lta <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = .multilpa_group_sample_size(object),
            class = "logLik")
}

#' @noRd
#' @export
nobs.multilpa_lta <- function(object, ...) .multilpa_group_sample_size(object)

#' @noRd
#' @export
parameter_inference.multilpa_lta <-function(x, data = NULL, level = 0.95,
                                             step = 1e-4,
                                             vcov_type = c("observed", "robust", "opg"),
                                             adjust = .multilpa_p_adjust_methods,
                                             method = c("wald", "bootstrap"),
                                             iter = 199L, n_starts = 10L,
                                             max_iter = 5000L, tol = 1e-10,
                                             seed = NULL,
                                             boundary = c("error", "fix")) {
  .latents_require("`level` must be a single finite number in (0, 1)" =
    is.numeric(level) && length(level) == 1L && is.finite(level) && level > 0 && level < 1)
  # Transition inference holds nothing at a bound, and a fit with an active
  # bound is refused (`latents_boundary_fit`) either way; so "fix" and "error"
  # agree wherever a table is returned.
  .latents_match_formal(boundary)
  adjust <- .latents_match_formal(adjust)
  if (identical(.latents_match_formal(method), "bootstrap")) {
    return(.lta_bootstrap_dispatch(x, data, level, iter, n_starts, max_iter, tol,
                                   seed, adjust))
  }
  if (!is.null(data)) .lta_check_data(x, data)
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), .latents_match_formal(vcov_type),
                                      missing(vcov_type))
  inference <- .lta_inference(x, vcov_type, step)
  table <- rbind(
    cbind(block = "transition", .lta_coefficient_table(x, "transition", level, vcov_type, inference)[,
      c("group_class", "from", "to", "term", "estimate", "standard_error",
        "statistic", "p_value", "conf_low", "conf_high", "estimated")]),
    cbind(block = "initial", transform(
      .lta_coefficient_table(x, "initial", level, vcov_type, inference),
      from = NA_character_,
      to = profile)[, c("group_class", "from", "to", "term", "estimate",
                        "standard_error", "statistic", "p_value", "conf_low",
                        "conf_high", "estimated")]))
  if ((x$order %||% 1L) >= 2L) {
    table$previous <- NA_character_
    second <- .lta_coefficient_table(x, "transition2", level, vcov_type, inference)
    second$block <- "transition2"
    table <- rbind(table, second[, names(table)])
  }
  # The family is the tests the table carries; p.adjust() leaves NA p-values
  # out of it.
  table$p_adjusted <- stats::p.adjust(table$p_value, method = adjust)
  ordered <- append(setdiff(names(table), "p_adjusted"), "p_adjusted",
                    after = match("p_value", names(table)))
  structure(table[, ordered], covariance_unconstrained = inference$vcov,
            group_scores = inference$group_scores, vcov_type = vcov_type, level = level,
            adjust = adjust)
}


#' @noRd
#' @export
plot.multilpa_lta <-function(x, main = NULL, subtitle = NULL, ...) {
  .latents_warn_unused_arguments(
    list(...), "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  .gg_require()
  frame <- .lta_table(x, "transitions")
  # One panel per group class (and occasion), in numeric occasion order --
  # as plain strings "occasion 10" would sort before "occasion 2" -- with the
  # occasion on its own line so narrow panels do not cut the label.
  panels <- if (all(is.na(frame$occasion))) frame$group_class else
    paste0(frame$group_class, "\noccasion ", frame$occasion)
  ordering <- order(frame$group_class, as.numeric(frame$occasion))
  frame$panel <- factor(panels, levels = unique(panels[ordering]))
  # The same white-to-blue probability ramp and legible ink as the profile
  # models' transition view (.gg_view_transitions()).
  frame$text_colour <- .gg_ink(.gg_probability_fill(frame$probability))
  ggplot2::ggplot(frame, ggplot2::aes(x = to, y = from, fill = probability)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", probability),
                                    colour = text_colour), size = 3.2) +
    ggplot2::scale_colour_identity() +
    .gg_probability_scale("Probability") +
    ggplot2::facet_wrap(~panel) +
    .gg_theme() +
    ggplot2::labs(x = "To", y = "From") +
    .gg_titles(main, subtitle, "Transition probabilities",
               "Model-implied, averaged over the groups observed", NULL)
}
