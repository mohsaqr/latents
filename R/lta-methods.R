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
    "group_classes", "assignments", "fit", "starts")
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
  theta <- .lta_pack(.lta_parameters(fit), fit$variance_model, fit$extra_data)
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
  rownames(table) <- NULL
  table
}

#' Model-implied probabilities averaged over the groups observed at each
#' occasion (transitions) or over groups (initial distribution)
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
    weights <- vapply(transition_occasions, function(t) {
      sum((unit_weights %||% 1) * fit$layout$within[, t + 1L])
    }, numeric(1))
    if (!by_occasion) per_occasion <- list(Reduce(`+`, Map(`*`, per_occasion, weights)) /
                                             sum(weights))
    do.call(rbind, lapply(seq_along(per_occasion), function(t) {
      cells <- expand.grid(from = seq_len(n_states), to = seq_len(n_states))
      data.frame(group_class = classes[h],
                 occasion = if (by_occasion) t + 1L else NA_integer_,
                 from = profiles[cells$from], to = profiles[cells$to],
                 probability = per_occasion[[t]][cbind(cells$from, cells$to)],
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
  do.call(rbind, lapply(seq_len(fit$n_group_classes), function(h) {
    occasions <- seq_along(fit$designs$transition)[-1L]
    weights <- vapply(occasions, function(t) {
      sum((unit_weights %||% 1) * fit$layout$within[, t + 1L])
    }, numeric(1))
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
    data.frame(group_class = classes[h],
               previous = profiles[(cells$pair - 1L) %/% n_states + 1L],
               from = profiles[(cells$pair - 1L) %% n_states + 1L],
               to = profiles[cells$to],
               probability = averaged[cbind(cells$pair, cells$to)],
               stringsAsFactors = FALSE)
  }))
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

#' Build one table of a general transition fit
#' @noRd
.lta_table <- function(fit, what, level = 0.95, vcov_type = "observed") {
  classes <- names(fit$group_probabilities)
  profiles <- paste0("profile_", seq_len(fit$n_profiles))
  switch(what,
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
                                count = unname(colSums(fit$group_posteriors *
                                                         (fit$sampling_weights %||% 1))),
                               stringsAsFactors = FALSE),
    assignments = {
      profile <- max.col(fit$subject_posteriors, ties.method = "first")
      group_class <- max.col(fit$group_posteriors, ties.method = "first")
      cbind(data.frame(row = seq_len(fit$n_observations)),
            stats::setNames(data.frame(fit$group_ids[fit$group_index],
                                       stringsAsFactors = FALSE), fit$id),
            data.frame(occasion = fit$occasion_of_row, profile = profiles[profile],
                       posterior = fit$subject_posteriors[cbind(seq_along(profile), profile)],
                       group_class = classes[group_class[fit$group_index]],
                       stringsAsFactors = FALSE))
    },
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
  stopifnot("`level` must be a single finite number in (0, 1)" =
    is.numeric(level) && length(level) == 1L && is.finite(level) && level > 0 && level < 1)
  what <- match.arg(what, c(.lta_tables(x), "all"))
  # No table of this fit takes a further argument. A released call that passed
  # one (`estimated = FALSE`, say) still gets its table, but is told the
  # argument changed nothing rather than left to read the table as filtered.
  .latents_warn_unused_arguments(
    list(...), "get_results()",
    "The tables of a general transition fit take only `what`, `level` and `vcov_type`.")
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), match.arg(vcov_type),
                                      missing(vcov_type))
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.lta_tables(x), function(name) {
      .lta_table(x, name, level, vcov_type)
    }), .lta_tables(x)))
  }
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
  cat(paste("Tables: get_results(x, what = ), e.g. \"transitions\",",
            "\"transition_coefficients\", \"profiles\"; summary(x).\n"))
  invisible(x)
}

#' @noRd
#' @param object A `multilpa_lta` fit.
#' @export
summary.multilpa_lta <- function(object, ...) {
  get_results.multilpa_lta(object, "all")
}

#' @noRd
#' @export
coef.multilpa_lta <- function(object, ...) {
  stats::setNames(.lta_pack(.lta_parameters(object), object$variance_model,
                            object$extra_data),
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
                                                match.arg(type), missing(type)))$vcov
}

#' @noRd
#' @export
logLik.multilpa_lta <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters, nobs = object$n_groups,
            class = "logLik")
}

#' @noRd
#' @export
nobs.multilpa_lta <- function(object, ...) object$n_groups

#' @details For a general transition fit (an [lta()] fit with covariate,
#'   occasion-varying or second-order transitions), the table has one row per
#'   transition, initial and second-order logit: `block`, `group_class`,
#'   `from`, `to`, `term`, `estimate`, `standard_error`, `statistic`,
#'   `p_value`, `p_adjusted`, `conf_low`, `conf_high` (and `previous` for a
#'   second-order fit). `method = "bootstrap"` resamples persons, refits, and
#'   matches every replicate's profiles (on the measurement) and group classes
#'   to the original's, rebasing the initial logits when the reference profile
#'   moves; it needs `data`, and works where Wald inference is refused (for
#'   example full covariance structures). Transition inference holds nothing
#'   at a bound: with no bound active `boundary = "fix"` gives the same table
#'   as `"error"`, and a fit with an active bound raises `latents_boundary_fit`
#'   under either.
#' @rdname parameter_inference
#' @export
parameter_inference.multilpa_lta <- function(x, data = NULL, level = 0.95,
                                             step = 1e-4,
                                             vcov_type = c("observed", "robust", "opg"),
                                             adjust = .multilpa_p_adjust_methods,
                                             method = c("wald", "bootstrap"),
                                             iter = 199L, n_starts = 10L,
                                             max_iter = 5000L, tol = 1e-10,
                                             seed = NULL,
                                             boundary = c("error", "fix")) {
  stopifnot("`level` must be a single finite number in (0, 1)" =
    is.numeric(level) && length(level) == 1L && is.finite(level) && level > 0 && level < 1)
  # Transition inference holds nothing at a bound, and a fit with an active
  # bound is refused (`latents_boundary_fit`) either way; so "fix" and "error"
  # agree wherever a table is returned.
  match.arg(boundary)
  adjust <- match.arg(adjust)
  if (identical(match.arg(method), "bootstrap")) {
    return(.lta_bootstrap_dispatch(x, data, level, iter, n_starts, max_iter, tol,
                                   seed, adjust))
  }
  if (!is.null(data)) .lta_check_data(x, data)
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), match.arg(vcov_type),
                                      missing(vcov_type))
  inference <- .lta_inference(x, vcov_type, step)
  table <- rbind(
    cbind(block = "transition", .lta_coefficient_table(x, "transition", level, vcov_type, inference)[,
      c("group_class", "from", "to", "term", "estimate", "standard_error",
        "statistic", "p_value", "conf_low", "conf_high")]),
    cbind(block = "initial", transform(
      .lta_coefficient_table(x, "initial", level, vcov_type, inference),
      from = NA_character_,
      to = profile)[, c("group_class", "from", "to", "term", "estimate",
                        "standard_error", "statistic", "p_value", "conf_low",
                        "conf_high")]))
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


#' @rdname latents-plots
#' @export
plot.multilpa_lta <- function(x, main = NULL, subtitle = NULL, ...) {
  .latents_warn_unused_arguments(
    list(...), "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  .gg_require()
  frame <- .lta_table(x, "transitions")
  frame$panel <- if (all(is.na(frame$occasion))) frame$group_class else
    paste(frame$group_class, "occasion", frame$occasion)
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
