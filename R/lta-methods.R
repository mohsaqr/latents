# Tables and S3 methods for the general transition model (multilpa_lta).

utils::globalVariables(c("from", "to", "probability", "occasion"))

.lta_tables <- function(fit = NULL) {
  c("transitions", "transition_coefficients",
    if (!is.null(fit) && isTRUE(fit$order >= 2L))
      c("second_order_transitions", "second_order_coefficients"),
    "initial", "initial_coefficients", "profiles",
    if (!is.null(fit$extra_data$ordinal)) "ordinal",
    if (!is.null(fit$extra_data$count)) "count_means",
    "group_classes", "assignments", "fit", "starts")
}

#' Coefficient table (log odds) with Wald columns when inference applies
#' @noRd
.lta_coefficient_table <- function(fit, block, level, vcov_type) {
  inference <- tryCatch(.lta_inference(fit, vcov_type), error = function(condition) {
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
  parts <- do.call(rbind, strsplit(names_all[selected], ".", fixed = TRUE))
  # name layout: block.group_class.profile(s).term ; a term may contain dots.
  pieces <- strsplit(names_all[selected], ".", fixed = TRUE)
  table <- data.frame(
    group_class = vapply(pieces, `[`, character(1), 2L),
    outcome = vapply(pieces, `[`, character(1), 3L),
    term = vapply(pieces, function(p) paste(p[-(1:3)], collapse = "."), character(1)),
    estimate = theta[selected], stringsAsFactors = FALSE)
  rm(parts)
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
  z <- stats::qnorm(1 - (1 - level) / 2)
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
  n_states <- fit$n_profiles
  profiles <- paste0("profile_", seq_len(n_states))
  classes <- names(fit$group_probabilities)
  if (identical(what, "initial")) {
    return(do.call(rbind, lapply(seq_len(fit$n_group_classes), function(h) {
      p <- exp(.lta_log_initial(fit$designs$initial,
                                matrix(parameters$initial[, , h], ncol(fit$designs$initial))))
      data.frame(group_class = classes[h], profile = profiles,
                 probability = colMeans(p), stringsAsFactors = FALSE)
    })))
  }
  by_occasion <- !identical(fit$transitions, "homogeneous") ||
    length(fit$transition_covariates) > 0L
  rows <- lapply(seq_len(fit$n_group_classes), function(h) {
    per_occasion <- lapply(seq_along(fit$designs$transition), function(t) {
      active <- fit$layout$within[, t + 1L]
      log_p <- if (isTRUE(parameters$stayer[h])) {
        .lta_identity_transition(sum(active), n_states)
      } else .lta_log_transition(fit$designs$transition[[t]][active, , drop = FALSE],
                                 array(parameters$transition[, , , h],
                                       dim(parameters$transition)[1:3]))
      apply(exp(log_p), c(2L, 3L), mean)
    })
    weights <- vapply(seq_along(fit$designs$transition), function(t) {
      sum(fit$layout$within[, t + 1L])
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
  n_states <- fit$n_profiles
  profiles <- paste0("profile_", seq_len(n_states))
  classes <- names(fit$group_probabilities)
  do.call(rbind, lapply(seq_len(fit$n_group_classes), function(h) {
    occasions <- seq_along(fit$designs$transition)[-1L]
    weights <- vapply(occasions, function(t) sum(fit$layout$within[, t + 1L]), numeric(1))
    averaged <- Reduce(`+`, Map(function(t, w) {
      active <- fit$layout$within[, t + 1L]
      log_p <- .lta_log_transition2(fit$designs$transition[[t]][active, , drop = FALSE],
                                    array(parameters$transition2[, , , h],
                                          dim(parameters$transition2)[1:3]))
      apply(exp(log_p), c(2L, 3L), mean) * w
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
      data.frame(occasion = if (length(fit$measurement) == 1L) NA_integer_ else b,
                 profile = profiles[cells$profile],
                 indicator = fit$continuous[cells$indicator],
                 mean = block$means[cbind(cells$profile, cells$indicator)],
                 variance = block$variances[cbind(cells$profile, cells$indicator)],
                 stringsAsFactors = FALSE)
    })),
    ordinal = .lta_extra_table(fit, "ordinal"),
    count_means = .lta_extra_table(fit, "count_means"),
    group_classes = data.frame(group_class = classes,
                               weight = unname(fit$group_probabilities),
                               count = unname(colSums(fit$group_posteriors)),
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

#' Tables of a latent transition fit with occasion- or covariate-dependent
#' transitions or occasion-specific measurement
#'
#' The result of [lta()] with `transitions = "occasion"`,
#' `transition_covariates`, `initial_covariates` or `measurement = "occasion"`.
#'
#' @param x A `multilpa_lta` fit.
#' @param what Which table; `"all"` for a named list of every table.
#' @param level Confidence level of the coefficient intervals.
#' @param vcov_type `"observed"`, `"robust"` (CR0 over groups) or `"opg"`.
#' @param ... Unused.
#' @section Tables:
#' \describe{
#'   \item{`transitions`}{Model-implied transition probabilities, one row per
#'     group class, origin `from` and destination `to`; with occasion-varying
#'     or covariate-dependent transitions, one set per destination
#'     `occasion`, averaged over the groups observed there.}
#'   \item{`transition_coefficients`}{One row per group class, move `from` ->
#'     `to` and design `term`: the log odds of that move rather than staying
#'     in `from`, its standard error, Wald statistic, p-value, interval and
#'     `odds_ratio`. Terms are `(Intercept)` or `occasion_t` and each
#'     covariate column.}
#'   \item{`second_order_transitions`, `second_order_coefficients`}{With
#'     `order = 2`: the probability of moving `to` given the profiles at the
#'     two previous occasions (`previous`, `from`), averaged over the groups
#'     observed from the third occasion on, and its logit coefficients
#'     (reference: staying in `from`). The `transitions` tables then describe
#'     the first move only.}
#'   \item{`initial`, `initial_coefficients`}{The initial distribution
#'     (averaged over groups) and its logit coefficients, reference: the last
#'     profile.}
#'   \item{`profiles`}{Means and variances per profile and indicator, per
#'     `occasion` when measurement is occasion-specific.}
#'   \item{`group_classes`, `assignments`, `fit`, `starts`}{As for other fits.}
#' }
#' @return A base `data.frame`, or a named list of them for `what = "all"`.
#' @examples
#' # Does the grade in the previous course shift the chance of switching
#' # engagement profile? `previous_grade` changes from course to course.
#' fit <- lta(course_engagement, c("browse", "lectures"), "student",
#'            n_profiles = 2, time = "sequence",
#'            transition_covariates = "previous_grade", n_starts = 2, seed = 1)
#' get_results(fit, "transition_coefficients")
#' get_results(fit, "transitions")
#' @export
get_results.multilpa_lta <- function(x, what = "transitions", level = 0.95,
                                     vcov_type = c("observed", "robust", "opg"), ...) {
  what <- match.arg(what, c(.lta_tables(x), "all"))
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), match.arg(vcov_type),
                                      missing(vcov_type))
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.lta_tables(x), function(name) {
      .lta_table(x, name, level, vcov_type)
    }), .lta_tables(x)))
  }
  .lta_table(x, what, level, vcov_type)
}

#' @rdname get_results.multilpa_lta
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.multilpa_lta <- function(x, row.names = NULL, optional = FALSE,
                                       what = "transitions", ...) {
  get_results.multilpa_lta(x, what = what, ...)
}

#' @rdname get_results.multilpa_lta
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
  cat(sprintf(paste0("Latent transition model: %d profiles, %d group classes (%s)\n",
                     "%d groups, %d observations, %d occasions; log-likelihood %.3f, ",
                     "%d parameters, BIC %.3f\nconverged: %s\n"),
              x$n_profiles, x$n_group_classes, paste(pieces, collapse = "; "),
              x$n_groups, x$n_observations, x$n_occasions, x$log_likelihood,
              x$n_parameters, x$bic, if (x$converged) "yes" else "no"))
  .latents_print_extra(x)
  .latents_print_weights(x)
  cat(paste("Tables: get_results(x, what = ), e.g. \"transitions\",",
            "\"transition_coefficients\", \"profiles\"; summary(x).\n"))
  invisible(x)
}

#' @rdname get_results.multilpa_lta
#' @param object A `multilpa_lta` fit.
#' @export
summary.multilpa_lta <- function(object, ...) {
  get_results.multilpa_lta(object, "all")
}

#' @rdname get_results.multilpa_lta
#' @export
coef.multilpa_lta <- function(object, ...) {
  stats::setNames(.lta_pack(.lta_parameters(object), object$variance_model,
                            object$extra_data),
                  .lta_names(object))
}

#' @rdname get_results.multilpa_lta
#' @param type Covariance type, as `vcov_type`.
#' @export
vcov.multilpa_lta <- function(object, type = c("observed", "robust", "opg"), ...) {
  .lta_inference(object, .latents_weighted_vcov(.latents_is_weighted(object),
                                                match.arg(type), missing(type)))$vcov
}

#' @rdname get_results.multilpa_lta
#' @export
logLik.multilpa_lta <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters, nobs = object$n_groups,
            class = "logLik")
}

#' @rdname get_results.multilpa_lta
#' @export
nobs.multilpa_lta <- function(object, ...) object$n_groups

#' @rdname parameter_inference
#' @export
parameter_inference.multilpa_lta <- function(x, data = NULL, level = 0.95,
                                             step = 1e-4,
                                             vcov_type = c("observed", "robust", "opg"),
                                             ...) {
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), match.arg(vcov_type),
                                      missing(vcov_type))
  inference <- .lta_inference(x, vcov_type, step)
  table <- rbind(
    cbind(block = "transition", .lta_coefficient_table(x, "transition", level, vcov_type)[,
      c("group_class", "from", "to", "term", "estimate", "standard_error",
        "statistic", "p_value", "conf_low", "conf_high")]),
    cbind(block = "initial", transform(
      .lta_coefficient_table(x, "initial", level, vcov_type), from = NA_character_,
      to = profile)[, c("group_class", "from", "to", "term", "estimate",
                        "standard_error", "statistic", "p_value", "conf_low",
                        "conf_high")]))
  structure(table, covariance_unconstrained = inference$vcov,
            group_scores = inference$group_scores, vcov_type = vcov_type, level = level)
}

#' Plot transition probabilities of a general transition fit
#'
#' One panel per group class (and occasion, when transitions vary): a tile
#' per origin and destination, labelled with its probability.
#' @param x A `multilpa_lta` fit.
#' @param main,subtitle Optional title and subtitle.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @export
plot.multilpa_lta <- function(x, main = NULL, subtitle = NULL, ...) {
  .gg_require()
  frame <- .lta_table(x, "transitions")
  frame$panel <- if (all(is.na(frame$occasion))) frame$group_class else
    paste(frame$group_class, "occasion", frame$occasion)
  ggplot2::ggplot(frame, ggplot2::aes(x = to, y = from, fill = probability)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", probability)), size = 3.2) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#4A6FE3", limits = c(0, 1)) +
    ggplot2::facet_wrap(~panel) +
    .gg_theme() +
    ggplot2::labs(x = "To", y = "From") +
    .gg_titles(main, subtitle, "Transition probabilities",
               "Model-implied, averaged over the groups observed", NULL)
}
