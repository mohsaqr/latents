# Comparing latent transition models: simulation from a fitted model, the
# nesting rule, the parametric bootstrap likelihood-ratio test, and
# enumeration over numbers of profiles and group classes.

#' The lta() arguments that define a transition fit, so it can be refitted
#' @noRd
.lta_arguments <- function(fit) {
  if (inherits(fit, "multilpa_lta")) return(fit$arguments)
  stopifnot(inherits(fit, "multilpa_transitions"))
  list(vars = fit$vars, id = fit$id, time = fit$time,
       n_profiles = fit$n_profiles, n_group_classes = fit$n_group_classes,
       variance_model = fit$variance_model, categorical = fit$categorical,
       occasions = fit$occasions, transitions = "homogeneous",
       transition_covariates = character(), initial_covariates = character(),
       measurement = "invariant", order = 1L, mover_stayer = FALSE,
       missing = fit$missing, covariance_model = fit$covariance_model,
       model = fit$covariance_structure,
       min_variance = fit$min_variance, min_probability = fit$min_probability)
}

#' A transition fit in the general engine's parameterization
#'
#' Homogeneous fits carry probabilities; they are turned into the same
#' logits the general engine uses (initial: reference last profile;
#' transitions: reference staying), on an intercept-only design.
#' @return A list with `parameters` (natural means), `designs`, `layout`,
#'   `occasion_of_row`, `n_profiles`, `vars`, `continuous`, `categorical`.
#' @noRd
.lta_general_view <- function(fit, data) {
  if (inherits(fit, "multilpa_lta")) {
    if (!identical(fit$covariance_model %||% "diagonal", "diagonal")) {
      stop(errorCondition("Simulation from a transition fit needs diagonal covariances.",
                          class = "latents_unsupported_inference", call = NULL))
    }
    return(list(parameters = list(measurement = fit$measurement,
                                  initial = fit$initial_coefficients,
                                  transition = fit$transition_coefficients,
                                  transition2 = fit$second_order_coefficients,
                                  group_probabilities = fit$group_probabilities,
                                  stayer = fit$stayer),
                designs = fit$designs, layout = fit$layout,
                occasion_of_row = fit$occasion_of_row, n_profiles = fit$n_profiles,
                vars = fit$vars, continuous = fit$continuous,
                categorical = fit$categorical, codes = fit$codes))
  }
  if (identical(fit$covariance_model, "full")) {
    stop(errorCondition("Simulation from a transition fit needs diagonal covariances.",
                        class = "latents_unsupported_inference", call = NULL))
  }
  K <- fit$n_profiles # nolint: object_name_linter. The number of profiles.
  layout <- .multilpa_sequence_layout(fit$group_index, fit$time_values, fit$n_groups,
                                      fit$occasions)
  occasion_of_row <- integer(length(fit$group_index))
  occasion_of_row[layout$slot[!is.na(layout$slot)]] <- col(layout$slot)[!is.na(layout$slot)]
  designs <- .lta_designs(layout, NULL, NULL, "homogeneous")
  initial <- array(vapply(seq_len(fit$n_group_classes), function(h) {
    p <- fit$initial_probabilities[h, ]
    log(p[-K] / p[K])
  }, numeric(K - 1L)), c(1L, K - 1L, fit$n_group_classes))
  transition <- array(vapply(seq_len(fit$n_group_classes), function(h) {
    vapply(seq_len(K), function(k) {
      row <- fit$transition_probabilities[k, , h]
      log(row[setdiff(seq_len(K), k)] / row[k])
    }, numeric(K - 1L))
  }, numeric((K - 1L) * K)), c(1L, K - 1L, K, fit$n_group_classes))
  block <- list(means = fit$means, variances = fit$variances)
  if (!is.null(fit$response_probabilities)) {
    block$response_probabilities <- fit$response_probabilities
  }
  list(parameters = list(measurement = list(block), initial = initial,
                         transition = transition, transition2 = NULL,
                         group_probabilities = fit$group_probabilities),
       designs = designs, layout = layout, occasion_of_row = occasion_of_row,
       n_profiles = K, vars = fit$vars, continuous = fit$continuous,
       categorical = fit$categorical, codes = fit$categorical_data)
}

#' Simulate indicators from a fitted transition model, keeping the data's
#' groups, occasions and covariates (the model conditions on them)
#' @param fit A `multilpa_transitions` or `multilpa_lta` fit.
#' @param data The data it was fitted to.
#' @return `data` with the indicator columns replaced.
#' @noRd
.lta_simulate <- function(fit, data) {
  view <- .lta_general_view(fit, data)
  p <- view$parameters
  K <- view$n_profiles # nolint: object_name_linter. The number of profiles.
  layout <- view$layout
  n_groups <- nrow(layout$slot)
  n_classes <- length(p$group_probabilities)
  classes <- sample.int(n_classes, n_groups, replace = TRUE,
                        prob = p$group_probabilities)
  states <- matrix(NA_integer_, n_groups, layout$n_occasions)
  # Each occasion's state depends on the previous one(s); the draw is
  # vectorized over groups and sequential over occasions.
  invisible(lapply(seq_len(n_classes), function(h) {
    members <- which(classes == h)
    if (length(members) == 0L) return(NULL)
    initial <- exp(.lta_log_initial(view$designs$initial[members, , drop = FALSE],
                                    matrix(p$initial[, , h], ncol(view$designs$initial))))
    states[members, 1L] <<- .multilpa_draw_rows(initial)
    lapply(seq_len(layout$n_occasions)[-1L], function(t) {
      if (isTRUE(p$stayer[h])) {
        states[members, t] <<- states[members, t - 1L]
        return(NULL)
      }
      second <- !is.null(p$transition2) && t >= 3L
      if (second) {
        log_p <- .lta_log_transition2(
          view$designs$transition[[t - 1L]][members, , drop = FALSE],
          array(p$transition2[, , , h], dim(p$transition2)[1:3]))
        pair <- (states[members, t - 2L] - 1L) * K + states[members, t - 1L]
        rows <- t(vapply(seq_along(members), function(i) log_p[i, pair[i], ],
                         numeric(K)))
      } else {
        log_p <- .lta_log_transition(
          view$designs$transition[[t - 1L]][members, , drop = FALSE],
          array(p$transition[, , , h], dim(p$transition)[1:3]))
        origin <- states[members, t - 1L]
        rows <- t(vapply(seq_along(members), function(i) log_p[i, origin[i], ],
                         numeric(K)))
      }
      states[members, t] <<- .multilpa_draw_rows(matrix(exp(rows), length(members)))
    })
  }))
  profile <- integer(nrow(data))
  observed <- !is.na(layout$slot)
  profile[layout$slot[observed]] <- states[observed]
  occasion <- view$occasion_of_row
  simulated <- data
  invisible(lapply(seq_along(p$measurement), function(b) {
    rows <- if (length(p$measurement) == 1L) seq_len(nrow(data)) else which(occasion == b)
    block <- p$measurement[[b]]
    lapply(view$continuous, function(v) {
      simulated[[v]][rows] <<- stats::rnorm(length(rows), block$means[profile[rows], v],
                                            sqrt(block$variances[profile[rows], v]))
    })
    lapply(seq_along(view$categorical), function(i) {
      v <- view$categorical[i]
      codes <- view$codes[, i]
      # Each code's value as it appears in the data, in the data's type.
      value_of_code <- data[[v]][match(seq_len(max(codes, na.rm = TRUE)), codes)]
      drawn <- .multilpa_draw_rows(
        block$response_probabilities[[i]][profile[rows], , drop = FALSE])
      simulated[[v]][rows] <<- value_of_code[drawn]
    })
  }))
  # Missing values stay where the data had them (FIML fits).
  .multilpa_carry_missingness(simulated, data, view$vars)
}

#' Is one transition model nested in another?
#'
#' Every option may only be relaxed from the null to the alternative: at most
#' as many profiles and group classes; homogeneous inside occasion-varying
#' transitions and invariant inside occasion-specific measurement; first
#' order inside second; covariate sets included; equal inside varying
#' variances.
#' @noRd
.lta_nested <- function(null, alternative) {
  a <- .lta_arguments(null)
  b <- .lta_arguments(alternative)
  rank <- function(value, order) match(value, order)
  identical(a$vars, b$vars) && identical(a$id, b$id) && identical(a$time, b$time) &&
    a$n_profiles <= b$n_profiles && a$n_group_classes <= b$n_group_classes &&
    rank(a$transitions, c("homogeneous", "occasion")) <=
      rank(b$transitions, c("homogeneous", "occasion")) &&
    rank(a$measurement, c("invariant", "occasion")) <=
      rank(b$measurement, c("invariant", "occasion")) &&
    a$order <= b$order &&
    isTRUE(a$mover_stayer) <= isTRUE(b$mover_stayer) &&
    all(a$transition_covariates %in% b$transition_covariates) &&
    all(a$initial_covariates %in% b$initial_covariates) &&
    rank(a$variance_model, c("equal", "varying")) <=
      rank(b$variance_model, c("equal", "varying"))
}

#' Refit a transition model's specification to (simulated) data
#' @noRd
.lta_refit <- function(fit, data, n_starts, max_iter, tol) {
  a <- .lta_arguments(fit)
  lta(data, a$vars, a$id, n_profiles = a$n_profiles, time = a$time,
      n_group_classes = a$n_group_classes, variance_model = a$variance_model,
      categorical = a$categorical %||% character(), occasions = a$occasions,
      transitions = a$transitions, transition_covariates = a$transition_covariates,
      initial_covariates = a$initial_covariates, measurement = a$measurement,
      order = a$order, mover_stayer = isTRUE(a$mover_stayer), model = a$model,
      missing = a$missing %||% "error", covariance_model = a$covariance_model %||% "diagonal",
      n_starts = n_starts, max_iter = max_iter,
      tol = tol, min_variance = a$min_variance, min_probability = a$min_probability)
}

#' Parametric bootstrap likelihood-ratio test between nested transition fits
#' @noRd
.lta_bootstrap_lrt <- function(null_model, alternative_model, data, iter, n_starts,
                               max_iter, tol, seed, call) {
  if (is.null(data)) {
    stop(errorCondition(paste(
      "bootstrap_lrt() for transition models needs `data`: the simulation keeps",
      "its groups, occasions and covariates."), class = "latents_bad_argument",
      call = NULL))
  }
  if (!.lta_nested(null_model, alternative_model)) {
    stop(errorCondition(paste(
      "The null transition model must be nested in the alternative: same",
      "indicators, groups and occasions, and every option at most as free."),
      class = "latents_bad_argument", call = NULL))
  }
  if (!isTRUE(null_model$converged) || !isTRUE(alternative_model$converged)) {
    stop(errorCondition("Both fits must have converged.",
                        class = "latents_no_converge", call = NULL))
  }
  check <- .lta_refit(null_model, data, 1L, 0L, tol)
  if (!isTRUE(all.equal(check$n_observations, null_model$n_observations))) {
    stop(errorCondition("`data` does not reproduce the data these fits were made from.",
                        class = "latents_bad_inference_data", call = NULL))
  }
  .multilpa_check_seed(seed)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)  # nolint: object_name_linter. R's name for the RNG state.
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  reversal_window <- 100 * tol * (1 + abs(null_model$log_likelihood))
  observed <- 2 * (alternative_model$log_likelihood - null_model$log_likelihood)
  if (observed < -reversal_window) {
    stop(errorCondition(sprintf(paste(
      "Alternative has lower likelihood by %.3g; improve its optimization first."),
      -observed), class = "latents_reversed_likelihood", call = NULL))
  }
  observed <- max(0, observed)
  replicates <- do.call(rbind, lapply(seq_len(iter), function(i) {
    warning_text <- character()
    tryCatch(withCallingHandlers({
      simulated <- .lta_simulate(null_model, data)
      models <- lapply(list(null_model, alternative_model), .lta_refit, simulated,
                       n_starts, max_iter, tol)
      statistic <- 2 * (models[[2L]]$log_likelihood - models[[1L]]$log_likelihood)
      valid <- all(vapply(models, function(m) isTRUE(m$converged), logical(1))) &&
        statistic >= -reversal_window
      data.frame(replicate = i, statistic = if (valid) max(0, statistic) else NA_real_,
                 valid = valid,
                 boundary = any(vapply(models, function(m) isTRUE(m$boundary), logical(1))),
                 logits_settled = NA, null_replications = NA_integer_,
                 alternative_replications = NA_integer_,
                 warnings = if (length(warning_text) == 0L) NA_character_ else
                   paste(unique(warning_text), collapse = "; "),
                 error = if (valid) NA_character_ else "Nonconvergence or reversed likelihood")
    }, warning = function(warning) {
      warning_text <<- c(warning_text, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(error) {
      data.frame(replicate = i, statistic = NA_real_, valid = FALSE, boundary = NA,
                 logits_settled = NA, null_replications = NA_integer_,
                 alternative_replications = NA_integer_,
                 warnings = paste(unique(warning_text), collapse = "; "),
                 error = conditionMessage(error))
    })
  }))
  valid <- all(replicates$valid)
  p_value <- if (valid) (1 + sum(replicates$statistic >= observed)) / (iter + 1) else NA_real_
  if (!valid) {
    warning(warningCondition(paste(
      "Some bootstrap fits failed validation, so p_value is NA.",
      "summary() reports every replicate; improve fitting and rerun."),
      class = "latents_failed_replicates"))
  }
  describe <- function(fit) {
    a <- .lta_arguments(fit)
    paste(c(sprintf("lta with %d profiles", a$n_profiles),
            if (isTRUE(a$mover_stayer)) "stayer class",
            if (a$order >= 2L) "second order",
            if (!identical(a$transitions, "homogeneous")) "occasion transitions",
            if (length(a$transition_covariates)) paste("transitions on",
                                                        paste(a$transition_covariates, collapse = "+")),
            if (length(a$initial_covariates)) paste("initial on",
                                                     paste(a$initial_covariates, collapse = "+")),
            if (!identical(a$measurement, "invariant")) "occasion measurement"),
          collapse = ", ")
  }
  result <- list(statistic = observed, p_value = p_value,
                 monte_carlo_se = if (valid) sqrt(p_value * (1 - p_value) / (iter + 1)) else NA_real_,
                 iter = iter, n_valid = sum(replicates$valid), replicates = replicates,
                 fixed = character(),
                 null_profiles = .lta_arguments(null_model)$n_profiles,
                 null_group_classes = .lta_arguments(null_model)$n_group_classes,
                 alternative_profiles = .lta_arguments(alternative_model)$n_profiles,
                 alternative_group_classes = .lta_arguments(alternative_model)$n_group_classes,
                 null_family = describe(null_model),
                 alternative_family = describe(alternative_model), call = call)
  class(result) <- "multilpa_bootstrap_lrt"
  result
}

#' Enumerate transition models over numbers of profiles and group classes
#' @noRd
.lta_enumerate <- function(data, vars, id, time, n_profiles, n_group_classes,
                           model, seed, extra, call) {
  unknown <- setdiff(names(extra), setdiff(names(formals(lta)),
                                           c("data", "vars", "id", "time", "n_profiles",
                                             "n_group_classes", "seed", "model")))
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf("%s %s not an argument of lta().",
                                paste(sprintf("`%s`", unknown), collapse = ", "),
                                if (length(unknown) == 1L) "is" else "are"),
                        class = "latents_bad_argument", call = NULL))
  }
  if (any(n_profiles < 2)) {
    stop(errorCondition("Transition models need at least two profiles.",
                        class = "latents_bad_transition", call = NULL))
  }
  models <- if (is.null(model)) NA_character_ else model
  grid <- expand.grid(n_profiles = unique(n_profiles),
                      n_group_classes = unique(n_group_classes), model = models,
                      stringsAsFactors = FALSE)
  runs <- lapply(seq_len(nrow(grid)), function(i) {
    warnings <- character()
    error_text <- NA_character_
    fit <- tryCatch(withCallingHandlers(do.call(lta, c(
      list(data = data, vars = vars, id = id, time = time,
           n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i],
           seed = seed, model = if (is.na(grid$model[i])) NULL else grid$model[i]),
      extra)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }), error = function(error) {
        error_text <<- conditionMessage(error)
        NULL
      })
    entropy <- if (is.null(fit)) NA_real_ else
      .multilpa_relative_entropy(fit$subject_posteriors)
    row <- data.frame(
      n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i],
      model = grid$model[i],
      log_likelihood = if (is.null(fit)) NA_real_ else fit$log_likelihood,
      n_parameters = if (is.null(fit)) NA_integer_ else as.integer(fit$n_parameters),
      aic = if (is.null(fit)) NA_real_ else fit$aic,
      bic = if (is.null(fit)) NA_real_ else fit$bic,
      bic_individual = if (is.null(fit)) NA_real_ else fit$bic_individual,
      entropy = entropy,
      converged = !is.null(fit) && isTRUE(fit$converged),
      boundary = if (is.null(fit)) NA else isTRUE(fit$boundary),
      n_best_replicated = if (is.null(fit)) NA_integer_ else fit$n_best_replicated,
      warnings = paste(unique(warnings), collapse = "; "), error = error_text,
      stringsAsFactors = FALSE)
    list(fit = fit, row = row)
  })
  table <- do.call(rbind, lapply(runs, `[[`, "row"))
  rownames(table) <- NULL
  structure(list(table = table, fits = lapply(runs, `[[`, "fit"), call = call),
            class = "latents_transition_enumeration")
}

#' Tables of a transition-model enumeration
#'
#' The result of `enumerate_classes(..., time = )`: one transition model per
#' number of profiles and group classes (and covariance structure, when
#' `model` names several).
#'
#' @param x A `latents_transition_enumeration` result.
#' @param what `"candidates"` (one row per candidate) or `"best"` (the lowest
#'   value of `criterion`).
#' @param criterion `"bic"` (penalized by the number of groups), `"aic"` or
#'   `"bic_individual"`.
#' @param ... Unused.
#' @return A base `data.frame`, one row per candidate: `n_profiles`,
#'   `n_group_classes`, `model`, `log_likelihood`, `n_parameters`, `aic`,
#'   `bic`, `bic_individual`, relative profile `entropy`, `converged`,
#'   `boundary`, `n_best_replicated`, `warnings` and `error`.
#' @examples
#' grid <- enumerate_classes(course_engagement, c("browse", "lectures"),
#'                           "student", time = "sequence", n_profiles = 2:3,
#'                           n_group_classes = 1, n_starts = 2, seed = 1)
#' get_results(grid)
#' get_results(grid, "best")
#' @export
get_results.latents_transition_enumeration <- function(
    x, what = c("candidates", "best"),
    criterion = c("bic", "aic", "bic_individual"), ...) {
  what <- match.arg(what)
  criterion <- match.arg(criterion)
  if (identical(what, "candidates")) return(x$table)
  usable <- which(is.finite(x$table[[criterion]]))
  if (length(usable) == 0L) {
    stop(errorCondition("No candidate was fitted.", class = "latents_failed_candidate",
                        call = NULL))
  }
  x$table[usable[which.min(x$table[[criterion]][usable])], , drop = FALSE]
}

#' @rdname get_results.latents_transition_enumeration
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_transition_enumeration <- function(x, row.names = NULL,
                                                         optional = FALSE, ...) {
  x$table
}

#' @rdname get_results.latents_transition_enumeration
#' @export
print.latents_transition_enumeration <- function(x, ...) {
  cat(sprintf("Transition-model enumeration: %d candidates\n", nrow(x$table)))
  print(x$table[, c("n_profiles", "n_group_classes", "model", "log_likelihood",
                    "n_parameters", "bic", "aic", "entropy", "converged", "boundary")],
        digits = 4, row.names = FALSE)
  cat("Pick one with candidate_fit(x, n_profiles = , n_group_classes = ).\n")
  invisible(x)
}

#' @rdname get_results.latents_transition_enumeration
#' @param object A `latents_transition_enumeration` result.
#' @export
summary.latents_transition_enumeration <- function(object, ...) object$table

#' Plot information criteria across transition models
#'
#' The criterion against the number of profiles, one line per number of group
#' classes (and structure), told apart by colour, shape and line type.
#' @param x A `latents_transition_enumeration` result.
#' @param criterion `"bic"`, `"aic"` or `"bic_individual"`.
#' @param main,subtitle Optional title and subtitle.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @export
plot.latents_transition_enumeration <- function(x, criterion = c("bic", "aic",
                                                                 "bic_individual"),
                                                main = NULL, subtitle = NULL, ...) {
  criterion <- match.arg(criterion)
  .gg_require()
  frame <- x$table
  frame$value <- frame[[criterion]]
  frame <- subset(frame, is.finite(value))
  frame$series <- factor(paste0(frame$n_group_classes, " group class",
                                ifelse(frame$n_group_classes == 1L, "", "es"),
                                ifelse(is.na(frame$model), "", paste0(", ", frame$model))))
  series <- levels(frame$series)
  ggplot2::ggplot(frame, ggplot2::aes(x = n_profiles, y = value, colour = series,
                                      fill = series, shape = series,
                                      linetype = series)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.8, colour = "grey20") +
    ggplot2::scale_colour_manual(values = rep_len(.gg_okabe_ito, length(series))) +
    ggplot2::scale_fill_manual(values = rep_len(.gg_okabe_ito, length(series))) +
    ggplot2::scale_shape_manual(values = rep_len(.gg_shapes, length(series))) +
    ggplot2::scale_linetype_manual(values = rep_len(c("solid", "dashed", "dotdash",
                                                      "dotted", "longdash"),
                                                    length(series))) +
    ggplot2::scale_x_continuous(breaks = sort(unique(frame$n_profiles))) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = "Profiles", y = toupper(criterion)) +
    .gg_titles(main, subtitle, "Transition models compared", "Lower is better", NULL)
}
