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
  if (identical(fit$covariance_model, "full")) block$covariances <- fit$covariances
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
    drawn <- .lta_draw_continuous(block, profile[rows], view$continuous)
    lapply(view$continuous, function(v) simulated[[v]][rows] <<- drawn[, v])
    lapply(seq_along(view$categorical), function(i) {
      v <- view$categorical[i]
      codes <- view$codes[, i]
      # Each code's value as it appears in the data, in the data's type.
      value_of_code <- data[[v]][match(seq_len(max(codes, na.rm = TRUE)), codes)]
      drawn <- .multilpa_draw_rows(
        block$response_probabilities[[i]][profile[rows], , drop = FALSE])
      simulated[[v]][rows] <<- value_of_code[drawn]
    })
    if (!is.null(fit$extra_data)) {
      drawn <- .latents_draw_extra(
        c(block, list(extra_data = .latents_extra_rows(fit$extra_data, rows))),
        profile[rows])
      lapply(names(drawn), function(v) {
        simulated[[v]][rows] <<- .latents_as_column_type(drawn[[v]], data[[v]])
      })
    }
  }))
  # Missing values stay where the data had them (FIML fits).
  .multilpa_carry_missingness(simulated, data, view$vars)
}

#' Draw continuous indicators given each row's profile
#'
#' Diagonal fits draw each indicator independently, in the order the package
#' always has, so seeded results are unchanged. A fit that carries full
#' covariance matrices draws each profile's rows as the profile mean plus
#' standard normals times the upper Cholesky factor of its covariance.
#'
#' @param block A measurement block: `means` and `variances` (profiles by
#'   indicators) and, for a full-covariance fit, `covariances` (indicators by
#'   indicators by profiles).
#' @param profile The profile of every row to draw.
#' @param continuous The continuous indicator names.
#' @return A rows-by-indicators matrix with those column names.
#' @noRd
.lta_draw_continuous <- function(block, profile, continuous) {
  n <- length(profile)
  if (is.null(block$covariances)) {
    return(vapply(continuous, function(v) {
      stats::rnorm(n, block$means[profile, v], sqrt(block$variances[profile, v]))
    }, numeric(n)) |> matrix(n, length(continuous), dimnames = list(NULL, continuous)))
  }
  noise <- matrix(stats::rnorm(n * length(continuous)), n, length(continuous))
  correlated <- Reduce(function(drawn, k) {
    rows <- which(profile == k)
    if (length(rows) == 0L) return(drawn)
    # General fits store the array unnamed, in `continuous` order.
    covariance <- block$covariances[, , k]
    if (!is.null(rownames(covariance))) covariance <- covariance[continuous, continuous]
    factor <- chol(covariance)
    drawn[rows, ] <- noise[rows, , drop = FALSE] %*% factor
    drawn
  }, sort(unique(profile)), noise)
  structure(correlated + block$means[profile, continuous, drop = FALSE],
            dimnames = list(NULL, continuous))
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
  structure_for <- function(a) a$model %||%
    .multilpa_resolve_structure(a$variance_model, a$covariance_model %||% "diagonal")
  covariance_nested <- all(vapply(seq_len(3L), function(position) {
    rank(substr(structure_for(a), position, position), c("I", "E", "V")) <=
      rank(substr(structure_for(b), position, position), c("I", "E", "V"))
  }, logical(1)))
  identical(a$vars, b$vars) && identical(a$id, b$id) && identical(a$time, b$time) &&
    identical(a$categorical %||% character(), b$categorical %||% character()) &&
    identical(a$occasions %||% "observed", b$occasions %||% "observed") &&
    identical(a$missing %||% "error", b$missing %||% "error") && covariance_nested &&
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
      rank(b$variance_model, c("equal", "varying")) &&
    identical(a$ordinal %||% character(), b$ordinal %||% character()) &&
    identical(a[["count"]] %||% character(), b[["count"]] %||% character()) &&
    identical(a$count_model %||% "poisson", b$count_model %||% "poisson") &&
    rank(a$count_dispersion %||% "varying", c("equal", "varying")) <=
      rank(b$count_dispersion %||% "varying", c("equal", "varying")) &&
    identical(null$extra_data$ordinal_levels, alternative$extra_data$ordinal_levels)
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
      ordinal = a$ordinal %||% character(), count = a$count %||% character(),
      count_model = a$count_model %||% "poisson",
      count_dispersion = a$count_dispersion %||% "varying", weights = fit$weights,
      n_starts = n_starts, max_iter = max_iter,
      tol = tol, min_variance = a$min_variance, min_probability = a$min_probability)
}

#' Verify original indicators, occasions and covariate designs of an LTA fit
#' @param fit A fitted transition model.
#' @param data Its original data frame.
#' @return `NULL`, invisibly; refuses mismatched fitting data.
#' @noRd
.lta_check_data <- function(fit, data) {
  a <- .lta_arguments(fit)
  required <- c(a$vars, a$id, a$time, a$transition_covariates, a$initial_covariates)
  refuse <- function() stop(errorCondition(
    "`data` must reproduce the indicators, groups, occasions and covariates of both fits.",
    class = "latents_bad_inference_data", call = NULL))
  if (!is.data.frame(data) || nrow(data) != fit$n_observations ||
      !all(required %in% names(data))) refuse()
  if (!inherits(fit, "multilpa_lta")) {
    .multilpa_check_transition_data(fit, data)
    if (length(fit$categorical) &&
        !identical(.multilpa_encode_categorical(data[fit$categorical])$levels,
                   fit$categorical_levels)) refuse()
    if (!identical(.multilpa_time_values(data, a$time, a$id, a$vars),
                   fit$time_values)) refuse()
    return(invisible(NULL))
  }
  groups <- .multilpa_prepare_groups(data[[a$id]])
  time_values <- .multilpa_time_values(data, a$time, a$id, a$vars)
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, a$occasions)
  designs <- .lta_designs(layout, .lta_covariate_matrix(data, a$transition_covariates),
                         .lta_covariate_matrix(data, a$initial_covariates), a$transitions)
  x <- sweep(as.matrix(data[fit$continuous]), 2L, fit$centers, "-")
  encoded <- if (!length(fit$categorical)) NULL else
    .multilpa_encode_categorical(data[fit$categorical])
  codes <- encoded$codes
  extra <- .latents_prepare_extra(data, fit$ordinal, fit$count, a$missing)
  same <- function(left, right) isTRUE(all.equal(left, right, check.attributes = FALSE,
                                                tolerance = 1e-12))
  if (!identical(groups$ids, fit$group_ids) || !identical(groups$index, fit$group_index) ||
      !identical(time_values, fit$time_values) || !identical(layout, fit$layout) ||
      !identical(encoded$levels, fit$categorical_levels) || !same(designs, fit$designs) ||
      !same(x, fit$x) || !same(codes, fit$codes) ||
      !identical(extra$ordinal_levels, fit$extra_data$ordinal_levels) ||
      !same(extra$ordinal, fit$extra_data$ordinal) ||
      !same(extra$count, fit$extra_data$count)) refuse()
  invisible(NULL)
}

#' Parametric bootstrap likelihood-ratio test between nested transition fits
#' @noRd
.lta_bootstrap_lrt <- function(null_model, alternative_model, data, iter, n_starts,
                               max_iter, tol, seed, call) {
  positive_integer <- function(value) is.numeric(value) && length(value) == 1L &&
    is.finite(value) && value >= 1 && value <= .Machine$integer.max && value == floor(value)
  stopifnot("`iter` must be a positive integer" = positive_integer(iter),
            "`n_starts` must be a positive integer" = positive_integer(n_starts),
            "`max_iter` must be a nonnegative integer" = positive_integer(max_iter + 1),
            "`tol` must be a positive finite number" = is.numeric(tol) &&
              length(tol) == 1L && is.finite(tol) && tol > 0)
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
  lapply(list(null_model, alternative_model), .lta_check_data, data = data)
  # Whether the null can be simulated from is a property of the fit, not of a
  # replicate: refuse here, or every replicate fails and the classed refusal
  # surfaces only as a generic failed-replicates warning.
  .lta_general_view(null_model, data)
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
