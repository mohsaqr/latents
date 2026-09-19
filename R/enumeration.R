#' Enumerate numbers of individual profiles and group classes
#'
#' Fits every requested combination and retains errors and warnings alongside
#' successful fits. Information criteria are descriptive: neither the smallest
#' BIC nor high entropy guarantees the correct number of classes. Unconverged,
#' boundary, and unreplicated fits are reported, not silently selected.
#' @param data Data frame.
#' @param indicators Continuous indicator names.
#' @param group Group identifier column name.
#' @param profiles Positive integer profile counts to try.
#' @param group_classes Positive integer group-class counts to try.
#' @param seed Optional reproducible seed for each fit.
#' @param ... Further arguments to [multilpa()].
#' @return A list with `table` (criteria/diagnostics) and `fits` (models or NULL).
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", profiles = 1:2,
#'                               group_classes = 1, n_starts = 2, seed = 1)
#' candidates$table
#' @export
enumerate_classes <- function(data, indicators, group, profiles = 1:4,
                              group_classes = 1:3, seed = NULL, ...) {
  stopifnot(is.data.frame(data), is.character(indicators), is.character(group),
            is.numeric(profiles), length(profiles) > 0L,
            all(is.finite(profiles)), all(profiles >= 1), all(profiles == as.integer(profiles)),
            is.numeric(group_classes), length(group_classes) > 0L,
            all(is.finite(group_classes)), all(group_classes >= 1),
            all(group_classes == as.integer(group_classes)))
  extra <- list(...)
  if (any(names(extra) %in% c("n_profiles", "n_group_classes", "data", "indicators", "group", "seed"))) {
    stop("Specify model counts through profiles and group_classes.")
  }
  grid <- expand.grid(n_profiles = unique(profiles), n_group_classes = unique(group_classes))
  runs <- lapply(seq_len(nrow(grid)), function(i) {
    warnings <- character()
    error_text <- NA_character_
    fit <- tryCatch(withCallingHandlers(do.call(multilpa,
      c(list(data = data, indicators = indicators, group = group,
             n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i], seed = seed), extra)),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
      }), error = function(error) {
        error_text <<- conditionMessage(error)
        NULL
      })
    row <- cbind(
      data.frame(n_profiles = grid$n_profiles[i], n_group_classes = grid$n_group_classes[i],
        log_likelihood = if (is.null(fit)) NA_real_ else fit$log_likelihood,
        n_parameters = if (is.null(fit)) NA_integer_ else fit$n_parameters),
      .multilpa_enumeration_indices(fit),
      data.frame(
        profile_entropy = if (is.null(fit)) NA_real_ else
          .multilpa_relative_entropy(fit$subject_posteriors),
        group_entropy = if (is.null(fit)) NA_real_ else
          .multilpa_relative_entropy(fit$group_posteriors),
        converged = !is.null(fit) && fit$converged,
        boundary = if (is.null(fit)) NA else fit$boundary,
        n_best_replicated = if (is.null(fit)) NA_integer_ else fit$n_best_replicated,
        warnings = paste(unique(warnings), collapse = "; "), error = error_text))
    list(fit = fit, row = row)
  })
  result <- list(table = do.call(rbind, lapply(runs, `[[`, "row")),
       fits = lapply(runs, `[[`, "fit"), call = match.call())
  class(result) <- "multilpa_enumeration"
  result
}

#' Generate observations from a fitted discrete multilevel model
#' @param object Fitted model.
#' @return Simulated data frame with the original group layout.
#' @noRd
.multilpa_simulate <- function(object) {
  stopifnot(inherits(object, "multilpa"), !inherits(object, "multilpa_covariates"))
  group_class <- sample.int(object$n_group_classes, object$n_groups,
                            replace = TRUE, prob = object$group_probabilities)
  profile <- .multilpa_draw_rows(
    object$profile_probabilities[group_class[object$group_index], , drop = FALSE])
  continuous <- .multilpa_continuous_names(object)
  dimension <- length(continuous)
  result <- data.frame(row.names = seq_len(object$n_observations))
  if (dimension > 0L) {
    noise <- matrix(stats::rnorm(object$n_observations * dimension),
                    object$n_observations)
    values <- if (identical(object$covariance_model, "full")) {
      factors <- lapply(seq_len(object$n_profiles),
                        function(k) chol(object$covariances[, , k]))
      t(matrix(vapply(seq_len(object$n_observations), function(i) {
        as.vector(noise[i, , drop = FALSE] %*% factors[[profile[i]]]) +
          object$means[profile[i], ]
      }, numeric(dimension)), nrow = dimension, ncol = object$n_observations))
    } else {
      object$means[profile, , drop = FALSE] +
        noise * sqrt(object$variances[profile, , drop = FALSE])
    }
    values <- as.data.frame(values)
    names(values) <- continuous
    result <- cbind(result, values)
  }
  ## Categorical indicators are returned as their original category values, not
  ## as codes, so that refitting re-encodes them to exactly the same levels.
  blocks <- object$response_probabilities
  if (!is.null(blocks) && length(blocks) > 0L) {
    drawn <- lapply(names(blocks), function(indicator) {
      probabilities <- blocks[[indicator]][profile, , drop = FALSE]
      ## Levels are stored as character. Returning them as such would make a
      ## refit sort "10" before "2", permuting the categories relative to the
      ## fit being bootstrapped, so a numeric level set goes back as numeric.
      levels_observed <- object$categorical_levels[[indicator]]
      numeric_levels <- suppressWarnings(as.numeric(levels_observed))
      if (!anyNA(numeric_levels)) levels_observed <- numeric_levels
      levels_observed[.multilpa_draw_rows(probabilities)]
    })
    names(drawn) <- names(blocks)
    result <- cbind(result, as.data.frame(drawn, stringsAsFactors = FALSE))
  }
  result <- result[, object$indicators, drop = FALSE]
  result[[object$group]] <- object$group_values[object$group_index]
  result
}

#' Draw one category per row from a matrix of row-wise probabilities
#'
#' Inverse-CDF sampling in one pass. The last cumulative probability is set to
#' exactly one so that a rounding shortfall cannot leave a draw unmatched and
#' fall back silently to the first category.
#'
#' @param probabilities Matrix whose rows each sum to one.
#' @return An integer vector, one column index per row.
#' @noRd
.multilpa_draw_rows <- function(probabilities) {
  stopifnot("`probabilities` must be a matrix" = is.matrix(probabilities),
            "`probabilities` must have at least one column" = ncol(probabilities) >= 1L)
  cumulative <- t(apply(probabilities, 1L, cumsum))
  if (ncol(probabilities) == 1L) cumulative <- matrix(1, nrow(probabilities), 1L)
  cumulative[, ncol(cumulative)] <- 1
  max.col(cumulative >= stats::runif(nrow(probabilities)), ties.method = "first")
}

#' Parametric bootstrap likelihood-ratio comparison
#'
#' Simulates complete indicators under the null model while preserving observed
#' group sizes, refits both models, and compares their likelihood differences.
#' Models must differ by exactly one individual profile or one group class, with
#' the other count and covariance specification fixed. This is a native
#' parametric bootstrap, not an implementation of Mplus TECH14. It does not use
#' a chi-square reference distribution. Any failed/nonconverged or reversed
#' replicate makes the p-value NA, avoiding silent deletion of difficult fits.
#' @param null_model Smaller, converged [multilpa()] model on complete data.
#' @param alternative_model Larger model fitted to exactly the same data.
#' @param data Original data, used to verify both fitted likelihoods.
#' @param n_boot Number of simulated datasets (at least two; use many for inference).
#' @param n_starts Number of starts for each simulated fit.
#' @param max_iter Maximum EM iterations for each simulated fit.
#' @param tol Relative likelihood convergence tolerance.
#' @param seed Optional seed, with caller RNG state restored.
#' @return Observed statistic, finite-simulation corrected p-value, Monte Carlo
#'   standard error, replicate diagnostics and counts. Inspect failed starts,
#'   boundary flags and likelihood replication before interpreting results.
#' @examples
#' # After fitting nested models on d:
#' # bootstrap_lrt(smaller, larger, d, n_boot = 199, seed = 1)
#' @export
bootstrap_lrt <- function(null_model, alternative_model, data,
                                 n_boot = 199L, n_starts = 10L, max_iter = 1000L,
                                 tol = 1e-8, seed = NULL) {
  stopifnot(inherits(null_model, "multilpa"), inherits(alternative_model, "multilpa"),
            is.data.frame(data), is.numeric(n_boot), length(n_boot) == 1L,
            is.finite(n_boot), n_boot >= 2L, n_boot == as.integer(n_boot),
            is.numeric(n_starts), length(n_starts) == 1L, is.finite(n_starts),
            n_starts >= 1, n_starts == as.integer(n_starts),
            is.numeric(max_iter), length(max_iter) == 1L, is.finite(max_iter),
            max_iter >= 1, max_iter == as.integer(max_iter),
            is.numeric(tol), length(tol) == 1L, is.finite(tol), tol > 0)
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed),
              seed >= 0, seed <= .Machine$integer.max, seed == as.integer(seed))
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
              rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  fields <- c("indicators", "group", "group_values", "group_index", "variance_model", "min_variance", "covariance_model")
  fields <- c(fields, "categorical", "categorical_levels", "min_probability")
  if (!all(vapply(fields, function(field) identical(null_model[[field]], alternative_model[[field]]), logical(1)))) {
    stop("Models must use the same observations, group layout and covariance specification.")
  }
  delta <- c(alternative_model$n_profiles - null_model$n_profiles,
             alternative_model$n_group_classes - null_model$n_group_classes)
  if (!all(delta >= 0) || sum(delta) != 1) stop("Models must differ by exactly one class at one level.")
  invisible(lapply(list(null_model, alternative_model), function(model) {
    if (!isTRUE(model$converged) || isTRUE(model$boundary)) {
      stop("Original models must be converged with inactive variance bounds.")
    }
    if (!all(c(model$indicators, model$group) %in% names(data)) ||
        nrow(data) != model$n_observations ||
        !identical(data[[model$group]], model$group_values[model$group_index])) {
      stop("data must preserve the original observations and group ordering.")
    }
    continuous <- .multilpa_continuous_names(model)
    x <- if (length(continuous) == 0L) matrix(numeric(0), nrow(data), 0L) else
      as.matrix(data[continuous])
    if (!is.numeric(x) || any(!is.finite(x))) stop("Bootstrap currently requires complete finite indicators.")
    if (!is.null(model$indicator_data) && length(continuous) > 0L &&
        !identical(x, model$indicator_data)) {
      stop("data must reproduce the original indicator data and row order.")
    }
    encoded <- if (length(model$categorical %||% character()) == 0L) NULL else
      .multilpa_encode_categorical(data[, model$categorical, drop = FALSE])
    codes <- encoded$codes
    if (anyNA(codes)) stop("Bootstrap currently requires complete finite indicators.")
    if (!is.null(encoded) && !identical(encoded$levels, model$categorical_levels)) {
      stop("data must reproduce the original categorical levels and coding.")
    }
    if (!is.null(codes) && !identical(unname(codes), unname(model$categorical_data))) {
      stop("data must reproduce the original categorical indicators and row order.")
    }
    likelihood <- .multilpa_expectation(x, model$group_index, model,
                                        codes)$log_likelihood
    if (abs(likelihood - model$log_likelihood) > 1e-7 * (1 + abs(likelihood))) {
      stop("data do not reproduce the fitted model likelihood.")
    }
  }))
  observed <- 2 * (alternative_model$log_likelihood - null_model$log_likelihood)
  if (observed < -1e-5) stop("Alternative has lower likelihood; improve its optimization first.")
  observed <- max(0, observed)
  replicates <- do.call(rbind, lapply(seq_len(n_boot), function(i) {
    warning_text <- character()
    tryCatch(withCallingHandlers({
      simulated <- .multilpa_simulate(null_model)
      models <- lapply(list(null_model, alternative_model), function(model) {
        multilpa(simulated, model$indicators, model$group, model$n_profiles,
                    model$n_group_classes, model$variance_model, n_starts = n_starts,
                    max_iter = max_iter, tol = tol, min_variance = model$min_variance,
                    covariance_model = if (is.null(model$covariance_model)) "diagonal" else model$covariance_model,
                    categorical = model$categorical %||% character(),
                    min_probability = model$min_probability %||% 1e-10)
      })
      statistic <- 2 * (models[[2L]]$log_likelihood - models[[1L]]$log_likelihood)
      valid <- all(vapply(models, `[[`, logical(1), "converged")) && statistic >= -1e-5
      data.frame(replicate = i, statistic = if (valid) max(0, statistic) else NA_real_,
        valid = valid, boundary = any(vapply(models, `[[`, logical(1), "boundary")),
        null_replications = models[[1L]]$n_best_replicated,
        alternative_replications = models[[2L]]$n_best_replicated,
        warnings = paste(unique(warning_text), collapse = "; "),
        error = if (valid) NA_character_ else "Nonconvergence or reversed likelihood")
    }, warning = function(warning) {
      warning_text <<- c(warning_text, conditionMessage(warning))
    }), error = function(error) {
      data.frame(replicate = i, statistic = NA_real_, valid = FALSE, boundary = NA,
        null_replications = NA_integer_, alternative_replications = NA_integer_,
        warnings = paste(unique(warning_text), collapse = "; "), error = conditionMessage(error))
    })
  }))
  valid <- all(replicates$valid)
  p_value <- if (valid) (1 + sum(replicates$statistic >= observed)) / (n_boot + 1) else NA_real_
  if (!valid) warning("Some bootstrap fits failed validation; p_value is NA. Inspect $replicates and improve fitting.")
  list(statistic = observed, p_value = p_value,
       monte_carlo_se = if (valid) sqrt(p_value * (1 - p_value) / (n_boot + 1)) else NA_real_,
       n_boot = n_boot, n_valid = sum(replicates$valid), replicates = replicates,
       call = match.call())
}

#' Spread the information criteria of one candidate into a single row
#' @param fit A fitted `multilpa` model, or `NULL` for a failed candidate.
#' @return A one-row `data.frame` of criteria, all `NA_real_` when `fit` is `NULL`.
#' @noRd
.multilpa_enumeration_indices <- function(fit) {
  names_wanted <- c("aic", "bic_groups", "bic_individual", "sabic_groups",
                    "sabic_individual", "caic_groups", "caic_individual",
                    "awe_groups", "awe_individual", "icl_groups", "icl_individual")
  if (is.null(fit)) {
    return(as.data.frame(stats::setNames(
      rep(list(NA_real_), length(names_wanted)), names_wanted)))
  }
  indices <- information_criteria(fit)
  labels <- ifelse(indices$convention == "none", indices$criterion,
                   paste(indices$criterion,
                         sub("individuals", "individual", indices$convention),
                         sep = "_"))
  values <- stats::setNames(indices$value, labels)
  as.data.frame(as.list(values[names_wanted]))
}
