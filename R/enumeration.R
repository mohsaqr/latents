#' Enumerate numbers of individual profiles and group classes
#'
#' Fits every requested combination and retains errors and warnings alongside
#' successful fits. Information criteria are descriptive: neither the smallest
#' BIC nor high entropy guarantees the correct number of classes. Unconverged,
#' boundary, and unreplicated fits are reported, not silently selected.
#' @param data Data frame.
#' @param vars Continuous indicator names.
#' @param id Group identifier column name.
#' @param n_profiles Positive integer profile counts to try.
#' @param n_group_classes Positive integer group-class counts to try.
#' @param seed Optional reproducible seed for each fit.
#' @param ... Further arguments to [multilpa()].
#' @return An object of class `multilpa_enumeration`. Read it with the verbs
#'   that describe it rather than by reaching into it: [as.data.frame()] gives
#'   one row per candidate model with every criterion and diagnostic,
#'   [summary()] gives one row per information criterion naming the candidate
#'   that minimises it, [plot()] draws one criterion across the grid, and
#'   [candidate_fit()] returns the fitted model for one cell of the grid.
#' @seealso [candidate_fit()] to take one fitted model out of the grid,
#'   [summary.multilpa_enumeration()] for the criterion-by-criterion comparison.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(candidates)
#' summary(candidates)
#' @export
enumerate_classes <- function(data, vars, id, n_profiles = 1:4,
                              n_group_classes = 1:3, seed = NULL, ...) {
  stopifnot(is.data.frame(data), is.character(vars), is.character(id),
            is.numeric(n_profiles), length(n_profiles) > 0L,
            all(is.finite(n_profiles)), all(n_profiles >= 1), all(n_profiles == as.integer(n_profiles)),
            is.numeric(n_group_classes), length(n_group_classes) > 0L,
            all(is.finite(n_group_classes)), all(n_group_classes >= 1),
            all(n_group_classes == as.integer(n_group_classes)))
  # Every name this guard used to reject is now a formal, so R refuses the call
  # with "matched by multiple actual arguments" before `...` is assembled.
  extra <- list(...)
  grid <- expand.grid(n_profiles = unique(n_profiles), n_group_classes = unique(n_group_classes))
  runs <- lapply(seq_len(nrow(grid)), function(i) {
    warnings <- character()
    error_text <- NA_character_
    fit <- tryCatch(withCallingHandlers(do.call(multilpa,
      c(list(data = data, vars = vars, id = id,
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

#' Take one fitted model out of an enumeration grid
#'
#' An enumeration keeps every candidate it fitted. This returns one of them, so
#' that a caller who wants to plot, summarise or test a particular candidate
#' names it by its class counts instead of indexing the grid by position.
#'
#' @param x An `multilpa_enumeration` result from [enumerate_classes()].
#' @param n_profiles Number of individual profiles identifying the candidate.
#' @param n_group_classes Number of group classes identifying the candidate.
#' @return The fitted [multilpa()] model for that cell of the grid.
#' @section Conditions:
#'   `multilpa_unknown_candidate` when the requested class counts are not in
#'   the grid, and `multilpa_failed_candidate` when they are in the grid but
#'   that fit raised an error, so no model exists to return. Neither situation
#'   returns `NULL`, because a `NULL` would flow silently into whatever the
#'   caller did next.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' candidate_fit(candidates, n_profiles = 2, n_group_classes = 1)
#' @export
candidate_fit <- function(x, n_profiles, n_group_classes = 1L) {
  stopifnot(
    "`x` must be an `multilpa_enumeration` result" =
      inherits(x, "multilpa_enumeration"),
    "`n_profiles` must be a single positive integer" =
      is.numeric(n_profiles) && length(n_profiles) == 1L &&
      is.finite(n_profiles) && n_profiles >= 1,
    "`n_group_classes` must be a single positive integer" =
      is.numeric(n_group_classes) && length(n_group_classes) == 1L &&
      is.finite(n_group_classes) && n_group_classes >= 1)
  grid <- x$table
  at <- which(grid$n_profiles == n_profiles &
                grid$n_group_classes == n_group_classes)
  if (length(at) != 1L) {
    stop(errorCondition(sprintf(
      "No candidate with %d profiles and %d group classes was enumerated.",
      as.integer(n_profiles), as.integer(n_group_classes)),
      class = "multilpa_unknown_candidate", call = NULL))
  }
  fit <- x$fits[[at]]
  if (is.null(fit)) {
    stop(errorCondition(sprintf(
      "The candidate with %d profiles and %d group classes could not be fitted: %s",
      as.integer(n_profiles), as.integer(n_group_classes), grid$error[at]),
      class = "multilpa_failed_candidate", call = NULL))
  }
  fit
}

#' Summarise a class-enumeration grid
#'
#' Reports, for every information criterion the grid carries, which candidate
#' minimises it. Criteria disagree by construction: they differ in how they
#' penalise parameters and in whether they count individuals or independent
#' groups. Laying the minima side by side shows that disagreement instead of
#' hiding it behind one default. Nothing here selects a model.
#'
#' @param object An `multilpa_enumeration` result from [enumerate_classes()].
#' @param ... Reserved for compatibility with `summary()`.
#' @return An object of class `summary_multilpa_enumeration`, with a `print`
#'   method and an [as.data.frame()] accessor. `as.data.frame()` returns one
#'   row per information criterion, with columns `criterion`, `convention`
#'   (`"groups"`, `"individuals"` or `NA` when the criterion uses no sample
#'   size), `n_profiles` and `n_group_classes` of the minimising candidate, and
#'   its `value`. Only converged candidates are eligible; a criterion with no
#'   converged candidate has `NA` in every remaining column. Pass
#'   `what = "candidates"` for the full grid instead.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' summary(candidates)
#' as.data.frame(summary(candidates))
#' @export
summary.multilpa_enumeration <- function(object, ...) {
  stopifnot("`object` must be an `multilpa_enumeration` result" =
              inherits(object, "multilpa_enumeration"))
  grid <- object$table
  result <- list(criteria = .multilpa_enumeration_minima(grid),
                 candidates = grid,
                 n_candidates = nrow(grid),
                 n_converged = sum(grid$converged),
                 n_failed = sum(is.na(grid$log_likelihood)),
                 n_boundary = sum(grid$boundary %in% TRUE),
                 call = object$call)
  class(result) <- "summary_multilpa_enumeration"
  result
}

#' Which candidate minimises each information criterion
#' @param grid The enumeration table.
#' @return One row per criterion carried by the grid.
#' @noRd
.multilpa_enumeration_minima <- function(grid) {
  present <- intersect(.multilpa_enumeration_criteria(), names(grid))
  eligible <- grid$converged %in% TRUE
  rows <- lapply(present, function(name) {
    values <- grid[[name]]
    values[!eligible] <- NA_real_
    best <- if (all(is.na(values))) NA_integer_ else which.min(values)
    data.frame(
      criterion = sub("_(groups|individual)$", "", name),
      convention = if (grepl("_groups$", name)) "groups" else
        if (grepl("_individual$", name)) "individuals" else NA_character_,
      n_profiles = if (is.na(best)) NA_integer_ else
        as.integer(grid$n_profiles[best]),
      n_group_classes = if (is.na(best)) NA_integer_ else
        as.integer(grid$n_group_classes[best]),
      value = if (is.na(best)) NA_real_ else values[best],
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Print an enumeration summary
#' @param x A `summary_multilpa_enumeration` object.
#' @param digits Number of printed significant digits.
#' @param ... Passed to the underlying `data.frame` printing.
#' @return The summary, invisibly.
#' @examples
#' # After enumerating: print(summary(candidates))
#' @export
print.summary_multilpa_enumeration <- function(x, digits = 4L, ...) {
  stopifnot("`x` must be a `summary_multilpa_enumeration` object" =
              inherits(x, "summary_multilpa_enumeration"),
            "`digits` must be a single number between 1 and 22" =
              is.numeric(digits) && length(digits) == 1L && is.finite(digits) &&
              digits >= 1 && digits <= 22)
  cat(sprintf("Class enumeration: %d candidates, %d converged, %d failed to fit\n",
              x$n_candidates, x$n_converged, x$n_failed))
  cat("\nCandidate minimising each information criterion:\n")
  print(x$criteria, digits = digits, row.names = FALSE, ...)
  disagreement <- unique(x$criteria[stats::complete.cases(
    x$criteria[c("n_profiles", "n_group_classes")]), c("n_profiles", "n_group_classes")])
  cat(sprintf("\n%d distinct candidate(s) are minimal under some criterion.\n",
              nrow(disagreement)))
  cat("No candidate is selected automatically. Choose one convention and keep it.\n")
  invisible(x)
}

#' Tidy an enumeration summary
#' @param x A `summary_multilpa_enumeration` object.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what `"criteria"` returns one row per information criterion with the
#'   candidate that minimises it; `"candidates"` returns the full grid, exactly
#'   as [as.data.frame.multilpa_enumeration()] returns it.
#' @param ... Ignored.
#' @return A base `data.frame`. For `"criteria"`, one row per information
#'   criterion, with columns `criterion`, `convention`, `n_profiles`,
#'   `n_group_classes` and `value`. For `"candidates"`, one row per candidate
#'   model in the grid.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
#' candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
#'                               n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(summary(candidates))
#' @export
as.data.frame.summary_multilpa_enumeration <- function(x, row.names = NULL,
                                                       optional = FALSE,
                                                       what = c("criteria",
                                                                "candidates"),
                                                       ...) {
  stopifnot("`x` must be a `summary_multilpa_enumeration` object" =
              inherits(x, "summary_multilpa_enumeration"))
  what <- match.arg(what)
  result <- switch(what, criteria = x$criteria, candidates = x$candidates)
  row.names(result) <- row.names
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
  result <- result[, object$vars, drop = FALSE]
  result[[object$id]] <- object$group_values[object$group_index]
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
#' @param data Optional. The data both models were fitted to, used to verify
#'   both fitted likelihoods; when omitted it is rebuilt from what the null
#'   model stores.
#' @param iter Number of simulated datasets (at least two; use many for inference).
#' @param n_starts Number of starts for each simulated fit.
#' @param max_iter Maximum EM iterations for each simulated fit.
#' @param tol Relative likelihood convergence tolerance.
#' @param seed Optional seed, with caller RNG state restored.
#' @return An object of class `multilpa_bootstrap_lrt`, carrying the observed
#'   statistic, the finite-simulation corrected p-value, its Monte Carlo
#'   standard error, and one record per replicate. Read it with the verbs that
#'   describe it: [as.data.frame()] gives one row per simulated replicate,
#'   `as.data.frame(what = "test")` the single-row test result, [summary()] the
#'   test beside the replicate diagnostics, and [plot()] the simulated null
#'   distribution with the observed statistic marked. Inspect failed starts,
#'   boundary flags and likelihood replication before interpreting results.
#' @examples
#' # After fitting nested models on d:
#' # bootstrap_lrt(smaller, larger, d, iter = 199, seed = 1)
#' @export
bootstrap_lrt <- function(null_model, alternative_model, data = NULL,
                                 iter = 199L, n_starts = 10L, max_iter = 1000L,
                                 tol = 1e-8, seed = NULL) {
  stopifnot(inherits(null_model, "multilpa"), inherits(alternative_model, "multilpa"))
  # The null model is the one being simulated from, so its own columns are the
  # ones that matter when the caller does not supply data.
  data <- .multilpa_resolve_data(null_model, data)
  stopifnot(is.data.frame(data), is.numeric(iter), length(iter) == 1L,
            is.finite(iter), iter >= 2L, iter == as.integer(iter),
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
  fields <- c("vars", "id", "group_values", "group_index", "variance_model", "min_variance", "covariance_model")
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
    if (!all(c(model$vars, model$id) %in% names(data)) ||
        nrow(data) != model$n_observations ||
        !identical(data[[model$id]], model$group_values[model$group_index])) {
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
  replicates <- do.call(rbind, lapply(seq_len(iter), function(i) {
    warning_text <- character()
    tryCatch(withCallingHandlers({
      simulated <- .multilpa_simulate(null_model)
      models <- lapply(list(null_model, alternative_model), function(model) {
        multilpa(simulated, model$vars, model$id, model$n_profiles,
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
  p_value <- if (valid) (1 + sum(replicates$statistic >= observed)) / (iter + 1) else NA_real_
  if (!valid) {
    warning(warningCondition(paste(
      "Some bootstrap fits failed validation, so p_value is NA.",
      "summary() reports every replicate; improve fitting and rerun."),
      class = "multilpa_failed_replicates"))
  }
  result <- list(statistic = observed, p_value = p_value,
       monte_carlo_se = if (valid) sqrt(p_value * (1 - p_value) / (iter + 1)) else NA_real_,
       iter = iter, n_valid = sum(replicates$valid), replicates = replicates,
       null_profiles = null_model$n_profiles,
       null_group_classes = null_model$n_group_classes,
       alternative_profiles = alternative_model$n_profiles,
       alternative_group_classes = alternative_model$n_group_classes,
       call = match.call())
  class(result) <- "multilpa_bootstrap_lrt"
  result
}

#' Print a parametric bootstrap likelihood-ratio comparison
#' @param x An `multilpa_bootstrap_lrt` result.
#' @param ... Reserved for compatibility with `print()`.
#' @return The input, invisibly.
#' @examples
#' # After bootstrapping: print(comparison)
#' @export
print.multilpa_bootstrap_lrt <- function(x, ...) {
  stopifnot("`x` must be an `multilpa_bootstrap_lrt` result" =
              inherits(x, "multilpa_bootstrap_lrt"))
  cat(sprintf("Parametric bootstrap likelihood-ratio comparison\n"))
  cat(sprintf("Null: %d profiles, %d group classes; alternative: %d profiles, %d group classes\n",
              x$null_profiles, x$null_group_classes,
              x$alternative_profiles, x$alternative_group_classes))
  cat(sprintf("Observed statistic: %.6f\n", x$statistic))
  cat(sprintf("p-value: %s (Monte Carlo SE %s) from %d of %d valid replicates\n",
              format(x$p_value, digits = 4L), format(x$monte_carlo_se, digits = 3L),
              x$n_valid, x$iter))
  if (is.na(x$p_value)) {
    cat("The p-value is withheld because not every replicate was valid; summary() lists them.\n")
  }
  invisible(x)
}

#' Summarise a parametric bootstrap likelihood-ratio comparison
#' @param object An `multilpa_bootstrap_lrt` result.
#' @param ... Reserved for compatibility with `summary()`.
#' @return An object of class `summary_multilpa_bootstrap_lrt`, with a `print`
#'   method and an [as.data.frame()] accessor. `as.data.frame()` returns the
#'   one-row test table by default and one row per replicate with
#'   `what = "replicates"`.
#' @examples
#' # After bootstrapping: summary(comparison)
#' @export
summary.multilpa_bootstrap_lrt <- function(object, ...) {
  stopifnot("`object` must be an `multilpa_bootstrap_lrt` result" =
              inherits(object, "multilpa_bootstrap_lrt"))
  result <- list(test = .multilpa_bootstrap_test_frame(object),
                 replicates = object$replicates,
                 n_boundary = sum(object$replicates$boundary %in% TRUE),
                 n_errors = sum(!is.na(object$replicates$error)),
                 call = object$call)
  class(result) <- "summary_multilpa_bootstrap_lrt"
  result
}

#' The bootstrap test as a single tidy row
#' @param x An `multilpa_bootstrap_lrt` result.
#' @return A one-row `data.frame` describing the comparison.
#' @noRd
.multilpa_bootstrap_test_frame <- function(x) {
  data.frame(null_profiles = x$null_profiles,
             null_group_classes = x$null_group_classes,
             alternative_profiles = x$alternative_profiles,
             alternative_group_classes = x$alternative_group_classes,
             statistic = x$statistic, p_value = x$p_value,
             monte_carlo_se = x$monte_carlo_se,
             iter = x$iter, n_valid = x$n_valid,
             row.names = NULL)
}

#' Print a bootstrap likelihood-ratio summary
#' @param x A `summary_multilpa_bootstrap_lrt` object.
#' @param digits Number of printed significant digits.
#' @param ... Passed to the underlying `data.frame` printing.
#' @return The summary, invisibly.
#' @examples
#' # After bootstrapping: print(summary(comparison))
#' @export
print.summary_multilpa_bootstrap_lrt <- function(x, digits = 4L, ...) {
  stopifnot("`x` must be a `summary_multilpa_bootstrap_lrt` object" =
              inherits(x, "summary_multilpa_bootstrap_lrt"),
            "`digits` must be a single number between 1 and 22" =
              is.numeric(digits) && length(digits) == 1L && is.finite(digits) &&
              digits >= 1 && digits <= 22)
  cat("Parametric bootstrap likelihood-ratio comparison\n\n")
  print(x$test, digits = digits, row.names = FALSE, ...)
  cat(sprintf("\n%d replicate(s) reached a parameter boundary; %d raised an error.\n",
              x$n_boundary, x$n_errors))
  cat("The reference distribution is simulated, not chi-square.\n")
  invisible(x)
}

#' Tidy a parametric bootstrap likelihood-ratio comparison
#' @param x An `multilpa_bootstrap_lrt` result, or its summary.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what `"test"` returns the single-row test result; `"replicates"`
#'   returns one row per simulated replicate.
#' @param ... Ignored.
#' @return A base `data.frame`. `"test"` has one row, with columns
#'   `null_profiles`, `null_group_classes`, `alternative_profiles`,
#'   `alternative_group_classes`, `statistic`, `p_value`, `monte_carlo_se`,
#'   `iter` and `n_valid`. `"replicates"` has one row per simulated dataset,
#'   with columns `replicate`, `statistic`, `valid`, `boundary`,
#'   `null_replications`, `alternative_replications`, `warnings` and `error`.
#' @examples
#' # After bootstrapping: as.data.frame(comparison, what = "replicates")
#' @export
as.data.frame.multilpa_bootstrap_lrt <- function(x, row.names = NULL,
                                                 optional = FALSE,
                                                 what = c("test", "replicates"),
                                                 ...) {
  stopifnot("`x` must be an `multilpa_bootstrap_lrt` result" =
              inherits(x, "multilpa_bootstrap_lrt"))
  what <- match.arg(what)
  result <- switch(what, test = .multilpa_bootstrap_test_frame(x),
                   replicates = x$replicates)
  row.names(result) <- row.names
  result
}

#' @rdname as.data.frame.multilpa_bootstrap_lrt
#' @export
as.data.frame.summary_multilpa_bootstrap_lrt <- function(x, row.names = NULL,
                                                         optional = FALSE,
                                                         what = c("test",
                                                                  "replicates"),
                                                         ...) {
  stopifnot("`x` must be a `summary_multilpa_bootstrap_lrt` object" =
              inherits(x, "summary_multilpa_bootstrap_lrt"))
  what <- match.arg(what)
  result <- switch(what, test = x$test, replicates = x$replicates)
  row.names(result) <- row.names
  result
}

#' Plot a simulated bootstrap null distribution
#'
#' Draws the replicate likelihood-ratio statistics as a histogram with the
#' observed statistic marked, so that the p-value can be read as a tail area of
#' the distribution that was actually simulated. Invalid replicates carry no
#' statistic and are counted in the subtitle rather than dropped silently.
#'
#' @param x An `multilpa_bootstrap_lrt` result.
#' @param main,subtitle Panel title and secondary line, or `NULL` for defaults.
#' @param palette A vector of at least two colours: the histogram fill and the
#'   observed-statistic marker. `NULL` uses the package palette.
#' @param style A list of visual constants, as built by `.multilpa_style()`.
#' @param ... Further named visual constants, merged into `style`.
#' @return The input, invisibly. Called for the side effect of drawing.
#' @section Conditions:
#'   `multilpa_nothing_to_plot` when no replicate produced a finite statistic.
#' @examples
#' # After bootstrapping: plot(comparison)
#' @export
plot.multilpa_bootstrap_lrt <- function(x, main = NULL, subtitle = NULL,
                                        palette = NULL,
                                        style = .multilpa_style(), ...) {
  stopifnot("`x` must be an `multilpa_bootstrap_lrt` result" =
              inherits(x, "multilpa_bootstrap_lrt"))
  values <- x$replicates$statistic
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    stop(errorCondition(
      "No replicate produced a finite statistic; there is nothing to plot.",
      class = "multilpa_nothing_to_plot", call = NULL))
  }
  style <- utils::modifyList(style, list(...))
  colours <- if (is.null(palette)) .multilpa_palette(2L) else
    rep(palette, length.out = 2L)
  previous <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA, mar = style$margins)
  breaks <- pretty(range(c(values, x$statistic, 0)), n = 12L)
  counts <- graphics::hist(values, breaks = breaks, plot = FALSE)
  .multilpa_panel(xlim = range(breaks), ylim = c(0, max(counts$counts) * 1.12),
    xlab = "Simulated likelihood-ratio statistic", ylab = "Replicates",
    main = if (is.null(main)) "Simulated null distribution" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d of %d replicates valid; observed %.3f, p = %s", x$n_valid, x$iter,
      x$statistic, format(x$p_value, digits = 3L)) else subtitle,
    style = style)
  graphics::rect(counts$breaks[-length(counts$breaks)], 0,
                 counts$breaks[-1L], counts$counts, col = colours[1L],
                 border = style$panel_fill, lwd = 1.2)
  graphics::abline(v = x$statistic, col = colours[2L], lwd = style$line_width,
                   lty = 2L)
  graphics::text(x$statistic, max(counts$counts) * 1.06, " observed",
                 adj = c(0, 0.5), col = colours[2L],
                 cex = style$label_text_size, font = 2L)
  invisible(x)
}

#' Spread the information criteria of one candidate into a single row
#' @param fit A fitted `multilpa` model, or `NULL` for a failed candidate.
#' @return A one-row `data.frame` of criteria, all `NA_real_` when `fit` is `NULL`.
#' @noRd
.multilpa_enumeration_indices <- function(fit) {
  names_wanted <- .multilpa_enumeration_criteria()
  if (is.null(fit)) {
    return(as.data.frame(stats::setNames(
      rep(list(NA_real_), length(names_wanted)), names_wanted)))
  }
  ## The same pivot `information_criteria(format = "wide")` performs, so the
  ## grid and a single fit cannot name the same quantity differently. The wide
  ## form leads with the likelihood and parameter count, which the grid already
  ## carries in its own columns.
  wide <- information_criteria(fit, format = "wide")
  as.data.frame(as.list(stats::setNames(
    lapply(names_wanted, function(name) wide[[name]]), names_wanted)))
}

#' The information criteria an enumeration grid carries
#'
#' Named once so that the failed-candidate row and the fitted row cannot drift
#' apart, and so that a renamed criterion is caught rather than silently
#' turned into a missing column.
#'
#' @return A character vector of grid column names, in grid order.
#' @noRd
.multilpa_enumeration_criteria <- function() {
  c("aic", "kic", "bic_groups", "bic_individual", "sabic_groups",
    "sabic_individual", "caic_groups", "caic_individual", "awe_groups",
    "awe_individual", "icl_groups", "icl_individual", "clc_groups",
    "clc_individual")
}
