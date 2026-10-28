# Class enumeration for mixture regressions.

#' Compare mixture regressions with different numbers of classes
#'
#' Fits [mixture_regression()] for every requested number of classes (crossed
#' with numbers of group classes, for the two-level model) and returns their
#' fit statistics in one table. Optionally adds a parametric bootstrap
#' likelihood-ratio test of `k - 1` against `k` classes, at the same number of
#' group classes (McLachlan 1987; Nylund, Asparouhov and Muthen 2007):
#' outcomes are simulated from the `k - 1` class fit, both models are refitted
#' to every replicate, and the p-value is `(1 + b) / (1 + V)`, where `b`
#' counts valid replicate statistics at least as large as the observed one and
#' `V` is the number of valid replicates.
#'
#' A replicate is invalid when a refit fails, does not converge, warns of
#' separation or a degenerate random-effect covariance, or gives a statistic
#' below zero beyond rounding. Under the null hypothesis the extra class of the
#' `k`-class refit is spurious, and a spurious class with membership
#' covariates is expected to drift towards a separating or slowly converging
#' solution: an invalid replicate is a property of the null distribution, not
#' a failure of the test. The p-value is therefore computed over the valid
#' replicates, their count and the invalid count are reported, and
#' `latents_failed_replicates` is warned naming the count; the p-value is `NA`
#' only when fewer than half of the requested replicates are valid. Excluding
#' invalid replicates conditions the null distribution on a valid refit, so a
#' large invalid share is worth reporting with the result.
#'
#' @inheritParams mixture_regression
#' @param n_classes Positive integer class counts to fit. Default `1:4`.
#' @param n_group_classes Positive integer group-class counts (two-level
#'   model; needs `id`), crossed with `n_classes`. Default 1.
#' @param bootstrap Number of bootstrap replicates for the likelihood-ratio
#'   test; `0` (the default) skips it. Not available with sampling weights.
#' @param bootstrap_starts Random starts for each bootstrap refit. Default 3.
#' @param ... Further arguments to [mixture_regression()], such as `family`,
#'   `id`, `class_level`, `common`, `membership` or `variance`. `n_classes`,
#'   `n_group_classes` and `vcov_type` are not accepted here.
#' @return An object of class `latents_regression_enumeration`. Its table
#'   (from [as.data.frame()] or `get_results(x, "fit")`) has one row per
#'   model, with columns `n_classes`, `n_group_classes`, `log_likelihood`,
#'   `n_parameters`, `aic`, `bic`, `bic_rows`, `sabic`, `icl`, `entropy`,
#'   `smallest_share`, `converged`, `n_best_replicated`, `note` (why a model
#'   could not be fitted), `best_bic` (`TRUE` for the minimum-BIC row), and,
#'   with `bootstrap > 0`, `blrt_statistic`, `blrt_p_value`,
#'   `blrt_replicates` (the valid replicates), `blrt_n_invalid` (the failed,
#'   unconverged, separated or reversed-likelihood replicates) and
#'   `blrt_flagged` (the same count, under its earlier name). The test
#'   columns are `NA` where the `k - 1` or `k` class model is missing or
#'   failed, and the p-value is `NA` when the observed fits fail validation or
#'   fewer than half of the requested replicates are valid.
#'   `get_results(x, "model", n_classes = , n_group_classes = )` returns one
#'   fitted model. `summary()` returns the same table as [as.data.frame()].
#'   `plot()` draws AIC, BIC, SABIC and ICL across the classes (see
#'   [plot_enumeration()] for its other views, and `plot_views(x)` for the
#'   views a regression grid offers).
#' @section Conditions:
#'   A model whose fit raises `latents_no_valid_start` or `latents_bad_data`
#'   at some count gets a row of `NA` with the condition message in `note`
#'   (and `print()` says how many); when every model fails with the same
#'   `latents_bad_data` error (a column that does not exist), that error is
#'   raised instead. Other errors stop the enumeration.
#'   `latents_bad_argument` for invalid `n_classes`, `n_group_classes`,
#'   `bootstrap` or `bootstrap_starts`, for an invalid `seed`, for an argument
#'   in `...` that [mixture_regression()] does not have, for `n_classes`,
#'   `n_group_classes` or `vcov_type` passed through `...`, for a bootstrap
#'   test whose outcome is not a plain column name, and from `get_results()`
#'   for a model that is not in the grid or could not be fitted.
#'   `latents_unsupported_weights` for a bootstrap test with sampling
#'   weights. `latents_failed_replicates` is warned, naming the count, when
#'   some replicates are invalid (the p-value then uses the valid ones), and
#'   when the observed fits fail validation or fewer than half of the
#'   replicates are valid (the p-value is then `NA`).
#' @references
#' McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test
#' statistic for the number of components in a normal mixture. *Applied
#' Statistics*, 36, 318--324.
#'
#' Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
#' number of classes in latent class analysis and growth mixture modeling: A
#' Monte Carlo simulation study. *Structural Equation Modeling*, 14,
#' 535--569.
#' @examples
#' classes <- enumerate_regressions(score ~ hours, data = study_hours,
#'                                  n_classes = 1:3, n_starts = 3, seed = 1)
#' classes
#' \donttest{
#' with_test <- enumerate_regressions(score ~ hours, data = study_hours,
#'                                    n_classes = 1:2, n_starts = 3, seed = 1,
#'                                    bootstrap = 19)
#' as.data.frame(with_test)
#' }
#' @export
enumerate_regressions <- function(formula, data, n_classes = 1:4,
                             n_group_classes = 1L, bootstrap = 0L,
                             bootstrap_starts = 3L, seed = NULL, ...) {
  data_expression <- substitute(data)
  .latents_require(
    "`n_classes` must be positive whole numbers" =
      is.numeric(n_classes) && length(n_classes) >= 1L &&
      all(vapply(n_classes, .mixture_is_count, logical(1))),
    "`n_group_classes` must be positive whole numbers" =
      is.numeric(n_group_classes) && length(n_group_classes) >= 1L &&
      all(vapply(n_group_classes, .mixture_is_count, logical(1))),
    "`bootstrap` must be a single non-negative whole number" =
      is.numeric(bootstrap) && length(bootstrap) == 1L && is.finite(bootstrap) &&
      bootstrap >= 0 && bootstrap == floor(bootstrap) && bootstrap <= .Machine$integer.max,
    "`bootstrap_starts` must be a single non-negative whole number" =
      is.numeric(bootstrap_starts) && length(bootstrap_starts) == 1L &&
      is.finite(bootstrap_starts) && bootstrap_starts >= 0 &&
      bootstrap_starts == floor(bootstrap_starts) && bootstrap_starts <= .Machine$integer.max
  )
  .multilpa_check_seed(seed)
  arguments <- list(...)
  # A misspelt argument would fail every cell with base R's "unused argument"
  # error, whose call held mixture_regression()'s body and the data by value
  # (API-1); name it back instead, as enumerate_classes() does.
  supplied <- names(arguments) %||% rep("", length(arguments))
  unknown <- !nzchar(supplied) | !supplied %in% names(formals(mixture_regression))
  if (any(unknown)) {
    .multilpa_reject_extra_arguments(arguments[unknown], "enumerate_regressions()",
      "Its `...` takes named arguments of mixture_regression(), such as `family`, `id` or `common`.")
  }
  if (bootstrap > 0 && !is.null(arguments$weights)) {
    .latents_refuse_weights(paste(
      "the bootstrap likelihood ratio test: it draws samples from the fitted",
      "model, not from the sampling design. Compare weighted fits on BIC"))
  }
  if (any(c("n_classes", "n_group_classes", "vcov_type") %in% names(arguments))) {
    stop(errorCondition(paste(
      "Pass the class counts as `n_classes` and `n_group_classes` of",
      "enumerate_regressions() itself; `vcov_type` is not used here."),
      class = "latents_bad_argument", call = NULL))
  }
  grid <- expand.grid(n_classes = sort(unique(as.integer(n_classes))),
                      n_group_classes = sort(unique(as.integer(n_group_classes))))
  .mixture_with_seed(seed, {
    fits <- lapply(seq_len(nrow(grid)), function(row) {
      fit <- .mixture_try_fit(formula, data, grid$n_classes[row],
                              grid$n_group_classes[row], arguments)
      # A failed cell keeps its condition; a fit, regression or growth mixture,
      # records a call that names the data rather than holding it by value.
      if (inherits(fit, "condition")) return(fit)
      fit$call <- .latents_candidate_call("mixture_regression", data_expression,
        c(list(formula = formula, n_classes = grid$n_classes[row],
               n_group_classes = grid$n_group_classes[row],
               vcov_type = "none"), arguments))
      fit
    })
    # A defect of the data (a column that does not exist) fails every cell for
    # one reason no class count changes: it is raised, not tabulated (API-2).
    .latents_refuse_shared_data_error(fits)
    table <- do.call(rbind, Map(function(fit, row) {
      .mixture_enumeration_row(fit, grid$n_classes[row],
                              grid$n_group_classes[row])
    }, fits, seq_len(nrow(grid))))
    # When every cell failed there is no lowest BIC, and no row is marked.
    lowest <- if (all(is.na(table$bic))) NA_real_ else min(table$bic, na.rm = TRUE)
    table$best_bic <- !is.na(table$bic) & table$bic == lowest
    if (bootstrap > 0) {
      tests <- do.call(rbind, lapply(seq_len(nrow(grid)), function(row) {
        .mixture_blrt(fits, grid, row, as.integer(bootstrap),
                     as.integer(bootstrap_starts), arguments)
      }))
      table <- cbind(table, tests)
    }
    structure(list(table = table, fits = fits, grid = grid),
              class = "latents_regression_enumeration")
  })
}

#' Fit one grid cell, keeping a failure as a note rather than an abort
#'
#' Only the conditions a fit legitimately raises for an unsupportable count
#' are caught; anything else is a defect and propagates. `do.call()` records
#' the arguments by value, so a propagating error is raised again with its
#' call bounded: it never carries the data (condition-calls.R).
#' @noRd
.mixture_try_fit <- function(formula, data, n_classes, n_group_classes,
                            arguments) {
  call_arguments <- c(list(formula = formula, data = data,
                           n_classes = n_classes,
                           n_group_classes = n_group_classes,
                           vcov_type = "none"), arguments)
  tryCatch(do.call(mixture_regression, call_arguments),
           latents_no_valid_start = .latents_bounded_condition,
           latents_bad_data = .latents_bounded_condition,
           error = function(condition) stop(.latents_bounded_condition(condition)))
}

#' One enumeration row
#' @noRd
.mixture_enumeration_row <- function(fit, n_classes, n_group_classes) {
  if (inherits(fit, "condition")) {
    return(data.frame(n_classes = n_classes, n_group_classes = n_group_classes,
                      log_likelihood = NA_real_, n_parameters = NA_integer_,
                      aic = NA_real_, bic = NA_real_, bic_rows = NA_real_,
                      sabic = NA_real_, icl = NA_real_, entropy = NA_real_,
                      smallest_share = NA_real_, converged = NA,
                      n_best_replicated = NA_integer_,
                      note = conditionMessage(fit)))
  }
  # A growth mixture fit has its own fit table with the same core columns.
  stats_row <- as.data.frame(get_results(fit, "fit"))
  data.frame(n_classes = n_classes, n_group_classes = n_group_classes,
             stats_row[c("log_likelihood", "n_parameters", "aic", "bic",
                         "bic_rows", "sabic", "icl", "entropy",
                         "smallest_share", "converged",
                         "n_best_replicated")],
             note = NA_character_)
}

#' Bootstrap likelihood-ratio test of the row against one class fewer
#' @noRd
.mixture_blrt <- function(fits, grid, row, replicates, starts, arguments) {
  empty <- data.frame(blrt_statistic = NA_real_, blrt_p_value = NA_real_,
                      blrt_replicates = NA_integer_, blrt_n_invalid = NA_integer_,
                      blrt_flagged = NA_integer_)
  k <- grid$n_classes[row]
  h <- grid$n_group_classes[row]
  smaller <- which(grid$n_classes == k - 1L & grid$n_group_classes == h)
  if (length(smaller) == 0L) return(empty)
  null_fit <- fits[[smaller]]
  alternative_fit <- fits[[row]]
  if (inherits(null_fit, "condition") ||
      inherits(alternative_fit, "condition")) return(empty)
  observed <- 2 * (alternative_fit$log_likelihood - null_fit$log_likelihood)
  response <- null_fit$spec$response_name
  base_data <- null_fit$spec$model_data
  if (!identical(make.names(response), response) ||
      !response %in% names(base_data)) {
    stop(errorCondition(paste(
      "The bootstrap test needs the outcome to be a plain column name, since",
      "it replaces that column with simulated values."),
      class = "latents_bad_argument", call = NULL))
  }
  refit_arguments <- arguments
  refit_arguments$n_starts <- starts
  refit_arguments$seed <- NULL
  refit_arguments$missing <- NULL
  tolerance <- 100 * (arguments$tol %||% 1e-8) *
    (1 + abs(null_fit$log_likelihood))
  if (!isTRUE(null_fit$converged) || !isTRUE(alternative_fit$converged) ||
      !is.finite(observed) || observed < -tolerance) {
    warning(warningCondition(
      "The observed fits failed validation, so blrt_p_value is NA.",
      class = "latents_failed_replicates", call = NULL))
    return(empty)
  }
  observed <- max(0, observed)
  # Replicate refits routinely meet unconverged or degenerate starts. Those
  # warnings mark the replicate invalid and are counted in `blrt_n_invalid`,
  # not repeated hundreds of times on the console.
  statistics <- vapply(seq_len(replicates), function(b) {
    simulated <- base_data
    simulated[[response]] <- if (inherits(null_fit, "latents_growth_mixture")) {
      .growth_draw(null_fit)
    } else .mixture_draw(null_fit)
    warned <- FALSE
    refit <- function(n) {
      tryCatch(withCallingHandlers(
        .mixture_try_fit(null_fit$spec$formula, simulated, n, h,
                        refit_arguments),
        latents_unconverged = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        },
        latents_degenerate_start = function(w) invokeRestart("muffleWarning"),
        latents_random_boundary = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        },
        latents_separation = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        }), error = identity)
    }
    null_refit <- refit(k - 1L)
    alternative_refit <- refit(k)
    invalid <- warned || inherits(null_refit, "condition") ||
      inherits(alternative_refit, "condition") ||
      !isTRUE(null_refit$converged) || !isTRUE(alternative_refit$converged)
    statistic <- if (invalid) NA_real_ else
      2 * (alternative_refit$log_likelihood - null_refit$log_likelihood)
    if (!is.finite(statistic) || statistic < -tolerance) return(NA_real_)
    max(0, statistic)
  }, numeric(1))
  # An invalid replicate is part of the null distribution's behaviour (a
  # spurious class separating or stalling), so the p-value is taken over the
  # valid replicates, and withheld only when fewer than half are valid.
  valid <- is.finite(statistics)
  n_valid <- sum(valid)
  n_invalid <- replicates - n_valid
  enough <- 2L * n_valid >= replicates && n_valid > 0L
  p_value <- if (enough) (1 + sum(statistics[valid] >= observed)) / (1 + n_valid) else
    NA_real_
  if (n_invalid > 0L) {
    warning(warningCondition(sprintf(paste(
      "%d of %d bootstrap replicates were invalid (a refit failed, did not",
      "converge or separated)%s"), n_invalid, replicates,
      if (enough) sprintf("; the p-value uses the %d valid replicates.", n_valid) else
        "; fewer than half are valid, so blrt_p_value is NA."),
      class = "latents_failed_replicates", call = NULL))
  }
  data.frame(blrt_statistic = observed, blrt_p_value = p_value,
             blrt_replicates = n_valid, blrt_n_invalid = n_invalid,
             blrt_flagged = n_invalid)
}

#' @rdname enumerate_regressions
#' @param x A `latents_regression_enumeration` object.
#' @param row.names,optional Not used; part of the generic.
#' @export
as.data.frame.latents_regression_enumeration <- function(x, row.names = NULL,
                                                     optional = FALSE, ...) {
  x$table
}

#' @rdname enumerate_regressions
#' @param what `"fit"` (the default) for the comparison table, `"model"` for
#'   one fitted model, `"tree"` for the class tree, or `"all"` for the named
#'   list of the two tables (`fit` and `tree`). Any other name raises
#'   `latents_bad_argument`.
#' @export
get_results.latents_regression_enumeration <- function(x, what = "fit",
                                                   n_classes = NULL,
                                                   n_group_classes = 1L, ...) {
  what <- .multilpa_match_table(what, c("fit", "model", "tree"))
  if (identical(what, "tree")) {
    # The tree counts its classes as `n_profiles`, as every grid's tree does.
    extra <- list(...)
    if (!missing(n_classes)) extra$n_profiles <- n_classes
    if (!missing(n_group_classes)) extra$n_group_classes <- n_group_classes
    return(.latents_tree_results(x, extra))
  }
  .multilpa_reject_extra_arguments(list(...), "get_results()", paste(
    "A regression grid's \"model\" takes `n_classes` and `n_group_classes`;",
    "its \"fit\" table takes no argument."))
  if (identical(what, "fit")) return(x$table)
  if (identical(what, "all")) {
    # The grid's tables; "model" is one fitted model, not a table.
    return(list(fit = x$table, tree = .latents_tree_results(x, list())))
  }
  .latents_require("`n_classes` must name one fitted class count" =
                     .mixture_is_count(n_classes %||% NA_real_))
  row <- which(x$grid$n_classes == n_classes &
                 x$grid$n_group_classes == n_group_classes)
  if (length(row) != 1L || inherits(x$fits[[row]], "condition")) {
    stop(errorCondition(sprintf(paste(
      "No fitted model with %s classes and %s group classes in this",
      "enumeration."), n_classes, n_group_classes),
      class = "latents_bad_argument", call = NULL))
  }
  fit <- x$fits[[row]]
  fit$inference <- if (inherits(fit, "latents_growth_mixture")) {
    .growth_resolve_inference(fit)
  } else .mixture_resolve_inference(fit, NULL)
  fit
}

#' @rdname enumerate_regressions
#' @param object A `latents_regression_enumeration` object.
#' @export
summary.latents_regression_enumeration <- function(object, ...) {
  stopifnot("`object` must be a `latents_regression_enumeration` grid" =
              inherits(object, "latents_regression_enumeration"))
  .latents_warn_unused_arguments(list(...), "summary()",
                                 "A regression grid's summary is its table.")
  .subset2(object, "table")
}

#' @rdname latents-plots
#' @export
plot.latents_regression_enumeration <- function(x, what = c("aic", "bic", "sabic",
                                                            "icl"),
                                                combine = TRUE, labels = TRUE,
                                                mark_minimum = TRUE, main = NULL,
                                                subtitle = NULL, cutoff = 0.8,
                                                using = "bic", ...) {
  stopifnot("`x` must be a `latents_regression_enumeration` grid" =
              inherits(x, "latents_regression_enumeration"))
  plot_enumeration(x, what = what, combine = combine, labels = labels,
                   mark_minimum = mark_minimum, main = main,
                   subtitle = subtitle, cutoff = cutoff, using = using, ...)
}

#' @rdname enumerate_regressions
#' @param digits Significant digits.
#' @export
print.latents_regression_enumeration <- function(x, digits = 4L, ...) {
  cat("Mixture regression class enumeration\n\n")
  shown <- x$table[intersect(c("n_classes", "n_group_classes",
                               "log_likelihood", "n_parameters", "bic",
                               "icl", "entropy", "smallest_share",
                               "converged", "blrt_p_value", "blrt_n_invalid",
                               "best_bic"),
                             names(x$table))]
  print(shown, digits = digits, row.names = FALSE)
  failed <- !is.na(x$table$note)
  if (any(failed)) {
    cat(sprintf(paste("\n%d model(s) could not be fitted; the `note` column of",
                      "as.data.frame(x) gives the reason.\n"), sum(failed)))
  }
  cat("\nOne model: get_results(x, \"model\", n_classes = ).\n")
  invisible(x)
}
