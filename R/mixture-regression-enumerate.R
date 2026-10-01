# Class enumeration for mixture regressions.

#' Compare mixture regressions with different numbers of classes
#'
#' Fits [mixture_regression()] for every requested number of classes (and of group
#'                          classes, for the two-level model) and returns their fit statistics in one
#' table. Optionally adds a parametric bootstrap likelihood-ratio test of
#' `k - 1` against `k` classes (McLachlan 1987; Nylund, Asparouhov and Muthen
#' 2007): outcomes are simulated from the `k - 1` class fit, both models are
#' refitted to every replicate, and the p-value is the share of replicate
#' statistics at least as large as the observed one, `(1 + b) / (1 + B)`.
#'
#' @inheritParams mixture_regression
#' @param n_classes Integer vector of class counts to fit.
#' @param n_group_classes Integer vector of group-class counts (two-level
#'   model; needs `id`). Crossed with `n_classes`.
#' @param bootstrap Number of bootstrap replicates for the likelihood-ratio
#'   test; `0` skips it.
#' @param bootstrap_starts Random starts for each bootstrap refit.
#' @param ... Further arguments to [mixture_regression()], such as `family`, `id`,
#'   `class_level`, `common`, `membership` or `variance`.
#' @return An object of class `latents_regression_enumeration`. Its table (from
#'   [as.data.frame()] or `get_results(x, "fit")`) has one row per model:
#'   `n_classes`, `n_group_classes`, `log_likelihood`, `n_parameters`,
#'   `aic`, `bic`, `bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share`,
#'   `converged`, `n_best_replicated`, `best_bic` (the minimum-BIC row), and,
#'   with `bootstrap > 0`, `blrt_statistic`, `blrt_p_value` and
#'   `blrt_replicates` (the successful replicates) and `blrt_flagged` (the
#'   replicates in which a refit did not converge or separated). `get_results(x, "model",
#'   n_classes = , n_group_classes = )` returns one fitted model.
#' @section Conditions:
#'   `latents_bad_argument` for an empty or invalid grid. A model that cannot
#'   be fitted at some count (for example every start degenerates) gets a row
#'   of `NA` with its condition message in `note`, rather than stopping the
#'   comparison.
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
  stopifnot(
    "`n_classes` must be positive whole numbers" =
      is.numeric(n_classes) && length(n_classes) >= 1L &&
      all(vapply(n_classes, .mixture_is_count, logical(1))),
    "`n_group_classes` must be positive whole numbers" =
      is.numeric(n_group_classes) && length(n_group_classes) >= 1L &&
      all(vapply(n_group_classes, .mixture_is_count, logical(1))),
    "`bootstrap` must be a single non-negative whole number" =
      length(bootstrap) == 1L && is.finite(bootstrap) && bootstrap >= 0,
    "`bootstrap_starts` must be a single non-negative whole number" =
      length(bootstrap_starts) == 1L && is.finite(bootstrap_starts) &&
      bootstrap_starts >= 0
  )
  arguments <- list(...)
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
      .mixture_try_fit(formula, data, grid$n_classes[row],
                      grid$n_group_classes[row], arguments)
    })
    table <- do.call(rbind, Map(function(fit, row) {
      .mixture_enumeration_row(fit, grid$n_classes[row],
                              grid$n_group_classes[row])
    }, fits, seq_len(nrow(grid))))
    table$best_bic <- !is.na(table$bic) &
      table$bic == min(table$bic, na.rm = TRUE)
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
#' are caught; anything else is a defect and propagates.
#' @noRd
.mixture_try_fit <- function(formula, data, n_classes, n_group_classes,
                            arguments) {
  call_arguments <- c(list(formula = formula, data = data,
                           n_classes = n_classes,
                           n_group_classes = n_group_classes,
                           vcov_type = "none"), arguments)
  tryCatch(do.call(mixture_regression, call_arguments),
           latents_no_valid_start = function(condition) condition,
           latents_bad_data = function(condition) condition)
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
  stats_row <- .mixture_fit_table(fit)
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
                      blrt_replicates = NA_integer_, blrt_flagged = NA_integer_)
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
  # Replicate refits routinely meet unconverged or degenerate starts. Those
  # warnings are counted per replicate and reported in `blrt_flagged`, not
  # repeated hundreds of times on the console.
  flagged <- 0L
  statistics <- vapply(seq_len(replicates), function(b) {
    simulated <- base_data
    simulated[[response]] <- .mixture_draw(null_fit)
    warned <- FALSE
    refit <- function(n) {
      withCallingHandlers(
        .mixture_try_fit(null_fit$spec$formula, simulated, n, h,
                        refit_arguments),
        latents_unconverged = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        },
        latents_degenerate_start = function(w) invokeRestart("muffleWarning"),
        latents_separation = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        })
    }
    null_refit <- refit(k - 1L)
    alternative_refit <- refit(k)
    if (warned) flagged <<- flagged + 1L
    if (inherits(null_refit, "condition") ||
        inherits(alternative_refit, "condition")) return(NA_real_)
    2 * (alternative_refit$log_likelihood - null_refit$log_likelihood)
  }, numeric(1))
  valid <- statistics[is.finite(statistics)]
  data.frame(blrt_statistic = observed,
             blrt_p_value = (1 + sum(valid >= observed)) / (1 + length(valid)),
             blrt_replicates = length(valid), blrt_flagged = flagged)
}

#' @rdname enumerate_regressions
#' @param x A `latents_regression_enumeration` object.
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_regression_enumeration <- function(x, row.names = NULL,
                                                     optional = FALSE, ...) {
  x$table
}

#' @rdname enumerate_regressions
#' @param what `"fit"` for the comparison table, `"model"` for one fitted
#'   model.
#' @export
get_results.latents_regression_enumeration <- function(x, what = c("fit", "model"),
                                                   n_classes = NULL,
                                                   n_group_classes = 1L, ...) {
  what <- match.arg(what)
  if (identical(what, "fit")) return(x$table)
  stopifnot("`n_classes` must name one fitted class count" =
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
  fit$inference <- .mixture_inference(fit, "observed")
  fit
}

#' @rdname enumerate_regressions
#' @param digits Significant digits.
#' @export
print.latents_regression_enumeration <- function(x, digits = 4L, ...) {
  cat("Mixture regression class enumeration\n\n")
  shown <- x$table[intersect(c("n_classes", "n_group_classes",
                               "log_likelihood", "n_parameters", "bic",
                               "icl", "entropy", "smallest_share",
                               "converged", "blrt_p_value", "best_bic"),
                             names(x$table))]
  print(shown, digits = digits, row.names = FALSE)
  cat("\nOne model: get_results(x, \"model\", n_classes = ).\n")
  invisible(x)
}
