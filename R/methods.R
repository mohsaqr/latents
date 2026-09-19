#' Print a fitted multilevel latent profile model
#' @param x An `multilpa` model.
#' @param ... Reserved for compatibility with `print()`.
#' @return The input model, invisibly.
#' @examples
#' # Printing a fitted model displays its convergence and information criteria.
#' @export
print.multilpa <- function(x, ...) {
  stopifnot(inherits(x, "multilpa"))
  cat(sprintf("Two-level latent profile analysis: %d profiles, %d group classes\n",
              x$n_profiles, x$n_group_classes))
  covariance_model <- x$covariance_model %||% "diagonal"
  cat(sprintf("%d individuals in %d groups; %s %s residual covariance\n",
              x$n_observations, x$n_groups, x$variance_model, covariance_model))
  cat(sprintf("Log likelihood: %.6f | AIC: %.3f | BIC (groups): %.3f\n",
              x$log_likelihood, x$aic, x$bic))
  cat(sprintf("Converged: %s | iterations: %d | best start: %d/%d\n",
              x$converged, x$iterations, x$best_start, nrow(x$starts)))
  if (!is.null(x$n_informative) && x$n_informative < x$n_observations) {
    cat(sprintf("%d row(s) carry no observed indicator; the individual-level BIC uses %d.\n",
                x$n_observations - x$n_informative, x$n_informative))
  }
  if (x$boundary) cat("A variance is at its specified lower bound.\n")
  if (x$small_classes) cat("An effective class membership is below one.\n")
  invisible(x)
}

#' Summarize a fitted multilevel latent profile model
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `summary()`.
#' @return A `summary_multilpa` list of estimates, effective class counts,
#'   information criteria, and restart diagnostics.
#' @examples
#' # After fitting: summary(fit)
#' @export
summary.multilpa <- function(object, ...) {
  stopifnot(inherits(object, "multilpa"))
  fields <- c("call", "n_observations", "n_informative", "n_groups",
              "n_profiles", "n_group_classes",
              "variance_model", "covariance_model", "missing", "covariances",
              "means", "variances", "standard_deviations",
              "profile_probabilities", "group_probabilities", "log_likelihood",
              "n_parameters", "aic", "bic", "bic_individual", "converged",
              "iterations", "boundary", "min_variance", "small_classes", "starts",
              "best_start", "n_best_replicated", "replication_tolerance")
  result <- object[fields]
  result$effective_profile_counts <- colSums(object$subject_posteriors)
  result$effective_group_counts <- colSums(object$group_posteriors)
  class(result) <- "summary_multilpa"
  result
}

#' Print a multilevel LPA summary
#' @param x A `summary_multilpa` object.
#' @param digits Number of printed significant digits.
#' @param ... Additional arguments passed to matrix printing.
#' @return The summary, invisibly.
#' @examples
#' # After fitting: print(summary(fit), digits = 4)
#' @export
print.summary_multilpa <- function(x, digits = 4L, ...) {
  stopifnot(inherits(x, "summary_multilpa"), is.numeric(digits),
            length(digits) == 1L, is.finite(digits), digits >= 1, digits <= 22)
  cat(sprintf("Multilevel LPA: %d profiles and %d group classes\n",
              x$n_profiles, x$n_group_classes))
  cat(sprintf("Individuals: %d; groups: %d; parameters: %d; converged: %s\n",
              x$n_observations, x$n_groups, x$n_parameters, x$converged))
  cat("\nProfile means:\n")
  print(x$means, digits = digits, ...)
  cat("\nProfile standard deviations:\n")
  print(x$standard_deviations, digits = digits, ...)
  if (identical(x$covariance_model, "full")) {
    cat("\nProfile residual covariance matrices:\n")
    print(x$covariances, digits = digits, ...)
  }
  cat("\nProfile probabilities within each group class:\n")
  print(x$profile_probabilities, digits = digits, ...)
  cat("\nGroup-class probabilities:\n")
  print(x$group_probabilities, digits = digits, ...)
  cat("\nEffective individual memberships:\n")
  print(x$effective_profile_counts, digits = digits, ...)
  cat("\nEffective group memberships:\n")
  print(x$effective_group_counts, digits = digits, ...)
  cat(sprintf("\nLog likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              x$log_likelihood, x$aic, x$bic, x$bic_individual))
  cat(sprintf("Best likelihood replicated in %d/%d starts (absolute tolerance %.3g).\n",
              x$n_best_replicated, nrow(x$starts), x$replication_tolerance))
  if (!x$converged) cat("WARNING: best start did not converge.\n")
  if (x$boundary) cat("WARNING: at least one variance is at min_variance.\n")
  if (x$small_classes) cat("WARNING: an effective class membership is below one.\n")
  cat("\nStart diagnostics:\n")
  print(x$starts, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Extract the multilevel model log likelihood
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `logLik()`.
#' @return A `logLik` object with parameter count `df` and the number of observed
#'   groups as `nobs`. Thus `stats::BIC()` uses group-count BIC. For the
#'   individual-count alternative, and every other criterion, call
#'   [information_criteria()], which reports both conventions side by side.
#' @examples
#' # After fitting: logLik(fit); AIC(fit); BIC(fit)
#' @export
#' @importFrom stats logLik
logLik.multilpa <- function(object, ...) {
  stopifnot(inherits(object, "multilpa"))
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' Extract the number of independent groups
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `nobs()`.
#' @return Number of observed groups, which are the independent units of the
#'   two-level likelihood. For the individual count alongside every other
#'   sample-size-dependent quantity, call [information_criteria()].
#' @examples
#' # After fitting: nobs(fit)
#' @export
#' @importFrom stats nobs
nobs.multilpa <- function(object, ...) {
  stopifnot(inherits(object, "multilpa"))
  object$n_groups
}
