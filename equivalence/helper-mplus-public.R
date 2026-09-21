#' Read the narrow two-class, four-indicator Mplus example output format
#' @param lines Character lines from official examples 7.9 or 7.10.
#' @return Published scalar fit statistics and class-specific estimates.
parse_public_mplus <- function(lines) {
  stopifnot(is.character(lines), !anyNA(lines))
  scalar_patterns <- c(n_observations = "^Number of observations",
                       n_parameters = "^Number of Free Parameters",
                       log_likelihood = "^\\s*H0 Value", aic = "^\\s*Akaike ",
                       bic_individual = "^\\s*Bayesian ", entropy = "^\\s*Entropy ")
  scalars <- vapply(scalar_patterns, function(pattern) {
    stopifnot(is.character(pattern))
    matched <- grep(pattern, lines, value = TRUE)
    stopifnot(length(matched) == 1L)
    as.numeric(sub(".*\\s(-?[0-9]+(?:\\.[0-9]+)?)\\s*$", "\\1", matched, perl = TRUE))
  }, numeric(1))
  class_indices <- grep("^Latent Class [12]$", trimws(lines))
  stopifnot(length(class_indices) == 2L)
  estimates <- lapply(class_indices, function(index) {
    stopifnot(is.integer(index), length(index) == 1L)
    rows <- grep("^\\s+Y[1-4]\\s", lines[seq.int(index + 1L, index + 13L)], value = TRUE)
    stopifnot(length(rows) == 8L)
    values <- as.numeric(sub("^\\s+Y[1-4]\\s+(-?[0-9.]+).*", "\\1", rows))
    list(means = values[seq_len(4L)], variances = values[4L + seq_len(4L)])
  })
  table_anchors <- c(model = "^BASED ON THE ESTIMATED MODEL$",
                     modal = "^BASED ON THEIR MOST LIKELY LATENT CLASS MEMBERSHIP$")
  counts <- lapply(table_anchors, function(pattern) {
    stopifnot(is.character(pattern))
    index <- grep(pattern, lines)
    stopifnot(length(index) == 1L)
    rows <- grep("^\\s+[12]\\s+[0-9.]+\\s+[0-9.]+\\s*$",
                 lines[seq.int(index + 1L, index + 11L)], value = TRUE)
    stopifnot(length(rows) == 2L)
    parsed <- utils::read.table(text = paste(rows, collapse = "\n"),
                                col.names = c("class", "count", "proportion"))
    stopifnot(!anyNA(parsed), identical(parsed$class, 1:2))
    parsed
  })
  means <- do.call(rbind, lapply(estimates, `[[`, "means"))
  variances <- do.call(rbind, lapply(estimates, `[[`, "variances"))
  ordering <- order(means[, 1L])
  c(as.list(scalars), list(means = means[ordering, , drop = FALSE],
                          variances = variances[ordering, , drop = FALSE],
                          proportions = counts$model$proportion[ordering],
                          effective_counts = counts$model$count[ordering],
                          modal_counts = counts$modal$count[ordering]))
}

#' Compare a fitted single-level limit against published Mplus estimates
#' @param fit An multilpa fit with two profiles and one group class.
#' @param expected Parsed published estimates ordered by first indicator mean.
#' @return A table of maximum absolute discrepancies and absolute tolerances.
compare_public_mplus <- function(fit, expected) {
  stopifnot(inherits(fit, "multilpa"), is.list(expected),
            fit$n_profiles == 2L, fit$n_group_classes == 1L)
  ordering <- order(fit$means[, "y1"])
  probabilities <- fit$subject_posteriors
  observed <- list(n_observations = fit$n_observations,
                   n_parameters = fit$n_parameters,
                   log_likelihood = fit$log_likelihood, aic = fit$aic,
                   bic_individual = fit$bic_individual,
                   entropy = 1 + sum(probabilities * log(pmax(probabilities, .Machine$double.xmin))) /
                     (fit$n_observations * log(fit$n_profiles)),
                   means = unname(fit$means[ordering, , drop = FALSE]),
                   variances = unname(fit$variances[ordering, , drop = FALSE]),
                   proportions = as.numeric(fit$profile_probabilities[1L, ordering]),
                   effective_counts = as.numeric(fit$effective_profile_counts[ordering]),
                   modal_counts = tabulate(fit$subject_profiles, nbins = 2L)[ordering])
  stopifnot(setequal(names(observed), names(expected)))
  errors <- vapply(names(observed), function(field) {
    stopifnot(is.character(field), length(observed[[field]]) == length(expected[[field]]))
    max(abs(observed[[field]] - expected[[field]]))
  }, numeric(1))
  # Counts combine five-decimal printing with accumulated optimizer differences:
  # 0.001 subjects is an explicit optimization allowance, not reported precision.
  tolerance <- c(n_observations = 0, n_parameters = 0, log_likelihood = 0.0005,
                 aic = 0.0005, bic_individual = 0.0005, entropy = 0.0005,
                 means = 0.0005, variances = 0.0005, proportions = 0.000005,
                 effective_counts = 0.001, modal_counts = 0)
  data.frame(metric = names(errors), max_abs_difference = unname(errors),
             tolerance = unname(tolerance[names(errors)]),
             passed = unname(errors <= tolerance[names(errors)]))
}
