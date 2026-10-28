# The checks of diagnostics(): one row per quantity, grouped in five domains.
# Classification follows each family's eligible posterior rows and weights.
# Values have no rule-of-thumb cutoffs; statuses mark facts about the fit.

#' A class label for a message: "class 2" for a number, the label otherwise
#' @noRd
.latents_class_label <- function(label) {
  label <- as.character(label)
  ifelse(grepl("^[0-9]+$", label), paste("class", label), label)
}

# The level of the bivariate-residual test, after Benjamini-Hochberg
# adjustment. The other quantities are reported without a verdict: their
# published cut-offs (AvePP .70, OCC 5, a 5% smallest class) are rules of
# thumb, not derived criteria, and the package does not impose them.
.latents_residual_alpha <- 0.05

.latents_check_references <- c(
  replication = "Masyn (2013)",
  avepp = "Nagin (2005)",
  occ = "Nagin (2005)",
  entropy = "Clark & Muthen (2009); Masyn (2013)",
  residuals = "Vermunt & Magidson (2005)")

#' One row of the checks table
#' @noRd
.latents_check_row <- function(domain, check, level = "model", value = NA_real_,
                               display = "", criterion = "", status = "info",
                               detail = "", reference = "") {
  data.frame(domain = domain, check = check, level = level, value = as.numeric(value),
             display = display, criterion = criterion, status = status, detail = detail,
             reference = reference, stringsAsFactors = FALSE)
}

#' The estimation checks: convergence, replication, starts, bounds
#' @noRd
.latents_estimation_checks <- function(x) {
  replication <- .multilpa_replication_table(x)
  converged <- isTRUE(x$converged)
  rows <- list(.latents_check_row(
    "Estimation", "best start converged", value = as.numeric(converged),
    display = if (converged) "yes" else "no",
    status = if (converged) "ok" else "fail",
    detail = if (converged) "" else "raise max_iter; the estimates are not at a maximum"))
  k <- replication$n_best_replicated
  n <- replication$n_starts
  if (!is.na(k) && !is.na(n)) {
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Estimation", "best log likelihood reached by", value = k,
      display = sprintf("%d of %d starts", k, n), status = "info",
      detail = if (n < 2L) "one start: the maximum is not replicated" else
        if (k == 1L) "one start only: the maximum is not replicated" else "",
      reference = .latents_check_references[["replication"]])
  }
  if (!is.na(replication$n_converged)) {
    unfinished <- n - replication$n_converged
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Estimation", "starts failed or unfinished", value = unfinished,
      display = sprintf("%d of %d%s", unfinished, n,
                        if (isTRUE(replication$n_screened > 0L)) {
                          sprintf(" (%d screened)", replication$n_screened)
                        } else ""),
      status = if (unfinished == 0L) "ok" else "check",
      detail = if (unfinished == 0L) "" else sprintf("%d failed, %d stopped at max_iter",
                                                     replication$n_failed,
                                                     unfinished - replication$n_failed))
  }
  boundary <- x$boundary %||% x$random_boundary
  if (!is.null(boundary)) {
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Estimation", "estimates at a bound", value = as.numeric(isTRUE(boundary)),
      display = if (isTRUE(boundary)) "yes" else "no",
      status = if (isTRUE(boundary)) "check" else "ok",
      detail = if (isTRUE(boundary)) {
        "a variance or probability sits at its lower bound; its standard error is not available"
      } else "")
  }
  if (isTRUE(x$small_classes)) {
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Estimation", "effective class size below one", value = 1, display = "yes",
       status = "fail",
      detail = "a class is estimated from less than one unit: fit fewer classes")
  }
  do.call(rbind, rows)
}

#' Eligible posterior rows and their sampling weights at a classification level
#'
#' The adapter keeps each family's actual classification units: observations,
#' groups, or persons and their clusters. Empty transition rows follow the
#' fitted drop_empty policy, and zero-weight units never enter any summary.
#' @param x A fitted model.
#' @param level The classification level, individuals or groups.
#' @return A list with posteriors and aligned weights; weights is NA only for
#'   an old weighted object that lacks a row-to-group map.
#' @noRd
.latents_classification_inputs <- function(x, level) {
  if (.multilpa_any_fit(x)) {
    posterior <- .multilpa_level_posteriors(x, level)
    return(list(posteriors = posterior, weights = if (is.null(posterior)) NULL else
      .multilpa_level_weights(x, level, nrow(posterior))))
  }
  weights <- NULL
  if (inherits(x, "latents_mixture_regression")) {
    grouped <- identical(x$spec$nesting, "group")
    posterior <- if (identical(level, "groups")) {
      if (grouped) x$expectation$group_tau else x$expectation$rho
    } else x$expectation$tau
    weights <- if (identical(level, "groups")) x$spec$sampling_weights else
      x$spec$row_weights
  } else if (inherits(x, "latents_growth_mixture")) {
    posterior <- if (identical(level, "groups")) x$expectation$structure$rho else
      x$expectation$posterior
    # Sampling weights are supported for persons, not the hierarchical model.
    weights <- if (identical(level, "individuals")) x$spec$sampling_weights
  } else {
    posterior <- .latents_level_posteriors(x, level)
    weights <- .latents_posterior_weights(x, level)
  }
  if (is.null(posterior)) return(list(posteriors = NULL, weights = NULL))
  if (identical(weights, NA)) return(list(posteriors = posterior, weights = NA))
  kept <- rowSums(is.na(posterior)) == 0L
  if (inherits(x, "multilpa_lta")) {
    kept <- kept & if (identical(level, "groups")) {
      .multilpa_informative_groups(x)
    } else .multilpa_informative_rows(x)
  }
  if (!is.null(weights)) {
    stopifnot("each posterior row has its sampling weight" =
                length(weights) == nrow(posterior) && all(is.finite(weights)))
    kept <- kept & weights > 0
  }
  list(posteriors = posterior[kept, , drop = FALSE],
       weights = if (is.null(weights)) NULL else unname(weights[kept]))
}

#' Per-class quantities, using the same eligible rows and weights as entropy
#' @noRd
.latents_level_classification <- function(x, level) {
  pieces <- .latents_classification_inputs(x, level)
  .latents_posterior_classification(pieces$posteriors, level, pieces$weights)
}

#' A stored posterior matrix at one classification level
#' @noRd
.latents_level_posteriors <- function(x, level) {
  posteriors <- if (identical(level, "groups")) x$group_posteriors else x$subject_posteriors
  if (!is.matrix(posteriors)) {
    stop(errorCondition(sprintf("This fit stores no %s posteriors to classify.", level),
                        class = "latents_bad_argument", call = NULL))
  }
  posteriors
}

#' Sampling weights aligned with a posterior-based family's rows
#'
#' Cross-level fits retain the row-to-group map in their sufficient statistics.
#' An older object without either map reports the individual level unavailable.
#' @noRd
.latents_posterior_weights <- function(x, level) {
  if (!.latents_is_weighted(x)) return(NULL)
  posteriors <- .latents_level_posteriors(x, level)
  weights <- unname(x$sampling_weights)
  if (identical(level, "groups")) return(as.numeric(weights))
  index <- x$group_index %||% x$sufficient_statistics$index
  if (!is.null(index) && length(index) == nrow(posteriors)) return(weights[index])
  if (length(weights) == nrow(posteriors)) return(as.numeric(weights))
  NA
}

#' Classification quantities from eligible posterior rows
#' @noRd
.latents_posterior_classification <- function(posteriors, level, weights = NULL) {
  stopifnot("posteriors must be a numeric matrix" = is.matrix(posteriors) &&
              is.numeric(posteriors))
  .latents_require_classification_weights(weights, level)
  kept <- rowSums(is.na(posteriors)) == 0L
  if (!is.null(weights)) kept <- kept & weights > 0
  posteriors <- posteriors[kept, , drop = FALSE]
  if (!is.null(weights)) weights <- weights[kept]
  table <- .multilpa_classification_summary(posteriors, level, weights)
  data.frame(level = table$level, class = as.character(table$class),
             estimated_n = table$estimated_n,
             estimated_proportion = table$estimated_proportion,
             n_assigned = table$n_modal, average_posterior = table$average_posterior,
             odds_correct_classification = table$odds_correct_classification,
             stringsAsFactors = FALSE)
}

#' Refuse a weighted classification level whose weights cannot be aligned
#' @noRd
.latents_require_classification_weights <- function(weights, level) {
  if (identical(weights, NA)) {
    stop(errorCondition(sprintf(paste(
      "The %s level of this weighted fit cannot be weighted: the fit stores its",
      "sampling weights per group and no map from rows to groups."), level),
      class = "latents_unweightable_level", call = NULL))
  }
  invisible(NULL)
}

#' The levels at which a fit classifies and their class counts
#' @noRd
.latents_classifying_levels <- function(x) {
  if (inherits(x, "multilpa_additive")) return(c(groups = ncol(x$group_posteriors)))
  if (inherits(x, "latents_mixture_regression")) {
    if (identical(x$spec$nesting, "group")) return(c(groups = x$spec$n_classes))
    return(c(individuals = x$spec$n_classes,
             if (identical(x$spec$nesting, "two-level")) c(groups = x$spec$n_group_classes)))
  }
  if (inherits(x, "latents_growth_mixture")) {
    return(c(individuals = x$spec$n_classes,
             if (!is.null(x$expectation$structure$rho))
               c(groups = ncol(x$expectation$structure$rho))))
  }
  if (inherits(x, "multilpa_cross_level")) {
    classes <- c(individuals = ncol(x$subject_posteriors), groups = ncol(x$group_posteriors))
    return(classes[classes > 1L | names(classes) == "individuals"])
  }
  # Classification includes the noise alternative, which n_profiles omits.
  # Read the same posterior columns that the entropy/classification tables use.
  individuals <- if (.multilpa_any_fit(x)) ncol(.multilpa_class_posteriors(x)) else
    x$n_profiles %||% x$n_classes %||% NA_integer_
  groups <- x$n_group_classes %||% 1L
  levels <- c(individuals = individuals)
  if (!isTRUE(x$single_level) && isTRUE(groups > 1L)) levels <- c(levels, groups = groups)
  levels
}

#' Classification tables for the newly supported model families
#'
#' Every diagnostics object has entropy as its primary table. A legacy
#' weighted object lacking a row map retains an entropy row with unavailable
#' values; its checks record the reason and unsupported classification rows
#' are omitted.
#' @noRd
.latents_family_diagnostic_tables <- function(x) {
  levels <- .latents_classifying_levels(x)
  pieces <- lapply(names(levels), function(level) .latents_classification_inputs(x, level))
  entropy <- do.call(rbind, Map(function(piece, level, k) {
    available <- !identical(piece$weights, NA)
    data.frame(level = level, n_classes = k,
      n_units = if (available) nrow(piece$posteriors) else NA_integer_,
      entropy_sum = if (available) .multilpa_entropy_sum(piece$posteriors, piece$weights) else NA_real_,
      relative_entropy = if (available) .multilpa_relative_entropy(piece$posteriors, piece$weights) else NA_real_)
  }, pieces, names(levels), unname(levels)))
  table_of <- function(fun) {
    tables <- Map(function(piece, level) {
      if (identical(piece$weights, NA)) return(NULL)
      fun(piece$posteriors, level, piece$weights)
    }, pieces, names(levels))
    combined <- do.call(rbind, tables)
    # Preserve the table contract even when an old object's only level has
    # no weight map. Its entropy/checks explain why no class rows are shown.
    if (is.null(combined)) fun(matrix(1, 1L, 1L), "individuals", NULL)[FALSE, , drop = FALSE] else
      combined
  }
  list(entropy = entropy,
       classification = table_of(.multilpa_classification_summary),
       average_posteriors = table_of(.multilpa_average_posterior_rows))
}

#' Class-size and classification checks of one level
#' @noRd
.latents_level_checks <- function(x, level, n_classes, entropy) {
  label <- if (identical(level, "groups")) "group" else "unit"
  if (isTRUE(n_classes == 1L)) {
    return(rbind(
      .latents_check_row("Classification", "relative entropy", level, status = "n/a",
                         detail = "one class: nothing to classify"),
      .latents_check_row("Classification", "lowest average posterior", level, status = "n/a",
                         detail = "one class: nothing to classify")))
  }
  table <- tryCatch(.latents_level_classification(x, level),
                    latents_unweightable_level = function(condition) conditionMessage(condition))
  if (is.character(table)) {
    return(.latents_check_row("Classification", "class sizes and classification", level,
                              status = "n/a", detail = table))
  }
  rows <- list()
  if (any(!is.na(table$estimated_proportion))) {
    smallest <- which.min(table$estimated_proportion)
    share <- table$estimated_proportion[[smallest]]
    count <- table$estimated_n[[smallest]]
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Class size", "smallest class", level, value = share,
      display = sprintf("%.1f%% (n = %s, %s)", 100 * share, format(round(count), big.mark = ","),
                        .latents_class_label(table$class[[smallest]])),
      status = "info")
    gap <- max(abs(table$n_assigned / sum(table$n_assigned) - table$estimated_proportion))
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Class size", "modal vs estimated proportions", level, value = gap,
      display = sprintf("largest gap %.3f", gap), status = "info")
  }
  share <- .latents_family_smallest_share(x)
  if (all(is.na(table$estimated_proportion)) && !is.na(share)) {
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Class size", "smallest class", level, value = share,
      display = sprintf("%.1f%%", 100 * share), status = "info")
  }
  empty <- table$class[table$n_assigned == 0L]
  rows[[length(rows) + 1L]] <- .latents_check_row(
    "Class size", sprintf("classes with no assigned %s", label), level, value = length(empty),
    display = sprintf("%d of %d", length(empty), nrow(table)),
    status = if (length(empty) == 0L) "ok" else "check",
    detail = if (length(empty) == 0L) "" else
      sprintf("%s: no %s has it as its most likely class",
              paste(.latents_class_label(empty), collapse = ", "), label))
  if (!is.na(entropy)) {
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Classification", "relative entropy", level, value = entropy,
      display = sprintf("%.3f", entropy), status = "info",
      reference = .latents_check_references[["entropy"]])
  }
  assigned <- table[table$n_assigned > 0L & !is.na(table$average_posterior), , drop = FALSE]
  if (nrow(assigned) > 0L) {
    lowest <- which.min(assigned$average_posterior)
    value <- assigned$average_posterior[[lowest]]
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Classification", "lowest average posterior", level, value = value,
      display = sprintf("%.3f (%s)", value, .latents_class_label(assigned$class[[lowest]])),
      status = "info", reference = .latents_check_references[["avepp"]])
  }
  odds <- assigned$odds_correct_classification
  if (any(is.finite(odds))) {
    lowest <- which.min(ifelse(is.finite(odds), odds, Inf))
    value <- odds[[lowest]]
    rows[[length(rows) + 1L]] <- .latents_check_row(
      "Classification", "lowest odds of correct classification", level, value = value,
      display = sprintf("%.1f (%s)", value, .latents_class_label(assigned$class[[lowest]])),
      status = "info", reference = .latents_check_references[["occ"]])
  }
  do.call(rbind, rows)
}

#' The smallest class share a family stores in its fit table, if any
#' @noRd
.latents_family_smallest_share <- function(x) {
  # Only these families store a smallest-class share in their fit table.
  if (!inherits(x, c("latents_mixture_regression", "multilpa_additive",
                     "latents_growth_mixture"))) {
    return(NA_real_)
  }
  fit <- get_results(x, "fit")
  share <- fit$smallest_share %||% fit$smallest_class %||% NA_real_
  if (isTRUE(share > 1)) NA_real_ else as.numeric(share)
}

#' Relative entropy at each actual classification level
#' @noRd
.latents_level_entropy <- function(x) {
  table <- if (.multilpa_any_fit(x)) .multilpa_entropy_table(x) else
    .latents_family_diagnostic_tables(x)$entropy
  stats::setNames(table$relative_entropy, table$level)
}

#' The local-independence check, from the bivariate residuals
#' @noRd
.latents_residual_checks <- function(residuals, reason, residual_note = NULL) {
  if (is.null(residuals) || nrow(residuals) == 0L) {
    return(.latents_check_row("Local independence", "bivariate residuals", status = "n/a",
                              detail = reason))
  }
  if (all(is.na(residuals$residual))) {
    return(.latents_check_row("Local independence", "bivariate residuals", status = "n/a",
                              detail = "no indicator pair has enough jointly observed units"))
  }
  adjusted <- residuals$p_adjusted
  if (all(is.na(adjusted))) {
    # Without a test (under sampling weights the residuals are descriptive),
    # a count of significant pairs would be a count of tests never run: the
    # residuals are reported, not judged.
    worst <- residuals[which.max(ifelse(is.na(residuals$residual), -Inf,
                                        abs(residuals$residual))), , drop = FALSE]
    return(.latents_check_row(
      "Local independence", "largest bivariate residual",
      value = abs(worst$residual),
      display = sprintf("%.3f (%d pairs)", worst$residual, nrow(residuals)),
      criterion = "not tested", status = "info",
      detail = sprintf("largest: %s x %s in %s; %s", worst$indicator_1, worst$indicator_2,
                       worst$profile, residual_note %||%
                         "no test is available for these residuals: read them as descriptive"),
      reference = .latents_check_references[["residuals"]]))
  }
  tested <- !is.na(adjusted)
  untested <- sum(!tested)
  flagged <- sum(adjusted[tested] < .latents_residual_alpha)
  worst <- residuals[which.min(ifelse(tested, adjusted, Inf)), , drop = FALSE]
  detail <- sprintf("largest: %s x %s in %s, residual %.3f, adjusted p %s", worst$indicator_1,
                    worst$indicator_2, worst$profile, worst$residual,
                    format.pval(worst$p_adjusted, digits = 2L, eps = 1e-4))
  .latents_check_row(
    "Local independence", "pairs dependent within a class", value = flagged,
    display = sprintf("%d of %d tested pairs%s", flagged, sum(tested),
                      if (untested > 0L) sprintf(" (%d not tested)", untested) else ""),
    criterion = sprintf("BH-adjusted p < %g", .latents_residual_alpha),
    # A pair without a test cannot pass: with untested pairs and none flagged,
    # the check is reported, not passed.
    status = if (flagged > 0L) "check" else if (untested > 0L) "info" else "ok",
    detail = if (flagged == 0L) {
      paste0(detail, if (untested > 0L) sprintf("; %d pairs have no test", untested) else "")
    } else paste0(detail, "; consider local_dependence = or a full covariance model"),
    reference = .latents_check_references[["residuals"]])
}

#' The information line: no verdict, a pointer to model comparison
#' @noRd
.latents_information_check <- function(x) {
  fit <- if (.multilpa_any_fit(x)) NULL else get_results(x, "fit")
  bic <- x$bic %||% fit$bic %||% NA_real_
  parameters <- x$n_parameters %||% fit$n_parameters %||% NA_integer_
  .latents_check_row(
    "Information", "log likelihood, parameters, BIC", value = x$log_likelihood,
    display = sprintf("%.2f, %s, %.2f", x$log_likelihood,
                      format(parameters), bic),
    criterion = "", status = "info",
    detail = "meaningful only against other models: compare_models()")
}

#' Every check of a fit
#' @param x A fitted model.
#' @param residuals The bivariate residual table, or `NULL`.
#' @param residual_reason Why there is none, for the `n/a` row.
#' @return The checks table.
#' @noRd
.latents_diagnostic_checks <- function(x, residuals, residual_reason) {
  levels <- .latents_classifying_levels(x)
  entropy <- .latents_level_entropy(x)
  level_rows <- lapply(names(levels), function(level) {
    .latents_level_checks(x, level, levels[[level]],
                          entropy[level] %||% NA_real_)
  })
  checks <- do.call(rbind, c(list(.latents_estimation_checks(x)), level_rows,
                             list(.latents_residual_checks(
                               residuals, residual_reason,
                               if (.latents_is_weighted(x)) {
                                 "not tested under sampling weights: read them as descriptive"
                               }),
                                  .latents_information_check(x))))
  domains <- c("Estimation", "Class size", "Classification", "Local independence",
               "Information")
  checks <- checks[order(match(checks$domain, domains), seq_len(nrow(checks))), , drop = FALSE]
  row.names(checks) <- NULL
  checks
}

#' Print the checks grouped by domain, attention first
#' @noRd
.latents_print_checks <- function(checks) {
  attention <- checks[checks$status %in% c("fail", "check"), , drop = FALSE]
  if (nrow(attention) == 0L) {
    cat("No check needs attention.\n")
  } else {
    cat(sprintf("%d %s attention: %s\n", nrow(attention),
                if (nrow(attention) == 1L) "check needs" else "checks need",
                paste(sprintf("%s%s", attention$check,
                              ifelse(attention$level %in% "groups", " (groups)", "")),
                      collapse = "; ")))
  }
  invisible(lapply(unique(checks$domain), function(domain) {
    rows <- checks[checks$domain == domain, , drop = FALSE]
    cat(sprintf("\n%s\n", domain))
    invisible(lapply(seq_len(nrow(rows)), function(index) {
      row <- rows[index, , drop = FALSE]
      level <- if (row$level %in% c("individuals", "groups")) row$level else ""
      cat(sprintf("  %-6s %-38s %-12s %-24s %s\n", row$status, row$check, level,
                  row$display, row$criterion))
      # The largest residual pair is worth reading even when none is flagged.
      if ((!row$status %in% c("ok", "info") || identical(row$domain, "Local independence")) &&
          nzchar(row$detail)) {
        cat(sprintf("         %s\n", row$detail))
      }
    }))
  }))
  invisible(NULL)
}
