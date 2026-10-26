# Every output, not one or two measures.
#
# Until 2026-10-26 each suite chose by hand what to compare, and the choices
# were narrow: the depmixS4 suite compared the first indicator's means and the
# transitions, mclust the means and mixing proportions, poLCA ran without
# standard errors. Nothing compared posteriors, classifications, the
# classification error matrix or the BCH weights, which is how the empty-row
# defect in the error matrix (fixed 2026-10-26) went unseen.
#
# This file makes the comparison complete and makes its completeness visible.
#
#   fit_outputs()       every value a fit reports, one row per table column,
#                       read from the package's own catalogue
#                       (`get_results(fit, "all")`), so a column added to the
#                       package appears here without anyone listing it.
#   reference_tables()  the same tables built from a reference program's
#                       estimates and posteriors, with every derived table
#                       (counts, average posteriors, odds of correct
#                       classification, error matrix, BCH weights, entropy,
#                       information criteria, residuals, Gaussian diagnostics)
#                       written here from its definition.
#   compare_tables()    one comparison row per output column, over every cell
#                       of it, with the worst cell reported.
#   compare_fit()       the three above, plus the ledger: for every output
#                       column, how many values were produced, how many were
#                       compared, and, for anything not compared, the reason.
#
# The ledger is the point. A column that is neither compared nor exempted with
# a stated reason is `not covered`, and `run.R` fails on it.

# ---------------------------------------------------------------------------
# What a fit reports
# ---------------------------------------------------------------------------

#' Columns that name a row of a table rather than carry a value
#'
#' `profile` and `group_class` are keys everywhere except `assignments`, where
#' they are the modal classes, so that table is keyed by position instead (see
#' [.output_keys()]).
.output_identifiers <- c(
  "row", "group", "unit", ".row", "profile", "group_class", "class",
  "assigned_class", "true_class", "indicator", "indicator_2", "indicator_1",
  "category", "level", "context", "kind", "n_classes", "outcome", "term",
  "parameter", "time", "occasion", "from", "to", "state", "criterion",
  "convention", "definition")

#' Outputs that are not estimates, and why they are never compared
#'
#' Keys are a table name or `table$column`. Anything else a fit reports must be
#' compared or exempted by the suite with its own reason.
.output_exemptions <- c(
  data = "the input data, returned as given",
  starts = "this program's random-start trace; another program's starts are not comparable",
  stages = "this program's search stages; another program's search is not comparable",
  `model$converged` = "this optimiser's own stopping flag",
  `model$iterations` = "this optimiser's iteration count",
  `model$best_start` = "which of this program's random starts won",
  `model$n_best_replicated` = "how many of this program's random starts reached the best value",
  `model$boundary` = "this package's warning flag, with its own threshold (unit-tested)",
  `model$small_classes` = "this package's warning flag, with its own threshold (unit-tested)")

#' Columns named like an identifier that are values in a particular table
#'
#' `time` is the occasion key of a transition fit's `sequences` table and the
#' random-slope prediction of a growth fit's `random_effects` table.
.output_value_columns <- list(random_effects = "time")

#' The key columns of one output table
#' @param table A table from `get_results()`.
#' @param name Its catalogue name.
#' @return A character vector of column names.
.output_keys <- function(table, name) {
  if (identical(name, "assignments")) return(".row")
  labels <- names(table)[!vapply(table, function(column) {
    is.numeric(column) || is.logical(column)
  }, logical(1))]
  setdiff(union(intersect(names(table), .output_identifiers), labels),
          .output_value_columns[[name]] %||% character())
}

#' Every table a fit reports, as compared here
#'
#' `assignments` gains a `.row` key, and loses the columns that echo the input
#' (the indicators and the identifier), which are data rather than output.
#'
#' @param fit A fitted model of the package.
#' @return A named list of `data.frame`s.
fit_tables <- function(fit) {
  tables <- get_results(fit, "all")
  if (!is.null(tables$assignments)) {
    # The indicators, the identifier and the sampling-weight column are input.
    weight_column <- if (is.character(fit$weights)) fit$weights else character()
    echoed <- intersect(names(tables$assignments), c(fit$vars, fit$id, weight_column))
    assignments <- tables$assignments[setdiff(names(tables$assignments), echoed)]
    tables$assignments <- data.frame(.row = seq_len(nrow(assignments)), assignments,
                                     check.names = FALSE)
  }
  tables
}

#' Every value a fit reports, one row per table column
#'
#' @param fit A fitted model, or the list [fit_tables()] returned for one.
#' @return A `data.frame` with columns `output` (`table$column`), `table`,
#'   `column` and `n_values` (cells, missing ones included: a missing value is
#'   an output too).
fit_outputs <- function(fit) {
  tables <- if (inherits(fit, "list") && !is.null(names(fit)) &&
                all(vapply(fit, is.data.frame, logical(1)))) fit else fit_tables(fit)
  rows <- lapply(names(tables), function(name) {
    table <- tables[[name]]
    values <- setdiff(names(table)[vapply(table, function(column) {
      is.numeric(column) || is.logical(column)
    }, logical(1))], .output_keys(table, name))
    if (length(values) == 0L || nrow(table) == 0L) return(NULL)
    data.frame(output = paste0(name, "$", values), table = name, column = values,
               n_values = nrow(table), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(output = character(), table = character(), column = character(),
                      n_values = integer(), stringsAsFactors = FALSE)
  }
  row.names(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Comparing a whole column
# ---------------------------------------------------------------------------

#' Compare every cell of one output column
#'
#' Missing values must be missing in both, and infinite values infinite with the
#' same sign in both; such a pair agrees exactly. Every other pair is compared
#' on the relative scale, `|r - o| / max(|r|, 1)`, which is the absolute
#' difference for a probability and the relative one for a large quantity.
#' Logical columns must match exactly.
#'
#' @param quantity Label for the row.
#' @param reference,obtained Vectors of equal length.
#' @param tolerance Positive tolerance.
#' @param output The `table$column` key, or `NA`.
#' @return A one-row `data.frame` shaped like [compare_values()]: `reference`,
#'   `obtained` and `difference` are those of the worst cell, and `n_values`
#'   is how many cells were compared.
compare_block <- function(quantity, reference, obtained, tolerance, output = NA_character_) {
  stopifnot(
    "`reference` and `obtained` must have the same length" =
      length(reference) == length(obtained),
    "`tolerance` must be one positive number" =
      is.numeric(tolerance) && length(tolerance) == 1L && is.finite(tolerance) &&
        tolerance > 0)
  exact <- is.logical(reference) || is.logical(obtained)
  reference <- as.numeric(reference)
  obtained <- as.numeric(obtained)
  n <- length(reference)
  both_missing <- is.na(reference) & is.na(obtained)
  both_infinite <- is.infinite(reference) & is.infinite(obtained) &
    sign(reference) == sign(obtained)
  matched <- both_missing | both_infinite
  difference <- abs(reference - obtained) / pmax(abs(reference), 1)
  difference[matched] <- 0
  difference[is.na(difference)] <- Inf
  if (exact) tolerance <- 0.5
  finite_pairs <- is.finite(reference) & is.finite(obtained)
  worst <- if (any(!is.finite(difference))) which(!is.finite(difference))[1L] else
    if (any(finite_pairs)) which(finite_pairs)[which.max(difference[finite_pairs])] else NA
  note <- if (n > 0L && all(matched)) ", all missing or infinite in both" else ""
  if (is.na(worst)) {
    picked <- c(0, 0, 0)
  } else {
    picked <- c(reference[[worst]], obtained[[worst]], difference[[worst]])
  }
  agrees <- n > 0L && all(difference <= tolerance)
  # A matched non-finite column carries no finite cell to show; zeros stand in
  # so the row is not voided by the finiteness check in `.evaluate_suite()`.
  if (agrees && !all(is.finite(picked))) picked <- c(0, 0, 0)
  data.frame(quantity = sprintf("%s (%d values%s)", quantity, n, note),
             reference = picked[[1L]], obtained = picked[[2L]],
             difference = picked[[3L]], tolerance = tolerance, agrees = agrees,
             output = output, n_values = n, stringsAsFactors = FALSE)
}

#' Look up the tolerance for one output column
#' @param tolerance Named numeric: `table$column`, then `table`, then `default`.
#' @param table,column The output.
#' @return A single tolerance.
.output_tolerance <- function(tolerance, table, column) {
  key <- paste0(table, "$", column)
  for_key <- tolerance[c(key, table, "default")]
  for_key <- for_key[!is.na(for_key)]
  stopifnot("a tolerance must be given for every output, or a default" =
              length(for_key) > 0L)
  unname(for_key[[1L]])
}

#' The tolerances `compare_fit()` uses unless told otherwise
#'
#' Two converged optimisers agree on estimates to about the square root of
#' their stopping tolerance; standard errors come from numerical second
#' derivatives in at least one of the two programs, and agree less closely.
.output_default_tolerances <- c(
  default = 1e-6,
  `profiles$mean_standard_error` = 1e-4,
  `profiles$variance_standard_error` = 1e-4,
  `responses$probability_standard_error` = 1e-4,
  `responses$threshold` = 1e-4,
  # OCC divides by 1 - AvePP: a gap d in the average posterior a becomes a
  # relative gap of d / (a (1 - a)) in the odds, about 2000 d at a = 0.9995.
  # Posteriors agreeing to 1e-8 therefore support OCC to about 1e-4.
  `classification$odds_correct_classification` = 1e-4)

#' Compare every expected table with what a fit reports
#'
#' Rows are joined on the key columns the two tables share. A reference row the
#' fit does not report is a disagreement; a fit row the reference does not
#' cover is left uncompared and shows up in the ledger as partial coverage.
#'
#' @param label Label for the rows.
#' @param obtained Named list of tables, from [fit_tables()].
#' @param expected Named list of tables in the same shapes, each holding its
#'   key columns and the value columns to compare.
#' @param tolerance Named tolerances, as for [.output_tolerance()].
#' @return A comparison `data.frame` with `output` and `n_values` columns.
compare_tables <- function(label, obtained, expected, tolerance = .output_default_tolerances) {
  rows <- lapply(names(expected), function(name) {
    want <- expected[[name]]
    have <- obtained[[name]]
    if (is.null(want) || nrow(want) == 0L) return(NULL)
    if (is.null(have)) {
      return(compare_values(sprintf("%s: %s is reported", label, name),
                            1, 0, tolerance = 0.5))
    }
    keys <- intersect(.output_keys(want, name), names(have))
    if (nrow(want) == 1L && nrow(have) == 1L) {
      # A one-row table (the criteria, the model row) is matched by position:
      # its text columns describe the fit, and one of them may be missing.
      want$.row <- 1L
      have$.row <- 1L
      keys <- ".row"
    }
    stopifnot("an expected table must share a key with the fitted one" =
                length(keys) > 0L)
    # Values are the numeric and logical columns; text columns label the row.
    values <- setdiff(intersect(names(want), names(have)),
                      union(keys, .output_keys(want, name)))
    values <- values[vapply(want[values], function(column) {
      is.numeric(column) || is.logical(column)
    }, logical(1))]
    as_key <- function(frame) {
      frame[keys] <- lapply(frame[keys], function(column) {
        if (is.numeric(column)) format(column, digits = 15, trim = TRUE) else
          as.character(column)
      })
      frame
    }
    joined <- merge(as_key(want[c(keys, values)]), as_key(have[c(keys, values)]),
                    by = keys, all.x = TRUE, suffixes = c(".reference", ".obtained"),
                    sort = TRUE)
    found <- stats::complete.cases(joined[keys]) &
      !is.na(match(do.call(paste, c(joined[keys], sep = "\r")),
                   do.call(paste, c(as_key(have)[keys], sep = "\r"))))
    present <- compare_values(sprintf("%s: %s rows the reference has, found in the fit",
                                      label, name),
                              nrow(want), sum(found), tolerance = 0.5)
    blocks <- lapply(values, function(column) {
      compare_block(sprintf("%s: %s$%s", label, name, column),
                    joined[[paste0(column, ".reference")]],
                    joined[[paste0(column, ".obtained")]],
                    tolerance = .output_tolerance(tolerance, name, column),
                    output = paste0(name, "$", column))
    })
    do.call(rbind, c(list(present), blocks))
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# The ledger
# ---------------------------------------------------------------------------

#' One case of a suite: its comparisons and its coverage ledger
#'
#' @param comparisons A comparison `data.frame`.
#' @param ledger A ledger `data.frame` from [output_ledger()], or `NULL`.
#' @return An object of class `harness_case`.
harness_case <- function(comparisons, ledger = NULL) {
  stopifnot("`comparisons` must be a data.frame" = is.data.frame(comparisons))
  structure(list(comparisons = comparisons, ledger = ledger), class = "harness_case")
}

#' Stack the cases of a suite
#' @param cases A list of `harness_case` objects or comparison frames.
#' @return One `harness_case`.
stack_cases <- function(cases) {
  as_case <- function(case) if (inherits(case, "harness_case")) case else harness_case(case)
  cases <- lapply(cases, as_case)
  comparisons <- do.call(rbind, lapply(cases, function(case) .with_output_columns(case$comparisons)))
  ledger <- do.call(rbind, lapply(cases, `[[`, "ledger"))
  harness_case(comparisons, ledger)
}

#' Give a comparison frame the `output` and `n_values` columns
#' @param comparisons A comparison `data.frame`.
#' @return The same frame with both columns present.
.with_output_columns <- function(comparisons) {
  if (is.null(comparisons$output)) comparisons$output <- rep(NA_character_, nrow(comparisons))
  if (is.null(comparisons$n_values)) comparisons$n_values <- rep(1L, nrow(comparisons))
  comparisons
}

#' Account for every output of one fit
#'
#' @param fit_label Label of the fit, as used in the comparisons.
#' @param inventory The frame from [fit_outputs()].
#' @param comparisons The comparisons made for this fit.
#' @param exempt Named character vector of reasons, keyed by `table$column` or
#'   by table, for outputs this reference cannot provide.
#' @return A `data.frame`, one row per output column: `fit`, `output`,
#'   `n_produced`, `n_compared`, `status` (`compared`, `partial`, `exempt` or
#'   `not covered`) and `reason`.
output_ledger <- function(fit_label, inventory, comparisons, exempt = character()) {
  stopifnot("`exempt` must be a named character vector of reasons" =
              is.character(exempt) && (length(exempt) == 0L || !is.null(names(exempt))))
  comparisons <- .with_output_columns(comparisons)
  named <- unique(stats::na.omit(comparisons$output))
  unknown <- setdiff(named, inventory$output)
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf("%s: compared outputs the fit does not report: %s",
                                fit_label, paste(unknown, collapse = ", ")),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  stale <- setdiff(names(exempt), c(inventory$output, inventory$table))
  if (length(stale) > 0L) {
    stop(errorCondition(sprintf("%s: exemptions for outputs the fit does not report: %s",
                                fit_label, paste(stale, collapse = ", ")),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  compared <- tapply(comparisons$n_values[!is.na(comparisons$output)],
                     comparisons$output[!is.na(comparisons$output)], sum)
  reasons <- c(exempt, .output_exemptions)
  reason_for <- function(output, table) {
    found <- reasons[c(output, table)]
    found <- found[!is.na(found)]
    if (length(found) == 0L) NA_character_ else unname(found[[1L]])
  }
  n_compared <- as.integer(ifelse(inventory$output %in% names(compared),
                                  compared[inventory$output], 0L))
  reason <- vapply(seq_len(nrow(inventory)), function(i) {
    reason_for(inventory$output[[i]], inventory$table[[i]])
  }, character(1))
  status <- ifelse(n_compared >= inventory$n_values & n_compared > 0L, "compared",
            ifelse(n_compared > 0L, "partial",
            ifelse(!is.na(reason), "exempt", "not covered")))
  data.frame(fit = fit_label, output = inventory$output,
             n_produced = inventory$n_values, n_compared = n_compared,
             status = status, reason = ifelse(status == "compared", NA_character_, reason),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Reference tables, written from their definitions
# ---------------------------------------------------------------------------

#' The reference a suite supplies, in one shape
#'
#' Every element is optional except the likelihood, the parameter count and the
#' posteriors; what is absent is not compared, and the ledger says so.
#' Classes must already be in the fit's order: see [align_reference()].
#'
#' @param log_likelihood,n_parameters The reference's maximum and its count.
#' @param posteriors Rows-by-profiles posterior matrix, rows in the fit's data
#'   order; the marginal profile posteriors for a two-level model.
#' @param rows Row numbers of `posteriors` in the fit's data, when the
#'   reference reports fewer rows than the data (default: every row).
#' @param groups Group identifier of each row of the data.
#' @param means,variances Profiles-by-indicators matrices with indicator names.
#' @param covariances Indicators-by-indicators-by-profiles array, for a full
#'   covariance model.
#' @param responses Named list of profiles-by-categories matrices whose column
#'   names are the fit's category labels.
#' @param standard_errors List with any of `mean`, `variance` (matrices like
#'   `means`) and `probability` (a list like `responses`).
#' @param profile_probabilities Group-classes-by-profiles matrix.
#' @param group_probabilities Group-class probabilities.
#' @param group_posteriors Groups-by-group-classes matrix, rows named by group.
#' @param group_log_likelihood Each group's log-likelihood contribution, named
#'   by group.
#' @param data The indicator data the reference was fitted to, in the fit's
#'   row order, for the residual and Gaussian diagnostic tables.
#' @param row_conditional_entropy One value per row of `posteriors`, in the
#'   same order: `E_C[H(X_i | C, y_i)] = -sum_h P(C = h | y_g) sum_k
#'   P(k | h, y_i) log P(k | h, y_i)`, the row's profile entropy within its
#'   group class, computed from the reference's own estimates. Needed for the
#'   entropy-charged criteria of a model with several group classes; with one
#'   group class it is the row posterior entropy and need not be given.
#' @return A list of class `harness_reference`.
harness_reference <- function(log_likelihood, n_parameters, posteriors, rows = NULL,
                              groups = NULL, means = NULL, variances = NULL,
                              covariances = NULL, responses = NULL,
                              standard_errors = NULL, profile_probabilities = NULL,
                              group_probabilities = NULL, group_posteriors = NULL,
                              group_log_likelihood = NULL, data = NULL,
                              row_conditional_entropy = NULL) {
  stopifnot(
    "`log_likelihood` must be one finite number" =
      is.numeric(log_likelihood) && length(log_likelihood) == 1L && is.finite(log_likelihood),
    "`n_parameters` must be one whole number" =
      is.numeric(n_parameters) && length(n_parameters) == 1L,
    "`posteriors` must be a numeric matrix" = is.matrix(posteriors) && is.numeric(posteriors),
    "`row_conditional_entropy` must be NULL or one finite value per row" =
      is.null(row_conditional_entropy) ||
      (is.numeric(row_conditional_entropy) &&
         length(row_conditional_entropy) == nrow(posteriors) &&
         all(is.finite(row_conditional_entropy))))
  structure(list(log_likelihood = log_likelihood, n_parameters = n_parameters,
                 posteriors = unname(posteriors),
                 rows = rows %||% seq_len(nrow(posteriors)), groups = groups,
                 means = means, variances = variances, covariances = covariances,
                 responses = responses, standard_errors = standard_errors,
                 profile_probabilities = profile_probabilities,
                 group_probabilities = group_probabilities,
                 group_posteriors = group_posteriors,
                 group_log_likelihood = group_log_likelihood, data = data,
                 row_conditional_entropy = row_conditional_entropy),
            class = "harness_reference")
}

#' Put a reference's classes in the fit's order
#'
#' Profiles are matched on the measurement parameters, then group classes on
#' the profile mix within each group class, each by exhaustive search.
#'
#' @param reference A `harness_reference`.
#' @param fit The fitted model.
#' @return The reference with every class-indexed element permuted.
align_reference <- function(reference, fit) {
  signature <- function(means, responses) {
    do.call(cbind, c(if (!is.null(means)) list(unname(as.matrix(means))),
                     lapply(responses, function(block) unname(as.matrix(block)))))
  }
  ours <- signature(if (!is.null(reference$means)) fit$means[, colnames(reference$means),
                                                              drop = FALSE],
                    fit$response_probabilities[names(reference$responses)])
  theirs <- signature(reference$means, reference$responses)
  order_profiles <- if (is.null(theirs)) seq_len(fit$n_profiles) else
    align_classes(ours, theirs)
  permute_rows <- function(m) if (is.null(m)) NULL else m[order_profiles, , drop = FALSE]
  reference$means <- permute_rows(reference$means)
  reference$variances <- permute_rows(reference$variances)
  if (!is.null(reference$covariances)) {
    reference$covariances <- reference$covariances[, , order_profiles, drop = FALSE]
  }
  reference$responses <- lapply(reference$responses, permute_rows)
  reference$standard_errors <- lapply(reference$standard_errors, function(block) {
    if (is.list(block)) lapply(block, permute_rows) else permute_rows(block)
  })
  reference$posteriors <- reference$posteriors[, order_profiles, drop = FALSE]
  if (!is.null(reference$profile_probabilities)) {
    reference$profile_probabilities <-
      reference$profile_probabilities[, order_profiles, drop = FALSE]
    if (nrow(reference$profile_probabilities) > 1L) {
      order_groups <- align_classes(unname(fit$profile_probabilities),
                                    unname(reference$profile_probabilities))
      reference$profile_probabilities <-
        reference$profile_probabilities[order_groups, , drop = FALSE]
      reference$group_probabilities <- reference$group_probabilities[order_groups]
      if (!is.null(reference$group_posteriors)) {
        reference$group_posteriors <- reference$group_posteriors[, order_groups, drop = FALSE]
      }
    }
  }
  reference
}

#' Classification entropy, `-sum p log p`, with `0 log 0 = 0`
#' @param p A posterior matrix.
#' @return A single number.
.reference_entropy <- function(p) {
  -sum(ifelse(p > 0, p * log(p), 0))
}

#' The classification tables of one level, from its posteriors
#'
#' Modal assignment takes the first class on a tie. Odds of correct
#' classification are `Inf` where the average posterior is one to within
#' `sqrt(.Machine$double.eps)`, the package's documented convention, and `NA`
#' where a class has no modal member or a share of zero or one.
#'
#' @param p Posteriors of the classified units.
#' @param level `"individuals"` or `"groups"`.
#' @param units Unit labels for the BCH table.
#' @return A named list of tables.
.reference_level_tables <- function(p, level, units) {
  k <- ncol(p)
  classes <- seq_len(k)
  modal <- max.col(p, ties.method = "first")
  n_modal <- tabulate(modal, nbins = k)
  average <- t(vapply(classes, function(assigned) {
    if (n_modal[[assigned]] == 0L) rep(NA_real_, k) else
      colMeans(p[modal == assigned, , drop = FALSE])
  }, numeric(k)))
  estimated <- colSums(p)
  share <- estimated / nrow(p)
  diagonal <- diag(average)
  odds <- ifelse(is.na(diagonal) | share <= 0 | share >= 1, NA_real_,
          ifelse(diagonal >= 1 - sqrt(.Machine$double.eps), Inf,
                 (diagonal / (1 - diagonal)) / (share / (1 - share))))
  errors <- t(vapply(classes, function(true_class) {
    vapply(classes, function(assigned) sum(p[modal == assigned, true_class]),
           numeric(1)) / sum(p[, true_class])
  }, numeric(k)))
  pairs <- expand.grid(assigned = classes, true = classes)
  inverse <- tryCatch(solve(errors), error = function(condition) NULL)
  list(
    counts = data.frame(level = level, class = classes, effective_count = estimated,
                        effective_proportion = share),
    classification = data.frame(level = level, class = classes, n_modal = n_modal,
                                proportion_modal = n_modal / nrow(p),
                                estimated_n = estimated, estimated_proportion = share,
                                average_posterior = diagonal,
                                odds_correct_classification = odds),
    average_posteriors = data.frame(level = level,
                                    assigned_class = rep(classes, each = k),
                                    class = rep(classes, times = k),
                                    n_assigned = rep(n_modal, each = k),
                                    average_posterior = as.vector(t(average))),
    classification_errors = data.frame(level = level, true_class = pairs$true,
                                       assigned_class = pairs$assigned,
                                       probability = errors[cbind(pairs$true, pairs$assigned)]),
    # A one-class level has the 1 x 1 error matrix [1], whose inverse weights
    # every unit by one.
    bch_weights = if (!is.null(inverse)) {
      data.frame(level = level, unit = rep(units, each = k),
                 assigned_class = rep(modal, each = k), class = rep(classes, length(modal)),
                 weight = as.vector(t(inverse[modal, , drop = FALSE])))
    },
    entropy = data.frame(level = level, n_classes = k, n_units = nrow(p),
                         entropy_sum = .reference_entropy(p),
                         relative_entropy = if (k > 1L)
                           1 - .reference_entropy(p) / (nrow(p) * log(k)) else NA_real_))
}

#' Every information criterion, from its definition
#'
#' Definitions as printed by `get_results(x, "information_criteria",
#' definitions = TRUE)`, and as in Nylund, Asparouhov and Muthen (2007),
#' Biernacki, Celeux and Govaert (2000) and Banfield and Raftery (1993). The
#' caller supplies each convention's entropy charge: the complete-data
#' entropy under both for a profile or class model, each level's own entropy
#' for a transition model (see [reference_tables()]).
#'
#' @return A one-row `data.frame`.
.reference_criteria <- function(log_likelihood, q, n_individuals, n_groups,
                                entropy_individuals, entropy_groups) {
  deviance <- -2 * log_likelihood
  by <- function(n, entropy) {
    c(bic = deviance + q * log(n), sabic = deviance + q * log((n + 2) / 24),
      caic = deviance + q * (log(n) + 1),
      awe = deviance + 2 * entropy + 2 * q * (1.5 + log(n)),
      icl = deviance + q * log(n) + 2 * entropy, clc = deviance + 2 * entropy)
  }
  groups <- by(n_groups, entropy_groups)
  individuals <- by(n_individuals, entropy_individuals)
  data.frame(log_likelihood = log_likelihood, n_parameters = q,
             aic = deviance + 2 * q, kic = deviance + 3 * (q + 1),
             bic_groups = groups[["bic"]], bic_individual = individuals[["bic"]],
             sabic_groups = groups[["sabic"]], sabic_individual = individuals[["sabic"]],
             caic_groups = groups[["caic"]], caic_individual = individuals[["caic"]],
             awe_groups = groups[["awe"]], awe_individual = individuals[["awe"]],
             icl_groups = groups[["icl"]], icl_individual = individuals[["icl"]],
             clc_groups = groups[["clc"]], clc_individual = individuals[["clc"]])
}

#' Within-profile residual association for every pair of indicators
#'
#' Gaussian pairs: the posterior-weighted residual correlation against the
#' model-implied one, Fisher-z statistic. Categorical pairs: the
#' posterior-weighted two-way table against the product of the profile's
#' margins, Pearson chi-square, both tables summarised by Cramer's V. Written
#' from the definitions in `?get_results`, not from the package's code.
#'
#' @return A `data.frame` keyed by `profile`, `indicator_1`, `indicator_2`.
.reference_residuals <- function(reference, continuous, categorical, levels) {
  p <- reference$posteriors
  data <- reference$data
  cramer <- function(table) {
    total <- sum(table)
    smaller <- min(dim(table)) - 1L
    if (total <= 0 || smaller < 1L) return(NA_real_)
    independent <- outer(rowSums(table), colSums(table)) / total
    positive <- independent > 0
    sqrt(sum((table[positive] - independent[positive])^2 / independent[positive]) /
           (total * smaller))
  }
  clamp <- function(r) pmin(pmax(r, -0.999999), 0.999999)
  rows <- lapply(seq_len(ncol(p)), function(k) {
    label <- sprintf("profile_%d", k)
    gaussian <- if (length(continuous) >= 2L) {
      pairs <- utils::combn(continuous, 2L)
      lapply(seq_len(ncol(pairs)), function(index) {
        pair <- pairs[, index]
        complete <- stats::complete.cases(data[pair])
        w <- p[complete, k]
        centred <- sweep(as.matrix(data[complete, pair]), 2L, reference$means[k, pair], "-")
        observed_cov <- crossprod(centred * w, centred)
        implied <- if (!is.null(reference$covariances))
          reference$covariances[pair, pair, k] else diag(reference$variances[k, pair])
        observed <- observed_cov[1L, 2L] / sqrt(prod(diag(observed_cov)))
        expected <- implied[1L, 2L] / sqrt(prod(diag(implied)))
        effective <- sum(w)
        statistic <- (atanh(clamp(observed)) - atanh(clamp(expected))) * sqrt(effective - 3)
        data.frame(profile = label, indicator_1 = pair[[1L]], indicator_2 = pair[[2L]],
                   kind = "gaussian", observed = observed, expected = expected,
                   residual = observed - expected, effective_n = effective,
                   statistic = statistic, df = NA_integer_,
                   p_value = 2 * stats::pnorm(-abs(statistic)))
      })
    }
    nominal <- if (length(categorical) >= 2L) {
      pairs <- utils::combn(categorical, 2L)
      lapply(seq_len(ncol(pairs)), function(index) {
        first <- pairs[1L, index]
        second <- pairs[2L, index]
        complete <- !is.na(data[[first]]) & !is.na(data[[second]])
        w <- p[complete, k]
        a <- factor(as.character(data[[first]][complete]), levels[[first]])
        b <- factor(as.character(data[[second]][complete]), levels[[second]])
        observed <- as.matrix(stats::xtabs(w ~ a + b))
        expected <- sum(w) * outer(reference$responses[[first]][k, levels[[first]]],
                                   reference$responses[[second]][k, levels[[second]]])
        positive <- expected > 0
        statistic <- if (any(!positive & observed > 0)) Inf else
          sum((observed[positive] - expected[positive])^2 / expected[positive])
        degrees <- (nrow(observed) - 1L) * (ncol(observed) - 1L)
        data.frame(profile = label, indicator_1 = first, indicator_2 = second,
                   kind = "categorical", observed = cramer(observed),
                   expected = cramer(expected),
                   residual = cramer(observed) - cramer(expected), effective_n = sum(w),
                   statistic = statistic, df = as.integer(degrees),
                   p_value = stats::pchisq(statistic, degrees, lower.tail = FALSE))
      })
    }
    pieces <- c(gaussian, nominal)
    if (length(pieces) == 0L) NULL else do.call(rbind, pieces)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(NULL)
  # The default adjustment is "none", so the adjusted p-value is the raw one.
  out$p_adjusted <- out$p_value
  out
}

#' Covariance and support diagnostics of each Gaussian profile, by definition
#'
#' Mass over rows with a Gaussian value; Kish's effective n
#' of the posteriors; the smallest variance and the smallest eigenvalue of the
#' covariance, each relative to the observed marginal variance (divisor n);
#' and whether that eigenvalue sits at the variance floor.
#'
#' `positive_rows` is the exact value: a Gaussian density is positive
#' everywhere, so every row with a Gaussian value has a positive posterior in
#' every component. Counting computed posteriors above zero instead counts the
#' rows whose posterior did not underflow, which differs between programs (on
#' thyroid VVV the package keeps a posterior of 1.2e-320 that mclust rounds to
#' zero). Disagreement here is a finding about that column, not tolerance to
#' tune.
#'
#' @return A `data.frame` keyed by `context` and `profile`.
.reference_gaussian_diagnostics <- function(reference, continuous, min_variance) {
  data <- as.matrix(reference$data[continuous])
  p <- reference$posteriors
  with_value <- rowSums(!is.na(data)) > 0L
  overall <- apply(data, 2L, function(column) {
    column <- column[!is.na(column)]
    mean((column - mean(column))^2)
  })
  do.call(rbind, lapply(seq_len(ncol(p)), function(k) {
    contribution <- p[, k] * with_value
    positive <- contribution > 0
    scaled <- contribution[positive] / max(contribution)
    covariance <- if (!is.null(reference$covariances))
      matrix(reference$covariances[continuous, continuous, k], length(continuous)) else
      diag(reference$variances[k, continuous], length(continuous))
    standardized <- covariance / sqrt(outer(overall, overall))
    data.frame(context = "profiles", profile = sprintf("profile_%d", k),
               gaussian_mass = sum(contribution), positive_rows = sum(with_value),
               weight_effective_n = sum(scaled)^2 / sum(scaled^2),
               min_variance_ratio = min(diag(covariance) / overall),
               min_standardized_eigenvalue = min(eigen(standardized, symmetric = TRUE,
                                                       only.values = TRUE)$values),
               boundary = min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <=
                 min_variance * (1 + 1e-8))
  }))
}

#' Build every table a reference can supply, in the fit's shapes
#'
#' @param reference An aligned `harness_reference`.
#' @param fit The fitted model, for its dimensions, labels and settings only;
#'   no estimate is read from it.
#' @return A named list of tables, with an `exempt` attribute naming (with
#'   reasons) the outputs left out because the reference cannot supply them.
reference_tables <- function(reference, fit) {
  stopifnot("`reference` must come from harness_reference()" =
              inherits(reference, "harness_reference"),
            "sampling weights are not handled by the reference tables" =
              is.null(fit$sampling_weights))
  p <- reference$posteriors
  k <- ncol(p)
  rows <- reference$rows
  n_rows <- nrow(p)
  data <- reference$data
  continuous <- intersect(fit$vars, colnames(reference$means) %||% character())
  categorical <- names(reference$responses) %||% character()
  informative <- if (is.null(data)) rep(TRUE, n_rows) else
    rowSums(!is.na(data[intersect(fit$vars, names(data))])) > 0L
  tables <- list()

  if (!is.null(reference$means)) {
    grid <- expand.grid(profile = seq_len(k), indicator = continuous,
                        stringsAsFactors = FALSE)
    at <- cbind(grid$profile, match(grid$indicator, colnames(reference$means)))
    variance <- if (!is.null(reference$covariances)) {
      vapply(seq_len(nrow(grid)), function(i) {
        reference$covariances[grid$indicator[[i]], grid$indicator[[i]], grid$profile[[i]]]
      }, numeric(1))
    } else reference$variances[, continuous, drop = FALSE][at]
    profiles <- data.frame(grid, mean = reference$means[, continuous, drop = FALSE][at],
                           variance = variance, standard_deviation = sqrt(variance))
    se <- reference$standard_errors
    if (!is.null(se$mean)) profiles$mean_standard_error <- se$mean[, continuous, drop = FALSE][at]
    if (!is.null(se$variance)) {
      profiles$variance_standard_error <- se$variance[, continuous, drop = FALSE][at]
    }
    tables$profiles <- profiles
    covariance_of <- function(a, b, profile) {
      if (!is.null(reference$covariances)) return(reference$covariances[a, b, profile])
      if (identical(a, b)) reference$variances[profile, a] else 0
    }
    pairs <- expand.grid(profile = seq_len(k), indicator = continuous,
                         indicator_2 = continuous, stringsAsFactors = FALSE)
    tables$covariances <- data.frame(pairs, covariance = vapply(seq_len(nrow(pairs)), function(i) {
      covariance_of(pairs$indicator[[i]], pairs$indicator_2[[i]], pairs$profile[[i]])
    }, numeric(1)))
  }

  if (length(categorical) > 0L) {
    tables$responses <- do.call(rbind, lapply(categorical, function(item) {
      block <- reference$responses[[item]]
      se <- reference$standard_errors$probability[[item]]
      cumulative <- t(apply(block, 1L, cumsum))
      threshold <- stats::qlogis(pmin(pmax(cumulative, 0), 1))
      threshold[, ncol(block)] <- NA_real_
      out <- data.frame(profile = rep(seq_len(k), ncol(block)), indicator = item,
                        category = rep(colnames(block), each = k),
                        probability = as.vector(block), threshold = as.vector(threshold),
                        stringsAsFactors = FALSE)
      if (!is.null(se)) out$probability_standard_error <- as.vector(se)
      out
    }))
  }

  gamma <- reference$profile_probabilities %||% matrix(colMeans(p), 1L)
  n_group_classes <- nrow(gamma)
  delta <- reference$group_probabilities %||% rep(1, n_group_classes)
  tables$profile_probabilities <- data.frame(
    group_class = rep(seq_len(n_group_classes), times = k),
    profile = rep(seq_len(k), each = n_group_classes),
    probability = as.vector(gamma), group_class_probability = rep(delta, times = k))

  group_of_row <- reference$groups %||% rows
  group_ids <- unique(group_of_row)
  zg <- reference$group_posteriors %||%
    matrix(1, length(group_ids), 1L, dimnames = list(as.character(group_ids), NULL))
  zg <- zg[as.character(group_ids), , drop = FALSE]
  informative_group <- as.character(group_ids) %in%
    as.character(unique(group_of_row[rows[informative]]))

  tables$posteriors <- data.frame(row = rep(rows, times = k), profile = rep(seq_len(k), each = n_rows),
                                  posterior = as.vector(p),
                                  modal = as.vector(col(p) == max.col(p, "first")))
  group_size <- as.vector(table(factor(as.character(group_of_row),
                                       as.character(group_ids))))
  tables$group_posteriors <- data.frame(
    group = rep(as.character(group_ids), times = ncol(zg)),
    group_class = rep(seq_len(ncol(zg)), each = nrow(zg)),
    group_size = rep(group_size, times = ncol(zg)),
    posterior = as.vector(zg), modal = as.vector(col(zg) == max.col(zg, "first")))
  if (!is.null(reference$group_log_likelihood)) {
    tables$group_posteriors$log_likelihood <-
      unname(reference$group_log_likelihood[tables$group_posteriors$group])
  }
  modal_group <- max.col(zg, "first")
  names(modal_group) <- as.character(group_ids)
  assignments <- data.frame(.row = rows, profile = max.col(p, "first"),
                            group_class = unname(modal_group[as.character(group_of_row[rows])]),
                            uncertainty = 1 - apply(p, 1L, max))
  assignments[sprintf("posterior_profile_%d", seq_len(k))] <- as.data.frame(p)
  tables$assignments <- assignments

  individual <- .reference_level_tables(p[informative, , drop = FALSE], "individuals",
                                        rows[informative])
  group <- .reference_level_tables(zg[informative_group, , drop = FALSE], "groups",
                                   as.character(group_ids)[informative_group])
  for (name in names(individual)) {
    tables[[name]] <- rbind(individual[[name]], group[[name]])
  }

  n_individuals <- sum(informative)
  n_groups <- sum(informative_group)
  # ICL, AWE and CLC charge the entropy of the complete-data classification,
  # H(group classes) + sum_i E_C[H(X_i | C, y_i)] over the classified groups
  # and rows, under both conventions; with one group class that is the row
  # posterior entropy. A transition model keeps each level's own entropy.
  # Without the rows' within-class entropies a model with several group
  # classes has no reference for the three, and they are exempted by name, so
  # they can never agree by default.
  exempt <- character()
  entropies <- if (inherits(fit, "multilpa_transitions")) {
    c(individuals = .reference_entropy(p[informative, , drop = FALSE]),
      groups = .reference_entropy(zg[informative_group, , drop = FALSE]))
  } else {
    complete <- if (ncol(zg) == 1L) {
      .reference_entropy(p[informative, , drop = FALSE])
    } else if (!is.null(reference$row_conditional_entropy)) {
      .reference_entropy(zg[informative_group, , drop = FALSE]) +
        sum(reference$row_conditional_entropy[informative])
    } else {
      charged <- paste0("information_criteria$",
                        c(outer(c("awe", "icl", "clc"), c("groups", "individual"),
                                paste, sep = "_")))
      exempt <- stats::setNames(rep(paste(
        "the reference supplies no row entropies within the group classes",
        "(`row_conditional_entropy`), which the complete-data entropy of a",
        "model with several group classes needs"), length(charged)), charged)
      NA_real_
    }
    c(individuals = complete, groups = complete)
  }
  tables$information_criteria <- .reference_criteria(
    reference$log_likelihood, reference$n_parameters, n_individuals, n_groups,
    entropies[["individuals"]], entropies[["groups"]])
  tables$model <- data.frame(n_observations = n_rows, n_informative = n_individuals,
                             n_groups = length(group_ids), n_profiles = k,
                             n_group_classes = n_group_classes,
                             n_parameters = reference$n_parameters,
                             # Nothing is held fixed in a reference fit, so its
                             # count is the count with the measurement model.
                             n_parameters_with_measurement = reference$n_parameters,
                             log_likelihood = reference$log_likelihood,
                             aic = tables$information_criteria$aic,
                             bic_groups = tables$information_criteria$bic_groups,
                             bic_individual = tables$information_criteria$bic_individual,
                             .row = 1L)

  if (!is.null(data)) {
    tables$residuals <- .reference_residuals(reference, continuous, categorical,
                                             fit$categorical_levels)
    if (length(continuous) > 0L) {
      tables$gaussian_diagnostics <- .reference_gaussian_diagnostics(
        reference, continuous, fit$min_variance %||% 1e-6)
    }
  }
  # Outputs this reference cannot speak to are left out here, and named in the
  # `exempt` attribute, which `compare_fit()` adds to the suite's exemptions
  # so that the ledger gives the reason.
  tables <- .drop_exempt(Filter(Negate(is.null), tables), exempt)
  attr(tables, "exempt") <- exempt
  tables
}

#' Leave exempted outputs out of the expected tables
#'
#' An exemption says the reference cannot speak to an output, so the output is
#' not compared at all; the ledger records the reason.
#'
#' @param expected Named list of expected tables.
#' @param exempt Named reasons, keyed by table or `table$column`.
#' @return The tables without the exempted ones and columns.
.drop_exempt <- function(expected, exempt) {
  keys <- names(exempt) %||% character()
  expected <- expected[setdiff(names(expected), keys)]
  lapply(stats::setNames(names(expected), names(expected)), function(name) {
    columns <- sub("^[^$]*[$]", "", keys[startsWith(keys, paste0(name, "$"))])
    expected[[name]][setdiff(names(expected[[name]]), columns)]
  })
}

#' Compare a fit with a reference on every output, and account for the rest
#'
#' @param label Label for the fit's rows.
#' @param fit The fitted model.
#' @param reference A `harness_reference`, in its own class order.
#' @param tolerance Named tolerances, merged over [.output_default_tolerances].
#' @param exempt Named reasons for outputs this reference cannot provide.
#' @param extra Further comparison rows for this fit (a frame from
#'   [compare_values()] or [compare_block()]), which count towards the ledger
#'   through their `output` column.
#' @return A `harness_case`.
compare_fit <- function(label, fit, reference, tolerance = numeric(),
                        exempt = character(), extra = NULL) {
  tolerances <- .output_default_tolerances
  tolerances[names(tolerance)] <- tolerance
  aligned <- align_reference(reference, fit)
  obtained <- fit_tables(fit)
  obtained$model$.row <- 1L
  built <- reference_tables(aligned, fit)
  # The reference's own exemptions join the suite's, so the ledger gives
  # their reason; a suite's own reason for the same output takes precedence.
  own <- attr(built, "exempt") %||% character()
  exempt <- c(exempt, own[setdiff(names(own), names(exempt))])
  expected <- .drop_exempt(built, exempt)
  comparisons <- compare_tables(label, obtained, expected, tolerances)
  if (!is.null(extra)) comparisons <- rbind(comparisons, .with_output_columns(extra))
  harness_case(comparisons, output_ledger(label, fit_outputs(fit), comparisons, exempt))
}
