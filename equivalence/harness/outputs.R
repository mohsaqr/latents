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
#   compare_fit_tables() the same for a fit whose expected tables the suite
#                       writes itself (families reference_tables() does not
#                       build: additive, cross-level, general transition).
#   compare_verb(),     the same ledger for tables a verb returns rather than
#   table_ledger()      a fit's catalogue (three_step, r3step,
#                       pool_imputations, predict, parameter_inference, the
#                       enumeration's lmr table). keyed_table() declares the
#                       keys of a table whose columns the catalogue rules
#                       would misread; every other column is then a value.
#
# reference_tables() handles sampling weights (the package's documented
# rescaling and positive-support convention), membership coefficients,
# ordinal and count indicators, and their residual kinds.
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
  `replication$n_starts` = "this program's completed-start budget (unit-tested)",
  `replication$n_screened` = "this program's screened-start budget (unit-tested)",
  `replication$n_converged` = "this optimiser's converged-start count (unit-tested)",
  `replication$n_failed` = "this optimiser's failed-start count (unit-tested)",
  `replication$n_best_replicated` = "this optimiser's repeated-best count (unit-tested)",
  `replication$share` = "this optimiser's repeated-best share (unit-tested)",
  `replication$tolerance` = "this program's threshold for counting repeated starts (unit-tested)",
  `replication$best_start` = "which of this program's starts won (unit-tested)",
  `model$converged` = "this optimiser's own stopping flag",
  `model$iterations` = "this optimiser's iteration count",
  `model$best_start` = "which of this program's random starts won",
  `model$n_best_replicated` = "how many of this program's random starts reached the best value",
  `model$boundary` = "this package's warning flag, with its own threshold (unit-tested)",
  `model$small_classes` = "this package's warning flag, with its own threshold (unit-tested)",
  # A membership-covariate fit's one-row table.
  `model$extreme_logits` = "this package's warning flag, with its own threshold (unit-tested)",
  `model$n_starts` = "this program's completed-start budget (unit-tested)",
  # The one-row `fit` table of the additive, cross-level and general
  # transition families carries the same bookkeeping under another name.
  `fit$converged` = "this optimiser's own stopping flag",
  `fit$iterations` = "this optimiser's iteration count",
  `fit$boundary` = "this package's warning flag, with its own threshold (unit-tested)",
  `fit$kkt` = "this optimiser's check of its own first-order conditions (unit-tested)",
  `fit$n_starts` = "this program's completed-start budget (unit-tested)",
  `fit$n_best_replicated` = "how many of this program's random starts reached the best value")

#' Columns named like an identifier that are values in a particular table
#'
#' `time` is the occasion key of a transition fit's `sequences` table and the
#' random-slope prediction of a growth fit's `random_effects` table.
.output_value_columns <- list(random_effects = "time")

#' The key columns of one output table
#'
#' A table that declares its keys ([keyed_table()]) is keyed by them alone.
#' Otherwise the keys are the identifier-named columns and the text columns.
#'
#' @param table A table from `get_results()`.
#' @param name Its catalogue name.
#' @return A character vector of column names.
.output_keys <- function(table, name) {
  declared <- attr(table, "harness_keys")
  if (!is.null(declared)) return(declared)
  if (identical(name, "assignments")) return(".row")
  labels <- names(table)[!vapply(table, function(column) {
    is.numeric(column) || is.logical(column)
  }, logical(1))]
  setdiff(union(intersect(names(table), .output_identifiers), labels),
          .output_value_columns[[name]] %||% character())
}

#' The value columns of one output table
#'
#' Without declared keys, a text column labels its row, so the values are the
#' numeric and logical columns that are not keys. A table with declared keys
#' has no labelling columns: every other column is a value, a text one
#' included (a modal class reported as `"profile_2"`, a note), and is compared
#' exactly.
#'
#' @param table A table.
#' @param name Its catalogue name.
#' @return A character vector of column names.
.output_values <- function(table, name) {
  candidates <- setdiff(names(table), c(.output_keys(table, name), ".row"))
  keep <- if (is.null(attr(table, "harness_keys"))) {
    vapply(table[candidates], function(column) is.numeric(column) || is.logical(column),
           logical(1))
  } else {
    !vapply(table[candidates], is.list, logical(1))
  }
  candidates[keep]
}

#' Declare the columns that identify a row of a table
#'
#' Verb results are not catalogue tables, and their columns do not follow the
#' catalogue's naming: `predict()` reports the modal profile in a `profile`
#' column, which the catalogue rules would take for a key, and
#' `pool_imputations()` keys its per-imputation rows by an integer
#' `imputation`. Declaring the keys makes every other column a value that has
#' to be compared or exempted.
#'
#' @param table A `data.frame`, possibly with a class of its own.
#' @param keys Names of the columns that identify a row, jointly unique.
#' @return A plain `data.frame` (classes and attributes dropped) carrying the
#'   keys in its `harness_keys` attribute.
keyed_table <- function(table, keys) {
  stopifnot("`table` must be a data.frame" = is.data.frame(table),
            "`keys` must name at least one column" =
              is.character(keys) && length(keys) > 0L && !anyNA(keys))
  absent <- setdiff(keys, names(table))
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf("keys not in the table: %s", paste(absent, collapse = ", ")),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  plain <- list2DF(lapply(as.list(table), function(column) {
    if (is.factor(column)) as.character(column) else column
  }), nrow = nrow(table))
  .check_unique_keys(plain, keys, "the table")
  attr(plain, "harness_keys") <- keys
  plain
}

#' Refuse keys that do not identify the rows
#'
#' A join on keys that repeat pairs rows arbitrarily, and a comparison over
#' it would mean nothing.
#'
#' @param table A table.
#' @param keys Its key columns.
#' @param what A description for the message.
#' @return `NULL`, invisibly.
.check_unique_keys <- function(table, keys, what) {
  absent <- setdiff(keys, names(table))
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf("%s lacks the key columns %s", what,
                                paste(absent, collapse = ", ")),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  if (nrow(table) > 0L && anyDuplicated(table[keys]) > 0L) {
    stop(errorCondition(sprintf("%s: the keys %s do not identify its rows", what,
                                paste(keys, collapse = ", ")),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  invisible(NULL)
}

#' Every table a fit reports, as compared here
#'
#' `assignments` gains a `.row` key, and loses the columns that echo the input
#' (the indicators and the identifier), which are data rather than output.
#'
#' The additive, cross-level and general transition families report the modal
#' class as a label (`"group_class_2"`), not a number. Under the catalogue
#' rules a text column labels its row, so those labels would never be
#' compared: in `assignments` (keyed by position) they would not be counted at
#' all. Their `assignments` and `groups` tables therefore declare their keys
#' (position; the group identifier), which makes the labels values that are
#' compared exactly. Every other family keeps the catalogue rules.
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
  if (inherits(fit, .labelled_class_families)) {
    if (!is.null(tables$assignments)) {
      tables$assignments <- keyed_table(tables$assignments, ".row")
    }
    if (!is.null(tables$groups)) tables$groups <- keyed_table(tables$groups, fit$id)
  }
  tables
}

#' Families whose unit tables report the modal class as a text label
.labelled_class_families <- c("multilpa_additive", "multilpa_cross_level", "multilpa_lta")

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
    values <- .output_values(table, name)
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
#' Logical columns must match exactly. Text columns (a modal class label, a
#' note) must match exactly too, missing in both counting as a match; their
#' row reports `reference = 1`, `obtained` 1 or 0 for the worst cell matching
#' or not, and quotes the first mismatching pair in `quantity`.
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
  textual <- function(x) is.character(x) || is.factor(x)
  if (textual(reference) || textual(obtained)) {
    return(.compare_text_block(quantity, as.character(reference), as.character(obtained),
                               output))
  }
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

#' Compare every cell of a text column exactly
#' @param quantity Label for the row.
#' @param reference,obtained Character vectors of equal length.
#' @param output The `table$column` key, or `NA`.
#' @return A one-row `data.frame` shaped like [compare_block()]'s.
.compare_text_block <- function(quantity, reference, obtained, output) {
  n <- length(reference)
  matched <- (is.na(reference) & is.na(obtained)) |
    (!is.na(reference) & !is.na(obtained) & reference == obtained)
  first <- which(!matched)[1L]
  note <- if (is.na(first)) "" else
    sprintf(", first mismatch %s against %s", encodeString(reference[[first]], quote = '"'),
            encodeString(obtained[[first]], quote = '"'))
  data.frame(quantity = sprintf("%s (%d values, text%s)", quantity, n, note),
             reference = 1, obtained = if (is.na(first)) 1 else 0,
             difference = if (is.na(first)) 0 else 1, tolerance = 0.5,
             agrees = n > 0L && is.na(first), output = output, n_values = n,
             stringsAsFactors = FALSE)
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
  # The replication table repeats the fitted likelihood. Check it against
  # the independent reference too; only search bookkeeping is exempt. Some
  # reference suites call their one-row model table `fit`, others `model`.
  if (!is.null(obtained$replication) && is.null(expected$replication)) {
    likelihood <- expected$model$log_likelihood %||% expected$fit$log_likelihood %||%
      expected$information_criteria$log_likelihood
    if (length(likelihood) == 1L) {
      expected$replication <- data.frame(log_likelihood = likelihood)
    }
  }
  rows <- lapply(names(expected), function(name) {
    want <- expected[[name]]
    have <- obtained[[name]]
    if (is.null(want) || nrow(want) == 0L) return(NULL)
    if (is.null(have)) {
      return(compare_values(sprintf("%s: %s is reported", label, name),
                            1, 0, tolerance = 0.5))
    }
    # Keys a table declares (keyed_table()) bind both sides, so the
    # comparison takes the same value columns the inventory counts.
    declared <- attr(have, "harness_keys") %||% attr(want, "harness_keys")
    if (!is.null(declared)) {
      attr(want, "harness_keys") <- declared
      .check_unique_keys(want, declared, sprintf("%s: the expected %s", label, name))
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
    # Values are the numeric and logical columns, text columns labelling the
    # row, unless the table declares its keys (see .output_values()).
    values <- setdiff(intersect(.output_values(want, name), names(have)), keys)
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
    # A missing key is a value like any other (parameter_inference() gives
    # group-class probabilities `term = NA`): merge() pairs NA with NA, and the
    # pasted keys match "NA" with "NA", so such a row counts as found.
    found <- !is.na(match(do.call(paste, c(joined[keys], sep = "\r")),
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
# Verb results: the same ledger for tables that are not a fit's catalogue
# ---------------------------------------------------------------------------
#
# three_step(), r3step(), pool_imputations(), predict(),
# parameter_inference() and get_results(enumeration, "lmr") return tables
# that get_results(fit, "all") never lists, so fit_outputs() cannot inventory
# them and their comparisons would sit outside the not-covered gate. A suite
# names those tables itself, and every value column of each is accounted for
# exactly as a fit's are: compared in full, exempted with a reason, or a
# failure of the run.

#' Account for every value of verb tables
#'
#' @param label Label of the verb call, as used in the comparisons (it fills
#'   the ledger's `fit` column).
#' @param tables Named list of the verb's tables, names being how the
#'   comparisons' `output` keys call them (`"three_step"` for
#'   `"three_step$estimate"`). Declare the keys of a table whose columns do
#'   not follow the catalogue naming with [keyed_table()].
#' @param comparisons The comparisons made for these tables.
#' @param exempt Named reasons, keyed by `table$column` or by table.
#' @return A ledger `data.frame` like [output_ledger()]'s.
table_ledger <- function(label, tables, comparisons, exempt = character()) {
  if (!is.list(tables) || is.data.frame(tables) || length(tables) == 0L ||
      is.null(names(tables)) || any(!nzchar(names(tables))) ||
      !all(vapply(tables, is.data.frame, logical(1)))) {
    stop(errorCondition(sprintf("%s: `tables` must be a named list of data.frames", label),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  # An empty result has nothing to account for, and would pass vacuously.
  empty <- names(tables)[vapply(tables, function(table) nrow(table) == 0L, logical(1))]
  if (length(empty) > 0L) {
    stop(errorCondition(sprintf("%s: verb tables with no rows cannot be accounted for: %s",
                                label, paste(empty, collapse = ", ")),
                        class = "multilpa_bad_ledger", call = NULL))
  }
  output_ledger(label, fit_outputs(tables), comparisons, exempt)
}

#' Compare verb tables with expected ones, and account for every value
#'
#' @param label Label for the rows.
#' @param obtained Named list of the verb's tables (see [table_ledger()]).
#' @param expected Named list of tables in the same shapes, built from the
#'   reference; a table's declared keys are taken from `obtained`.
#' @param tolerance Named tolerances (`table$column`, `table`, `default`),
#'   merged over [.output_default_tolerances].
#' @param exempt Named reasons for outputs the reference cannot provide.
#' @param extra Further comparison rows, counted through their `output`.
#' @return A `harness_case`.
compare_verb <- function(label, obtained, expected, tolerance = numeric(),
                         exempt = character(), extra = NULL) {
  .compare_case(label, obtained, expected, tolerance, exempt, extra,
                inventory = NULL)
}

#' Compare a fit of any family with tables a suite built itself
#'
#' The path for the families [reference_tables()] does not build (the
#' additive, cross-level and general transition families, growth and
#' regression mixtures): the suite writes the expected tables, and every
#' output of `fit_tables(fit)` is accounted for as [compare_fit()] does.
#'
#' @inheritParams compare_verb
#' @param fit The fitted model.
#' @return A `harness_case`.
compare_fit_tables <- function(label, fit, expected, tolerance = numeric(),
                               exempt = character(), extra = NULL) {
  obtained <- fit_tables(fit)
  .compare_case(label, obtained, expected, tolerance, exempt, extra,
                inventory = fit_outputs(obtained))
}

#' The comparison and ledger shared by the verb and fit paths
#' @param inventory `NULL` for verb tables (inventoried by [table_ledger()]),
#'   or a fit's inventory from [fit_outputs()].
#' @return A `harness_case`.
.compare_case <- function(label, obtained, expected, tolerance, exempt, extra, inventory) {
  tolerances <- .output_default_tolerances
  tolerances[names(tolerance)] <- tolerance
  expected <- .drop_exempt(expected, exempt)
  comparisons <- compare_tables(label, obtained, expected, tolerances)
  if (!is.null(extra)) comparisons <- rbind(comparisons, .with_output_columns(extra))
  ledger <- if (is.null(inventory)) table_ledger(label, obtained, comparisons, exempt) else
    output_ledger(label, inventory, comparisons, exempt)
  harness_case(comparisons, ledger)
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
#'   `means`), `probability` (a list like `responses`), `location` (a
#'   profiles-by-ordinal-indicators matrix), `count_mean` and
#'   `count_dispersion` (matrices like `count_means`).
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
#' @param weights For a fit with sampling weights: one raw sampling weight per
#'   independent unit (each group of a two-level fit, each row of a
#'   single-level one), named by the unit's group label or in the order the
#'   units first appear. [reference_tables()] applies the package's documented
#'   convention (`?multilpa`, `weights`): weights rescaled to sum to the number
#'   of units of positive weight, and a unit of weight zero in no sample size,
#'   entropy or classification table, though its posteriors are reported.
#'   `log_likelihood` must then be the pseudo log likelihood under those
#'   rescaled weights, `sum_j w_j log L_j`.
#' @param coefficients For a membership-covariate fit: a list with `profile`
#'   (terms by profiles) and, for group covariates, `group` (terms by group
#'   classes) matrices of multinomial logits, one column per class in the
#'   reference's class order and against any baseline; rows named by the fit's
#'   term labels, where a term naming a group class (`"x:group_class_2"`,
#'   `profile_slopes = "group_class"`) uses the reference's group-class
#'   numbering. Logits are identified up to a constant per term, so
#'   [reference_tables()] re-expresses them against the fit's baseline, its
#'   last class.
#' @param ordinal Named list of profiles-by-categories probability matrices of
#'   the ordinal indicators, column names the fit's category labels in order.
#' @param count_means Profiles-by-indicators matrix of count means, with
#'   indicator names.
#' @param count_dispersion Matrix like `count_means` of negative-binomial
#'   dispersions (`alpha` of `mu + alpha mu^2`).
#' @return A list of class `harness_reference`.
harness_reference <- function(log_likelihood, n_parameters, posteriors, rows = NULL,
                              groups = NULL, means = NULL, variances = NULL,
                              covariances = NULL, responses = NULL,
                              standard_errors = NULL, profile_probabilities = NULL,
                              group_probabilities = NULL, group_posteriors = NULL,
                              group_log_likelihood = NULL, data = NULL,
                              row_conditional_entropy = NULL, weights = NULL,
                              coefficients = NULL, ordinal = NULL, count_means = NULL,
                              count_dispersion = NULL) {
  matrix_list <- function(x) {
    is.null(x) || (is.list(x) && !is.null(names(x)) &&
                     all(vapply(x, function(block) is.null(block) ||
                                  (is.matrix(block) && is.numeric(block)), logical(1))))
  }
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
         all(is.finite(row_conditional_entropy))),
    "`weights` must be NULL or finite non-negative weights, some positive" =
      is.null(weights) || (is.numeric(weights) && all(is.finite(weights)) &&
                             all(weights >= 0) && any(weights > 0)),
    "`coefficients` must be a named list of numeric matrices with row names" =
      matrix_list(coefficients) &&
      all(vapply(coefficients, function(block) !is.null(rownames(block)), logical(1))),
    "`ordinal` must be a named list of numeric matrices with category names" =
      matrix_list(ordinal) &&
      all(vapply(ordinal, function(block) !is.null(colnames(block)), logical(1))),
    "`count_means` must be a numeric matrix with indicator names" =
      is.null(count_means) || (is.matrix(count_means) && !is.null(colnames(count_means))),
    "`count_dispersion` must match `count_means`" =
      is.null(count_dispersion) ||
      (!is.null(count_means) && identical(dim(count_dispersion), dim(count_means))))
  structure(list(log_likelihood = log_likelihood, n_parameters = n_parameters,
                 posteriors = unname(posteriors),
                 rows = rows %||% seq_len(nrow(posteriors)), groups = groups,
                 means = means, variances = variances, covariances = covariances,
                 responses = responses, standard_errors = standard_errors,
                 profile_probabilities = profile_probabilities,
                 group_probabilities = group_probabilities,
                 group_posteriors = group_posteriors,
                 group_log_likelihood = group_log_likelihood, data = data,
                 row_conditional_entropy = row_conditional_entropy, weights = weights,
                 coefficients = coefficients, ordinal = ordinal, count_means = count_means,
                 count_dispersion = count_dispersion),
            class = "harness_reference")
}

#' Put a reference's classes in the fit's order
#'
#' Profiles are matched on the measurement parameters (means, response,
#' ordinal category and count-mean blocks), then group classes on the profile
#' mix within each group class, each by exhaustive search. A
#' membership-covariate fit has no profile mix, and its group classes are
#' matched on the group posteriors instead. The fit's estimates are read for
#' this matching only.
#'
#' @param reference A `harness_reference`.
#' @param fit The fitted model.
#' @return The reference with every class-indexed element permuted, and the
#'   group-class labels inside coefficient term names renumbered.
align_reference <- function(reference, fit) {
  signature <- function(blocks) {
    blocks <- Filter(Negate(is.null), blocks)
    if (length(blocks) == 0L) NULL else
      do.call(cbind, lapply(blocks, function(block) unname(as.matrix(block))))
  }
  fit_ordinal <- lapply(stats::setNames(nm = names(reference$ordinal)), function(item) {
    .ordinal_probabilities(fit$ordinal_intercepts[[item]], fit$ordinal_locations[, item])
  })
  ours <- signature(c(
    list(if (!is.null(reference$means)) fit$means[, colnames(reference$means), drop = FALSE]),
    fit$response_probabilities[names(reference$responses)], fit_ordinal,
    list(if (!is.null(reference$count_means))
      log(fit$count_means[, colnames(reference$count_means), drop = FALSE]))))
  theirs <- signature(c(list(reference$means), reference$responses, reference$ordinal,
                        list(if (!is.null(reference$count_means)) log(reference$count_means))))
  order_profiles <- if (is.null(theirs)) seq_len(fit$n_profiles) else
    align_classes(ours, theirs)
  permute_rows <- function(m) if (is.null(m)) NULL else m[order_profiles, , drop = FALSE]
  reference$means <- permute_rows(reference$means)
  reference$variances <- permute_rows(reference$variances)
  if (!is.null(reference$covariances)) {
    reference$covariances <- reference$covariances[, , order_profiles, drop = FALSE]
  }
  reference$responses <- lapply(reference$responses, permute_rows)
  reference$ordinal <- lapply(reference$ordinal, permute_rows)
  reference$count_means <- permute_rows(reference$count_means)
  reference$count_dispersion <- permute_rows(reference$count_dispersion)
  reference$standard_errors <- lapply(reference$standard_errors, function(block) {
    if (is.list(block)) lapply(block, permute_rows) else permute_rows(block)
  })
  reference$posteriors <- reference$posteriors[, order_profiles, drop = FALSE]
  if (!is.null(reference$coefficients$profile)) {
    reference$coefficients$profile <-
      reference$coefficients$profile[, order_profiles, drop = FALSE]
  }
  order_groups <- NULL
  if (!is.null(reference$profile_probabilities)) {
    reference$profile_probabilities <-
      reference$profile_probabilities[, order_profiles, drop = FALSE]
    if (nrow(reference$profile_probabilities) > 1L) {
      order_groups <- align_classes(unname(fit$profile_probabilities),
                                    unname(reference$profile_probabilities))
      reference$profile_probabilities <-
        reference$profile_probabilities[order_groups, , drop = FALSE]
    }
  } else if (inherits(fit, "multilpa_covariates") &&
               ncol(reference$group_posteriors %||% matrix(0, 0L, 1L)) > 1L &&
               !is.null(fit$group_posteriors)) {
    # Both in the order the groups first appear in the data, as
    # reference_tables() orders them (the fit names its rows by each group's
    # first data row, not by the group label).
    theirs_groups <- reference$group_posteriors
    if (!is.null(reference$groups) && !is.null(rownames(theirs_groups))) {
      theirs_groups <- theirs_groups[as.character(unique(reference$groups)), , drop = FALSE]
    }
    order_groups <- align_classes(t(unname(fit$group_posteriors)), t(unname(theirs_groups)))
  }
  if (!is.null(order_groups)) {
    reference$group_probabilities <- reference$group_probabilities[order_groups]
    if (!is.null(reference$group_posteriors)) {
      reference$group_posteriors <- reference$group_posteriors[, order_groups, drop = FALSE]
    }
    if (!is.null(reference$coefficients$group)) {
      reference$coefficients$group <-
        reference$coefficients$group[, order_groups, drop = FALSE]
    }
    if (!is.null(reference$coefficients$profile)) {
      rownames(reference$coefficients$profile) <- .renumber_group_classes(
        rownames(reference$coefficients$profile), order_groups)
    }
  }
  reference
}

#' Renumber the group classes named inside coefficient terms
#'
#' Under `profile_slopes = "group_class"` a profile logit's terms are
#' `group_class_h` and `x:group_class_h`. When the reference's group class `i`
#' is the fit's class `j` (`order[j] == i`), its terms are renamed to `j`.
#'
#' @param terms Term labels in the reference's numbering.
#' @param order The group-class permutation from [align_classes()].
#' @return The labels in the fit's numbering.
.renumber_group_classes <- function(terms, order) {
  pattern <- "(^|:)group_class_([0-9]+)($|:)"
  vapply(terms, function(term) {
    found <- regmatches(term, regexec(pattern, term))[[1L]]
    if (length(found) == 0L) return(term)
    renumbered <- match(as.integer(found[[3L]]), order)
    sub(pattern, sprintf("\\1group_class_%d\\3", renumbered), term)
  }, character(1), USE.NAMES = FALSE)
}

#' Ordinal category probabilities, from the adjacent-category definition
#'
#' `log P(k | c) = a_k + eta_c (k - 1) - log sum_m exp(a_m + eta_c (m - 1))`,
#' with `a_1 = 0` (`?multilpa`, `ordinal`).
#'
#' @param intercepts The `a_2, ..., a_K`.
#' @param locations One `eta_c` per profile.
#' @return A profiles-by-categories matrix.
.ordinal_probabilities <- function(intercepts, locations) {
  scores <- seq_len(length(intercepts) + 1L) - 1
  logits <- outer(locations, scores) +
    matrix(c(0, intercepts), length(locations), length(scores), byrow = TRUE)
  top <- apply(logits, 1L, max)
  shifted <- exp(logits - top)
  shifted / rowSums(shifted)
}

#' Classification entropy, `-sum p log p`, with `0 log 0 = 0`
#' @param p A posterior matrix.
#' @return A single number.
.reference_entropy <- function(p, weights = NULL) {
  if (is.null(weights)) return(-sum(ifelse(p > 0, p * log(p), 0)))
  -sum(weights * rowSums(ifelse(p > 0, p * log(p), 0)))
}

#' The classification tables of one level, from its posteriors
#'
#' Modal assignment takes the first class on a tie. Odds of correct
#' classification are `Inf` where the average posterior is one to within
#' `sqrt(.Machine$double.eps)`, the package's documented convention, and `NA`
#' where a class has no modal member or a share of zero or one.
#'
#' Under sampling weights (`?get_results`, "classification" to "entropy")
#' every quantity that estimates the population is weighted: class sizes and
#' shares are weighted posterior totals, average posteriors weighted means,
#' the error matrix `sum_i w_i p_it 1(W_i = s) / sum_i w_i p_it`, the BCH
#' weights its inverse, the entropy `sum_i w_i e_i` and the relative entropy
#' `1 - sum_i w_i e_i / (sum_i w_i log K)`; modal counts and `n_units` stay
#' counts.
#'
#' @param p Posteriors of the classified units.
#' @param level `"individuals"` or `"groups"`.
#' @param units Unit labels for the BCH table.
#' @param weights `NULL`, or the rescaled sampling weight of each row of `p`.
#' @return A named list of tables.
.reference_level_tables <- function(p, level, units, weights = NULL) {
  k <- ncol(p)
  classes <- seq_len(k)
  w <- weights %||% rep(1, nrow(p))
  modal <- max.col(p, ties.method = "first")
  n_modal <- tabulate(modal, nbins = k)
  average <- t(vapply(classes, function(assigned) {
    rows <- modal == assigned
    total <- sum(w[rows])
    if (n_modal[[assigned]] == 0L || total <= 0) rep(NA_real_, k) else
      colSums(p[rows, , drop = FALSE] * w[rows]) / total
  }, numeric(k)))
  estimated <- colSums(p * w)
  share <- estimated / sum(w)
  diagonal <- diag(average)
  odds <- ifelse(is.na(diagonal) | share <= 0 | share >= 1, NA_real_,
          ifelse(diagonal >= 1 - sqrt(.Machine$double.eps), Inf,
                 (diagonal / (1 - diagonal)) / (share / (1 - share))))
  errors <- t(vapply(classes, function(true_class) {
    vapply(classes, function(assigned) sum((w * p[, true_class])[modal == assigned]),
           numeric(1)) / sum(w * p[, true_class])
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
                         entropy_sum = .reference_entropy(p, weights),
                         relative_entropy = if (k > 1L)
                           1 - .reference_entropy(p, weights) / (sum(w) * log(k)) else
                             NA_real_))
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
#' model-implied one, Fisher-z statistic. Categorical pairs (categorical or
#' ordinal indicators, an ordinal one's categories being its levels): the
#' posterior-weighted two-way table against the product of the profile's
#' margins, Pearson chi-square, both tables summarised by Cramer's V. Mixed
#' Gaussian/category pairs, and `count` pairs (any pair with a count
#' indicator), use the posterior-weighted Pearson correlation of the two
#' observed values, a category entering as its number, centred on their
#' weighted means. Local independence implies zero correlation; there is no
#' reference test distribution for these pairs. Under sampling weights every
#' posterior is multiplied by its row's weight, and `statistic`, `df`,
#' `p_value` and `p_adjusted` are `NA`: design-based inference is not
#' implemented. Written from the definitions in `?get_results`, not from the
#' package's code.
#'
#' @param reference The reference, for its posteriors, data and Gaussian
#'   estimates.
#' @param continuous,categorical,counts Indicator names of each kind;
#'   `categorical` includes the ordinal indicators.
#' @param levels Named list of the category labels of `categorical`.
#' @param blocks Named list of profiles-by-categories probability matrices of
#'   `categorical`.
#' @param row_weights `NULL`, or the rescaled sampling weight of every row.
#' @return A `data.frame` keyed by `profile`, `indicator_1`, `indicator_2`.
.reference_residuals <- function(reference, continuous, categorical, levels,
                                 blocks = reference$responses, counts = character(),
                                 row_weights = NULL) {
  weighted <- !is.null(row_weights)
  p <- if (weighted) reference$posteriors * row_weights else reference$posteriors
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
  # A category enters a correlation as its number, in category order.
  numeric_values <- function(pair) {
    vapply(pair, function(name) {
      if (name %in% categorical) as.numeric(match(as.character(data[[name]]), levels[[name]])) else
        as.numeric(data[[name]])
    }, numeric(nrow(data)))
  }
  correlation_rows <- function(pairs, kind, k, label) {
    lapply(seq_len(ncol(pairs)), function(index) {
      pair <- pairs[, index]
      values <- numeric_values(pair)
      complete <- stats::complete.cases(values)
      weights <- p[complete, k]
      effective <- sum(weights)
      # cov.wt normalizes the weights and centres on weighted means. The
      # normalization cancels in the correlation; no fitted Gaussian mean
      # is used for the descriptive mixed-pair association.
      covariance <- if (effective > 0) {
        stats::cov.wt(values[complete, , drop = FALSE], wt = weights, method = "ML")$cov
      } else matrix(0, 2L, 2L)
      observed <- if (all(diag(covariance) > 0)) {
        stats::cov2cor(covariance)[1L, 2L]
      } else NA_real_
      data.frame(profile = label, indicator_1 = pair[[1L]], indicator_2 = pair[[2L]],
                 kind = kind, observed = observed, expected = 0,
                 residual = observed, effective_n = effective,
                 statistic = NA_real_, df = NA_integer_, p_value = NA_real_)
    })
  }
  in_data_order <- function(names_wanted) names(data)[names(data) %in% names_wanted]
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
        statistic <- if (weighted) NA_real_ else
          (atanh(clamp(observed)) - atanh(clamp(expected))) * sqrt(effective - 3)
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
        expected <- sum(w) * outer(blocks[[first]][k, levels[[first]]],
                                   blocks[[second]][k, levels[[second]]])
        positive <- expected > 0
        statistic <- if (weighted) NA_real_ else if (any(!positive & observed > 0)) Inf else
          sum((observed[positive] - expected[positive])^2 / expected[positive])
        degrees <- if (weighted) NA_integer_ else
          as.integer((nrow(observed) - 1L) * (ncol(observed) - 1L))
        data.frame(profile = label, indicator_1 = first, indicator_2 = second,
                   kind = "categorical", observed = cramer(observed),
                   expected = cramer(expected),
                   residual = cramer(observed) - cramer(expected), effective_n = sum(w),
                   statistic = statistic, df = degrees,
                   p_value = stats::pchisq(statistic, degrees, lower.tail = FALSE))
      })
    }
    mixed <- if (length(continuous) > 0L && length(categorical) > 0L) {
      pairs <- utils::combn(in_data_order(c(continuous, categorical)), 2L)
      cross_kind <- xor(pairs[1L, ] %in% continuous, pairs[2L, ] %in% continuous)
      correlation_rows(pairs[, cross_kind, drop = FALSE], "mixed", k, label)
    }
    counted <- in_data_order(c(continuous, categorical, counts))
    count <- if (length(counts) > 0L && length(counted) >= 2L) {
      pairs <- utils::combn(counted, 2L)
      with_count <- pairs[1L, ] %in% counts | pairs[2L, ] %in% counts
      correlation_rows(pairs[, with_count, drop = FALSE], "count", k, label)
    }
    pieces <- c(gaussian, nominal, mixed, count)
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
#' and whether that eigenvalue sits at the variance floor. Under sampling
#' weights (`?get_results`, `gaussian_diagnostics`) each row's contribution is
#' its posterior times its rescaled weight, the marginal variance is the
#' sampling-weighted one (divisor the weight total), and a row of weight zero
#' is not counted.
#'
#' `positive_rows` is the exact value: a Gaussian density is positive
#' everywhere, so every row with a Gaussian value has a positive posterior in
#' every component. Counting computed posteriors above zero instead counts the
#' rows whose posterior did not underflow, which differs between programs (on
#' thyroid VVV the package keeps a posterior of 1.2e-320 that mclust rounds to
#' zero). Disagreement here is a finding about that column, not tolerance to
#' tune.
#'
#' @param row_weights `NULL`, or the rescaled sampling weight of every row.
#' @return A `data.frame` keyed by `context` and `profile`.
.reference_gaussian_diagnostics <- function(reference, continuous, min_variance,
                                            row_weights = NULL) {
  data <- as.matrix(reference$data[continuous])
  p <- reference$posteriors
  with_value <- rowSums(!is.na(data)) > 0L
  overall <- if (is.null(row_weights)) {
    apply(data, 2L, function(column) {
      column <- column[!is.na(column)]
      mean((column - mean(column))^2)
    })
  } else {
    with_value <- with_value & row_weights > 0
    vapply(seq_len(ncol(data)), function(j) {
      used <- !is.na(data[, j]) & row_weights > 0
      w <- row_weights[used]
      centre <- sum(w * data[used, j]) / sum(w)
      sum(w * (data[used, j] - centre)^2) / sum(w)
    }, numeric(1))
  }
  do.call(rbind, lapply(seq_len(ncol(p)), function(k) {
    contribution <- p[, k] * with_value * (row_weights %||% 1)
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
#' Beyond the measurement tables, a reference can supply the `ordinal` table
#' (category probabilities; each profile's location derived from them by the
#' adjacent-category definition, `eta_c = (d_K - d_1) / (K - 1)` with
#' `d_k = log P(k | c) - log P(k | last profile)`), the `count_means` table,
#' and for a membership-covariate fit the `coefficients` table in place of
#' `profile_probabilities`. A weighted fit takes the reference's weights
#' through the package's documented convention (see [harness_reference()]).
#'
#' @param reference An aligned `harness_reference`.
#' @param fit The fitted model, for its dimensions, labels and settings only;
#'   no estimate is read from it.
#' @return A named list of tables, with an `exempt` attribute naming (with
#'   reasons) the outputs left out because the reference cannot supply them.
reference_tables <- function(reference, fit) {
  stopifnot("`reference` must come from harness_reference()" =
              inherits(reference, "harness_reference"))
  weighted <- !is.null(fit$sampling_weights)
  if (!identical(weighted, !is.null(reference$weights))) {
    stop(errorCondition(paste(
      "a weighted fit needs the reference's sampling weights, and an unweighted",
      "fit must not be given any"), class = "multilpa_bad_reference", call = NULL))
  }
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
  exempt <- character()

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
  tables$ordinal <- .reference_ordinal(reference, k)
  tables$count_means <- .reference_counts(reference, k)

  group_of_row <- reference$groups %||% rows
  group_ids <- unique(group_of_row)
  zg <- reference$group_posteriors %||%
    matrix(1, length(group_ids), 1L, dimnames = list(as.character(group_ids), NULL))
  zg <- zg[as.character(group_ids), , drop = FALSE]
  # Sampling weights, rescaled as the package rescales them; a unit of weight
  # zero keeps its posteriors but is in no count, entropy or sample size.
  unit_weights <- if (weighted) .reference_unit_weights(reference$weights, group_ids)
  row_weights <- if (weighted) unname(unit_weights[as.character(group_of_row[rows])])
  classified_row <- informative & (if (weighted) row_weights > 0 else TRUE)
  informative_group <- as.character(group_ids) %in%
    as.character(unique(group_of_row[rows[classified_row]]))
  if (weighted) informative_group <- informative_group & unit_weights > 0

  if (inherits(fit, "multilpa_covariates")) {
    # Class membership is a regression here, so there is no profile mix.
    n_group_classes <- ncol(zg)
    if (is.null(reference$coefficients)) {
      exempt[["coefficients"]] <- "the reference supplies no membership coefficients"
    } else {
      tables$coefficients <- .reference_coefficients(reference$coefficients)
    }
  } else {
    gamma <- reference$profile_probabilities %||% matrix(colMeans(p), 1L)
    n_group_classes <- nrow(gamma)
    delta <- reference$group_probabilities %||% rep(1, n_group_classes)
    tables$profile_probabilities <- data.frame(
      group_class = rep(seq_len(n_group_classes), times = k),
      profile = rep(seq_len(k), each = n_group_classes),
      probability = as.vector(gamma), group_class_probability = rep(delta, times = k))
  }

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

  individual <- .reference_level_tables(p[classified_row, , drop = FALSE], "individuals",
                                        rows[classified_row], row_weights[classified_row])
  group <- .reference_level_tables(zg[informative_group, , drop = FALSE], "groups",
                                   as.character(group_ids)[informative_group],
                                   unname(unit_weights[informative_group]))
  # One group class classifies nothing at the group level: the fit leaves its
  # group posteriors and the group rows of the classification tables out of
  # `get_results(x, "all")`, keeping the group rows of `counts` and `entropy`
  # (`?get_results`; package review API-9).
  group_level <- if (ncol(zg) > 1L) names(group) else c("counts", "entropy")
  if (ncol(zg) == 1L) tables$group_posteriors <- NULL
  tables[names(individual)] <- lapply(names(individual), function(name) {
    rbind(individual[[name]], if (name %in% group_level) group[[name]])
  })

  n_individuals <- sum(classified_row)
  n_groups <- sum(informative_group)
  # ICL, AWE and CLC charge the entropy of the complete-data classification,
  # H(group classes) + sum_i E_C[H(X_i | C, y_i)] over the classified groups
  # and rows, under both conventions; with one group class that is the row
  # posterior entropy. A transition model keeps each level's own entropy.
  # Without the rows' within-class entropies a model with several group
  # classes has no reference for the three, and they are exempted by name, so
  # they can never agree by default. Under sampling weights each unit's term
  # is multiplied by its weight.
  row_entropy <- function() {
    .reference_entropy(p[classified_row, , drop = FALSE], row_weights[classified_row])
  }
  group_entropy <- function() {
    .reference_entropy(zg[informative_group, , drop = FALSE],
                       unname(unit_weights[informative_group]))
  }
  entropies <- if (inherits(fit, "multilpa_transitions")) {
    c(individuals = row_entropy(), groups = group_entropy())
  } else {
    complete <- if (ncol(zg) == 1L) {
      row_entropy()
    } else if (!is.null(reference$row_conditional_entropy)) {
      group_entropy() + sum((row_weights %||% 1) * reference$row_conditional_entropy *
                              classified_row)
    } else {
      charged <- paste0("information_criteria$",
                        c(outer(c("awe", "icl", "clc"), c("groups", "individual"),
                                paste, sep = "_")))
      exempt <- c(exempt, stats::setNames(rep(paste(
        "the reference supplies no row entropies within the group classes",
        "(`row_conditional_entropy`), which the complete-data entropy of a",
        "model with several group classes needs"), length(charged)), charged))
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
  if (!is.null(tables$coefficients)) {
    tables$model$n_profile_covariates <- .reference_covariate_count(reference$coefficients$profile)
    tables$model$n_group_covariates <- .reference_covariate_count(reference$coefficients$group)
  }

  if (!is.null(data)) {
    # An ordinal indicator pairs as a categorical one, its levels its
    # categories. The package pairs the categorical indicators first and the
    # ordinal ones after them, each in the fit's order
    # (.multilpa_discrete_items()), whatever their order in `vars`.
    ordinal <- names(reference$ordinal)
    paired <- c(categorical, ordinal[order(match(ordinal, fit$ordinal %||% ordinal))])
    tables$residuals <- .reference_residuals(
      reference, continuous, paired,
      c(fit$categorical_levels, lapply(reference$ordinal, colnames)),
      blocks = c(reference$responses, reference$ordinal),
      counts = colnames(reference$count_means) %||% character(), row_weights = row_weights)
    if (length(continuous) > 0L) {
      tables$gaussian_diagnostics <- .reference_gaussian_diagnostics(
        reference, continuous, fit$min_variance %||% 1e-6, row_weights)
    }
  }
  # Outputs this reference cannot speak to are left out here, and named in the
  # `exempt` attribute, which `compare_fit()` adds to the suite's exemptions
  # so that the ledger gives the reason.
  tables <- .drop_exempt(Filter(Negate(is.null), tables), exempt)
  attr(tables, "exempt") <- exempt
  tables
}

#' One rescaled sampling weight per unit
#'
#' The package rescales the weights to sum to the number of units of positive
#' weight (`?multilpa`, `weights`).
#'
#' @param weights Raw weights, named by unit or in unit order.
#' @param group_ids The units, in the reference's order.
#' @return The rescaled weights, named by unit.
.reference_unit_weights <- function(weights, group_ids) {
  labels <- as.character(group_ids)
  if (!is.null(names(weights))) weights <- weights[labels]
  if (length(weights) != length(labels) || anyNA(weights)) {
    stop(errorCondition("the reference's weights must cover every unit once",
                        class = "multilpa_bad_reference", call = NULL))
  }
  stats::setNames(unname(weights) * sum(weights > 0) / sum(weights), labels)
}

#' The membership coefficient table, against the fit's baseline
#'
#' Multinomial logits are identified up to a constant per term, so each
#' block is re-expressed against its last class (the package's baseline) by
#' subtracting that column, which is then dropped.
#'
#' @param coefficients The aligned `coefficients` of a `harness_reference`.
#' @return A `data.frame` with `level`, `outcome`, `term`, `parameter` and
#'   `estimate`.
.reference_coefficients <- function(coefficients) {
  prefixes <- c(profile = "profile", group = "group_class")
  blocks <- Filter(Negate(is.null), coefficients[intersect(names(prefixes), names(coefficients))])
  do.call(rbind, lapply(names(blocks), function(level) {
    block <- blocks[[level]]
    against_last <- block[, -ncol(block), drop = FALSE] - block[, ncol(block)]
    data.frame(level = level,
               outcome = rep(sprintf("%s_%d", prefixes[[level]], seq_len(ncol(against_last))),
                             each = nrow(against_last)),
               term = rep(rownames(block), times = ncol(against_last)),
               parameter = "logit", estimate = as.vector(against_last),
               stringsAsFactors = FALSE)
  }))
}

#' How many covariates a coefficient block's terms name
#' @param block Terms-by-classes matrix, or `NULL`.
#' @return A whole number: the distinct covariates, intercepts (and the
#'   group-class intercepts of `profile_slopes = "group_class"`) not counted.
.reference_covariate_count <- function(block) {
  if (is.null(block)) return(0L)
  covariates <- sub("(^|:)group_class_[0-9]+$", "", rownames(block))
  length(setdiff(unique(covariates), c("(Intercept)", "")))
}

#' The ordinal table, from the reference's category probabilities
#' @param reference An aligned `harness_reference`.
#' @param k Number of profiles.
#' @return A `data.frame`, or `NULL` when the reference has no ordinal block.
.reference_ordinal <- function(reference, k) {
  if (length(reference$ordinal) == 0L) return(NULL)
  do.call(rbind, lapply(names(reference$ordinal), function(item) {
    block <- reference$ordinal[[item]]
    categories <- ncol(block)
    # log P(k | c) - log P(k | last) = eta_c (k - 1) - (normalisers), so the
    # location is the slope across categories; exactly zero for the last.
    against_last <- log(block) - matrix(log(block[k, ]), k, categories, byrow = TRUE)
    location <- (against_last[, categories] - against_last[, 1L]) / (categories - 1L)
    out <- data.frame(profile = rep(seq_len(k), times = categories), indicator = item,
                      category = rep(colnames(block), each = k),
                      probability = as.vector(block),
                      location = rep(unname(location), times = categories),
                      stringsAsFactors = FALSE)
    se <- reference$standard_errors$location
    if (!is.null(se)) out$location_standard_error <- rep(unname(se[, item]), times = categories)
    out
  }))
}

#' The count table, from the reference's means and dispersions
#' @param reference An aligned `harness_reference`.
#' @param k Number of profiles.
#' @return A `data.frame`, or `NULL` when the reference has no counts.
.reference_counts <- function(reference, k) {
  means <- reference$count_means
  if (is.null(means)) return(NULL)
  items <- colnames(means)
  column <- function(m) as.vector(m[, items, drop = FALSE])
  out <- data.frame(profile = rep(seq_len(k), times = length(items)),
                    indicator = rep(items, each = k), mean = column(means),
                    stringsAsFactors = FALSE)
  if (!is.null(reference$count_dispersion)) out$dispersion <- column(reference$count_dispersion)
  se <- reference$standard_errors
  if (!is.null(se$count_mean)) out$standard_error <- column(se$count_mean)
  if (!is.null(se$count_dispersion)) out$dispersion_standard_error <- column(se$count_dispersion)
  out
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
