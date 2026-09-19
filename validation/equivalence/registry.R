# Equivalence harness for multilpa.
#
# Every external comparison in this package reports through one tidy frame, so
# that the evidence can be regenerated and read as a single table rather than
# recovered by reading a pile of scripts. A suite is a file in suites/ that
# defines `suite_<name>()` returning the frame built by `compare_values()`.
#
# Run every suite from the project root:
#   Rscript validation/equivalence/run.R

#' Compare reference values against values this package produced
#'
#' @param quantity Character vector naming each compared quantity.
#' @param reference Numeric vector of reference values from the external source.
#' @param obtained Numeric vector of values multilpa produced.
#' @param tolerance Numeric tolerance, recycled across quantities.
#' @param scale `"absolute"` compares `|reference - obtained|`; `"relative"`
#'   divides that difference by `max(|reference|, 1)`, for quantities whose
#'   magnitude makes an absolute tolerance meaningless.
#' @return A `data.frame` with one row per quantity and columns `quantity`,
#'   `reference`, `obtained`, `difference`, `tolerance`, `agrees`.
compare_values <- function(quantity, reference, obtained, tolerance = 1e-8,
                           scale = c("absolute", "relative")) {
  scale <- match.arg(scale)
  stopifnot(
    "`quantity` must be character" = is.character(quantity),
    "`reference` must be numeric" = is.numeric(reference),
    "`obtained` must be numeric" = is.numeric(obtained),
    "`reference` and `obtained` must have one value per quantity" =
      length(reference) == length(quantity) && length(obtained) == length(quantity),
    "`tolerance` must be positive" = is.numeric(tolerance) && all(tolerance > 0)
  )
  divisor <- if (identical(scale, "relative")) pmax(abs(reference), 1) else 1
  difference <- abs(reference - obtained) / divisor
  data.frame(quantity = quantity, reference = reference, obtained = obtained,
             difference = difference, tolerance = tolerance,
             agrees = difference <= tolerance, stringsAsFactors = FALSE)
}

#' Compare against a value read from published digits
#'
#' A value printed to `digits` decimals states its quantity only to within half
#' a unit in the last place. Reporting the raw gap against such a reference
#' measures how the source was typeset, not how the implementation behaves: a
#' G-squared printed as `15.3` sits up to 0.05 from whatever was computed, so an
#' exact reproduction still shows a difference of that order.
#'
#' This expresses the gap as a fraction of the rounding the printed value
#' admits. A value of 0.77 means the computed quantity is well inside the
#' printed precision and therefore reproduces the published digits; anything
#' above 1 means it does not, and no amount of rounding explains it.
#'
#' @param quantity Character vector naming each compared quantity.
#' @param reference Numeric vector of published values.
#' @param obtained Numeric vector of values multilpa produced.
#' @param digits Number of decimals each reference was printed to, recycled
#'   across quantities.
#' @return A `data.frame` shaped like [compare_values()], whose `difference` is
#'   in units of the half-unit rounding bound and whose `tolerance` is 1.
compare_printed <- function(quantity, reference, obtained, digits) {
  stopifnot(
    "`digits` must be a nonnegative whole number" =
      is.numeric(digits) && all(digits >= 0) && all(digits == floor(digits)),
    "`reference` and `obtained` must have one value per quantity" =
      length(reference) == length(quantity) && length(obtained) == length(quantity))
  bound <- 0.5 * 10^(-digits)
  data.frame(quantity = quantity, reference = reference, obtained = obtained,
             difference = abs(reference - obtained) / bound,
             tolerance = 1, agrees = abs(reference - obtained) <= bound,
             stringsAsFactors = FALSE)
}

#' Load and run every registered equivalence suite
#'
#' @param suites Character vector of suite names, or `NULL` for every suite
#'   found in `validation/equivalence/suites`.
#' @param root Project root directory.
#' @return A `data.frame` with one row per compared quantity and columns
#'   `suite`, `source`, `precision`, `dataset`, `quantity`, `reference`,
#'   `obtained`, `difference`, `tolerance`, `agrees`, `seconds`. `precision` is
#'   `"machine"` when the reference was computed in this session at full double
#'   precision and `"printed"` when it was read from published digits, which
#'   bounds how closely any implementation can agree with it.
run_equivalence <- function(suites = NULL, root = ".") {
  suite_dir <- file.path(root, "validation", "equivalence", "suites")
  files <- sort(list.files(suite_dir, pattern = "[.]R$", full.names = TRUE))
  names(files) <- sub("[.]R$", "", basename(files))
  if (!is.null(suites)) {
    unknown <- setdiff(suites, names(files))
    if (length(unknown) > 0L) {
      stop(sprintf("Unknown suite(s): %s", paste(unknown, collapse = ", ")))
    }
    files <- files[suites]
  }
  rows <- lapply(names(files), function(name) {
    environment_for_suite <- new.env(parent = globalenv())
    sys.source(files[[name]], envir = environment_for_suite)
    entry <- get(paste0("suite_", gsub("-", "_", name)), envir = environment_for_suite)
    message(sprintf("running suite: %s", name))
    started <- proc.time()[["elapsed"]]
    result <- entry()
    elapsed <- proc.time()[["elapsed"]] - started
    stopifnot("a suite must return a data.frame with a `source` column" =
                is.data.frame(result) && "source" %in% names(result))
    if (is.null(result$dataset)) result$dataset <- NA_character_
    if (is.null(result$precision)) result$precision <- "machine"
    data.frame(suite = name, source = result$source,
               precision = result$precision, dataset = result$dataset,
               result[, c("quantity", "reference", "obtained", "difference",
                          "tolerance", "agrees")],
               seconds = round(elapsed, 2), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Summarise an equivalence report by source
#'
#' @param report The frame returned by [run_equivalence()].
#' @return A `data.frame` with one row per suite and source, giving the number
#'   of quantities compared, how many agreed, and the largest difference seen.
equivalence_summary <- function(report) {
  stopifnot("`report` must be an equivalence report" =
              is.data.frame(report) && all(c("suite", "source", "agrees") %in% names(report)))
  parts <- split(report, list(report$suite, report$source, report$precision),
                 drop = TRUE)
  summaries <- lapply(parts, function(part) {
    data.frame(suite = part$suite[[1L]], source = part$source[[1L]],
               precision = part$precision[[1L]],
               compared = nrow(part), agreed = sum(part$agrees),
               worst_difference = max(part$difference),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, summaries)
  row.names(out) <- NULL
  out[order(out$precision, out$suite, out$source), ]
}

#' Stack a per-indicator list of class-by-category matrices into one matrix
#'
#' @param blocks List of matrices sharing a row count, one per indicator.
#' @return A matrix with one row per latent class and every category in columns.
flatten_blocks <- function(blocks) {
  stopifnot("`blocks` must be a list of matrices" =
              is.list(blocks) && all(vapply(blocks, is.matrix, logical(1))))
  do.call(cbind, lapply(blocks, function(block) unname(as.matrix(block))))
}

#' Align two labelings of the same latent classes
#'
#' Latent class labels are arbitrary, so two programs fitting the same model
#' recover the same classes under some permutation. This finds the permutation
#' minimising total absolute distance, which is exact for the class counts used
#' here rather than a greedy match.
#'
#' @param target Matrix with one row per class, in the order to match.
#' @param candidate Matrix of the same dimensions in an unknown class order.
#' @return An integer permutation `p` such that `candidate[p, ]` aligns with
#'   `target`.
align_classes <- function(target, candidate) {
  target <- unname(as.matrix(target))
  candidate <- unname(as.matrix(candidate))
  stopifnot("`target` and `candidate` must have the same dimensions" =
              identical(dim(target), dim(candidate)))
  orders <- .class_permutations(nrow(target))
  distances <- vapply(seq_len(nrow(orders)), function(row) {
    sum(abs(target - candidate[orders[row, ], , drop = FALSE]))
  }, numeric(1))
  orders[which.min(distances), ]
}

#' Every permutation of seq_len(n), one per row
#' @param n Number of elements; kept small because this enumerates n! rows.
#' @return An integer matrix with `factorial(n)` rows.
.class_permutations <- function(n) {
  stopifnot("permutation search is only used for small class counts" = n <= 8L)
  if (n == 1L) return(matrix(1L, 1L, 1L))
  smaller <- .class_permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(value) {
    cbind(value, matrix(setdiff(seq_len(n), value)[smaller], nrow(smaller)),
          deparse.level = 0)
  }))
}

#' Render an equivalence report as markdown
#'
#' @param report The frame returned by [run_equivalence()].
#' @param summary_table The frame returned by [equivalence_summary()].
#' @return A character vector of markdown lines.
.equivalence_markdown <- function(report, summary_table) {
  row_text <- function(values) paste0("| ", paste(values, collapse = " | "), " |")
  header <- c(
    "# Equivalence report",
    "",
    sprintf("Generated by `Rscript validation/equivalence/run.R` under %s on %s.",
            R.version.string, format(Sys.Date())),
    "",
    sprintf("**%d quantities compared across %d sources; %d agreed within tolerance.**",
            nrow(report), length(unique(report$source)), sum(report$agrees)),
    "",
    sprintf("Largest difference against a reference computed here at full double precision: **%.3e**.",
            max(c(report$difference[report$precision == "machine"], 0))),
    sprintf("Largest gap against a reference read from published digits: **%.3g** of the rounding that value's printed precision admits, where 1 would be the whole of it. Published references are compared this way because the raw gap against a statistic printed to one decimal measures the typesetting, not the implementation.",
            max(c(report$difference[report$precision == "printed"], 0))),
    "",
    "## By source",
    "",
    row_text(c("Suite", "Reference source", "Reference precision", "Dataset(s)",
               "Compared", "Agreed", "Largest difference")),
    row_text(rep("---", 7L)))
  datasets <- vapply(seq_len(nrow(summary_table)), function(row) {
    part <- report[report$suite == summary_table$suite[[row]] &
                     report$source == summary_table$source[[row]] &
                     report$precision == summary_table$precision[[row]], , drop = FALSE]
    paste(sort(unique(stats::na.omit(part$dataset))), collapse = ", ")
  }, character(1))
  body <- vapply(seq_len(nrow(summary_table)), function(row) {
    row_text(c(summary_table$suite[[row]], summary_table$source[[row]],
               summary_table$precision[[row]], datasets[[row]],
               summary_table$compared[[row]], summary_table$agreed[[row]],
               sprintf("%.2e", summary_table$worst_difference[[row]])))
  }, character(1))
  failures <- report[!report$agrees, , drop = FALSE]
  tail_lines <- if (nrow(failures) == 0L) {
    c("", "Every compared quantity agreed within its stated tolerance.")
  } else {
    c("", "## Disagreements", "",
      row_text(c("Suite", "Quantity", "Reference", "Obtained", "Difference",
                 "Tolerance")),
      row_text(rep("---", 6L)),
      vapply(seq_len(nrow(failures)), function(row) {
        row_text(c(failures$suite[[row]], failures$quantity[[row]],
                   sprintf("%.6g", failures$reference[[row]]),
                   sprintf("%.6g", failures$obtained[[row]]),
                   sprintf("%.2e", failures$difference[[row]]),
                   sprintf("%.2e", failures$tolerance[[row]])))
      }, character(1)))
  }
  c(header, body, tail_lines, "",
    "Per-quantity detail, including every tolerance, is in `report.csv`.")
}
