# Reading Latent GOLD output and comparing it with this package's targets.
#
# Sourced by compare.R and by test-compare.R. Nothing here runs a model.
#
# What is compared, and how robustly it can be read:
#
# - The log likelihood and the parameter count, from the listing. Both are
#   invariant to how either program parameterises the model.
# - The posterior class probabilities, from the posterior file each model
#   writes. Also parameterisation-invariant, and the strongest single check:
#   two fits that agree on every posterior agree on the model.
# - Profile means and variances, recomputed here from Latent GOLD's own
#   posteriors. At convergence EM is at a fixed point, so its estimates are the
#   posterior-weighted means and variances; this recovers Latent GOLD's values
#   without scraping the parameter tables out of the listing.
# - Group-class proportions, as the mean of Latent GOLD's group posteriors.
#
# Standard errors, bivariate residuals and the Step-3 results exist only in the
# listing, whose table layout is not known until a real listing exists. Those
# rows are reported as `awaiting parser`, with this package's value beside
# them, rather than parsed by a guess that could misread a column.

#' The tolerances this comparison uses
#'
#' Every agreement is judged at the rounding Latent GOLD's printed value
#' admits plus an optimiser allowance, never at a bare number chosen after
#' seeing the result.
#' - `likelihood`: 1e-5, two independent searches stopped at their own
#'   convergence criteria (the equivalence harness uses the same value).
#' - `parameter`: 1e-4, for a posterior or a parameter from two optimisers.
#' @format A named numeric vector.
.lg_tolerances <- c(likelihood = 1e-5, parameter = 1e-4)

#' The latent variable names the kit's syntax uses, by level
#' @format A named character vector.
.lg_latent_names <- c(individual = "Cluster", group = "GClass", state = "State")

#' Half a unit in the last printed place of a set of numeric strings
#'
#' The file's precision is the most decimals any value carries. A writer that
#' drops trailing zeros prints `1.0000` as `1`; taking the fewest decimals
#' would read that as a file printed to whole numbers and widen the tolerance
#' to 0.5, large enough to call any posterior agreement.
#'
#' @param text Character vector of numbers as printed.
#' @return A single number: half a unit in the file's last printed place.
.lg_half_unit <- function(text) {
  text <- trimws(text[!is.na(text) & nzchar(trimws(text))])
  decimals <- vapply(strsplit(sub("[eE].*$", "", text), ".", fixed = TRUE),
                     function(parts) if (length(parts) > 1L) nchar(parts[[2L]]) else 0L,
                     integer(1))
  if (length(decimals) == 0L) return(NA_real_)
  0.5 * 10^(-max(decimals))
}

#' Read a delimited Latent GOLD text file whatever its separator
#'
#' The posterior file's separator depends on the build and the file extension,
#' so it is detected from the header rather than assumed. Every column is read
#' as text first, so that the printed precision is known before conversion.
#'
#' @param path File path.
#' @return A data frame of character columns with the original header names,
#'   or `NULL` when the file cannot be read as a table.
.lg_read_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2L) return(NULL)
  separator <- if (grepl("\t", lines[[1L]], fixed = TRUE)) "\t" else
    if (grepl(",", lines[[1L]], fixed = TRUE)) "," else ""
  table <- tryCatch(
    utils::read.table(text = lines, header = TRUE, sep = separator,
                      colClasses = "character", check.names = FALSE,
                      na.strings = character(), strip.white = TRUE,
                      comment.char = "", quote = "\""),
    # A malformed file is reported as unparsed by the caller, which is the
    # visible outcome this handler exists to produce; other conditions pass.
    error = function(condition) NULL)
  if (is.null(table) || ncol(table) < 2L) NULL else table
}

#' The posterior columns of one latent variable, as a numeric matrix
#'
#' @param table From [.lg_read_table()].
#' @param latent The latent variable's name, e.g. `"Cluster"`.
#' @return A list with `values` (rows as in `table`, one column per class in
#'   class order) and `half_unit`, or `NULL` if the columns are absent.
.lg_posterior_matrix <- function(table, latent) {
  pattern <- sprintf("^%s#([0-9]+)$", latent)
  columns <- grep(pattern, names(table), value = TRUE)
  if (length(columns) < 2L) return(NULL)
  columns <- columns[order(as.integer(sub(pattern, "\\1", columns)))]
  text <- unlist(table[columns], use.names = FALSE)
  values <- vapply(table[columns], function(column) as.numeric(gsub(",", ".", column)),
                   numeric(nrow(table)))
  if (anyNA(values)) return(NULL)
  list(values = matrix(values, nrow(table)), half_unit = .lg_half_unit(text))
}

#' Every ordering of `seq_len(k)`, one per row
#' @param k Number of classes.
#' @return An integer matrix with `factorial(k)` rows.
.lg_permutations <- function(k) {
  if (k <= 1L) return(matrix(seq_len(k), 1L))
  smaller <- .lg_permutations(k - 1L)
  do.call(rbind, lapply(seq_len(k), function(first) {
    rest <- matrix(setdiff(seq_len(k), first)[smaller], nrow = nrow(smaller))
    cbind(first, rest, deparse.level = 0L)
  }))
}

#' Match Latent GOLD's class labels to this package's
#'
#' Labels are arbitrary in both programs. The ordering of Latent GOLD's columns
#' that brings its posteriors closest to this package's is the alignment.
#'
#' @param target This package's posterior matrix.
#' @param reference Latent GOLD's, rows in the same order.
#' @return A list with `order` (Latent GOLD column for each package class) and
#'   `difference` (the largest absolute posterior disagreement under it).
.lg_align <- function(target, reference) {
  stopifnot("posterior matrices must have the same shape" =
              identical(dim(target), dim(reference)))
  orders <- .lg_permutations(ncol(target))
  gaps <- vapply(seq_len(nrow(orders)), function(index) {
    max(abs(target - reference[, orders[index, ], drop = FALSE]))
  }, numeric(1))
  best <- which.min(gaps)
  list(order = orders[best, ], difference = gaps[[best]])
}

#' The lines of a Latent GOLD listing
#'
#' The listing is Latin-1, not UTF-8: it prints `R\u00b2` as a single byte, and
#' reading it as UTF-8 yields invalid strings that no pattern matches. It is
#' tab-separated, with each label in the first field.
#'
#' @param path Path to a `.lst` file.
#' @return Character vector of lines, converted to UTF-8.
.lg_listing_lines <- function(path) {
  iconv(readLines(path, warn = FALSE), from = "latin1", to = "UTF-8", sub = "?")
}

#' One labelled number from a Latent GOLD listing
#'
#' @param lines The listing's lines.
#' @param label A regular expression matching the start of the line's label.
#' @return A list with `value` and `half_unit`, or `NULL` when absent.
.lg_listing_value <- function(lines, label) {
  hit <- grep(label, lines, value = TRUE)
  if (length(hit) == 0L) return(NULL)
  tokens <- regmatches(hit[[1L]], gregexpr("-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?",
                                           hit[[1L]]))[[1L]]
  if (length(tokens) == 0L) return(NULL)
  token <- tokens[[length(tokens)]]
  list(value = as.numeric(token), half_unit = .lg_half_unit(token))
}

#' A block of comparison rows
#'
#' @param case,quantity Labels.
#' @param reference Latent GOLD's value, `NA` when not available.
#' @param obtained This package's value.
#' @param tolerance The agreement bound, `NA` when nothing was compared.
#' @param status `NULL` to judge from the numbers, or a fixed status.
#' @param note Free text.
#' @return A data frame, one row per quantity.
.lg_rows <- function(case, quantity, reference, obtained, tolerance,
                     status = NULL, note = "") {
  difference <- abs(reference - obtained)
  judged <- ifelse(is.na(difference), "unparsed",
                   ifelse(difference <= tolerance, "agree", "disagree"))
  data.frame(case = case, quantity = quantity, reference = reference,
             obtained = obtained, difference = difference, tolerance = tolerance,
             status = status %||% judged, note = note, stringsAsFactors = FALSE)
}

#' Rows for quantities that exist only in the listing's tables
#'
#' @param target One case's targets.
#' @return Rows with status `awaiting parser`, or `NULL`.
.lg_pending_rows <- function(target) {
  case <- target$case
  parameters <- target$parameters
  standard_errors <- if (!is.null(parameters)) {
    .lg_rows(case, sprintf("standard error: %s %s %s %s", parameters$level,
                           parameters$outcome, parameters$term, parameters$parameter),
             NA_real_, parameters$standard_error, NA_real_,
             status = "awaiting parser",
             note = "Latent GOLD prints standard errors in the parameter tables of the listing.")
  }
  residuals <- if (!is.null(target$residuals)) {
    .lg_rows(case, sprintf("bivariate residual: %s %s-%s", target$residuals$profile,
                           target$residuals$indicator_1, target$residuals$indicator_2),
             NA_real_, target$residuals$statistic, NA_real_, status = "awaiting parser",
             note = "Different statistics; compare the ranking of pairs, not the values.")
  }
  steps <- if (!is.null(target$steps)) {
    means <- lapply(c("distal_bch", "distal_modal", "distal_proportional"), function(name) {
      table <- target$steps[[name]]
      .lg_rows(case, sprintf("%s: class %d mean", name, table$class), NA_real_,
               table$estimate, NA_real_, status = "awaiting parser",
               note = sprintf("Listing c09_%s.lst.", name))
    })
    slope <- target$steps$covariate_ml
    slope <- slope[slope$term != "(Intercept)", , drop = FALSE]
    c(means, list(.lg_rows(case, sprintf("covariate_ml: %s slope on %s", slope$outcome, slope$term),
                           NA_real_, slope$estimate, NA_real_, status = "awaiting parser",
                           note = "Log odds against the last class; convert Latent GOLD's coding before comparing.")))
  }
  do.call(rbind, c(list(standard_errors, residuals), steps))
}

#' Profile means and variances recovered from Latent GOLD's posteriors
#'
#' Each estimate is a posterior-weighted moment. Latent GOLD prints each
#' posterior to within `half_unit`, and that rounding perturbs a weighted
#' mean by at most `half_unit * sum(|y - m|) / (sum(w) - n * half_unit)`; the
#' tolerance is that bound plus the optimiser allowance, so it is derived from
#' the file rather than chosen.
#'
#' @param target One case's targets.
#' @param weights Latent GOLD's aligned individual posteriors.
#' @param half_unit Their printed rounding.
#' @return Comparison rows.
.lg_moment_rows <- function(target, weights, half_unit) {
  data <- target$data[match(target$posteriors$id, target$data$id), , drop = FALSE]
  profiles <- target$profiles
  n <- nrow(weights)
  denominators <- colSums(weights) - n * half_unit
  rows <- lapply(target$continuous, function(indicator) {
    y <- data[[indicator]]
    means <- colSums(weights * y) / colSums(weights)
    squared <- outer(y, means, "-")^2
    mean_bounds <- half_unit * colSums(abs(outer(y, means, "-"))) / denominators
    variances <- if (identical(target$variance_model, "equal")) {
      rep(sum(weights * squared) / n, ncol(weights))
    } else colSums(weights * squared) / colSums(weights)
    # A variance's derivative in its mean is zero at the weighted mean, so an
    # error in the mean enters only at second order, as its square.
    variance_bounds <- if (identical(target$variance_model, "equal")) {
      rep(half_unit * sum(squared) / n + sum(mean_bounds^2 * colSums(weights)) / n,
          ncol(weights))
    } else {
      half_unit * colSums(abs(sweep(squared, 2L, variances))) / denominators +
        mean_bounds^2
    }
    obtained <- profiles[profiles$indicator == indicator, , drop = FALSE]
    obtained <- obtained[order(obtained$profile), , drop = FALSE]
    rbind(
      .lg_rows(target$case, sprintf("mean: profile %d %s", obtained$profile, indicator),
               means, obtained$mean, mean_bounds + .lg_tolerances[["parameter"]],
               note = "Recomputed from Latent GOLD's posteriors."),
      .lg_rows(target$case, sprintf("variance: profile %d %s", obtained$profile, indicator),
               variances, obtained$variance,
               variance_bounds + .lg_tolerances[["parameter"]],
               note = "Recomputed from Latent GOLD's posteriors."))
  })
  do.call(rbind, rows)
}

#' Compare every quantity of one case
#'
#' @param target One case's targets, as written by make-kit.R.
#' @param returned Directory holding the files Latent GOLD wrote.
#' @return A data frame of comparison rows, one per quantity.
lg_compare_case <- function(target, returned) {
  case <- target$case
  stem <- if (identical(case, "c09_three_step")) "c09_step1" else case
  listing_path <- file.path(returned, sprintf("%s.lst", stem))
  posterior_path <- file.path(returned, sprintf("%s_posteriors.txt", stem))

  scalars <- if (!file.exists(listing_path)) {
    .lg_rows(case, c("log likelihood", "parameters"), NA_real_,
             c(target$log_likelihood, target$n_parameters), NA_real_,
             status = "missing", note = sprintf("No %s.", basename(listing_path)))
  } else {
    lines <- .lg_listing_lines(listing_path)
    # The labels are as Latent GOLD 6.1 writes them, confirmed against a real
    # listing: "Log-likelihood (LL)" and "Number of parameters (Npar)", each
    # followed by a tab and the value. "Npar" alone appears only as a column
    # heading in the summary table, with no value on that line.
    likelihood <- .lg_listing_value(lines, "^Log-likelihood \\(LL\\)")
    count <- .lg_listing_value(lines, "^Number of parameters \\(Npar\\)")
    rbind(
      .lg_rows(case, "log likelihood", likelihood$value %||% NA_real_, target$log_likelihood,
               (likelihood$half_unit %||% NA_real_) + .lg_tolerances[["likelihood"]],
               note = if (is.null(likelihood)) "No `Log-likelihood (LL)` line in the listing." else ""),
      .lg_rows(case, "parameters", count$value %||% NA_real_, target$n_parameters, 0,
               note = if (is.null(count)) "No `Number of parameters (Npar)` line in the listing." else ""))
  }

  individual_latent <- if (grepl("^c1[01]_", case)) .lg_latent_names[["state"]] else
    .lg_latent_names[["individual"]]
  table <- if (file.exists(posterior_path)) .lg_read_table(posterior_path)
  individual <- if (!is.null(table)) .lg_posterior_matrix(table, individual_latent)
  target_matrix <- as.matrix(target$posteriors[grep("^class_", names(target$posteriors))])

  posterior_rows <- if (!file.exists(posterior_path)) {
    .lg_rows(case, "posteriors: individuals", NA_real_, NA_real_, NA_real_,
             status = "missing", note = sprintf("No %s.", basename(posterior_path)))
  } else if (is.null(individual) || !"id" %in% names(table)) {
    .lg_rows(case, "posteriors: individuals", NA_real_, NA_real_, NA_real_,
             status = "unparsed",
             note = sprintf("Expected `id` and `%s#1`, ... columns; found: %s.", individual_latent,
                            if (is.null(table)) "no readable table" else
                              paste(names(table), collapse = ", ")))
  } else {
    rows <- match(target$posteriors$id, as.integer(table$id))
    if (anyNA(rows)) {
      .lg_rows(case, "posteriors: individuals", NA_real_, NA_real_, NA_real_,
               status = "unparsed", note = "Some ids of the data are absent from the posterior file.")
    } else {
      reference <- individual$values[rows, , drop = FALSE]
      alignment <- .lg_align(target_matrix, reference)
      aligned <- reference[, alignment$order, drop = FALSE]
      tolerance <- individual$half_unit + .lg_tolerances[["parameter"]]
      # The moments are the M-step only for complete data and a diagonal
      # covariance; with FIML or a full covariance they are not, and are not
      # recomputed.
      moments <- if (!is.null(target$profiles) && length(target$continuous) > 0L &&
                     identical(target$missing, "error") &&
                     identical(target$covariance_model, "diagonal")) {
        .lg_moment_rows(target, aligned, individual$half_unit)
      }
      rbind(
        .lg_rows(case, "posteriors: individuals (largest difference)", 0,
                 alignment$difference, tolerance,
                 note = sprintf("%d rows; Latent GOLD class order %s.", nrow(aligned),
                                paste(alignment$order, collapse = ","))),
        moments)
    }
  }

  group_rows <- if (!is.null(target$group_posteriors)) {
    group_latent <- .lg_latent_names[["group"]]
    groups <- if (!is.null(table)) .lg_posterior_matrix(table, group_latent)
    if (!file.exists(posterior_path)) {
      .lg_rows(case, "posteriors: groups", NA_real_, NA_real_, NA_real_, status = "missing",
               note = sprintf("No %s.", basename(posterior_path)))
    } else if (is.null(groups)) {
      .lg_rows(case, "posteriors: groups", NA_real_, NA_real_, NA_real_, status = "unparsed",
               note = sprintf("Expected `%s#1`, ... columns.", group_latent))
    } else {
      # Group posteriors repeat on every row of a group; the group of each row
      # comes from the data the kit wrote, so the file need not carry it.
      group_of_row <- target$data$g[match(as.integer(table$id), target$data$id)]
      first_rows <- match(target$group_posteriors$g, group_of_row)
      if (anyNA(first_rows)) {
        .lg_rows(case, "posteriors: groups", NA_real_, NA_real_, NA_real_, status = "unparsed",
                 note = "Some groups are absent from the posterior file.")
      } else {
        reference <- groups$values[first_rows, , drop = FALSE]
        target_groups <- as.matrix(target$group_posteriors[grep("^class_",
                                                                names(target$group_posteriors))])
        alignment <- .lg_align(target_groups, reference)
        aligned <- reference[, alignment$order, drop = FALSE]
        tolerance <- groups$half_unit + .lg_tolerances[["parameter"]]
        proportions <- if (!is.null(target$group_class_probabilities)) {
          .lg_rows(case, sprintf("group class %d proportion", seq_len(ncol(aligned))),
                   colMeans(aligned), unname(target$group_class_probabilities), tolerance,
                   note = "Mean of Latent GOLD's group posteriors.")
        }
        rbind(.lg_rows(case, "posteriors: groups (largest difference)", 0,
                       alignment$difference, tolerance,
                       note = sprintf("%d groups; Latent GOLD class order %s.", nrow(aligned),
                                      paste(alignment$order, collapse = ","))),
              proportions)
      }
    }
  }

  rbind(scalars, posterior_rows, group_rows, .lg_pending_rows(target))
}

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

#' The comparison statuses, in reporting order
#' @format A character vector.
.lg_statuses <- c("agree", "disagree", "missing", "unparsed", "awaiting parser")

#' The multilpa versions a set of targets was built with
#' @param targets A list of targets, as written by make-kit.R.
#' @return A character vector of distinct versions.
lg_target_versions <- function(targets) {
  unique(vapply(targets, function(target) target$multilpa_version, character(1)))
}

#' Compare every case
#'
#' @param targets A list of targets, as written by make-kit.R.
#' @param returned Directory holding the files Latent GOLD wrote.
#' @return A data frame, one row per compared quantity across every case.
lg_compare_all <- function(targets, returned) {
  comparison <- do.call(rbind, lapply(targets, lg_compare_case, returned = returned))
  row.names(comparison) <- NULL
  comparison
}

#' How many rows of each case have each status
#'
#' @param comparison From [lg_compare_all()].
#' @return A data frame, one row per case and one count column per status.
lg_status_counts <- function(comparison) {
  counts <- as.data.frame.matrix(table(factor(comparison$case),
                                       factor(comparison$status, levels = .lg_statuses)))
  data.frame(case = row.names(counts), counts, check.names = FALSE, row.names = NULL)
}

#' The rows where a Latent GOLD value was actually compared
#'
#' @param comparison From [lg_compare_all()].
#' @return The `agree` and `disagree` rows, without the note column.
lg_compared <- function(comparison) {
  compared <- comparison[comparison$status %in% c("agree", "disagree"), , drop = FALSE]
  compared[setdiff(names(compared), "note")]
}

#' One line summarising the whole comparison
#'
#' @param comparison From [lg_compare_all()].
#' @return A single string.
lg_totals <- function(comparison) {
  totals <- vapply(.lg_statuses, function(status) sum(comparison$status == status),
                   integer(1))
  paste(sprintf("%d %s", totals, names(totals)), collapse = ", ")
}

#' A data frame as the lines of a Markdown table
#' @param frame A data frame.
#' @return Character vector of table lines.
.lg_markdown_table <- function(frame) {
  cells <- vapply(frame, function(column) {
    if (is.numeric(column)) trimws(formatC(column, digits = 6, format = "g")) else
      as.character(column)
  }, character(nrow(frame)))
  cells <- matrix(cells, nrow(frame))
  cells[is.na(cells) | cells == "NA"] <- ""
  c(sprintf("| %s |", paste(names(frame), collapse = " | ")),
    sprintf("|%s|", paste(rep("---", ncol(frame)), collapse = "|")),
    sprintf("| %s |", apply(cells, 1L, paste, collapse = " | ")))
}

#' Write the comparison as CSV and as a Markdown report
#'
#' @param comparison From [lg_compare_all()].
#' @param directory Where to write `comparison.csv` and `REPORT.md`.
#' @param version The multilpa version the comparison was run under.
#' @return The two paths, invisibly.
lg_write_report <- function(comparison, directory, version) {
  csv <- file.path(directory, "comparison.csv")
  report <- file.path(directory, "REPORT.md")
  utils::write.csv(comparison, csv, row.names = FALSE)
  writeLines(c(
    "# Latent GOLD comparison",
    "",
    sprintf("multilpa %s against the output in `returned/`. Generated by `compare.R`;",
            version),
    "do not edit by hand.",
    "",
    sprintf("Totals: %s.", lg_totals(comparison)),
    "",
    "## Status by case",
    "",
    .lg_markdown_table(lg_status_counts(comparison)),
    "",
    "## Every row",
    "",
    .lg_markdown_table(comparison)), report)
  invisible(c(csv, report))
}

# ---------------------------------------------------------------------------
# Provenance of the retained Latent GOLD output
# ---------------------------------------------------------------------------

#' Fingerprint of the inputs Latent GOLD was given
#'
#' The retained output in `returned/` is evidence only for the data and syntax
#' it was produced from. This is the md5 of every `.dat` and `.lgs` of the kit,
#' so a rebuilt kit that differs can be detected instead of being compared
#' against stale output.
#'
#' @param kit The kit directory.
#' @return A named character vector, md5 by file name, sorted by name.
lg_fingerprint <- function(kit) {
  files <- sort(list.files(kit, pattern = "[.](dat|lgs)$"))
  stats::setNames(unname(tools::md5sum(file.path(kit, files))), files)
}

#' The fingerprint recorded beside the retained output
#'
#' @param returned The directory holding Latent GOLD's output.
#' @return The fingerprint as [lg_fingerprint()] returns it, or `NULL` when no
#'   record is there.
lg_recorded_fingerprint <- function(returned) {
  path <- file.path(returned, "INPUTS.md5")
  if (!file.exists(path)) return(NULL)
  recorded <- utils::read.table(path, header = FALSE, col.names = c("md5", "file"),
                                stringsAsFactors = FALSE)
  stats::setNames(recorded$md5, recorded$file)
}

#' Check the retained output against the kit it was produced from
#'
#' @param kit,returned The kit and output directories.
#' @return Invisibly `TRUE` when they match; otherwise a
#'   `multilpa_stale_latentgold_output` warning naming the files that differ.
lg_check_fingerprint <- function(kit, returned) {
  recorded <- lg_recorded_fingerprint(returned)
  if (is.null(recorded)) return(invisible(TRUE))
  current <- lg_fingerprint(kit)
  changed <- c(setdiff(names(current), names(recorded)),
               setdiff(names(recorded), names(current)),
               intersect(names(current), names(recorded))[
                 current[intersect(names(current), names(recorded))] !=
                   recorded[intersect(names(current), names(recorded))]])
  if (length(changed) == 0L) return(invisible(TRUE))
  warning(warningCondition(sprintf(
    "The kit has changed since Latent GOLD was run (%s). Rerun Latent GOLD on the current kit; the output in `returned/` is evidence for the previous one.",
    paste(sort(changed), collapse = ", ")),
    class = "multilpa_stale_latentgold_output"))
  invisible(FALSE)
}
