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
# - Standard errors, the natural-scale estimates they belong to, the ranking of
#   bivariate residuals and the Step-3 results, read from the listing's tables.
#   A standard error the listing cannot supply, or one that measures something
#   else, is reported as `not comparable` with the reason, never as agreement.

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
#' The listing is tab-separated with the label in the first field, so the label
#' is matched as that whole field rather than as a pattern: `AIC (based on LL)`
#' and `AIC3 (based on LL)` are different labels but one is a prefix of the
#' other, and a label ending in `)` has no word boundary after it.
#'
#' @param lines The listing's lines.
#' @param label The label, exactly as the listing writes it.
#' @return A list with `value` and `half_unit`, or `NULL` when the label is
#'   absent or carries no number (`Npar` heads a column of the summary table).
.lg_listing_value <- function(lines, label) {
  fields <- strsplit(lines, "\t", fixed = TRUE)
  hit <- Find(function(row) length(row) > 1L && identical(trimws(row[[1L]]), label),
              fields)
  if (is.null(hit)) return(NULL)
  numbers <- trimws(hit[-1L])
  numbers <- numbers[grepl("^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$", numbers)]
  if (length(numbers) == 0L) return(NULL)
  list(value = as.numeric(numbers[[1L]]), half_unit = .lg_half_unit(numbers[[1L]]))
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

# ---------------------------------------------------------------------------
# Tables in the listing: standard errors, bivariate residuals, Step-3
# ---------------------------------------------------------------------------
#
# Written against the listings Latent GOLD 6.1 produced on 2026-09-22, kept in
# returned/. Three facts about their layout decide how they are read:
#
# - A section is a line holding its title and nothing else (no tab), and runs
#   to the next such line. "Variances" appears twice, once among the
#   parameters and once among the paired comparisons, so the variance block is
#   looked for only inside "Parameters".
# - A parameter's name spans a varying number of fields: three in one listing
#   (`y1`, `<-`, `1`) and five in another (`y1`, `<-`, `1`, ``, ``), because a
#   model with class-specific variances widens the column for `| Cluster(1)`.
#   So the coefficient is found by content, not by column position.
# - Natural-scale values and their standard errors (class sizes, profile
#   probabilities, means) are printed under "EstimatedValues-Model" as small
#   blocks: a line naming the response, a line of column labels in which every
#   estimate is followed by `s.e.`, then one row per conditioning value.

#' Is each string a number as Latent GOLD prints one?
#' @param text Character vector.
#' @return A logical vector.
.lg_is_number <- function(text) {
  grepl("^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$", trimws(text))
}

#' Printed numbers as numbers, `NA` where a field is not one
#'
#' Latent GOLD prints `.` for a quantity it does not estimate. Converting only
#' the fields that are numbers keeps that apart from a conversion failure.
#'
#' @param text Character vector.
#' @return A numeric vector.
.lg_as_number <- function(text) {
  numbers <- rep(NA_real_, length(text))
  printed <- .lg_is_number(text)
  numbers[printed] <- as.numeric(trimws(text[printed]))
  numbers
}

#' The lines of one section of a listing
#'
#' @param lines The listing's lines.
#' @param title The section's title, exactly as printed on its own line.
#' @return The lines after the title up to the next title, or `character()`.
.lg_section <- function(lines, title) {
  titles <- which(!grepl("\t", lines, fixed = TRUE))
  first <- titles[trimws(lines[titles]) == title]
  if (length(first) == 0L) return(character())
  following <- titles[titles > first[[1L]]]
  last <- if (length(following) > 0L) following[[1L]] - 1L else length(lines)
  if (last <= first[[1L]]) character() else lines[(first[[1L]] + 1L):last]
}

#' Tab-separated fields of each line, trimmed
#' @param lines Character vector.
#' @return A list of character vectors.
.lg_fields <- function(lines) {
  lapply(strsplit(lines, "\t", fixed = TRUE), trimws)
}

#' The parameter tables of a listing, one row per printed parameter
#'
#' A reference category is printed as `0.0000` with `.` for its standard
#' error: it is fixed, not estimated, and its `standard_error` is `NA`.
#'
#' @param lines The listing's lines.
#' @return A data frame with `block` (`"regression"` or `"variance"`), `lhs`,
#'   `rhs` (`NA` for a variance), `condition` (e.g. `"Cluster(1)"`, or `NA`),
#'   `estimate`, `standard_error`, `estimate_text` and `standard_error_text`;
#'   `NULL` when the listing has no parameter section.
.lg_parameter_table <- function(lines) {
  section <- .lg_section(lines, "Parameters")
  if (length(section) == 0L) return(NULL)
  fields <- .lg_fields(section)
  heads <- vapply(fields, function(row) if (length(row) > 0L) row[[1L]] else "",
                  character(1))
  block <- cumsum(heads %in% c("Regression Parameters", "Variances"))
  block_name <- c("regression", "variance")[match(
    heads[heads %in% c("Regression Parameters", "Variances")], c("Regression Parameters",
                                                                "Variances"))]
  rows <- lapply(seq_along(fields), function(index) {
    row <- fields[[index]]
    if (block[[index]] == 0L || length(row) < 3L || !nzchar(row[[1L]]) ||
        row[[1L]] %in% c("term", "Regression Parameters", "Variances")) return(NULL)
    regression <- identical(row[[2L]], "<-")
    # The right-hand side can itself be `1`, the intercept, so it is taken by
    # position before the coefficient is looked for by content.
    rest <- if (regression) row[-(1:3)] else row[-1L]
    at <- which(.lg_is_number(rest))
    if (length(at) == 0L) return(NULL)
    label <- rest[seq_len(at[[1L]] - 1L)]
    label <- label[nzchar(label) & label != "|"]
    error_text <- if (length(rest) > at[[1L]]) rest[[at[[1L]] + 1L]] else "."
    data.frame(block = block_name[[block[[index]]]], lhs = row[[1L]],
               rhs = if (regression) row[[3L]] else NA_character_,
               condition = if (length(label) > 0L) paste(label, collapse = " ") else
                 NA_character_,
               estimate = .lg_as_number(rest[[at[[1L]]]]),
               standard_error = .lg_as_number(error_text),
               estimate_text = rest[[at[[1L]]]], standard_error_text = error_text,
               stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  if (is.null(table)) NULL else table
}

#' The natural-scale estimates of a listing, one row per printed value
#'
#' @param lines The listing's lines.
#' @return A data frame with `response`, `condition` (the conditioning
#'   variable, or `NA`), `condition_value`, `label` (a class or `"Mean"`),
#'   `estimate`, `standard_error` and their printed text; `NULL` when absent.
.lg_estimated_values <- function(lines) {
  fields <- .lg_fields(.lg_section(lines, "EstimatedValues-Model"))
  if (length(fields) < 3L) return(NULL)
  filled <- lapply(fields, function(row) which(nzchar(row)))
  # A block opens with a line naming only its response, followed by a line of
  # column labels holding at least one `s.e.`.
  opens <- which(vapply(seq_along(fields), function(index) {
    index < length(fields) && length(filled[[index]]) == 1L &&
      !.lg_is_number(fields[[index]][filled[[index]]]) &&
      "s.e." %in% fields[[index + 1L]]
  }, logical(1)))
  blocks <- lapply(opens, function(open) {
    response <- fields[[open]][filled[[open]]]
    header <- fields[[open + 1L]]
    errors <- which(header == "s.e.")
    estimates <- errors - 1L
    conditioning <- seq_len(estimates[[1L]] - 1L)
    condition_at <- conditioning[nzchar(header[conditioning])]
    condition <- if (length(condition_at) == 1L) header[[condition_at]] else NA_character_
    body <- seq(open + 2L, length.out = max(0L, length(fields) - open - 1L))
    ends <- which(vapply(body, function(index) {
      row <- fields[[index]]
      index %in% opens || !any(nzchar(row)) || length(row) < max(errors)
    }, logical(1)))
    body <- body[seq_len(if (length(ends) > 0L) ends[[1L]] - 1L else length(body))]
    do.call(rbind, lapply(body, function(index) {
      row <- fields[[index]]
      value <- if (length(condition_at) == 1L) row[[condition_at]] else NA_character_
      data.frame(response = response, condition = condition, condition_value = value,
                 label = header[estimates], estimate_text = row[estimates],
                 standard_error_text = row[errors],
                 estimate = .lg_as_number(row[estimates]),
                 standard_error = .lg_as_number(row[errors]),
                 stringsAsFactors = FALSE)
    }))
  })
  table <- do.call(rbind, blocks)
  if (is.null(table)) NULL else table
}

#' The bivariate residuals of the dependent variables
#'
#' Latent GOLD prints them as a lower triangle headed `Dependent`, with `.` on
#' the diagonal; the `Independent` and `Twolevel` blocks that follow are not
#' pairs of indicators and are not read.
#'
#' @param lines The listing's lines.
#' @return A data frame with `indicator_1`, `indicator_2` (in header order),
#'   `residual` and `residual_text`, or `NULL`.
.lg_bivariate_table <- function(lines) {
  fields <- .lg_fields(.lg_section(lines, "Bivariate Residuals"))
  header_at <- which(vapply(fields, function(row) identical(row[1L], "Dependent"),
                            logical(1)))
  if (length(header_at) == 0L) return(NULL)
  names_in <- fields[[header_at[[1L]]]][-1L]
  names_in <- names_in[nzchar(names_in)]
  rows <- fields[seq(header_at[[1L]] + 1L, length.out = length(names_in))]
  pairs <- do.call(rbind, lapply(seq_along(rows), function(i) {
    row <- rows[[i]]
    if (!identical(row[1L], names_in[[i]]) || i == 1L) return(NULL)
    data.frame(indicator_1 = names_in[seq_len(i - 1L)], indicator_2 = names_in[[i]],
               residual_text = row[1L + seq_len(i - 1L)], stringsAsFactors = FALSE)
  }))
  if (is.null(pairs) || !all(.lg_is_number(pairs$residual_text))) return(NULL)
  pairs$residual <- .lg_as_number(pairs$residual_text)
  pairs
}

#' A class number from a label such as `profile_2` or `group_class_1`
#' @param label Character vector.
#' @return An integer vector, `NA` where the label names no class.
.lg_class_number <- function(label) {
  numbered <- grepl("_[0-9]+$", label)
  classes <- rep(NA_integer_, length(label))
  classes[numbered] <- as.integer(sub("^.*_([0-9]+)$", "\\1", label[numbered]))
  classes
}

#' One of this package's logits, as Latent GOLD codes it
#'
#' Both programs model class membership as multinomial logits, but code them
#' against different classes: Latent GOLD with `parameters=first` against its
#' first class, this package against its last. The log odds between two
#' classes is the same in any coding, so this package's coefficient for class
#' `k` is Latent GOLD's for the class aligned to `k` minus Latent GOLD's for
#' the class aligned to the reference, summed over the right-hand-side terms
#' that make it up (a group-class intercept is Latent GOLD's intercept plus its
#' group-class effect). Its standard error is printed only when exactly one of
#' those terms is estimated rather than fixed at zero; otherwise it would need
#' their covariance, which the listing does not print.
#'
#' @param table From [.lg_parameter_table()].
#' @param latent The latent variable, e.g. `"Cluster"`.
#' @param class,reference Latent GOLD's class numbers for this package's class
#'   and its reference class.
#' @param rhs Latent GOLD's right-hand-side terms, e.g. `c("1", "GClass(2)")`.
#' @return A list with `estimate`, `half_unit`, `standard_error` (`NA` when
#'   not printed) and `standard_error_half_unit`, or `NULL` when a term is
#'   absent from the listing.
.lg_logit <- function(table, latent, class, reference, rhs) {
  wanted <- expand.grid(class = c(class, reference), rhs = rhs, stringsAsFactors = FALSE)
  wanted$sign <- ifelse(wanted$class == class, 1, -1)
  hits <- match(paste(sprintf("%s(%d)", latent, wanted$class), wanted$rhs),
                paste(table$lhs, table$rhs)[table$block == "regression"])
  if (anyNA(hits)) return(NULL)
  found <- table[table$block == "regression", , drop = FALSE][hits, , drop = FALSE]
  free <- !is.na(found$standard_error)
  list(estimate = sum(wanted$sign * found$estimate),
       half_unit = sum(vapply(found$estimate_text, .lg_half_unit, numeric(1))),
       standard_error = if (sum(free) == 1L) found$standard_error[free] else NA_real_,
       standard_error_half_unit = if (sum(free) == 1L)
         .lg_half_unit(found$standard_error_text[free]) else NA_real_)
}

#' Latent GOLD's value and standard error for one of this package's parameters
#'
#' @param parameter One row of a target's `parameters`.
#' @param table,values From [.lg_parameter_table()] and
#'   [.lg_estimated_values()].
#' @param alignment A list with `individual` and `group`, each giving Latent
#'   GOLD's class for each of this package's classes.
#' @param n_classes A list with `individual` and `group` class counts.
#' @return A list as [.lg_logit()] returns, or a string saying why the
#'   parameter could not be found.
.lg_locate_parameter <- function(parameter, table, values, alignment, n_classes) {
  individual <- .lg_latent_names[["individual"]]
  group <- .lg_latent_names[["group"]]
  class <- .lg_class_number(parameter$outcome)
  from_values <- function(hit) {
    if (length(hit) != 1L) return("Not found under EstimatedValues-Model.")
    list(estimate = values$estimate[[hit]],
         half_unit = .lg_half_unit(values$estimate_text[[hit]]),
         standard_error = values$standard_error[[hit]],
         standard_error_half_unit = .lg_half_unit(values$standard_error_text[[hit]]))
  }
  kind <- paste(parameter$level, parameter$parameter)
  if (identical(kind, "measurement mean")) {
    return(from_values(which(values$response == parameter$term &
                               values$condition %in% individual &
                               values$condition_value == alignment$individual[[class]] &
                               values$label == "Mean")))
  }
  if (identical(kind, "measurement variance")) {
    condition <- if (is.na(class)) NA_character_ else
      sprintf("%s(%d)", individual, alignment$individual[[class]])
    hit <- which(table$block == "variance" & table$lhs == parameter$term &
                   (table$condition %in% condition |
                      (is.na(condition) & is.na(table$condition))))
    if (length(hit) != 1L) return("Not found among the listing's variances.")
    return(list(estimate = table$estimate[[hit]],
                half_unit = .lg_half_unit(table$estimate_text[[hit]]),
                standard_error = table$standard_error[[hit]],
                standard_error_half_unit = .lg_half_unit(table$standard_error_text[[hit]])))
  }
  if (identical(kind, "profile probability")) {
    group_class <- .lg_class_number(parameter$term)
    return(from_values(which(values$response == individual & values$condition %in% group &
                               values$condition_value == alignment$group[[group_class]] &
                               values$label == alignment$individual[[class]])))
  }
  if (identical(kind, "group probability")) {
    return(from_values(which(values$response == group & is.na(values$condition) &
                               values$label == alignment$group[[class]])))
  }
  if (identical(parameter$parameter, "coefficient")) {
    level <- if (identical(parameter$level, "group")) "group" else "individual"
    latent <- if (identical(level, "group")) group else individual
    group_class <- if (grepl("^group_class_[0-9]+$", parameter$term))
      .lg_class_number(parameter$term)
    rhs <- if (!is.null(group_class)) {
      c("1", sprintf("%s(%d)", group, alignment$group[[group_class]]))
    } else if (identical(parameter$term, "(Intercept)")) "1" else parameter$term
    located <- .lg_logit(table, latent, alignment[[level]][[class]],
                         alignment[[level]][[n_classes[[level]]]], rhs)
    return(located %||% "Not found among the listing's regression parameters.")
  }
  sprintf("No rule maps a %s %s to the listing.", parameter$level, parameter$parameter)
}

#' An estimate row and a standard-error row for each located quantity
#'
#' @param case,labels Case name and one label per quantity.
#' @param located A list, one element per quantity, from
#'   [.lg_locate_parameter()] or [.lg_logit()].
#' @param estimates,errors This package's estimates and standard errors.
#' @param error_note Why a standard error that is printed cannot be compared,
#'   or `NULL` when it can.
#' @return Comparison rows.
.lg_located_rows <- function(case, labels, located, estimates, errors, error_note = NULL) {
  do.call(rbind, lapply(seq_along(located), function(index) {
    found <- located[[index]]
    if (is.character(found)) {
      return(.lg_rows(case, paste(c("estimate:", "standard error:"), labels[[index]]),
                      NA_real_, c(estimates[[index]], errors[[index]]), NA_real_,
                      status = "unparsed", note = found))
    }
    estimate <- .lg_rows(case, paste("estimate:", labels[[index]]), found$estimate,
                         estimates[[index]],
                         found$half_unit + .lg_tolerances[["parameter"]],
                         note = "Read from the listing.")
    error <- if (!is.null(error_note)) {
      .lg_rows(case, paste("standard error:", labels[[index]]), found$standard_error,
               errors[[index]], NA_real_, status = "not comparable", note = error_note)
    } else if (is.na(found$standard_error)) {
      .lg_rows(case, paste("standard error:", labels[[index]]), NA_real_, errors[[index]],
               NA_real_, status = "not comparable",
               note = "Latent GOLD codes this logit as a sum of two estimated terms; its standard error needs their covariance, which the listing does not print.")
    } else {
      .lg_rows(case, paste("standard error:", labels[[index]]), found$standard_error,
               errors[[index]],
               found$standard_error_half_unit + .lg_tolerances[["parameter"]],
               note = "Read from the listing.")
    }
    rbind(estimate, error)
  }))
}

#' Rows for the parameters and their standard errors
#'
#' @param target One case's targets.
#' @param lines The listing's lines.
#' @param alignment,n_classes As for [.lg_locate_parameter()].
#' @return Comparison rows, or `NULL` when the case has no parameter targets.
.lg_parameter_rows <- function(target, lines, alignment, n_classes) {
  parameters <- target$parameters
  if (is.null(parameters)) return(NULL)
  labels <- trimws(gsub(" +", " ", paste(parameters$level, parameters$outcome,
                                         ifelse(is.na(parameters$term), "", parameters$term),
                                         parameters$parameter)))
  table <- .lg_parameter_table(lines)
  values <- .lg_estimated_values(lines)
  located <- if (is.null(table) || is.null(values)) {
    rep(list("The listing has no parameter or EstimatedValues-Model section."),
        nrow(parameters))
  } else if (is.null(alignment$individual) ||
             (any(parameters$level %in% c("profile", "group")) && is.null(alignment$group))) {
    rep(list("Class labels could not be aligned, because the posteriors were not compared."),
        nrow(parameters))
  } else {
    lapply(seq_len(nrow(parameters)), function(index) {
      .lg_locate_parameter(parameters[index, , drop = FALSE], table, values, alignment,
                           n_classes)
    })
  }
  .lg_located_rows(target$case, labels, located, parameters$estimate,
                   parameters$standard_error)
}

#' Rows comparing how the two programs rank local dependence
#'
#' The statistics differ. Latent GOLD prints one residual per indicator pair,
#' pooled over classes; this package tests each pair within each profile. So
#' the values are not compared. What is compared is the order: each pair is
#' ranked by Latent GOLD's residual and by this package's largest absolute
#' within-profile statistic, and the two rankings must agree.
#'
#' @param target One case's targets.
#' @param lines The listing's lines.
#' @return Comparison rows, or `NULL` when the case has no residual targets.
.lg_residual_rows <- function(target, lines) {
  residuals <- target$residuals
  if (is.null(residuals)) return(NULL)
  ours <- tapply(abs(residuals$statistic),
                 paste(pmin(residuals$indicator_1, residuals$indicator_2),
                       pmax(residuals$indicator_1, residuals$indicator_2), sep = "-"),
                 max)
  theirs <- .lg_bivariate_table(lines)
  labels <- sprintf("bivariate residual rank: %s", names(ours))
  if (is.null(theirs)) {
    return(.lg_rows(target$case, labels, NA_real_, as.numeric(ours), NA_real_,
                    status = "unparsed", note = "No Dependent block under Bivariate Residuals."))
  }
  pairs <- paste(pmin(theirs$indicator_1, theirs$indicator_2),
                 pmax(theirs$indicator_1, theirs$indicator_2), sep = "-")
  printed <- theirs$residual[match(names(ours), pairs)]
  if (anyNA(printed)) {
    return(.lg_rows(target$case, labels, NA_real_, as.numeric(ours), NA_real_,
                    status = "unparsed", note = "A pair of this package's is absent from the listing."))
  }
  .lg_rows(target$case, labels, rank(-printed, ties.method = "min"),
           rank(-as.numeric(ours), ties.method = "min"), 0,
           note = sprintf("Rank 1 is the largest. Latent GOLD's residual %s; this package's largest |statistic| %.4f.",
                          theirs$residual_text[match(names(ours), pairs)], as.numeric(ours)))
}

#' The distal outcome's name in the three-step case, as `cases.R` writes it
#' @format A string.
.lg_distal_outcome <- "distal"

#' Rows for the Step-3 analyses
#'
#' The class means of a distal outcome are compared as estimates only. Their
#' standard errors are shown but not judged, because the two programs compute
#' different ones: this package's are cluster-robust, which on this
#' single-level case is the Huber-White sandwich, while Latent GOLD's modal and
#' proportional models assume one variance shared by the classes.
#' The R3STEP logits are compared with their standard errors, which are the
#' same model-based quantity in both programs.
#'
#' @param target One case's targets.
#' @param returned Directory holding the files Latent GOLD wrote.
#' @param alignment As for [.lg_locate_parameter()].
#' @return Comparison rows, or `NULL` when the case has no Step-3 targets.
.lg_step_rows <- function(target, returned, alignment) {
  steps <- target$steps
  if (is.null(steps)) return(NULL)
  case <- target$case
  read <- function(name) {
    path <- file.path(returned, sprintf("c09_%s.lst", name))
    if (file.exists(path)) .lg_listing_lines(path)
  }
  unavailable <- function(labels, estimates, errors, status, note) {
    .lg_rows(case, c(paste("estimate:", labels), paste("standard error:", labels)),
             NA_real_, c(estimates, errors), NA_real_, status = status, note = note)
  }
  means <- lapply(c("distal_bch", "distal_modal", "distal_proportional"), function(name) {
    table <- steps[[name]]
    labels <- sprintf("%s: class %d mean", name, table$class)
    lines <- read(name)
    if (is.null(lines)) {
      return(unavailable(labels, table$estimate, table$standard_error, "missing",
                         sprintf("No c09_%s.lst.", name)))
    }
    if (is.null(alignment$individual)) {
      return(unavailable(labels, table$estimate, table$standard_error, "unparsed",
                         "Class labels could not be aligned, because the posteriors were not compared."))
    }
    values <- .lg_estimated_values(lines)
    located <- lapply(table$class, function(class) {
      hit <- which(values$response == .lg_distal_outcome &
                     values$condition %in% .lg_latent_names[["individual"]] &
                     values$condition_value == alignment$individual[[class]] &
                     values$label == "Mean")
      if (length(hit) != 1L) return("Not found under EstimatedValues-Model.")
      list(estimate = values$estimate[[hit]],
           half_unit = .lg_half_unit(values$estimate_text[[hit]]),
           standard_error = values$standard_error[[hit]],
           standard_error_half_unit = .lg_half_unit(values$standard_error_text[[hit]]))
    })
    .lg_located_rows(case, labels, located, table$estimate, table$standard_error,
                     error_note = sprintf("This package's is %s; Latent GOLD's is model-based%s.",
                                          attr(table, "vcov_type") %||% "not recorded",
                                          if (identical(name, "distal_bch")) "" else
                                            ", with one variance shared by the classes"))
  })
  slopes <- steps$covariate_ml
  labels <- sprintf("covariate_ml: %s %s", slopes$outcome, slopes$term)
  lines <- read("covariate_ml")
  logits <- if (is.null(lines)) {
    unavailable(labels, slopes$estimate, slopes$standard_error, "missing",
                "No c09_covariate_ml.lst.")
  } else if (is.null(alignment$individual)) {
    unavailable(labels, slopes$estimate, slopes$standard_error, "unparsed",
                "Class labels could not be aligned, because the posteriors were not compared.")
  } else {
    table <- .lg_parameter_table(lines)
    reference <- attr(slopes, "reference_class") %||% length(alignment$individual)
    located <- lapply(seq_len(nrow(slopes)), function(index) {
      if (is.null(table)) return("The listing has no parameter section.")
      rhs <- if (identical(slopes$term[[index]], "(Intercept)")) "1" else slopes$term[[index]]
      .lg_logit(table, .lg_latent_names[["individual"]],
                alignment$individual[[.lg_class_number(slopes$outcome[[index]])]],
                alignment$individual[[reference]], rhs) %||%
        "Not found among the listing's regression parameters."
    })
    .lg_located_rows(case, labels, located, slopes$estimate, slopes$standard_error)
  }
  do.call(rbind, c(means, list(logits)))
}

#' Every row read from the listing's tables
#'
#' @param target One case's targets.
#' @param returned Directory holding the files Latent GOLD wrote.
#' @param lines The case's listing lines, or `NULL` when there is no listing.
#' @param alignment,n_classes As for [.lg_locate_parameter()].
#' @return Comparison rows, or `NULL`.
.lg_table_rows <- function(target, returned, lines, alignment, n_classes) {
  listing_rows <- if (is.null(lines)) {
    labels <- c(if (!is.null(target$parameters)) "parameters and standard errors",
                if (!is.null(target$residuals)) "bivariate residuals")
    if (length(labels) > 0L) {
      .lg_rows(target$case, labels, NA_real_, NA_real_, NA_real_, status = "missing",
               note = "No listing.")
    }
  } else {
    rbind(.lg_parameter_rows(target, lines, alignment, n_classes),
          .lg_residual_rows(target, lines))
  }
  rbind(listing_rows, .lg_step_rows(target, returned, alignment))
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

#' Information criteria, by the sample size each is computed on
#'
#' Latent GOLD labels a criterion `(based on LL)` when it uses the number of
#' cases and `(based on LL,Ngroups)` when it uses the number of groups. This
#' package draws the same distinction as `*_individual` and `*_groups`. The two
#' vocabularies do not line up by name, because what Latent GOLD calls a case
#' depends on the model: for a two-level model it is an observation, but for a
#' sequence model declared with `caseid` it is a whole sequence, which is this
#' package's *group*. So each criterion is matched to the column computed on
#' the same sample size, read from the listing, rather than by its label.
#'
#' Deliberately not compared, because the definitions differ rather than the
#' values: Latent GOLD's `AIC3` (this package reports `kic`, which is not
#' `-2LL + 3p`), and its `CLC`, `AWE`, entropy R-squared and classification
#' errors, which use different entropy conventions.
#'
#' @format A named character vector: the listing's label, and the stem of this
#'   package's column.
.lg_criteria <- c("AIC (based on LL)" = "aic",
                  "BIC (based on LL)" = "bic",
                  "CAIC (based on LL)" = "caic",
                  "SABIC (based on LL)" = "sabic")

#' Which convention of this package a Latent GOLD sample size is
#'
#' @param size The sample size Latent GOLD reports, or `NULL`.
#' @param target One case's targets.
#' @return `"individual"`, `"groups"`, or `NA` when it is neither.
.lg_convention <- function(size, target) {
  if (is.null(size) || is.null(target$n_observations)) return(NA_character_)
  if (isTRUE(all.equal(size, as.numeric(target$n_observations)))) return("individual")
  if (isTRUE(all.equal(size, as.numeric(target$n_groups)))) return("groups")
  NA_character_
}

#' Comparison rows for the information criteria and the sample sizes
#'
#' Each criterion is `-2 * LL` plus a penalty in the parameter count and the
#' sample size, both compared separately and exactly, so the only slack is
#' twice the likelihood's plus the criterion's own printed rounding.
#'
#' @param lines The listing's lines.
#' @param target One case's targets.
#' @return Comparison rows, or `NULL` when the targets predate the criteria.
.lg_criteria_rows <- function(lines, target) {
  if (is.null(target$criteria)) return(NULL)
  likelihood <- .lg_listing_value(lines, "Log-likelihood (LL)")
  slack <- 2 * ((likelihood$half_unit %||% 0) + .lg_tolerances[["likelihood"]])
  cases <- .lg_listing_value(lines, "Number of cases")
  groups <- .lg_listing_value(lines, "Number of groups")
  families <- list(list(suffix = "", size = cases, label = "cases"),
                   list(suffix = ",Ngroups", size = groups, label = "groups"))
  rows <- lapply(families, function(family) {
    # Latent GOLD prints the group-based family only for a model that has a
    # group level; its absence is the model's shape, not a parse failure.
    if (is.null(family$size)) return(NULL)
    convention <- .lg_convention(family$size$value, target)
    if (is.na(convention)) {
      return(.lg_rows(target$case, sprintf("sample size (%s)", family$label),
                      family$size$value, target$n_observations, 0,
                      note = "Latent GOLD's sample size is neither this fit's observations nor its groups."))
    }
    criteria <- lapply(names(.lg_criteria), function(label) {
      printed <- .lg_listing_value(lines, sub("\\)$", paste0(family$suffix, ")"), label))
      if (is.null(printed)) return(NULL)
      column <- if (identical(.lg_criteria[[label]], "aic")) "aic" else
        paste0(.lg_criteria[[label]], "_", convention)
      if (is.null(target$criteria[[column]])) return(NULL)
      .lg_rows(target$case, column, printed$value, target$criteria[[column]],
               printed$half_unit + slack,
               note = sprintf("Latent GOLD's `%s`, on %d %s.",
                              sub("\\)$", paste0(family$suffix, ")"), label),
                              as.integer(family$size$value), family$label))
    })
    size_row <- .lg_rows(target$case,
                         sprintf("n_%s", if (identical(convention, "individual"))
                           "observations" else "groups"),
                         family$size$value,
                         if (identical(convention, "individual")) target$n_observations else
                           target$n_groups, 0,
                         note = sprintf("Latent GOLD's `Number of %s`.", family$label))
    do.call(rbind, c(criteria, list(size_row)))
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

  lines <- NULL
  alignment <- list(individual = NULL, group = NULL)
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
    likelihood <- .lg_listing_value(lines, "Log-likelihood (LL)")
    count <- .lg_listing_value(lines, "Number of parameters (Npar)")
    rbind(
      .lg_criteria_rows(lines, target),
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
      matched <- .lg_align(target_matrix, reference)
      alignment$individual <- matched$order
      aligned <- reference[, matched$order, drop = FALSE]
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
                 matched$difference, tolerance,
                 note = sprintf("%d rows; Latent GOLD class order %s.", nrow(aligned),
                                paste(matched$order, collapse = ","))),
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
        matched <- .lg_align(target_groups, reference)
        alignment$group <- matched$order
        aligned <- reference[, matched$order, drop = FALSE]
        tolerance <- groups$half_unit + .lg_tolerances[["parameter"]]
        proportions <- if (!is.null(target$group_class_probabilities)) {
          .lg_rows(case, sprintf("group class %d proportion", seq_len(ncol(aligned))),
                   colMeans(aligned), unname(target$group_class_probabilities), tolerance,
                   note = "Mean of Latent GOLD's group posteriors.")
        }
        rbind(.lg_rows(case, "posteriors: groups (largest difference)", 0,
                       matched$difference, tolerance,
                       note = sprintf("%d groups; Latent GOLD class order %s.", nrow(aligned),
                                      paste(matched$order, collapse = ","))),
              proportions)
      }
    }
  }

  n_classes <- list(individual = ncol(target_matrix),
                    group = if (!is.null(target$group_posteriors))
                      length(grep("^class_", names(target$group_posteriors))))
  rbind(scalars, posterior_rows, group_rows,
        .lg_table_rows(target, returned, lines, alignment, n_classes))
}

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

#' The comparison statuses, in reporting order
#' @format A character vector.
.lg_statuses <- c("agree", "disagree", "not comparable", "missing", "unparsed")

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
