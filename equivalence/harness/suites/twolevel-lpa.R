# The default model of the package, the two-level latent profile model with
# continuous indicators (Vermunt 2003; profiles within discrete group
# classes), against Latent GOLD 6.1 and against a likelihood written here from
# its definition.
#
# Reference
# ---------
# Latent GOLD 6.1 is the reference implementation of this model. Its runs are
# retained in equivalence/latentgold/returned/ (see PROVENANCE.md there): for
# each case a listing (`.lst`, estimates and standard errors printed to four
# decimals, likelihood and criteria to four) and a posterior file (both
# levels' posteriors to six significant digits, with the data). The data are
# the kit files Latent GOLD read (equivalence/latentgold/kit/), checked against
# the md5 fingerprint recorded with the output, so the package is fitted to
# exactly the rows Latent GOLD saw. The listings are read with the kit's own
# parsers (equivalence/latentgold/compare-functions.R).
#
# Five cases, two group classes and two profiles each:
#
#   c01_continuous_varying  three indicators, diagonal, class-specific variances
#   c02_continuous_equal    as c01, one variance per indicator (equal)
#   c03_full_covariance     two correlated indicators, full class-specific
#                           covariance
#   c06_missing_fiml        c01 with 10% of values missing at random, FIML
#   c12_course_engagement   the bundled `course_engagement` data (1,422 rows in
#                           106 students, five indicators), the README model
#
# Correspondence of parameterisations
# -----------------------------------
#   Latent GOLD GClass (group, nominal)   <-> group classes; `GClass <- 1` is
#       P(group class), the package's group_probabilities (logits coded
#       against the first class with `parameters=first`; the natural values
#       are read from `EstimatedValues-Model`)
#   Latent GOLD Cluster (case, nominal), `Cluster <- 1 + GClass` <->
#       profiles with P(profile | group class), the package's
#       profile_probabilities
#   `y <- 1 + Cluster` <-> profile means (read from the `Profile` section,
#       which prints each profile's mean with its standard error)
#   `y | Cluster` <-> class-specific variances (`variance_model = "varying"`);
#       no such line <-> one variance per indicator (`"equal"`);
#       `y1 <-> y2 | Cluster` <-> `covariance_model = "full"`
#   `missing includeall` <-> `missing = "fiml"`; `bayes ... = 0` makes Latent
#       GOLD's estimates maximum likelihood, as the package's are
#   `standarderrors` (Latent GOLD's "standard" errors) <-> the inverse observed
#       information, `parameter_inference(vcov_type = "observed")`
#   Latent GOLD's `Entropy` is the entropy of the joint (group class, profile)
#       classification, which the package charges in CLC, ICL and AWE; its
#       `Entropy R-squared` divides by the entropy of the estimated class
#       sizes rather than by log K, so it is compared after rebuilding it from
#       the package's `entropy_sum` and class proportions.
#
# Two cases per dataset, as the Mplus two-level suite does:
#
#   1. At Latent GOLD's printed estimates. The package is scored there
#      (`max_iter = 0`, started at the printed values) and every output it
#      reports is compared, to 1e-10, with the from-scratch two-level Gaussian
#      likelihood (twolevel-lpa-reference.R) at the same parameters. Latent
#      GOLD's printed log-likelihood is compared with the package's at its
#      printed estimates, allowing for how far rounding the estimates to four
#      decimals can move it. Standard errors are withheld at a point the
#      package did not converge to, and are compared in case 2.
#   2. Converged. Latent GOLD's printed estimates are carried to the maximum
#      by Newton steps on the from-scratch likelihood and then by the
#      from-scratch EM to its fixed point. The package's own fit
#      (five random starts) is compared with that reference through
#      compare_fit(), every output, with standard errors from a numDeriv
#      Hessian of the same likelihood. `parameter_inference()` is compared in
#      full, observed and robust (the sandwich of numDeriv per-group scores
#      times G / (G - 1), as `?parameter_inference` documents). Everything
#      Latent GOLD printed that the package reports -- likelihood, criteria,
#      classification statistics and modal classification tables at both
#      levels, class sizes, estimates and their standard errors, and every
#      posterior at both levels -- is compared with the converged fit at half a
#      unit in the printed place, plus the package's measured stopping gap:
#      each printed quantity computed from scratch at the package's own
#      estimates and at the maximum, cell by cell.
#
# Reference checks. The from-scratch likelihood must reproduce, at the refined
# maximum, Latent GOLD's printed likelihood, its posteriors at both levels,
# its within-group-class row posteriors (the `joint2` columns, from which the
# row entropies of the complete-data entropy are built), its complete-data
# entropy and its standard errors, each to its printed digits. These rows
# carry no output: they show that the reference is Latent GOLD's model and
# Latent GOLD's maximum, not one of this file's own making.
#
# Vermunt, J. K. (2003). Multilevel latent class models. Sociological
# Methodology, 33(1), 213-239.
# Vermunt, J. K., & Magidson, J. (2021). Upgrade manual for Latent GOLD basic,
# advanced, syntax, and choice version 6.0. Statistical Innovations.

source(file.path("equivalence", "harness", "two-level-lca-reference.R"), local = environment())
source(file.path("equivalence", "harness", "twolevel-lpa-reference.R"), local = environment())

#' The Latent GOLD cases this suite compares
#' @return A named list, one element per case: its kit data file, indicators,
#'   covariance structure and missing-data handling.
.twolevel_lpa_cases <- function() {
  three <- c("y1", "y2", "y3")
  list(
    c01_continuous_varying = list(data_file = "c01.dat", indicators = three,
                                  structure = "varying", missing = "error"),
    c02_continuous_equal = list(data_file = "c01.dat", indicators = three,
                                structure = "equal", missing = "error"),
    c03_full_covariance = list(data_file = "c03.dat", indicators = c("y1", "y2"),
                               structure = "full", missing = "error"),
    c06_missing_fiml = list(data_file = "c06.dat", indicators = three,
                            structure = "varying", missing = "fiml"),
    c12_course_engagement = list(
      data_file = "c12.dat",
      indicators = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
      structure = "varying", missing = "error"))
}

suite_twolevel_lpa <- function() {
  require_suite_packages("numDeriv",
                         reason = "the reference observed information and per-group scores")
  directory <- file.path("equivalence", "latentgold")
  cases <- .twolevel_lpa_cases()
  data_files <- unique(vapply(cases, `[[`, character(1), "data_file"))
  require_suite_files(
    file.path(directory, "compare-functions.R"),
    file.path(directory, "kit", data_files),
    file.path(directory, "returned", c("INPUTS.md5", paste0(names(cases), ".lst"),
                                       paste0(names(cases), "_posteriors.txt"))),
    reason = "Latent GOLD 6.1's retained listings and posterior files, the kit data they were estimated from, and the kit's parsers")
  # The kit's listing parsers, into the environment this suite was sourced
  # into, where the functions below look them up.
  sys.source(file.path(directory, "compare-functions.R"), envir = environment(sys.function()))
  # The retained output is evidence for the data it was produced from and no
  # other: every data file read here must be the one Latent GOLD read.
  recorded <- lg_recorded_fingerprint(file.path(directory, "returned"))
  current <- lg_fingerprint(file.path(directory, "kit"))
  stopifnot("the kit data no longer match the data Latent GOLD was run on" =
              identical(unname(current[data_files]), unname(recorded[data_files])))
  stack_cases(lapply(names(cases), function(name) {
    do.call(.twolevel_lpa_case, c(list(case = name, directory = directory), cases[[name]]))
  }))
}

# ---------------------------------------------------------------------------
# Reading Latent GOLD's listing and posterior file
# ---------------------------------------------------------------------------

#' Half a unit in the last printed place of each printed number
#'
#' The posterior file prints six significant digits (`0.998227`,
#' `1.47299e-05`), so each value has its own resolution; the exponent counts.
#'
#' @param text Character vector of numbers as printed.
#' @return One bound per value.
.twolevel_half_units <- function(text) {
  text <- trimws(text)
  stopifnot("every printed value must be a number" = all(.lg_is_number(text)))
  mantissa <- sub("[eE].*$", "", text)
  scientific <- grepl("[eE]", text)
  exponent <- integer(length(text))
  exponent[scientific] <- as.integer(sub("^.*[eE]", "", text[scientific]))
  decimals <- ifelse(grepl(".", mantissa, fixed = TRUE), nchar(sub("^[^.]*[.]", "", mantissa)), 0L)
  0.5 * 10^(exponent - decimals)
}

#' The printed text of one labelled row of the listing, at given fields
#' @param lines The listing's lines.
#' @param label The row's first field.
#' @param columns Field positions to return.
#' @param after Line index to search after.
#' @return Character vector of the printed fields.
.twolevel_lg_row <- function(lines, label, columns, after = 0L) {
  fields <- .lg_fields(lines)
  hits <- which(vapply(fields, function(row) length(row) > 0L && identical(row[[1L]], label),
                       logical(1)))
  hits <- hits[hits > after]
  stopifnot("the listing prints the labelled row" = length(hits) >= 1L)
  row <- fields[[hits[[1L]]]]
  stopifnot("the labelled row has the requested fields" = max(columns) <= length(row))
  row[columns]
}

#' Latent GOLD's printed estimates and standard errors, labelled as the
#' package labels its parameters
#'
#' Group-class probabilities and profile probabilities within group class come
#' from `EstimatedValues-Model`, means from the `Profile` section, variances
#' and covariances from the parameter table, all as printed (four decimals).
#' Probabilities are renormalised to sum to one, after checking that their
#' printed sums depart from one by no more than rounding explains; the shift is
#' returned so that the rounding allowance can include it.
#'
#' @param lines The listing's lines.
#' @param layout From [twolevel_lpa_layout()].
#' @return A list with `parameters` (natural, in Latent GOLD's class
#'   numbering), `shift` (renormalisation shift of each free coordinate) and
#'   `printed`, a `data.frame` with `level`, `outcome`, `term`, `parameter`,
#'   `estimate_text` and `standard_error_text`.
.twolevel_lg_estimates <- function(lines, layout) {
  k <- layout$k
  g <- layout$g
  indicators <- layout$indicators
  values <- .lg_estimated_values(lines)
  stopifnot("the listing has an EstimatedValues-Model section" = !is.null(values))
  pick <- function(response, condition_value, label) {
    hit <- which(values$response == response & values$label == label &
                   (if (is.na(condition_value)) is.na(values$condition_value) else
                     values$condition_value %in% condition_value))
    stopifnot("each estimated value is printed once" = length(hit) == 1L)
    values[hit, c("estimate_text", "standard_error_text")]
  }
  delta <- do.call(rbind, lapply(seq_len(g), function(w) {
    data.frame(level = "group", outcome = sprintf("group_class_%d", w), term = "",
               parameter = "probability", pick("GClass", NA, as.character(w)))
  }))
  gamma <- do.call(rbind, lapply(seq_len(g), function(w) {
    do.call(rbind, lapply(seq_len(k), function(profile) {
      data.frame(level = "profile", outcome = sprintf("profile_%d", profile),
                 term = sprintf("group_class_%d", w), parameter = "probability",
                 pick("Cluster", as.character(w), as.character(profile)))
    }))
  }))
  # Means: the Profile section prints, per indicator, a `Mean` row holding
  # each GClass's, each Cluster's and the overall mean, each with its s.e.
  profile_fields <- .lg_fields(.lg_section(lines, "Profile"))
  cluster_at <- which(profile_fields[[1L]] == "Cluster")[[1L]]
  columns <- cluster_at + 2L * (seq_len(k) - 1L)
  stopifnot("the Profile section's Cluster columns are the classes with their s.e." =
              identical(profile_fields[[2L]][columns], as.character(seq_len(k))) &&
              all(profile_fields[[2L]][columns + 1L] == "s.e."))
  means <- do.call(rbind, lapply(indicators, function(indicator) {
    at <- which(vapply(profile_fields, function(row) identical(row[1L], indicator), logical(1)))
    stopifnot("each indicator has one Mean row in the Profile section" =
                length(at) == 1L && identical(profile_fields[[at + 1L]][1L], "Mean"))
    row <- profile_fields[[at + 1L]]
    data.frame(level = "measurement", outcome = sprintf("profile_%d", seq_len(k)),
               term = indicator, parameter = "mean", estimate_text = row[columns],
               standard_error_text = row[columns + 1L])
  }))
  table <- .lg_parameter_table(lines)
  variance_rows <- table[table$block == "variance", , drop = FALSE]
  variance_pattern <- "^Cluster\\(([0-9]+)\\)$"
  covariance_pattern <- "^<-> (\\S+) Cluster\\(([0-9]+)\\)$"
  covariance <- do.call(rbind, lapply(seq_len(nrow(variance_rows)), function(index) {
    row <- variance_rows[index, , drop = FALSE]
    condition <- row$condition
    off_diagonal <- !is.na(condition) && grepl(covariance_pattern, condition)
    profile <- if (is.na(condition)) NA_integer_ else if (off_diagonal)
      as.integer(sub(covariance_pattern, "\\2", condition)) else
        as.integer(sub(variance_pattern, "\\1", condition))
    other <- if (off_diagonal) sub(covariance_pattern, "\\1", condition) else row$lhs
    # The package names a covariance element `later:earlier` in indicator order.
    pair <- c(row$lhs, other)[order(match(c(row$lhs, other), indicators), decreasing = TRUE)]
    full <- identical(layout$structure, "full")
    data.frame(level = "measurement",
               outcome = if (is.na(profile)) "shared" else sprintf("profile_%d", profile),
               term = if (full) paste(pair, collapse = ":") else row$lhs,
               parameter = if (full) "covariance" else "variance",
               estimate_text = row$estimate_text, standard_error_text = row$standard_error_text)
  }))
  printed <- rbind(means, covariance, gamma, delta)
  printed$estimate <- as.numeric(printed$estimate_text)
  stopifnot("Latent GOLD printed every parameter of the model once" =
              nrow(printed) == nrow(twolevel_lpa_natural(numeric(layout$q), layout)$labels) &&
              !anyDuplicated(printed[c("level", "outcome", "term", "parameter")]))

  # The natural parameters, in Latent GOLD's class numbering.
  value_of <- function(level, outcome, term, parameter) {
    printed$estimate[printed$level == level & printed$outcome == outcome &
                       printed$term == term & printed$parameter == parameter]
  }
  renormalise <- function(p) {
    stopifnot("printed probabilities sum to one within their rounding" =
                abs(sum(p) - 1) <= length(p) * 5e-5)
    p / sum(p)
  }
  delta_printed <- vapply(seq_len(g), function(w) {
    value_of("group", sprintf("group_class_%d", w), "", "probability")
  }, numeric(1))
  gamma_printed <- t(vapply(seq_len(g), function(w) {
    vapply(seq_len(k), function(profile) {
      value_of("profile", sprintf("profile_%d", profile), sprintf("group_class_%d", w),
               "probability")
    }, numeric(1))
  }, numeric(k)))
  gamma_printed <- matrix(gamma_printed, g, k)
  natural <- twolevel_lpa_natural(numeric(layout$q), layout)$labels
  measurement <- natural[natural$level == "measurement", , drop = FALSE]
  measured <- vapply(seq_len(nrow(measurement)), function(index) {
    value_of("measurement", measurement$outcome[[index]], measurement$term[[index]],
             measurement$parameter[[index]])
  }, numeric(1))
  delta <- renormalise(delta_printed)
  gamma <- t(apply(gamma_printed, 1L, renormalise))
  theta <- c(measured, as.vector(t(gamma[, -k, drop = FALSE])), delta[-g])
  printed_theta <- c(measured, as.vector(t(gamma_printed[, -k, drop = FALSE])),
                     delta_printed[-g])
  list(parameters = twolevel_lpa_unpack(theta, layout), theta = theta,
       shift = abs(theta - printed_theta), printed = printed)
}

#' The columns of a two-level statistic in a listing block
#' @param header The block's header fields.
#' @param latent `"GClass"` or `"Cluster"`.
#' @param n Number of classes, or 1 for a single column.
#' @return Field positions.
.twolevel_lg_columns <- function(header, latent, n = 1L) {
  at <- which(header == latent)
  stopifnot("the block has one column group for the latent variable" = length(at) == 1L)
  at + seq_len(n) - 1L
}

#' Latent GOLD's modal classification table of one latent variable
#' @param lines The listing's lines.
#' @param latent `"GClass"` or `"Cluster"`.
#' @param n Number of classes.
#' @return A list with `cells` (latent by modal, printed text), `latent_total`
#'   and `modal_total`.
.twolevel_lg_modal_table <- function(lines, latent, n) {
  fields <- .lg_fields(lines)
  at <- which(vapply(fields, function(row) {
    length(row) >= 2L && grepl(sprintf("^%s +Classification Table$", latent), row[[1L]]) &&
      identical(row[[2L]], "Modal")
  }, logical(1)))
  stopifnot("the listing prints one modal classification table per latent variable" =
              length(at) == 1L &&
              identical(fields[[at + 1L]][seq_len(n + 2L)],
                        c("Latent", as.character(seq_len(n)), "Total")))
  rows <- fields[at + 1L + seq_len(n)]
  stopifnot("the rows are the latent classes in order" =
              identical(vapply(rows, `[[`, character(1), 1L), as.character(seq_len(n))))
  total <- fields[[at + 2L + n]]
  stopifnot("the table ends with its totals" = identical(total[[1L]], "Total"))
  list(cells = t(vapply(rows, function(row) row[1L + seq_len(n)], character(n))),
       latent_total = vapply(rows, function(row) row[[n + 2L]], character(1)),
       modal_total = total[1L + seq_len(n)])
}

# ---------------------------------------------------------------------------
# Printed comparisons
# ---------------------------------------------------------------------------

#' Compare printed values with obtained ones, at half a unit in each printed
#' place
#'
#' One row for the whole block, as [compare_block()] gives, with the
#' difference in units of each value's own rounding bound: 1 is the bound.
#'
#' @param quantity Label for the row.
#' @param text The printed values, as text.
#' @param obtained The values the package produced, one per printed value.
#' @param output The `table$column` compared, or `NA`.
#' @param allowance Further allowance per value, added to the bound (zero for
#'   a value computed at the point the program printed it at).
#' @return A one-row comparison `data.frame` with `n_values`.
.twolevel_printed <- function(quantity, text, obtained, output = NA_character_,
                              allowance = 0) {
  stopifnot("one obtained value per printed value" = length(text) == length(obtained),
            "something to compare" = length(text) > 0L)
  reference <- as.numeric(trimws(text))
  bound <- .twolevel_half_units(text) + allowance
  ratio <- abs(reference - obtained) / bound
  ratio[is.na(ratio)] <- Inf
  worst <- which.max(ratio)
  extra <- if (any(allowance > 0)) sprintf(" plus an allowance of up to %.2g", max(allowance)) else ""
  data.frame(quantity = sprintf("%s (%d values, half a unit of the printed place%s)",
                                quantity, length(text), extra),
             reference = reference[[worst]], obtained = obtained[[worst]],
             difference = ratio[[worst]], tolerance = 1, agrees = all(ratio <= 1),
             output = output, n_values = length(text), stringsAsFactors = FALSE)
}

#' The package's value for each printed estimate and standard error
#'
#' Every parameter Latent GOLD prints is matched to the column of the fit's
#' own tables that reports it; a standard error the tables do not report (a
#' class probability's, a covariance's) to `parameter_inference()`.
#'
#' @param tables From [fit_tables()].
#' @param inference The fit's `parameter_inference()` table, `term` filled.
#' @param structure The covariance structure.
#' @return A `data.frame` with the parameter labels, `what` (`estimate` or
#'   `standard_error`), `output` and `value`; a label can have several rows
#'   (an equal variance is reported once per profile).
.twolevel_fit_values <- function(tables, inference, structure) {
  profiles <- tables$profiles
  full <- identical(structure, "full")
  profile_label <- sprintf("profile_%d", profiles$profile)
  variance_label <- data.frame(
    level = "measurement",
    outcome = if (identical(structure, "equal")) "shared" else profile_label,
    term = if (full) paste(profiles$indicator, profiles$indicator, sep = ":") else
      profiles$indicator,
    parameter = if (full) "covariance" else "variance", stringsAsFactors = FALSE)
  mean_label <- data.frame(level = "measurement", outcome = profile_label,
                           term = profiles$indicator, parameter = "mean",
                           stringsAsFactors = FALSE)
  rows <- list(
    data.frame(mean_label, what = "estimate", output = "profiles$mean", value = profiles$mean),
    data.frame(mean_label, what = "standard_error", output = "profiles$mean_standard_error",
               value = profiles$mean_standard_error),
    data.frame(variance_label, what = "estimate", output = "profiles$variance",
               value = profiles$variance),
    data.frame(variance_label, what = "standard_error",
               output = "profiles$variance_standard_error",
               value = profiles$variance_standard_error))
  if (full) {
    covariances <- tables$covariances
    off <- covariances$indicator != covariances$indicator_2
    covariances <- covariances[off, , drop = FALSE]
    later <- ifelse(match(covariances$indicator, unique(profiles$indicator)) >
                      match(covariances$indicator_2, unique(profiles$indicator)),
                    covariances$indicator, covariances$indicator_2)
    earlier <- ifelse(later == covariances$indicator, covariances$indicator_2,
                      covariances$indicator)
    rows <- c(rows, list(data.frame(
      level = "measurement", outcome = sprintf("profile_%d", covariances$profile),
      term = paste(later, earlier, sep = ":"), parameter = "covariance",
      what = "estimate", output = "covariances$covariance", value = covariances$covariance)))
  }
  probabilities <- tables$profile_probabilities
  rows <- c(rows, list(
    data.frame(level = "profile", outcome = sprintf("profile_%d", probabilities$profile),
               term = sprintf("group_class_%d", probabilities$group_class),
               parameter = "probability", what = "estimate",
               output = "profile_probabilities$probability", value = probabilities$probability),
    data.frame(level = "group", outcome = sprintf("group_class_%d", probabilities$group_class),
               term = "", parameter = "probability", what = "estimate",
               output = "profile_probabilities$group_class_probability",
               value = probabilities$group_class_probability)))
  # Standard errors no fit table reports.
  elsewhere <- inference$level != "measurement" |
    (inference$parameter == "covariance" &
       vapply(strsplit(inference$term, ":", fixed = TRUE), function(pair) {
         length(pair) == 2L && !identical(pair[[1L]], pair[[2L]])
       }, logical(1)))
  rows <- c(rows, list(data.frame(
    inference[elsewhere, c("level", "outcome", "term", "parameter")],
    what = "standard_error", output = "parameter_inference$standard_error",
    value = inference$standard_error[elsewhere])))
  out <- do.call(rbind, lapply(rows, function(frame) {
    frame[] <- lapply(frame, function(column) if (is.factor(column)) as.character(column) else column)
    frame
  }))
  row.names(out) <- NULL
  out
}

#' Everything Latent GOLD printed that the package reports, as blocks
#'
#' Each block pairs printed text with the package's values, read from a set of
#' tables shaped like `fit_tables()`: the fit's own, or `reference_tables()`
#' built from scratch at a given point, which is how the stopping allowance is
#' measured (see `.twolevel_lpa_case()`). Values are in a fixed order whatever
#' the row order of the tables.
#'
#' @param label Label of the fit.
#' @param lines The listing's lines.
#' @param posterior_file The posterior file, from `.lg_read_table()`.
#' @param tables The tables to read.
#' @param inference The fit's `parameter_inference()` table, `term` filled.
#' @param structure The covariance structure.
#' @param printed Latent GOLD's printed parameters, in the fit's numbering.
#' @param group_ids,first_row The groups and each one's first data row.
#' @param fit_profile,fit_group The fit's class for each Latent GOLD class.
#' @param lg_profile,lg_group Latent GOLD's class for each of the fit's.
#' @return A list of blocks, each a list with `quantity`, `text`, `values` and
#'   `output`.
.twolevel_printed_blocks <- function(label, lines, posterior_file, tables, inference, structure,
                                     printed, group_ids, first_row, fit_profile, fit_group,
                                     lg_profile, lg_group) {
  block <- function(quantity, text, values, output) {
    stopifnot("one value per printed value" = length(text) == length(values))
    list(quantity = sprintf("%s: Latent GOLD's %s", label, quantity), text = text,
         values = values, output = output)
  }
  key <- function(frame) paste(frame$level, frame$outcome, frame$term, frame$parameter, sep = "|")

  # Estimates and standard errors.
  values <- .twolevel_fit_values(tables, inference, structure)
  values <- values[order(values$what, values$output, key(values), values$value), , drop = FALSE]
  parameters <- lapply(unique(paste(values$what, values$output)), function(pair) {
    rows <- values[paste(values$what, values$output) == pair, , drop = FALSE]
    at <- match(key(rows), key(printed))
    stopifnot("every reported value has a printed counterpart" = !anyNA(at))
    text <- if (identical(rows$what[[1L]], "estimate")) printed$estimate_text else
      printed$standard_error_text
    block(sprintf("printed %s", rows$output[[1L]]), text[at], rows$value, rows$output[[1L]])
  })

  # Likelihood, sizes and criteria.
  criteria <- tables$information_criteria
  criterion_labels <- c(
    "Log-likelihood (LL)" = "log_likelihood", "Number of parameters (Npar)" = "n_parameters",
    "BIC (based on LL)" = "bic_individual", "AIC (based on LL)" = "aic",
    "CAIC (based on LL)" = "caic_individual", "SABIC (based on LL)" = "sabic_individual",
    "BIC (based on LL,Ngroups)" = "bic_groups", "CAIC (based on LL,Ngroups)" = "caic_groups",
    "SABIC (based on LL,Ngroups)" = "sabic_groups", "CLC" = "clc_individual",
    "AWE" = "awe_individual", "ICL-BIC" = "icl_individual")
  scalars <- c(
    lapply(names(criterion_labels), function(name) {
      column <- criterion_labels[[name]]
      block(name, .twolevel_lg_row(lines, name, 2L), criteria[[column]],
            paste0("information_criteria$", column))
    }),
    # Latent GOLD's `Entropy` is the joint classification's: (CLC + 2 LL) / 2.
    list(block("complete-data Entropy, from the package's CLC",
               .twolevel_lg_row(lines, "Entropy", 2L),
               (criteria$clc_individual + 2 * criteria$log_likelihood) / 2,
               "information_criteria$clc_individual"),
         block("number of cases", .twolevel_lg_row(lines, "Number of cases", 2L),
               tables$model$n_observations, "model$n_observations"),
         block("number of groups", .twolevel_lg_row(lines, "Number of groups", 2L),
               tables$model$n_groups, "model$n_groups")))

  # Classification statistics and modal tables of each level.
  statistics_at <- which(vapply(.lg_fields(lines), function(row) {
    length(row) > 0L && identical(row[[1L]], "Classification Statistics")
  }, logical(1)))
  stopifnot("one Classification Statistics block" = length(statistics_at) == 1L)
  statistics_header <- .lg_fields(lines[[statistics_at]])[[1L]]
  overall_fields <- .lg_fields(.lg_section(lines, "ProbMeans-Posterior"))
  levels_of <- list(
    individuals = list(latent = "Cluster", n = length(fit_profile), fit_class = fit_profile),
    groups = list(latent = "GClass", n = length(fit_group), fit_class = fit_group))
  classification <- do.call(c, lapply(names(levels_of), function(level) {
    spec <- levels_of[[level]]
    on_level <- function(table) table[table$level == level, , drop = FALSE]
    classes <- on_level(tables$classification)
    classes <- classes[match(seq_len(spec$n), classes$class), , drop = FALSE]
    counts <- on_level(tables$counts)
    counts <- counts[match(seq_len(spec$n), counts$class), , drop = FALSE]
    entropy <- on_level(tables$entropy)
    averages <- on_level(tables$average_posteriors)
    column <- .twolevel_lg_columns(statistics_header, spec$latent)
    n_units <- sum(classes$n_modal)
    proportions <- counts$effective_proportion
    modal <- .twolevel_lg_modal_table(lines, spec$latent, spec$n)
    cell_grid <- expand.grid(latent = seq_len(spec$n), modal = seq_len(spec$n))
    cell_value <- vapply(seq_len(nrow(cell_grid)), function(index) {
      hit <- averages$assigned_class == spec$fit_class[[cell_grid$modal[[index]]]] &
        averages$class == spec$fit_class[[cell_grid$latent[[index]]]]
      averages$n_assigned[hit] * averages$average_posterior[hit]
    }, numeric(1))
    overall_columns <- .twolevel_lg_columns(overall_fields[[1L]], spec$latent, spec$n)
    stopifnot("ProbMeans-Posterior prints the classes in order" =
                identical(overall_fields[[2L]][overall_columns], as.character(seq_len(spec$n))) &&
                identical(overall_fields[[3L]][[1L]], "Overall"))
    list(
      block(sprintf("%s classification errors", spec$latent),
            .twolevel_lg_row(lines, "Classification errors", column, statistics_at),
            1 - sum(classes$n_modal * classes$average_posterior) / n_units,
            "classification$average_posterior"),
      # Latent GOLD divides by the entropy of the estimated class sizes.
      block(sprintf("%s Entropy R-squared, rebuilt from entropy_sum and the class proportions",
                    spec$latent),
            .twolevel_lg_row(lines, "Entropy R-squared", column, statistics_at),
            1 - entropy$entropy_sum / (n_units * -sum(proportions * log(proportions))),
            "entropy$entropy_sum"),
      block(sprintf("%s modal table cells (posterior mass of each class among each modal class)",
                    spec$latent),
            modal$cells[cbind(cell_grid$latent, cell_grid$modal)], cell_value,
            "average_posteriors$average_posterior"),
      block(sprintf("%s modal table latent totals", spec$latent), modal$latent_total,
            classes$estimated_n[spec$fit_class], "classification$estimated_n"),
      block(sprintf("%s modal table modal totals", spec$latent), modal$modal_total,
            classes$n_modal[spec$fit_class], "classification$n_modal"),
      block(sprintf("%s posterior class proportions (ProbMeans-Posterior)", spec$latent),
            overall_fields[[3L]][overall_columns], proportions[spec$fit_class],
            "counts$effective_proportion"))
  }))

  # Every posterior, at both levels.
  posteriors <- tables$posteriors
  posteriors <- posteriors[order(posteriors$row, posteriors$profile), , drop = FALSE]
  group_posteriors <- tables$group_posteriors
  group_index <- match(as.character(group_posteriors$group), as.character(group_ids))
  group_posteriors <- group_posteriors[order(group_index, group_posteriors$group_class), ,
                                       drop = FALSE]
  group_index <- match(as.character(group_posteriors$group), as.character(group_ids))
  row_text <- as.matrix(posterior_file[sprintf("Cluster#%d", seq_along(fit_profile))])
  group_text <- as.matrix(posterior_file[sprintf("GClass#%d", seq_along(fit_group))])[
    first_row, , drop = FALSE]
  posterior_blocks <- list(
    block("row posteriors", row_text[cbind(posteriors$row, lg_profile[posteriors$profile])],
          posteriors$posterior, "posteriors$posterior"),
    block("group posteriors", group_text[cbind(group_index, lg_group[group_posteriors$group_class])],
          group_posteriors$posterior, "group_posteriors$posterior"))
  c(parameters, scalars, classification, posterior_blocks)
}

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

#' Compare one Latent GOLD case
#' @param case Case name, as the files in `returned/` are named.
#' @param directory The Latent GOLD comparison directory.
#' @param data_file The kit data file Latent GOLD read.
#' @param indicators Indicator names.
#' @param structure `"varying"`, `"equal"` or `"full"`.
#' @param missing The package's `missing` argument.
#' @return A `harness_case`.
.twolevel_lpa_case <- function(case, directory, data_file, indicators, structure, missing) {
  k <- 2L
  g <- 2L
  data <- utils::read.delim(file.path(directory, "kit", data_file), na.strings = ".")
  data <- data[c("id", "g", indicators)]
  x <- as.matrix(data[indicators])
  groups <- data$g
  group_ids <- unique(groups)
  first_row <- match(group_ids, groups)
  layout <- twolevel_lpa_layout(indicators, k, g, structure)
  lines <- .lg_listing_lines(file.path(directory, "returned", paste0(case, ".lst")))
  posterior_file <- .lg_read_table(file.path(directory, "returned",
                                             paste0(case, "_posteriors.txt")))
  posterior_ids <- as.integer(posterior_file[[match("id", names(posterior_file))]])
  stopifnot(
    "Latent GOLD read the kit's rows, in the kit's order" =
      identical(posterior_ids, as.integer(data$id)),
    "Latent GOLD's sample sizes are the kit's" =
      identical(.lg_listing_value(lines, "Number of cases")$value, as.numeric(nrow(data))) &&
      identical(.lg_listing_value(lines, "Number of groups")$value,
                as.numeric(length(group_ids))),
    "Latent GOLD estimated the parameters this layout counts" =
      identical(.lg_listing_value(lines, "Number of parameters (Npar)")$value,
                as.numeric(layout$q)))
  lg <- .twolevel_lg_estimates(lines, layout)
  fit_arguments <- list(data = data[c("g", indicators)], vars = indicators, id = "g",
                        n_profiles = k, n_group_classes = g, missing = missing,
                        variance_model = if (identical(structure, "equal")) "equal" else "varying",
                        covariance_model = if (identical(structure, "full")) "full" else "diagonal")
  source_label <- "Latent GOLD 6.1"

  # ---- Case 1: at Latent GOLD's printed estimates ---------------------------
  label_printed <- sprintf("%s at Latent GOLD's printed estimates", case)
  parameters <- lg$parameters
  variances <- t(vapply(seq_len(k), function(profile) diag(parameters$covariances[, , profile]),
                        numeric(layout$j)))
  start <- starting_values(c(
    list(means = parameters$means, profile_probabilities = parameters$gamma,
         group_probabilities = parameters$delta),
    if (identical(structure, "full")) list(covariances = parameters$covariances) else
      list(variances = matrix(variances, k))))
  at_printed <- do.call(multilpa, c(fit_arguments,
                                    list(n_starts = 1L, max_iter = 0L, start = start)))
  reference_printed <- twolevel_lpa_reference(lg$theta, layout, x, groups, data[indicators])
  exempt_printed <- c(
    `profiles$mean_standard_error` = paste(
      "the package withholds standard errors at a point it did not converge to",
      "(max_iter = 0; 'Inference requires a converged fit'); they are compared in",
      "the converged case"),
    `profiles$variance_standard_error` = paste(
      "the package withholds standard errors at a point it did not converge to",
      "(max_iter = 0); they are compared in the converged case"))
  # Both sides evaluate one formula at the same parameters, so they agree to
  # rounding error; 1e-10 is the claim (observed: 5.5e-14 at worst, 2026-10-27).
  exact <- c(default = 1e-10, `classification$odds_correct_classification` = 1e-10)
  matched_printed <- compare_fit(label_printed, at_printed, reference_printed,
                                 tolerance = exact, exempt = exempt_printed)

  # Latent GOLD computed its log-likelihood at its unrounded estimates, which lie
  # within half a unit (5e-5) of each printed one, plus the renormalisation
  # shift of the probabilities. To second order the log-likelihood moves by at
  # most sum h |g| + h' |H| h / 2 over that box, with g and H the gradient and
  # Hessian of the from-scratch likelihood at the printed estimates.
  refined <- twolevel_lpa_refine(lg$theta, layout, x, groups)
  box <- 5e-5 + lg$shift
  likelihood_reach <- sum(box * abs(refined$start_gradient)) +
    as.numeric(crossprod(box, abs(refined$start_hessian) %*% box)) / 2
  likelihood_text <- .twolevel_lg_row(lines, "Log-likelihood (LL)", 2L)
  printed_at_lg <- .twolevel_printed(
    sprintf("%s: Latent GOLD's log-likelihood, the package scoring its printed estimates (allowance %.2g for their rounding)",
            label_printed, likelihood_reach),
    likelihood_text, as.numeric(logLik(at_printed)),
    output = "information_criteria$log_likelihood", allowance = likelihood_reach)
  comparisons_printed <- rbind(
    transform(matched_printed$comparisons, precision = "machine"),
    transform(printed_at_lg, precision = "printed"))
  case_printed <- harness_case(
    comparisons_printed,
    output_ledger(label_printed, fit_outputs(at_printed), comparisons_printed, exempt_printed))

  # ---- Case 2: converged ------------------------------------------------------
  label_converged <- sprintf("%s converged", case)
  # Newton steps on a numerical gradient resolve the maximum to about 1e-8 in
  # the parameters, which leaves some posteriors Latent GOLD prints to eight or
  # more significant places a hair outside their rounding (1.07 half-units at
  # worst, c03). The reference's own EM, continued to its fixed point (a step
  # below 1e-13), resolves it to rounding error and reproduces every printed
  # posterior (measured 2026-10-27).
  maximum <- twolevel_lpa_fixed_point(refined$theta, layout, x, groups)
  covariances <- twolevel_lpa_covariances(maximum$theta, layout, x, groups)
  reference_max <- twolevel_lpa_reference(maximum$theta, layout, x, groups, data[indicators],
                                          covariances)
  # The package's own search, at tol = 1e-15 (a log-likelihood change below
  # about 3e-12 per iteration): at the default 1e-12 the flat c03 likelihood
  # left its estimates 3e-6 from the maximum, outside the 1e-6 claimed for
  # posteriors; at 1e-15 they are within about 2e-7 (measured 2026-10-27).
  fit <- do.call(multilpa, c(fit_arguments, list(n_starts = 5L, seed = 20261027L,
                                                 tol = 1e-15, max_iter = 200000L)))
  # The reference's class for each of the fit's classes, and back.
  lg_profile <- align_classes(fit$means[, indicators, drop = FALSE], reference_max$means)
  lg_group <- align_classes(unname(fit$profile_probabilities),
                            reference_max$profile_probabilities[, lg_profile, drop = FALSE])
  fit_profile <- match(seq_len(k), lg_profile)
  fit_group <- match(seq_len(g), lg_group)

  # Standard errors: both sides invert an observed information at points at
  # most about 2e-7 apart, the package from its analytic score differenced
  # with step 1e-4 and the reference from a Richardson-extrapolated Hessian;
  # they agree to 5e-8 at worst (c03, measured 2026-10-27). 1e-6 is below
  # 1e-4 of the smallest standard error in any case (0.0196, c12), the
  # relative precision claimed.
  error_tolerance <- c(`profiles$mean_standard_error` = 1e-6,
                       `profiles$variance_standard_error` = 1e-6)
  tables <- fit_tables(fit)
  observed <- parameter_inference(fit, vcov_type = "observed")
  inference_observed <- transform(observed, term = ifelse(is.na(term), "", term))
  key <- function(frame) paste(frame$level, frame$outcome, frame$term, frame$parameter, sep = "|")

  # Where the package stopped. On a flat likelihood its stopping rule (a
  # log-likelihood change below tol * (1 + |LL|)) leaves its estimates up to
  # about 2e-7 from the maximum (c03, measured 2026-10-27), which moves a
  # posterior by more than the 5e-8 to which Latent GOLD prints some of them.
  # That gap is measured, not assumed: every printed quantity is computed from
  # scratch at the package's own estimates and at the refined maximum, and the
  # difference, cell by cell, is added to the half-unit bound. Standard errors
  # are taken at the maximum in both (their change over 2e-7 is far below the
  # printed 5e-5), so their allowance is zero.
  package_covariances <- if (identical(structure, "full")) {
    unname(fit$covariances[indicators, indicators, , drop = FALSE])
  } else {
    vapply(seq_len(k), function(profile) diag(fit$variances[profile, indicators], layout$j),
           matrix(0, layout$j, layout$j))
  }
  theta_package <- twolevel_lpa_pack(list(
    delta = unname(fit$group_probabilities)[fit_group],
    gamma = unname(fit$profile_probabilities)[fit_group, fit_profile, drop = FALSE],
    means = unname(fit$means[, indicators, drop = FALSE])[fit_profile, , drop = FALSE],
    covariances = array(package_covariances, c(layout$j, layout$j, k))[, , fit_profile,
                                                                       drop = FALSE]),
    layout)
  reference_package <- twolevel_lpa_reference(theta_package, layout, x, groups,
                                              data[indicators], covariances)
  block_sets <- lapply(
    list(fit = tables,
         package = reference_tables(align_reference(reference_package, fit), fit),
         maximum = reference_tables(align_reference(reference_max, fit), fit)),
    function(source_tables) {
      .twolevel_printed_blocks(label_converged, lines, posterior_file, source_tables,
                               inference_observed, structure,
                               printed = twolevel_lpa_relabel(lg$printed, fit_profile, fit_group),
                               group_ids = group_ids, first_row = first_row,
                               fit_profile = fit_profile, fit_group = fit_group,
                               lg_profile = lg_profile, lg_group = lg_group)
    })
  printed_rows <- do.call(rbind, Map(function(obtained, at_package, at_maximum) {
    stopifnot("the three block sets line up" =
                identical(obtained$quantity, at_package$quantity) &&
                length(obtained$values) == length(at_maximum$values))
    .twolevel_printed(sprintf("%s, against the package's stopping gap", obtained$quantity),
                      obtained$text, obtained$values, obtained$output,
                      allowance = abs(at_package$values - at_maximum$values))
  }, block_sets$fit, block_sets$package, block_sets$maximum))
  modal_rows <- rbind(
    compare_block(sprintf("%s: Latent GOLD's modal profile of each row", label_converged),
                  as.character(fit_profile[as.integer(posterior_file[["Cluster#"]])]),
                  as.character(tables$assignments$profile), tolerance = 0.5,
                  output = "assignments$profile"),
    compare_block(sprintf("%s: Latent GOLD's modal group class of each row", label_converged),
                  as.character(fit_group[as.integer(posterior_file[["GClass#"]])]),
                  as.character(tables$assignments$group_class), tolerance = 0.5,
                  output = "assignments$group_class"))
  posterior_text <- function(latent, n) {
    as.matrix(posterior_file[sprintf("%s#%d", latent, seq_len(n))])
  }
  row_text <- posterior_text("Cluster", k)
  group_text <- posterior_text("GClass", g)[first_row, , drop = FALSE]
  likelihood_text <- .twolevel_lg_row(lines, "Log-likelihood (LL)", 2L)

  # The reference against Latent GOLD: its maximum, posteriors, within-class
  # row posteriors, complete-data entropy and standard errors.
  evaluated <- twolevel_lpa_evaluate(twolevel_lpa_unpack(maximum$theta, layout), x, groups)
  conditional_text <- do.call(cbind, lapply(seq_len(g), function(w) {
    as.matrix(posterior_file[sprintf("joint2_%d_%d", w, seq_len(k))])
  }))
  group_posterior <- evaluated$group_posteriors
  complete_entropy <- -sum(ifelse(group_posterior > 0, group_posterior * log(group_posterior), 0)) +
    sum(evaluated$row_conditional_entropy)
  reference_errors <- sqrt(diag(covariances$observed))
  error_at <- match(key(lg$printed), key(covariances$labels))
  reference_checks <- rbind(
    .twolevel_printed(sprintf("%s: Latent GOLD's log-likelihood, from-scratch likelihood at Latent GOLD's estimates refined to the maximum",
                              label_converged), likelihood_text, maximum$log_likelihood),
    .twolevel_printed(sprintf("%s: Latent GOLD's row posteriors, from-scratch at the refined maximum",
                              label_converged), row_text, as.vector(evaluated$posteriors)),
    .twolevel_printed(sprintf("%s: Latent GOLD's group posteriors, from-scratch at the refined maximum",
                              label_converged), group_text, as.vector(evaluated$group_posteriors)),
    .twolevel_printed(sprintf("%s: Latent GOLD's row posteriors within each group class (joint2), from-scratch at the refined maximum",
                              label_converged),
                      conditional_text, as.vector(do.call(cbind, evaluated$conditional))),
    .twolevel_printed(sprintf("%s: Latent GOLD's complete-data Entropy, from-scratch group entropy plus row_conditional_entropy",
                              label_converged),
                      .twolevel_lg_row(lines, "Entropy", 2L), complete_entropy),
    .twolevel_printed(sprintf("%s: Latent GOLD's standard errors, from-scratch observed information",
                              label_converged),
                      lg$printed$standard_error_text, reference_errors[error_at]))

  matched_converged <- compare_fit(
    label_converged, fit, reference_max, tolerance = error_tolerance,
    extra = rbind(
      compare_values(sprintf("%s: log-likelihood of Latent GOLD's estimates refined to the maximum",
                             label_converged),
                     maximum$log_likelihood, as.numeric(logLik(fit)), tolerance = 1e-10,
                     scale = "relative", output = "information_criteria$log_likelihood")))
  on_inference <- printed_rows$output == "parameter_inference$standard_error"
  converged_printed <- rbind(printed_rows[!on_inference, , drop = FALSE], modal_rows)
  comparisons_converged <- rbind(
    transform(matched_converged$comparisons, precision = "machine"),
    transform(converged_printed, precision = "printed"),
    transform(reference_checks, precision = "printed"))
  case_converged <- harness_case(
    comparisons_converged,
    output_ledger(label_converged, fit_outputs(fit), comparisons_converged))

  # ---- parameter_inference(), observed and robust, in full ---------------------
  labels_fit <- twolevel_lpa_relabel(covariances$labels, fit_profile, fit_group)
  inference_case <- function(vcov_type, table, extra = NULL) {
    expected <- twolevel_lpa_inference_table(labels_fit,
                                             sqrt(diag(covariances[[vcov_type]])))
    compare_verb(sprintf("%s parameter_inference(vcov_type = \"%s\")", case, vcov_type),
                 list(parameter_inference = keyed_table(
                   table, c("level", "outcome", "term", "parameter"))),
                 list(parameter_inference = expected),
                 # The estimates sit where the fit stopped (up to about 2e-7
                 # from the maximum, c03; worst difference 7.2e-7 in a Wald
                 # statistic); the errors as for the fit's tables above.
                 tolerance = c(default = 1e-6), extra = extra)
  }
  robust <- parameter_inference(fit, vcov_type = "robust")
  # Latent GOLD's printed standard errors that no fit table reports (class
  # probabilities, covariances) are counted against this table.
  case_observed <- inference_case(
    "observed", inference_observed,
    extra = printed_rows[on_inference, , drop = FALSE])
  case_robust <- inference_case("robust", transform(robust, term = ifelse(is.na(term), "", term)))
  case_observed$comparisons$precision <- ifelse(
    grepl("half a unit of the printed place", case_observed$comparisons$quantity, fixed = TRUE),
    "printed", "machine")
  case_robust$comparisons$precision <- "machine"

  combined <- stack_cases(list(case_printed, case_converged, case_observed, case_robust))
  combined$comparisons <- transform(combined$comparisons, source = source_label,
                                    dataset = case)
  combined
}
