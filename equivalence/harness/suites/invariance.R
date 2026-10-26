# Invariance: transformations whose effect on every output is known exactly.
#
# The other suites ask whether this package agrees with another program on the
# data that program was shown. None of them asks whether the answer survives a
# change that should not matter, and every defect of the 2026-10-26 review that
# returned a wrong answer silently was of that kind: a unit change (G1, MR-F3,
# MR-F5), a relabelling (R2, CAT-F3), zero-weight units (ML-1). Each case here
# refits transformed data and compares EVERY output table of the new fit with
# the original fit's tables carried through the same transformation by hand:
#
#   rows permuted         row keys move; nothing else changes
#   groups renamed        group keys change, to labels with `: , [ ]` in them
#   categories renamed    category keys change, likewise
#   indicators renamed    indicator keys change, likewise
#   y -> a + c y          means, variances, covariances and their standard
#                         errors transform; the log-likelihood shifts by
#                         -sum(observed cells x log c); posteriors and every
#                         classification table are unchanged
#   data duplicated       estimates unchanged, standard errors divided by
#                         sqrt(2), likelihood and counts doubled
#   zero-weight units     nothing changes
#   predictor + 2000      slopes and their errors unchanged; the intercept
#                         moves by -2000 x slope, with the delta-method error
#
# The transformed fit's classes are put in the original's order by refitting
# from its own solution relabelled, so the comparison sees the transformed
# search's answer, not a label permutation.

suite_invariance <- function() {
  stack_cases(c(
    .invariance_gaussian_single(),
    .invariance_gaussian_two_level(),
    .invariance_categorical(),
    .invariance_weights(),
    .invariance_regression()))
}

# ---------------------------------------------------------------------------
# Fitting and relabelling
# ---------------------------------------------------------------------------

#' Fit a model of the invariance cases
#' @param spec A list of `multilpa()` arguments other than `data`.
#' @param data The data.
#' @param start Optional start list.
#' @return A fitted model.
.invariance_fit <- function(spec, data, start = NULL) {
  arguments <- c(list(data = data), spec)
  if (!is.null(start)) {
    arguments$start <- start
    arguments$n_starts <- 1L
  }
  do.call(multilpa, arguments)
}

#' The fitted model's own solution as a start list, classes permuted
#' @param fit A fitted model.
#' @param profiles,groups Permutations of the profiles and group classes.
#' @return A list for `multilpa(start = )`.
.invariance_start <- function(fit, profiles, groups) {
  start <- unclass(starting_values(fit))
  if (!is.null(start$means)) start$means <- start$means[profiles, , drop = FALSE]
  if (!is.null(start$variances)) start$variances <- start$variances[profiles, , drop = FALSE]
  if (!is.null(start$covariances)) start$covariances <- start$covariances[, , profiles, drop = FALSE]
  start$profile_probabilities <- start$profile_probabilities[groups, profiles, drop = FALSE]
  start$group_probabilities <- start$group_probabilities[groups]
  if (!is.null(start$response_probabilities)) {
    start$response_probabilities <- lapply(start$response_probabilities,
                                           function(block) block[profiles, , drop = FALSE])
  }
  start
}

#' Refit transformed data and put its classes in the expected order
#'
#' The transformed data are searched from scratch with the original's settings.
#' Its classes are then matched to the expected tables (profiles on the
#' expected means and response probabilities, group classes on the profile mix)
#' and, when the order differs, the model is refitted from its own solution
#' relabelled, which converges where it starts and keeps the order.
#'
#' @param spec The model's arguments.
#' @param data The transformed data.
#' @param expected The expected tables.
#' @return A fitted model in the expected class order.
.invariance_refit <- function(spec, data, expected) {
  searched <- .invariance_fit(spec, data)
  signature <- function(profiles, responses, k) {
    parts <- list()
    if (!is.null(profiles) && nrow(profiles) > 0L) {
      parts$means <- matrix(profiles$mean[order(profiles$indicator, profiles$profile)], k)
    }
    if (!is.null(responses) && nrow(responses) > 0L) {
      ordered <- responses[order(responses$indicator, responses$category, responses$profile), ]
      parts$responses <- matrix(ordered$probability, k)
    }
    do.call(cbind, parts)
  }
  ours <- fit_tables(searched)
  k <- searched$n_profiles
  # Means are compared on the scale of the expected tables, which is the
  # transformed scale, so a rescaled fit aligns on its own units.
  profiles <- align_classes(signature(expected$profiles, expected$responses, k),
                            signature(ours$profiles, ours$responses, k))
  mix <- unname(searched$profile_probabilities)[, profiles, drop = FALSE]
  wanted <- matrix(expected$profile_probabilities$probability,
                   nrow = max(expected$profile_probabilities$group_class))
  groups <- if (nrow(mix) > 1L) align_classes(wanted, mix) else 1L
  if (identical(profiles, seq_len(k)) && identical(groups, seq_len(nrow(mix)))) {
    return(searched)
  }
  .invariance_fit(spec, data, .invariance_start(searched, profiles, groups))
}

#' Compare a transformed fit with its expected tables and keep the ledger
#' @param label Case label.
#' @param fit The transformed, relabelled fit.
#' @param expected Expected tables.
#' @param tolerance Named tolerances.
#' @param exempt Named reasons.
#' @return A `harness_case`.
.invariance_compare <- function(label, fit, expected, tolerance = c(default = 1e-6),
                                exempt = character()) {
  tolerances <- .output_default_tolerances
  tolerances[names(tolerance)] <- tolerance
  obtained <- fit_tables(fit)
  # Search traces, the input echo and the optimiser's bookkeeping are not
  # outputs of the model, and the search on rearranged data may take another
  # path to the same maximum.
  expected <- expected[setdiff(names(expected), c("data", "starts", "stages"))]
  bookkeeping <- grep("[$]", names(.output_exemptions), value = TRUE)
  expected <- lapply(stats::setNames(names(expected), names(expected)), function(name) {
    dropped <- sub("^[^$]*[$]", "", bookkeeping[startsWith(bookkeeping, paste0(name, "$"))])
    expected[[name]][setdiff(names(expected[[name]]), dropped)]
  })
  comparisons <- compare_tables(label, obtained, expected, tolerances)
  comparisons$source <- "invariance"
  harness_case(comparisons, output_ledger(label, fit_outputs(fit), comparisons, exempt))
}

# ---------------------------------------------------------------------------
# Carrying tables through a transformation
# ---------------------------------------------------------------------------

#' Replace key values in every table that carries the key
#' @param tables Named list of tables.
#' @param columns Columns holding the key.
#' @param mapping Named character vector, old value to new value.
#' @param where Optional function of a table returning the rows to change.
#' @return The tables.
.rekey <- function(tables, columns, mapping, where = NULL) {
  lapply(tables, function(table) {
    present <- intersect(columns, names(table))
    if (length(present) == 0L || nrow(table) == 0L) return(table)
    rows <- if (is.null(where)) rep(TRUE, nrow(table)) else where(table)
    table[present] <- lapply(table[present], function(column) {
      text <- as.character(column)
      hit <- rows & text %in% names(mapping)
      text[hit] <- unname(mapping[text[hit]])
      text
    })
    table
  })
}

#' Move row keys after the rows of the data were permuted
#' @param tables Expected tables.
#' @param old_of_new `old_of_new[r]` is the original row now at position `r`.
#' @return The tables.
.expect_rows <- function(tables, old_of_new) {
  mapping <- stats::setNames(as.character(seq_along(old_of_new)), as.character(old_of_new))
  tables <- .rekey(tables, c("row", ".row"), mapping)
  .rekey(tables, "unit", mapping, where = function(table) {
    if (is.null(table$level)) rep(TRUE, nrow(table)) else table$level == "individuals"
  })
}

#' Rename groups
#' @param tables Expected tables.
#' @param mapping Old group label to new.
#' @return The tables.
.expect_groups <- function(tables, mapping) {
  tables <- .rekey(tables, "group", mapping)
  .rekey(tables, "unit", mapping, where = function(table) {
    !is.null(table$level) & table$level == "groups"
  })
}

#' Carry the tables through y -> shift + scale * y, indicator by indicator
#'
#' @param tables Expected tables.
#' @param scale,shift Named by indicator.
#' @param data The original data, to count observed cells.
#' @param group The group column.
#' @return The tables.
.expect_affine <- function(tables, scale, shift, data, group) {
  indicators <- names(scale)
  observed <- !is.na(data[indicators])
  log_scale <- log(abs(scale))
  by_row <- as.vector(observed %*% log_scale)
  delta <- -sum(by_row)
  by_group <- tapply(by_row, as.character(data[[group]]), sum)
  tables$profiles <- within(tables$profiles, {
    c_ <- scale[indicator]
    mean <- shift[indicator] + c_ * mean
    variance <- c_^2 * variance
    standard_deviation <- abs(c_) * standard_deviation
    if (exists("mean_standard_error")) mean_standard_error <- abs(c_) * mean_standard_error
    if (exists("variance_standard_error")) variance_standard_error <- c_^2 * variance_standard_error
    rm(c_)
  })
  if (!is.null(tables$covariances) && nrow(tables$covariances) > 0L) {
    tables$covariances$covariance <- tables$covariances$covariance *
      scale[tables$covariances$indicator] * scale[tables$covariances$indicator_2]
  }
  if (!is.null(tables$group_posteriors$log_likelihood)) {
    tables$group_posteriors$log_likelihood <- tables$group_posteriors$log_likelihood -
      unname(by_group[as.character(tables$group_posteriors$group)])
  }
  criteria <- setdiff(names(tables$information_criteria), "n_parameters")
  tables$information_criteria[criteria] <- lapply(criteria, function(name) {
    value <- tables$information_criteria[[name]]
    if (identical(name, "log_likelihood")) value + delta else value - 2 * delta
  })
  tables$model <- within(tables$model, {
    log_likelihood <- log_likelihood + delta
    aic <- aic - 2 * delta
    bic_groups <- bic_groups - 2 * delta
    bic_individual <- bic_individual - 2 * delta
  })
  tables
}

#' Carry the tables through duplicating every unit under a new label
#'
#' Single-level fits only (each unit is its own group). Estimates, posteriors
#' and every share are unchanged; counts, the likelihood and the observed
#' information double, so standard errors fall by `sqrt(2)`; the criteria are
#' recomputed from their definitions with everything doubled.
#'
#' @param tables Expected tables of the original fit.
#' @param n Number of original rows.
#' @return The tables.
.expect_duplicated <- function(tables, n) {
  copy_rows <- function(table, column) {
    if (is.null(table) || nrow(table) == 0L) return(table)
    copied <- table
    copied[[column]] <- as.character(as.integer(as.character(copied[[column]])) + n)
    table[[column]] <- as.character(table[[column]])
    rbind(table, copied)
  }
  tables$profiles <- within(tables$profiles, {
    mean_standard_error <- mean_standard_error / sqrt(2)
    variance_standard_error <- variance_standard_error / sqrt(2)
  })
  tables$counts$effective_count <- 2 * tables$counts$effective_count
  tables$posteriors <- copy_rows(tables$posteriors, "row")
  tables$posteriors$group <- NULL
  tables$group_posteriors <- copy_rows(tables$group_posteriors, "group")
  tables$assignments <- copy_rows(tables$assignments, ".row")
  tables$classification <- within(tables$classification, {
    n_modal <- 2L * n_modal
    estimated_n <- 2 * estimated_n
  })
  tables$average_posteriors$n_assigned <- 2L * tables$average_posteriors$n_assigned
  individual_bch <- subset(tables$bch_weights, level == "individuals")
  group_bch <- subset(tables$bch_weights, level == "groups")
  tables$bch_weights <- rbind(copy_rows(individual_bch, "unit"), copy_rows(group_bch, "unit"))
  tables$entropy <- within(tables$entropy, {
    n_units <- 2L * n_units
    entropy_sum <- 2 * entropy_sum
  })
  if (!is.null(tables$residuals) && nrow(tables$residuals) > 0L) {
    tables$residuals <- within(tables$residuals, {
      gaussian <- kind == "gaussian"
      statistic[gaussian] <- statistic[gaussian] *
        sqrt(2 * effective_n[gaussian] - 3) / sqrt(effective_n[gaussian] - 3)
      statistic[!gaussian] <- 2 * statistic[!gaussian]
      effective_n <- 2 * effective_n
      p_value <- ifelse(gaussian, 2 * stats::pnorm(-abs(statistic)),
                        stats::pchisq(statistic, df, lower.tail = FALSE))
      p_adjusted <- p_value
      rm(gaussian)
    })
  }
  if (!is.null(tables$gaussian_diagnostics)) {
    tables$gaussian_diagnostics <- within(tables$gaussian_diagnostics, {
      gaussian_mass <- 2 * gaussian_mass
      positive_rows <- 2L * positive_rows
      weight_effective_n <- 2 * weight_effective_n
    })
  }
  ic <- tables$information_criteria
  entropy <- stats::setNames(tables$entropy$entropy_sum, tables$entropy$level)
  # A single-level fit has one group class, so the complete-data entropy both
  # conventions charge is the (doubled) row posterior entropy.
  stopifnot("duplication is carried through for single-level fits only" =
              all(tables$profile_probabilities$group_class == 1L))
  tables$information_criteria <- .reference_criteria(
    2 * ic$log_likelihood, ic$n_parameters, 2 * tables$model$n_informative,
    2 * tables$model$n_groups, entropy[["individuals"]], entropy[["individuals"]])
  tables$model <- within(tables$model, {
    n_observations <- 2L * n_observations
    n_informative <- 2L * n_informative
    n_groups <- 2L * n_groups
    log_likelihood <- 2 * log_likelihood
    aic <- tables$information_criteria$aic
    bic_groups <- tables$information_criteria$bic_groups
    bic_individual <- tables$information_criteria$bic_individual
  })
  tables
}

# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------

#' A label with every character the 2026-10-26 review found mishandled
#' @param text The base text.
#' @return The text decorated with `: , [ ]`.
.awkward <- function(text) sprintf("%s: x, [%s]", text, text)

#' Row permutation, rescaling and duplication on a single-level Gaussian fit
.invariance_gaussian_single <- function() {
  data <- datasets::faithful
  data$unit <- seq_len(nrow(data))
  spec <- list(vars = c("eruptions", "waiting"), id = "unit", n_profiles = 2L,
               n_group_classes = 1L, n_starts = 10L, seed = 1L, tol = 1e-12,
               max_iter = 5000L, covariance_model = "full", variance_model = "varying")
  original <- .invariance_fit(spec, data)
  tables <- fit_tables(original)
  label <- "faithful VVV 2"

  set.seed(20261026L)
  old_of_new <- sample.int(nrow(data))
  permuted <- data[old_of_new, ]
  permuted$unit <- seq_len(nrow(permuted))
  # Unit labels are row numbers, so a moved row's unit label moves with it.
  by_row <- .expect_rows(tables, old_of_new)
  by_row <- .expect_groups(by_row, stats::setNames(as.character(seq_along(old_of_new)),
                                                   as.character(old_of_new)))
  rows_case <- .invariance_compare(sprintf("%s, rows permuted", label),
                                   .invariance_refit(spec, permuted, by_row), by_row)

  affine <- function(scale, shift, tolerance) {
    # The default variance bound (`min_variance = NULL`) is 1e-6 times the
    # smallest indicator variance in the units fitted, so it moves with the
    # data and the suite tests the default. Until review 2026-10-26 round 2
    # (K2-2) it was an absolute 1e-6 and was carried into the new units here by
    # hand, as 1e-6 * min(c^2).
    moved <- data
    moved[names(scale)] <- Map(function(column, c_, a_) a_ + c_ * column,
                               data[names(scale)], scale, shift)
    expected <- .expect_affine(tables, scale, shift, data, "unit")
    .invariance_compare(sprintf("%s, y -> a + c y with c = (%s)", label,
                                paste(format(scale), collapse = ", ")),
                        .invariance_refit(spec, moved, expected), expected,
                        tolerance = tolerance)
  }
  moderate <- affine(c(eruptions = 1e-3, waiting = 1e3), c(eruptions = 50, waiting = -2000),
                     c(default = 1e-6))
  extreme <- affine(c(eruptions = 1e-5, waiting = 1e6), c(eruptions = 0, waiting = 0),
                    c(default = 1e-6))

  doubled <- rbind(data, transform(data, unit = unit + nrow(data)))
  duplicated_case <- .invariance_compare(
    sprintf("%s, every unit duplicated", label),
    .invariance_refit(spec, doubled, .expect_duplicated(tables, nrow(data))),
    .expect_duplicated(tables, nrow(data)))

  list(rows_case, moderate, extreme, duplicated_case)
}

#' Row and group permutation, group renaming and rescaling on a two-level fit
.invariance_gaussian_two_level <- function() {
  set.seed(7L)
  groups <- 40L
  data <- data.frame(school = rep(seq_len(groups), each = 8L))
  type <- 1L + (data$school %% 2L)
  profile <- ifelse(stats::runif(nrow(data)) < c(0.8, 0.25)[type], 1L, 2L)
  data$y1 <- stats::rnorm(nrow(data), c(-1.5, 1.5)[profile])
  data$y2 <- stats::rnorm(nrow(data), c(1, -1)[profile], 0.8)
  data$y3 <- stats::rnorm(nrow(data), c(0, 2)[profile], 1.2)
  data$y2[c(3L, 50L, 101L)] <- NA
  spec <- list(vars = c("y1", "y2", "y3"), id = "school", n_profiles = 2L,
               n_group_classes = 2L, n_starts = 20L, seed = 1L, tol = 1e-12,
               max_iter = 5000L, missing = "fiml")
  original <- .invariance_fit(spec, data)
  tables <- fit_tables(original)
  label <- "two-level VVI 2 x 2"

  set.seed(20261027L)
  old_of_new <- sample.int(nrow(data))
  renamed <- stats::setNames(.awkward(sprintf("school %02d", seq_len(groups))),
                             as.character(seq_len(groups)))
  moved <- data[old_of_new, ]
  moved$school <- unname(renamed[as.character(moved$school)])
  expected <- .expect_groups(.expect_rows(tables, old_of_new), renamed)
  shuffled <- .invariance_compare(sprintf("%s, rows permuted and groups renamed", label),
                                  .invariance_refit(spec, moved, expected), expected)

  scale <- c(y1 = 1e-4, y2 = 1, y3 = 1e4)
  shift <- c(y1 = 3, y2 = 0, y3 = 2000)
  # The default bound is relative to the data, so nothing is carried by hand
  # (review 2026-10-26 round 2, K2-2).
  scaled <- data
  scaled[names(scale)] <- Map(function(column, c_, a_) a_ + c_ * column,
                              data[names(scale)], scale, shift)
  expected_scaled <- .expect_affine(tables, scale, shift, data, "school")
  rescaled <- .invariance_compare(sprintf("%s, y -> a + c y", label),
                                  .invariance_refit(spec, scaled, expected_scaled),
                                  expected_scaled)
  list(shuffled, rescaled)
}

#' Category and indicator renaming on single- and two-level categorical fits
.invariance_categorical <- function() {
  require_suite_packages("poLCA", reason = "the carcinoma data")
  carcinoma <- get(utils::data(list = "carcinoma", package = "poLCA",
                               envir = environment()), envir = environment())
  items <- c("A", "B", "C", "D", "E", "F", "G")
  data <- carcinoma
  data$unit <- seq_len(nrow(data))
  spec <- list(vars = items, id = "unit", n_profiles = 3L, n_group_classes = 1L,
               categorical = items, n_starts = 20L, seed = 1L, tol = 1e-12,
               max_iter = 5000L)
  original <- .invariance_fit(spec, data)
  tables <- fit_tables(original)
  label <- "carcinoma 3 classes"

  categories <- c(`1` = .awkward("no"), `2` = .awkward("yes"))
  recoded <- data
  recoded[items] <- lapply(data[items], function(column) {
    factor(unname(categories[as.character(column)]), levels = unname(categories))
  })
  by_category <- .rekey(tables, "category", categories)
  category_case <- .invariance_compare(sprintf("%s, categories renamed", label),
                                       .invariance_refit(spec, recoded, by_category),
                                       by_category)

  names_new <- stats::setNames(.awkward(sprintf("item %s", items)), items)
  renamed <- data
  names(renamed)[match(items, names(renamed))] <- unname(names_new)
  renamed_spec <- spec
  renamed_spec$vars <- unname(names_new)
  renamed_spec$categorical <- unname(names_new)
  by_item <- .rekey(tables, c("indicator", "indicator_1", "indicator_2"), names_new)
  item_case <- .invariance_compare(sprintf("%s, indicators renamed", label),
                                   .invariance_refit(renamed_spec, renamed, by_item),
                                   by_item)
  list(category_case, item_case)
}

#' Zero-weight units change nothing (review item ML-1)
.invariance_weights <- function() {
  set.seed(11L)
  groups <- 15L
  data <- data.frame(school = rep(seq_len(groups), each = 10L))
  profile <- 1L + (stats::runif(nrow(data)) < c(0.2, 0.8)[1L + data$school %% 2L])
  data$y1 <- stats::rnorm(nrow(data), c(-1.5, 1.5)[profile])
  data$y2 <- stats::rnorm(nrow(data), c(1, -1)[profile])
  data$w <- 1
  spec <- list(vars = c("y1", "y2"), id = "school", n_profiles = 2L,
               n_group_classes = 2L, n_starts = 20L, seed = 1L, tol = 1e-12,
               max_iter = 5000L, weights = "w")
  original <- .invariance_fit(spec, data)
  tables <- fit_tables(original)
  extra <- data.frame(school = rep(groups + seq_len(4L), each = 10L))
  extra$y1 <- stats::rnorm(nrow(extra), 8)
  extra$y2 <- stats::rnorm(nrow(extra), -8)
  extra$w <- 0
  padded <- rbind(data, extra)
  # `n_observations` and `n_groups` count the rows and groups in the data, as
  # they count rows with nothing observed; the sample the model is fitted to
  # (`n_informative`, the group sample size of the criteria) is unchanged.
  expected <- tables
  expected$model$n_observations <- nrow(padded)
  expected$model$n_groups <- length(unique(padded$school))
  list(.invariance_compare("two-level 2 x 2, four zero-weight groups added",
                           .invariance_refit(spec, padded, expected), expected))
}

#' Shifting a predictor's origin (review item MR-F5)
.invariance_regression <- function() {
  data <- study_hours
  spec <- list(formula = score ~ hours, n_classes = 2L, n_starts = 10L, seed = 1L,
               tol = 1e-12, max_iter = 5000L)
  original <- do.call(mixture_regression, c(list(data = data), spec))
  shifted_data <- transform(data, hours = hours + 2000)
  shifted <- do.call(mixture_regression, c(list(data = shifted_data), spec))
  ours <- as.data.frame(get_results(original, "coefficients"))
  theirs <- as.data.frame(get_results(shifted, "coefficients"))
  # Classes by their slopes, which the shift does not move.
  slope <- function(table) {
    rows <- table[table$term == "hours", , drop = FALSE]
    rows$estimate[order(rows$class)]
  }
  order_ours <- order(slope(ours))
  order_theirs <- order(slope(theirs))
  # Rows of one term, in class order, then in the order given.
  pick <- function(table, term, classes) {
    rows <- table[table$term == term, , drop = FALSE]
    rows[order(rows$class), , drop = FALSE][classes, , drop = FALSE]
  }
  covariance <- stats::vcov(original)
  coefficient_names <- names(stats::coef(original))
  intercept_se <- vapply(order_ours, function(class) {
    intercept <- match(sprintf("coefficient.class_%d.(Intercept)", class), coefficient_names)
    hours <- match(sprintf("coefficient.class_%d.hours", class), coefficient_names)
    stopifnot("one intercept and one slope per class" =
                !is.na(intercept) && !is.na(hours))
    gradient <- c(1, -2000)
    block <- covariance[c(intercept, hours), c(intercept, hours)]
    sqrt(drop(t(gradient) %*% block %*% gradient))
  }, numeric(1))
  label <- "mixture regression, hours + 2000"
  expected_intercept <- pick(ours, "(Intercept)", order_ours)$estimate -
    2000 * pick(ours, "hours", order_ours)$estimate
  rows <- rbind(
    compare_values(sprintf("%s: log-likelihood", label), as.numeric(stats::logLik(original)),
                   as.numeric(stats::logLik(shifted)), tolerance = 1e-8, scale = "relative"),
    compare_values(sprintf("%s: slope, class %d", label, seq_along(order_ours)),
                   pick(ours, "hours", order_ours)$estimate,
                   pick(theirs, "hours", order_theirs)$estimate, tolerance = 1e-6,
                   scale = "relative"),
    compare_values(sprintf("%s: slope standard error, class %d", label, seq_along(order_ours)),
                   pick(ours, "hours", order_ours)$std_error,
                   pick(theirs, "hours", order_theirs)$std_error,
                   tolerance = 1e-4, scale = "relative"),
    compare_values(sprintf("%s: intercept, class %d", label, seq_along(order_ours)),
                   expected_intercept, pick(theirs, "(Intercept)", order_theirs)$estimate,
                   tolerance = 1e-6, scale = "relative"),
    compare_values(sprintf("%s: intercept standard error, class %d", label,
                           seq_along(order_ours)),
                   intercept_se,
                   pick(theirs, "(Intercept)", order_theirs)$std_error,
                   tolerance = 1e-4, scale = "relative"))
  rows$source <- "invariance"
  list(harness_case(rows))
}
