# Holding measurement parameters at supplied values.
#
# A staged analysis estimates the measurement model once and then asks what the
# membership structure looks like with that measurement treated as known. Until
# now `start` only chose where EM began, so a supplied measurement solution
# drifted away from the values it was given. `fixed` holds the named blocks at
# those values for the whole fit, which is a different model with fewer free
# parameters, not a different starting point for the same one.

#' Names of the measurement blocks a fit can hold fixed
#' @return Character vector of the accepted `fixed` values, excluding the
#'   `"measurement"` shorthand.
#' @noRd
.multilpa_fixed_blocks <- function() {
  c("means", "variances", "response_probabilities")
}

#' Validate a `fixed` specification against the supplied start
#'
#' @param fixed Character vector naming measurement blocks to hold, possibly
#'   containing the shorthand `"measurement"`, or `character()` for none.
#' @param start Validated start list, or `NULL`.
#' @param covariance_model Diagonal or full residual covariance.
#' @param n_continuous Number of continuous indicators.
#' @param n_categories Category counts per categorical indicator, or `NULL`.
#' @return The expanded character vector of held block names, in a fixed order.
#'   Raises `multilpa_bad_fixed` when the request cannot be met.
#' @noRd
.multilpa_validate_fixed <- function(fixed, start, covariance_model,
                                     n_continuous, n_categories) {
  if (length(fixed) == 0L) return(character(0))
  if (!is.character(fixed) || anyNA(fixed)) {
    stop(errorCondition(
      "`fixed` must be a character vector of measurement block names.",
      class = "multilpa_bad_fixed", call = NULL))
  }
  # The shorthand names every block the model actually has, so that a caller
  # holding "measurement" does not have to know which blocks a mixed-mode or
  # full-covariance fit carries.
  available <- c(if (n_continuous > 0L) c("means", "variances"),
                 if (!is.null(n_categories)) "response_probabilities")
  if ("measurement" %in% fixed) fixed <- c(setdiff(fixed, "measurement"), available)
  fixed <- unique(fixed)
  unknown <- setdiff(fixed, .multilpa_fixed_blocks())
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "`fixed` may name %s or \"measurement\"; received %s.",
      paste(sprintf("\"%s\"", .multilpa_fixed_blocks()), collapse = ", "),
      paste(sprintf("\"%s\"", unknown), collapse = ", ")),
      class = "multilpa_bad_fixed", call = NULL))
  }
  absent <- setdiff(fixed, available)
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf(
      "This model has no %s to hold fixed.",
      paste(sprintf("`%s`", absent), collapse = " or ")),
      class = "multilpa_bad_fixed", call = NULL))
  }
  if (is.null(start)) {
    stop(errorCondition(
      "`fixed` needs `start` to supply the values to hold. Pass starting_values(fit).",
      class = "multilpa_bad_fixed", call = NULL))
  }
  supplied <- vapply(fixed, function(block) {
    if (identical(block, "variances") && identical(covariance_model, "full")) {
      !is.null(start$covariances)
    } else !is.null(start[[block]])
  }, logical(1))
  if (!all(supplied)) {
    stop(errorCondition(sprintf(
      "`start` carries no %s to hold fixed.",
      paste(sprintf("`%s`", fixed[!supplied]), collapse = " or ")),
      class = "multilpa_bad_fixed", call = NULL))
  }
  # Holding everything leaves no measurement parameter free, which is allowed
  # and is exactly the staged case; holding nothing is the ordinary fit.
  intersect(.multilpa_fixed_blocks(), fixed)
}

#' Extract the held values from a validated start
#'
#' @param start Validated start list.
#' @param fixed Expanded character vector of held block names.
#' @param covariance_model Diagonal or full residual covariance.
#' @return A list carrying only the held blocks, named as the maximization step
#'   reads them, or `NULL` when nothing is held.
#' @noRd
.multilpa_held_parameters <- function(start, fixed, covariance_model) {
  if (length(fixed) == 0L) return(NULL)
  held <- list()
  if ("means" %in% fixed) held$means <- start$means
  if ("variances" %in% fixed) {
    # Under a full covariance model the spread is the covariance array; under a
    # diagonal one it is the variance matrix. `fixed = "variances"` names the
    # spread either way, so the caller does not track which is in use.
    if (identical(covariance_model, "full")) {
      held$covariances <- start$covariances
    } else held$variances <- start$variances
  }
  if ("response_probabilities" %in% fixed) {
    held$response_probabilities <- start$response_probabilities
  }
  held
}

#' Apply the held values to a set of initial parameters
#'
#' Every restart must begin from the held values, otherwise a random start
#' would report a measurement solution that was never estimated and never
#' supplied.
#'
#' @param parameters Initial parameters from the initializer or the start.
#' @param held Held values from [.multilpa_held_parameters()], or `NULL`.
#' @return The parameters with every held block replaced.
#' @noRd
.multilpa_apply_held <- function(parameters, held) {
  if (is.null(held)) return(parameters)
  parameters[names(held)] <- held
  # A held covariance array fixes the diagonal it implies, so the variance
  # matrix cannot be left at whatever the initializer produced.
  if (!is.null(held$covariances)) {
    dimension <- dim(held$covariances)[1L]
    parameters$variances <- t(matrix(vapply(
      seq_len(dim(held$covariances)[3L]), function(profile) {
        diag(matrix(held$covariances[, , profile], dimension, dimension))
      }, numeric(dimension)), dimension, dim(held$covariances)[3L]))
  }
  parameters
}

#' Count the measurement parameters a `fixed` specification removes
#'
#' @param fixed Expanded character vector of held block names.
#' @param n_profiles Number of individual profiles.
#' @param n_continuous Number of continuous indicators.
#' @param n_categories Category counts per categorical indicator, or `NULL`.
#' @param variance_model Varying or equal across profiles.
#' @param covariance_model Diagonal or full residual covariance.
#' @return Number of parameters that are no longer free.
#' @noRd
.multilpa_fixed_parameters <- function(fixed, n_profiles, n_continuous,
                                       n_categories, variance_model,
                                       covariance_model) {
  if (length(fixed) == 0L) return(0L)
  covariance_parameters <- if (identical(covariance_model, "full")) {
    n_continuous * (n_continuous + 1) / 2
  } else n_continuous
  as.integer(
    (if ("means" %in% fixed) n_profiles * n_continuous else 0L) +
    (if ("variances" %in% fixed) {
      if (identical(variance_model, "varying")) n_profiles * covariance_parameters
      else covariance_parameters
    } else 0L) +
    (if ("response_probabilities" %in% fixed && !is.null(n_categories)) {
      .multilpa_categorical_parameters(n_profiles, n_categories)
    } else 0L))
}

#' Complete a measurement-only start with free mixing values
#'
#' `starting_values(what = "measurement")` deliberately drops the mixing
#' blocks, because they are sized for the model that produced them and a staged
#' fit usually changes the number of group classes. They are free parameters
#' here, so only a starting point is needed, and a uniform one is used.
#'
#' @param start Start list, possibly without mixing blocks, or `NULL`.
#' @param n_profiles Number of individual profiles.
#' @param n_group_classes Number of latent group classes.
#' @return The start with `profile_probabilities` and `group_probabilities`
#'   present. Blocks the caller supplied are left exactly as given, so a
#'   mismatched one still fails validation rather than being silently replaced.
#' @noRd
.multilpa_complete_start <- function(start, n_profiles, n_group_classes) {
  if (is.null(start)) return(NULL)
  if (is.null(start$profile_probabilities)) {
    start$profile_probabilities <- matrix(1 / n_profiles, n_group_classes,
                                          n_profiles)
  }
  if (is.null(start$group_probabilities)) {
    start$group_probabilities <- rep(1 / n_group_classes, n_group_classes)
  }
  start
}

#' Fit a two-level model in stages, holding the measurement solution
#'
#' Estimates the measurement model once, on its own, and then estimates the
#' group-class structure with that measurement held fixed. This is the staged
#' workflow used when the measurement solution is meant to be decided before
#' the grouping question is asked, so that adding group classes cannot move the
#' profiles they are meant to describe.
#'
#' @inheritParams multilpa
#' @param id Name of the observed group identifier column. Character, factor, or
#'   numeric identifiers are supported; missing identifiers are not. Required
#'   here: this model has a second level by construction, so it has no
#'   single-level form and does not take `multilpa()`'s `id = NULL`.
#' @param n_group_classes Number of latent group classes to estimate in the
#'   second stage. Must be at least two; one would leave nothing to estimate.
#' @param measurement Optional fitted `multilpa` model to use as the first
#'   stage, so that a measurement solution already chosen and inspected is
#'   carried forward rather than refitted. Its profiles, indicators and
#'   measurement options must match the ones requested here. When `NULL`, the
#'   first stage is fitted here with one group class.
#' @param n_starts Number of EM starts, used for both stages.
#' @return A `multilpa` object for the second stage, so that every accessor,
#'   diagnostic and method works on it unchanged. Its measurement parameters
#'   are exactly the first stage's. It additionally carries `fixed`, naming the
#'   held blocks; `n_parameters`, counting only the parameters this stage
#'   estimated; and `n_parameters_with_measurement`, which adds the held
#'   measurement back. Use [get_results()] with `what = "stages"` for a tidy
#'   two-row summary of both stages.
#' @details First-stage uncertainty is **not** propagated. The second stage
#'   treats the measurement solution as known, so its standard errors,
#'   information criteria and likelihood-ratio comparisons are conditional on
#'   that solution and are narrower than they would be if the measurement had
#'   been estimated jointly. This is a property of staging itself, not of this
#'   implementation, and it is the reason both parameter counts are reported:
#'   compare staged fits with one another using `n_parameters`, and compare a
#'   staged fit with a jointly estimated one using
#'   `n_parameters_with_measurement`, remembering that the staged likelihood is
#'   not the joint maximum and the comparison is descriptive.
#' @seealso [multilpa()] with `fixed` for finer control over which blocks are
#'   held, and [starting_values()] with `what = "measurement"` for the values a
#'   stage hands on.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(30), each = 8),
#'   score_a = stats::rnorm(240), score_b = stats::rnorm(240)
#' )
#' staged <- fit_staged(example_data, c("score_a", "score_b"), "school",
#'                      n_profiles = 2, n_group_classes = 2, n_starts = 2,
#'                      seed = 1)
#' get_results(staged, what = "stages")
#' get_results(staged, what = "profile_probabilities")
#' @export
fit_staged <- function(data, vars, id, n_profiles,
                       n_group_classes = 2L, measurement = NULL,
                       variance_model = c("varying", "equal"), n_starts = 10L,
                       max_iter = 1000L, tol = 1e-8, min_variance = 1e-6,
                       seed = NULL, missing = c("error", "fiml"),
                       covariance_model = c("diagonal", "full"),
                       categorical = character(), min_probability = 1e-10,
                       time = NULL) {
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  missing <- match.arg(missing)
  # The data contract itself is checked once, by multilpa(), which both stages
  # call; only what is specific to staging is checked here. `is.numeric` is
  # explicit because `"3" >= 2L` compares as strings and would pass silently.
  stopifnot(
    "`n_group_classes` must be at least two; staging one leaves nothing to estimate" =
      is.numeric(n_group_classes) && length(n_group_classes) == 1L &&
      is.finite(n_group_classes) && n_group_classes >= 2L &&
      n_group_classes == floor(n_group_classes),
    "`measurement` must be NULL or a fitted `multilpa` model" =
      is.null(measurement) || inherits(measurement, "multilpa"))
  ## A staged fit holds the whole measurement, and `multilpa()` refuses to hold
  ## the variances of a structure that is maximized across every profile at
  ## once. This stage cannot name those structures either -- it passes
  ## `variance_model` and `covariance_model`, which reach only four of the
  ## fourteen -- so a constrained measurement used to be fitted as the nearest
  ## structure those two can express: a VEI measurement came back recorded as
  ## VVI, holding the right numbers under the wrong model and counting one
  ## variance parameter too many. It is refused here instead, where the reason
  ## can name the structure the caller actually fitted.
  if (!is.null(measurement)) {
    held_structure <- measurement$covariance_structure %||% NA_character_
    if (!is.na(held_structure) &&
        !held_structure %in% .multilpa_inferable_structures()) {
      stop(errorCondition(sprintf(paste(
        "`measurement` was fitted with the %s covariance structure, which is",
        "maximized across every profile at once, so holding it would not leave",
        "each profile at its maximum. Refit the measurement as EEI, VVI, EEE or",
        "VVV, or fit both levels together with `multilpa()`, which can estimate",
        "%s directly."), held_structure, held_structure),
        class = "multilpa_bad_stage", call = NULL))
    }
  }
  call <- match.call()
  shared <- list(data = data, vars = vars, id = id,
                 n_profiles = n_profiles, variance_model = variance_model,
                 n_starts = n_starts, max_iter = max_iter, tol = tol,
                 min_variance = min_variance, seed = seed, missing = missing,
                 covariance_model = covariance_model, categorical = categorical,
                 min_probability = min_probability, time = time)
  if (is.null(measurement)) {
    measurement <- do.call(multilpa, c(shared, list(n_group_classes = 1L)))
  } else {
    .multilpa_check_measurement(measurement, n_profiles, vars, categorical,
                                covariance_model, variance_model,
                                .multilpa_requested_levels(data, categorical))
  }
  result <- do.call(multilpa, c(shared, list(
    n_group_classes = n_group_classes,
    start = .multilpa_stage_start(measurement),
    fixed = "measurement")))
  result$call <- call
  result$staged <- TRUE
  result$n_parameters_with_measurement <- result$n_parameters +
    .multilpa_fixed_parameters(result$fixed, n_profiles,
                               length(setdiff(vars, categorical)),
                               .multilpa_category_counts(measurement),
                               variance_model, covariance_model)
  result$stage_one <- measurement
  result
}

#' The measurement values a stage hands on, with their labels kept
#'
#' [starting_values()] drops dimension names, because [multilpa()] reads a start
#' by position. A staged fit cannot afford that for the categorical block: the
#' second stage may name `categorical` in a different order, or see the same
#' factor with its levels declared in a different order, and a positional hand
#' over would then attach a first-stage response distribution to the wrong item
#' or the wrong category and report a likelihood for a model nobody asked for.
#' The labelled blocks are handed over instead, so [multilpa()] aligns them by
#' indicator name and category label and refuses an encoding it cannot align.
#'
#' @param measurement The first-stage fitted `multilpa` model.
#' @return A `multilpa_start` carrying the measurement blocks, with the
#'   categorical blocks named by indicator and their columns named by category.
#' @noRd
.multilpa_stage_start <- function(measurement) {
  start <- starting_values(measurement, what = "measurement")
  if (!is.null(measurement$response_probabilities)) {
    start$response_probabilities <- measurement$response_probabilities
  }
  start
}

#' The category labels a requested stage will observe
#'
#' The second stage encodes its own categorical columns, so the labels it will
#' use come from the data it is about to be given, not from the first stage.
#'
#' @param data The data frame the stage will be fitted to.
#' @param categorical Requested categorical indicator names.
#' @return A named list of category labels per indicator, in code order, or
#'   `NULL` when the request names no categorical indicator or names a column
#'   `data` does not have, in which case [multilpa()] reports the broken data
#'   contract in its own words.
#' @noRd
.multilpa_requested_levels <- function(data, categorical) {
  if (length(categorical) == 0L || !is.data.frame(data) ||
      !all(categorical %in% names(data))) {
    return(NULL)
  }
  .multilpa_encode_categorical(data[, categorical, drop = FALSE])$levels
}

#' Category counts of a fitted model, or NULL when it has no categorical block
#' @param object A fitted `multilpa` model.
#' @return Integer vector of category counts, or `NULL`.
#' @noRd
.multilpa_category_counts <- function(object) {
  if (is.null(object$response_probabilities)) return(NULL)
  vapply(object$response_probabilities, ncol, integer(1))
}

#' Check that a supplied first stage matches the second stage being requested
#'
#' A measurement solution fitted to different indicators, a different number of
#' profiles or a different covariance structure cannot be held fixed here; the
#' mismatch is reported rather than reshaped.
#'
#' @param measurement Fitted `multilpa` model offered as the first stage.
#' @param n_profiles Requested number of profiles.
#' @param vars Requested indicator names.
#' @param categorical Requested categorical indicator names.
#' @param covariance_model Requested residual covariance structure.
#' @param variance_model Requested variance constraint.
#' @param levels Named list of the category labels the requested stage will
#'   observe, in code order, or `NULL` when there are no categorical
#'   indicators. The labels, not their order, are what must match: a stage that
#'   sees the same categories in a different order is aligned by label, while a
#'   stage that sees a different category cannot be aligned at all.
#' @return `NULL`, invisibly; raises `multilpa_bad_stage` on the first mismatch.
#' @noRd
.multilpa_check_measurement <- function(measurement, n_profiles, vars,
                                        categorical, covariance_model,
                                        variance_model, levels = NULL) {
  mismatch <- function(what, expected, received) {
    stop(errorCondition(sprintf(
      "`measurement` was fitted with %s %s, but %s was requested.",
      what, received, expected), class = "multilpa_bad_stage", call = NULL))
  }
  if (!identical(as.integer(measurement$n_profiles), as.integer(n_profiles))) {
    mismatch("n_profiles", n_profiles, measurement$n_profiles)
  }
  if (!identical(measurement$vars, vars)) {
    mismatch("vars", paste(vars, collapse = ", "),
             paste(measurement$vars, collapse = ", "))
  }
  if (!identical(sort(measurement$categorical %||% character()), sort(categorical))) {
    mismatch("categorical indicators",
             paste(categorical, collapse = ", ") ,
             paste(measurement$categorical %||% character(), collapse = ", "))
  }
  # The indicators are the same ones; the categories within them must be the
  # same values too, or the held distribution has no meaning on this scale.
  # Order is not part of the requirement, because the hand over is by label.
  if (!is.null(levels)) {
    invisible(lapply(categorical, function(indicator) {
      held <- measurement$categorical_levels[[indicator]]
      wanted <- levels[[indicator]]
      if (is.null(held) || !setequal(held, wanted)) {
        stop(errorCondition(sprintf(
          "`measurement` holds categories %s for `%s`, but this stage observes %s.",
          paste(sprintf("\"%s\"", held %||% character()), collapse = ", "),
          indicator,
          paste(sprintf("\"%s\"", wanted), collapse = ", ")),
          class = "multilpa_bad_stage", call = NULL))
      }
    }))
  }
  if (!identical(measurement$covariance_model %||% "diagonal", covariance_model)) {
    mismatch("covariance_model", covariance_model, measurement$covariance_model)
  }
  if (!identical(measurement$variance_model, variance_model)) {
    mismatch("variance_model", variance_model, measurement$variance_model)
  }
  invisible(NULL)
}
