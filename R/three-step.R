#' Posteriors and modal assignments at a chosen level
#'
#' Only the classified units: rows with at least one observed indicator, and
#' groups with at least one such row, as for entropy, the average posteriors
#' and the information criteria. A row with nothing observed has the class
#' prior for its posterior, so counting it in the classification error matrix
#' charges the larger class with errors no observation made, and carrying it
#' into a three-step model relates its outcome to a class it was never
#' assigned by its data. `drop_empty = FALSE` keeps every unit, as before.
#'
#' @param object A fitted model of this package.
#' @param level `"individuals"` or `"groups"`.
#' @return A list with `posteriors`, `modal`, `units` (row numbers or group
#'   identifiers), `group_index` (the clusters, numbered densely in order of
#'   first appearance), `n_classes` and `kept`, the positions of the kept
#'   units among all rows or all groups, which a caller uses to subset any
#'   per-unit data it pairs with them, and `weights`, the kept units' sampling
#'   weights (`NULL` for an unweighted fit).
#' @noRd
.multilpa_level_assignments <- function(object, level) {
  stopifnot("`object` must be a fitted model of this package" =
              .multilpa_any_fit(object))
  level <- match.arg(level, c("individuals", "groups"))
  if (identical(level, "individuals")) {
    kept <- which(.multilpa_informative_rows(object))
    clusters <- object$group_index[kept]
    list(posteriors = object$subject_posteriors[kept, , drop = FALSE],
         modal = object$subject_profiles[kept], units = kept,
         group_index = match(clusters, unique(clusters)),
         n_classes = object$n_profiles, kept = kept,
         weights = .multilpa_level_weights(object, "individuals", length(kept)))
  } else {
    kept <- which(.multilpa_informative_groups(object))
    list(posteriors = object$group_posteriors[kept, , drop = FALSE],
         modal = object$group_classes[kept], units = object$group_values[kept],
         group_index = seq_along(kept), n_classes = object$n_group_classes,
         kept = kept,
         weights = .multilpa_level_weights(object, "groups", length(kept)))
  }
}

#' Require a first-stage fit whose classification error is not predictor-dependent
#'
#' The profile-level three-step verbs take a covariate-free [multilpa()] fit
#' (class `multilpa`) or an [lta()] fit with homogeneous first-order
#' transitions (class `multilpa_transitions`); trajectory models are checked
#' by `.trajectory_check_three_step()` before this is reached. Every other
#' fitted model of the package is named and refused with
#' `latents_unsupported_three_step`, so a fit of the package is never told it
#' is not one (review 2026-10-27, INF-2); anything else is not a fit and gets
#' `latents_not_a_fit`.
#'
#' @param object The `x` passed to [three_step()] or [r3step()].
#' @return `NULL`, invisibly, or `latents_unsupported_three_step`.
#' @noRd
.multilpa_check_three_step_fit <- function(object) {
  unsupported <- function(...) {
    stop(errorCondition(paste(...), class = "latents_unsupported_three_step",
                        call = NULL))
  }
  if (inherits(object, "multilpa_covariates")) {
    unsupported(
      "Three-step correction is unavailable for a fit that already uses",
      "membership covariates: this method has no classification-error",
      "adjustment conditional on those covariates. Fit the measurement model",
      "without membership covariates before using three_step() or r3step().")
  }
  if (inherits(object, "multilpa_lta")) {
    arguments <- object$arguments
    if (length(arguments$initial_covariates) > 0L ||
        length(arguments$transition_covariates) > 0L || isTRUE(object$mover_stayer)) {
      unsupported(
        "The step-one transition model has covariates (or a stayer class): its",
        "classification errors are conditional on them, which the fixed error",
        "matrix of the three-step correction does not allow. Fit the step-one",
        "model without them.")
    }
    unsupported(
      "three_step() and r3step(covariates = ) are implemented for the transition",
      "model lta() fits with its default homogeneous first-order transitions,",
      "invariant measurement and continuous or categorical indicators; this fit",
      "comes from the general transition engine (occasion-varying or",
      "second-order transitions, occasion-specific measurement, or ordinal or",
      "count indicators). Use r3step(initial_covariates = ,",
      "transition_covariates = ) for covariates of the initial profile and the",
      "transitions, or fit the default transition model for a profile-level",
      "analysis.")
  }
  if (inherits(object, "multilpa_additive")) {
    unsupported(
      "three_step() and r3step() are not implemented for a group-class family",
      "without individual profiles (multilpa(family = \"additive\"),",
      "\"dispersion\" or \"additive_dispersion\"). Fit the model with",
      "`family = \"profiles\"` for a three-step analysis.")
  }
  if (inherits(object, "multilpa_cross_level")) {
    unsupported(
      "three_step() and r3step() are not implemented for a cross-level model",
      "(multilpa(family = \"restricted_cross_level\") or",
      "\"full_cross_level\"). Fit the model with `family = \"profiles\"` for",
      "a three-step analysis.")
  }
  if (inherits(object, "latents_mixture_regression")) {
    unsupported(
      "three_step() and r3step() relate an outcome or covariates to person",
      "classes, and this mixture regression classifies observations. Fit it",
      "with `id` and `class_level = \"group\"` for a three-step analysis.")
  }
  .latents_require_fit(object, accepted = c("multilpa", "multilpa_transitions"),
                       expected = paste("a covariate-free fit from multilpa() or",
                                        "lta(), or a trajectory model from",
                                        "mixture_regression()"))
  invisible(NULL)
}

#' How often each true class is assigned to another
#'
#' Documented on `?get_results`, which is where a caller reaches this table from,
#' and consumed by [three_step()] and the BCH weights.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A base `data.frame`, one row per level and ordered pair of classes,
#'   whose probabilities sum to one within each `true_class`.
#' @noRd
.multilpa_classification_errors <- function(x,
                                            level = c("individuals", "groups",
                                                      "both")) {
  level <- match.arg(level)
  result <- do.call(rbind, lapply(.multilpa_classification_levels(x, level),
                                  .multilpa_error_frame, x = x))
  row.names(result) <- NULL
  result
}

#' The classification error table for one level
#' @param x A fitted model of this package.
#' @param level A single level name.
#' @return One row per ordered pair of classes at that level.
#' @noRd
.multilpa_error_frame <- function(x, level) {
  # A noise unit's modal class is 0 and its posteriors lack the noise column,
  # so the error matrix would not be a matrix of probabilities.
  .multilpa_refuse_noise(x, "The classification-error table")
  pieces <- .multilpa_level_assignments(x, level)
  matrix_form <- .multilpa_error_matrix(pieces)
  classes <- seq_len(pieces$n_classes)
  index <- expand.grid(assigned_class = classes, true_class = classes)
  data.frame(level = level, true_class = index$true_class,
             assigned_class = index$assigned_class,
             probability = matrix_form[cbind(index$true_class, index$assigned_class)],
             row.names = NULL)
}

#' The classification error matrix, rows indexed by the true class
#'
#' `D[t, s] = sum_i w_i p_it 1(W_i = s) / sum_i w_i p_it`, with
#' `w_i` the unit's sampling weight (1 without weights), so under weights it
#' is the population's error matrix (review 2026-10-26 round 2, S-4).
#'
#' @param pieces The list from `.multilpa_level_assignments()`.
#' @return A square matrix whose rows sum to one.
#' @noRd
.multilpa_error_matrix <- function(pieces) {
  classes <- seq_len(pieces$n_classes)
  weighted <- pieces$posteriors * (pieces$weights %||% 1)
  totals <- colSums(weighted)
  ## vapply() drops to a vector when there is only one class, so the shape is
  ## forced; a one-class model is degenerate but must not error here.
  result <- matrix(vapply(classes, function(assigned) {
    colSums(weighted * (pieces$modal == assigned)) / totals
  }, numeric(pieces$n_classes)), nrow = pieces$n_classes,
  dimnames = list(paste0("true_", classes), paste0("assigned_", classes)))
  result
}

#' The BCH inverse-classification-error weights
#'
#' Documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @param level `"individuals"`, `"groups"` or `"both"`.
#' @return A base `data.frame`, one row per level, unit and class.
#' @noRd
.multilpa_bch_weights <- function(x, level = c("individuals", "groups",
                                               "both")) {
  level <- match.arg(level)
  result <- do.call(rbind, lapply(.multilpa_classification_levels(x, level),
                                  .multilpa_weight_frame, x = x))
  row.names(result) <- NULL
  result
}

#' The BCH weights for one level
#' @param x A fitted model of this package.
#' @param level A single level name.
#' @return One row per unit and class at that level.
#' @noRd
.multilpa_weight_frame <- function(x, level) {
  .multilpa_refuse_noise(x, "The BCH weights table")
  pieces <- .multilpa_level_assignments(x, level)
  inverse <- .multilpa_invert_errors(.multilpa_error_matrix(pieces))
  weights <- inverse[pieces$modal, , drop = FALSE]
  classes <- seq_len(pieces$n_classes)
  ## t() first, so that transposing the matrix into a vector walks a unit's
  ## classes before moving to the next unit and the frame reads unit by unit.
  data.frame(level = level,
             unit = rep(pieces$units, each = pieces$n_classes),
             assigned_class = rep(pieces$modal, each = pieces$n_classes),
             class = rep(classes, times = length(pieces$modal)),
             weight = as.vector(t(weights)),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Invert a classification error matrix, refusing a singular one
#' @param errors The matrix from `.multilpa_error_matrix()`.
#' @return Its inverse.
#' @noRd
.multilpa_invert_errors <- function(errors) {
  inverse <- tryCatch(solve(errors), error = function(condition) {
    stop(errorCondition(
      paste("The classification error matrix cannot be inverted, so the classes",
            "are not separated well enough for a three-step correction.",
            "Inspect get_results(x, \"classification_errors\")."),
      class = "latents_inseparable_classes", call = NULL))
  })
  if (min(abs(eigen(as.matrix(errors), only.values = TRUE)$values)) < 1e-8) {
    stop(errorCondition(
      "The classification error matrix is numerically singular; the classes are not separable.",
      class = "latents_inseparable_classes", call = NULL))
  }
  inverse
}

#' Compare latent classes on an external outcome
#'
#' @description
#' `three_step()` estimates how an outcome that was not used to define the
#' classes differs between them, holding the measurement model fixed so the
#' classes keep their meaning. Comparing classes by modal assignment treats
#' uncertain assignments as certain and biases the differences; the
#' three-step methods correct for the classification error of the fitted
#' model. The default is the Bolck-Croon-Hagenaars (`"bch"`) weighting; the
#' maximum-likelihood (`"ml"`) method models the outcome given the true
#' class instead. Both assume that the outcome and the assigned class are
#' independent given the true class; the differences are associations, not
#' effects of class membership.
#'
#' @param x A covariate-free fit: from [multilpa()] (the profile family),
#'   from [lta()] with its default homogeneous first-order transitions and
#'   invariant measurement, or a trajectory model that classifies persons
#'   ([mixture_regression()] with `class_level = "group"`, with or without
#'   `random`) fitted without `membership` covariates. Refused: fits with
#'   membership covariates, multilevel growth mixtures (`cluster`), the
#'   other `multilpa()` families, observation-level mixture regressions,
#'   fits with weights or a noise component, and `lta()` fits with
#'   occasion-specific or covariate transitions, second order, occasion-
#'   specific measurement, a stayer class or ordinal or count indicators.
#' @param data The data frame carrying the outcome, in the fit's row order;
#'   columns shared with the fit are checked row by row. For a trajectory
#'   model, the data the model was fitted to.
#' @param outcome Name of the outcome column, without missing or infinite
#'   values and not an indicator of `x`: numeric for class means, or a
#'   factor, character or logical column for the probability of each
#'   category in each class. For
#'   `level = "groups"` it must be constant within each group; for a
#'   trajectory model, within each person.
#' @param level `"individuals"` (default) relates the outcome to the
#'   profiles, `"groups"` to the group classes. A trajectory model relates
#'   it to its person classes, reported as `level = "persons"`.
#' @param method `"bch"` (default): Bolck-Croon-Hagenaars weights. `"ml"`:
#'   the maximum-likelihood step-three model. `"proportional"`: posterior
#'   weights without correction. `"modal"`: modal assignment without
#'   correction. See *Methods*.
#' @param ci_level Confidence level of the intervals (default 0.95).
#' @param contrast `"none"` (default): one outcome mean per class (one
#'   probability per class and category). `"pairs"`: the difference between
#'   every pair of classes, tested against zero. `"overall"`: one Wald test
#'   that every class has the same mean (the same category distribution).
#' @param adjust Multiplicity correction across the pairwise differences,
#'   a `stats::p.adjust()` method; `"BH"` by default. Ignored for
#'   `contrast = "none"` or `"overall"`. With a binary outcome, its two
#'   complementary category rows count as one test per class pair; with
#'   more categories, every reported category difference enters the family.
#' @param vcov_type `"cluster"` (default) treats the fit's groups as the
#'   independent units of the sandwich variance, as the nesting requires;
#'   `"independent"` treats every observation as independent.
#' @param variance For `method = "ml"` with a numeric outcome (refused
#'   otherwise): `"varying"` (default) estimates an outcome variance in each
#'   class, `"equal"` one shared by the classes.
#' @param step_one `"fixed"` (default) conditions on the step-one estimates.
#'   `"estimated"` adds the uncertainty of the classification-error
#'   probabilities (and, under `"ml"`, of the class sizes) to the step-three
#'   covariance; available for [lpa()], [lca()] and [multilpa()] profile
#'   fits at either level when their Wald inference is available and their
#'   measurement was estimated jointly. Staged fits and fits holding a
#'   measurement block fixed are refused. See
#'   *Step-one uncertainty*.
#'
#' @return A data frame. For `contrast = "none"`, one row per class (per
#'   class and category for a categorical outcome, with `category` after
#'   `class`): `level`, `method`, `class`, `estimate`, `standard_error`,
#'   `conf_low`, `conf_high` and `effective_n`; means are not tested, and a
#'   probability's interval is formed on the logit scale. Under
#'   `method = "ml"` with a numeric outcome a `variance` column precedes
#'   `effective_n`, and under `"ml"` `effective_n` is the expected number of
#'   units in the class. For `contrast =
#'   "pairs"`, one row per pair (per pair and category): `level`, `method`,
#'   `class`, `reference_class`, `estimate` (class minus reference),
#'   `standard_error`, `statistic`, `p_value`, `p_adjusted`, `conf_low` and
#'   `conf_high`, as in [r3step()]. For `contrast = "overall"`, one row:
#'   `level`, `method`, `statistic`, `df` (`K - 1` for a mean,
#'   `(K - 1)(C - 1)` for `C` categories) and `p_value`. The attributes
#'   `vcov_type` and `step_one` record the choices, `adjust` on a pairs
#'   result and `variance` on a numeric ML result with `contrast = "none"`.
#'
#' @section Methods:
#' BCH gives each observation one weight per class, the row of the inverse
#' classification-error matrix for its assigned class, and takes weighted
#' means. With `method = "ml"` the assigned
#' class is an
#' indicator of the true class with the step-one error probabilities, and
#' the outcome is normal within each true class (a multinomial distribution
#' for a categorical outcome), with the class sizes fixed at the average
#' posteriors, as in Latent GOLD's step-three model. The likelihood is
#' maximised from a fixed set of starts plus twenty seeded random ones, so
#' the result does not depend on the caller's seed; the global maximum is
#' not guaranteed.
#' ML depends on the within-class outcome distribution: skewness or a wrong
#' equal-variance constraint can bias the means, and BCH is less sensitive to
#' that.
#'
#' A categorical outcome is described by its probability in each class. BCH
#' weights can be negative, so with poorly separated classes a share can
#' fall outside 0 to 1; it is then reported without an interval or test,
#' with a warning. The ML probabilities stay in range; where the BCH shares
#' do not, the ML maximum usually sits on the bound and is reported without
#' standard errors, with a warning.
#'
#' On an [lta()] fit, `level = "individuals"` uses one row per occasion and
#' relates the outcome to the profile at that occasion, with the
#' group-clustered variance accounting for the repeated occasions. It is a
#' marginal profile-outcome analysis, not a model of outcome effects on
#' transitions.
#'
#' @section Standard errors:
#' The variance is a sandwich over the `G` independent units (the groups,
#' or the observations under `vcov_type = "independent"`) with the
#' `G / (G - 1)` small-sample factor. It needs more units than classes;
#' with as many or fewer, the analysis is refused. A class carried by one unit (Kish's
#' effective number of units carrying it below 1.5) has no standard error:
#' its interval and the contrasts involving it are `NA`, with a warning.
#'
#' @section Step-one uncertainty:
#' With `step_one = "fixed"` the intervals ignore the uncertainty of the
#' classification-error probabilities, which lowers their coverage when
#' classification is uncertain. `"estimated"` adds the first-order
#' contribution of the step-one covariance, with the data and modal
#' assignments held fixed; the estimates
#' are unchanged, the standard errors and intervals grow, and the tests use
#' the corrected covariance. The first-stage contribution uses its
#' model-based observed information, even when the step-three sandwich is
#' clustered; it does not add a robust first-stage covariance for residual
#' dependence outside the fitted model. The correction is a large-sample one: in a
#' small, poorly separated sample the intervals still cover less than the
#' nominal level. The modal method is unaffected, as it uses no posteriors.
#' A staged fit or a fit with a held measurement block does not carry the
#' measurement covariance needed by this correction; use `step_one = "fixed"`
#' to report inference conditional on that measurement.
#' Neither option includes the uncertainty of having selected the class
#' solution.
#'
#' @seealso [r3step()] for covariates of class membership; [fit_staged()].
#' @examples
#' set.seed(11)
#' group <- rep(seq_len(40), each = 10)
#' truth <- 1L + as.integer(runif(400) > 0.5)
#' example_data <- data.frame(
#'   g = group,
#'   a = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
#'   b = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
#'   y = rnorm(400, ifelse(truth == 2L, 10, 0))
#' )
#' fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 4, seed = 1)
#' three_step(fit, example_data, "y")
#' three_step(fit, example_data, "y", contrast = "pairs")
#' # A categorical outcome: its probability in each class.
#' example_data <- transform(example_data, passed = y + rnorm(400, sd = 6) > 5)
#' three_step(fit, example_data, "passed")
#' three_step(fit, example_data, "passed", contrast = "overall")
#' @export
three_step <- function(x, data, outcome,
                       level = c("individuals", "groups"),
                       method = c("bch", "proportional", "modal", "ml"),
                       ci_level = 0.95, contrast = c("none", "pairs", "overall"),
                       adjust = c("BH", "holm", "hochberg", "hommel",
                                    "bonferroni", "BY", "none"),
                       vcov_type = c("cluster", "independent"),
                       variance = c("varying", "equal"),
                       step_one = c("fixed", "estimated")) {
  # Whether `variance` was given is read before it is matched: missing() is
  # FALSE once an argument has been reassigned.
  variance_given <- !missing(variance)
  step_one <- .latents_match_formal(step_one)
  finish <- function(result) {
    attr(result, "step_one") <- step_one
    result
  }
  if (identical(step_one, "estimated") &&
      (.trajectory_is_fit(x) || inherits(x, c("multilpa_transitions", "multilpa_lta")))) {
    stop(errorCondition(paste(
      "`step_one = \"estimated\"` is unavailable for transition and trajectory fits;",
      "use `step_one = \"fixed\"`, which still works for supported three-step fits."),
      class = "latents_unsupported_three_step", call = NULL))
  }
  level <- match.arg(level)
  # The released methods keep their abbreviations ("m" is "modal", as it was
  # before "ml" existed); anything else is matched against all four, so "ml"
  # is taken by its full name and an unknown method's error lists every one.
  methods <- eval(formals(three_step)$method)
  released <- pmatch(method[1L], c("bch", "proportional", "modal"))
  method <- if (identical(method, methods)) "bch" else
    if (length(method) == 1L && !is.na(released)) methods[released] else
      match.arg(method, methods)
  variance <- match.arg(variance)
  if (variance_given && !identical(method, "ml")) {
    stop(errorCondition(paste(
      "`variance` sets the class variances of the outcome in the ML step-three",
      "model; the other methods weight the outcome and model no variance. Use",
      "`method = \"ml\"`."), class = "latents_bad_argument", call = NULL))
  }
  contrast <- match.arg(contrast)
  adjust <- match.arg(adjust)
  vcov_type <- match.arg(vcov_type)
  .multilpa_refuse_noise(x, "three_step()")
  if (.latents_is_weighted(x)) .latents_refuse_weights("`three_step()` yet")
  .latents_check_data_frame(data, substitute(data))
  if (is.character(outcome) && length(outcome) == 1L && !outcome %in% names(data)) {
    stop(errorCondition(sprintf("`outcome` names `%s`, which `data` does not have.%s",
                                outcome, .latents_did_you_mean(outcome, names(data))),
                        class = "latents_bad_data", call = NULL))
  }
  .latents_require(
    "`data` must be a data frame" = is.data.frame(data),
    "`outcome` must name a single column of `data`" =
      is.character(outcome) && length(outcome) == 1L && outcome %in% names(data),
    "`outcome` must be numeric, a factor, character or logical" =
      is.numeric(data[[outcome]]) || is.factor(data[[outcome]]) ||
      is.character(data[[outcome]]) || is.logical(data[[outcome]]),
    "`outcome` must not be missing" = !anyNA(data[[outcome]]),
    "`ci_level` must be a single number in (0, 1)" =
      is.numeric(ci_level) && length(ci_level) == 1L && ci_level > 0 && ci_level < 1
  )
  categorical <- !is.numeric(data[[outcome]])
  if (categorical) {
    if (variance_given) {
      stop(errorCondition(paste(
        "`variance` sets the class variances of a numeric outcome; a categorical",
        "outcome has a probability for each category instead."),
        class = "latents_bad_argument", call = NULL))
    }
    # Categories in the factor's level order (unused levels dropped), or
    # sorted the same way in every locale. They code the outcome here; the
    # categories analysed are those of the units that enter step three.
    observed <- data[[outcome]]
    categories <- if (is.factor(observed)) levels(droplevels(observed)) else
      sort(unique(as.character(observed)), method = "radix")
    data[[outcome]] <- match(as.character(observed), categories)
  } else {
    .latents_require("`outcome` must be finite" = all(is.finite(data[[outcome]])))
  }
  if (.trajectory_is_fit(x)) {
    # Trajectory models classify persons: one outcome value per person.
    if (identical(level, "groups")) {
      stop(errorCondition(paste(
        "A trajectory model classifies persons and has no group classes, so",
        "`level = \"groups\"` does not apply; its outcome is related to the",
        "person classes."), class = "latents_bad_argument", call = NULL))
    }
    prepared <- .trajectory_three_step_inputs(x, data, outcome)
    pieces <- prepared$pieces
    values <- prepared$values
    level <- "persons"
  } else {
    .multilpa_check_three_step_fit(x)
    if (outcome %in% .multilpa_measured_items(x)) {
      stop(errorCondition(sprintf(
        "`%s` helped define the fitted classes, so it is not a distal outcome.",
        outcome), class = c("multilpa_indicator_reused", "latents_bad_outcome"),
        call = NULL))
    }
    pieces <- .multilpa_level_assignments(x, level)
    .multilpa_check_row_count(x, data)
    .multilpa_check_alignment(x, data)
    # Only the classified units enter step three, as they do the error matrix.
    values <- .multilpa_outcome_values(x, data, outcome, level)[pieces$kept]
  }
  if (categorical) {
    # A category carried only by units left out of step three (unclassified
    # rows, wholly excluded groups) is not part of the analysed outcome: kept,
    # it would be an empty column whose zero probability puts the ML model on
    # its bound and withholds every standard error. The codes follow the
    # category order, so sorting them keeps that order.
    present <- sort(unique(values))
    categories <- categories[present]
    values <- match(values, present)
    if (length(categories) < 2L) {
      stop(errorCondition(sprintf(paste(
        "`%s` takes one category among the units entering step three, so it",
        "has no distribution to compare across classes."), outcome),
        class = "latents_bad_outcome", call = NULL))
    }
  }
  correction <- if (identical(step_one, "estimated")) {
    .multilpa_step_one_correction(x, pieces, level, method)
  }
  if (categorical) {
    return(finish(.multilpa_three_step_categorical(
      pieces, as.integer(values), categories, level, method, ci_level, contrast,
      adjust, vcov_type, correction)))
  }
  if (identical(method, "ml")) {
    return(finish(.multilpa_three_step_ml(pieces, values, level, variance, ci_level,
                                         contrast, adjust, vcov_type, correction)))
  }
  weights <- .multilpa_step_weights(pieces, method)
  classes <- seq_len(pieces$n_classes)
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)

  totals <- colSums(weights)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop(errorCondition("A class has no positive total outcome weight.",
                        class = "latents_inseparable_classes", call = NULL))
  }
  # BCH weights may be negative. Multiplying and summing a constant outcome
  # directly can then leave a few ulps of cancellation, which looks like a
  # nonzero class contrast with an equally tiny standard error and can produce
  # a false significant p-value. Work on differences from an observed value.
  baseline <- values[[1L]]
  centered <- values - baseline
  centered_estimates <- as.vector(crossprod(weights, centered)) / totals
  estimates <- baseline + centered_estimates
  ## Cluster-robust variance of a weighted mean: the per-unit influence
  ## contributions are summed within each independent group, then across
  ## groups. Holding them as a matrix lets a difference between two classes be
  ## given a standard error from the same quantities, rather than assuming the
  ## two means are independent when they share every unit.
  influence <- sweep(weights * outer(centered, centered_estimates, "-"),
                     2L, totals, "/")
  clustered <- .multilpa_step_clusters(influence, pieces, vcov_type)
  # A class whose weight sits in one independent unit has a summed influence
  # of zero in that unit by the definition of its mean, so its sandwich is 0
  # by construction, not an estimate; the ML path already says so (INF-1).
  clustered <- .multilpa_step_thin_classes(clustered, weights, pieces, vcov_type)
  if (!is.null(correction)) {
    clustered <- correction(clustered, function(updated) {
      weights <- .multilpa_step_weights(updated, method)
      as.vector(crossprod(weights, centered)) / colSums(weights)
    })
  }

  if (identical(contrast, "overall")) {
    return(finish(.multilpa_step_overall(centered_estimates, clustered, pieces$n_classes,
                                        level, method, vcov_type)))
  }
  if (identical(contrast, "pairs")) {
    return(finish(.multilpa_step_pairs(centered_estimates, clustered, level, method, quantile,
                                      adjust, vcov_type)))
  }
  errors <- sqrt(colSums(clustered^2))
  result <- data.frame(level = level, method = method, class = classes,
                       estimate = estimates, standard_error = errors,
                       conf_low = estimates - quantile * errors,
                       conf_high = estimates + quantile * errors,
                       effective_n = totals^2 / colSums(weights^2),
                       row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "vcov_type") <- vcov_type
  column <- .latents_column
  finish(.latents_table(result, ci_level,
                        sprintf("Distal outcome by class (%s, %s%% CI)", method,
                                format(100 * ci_level)),
                        list(column("Class", "class", "integer"),
                             column("Mean", "estimate"),
                             column(sprintf("%s%% CI", format(100 * ci_level)),
                                    c("conf_low", "conf_high"), "ci", "estimate"),
                             column("SE", "standard_error"),
                             column("Effective n", "effective_n"))))
}

#' Propagate step-one uncertainty through fixed modal assignments
#'
#' The step-one data are the fit's own, as for [vcov()]: the `data` of
#' [three_step()] need carry only the outcome.
#' @param x,pieces,level,method The fit and analysed units of [three_step()].
#' @return A function appending covariance pseudo-rows to scaled contributions.
#'   Its estimator returns step-three estimates, or scores with a supplied
#'   implicit derivative transform.
#' @noRd
.multilpa_step_one_correction <- function(x, pieces, level, method) {
  stopifnot("`x` must be a multilpa fit" = inherits(x, "multilpa"),
            "`pieces` must describe the analysed units" = is.list(pieces))
  if (isTRUE(x$staged) || length(x$fixed %||% character()) > 0L) {
    stop(errorCondition(paste(
      "`step_one = \"estimated\"` needs the covariance of the estimated",
      "measurement parameters; a staged fit or one holding measurement",
      "blocks fixed does not carry that covariance. Use",
      "`step_one = \"fixed\"` for inference conditional on the measurement."),
      class = "latents_unsupported_three_step", call = NULL))
  }
  if (identical(method, "modal")) return(NULL)
  covariance <- tryCatch(stats::vcov(x, scale = "unconstrained"),
                         error = function(condition) {
    condition$message <- paste(conditionMessage(condition),
                                '`step_one = "fixed"` still works without step-one inference.')
    stop(condition)
  })
  root <- .multilpa_covariance_root(covariance)
  prepared <- .multilpa_inference_matrix(x, .multilpa_resolve_data(x, NULL))
  centered <- x
  centered$means <- sweep(x$means, 2L, prepared$centers, "-")
  full <- .multilpa_coefficients(centered, "unconstrained")
  free <- .multilpa_free_index(x, "unconstrained")
  restore <- .multilpa_restore_held(full, free)
  theta <- unname(full[free])
  # Each coordinate is stepped on its own scale, as the Wald information is,
  # so the correction does not depend on the units of the data.
  direct <- if (is.null(x$extra_data$direct)) NULL else .latents_center_direct(
    x$extra_data, prepared$x, x$sampling_weights[x$group_index])$direct
  steps <- 1e-5 * .multilpa_coordinate_scale(x, direct)[free]
  function(clustered, estimator, transform = NULL) {
    at <- function(point) {
      parameters <- .multilpa_decode(restore(point), x)
      posterior <- .multilpa_expectation(prepared$x, x$group_index, parameters,
                                         prepared$codes, unname(x$sampling_weights),
                                         x$extra_data)
      updated <- pieces
      updated$posteriors <- if (identical(level, "individuals")) {
        posterior$subject_posteriors[pieces$kept, , drop = FALSE]
      } else posterior$group_posteriors[pieces$kept, , drop = FALSE]
      estimator(updated)
    }
    size <- length(at(theta))
    derivative <- vapply(seq_along(theta), function(j) {
      up <- down <- theta
      up[j] <- up[j] + steps[j]
      down[j] <- down[j] - steps[j]
      (at(up) - at(down)) / (2 * steps[j])
    }, numeric(size))
    derivative <- matrix(derivative, nrow = size)
    if (!is.null(transform)) derivative <- transform %*% derivative
    # Appended after CR1 and thin-class masking: these rows are covariance
    # contributions, not additional observations or independent clusters.
    rbind(clustered, t(derivative %*% root))
  }
}

#' A square root of a step-one covariance, positive semidefinite allowed
#'
#' A sandwich or a constrained model (person centring makes the profile-
#' weighted means sum to zero within every unit) gives a covariance of
#' deficient rank, which a Cholesky factor refuses; its zero directions
#' carry no variance and contribute nothing.
#' @param covariance A symmetric covariance matrix.
#' @return `L` with `L %*% t(L)` equal to `covariance`, or
#'   `latents_singular_information` for a materially indefinite matrix.
#' @noRd
.multilpa_covariance_root <- function(covariance) {
  decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  values <- decomposition$values
  if (any(values < -sqrt(.Machine$double.eps) * max(abs(values), 1e-300))) {
    stop(errorCondition(paste(
      "The step-one covariance is not positive semidefinite, so its uncertainty",
      "cannot be added; `step_one = \"fixed\"` still works."),
      class = "latents_singular_information", call = NULL))
  }
  decomposition$vectors %*% diag(sqrt(pmax(values, 0)), length(values))
}

#' Sum influence contributions over the independent units, scaled by CR1
#'
#' Shared by every three-step method: the contributions are summed within
#' each group (`vcov_type = "cluster"`) or kept per observation, and scaled by
#' `sqrt(G / (G - 1))`, which scales every variance and covariance formed from
#' them, class means and differences alike.
#' @param influence Units-by-classes influence contributions.
#' @param pieces The list from `.multilpa_level_assignments()`.
#' @param vcov_type `"cluster"` or `"independent"`.
#' @return The scaled contributions, one row per independent unit.
#' @noRd
.multilpa_step_clusters <- function(influence, pieces, vcov_type) {
  units <- if (identical(vcov_type, "cluster")) {
    pieces$group_index
  } else seq_len(nrow(influence))
  clustered <- rowsum(influence, units, reorder = FALSE)
  .multilpa_require_clusters(
    nrow(clustered), pieces$n_classes,
    what = "A cluster-robust standard error for the class outcome means",
    unit = if (identical(vcov_type, "cluster")) "independent groups" else "observations",
    alternative = if (identical(vcov_type, "cluster")) {
      paste("`vcov_type = \"independent\"` treats every observation as its own",
            "unit instead, which ignores the nesting and must be reported as",
            "having done so.")
    })
  clustered * sqrt(.inference_cluster_factor(nrow(clustered)))
}

#' Withhold the standard error of a class carried by a single unit
#'
#' The sandwich of a class mean sums each independent unit's influence,
#' `sum_i w_ik (y_i - mu_k) / T_k` over the unit's observations. When every
#' observation with weight in class `k` lies in one unit, that unit's sum is
#' zero by the definition of `mu_k`, so the variance is zero by construction
#' rather than estimated, and a contrast with that class would borrow the
#' other class's variance alone. The units carrying a class are counted by
#' Kish's effective number, `(sum r)^2 / sum r^2` of the absolute weights
#' summed within each unit (absolute because BCH weights can be negative):
#' 1 for one unit, 2 for two equal ones. Below 1.5 the class's column is set
#' to `NA`, with a `latents_small_classes` warning, the rule the ML step-three
#' model applies to its responsibilities.
#'
#' @param clustered The scaled unit contributions from
#'   `.multilpa_step_clusters()`, one column per class.
#' @param carrying Units-by-classes weights saying how much of each class an
#'   observation carries: the method's weights, or the ML responsibilities.
#' @param pieces The list from `.multilpa_level_assignments()`.
#' @param vcov_type `"cluster"` or `"independent"`.
#' @return `clustered`, with the columns of thinly carried classes `NA`.
#' @noRd
.multilpa_step_thin_classes <- function(clustered, carrying, pieces, vcov_type) {
  clustered[, .multilpa_step_thin(carrying, pieces, vcov_type)] <- NA_real_
  clustered
}

#' Which classes are carried by a single independent unit
#'
#' The rule of `.multilpa_step_thin_classes()`, with its warning, for a
#' caller whose contributions have several columns per class.
#' @param carrying,pieces,vcov_type As in `.multilpa_step_thin_classes()`.
#' @return One logical per class, `TRUE` for a thinly carried class.
#' @noRd
.multilpa_step_thin <- function(carrying, pieces, vcov_type) {
  carried <- abs(carrying)
  if (identical(vcov_type, "cluster")) {
    carried <- rowsum(carried, pieces$group_index, reorder = FALSE)
  }
  carriers <- colSums(carried)^2 / colSums(carried^2)
  thin <- !(carriers >= 1.5)
  if (any(thin)) {
    several <- sum(thin) > 1L
    warning(warningCondition(sprintf(paste(
      "%s %s %s carried by a single %s, so %s standard error%s cannot be",
      "estimated and %s reported as NA."),
      if (several) "Classes" else "Class", paste(which(thin), collapse = ", "),
      if (several) "are each" else "is",
      if (identical(vcov_type, "cluster")) "group" else "observation",
      if (several) "their" else "its", if (several) "s" else "",
      if (several) "are" else "is"),
      class = "latents_small_classes", call = NULL))
  }
  thin
}

#' Pairwise differences between class outcome means
#'
#' @param estimates The per-class weighted means.
#' @param clustered Group-summed influence contributions, one column per class.
#' @param level,method The call's level and method, carried into the result.
#' @param quantile The normal quantile for the interval.
#' @param adjust A [stats::p.adjust()] method applied across the pairs.
#' @param vcov_type The variance requested, carried into the result's attribute.
#' @return One row per unordered pair of classes.
#' @noRd
.multilpa_step_pairs <- function(estimates, clustered, level, method, quantile,
                                 adjust, vcov_type) {
  n_classes <- length(estimates)
  if (n_classes < 2L) {
    stop(errorCondition(
      "A single class has no other class to be compared with.",
      class = "latents_inseparable_classes", call = NULL))
  }
  pairs <- utils::combn(n_classes, 2L)
  difference <- estimates[pairs[2L, ]] - estimates[pairs[1L, ]]
  errors <- sqrt(colSums((clustered[, pairs[2L, ], drop = FALSE] -
                            clustered[, pairs[1L, ], drop = FALSE])^2))
  statistic <- rep(NA_real_, length(errors))
  valid <- is.finite(errors) & errors > 0
  statistic[valid] <- difference[valid] / errors[valid]
  raw <- rep(NA_real_, length(errors))
  raw[valid] <- 2 * stats::pnorm(-abs(statistic[valid]))
  result <- data.frame(
    level = level, method = method, class = pairs[2L, ],
    reference_class = pairs[1L, ], estimate = difference,
    standard_error = errors, statistic = statistic, p_value = raw,
    p_adjusted = stats::p.adjust(raw, method = adjust),
    conf_low = difference - quantile * errors,
    conf_high = difference + quantile * errors,
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "adjust") <- adjust
  attr(result, "vcov_type") <- vcov_type
  column <- .latents_column
  percent <- format(100 * (2 * stats::pnorm(quantile) - 1))
  .latents_table(result, NULL,
                 sprintf("Distal outcome: class differences (%s, %s%% CI)", method, percent),
                 list(column("Class", "class", "integer"),
                      column("vs", "reference_class", "integer"),
                      column("Difference", "estimate"),
                      column(sprintf("%s%% CI", percent), c("conf_low", "conf_high"), "ci",
                             "estimate"),
                      column("p", "p_value", "p"),
                      column(sprintf("p (%s)", adjust), "p_adjusted", "p")))
}

#' The outcome, at the level being analysed
#' @return A numeric vector, one value per unit at that level.
#' @noRd
.multilpa_outcome_values <- function(object, data, outcome, level) {
  values <- data[[outcome]]
  .multilpa_check_row_count(object, data)
  if (identical(level, "individuals")) return(values)
  by_group <- split(values, object$group_index)
  constant <- vapply(by_group, function(v) length(unique(v)) == 1L, logical(1))
  if (!all(constant)) {
    stop(errorCondition(
      "For `level = \"groups\"` the outcome must be constant within each group.",
      class = "latents_bad_outcome", call = NULL))
  }
  vapply(by_group, function(v) v[[1L]], numeric(1))
}

#' Weights for the requested three-step method
#' @return A matrix with one row per unit and one column per class.
#' @noRd
.multilpa_step_weights <- function(pieces, method) {
  classes <- seq_len(pieces$n_classes)
  switch(method,
    bch = {
      inverse <- .multilpa_invert_errors(.multilpa_error_matrix(pieces))
      inverse[pieces$modal, , drop = FALSE]
    },
    proportional = pieces$posteriors,
    modal = vapply(classes, function(class) as.numeric(pieces$modal == class),
                   numeric(length(pieces$modal))))
}

#' Class-membership priors from a multinomial logit
#' @param gamma Coefficient matrix, predictors by non-reference classes.
#' @param design The model matrix.
#' @return A matrix of class probabilities, one row per unit.
#' @noRd
.multilpa_logit_prior <- function(gamma, design) {
  eta <- cbind(design %*% gamma, 0)
  scaled <- exp(eta - apply(eta, 1L, max))
  scaled / rowSums(scaled)
}

#' Relate latent class membership to external covariates
#'
#' @description
#' `r3step()` regresses membership in the classes of a fitted model on
#' covariates that were not used to define them, holding the measurement
#' model fixed. Treating the assigned classes as known outcomes distorts the
#' associations; the three-step approach instead treats each assignment as
#' an imperfect indicator of the true class, with the classification-error
#' probabilities of the fitted model built into a multinomial logit.
#' Coefficients are changes in the log odds of a class
#' against the reference class, the last one, holding the other covariates
#' constant; their exponentials are odds ratios. A transition fit can also
#' have its initial state and transitions regressed on covariates, a
#' separate specification. The correction assumes the classification-error
#' probabilities do not vary with the covariates within a class; the
#' standard errors treat them as fixed and ignore the choice of class
#' solution, so interval coverage can fall below the nominal level when the
#' classes are poorly separated; inspect the error matrix with
#' `get_results(x, "classification_errors")`. Associations are not causal
#' effects.
#'
#' @param x A covariate-free fit without weights or a noise component: from
#'   [multilpa()] (the profile family), from [lta()], or a trajectory model
#'   from [mixture_regression()] that classifies persons. A default `lta()`
#'   fit takes `covariates` or the Markov covariates; one with
#'   occasion-specific or second-order transitions, occasion-specific
#'   measurement, or ordinal or count indicators takes the Markov
#'   covariates only. Refused: fits with membership covariates, multilevel
#'   growth mixtures (`cluster`), the other `multilpa()` families,
#'   observation-level mixture regressions and stayer-class fits.
#' @param data The data frame carrying the covariates, in the fit's row
#'   order; columns shared with the fit are checked row by row.
#' @param covariates Names of numeric columns, complete and finite, that
#'   are not indicators of `x`. For `level = "groups"` each must be
#'   constant within a group; for a trajectory model, constant within a
#'   person and not part of the trajectory regression.
#' @param level `"individuals"` (default) predicts profile membership,
#'   `"groups"` group-class membership. A trajectory model predicts its
#'   person classes, reported as `level = "persons"` with outcomes
#'   `class_1`, `class_2` and so on.
#' @param ci_level Confidence level of the intervals (default 0.95).
#' @param vcov_type `"robust"` (default): the sandwich clustered on the
#'   independent units (the fit's groups; the persons of a trajectory model
#'   or a transition chain), which needs more units than the regression has
#'   coefficients, `"observed"` remaining available otherwise.
#'   `"observed"`: the observed information, treating every observation as
#'   independent in the pooled regression; with `by_group_class = TRUE`,
#'   the model-based variance of the two-level likelihood.
#' @param adjust Multiplicity correction across the covariate terms of every
#'   non-reference class (not the intercepts), a `stats::p.adjust()` method;
#'   `"BH"` by default.
#' @param by_group_class `FALSE` (default): one pooled multinomial logit with
#'   one intercept per profile. `TRUE`: the two-level form of the one-step
#'   model, with profile intercepts per group class, slopes shared across
#'   group classes, and the group-class proportions estimated jointly; the
#'   group classes then capture between-group differences in membership the
#'   covariates do not explain. Only for `level = "individuals"` on a
#'   `multilpa()` fit.
#' @param initial_covariates,transition_covariates For an [lta()] fit with
#'   one group class: covariates of the profile at the first occasion
#'   (constant within a person) and of every transition (read at the
#'   destination occasion, so they may vary over occasions). Either one fits
#'   the three-step latent Markov model (see *Transition fits*) instead of
#'   the occasion-by-occasion profile regression; `covariates`, `level` and
#'   `by_group_class` are refused then.
#'
#' @return A data frame with one row per non-reference class and term:
#'   `level`, `outcome`, `term`, `estimate`, `standard_error`, `statistic`,
#'   `p_value`, `p_adjusted` (`NA` on the intercept rows), `conf_low` and
#'   `conf_high`. Coefficients are log odds against the last class, named
#'   in the `reference_class` attribute; `vcov_type` and `adjust` are
#'   attributes too. With `by_group_class = TRUE` there is one intercept row
#'   per group class and non-reference profile, then the shared slopes and,
#'   at `level = "groups"`, the group-class logits, with the attributes
#'   `by_group_class` and `log_likelihood`.
#'
#'   With `initial_covariates` or `transition_covariates`, one row per
#'   coefficient of the step-three chain: `block` (`"initial"` or
#'   `"transition"`), `from` (the origin profile, `NA` for an initial
#'   logit), `to`, `term` and the estimate columns (`p_adjusted` on the
#'   covariate terms), with the attributes `log_likelihood`,
#'   `classification_errors` (the step-two error matrix, rows the true
#'   profile), `vcov_type` and `adjust`. Initial logits are
#'   against the last profile, transition logits against staying, as in
#'   `lta(initial_covariates =, transition_covariates =)`.
#'
#' @section Transition fits:
#' With `covariates` and `level = "individuals"`, an [lta()] fit is treated
#' occasion by occasion: the covariates predict the marginal profile at each
#' occasion, with the scores clustered by group under `"robust"`. That is
#' not a model of the initial state or the transitions.
#'
#' With `initial_covariates` or `transition_covariates`, each occasion is
#' classified from its own indicators by the transition fit's measurement
#' model (with the profile shares pooled over occasions) and assigned to its
#' modal profile; step three is a hidden Markov model for those assignments
#' with the classification errors fixed and the initial and transition
#' logits regressed on the covariates.
#' The transition fit's own posteriors are not used, because they draw on a
#' person's other occasions and their errors would not be independent across
#' occasions as step three assumes. The step-one fit must have one group
#' class and no covariates.
#' The step-three chain is homogeneous and first order, even when the
#' step-one transition structure varies over occasions or has higher order.
#' It uses one classification-error matrix pooled over occasions. Step-one
#' occasion-specific measurement is refused because its error probabilities
#' can change over occasions.
#'
#' @section Boundary estimates:
#' A step-three logit with a coefficient exceeding 20 in centred, scaled
#' design coordinates, or an average class probability at or below 1e-6
#' over the contributing design rows, is treated as separated or on a boundary. Wald
#' inference is refused for that regression, under both covariance types;
#' changing the covariance estimator does not make the boundary regular.
#'
#' @seealso [three_step()] for a distal outcome;
#'   `multilpa(profile_covariates = )` for the one-step alternative.
#' @examples
#' set.seed(21)
#' g <- rep(seq_len(50), each = 12)
#' x <- rnorm(600)
#' truth <- 1L + as.integer(runif(600) < plogis(-0.3 + 1.2 * x))
#' example_data <- data.frame(
#'   g = g, x = x,
#'   a = rnorm(600, ifelse(truth == 2L, 1.2, -1.2)),
#'   b = rnorm(600, ifelse(truth == 2L, 1.2, -1.2))
#' )
#' fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 4, seed = 1)
#' r3step(fit, example_data, "x")
#'
#' # Three-step regression of transition probabilities on covariates.
#' moves <- lta(course_engagement,
#'              vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'              id = "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
#' r3step(moves, course_engagement, transition_covariates = "sequence")
#' @export
r3step <- function(x, data, covariates = character(),
                   level = c("individuals", "groups"), ci_level = 0.95,
                   vcov_type = c("robust", "observed"),
                   adjust = c("BH", "holm", "hochberg", "hommel",
                                "bonferroni", "BY", "none"),
                   by_group_class = FALSE, initial_covariates = character(),
                   transition_covariates = character()) {
  level <- match.arg(level)
  vcov_type <- match.arg(vcov_type)
  adjust <- match.arg(adjust)
  if (length(initial_covariates) > 0L || length(transition_covariates) > 0L) {
    return(.r3step_transitions(x, data, covariates, level, ci_level, vcov_type,
                               adjust, by_group_class, initial_covariates,
                               transition_covariates))
  }
  trajectory <- .trajectory_is_fit(x)
  if (!trajectory) .multilpa_refuse_noise(x, "r3step()")
  if (.latents_is_weighted(x)) .latents_refuse_weights("`r3step()` yet")
  # Classed as three_step() classes the same mistakes (review 2026-10-27, INF-6).
  .latents_check_data_frame(data, substitute(data))
  .latents_require(
    "`covariates` must be a character vector naming at least one column of `data`" =
      is.character(covariates) && length(covariates) >= 1L && !anyNA(covariates))
  .r3step_check_columns(data, covariates, "covariates")
  .latents_require(
    "`covariates` must be numeric" =
      all(vapply(data[covariates], is.numeric, logical(1))),
    "`covariates` must not be missing" = !anyNA(data[covariates]),
    "`covariates` must be finite" =
      all(vapply(data[covariates], function(value) all(is.finite(value)), logical(1))),
    "`ci_level` must be a single number in (0, 1)" =
      is.numeric(ci_level) && length(ci_level) == 1L && is.finite(ci_level) &&
      ci_level > 0 && ci_level < 1,
    "`by_group_class` must be TRUE or FALSE" =
      is.logical(by_group_class) && length(by_group_class) == 1L &&
      !is.na(by_group_class)
  )
  if (trajectory) {
    # Trajectory models classify persons: one covariate value per person.
    if (identical(level, "groups") || by_group_class) {
      stop(errorCondition(paste(
        "A trajectory model classifies persons and has no group classes, so",
        "`level = \"groups\"` and `by_group_class = TRUE` do not apply."),
        class = "latents_bad_argument", call = NULL))
    }
    prepared <- .trajectory_r3step_inputs(x, data, covariates)
    pieces <- prepared$pieces
    design <- prepared$design
    level <- "persons"
  } else {
    .multilpa_check_three_step_fit(x)
    if (by_group_class && !identical(level, "individuals")) {
      stop(errorCondition(paste(
        "`by_group_class = TRUE` gives each group class its own profile",
        "intercepts, so it applies to `level = \"individuals\"` only."),
        class = "latents_bad_argument", call = NULL))
    }
    if (by_group_class && inherits(x, "multilpa_transitions")) {
      stop(errorCondition(paste(
        "`by_group_class = TRUE` is not available for a transition fit, whose",
        "group classes govern transitions rather than a profile mix."),
        class = "latents_bad_argument", call = NULL))
    }
    reused <- intersect(covariates, .multilpa_measured_items(x))
    if (length(reused) > 0L) {
      stop(errorCondition(sprintf(
        "Measurement indicator(s) %s already helped define the fitted classes and cannot be used as external predictors.",
        paste(sprintf("`%s`", reused), collapse = ", ")),
        class = c("multilpa_indicator_reused", "latents_bad_covariate"),
        call = NULL))
    }
    pieces <- .multilpa_level_assignments(x, level)
    .multilpa_check_row_count(x, data)
    .multilpa_check_alignment(x, data)
    .r3step_require_classes(pieces)
    design <- .multilpa_r3step_design(x, data, covariates, level)[pieces$kept, ,
                                                                   drop = FALSE]
  }
  .r3step_require_classes(pieces)
  errors <- .multilpa_error_matrix(pieces)
  if (by_group_class) {
    if (identical(vcov_type, "robust")) {
      .multilpa_require_clusters(
        length(unique(pieces$group_index)),
        (x$n_group_classes + length(covariates)) * (pieces$n_classes - 1L) +
          x$n_group_classes - 1L,
        what = "A cluster-robust covariance for the membership regression",
        alternative = paste("`vcov_type = \"observed\"` uses the model-based",
                            "information instead."))
    }
    return(.multilpa_r3step_by_group_class(
      x, design[, -1L, drop = FALSE], pieces, errors, ci_level, vcov_type,
      adjust))
  }
  ## The search runs on orthonormal coordinates of the design, which move with
  ## neither the units nor the origin of a covariate; a covariate held far from
  ## zero relative to its spread is otherwise nearly collinear with the
  ## intercept, and BFGS then stops far from the maximum or not at all.
  basis <- .r3step_design_basis(design) %||%
    diag(.multilpa_design_scale(design), ncol(design))
  terms <- colnames(design)
  design <- t(backsolve(basis, t(design), transpose = TRUE))
  ## The likelihood of the assigned class, given the covariates and an error
  ## matrix fixed at its step-one value.
  by_assigned <- t(errors)[pieces$modal, , drop = FALSE]
  n_free <- pieces$n_classes - 1L
  shape <- function(theta) matrix(theta, ncol(design), n_free)
  ## Refuse a sandwich the groups cannot support before paying for the fit.
  if (identical(vcov_type, "robust")) {
    .multilpa_require_clusters(
      length(unique(pieces$group_index)), ncol(design) * n_free,
      what = "A cluster-robust covariance for the membership regression",
      alternative = paste("`vcov_type = \"observed\"` uses the model-based",
                          "information instead, which does not allow for the",
                          "nesting."))
  }

  objective <- function(theta) {
    prior <- .multilpa_logit_prior(shape(theta), design)
    -sum(log(pmax(rowSums(prior * by_assigned), 1e-300)))
  }
  score_rows <- function(theta) {
    prior <- .multilpa_logit_prior(shape(theta), design)
    joint <- prior * by_assigned
    posterior <- joint / rowSums(joint)
    residual <- (posterior - prior)[, seq_len(n_free), drop = FALSE]
    do.call(cbind, lapply(seq_len(n_free), function(k) design * residual[, k]))
  }
  gradient <- function(theta) -colSums(score_rows(theta))

  ## The likelihood with the error matrix held fixed is not concave. When a
  ## covariate nearly separates the classes it can hold a ridge whose basin
  ## BFGS enters from zero and settles on, far from the finite maximum (Mplus
  ## UG ex7.3: -log L 81.7 on the ridge against 62.6 at the maximum). The
  ## proportional-assignment estimate, a multinomial logit weighted by the
  ## posteriors, is concave and finite and lies in the maximum's basin; the
  ## search runs from it and from zero and keeps the better converged end.
  proportional_objective <- function(theta) {
    prior <- .multilpa_logit_prior(shape(theta), design)
    -sum(pieces$posteriors * log(pmax(prior, 1e-300)))
  }
  proportional_gradient <- function(theta) {
    prior <- .multilpa_logit_prior(shape(theta), design)
    residual <- (pieces$posteriors - prior)[, seq_len(n_free), drop = FALSE]
    -colSums(do.call(cbind, lapply(seq_len(n_free), function(k) design * residual[, k])))
  }
  control <- list(maxit = 1000L, reltol = 1e-12)
  zero <- rep(0, ncol(design) * n_free)
  starts <- list(zero, stats::optim(zero, proportional_objective, proportional_gradient,
                                    method = "BFGS", control = control)$par)
  attempts <- lapply(starts, function(start)
    stats::optim(start, objective, gradient, method = "BFGS", control = control))
  converged <- vapply(attempts, function(attempt) identical(attempt$convergence, 0L),
                      logical(1))
  if (!any(converged)) {
    stop(errorCondition(
      sprintf("The membership regression did not converge (optim code %d).",
              attempts[[1L]]$convergence),
      class = "latents_no_converge", call = NULL))
  }
  attempts <- attempts[converged]
  fitted <- attempts[[which.min(vapply(attempts, `[[`, numeric(1), "value"))]]
  .r3step_require_interior(shape(fitted$par), design)
  information <- .multilpa_observed_hessian(
    function(step) objective(fitted$par + step),
    function(step) gradient(fitted$par + step),
    rep(1, length(fitted$par)), 1e-4)
  covariance <- information$inverse
  if (identical(vcov_type, "robust")) {
    scores <- rowsum(score_rows(fitted$par), pieces$group_index, reorder = FALSE)
    covariance <- .inference_sandwich(covariance, .multilpa_cross_product(scores),
                                      nrow(scores))
  }
  ## Back to the covariates' own units and origin, one logit column at a time.
  natural <- kronecker(diag(n_free), backsolve(basis, diag(ncol(design))))
  .multilpa_r3step_frame(as.vector(natural %*% fitted$par),
                         natural %*% covariance %*% t(natural), terms, n_free,
                         pieces$n_classes, level, ci_level, vcov_type, adjust)
}

#' Refuse a membership regression on a single class
#' @param pieces The level's assignments.
#' @return `NULL`, invisibly, or `latents_inseparable_classes`.
#' @noRd
.r3step_require_classes <- function(pieces) {
  if (pieces$n_classes < 2L) {
    stop(errorCondition(
      "A single class has no membership to predict.",
      class = "latents_inseparable_classes", call = NULL))
  }
  invisible(NULL)
}

#' Refuse covariate names that `data` does not have
#' @param data A data frame.
#' @param columns The names given.
#' @param argument The argument's name, for the message.
#' @return `NULL`, invisibly, or `latents_bad_data`.
#' @noRd
.r3step_check_columns <- function(data, columns, argument) {
  absent <- setdiff(columns, names(data))
  if (length(absent) == 0L) return(invisible(NULL))
  stop(errorCondition(sprintf(
    "`%s` names %s, which `data` does not have.%s", argument,
    .latents_name_list(absent), .latents_did_you_mean(absent, names(data))),
    class = "latents_bad_data", call = NULL))
}

#' r3step() with covariates of the initial profile and of the transitions
#' @return The table of `.lta_three_step()`.
#' @noRd
.r3step_transitions <- function(x, data, covariates, level, ci_level, vcov_type,
                                adjust, by_group_class, initial_covariates,
                                transition_covariates) {
  bad_argument <- function(message) {
    stop(errorCondition(message, class = "latents_bad_argument", call = NULL))
  }
  if (!inherits(x, c("multilpa_transitions", "multilpa_lta"))) {
    bad_argument(paste(
      "`initial_covariates` and `transition_covariates` apply to a transition",
      "model fitted by lta(); for profile membership use `covariates`."))
  }
  if (length(covariates) > 0L || !identical(level, "individuals") || by_group_class) {
    bad_argument(paste(
      "A three-step transition model takes `initial_covariates` and",
      "`transition_covariates` only: not `covariates`, `level` or",
      "`by_group_class`."))
  }
  if (.latents_is_weighted(x)) .latents_refuse_weights("`r3step()` yet")
  named <- function(value) is.character(value) && !anyNA(value) && !anyDuplicated(value)
  .latents_check_data_frame(data)
  .latents_require(
    "`initial_covariates` must be distinct column names of `data`" =
      named(initial_covariates),
    "`transition_covariates` must be distinct column names of `data`" =
      named(transition_covariates),
    "`ci_level` must be a single number in (0, 1)" =
      is.numeric(ci_level) && length(ci_level) == 1L && is.finite(ci_level) &&
      ci_level > 0 && ci_level < 1)
  .r3step_check_columns(data, initial_covariates, "initial_covariates")
  .r3step_check_columns(data, transition_covariates, "transition_covariates")
  .lta_three_step(x, data, initial_covariates, transition_covariates, ci_level,
                  vcov_type, adjust)
}

#' The covariate design at the level being analysed
#' @return A model matrix with an intercept.
#' @noRd
.multilpa_r3step_design <- function(object, data, covariates, level) {
  .multilpa_check_row_count(object, data)
  values <- data[covariates]
  if (identical(level, "groups")) {
    constant <- vapply(values, function(column) {
      all(vapply(split(column, object$group_index),
                 function(v) length(unique(v)) == 1L, logical(1)))
    }, logical(1))
    if (!all(constant)) {
      stop(errorCondition(
        "For `level = \"groups\"` every covariate must be constant within a group.",
        class = "latents_bad_covariate", call = NULL))
    }
    values <- values[!duplicated(object$group_index), , drop = FALSE]
  }
  .multilpa_r3step_check_rank(cbind(`(Intercept)` = 1, as.matrix(values)))
}

#' Refuse a rank-deficient R3STEP design
#' @param design A model matrix with an intercept.
#' @return `design`, or `latents_bad_covariate`.
#' @noRd
.multilpa_r3step_check_rank <- function(design) {
  if (is.null(.r3step_design_basis(design))) {
    stop(errorCondition(
      "The covariate design is rank deficient; remove constant or collinear predictors.",
      class = "latents_bad_covariate", call = NULL))
  }
  design
}

#' Orthonormal logit coordinates without a covariate-origin rank artefact
#'
#' Factor the centred, standardized design first and map its upper-triangular
#' factor back to the original coordinates. Direct QR of an intercept and a
#' predictor near 1e8 can lose rank even when their centred design is regular.
#' @param design A numeric design matrix with a leading column of ones.
#' @return An upper-triangular coordinate basis, or `NULL` for a deficient design.
#' @noRd
.r3step_design_basis <- function(design) {
  stopifnot(is.matrix(design), is.numeric(design), ncol(design) >= 1L,
            all(is.finite(design)), all(design[, 1L] == 1))
  if (ncol(design) == 1L) return(.multilpa_design_basis(design))
  centre <- colMeans(design[, -1L, drop = FALSE])
  centered <- sweep(design[, -1L, drop = FALSE], 2L, centre, "-")
  spread <- sqrt(colMeans(centered^2))
  if (any(!is.finite(spread)) || any(spread <= 0)) return(NULL)
  standardized <- cbind(1, sweep(centered, 2L, spread, "/"))
  basis <- .multilpa_design_basis(standardized)
  if (is.null(basis)) return(NULL)
  restore <- diag(c(1, spread), ncol(design))
  restore[1L, -1L] <- centre
  basis %*% restore
}

#' Refuse Wald inference for separated or empty step-three logit classes
#'
#' Coefficients are checked in orthonormal design coordinates, independently
#' of the units and origin of the supplied covariates. A class that has
#' negligible average probability over the design is a boundary even when its
#' intercept has not reached the coefficient limit.
#' @param coefficients A matrix with one row per design term and one column
#'   per non-reference logit.
#' @param design The contributing rows of that logit's design.
#' @param limit Largest absolute coefficient in orthonormal coordinates.
#' @param min_probability Largest probability treated as zero.
#' @return `NULL`, invisibly, or `latents_boundary_fit`.
#' @noRd
.r3step_require_interior <- function(coefficients, design, limit = 20,
                                      min_probability = 1e-6) {
  stopifnot(is.matrix(coefficients), is.matrix(design),
            nrow(coefficients) == ncol(design), ncol(coefficients) >= 1L,
            nrow(design) >= 1L)
  basis <- .r3step_design_basis(design)
  scaled <- if (is.null(basis)) {
    coefficients * .multilpa_design_scale(design)
  } else basis %*% coefficients
  probabilities <- .multilpa_logit_prior(coefficients, design)
  if (any(!is.finite(scaled)) ||
      any(abs(scaled) > limit) ||
      any(colMeans(probabilities) <= min_probability)) {
    stop(errorCondition(paste(
      "The step-three regression has separated or boundary logits, so finite",
      "Wald standard errors and tests are unavailable. Check sparse classes",
      "and transitions or reduce the membership covariates."),
      class = "latents_boundary_fit", call = NULL))
  }
  invisible(NULL)
}

#' Wald standard errors of a membership regression, refusing negative variances
#'
#' A negative diagonal says the information is not positive definite there:
#' that coefficient has no Wald standard error, which used to be reported as
#' 0. It is `NA`, with a `latents_singular_information` warning. Shared by the
#' pooled regression and the one with group-class intercepts.
#'
#' @param covariance The coefficients' covariance matrix.
#' @return One standard error per coefficient, `NA` where the variance is
#'   negative or not finite.
#' @noRd
.multilpa_r3step_standard_errors <- function(covariance) {
  variances <- diag(covariance)
  if (any(is.finite(variances) & variances < 0)) {
    warning(warningCondition(paste(
      "The covariance of the membership regression is not positive definite:",
      "some coefficients are not identified at this estimate, and their",
      "standard errors are NA."), class = "latents_singular_information",
      call = NULL))
  }
  valid <- is.finite(variances) & variances >= 0
  errors <- rep(NA_real_, length(variances))
  errors[valid] <- sqrt(variances[valid])
  errors
}

#' Assemble the R3STEP coefficient table
#' @return One row per non-reference class and term.
#' @noRd
.multilpa_r3step_frame <- function(estimates, covariance, terms, n_free,
                                   n_classes, level, ci_level, vcov_type,
                                   adjust) {
  errors <- .multilpa_r3step_standard_errors(covariance)
  statistic <- rep(NA_real_, length(errors))
  valid <- is.finite(errors) & errors > 0
  statistic[valid] <- estimates[valid] / errors[valid]
  quantile <- .inference_critical(ci_level, "one_minus")
  # Named as parameter_inference() names the same classes.
  prefix <- switch(level, groups = "group_class_", persons = "class_", "profile_")
  labels <- expand.grid(term = terms, outcome = paste0(prefix, seq_len(n_free)),
                        stringsAsFactors = FALSE)
  raw <- rep(NA_real_, length(errors))
  raw[valid] <- 2 * stats::pnorm(-abs(statistic[valid]))
  ## The intercepts are estimated, not tested: correcting across them would
  ## enlarge the family with hypotheses nobody asked about and make every
  ## covariate look less significant than the data say.
  tested <- labels$term != "(Intercept)"
  adjusted <- rep(NA_real_, length(raw))
  adjusted[tested] <- stats::p.adjust(raw[tested], method = adjust)
  result <- data.frame(
    level = level, outcome = labels$outcome, term = labels$term,
    estimate = estimates, standard_error = errors, statistic = statistic,
    p_value = raw, p_adjusted = adjusted,
    conf_low = estimates - quantile * errors,
    conf_high = estimates + quantile * errors,
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "reference_class") <- n_classes
  attr(result, "vcov_type") <- vcov_type
  attr(result, "adjust") <- adjust
  result
}

#' R3STEP with an intercept for every group class
#'
#' The membership regression of [r3step()] for profiles, with the two-level
#' structure the one-step model has: each group class has its own profile
#' intercepts, the covariate slopes are shared across group classes, and a
#' group's observations share one latent group class. With the step-one error
#' matrix `E` held fixed, a group's likelihood is
#' `sum_h eta_h prod_i sum_k pi_hk(z_i) E[k, a_i]`, where `a_i` is the assigned
#' profile. Its scores are posterior-weighted residuals, summed within group,
#' so the observed information and the group-clustered sandwich both come from
#' the same per-group rows.
#'
#' @param x The step-one fit.
#' @param z The covariate matrix, one row per observation, no intercept.
#' @param pieces The result of `.multilpa_level_assignments()` for individuals.
#' @param errors The step-one classification-error matrix.
#' @param ci_level,vcov_type,adjust As for [r3step()].
#' @return The tidy table [r3step()] returns, with group-class intercept rows
#'   and the group-class logits.
#' @noRd
.multilpa_r3step_by_group_class <- function(x, z, pieces, errors, ci_level,
                                            vcov_type, adjust) {
  model <- .multilpa_r3step_group_model(x, z, pieces, errors)
  n_classes <- pieces$n_classes
  n_free <- n_classes - 1L
  n_group_classes <- x$n_group_classes
  fitted <- stats::optim(model$start, model$objective, model$gradient,
                         method = "BFGS",
                         control = list(maxit = 1000L, reltol = 1e-12))
  if (!identical(fitted$convergence, 0L)) {
    stop(errorCondition(
      sprintf("The membership regression did not converge (optim code %d).",
              fitted$convergence),
      class = "latents_no_converge", call = NULL))
  }
  model$check_interior(fitted$par)
  information <- .multilpa_observed_hessian(
    function(step) model$objective(fitted$par + step),
    function(step) model$gradient(fitted$par + step),
    rep(1, length(fitted$par)), 1e-4)
  covariance <- information$inverse
  if (identical(vcov_type, "robust")) {
    scores <- model$group_scores(fitted$par)
    covariance <- .inference_sandwich(covariance, .multilpa_cross_product(scores),
                                      nrow(scores))
  }
  ## Back to the units and origin the covariates arrived in.
  estimates <- as.vector(model$transform %*% fitted$par)
  covariance <- model$transform %*% covariance %*% t(model$transform)
  .multilpa_r3step_group_frame(estimates, covariance, colnames(z), n_free,
                               n_group_classes, n_classes, ci_level, vcov_type,
                               adjust, -fitted$value)
}

#' The likelihood of R3STEP with group-class intercepts
#'
#' @param x The step-one fit.
#' @param z The covariate matrix, one row per observation, no intercept.
#' @param pieces The result of `.multilpa_level_assignments()` for individuals.
#' @param errors The step-one classification-error matrix.
#' @return A list with `objective` (negative log likelihood), `gradient`,
#'   `group_scores` (one row per group), `start`, and `transform`, the matrix
#'   that returns the search coordinates to the covariates' own units and
#'   origin.
#' @noRd
.multilpa_r3step_group_model <- function(x, z, pieces, errors) {
  n_classes <- pieces$n_classes
  n_free <- n_classes - 1L
  n_group_classes <- x$n_group_classes
  n_covariates <- ncol(z)
  group_index <- pieces$group_index
  n_groups <- max(group_index)
  by_assigned <- t(errors)[pieces$modal, , drop = FALSE]
  ## The search runs on the orthonormal coordinates of `cbind(1, z)` (see
  ## `.multilpa_design_basis()`), which move with neither the units nor the
  ## origin of a covariate. With `cbind(1, z) = Q R`, `z = r' + Q_z R_z`, so
  ## `alpha_h + z beta = (alpha_h + r' beta) + Q_z (R_z beta)`: the search
  ## estimates the bracketed intercepts and `R_z beta`, and `transform` maps
  ## them back to the covariates' own units and origin.
  basis <- .r3step_design_basis(cbind(1, z))
  shift <- if (is.null(basis)) rep(0, n_covariates) else basis[1L, -1L]
  slope_basis <- if (is.null(basis)) {
    diag(.multilpa_design_scale(z), n_covariates)
  } else basis[-1L, -1L, drop = FALSE]
  scaled_z <- t(backsolve(slope_basis, t(sweep(z, 2L, shift, "-")),
                          transpose = TRUE))
  n_alpha <- n_group_classes * n_free
  n_beta <- n_covariates * n_free
  slope_inverse <- backsolve(slope_basis, diag(n_covariates))
  ## beta = R_z^-1 beta'; alpha_h = alpha'_h - r' R_z^-1 beta'.
  transform <- diag(n_alpha + n_beta + n_group_classes - 1L)
  transform[n_alpha + seq_len(n_beta), n_alpha + seq_len(n_beta)] <-
    kronecker(diag(n_free), slope_inverse)
  transform[seq_len(n_alpha), n_alpha + seq_len(n_beta)] <-
    kronecker(diag(n_free), -matrix(1, n_group_classes, 1L) %*%
                (shift %*% slope_inverse))
  unpack <- function(theta) {
    list(alpha = matrix(theta[seq_len(n_alpha)], n_group_classes, n_free),
         beta = matrix(theta[n_alpha + seq_len(n_beta)], n_covariates, n_free),
         gamma = theta[n_alpha + n_beta + seq_len(n_group_classes - 1L)])
  }
  pieces_at <- function(theta) {
    parameters <- unpack(theta)
    slope_part <- scaled_z %*% parameters$beta
    per_class <- lapply(seq_len(n_group_classes), function(h) {
      linear <- cbind(sweep(slope_part, 2L, parameters$alpha[h, ], "+"), 0)
      prior <- exp(linear - apply(linear, 1L, max))
      prior <- prior / rowSums(prior)
      joint <- prior * by_assigned
      marginal <- pmax(rowSums(joint), 1e-300)
      list(prior = prior, posterior = joint / marginal,
           evidence = as.vector(rowsum(log(marginal), group_index, reorder = FALSE)))
    })
    log_eta <- c(parameters$gamma, 0)
    log_eta <- log_eta - .multilpa_log_sum_exp(matrix(log_eta, 1L))
    scores <- sweep(vapply(per_class, `[[`, numeric(n_groups), "evidence"),
                    2L, log_eta, "+")
    log_group <- .multilpa_log_sum_exp(scores)
    list(per_class = per_class, eta = exp(log_eta),
         group_posteriors = exp(sweep(scores, 1L, log_group, "-")),
         log_likelihood = sum(log_group))
  }
  group_scores <- function(theta) {
    state <- pieces_at(theta)
    residuals <- lapply(state$per_class, function(piece) {
      (piece$posterior - piece$prior)[, seq_len(n_free), drop = FALSE]
    })
    alpha_block <- do.call(cbind, lapply(seq_len(n_free), function(k) {
      vapply(seq_len(n_group_classes), function(h) {
        state$group_posteriors[, h] *
          as.vector(rowsum(residuals[[h]][, k], group_index, reorder = FALSE))
      }, numeric(n_groups))
    }))
    beta_block <- do.call(cbind, lapply(seq_len(n_free), function(k) {
      Reduce(`+`, lapply(seq_len(n_group_classes), function(h) {
        state$group_posteriors[, h] *
          rowsum(scaled_z * residuals[[h]][, k], group_index, reorder = FALSE)
      }))
    }))
    gamma_block <- sweep(state$group_posteriors[, seq_len(n_group_classes - 1L),
                                                drop = FALSE],
                         2L, state$eta[seq_len(n_group_classes - 1L)], "-")
    cbind(alpha_block, beta_block, gamma_block)
  }
  objective <- function(theta) -pieces_at(theta)$log_likelihood
  gradient <- function(theta) -colSums(group_scores(theta))

  ## Started from the step-one solution: the profile mix within each group
  ## class, no slopes, and the group-class proportions.
  mix <- pmax(x$profile_probabilities, 1e-8)
  start <- c(as.vector(log(mix[, seq_len(n_free), drop = FALSE] / mix[, n_classes])),
             rep(0, n_beta),
             log(pmax(x$group_probabilities[seq_len(n_group_classes - 1L)], 1e-8) /
                   max(x$group_probabilities[n_group_classes], 1e-8)))
  check_interior <- function(theta) {
    parameters <- unpack(theta)
    invisible(lapply(seq_len(n_group_classes), function(h) {
      .r3step_require_interior(rbind(parameters$alpha[h, ], parameters$beta),
                               cbind(1, scaled_z))
    }))
    if (length(parameters$gamma) > 0L) {
      .r3step_require_interior(matrix(parameters$gamma, 1L),
                               matrix(1, n_groups, 1L))
    }
    invisible(NULL)
  }
  list(objective = objective, gradient = gradient, group_scores = group_scores,
       start = start, transform = transform, check_interior = check_interior)
}

#' Assemble the table of R3STEP with group-class intercepts
#' @return One row per coefficient: profile intercepts for every group class,
#'   the shared slopes, then the group-class logits.
#' @noRd
.multilpa_r3step_group_frame <- function(estimates, covariance, terms, n_free,
                                         n_group_classes, n_classes, ci_level,
                                         vcov_type, adjust, log_likelihood) {
  errors_se <- .multilpa_r3step_standard_errors(covariance)
  ## Alpha is group classes by free classes and beta covariates by free
  ## classes, both stored column-major, so `expand.grid()` with the row index
  ## varying fastest labels them in the order they are packed.
  alpha_labels <- expand.grid(term = paste0("group_class_", seq_len(n_group_classes)),
                              outcome = paste0("profile_", seq_len(n_free)),
                              stringsAsFactors = FALSE)
  beta_labels <- expand.grid(term = terms,
                             outcome = paste0("profile_", seq_len(n_free)),
                             stringsAsFactors = FALSE)
  labels <- rbind(
    data.frame(level = "individuals", outcome = alpha_labels$outcome,
               term = alpha_labels$term, stringsAsFactors = FALSE),
    data.frame(level = "individuals", outcome = beta_labels$outcome,
               term = beta_labels$term, stringsAsFactors = FALSE),
    data.frame(level = rep("groups", n_group_classes - 1L),
               outcome = sprintf("group_class_%d", seq_len(n_group_classes - 1L)),
               term = rep("(Intercept)", n_group_classes - 1L),
               stringsAsFactors = FALSE))
  statistic <- rep(NA_real_, length(estimates))
  valid <- is.finite(errors_se) & errors_se > 0
  statistic[valid] <- estimates[valid] / errors_se[valid]
  raw <- rep(NA_real_, length(estimates))
  raw[valid] <- 2 * stats::pnorm(-abs(statistic[valid]))
  ## Only the covariate slopes are tested; intercepts and the group-class
  ## logits are estimated, not hypotheses, as in the pooled regression.
  tested <- labels$level == "individuals" & labels$term %in% terms
  adjusted <- rep(NA_real_, length(raw))
  adjusted[tested] <- stats::p.adjust(raw[tested], method = adjust)
  quantile <- .inference_critical(ci_level, "one_minus")
  result <- data.frame(
    labels, estimate = estimates, standard_error = errors_se,
    statistic = statistic, p_value = raw, p_adjusted = adjusted,
    conf_low = estimates - quantile * errors_se,
    conf_high = estimates + quantile * errors_se,
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "reference_class") <- n_classes
  attr(result, "vcov_type") <- vcov_type
  attr(result, "adjust") <- adjust
  attr(result, "by_group_class") <- TRUE
  attr(result, "log_likelihood") <- log_likelihood
  result
}

#' Every item that entered a fit's measurement model
#'
#' The fit's indicators and, for a fit with `local_dependence`, the items each
#' joint indicator was made from: `"a:b"` stands in `vars`, and `"a"` and
#' `"b"` helped define the classes just as much.
#' @param x A fitted model.
#' @return A character vector of column names.
#' @noRd
.multilpa_measured_items <- function(x) {
  unique(c(x$vars, x$original_vars))
}
