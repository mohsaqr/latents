#' Fit a finite mixture of regressions
#'
#' Fits a regression whose coefficients differ across latent classes: the
#' data are assumed to come from `n_classes` subpopulations, each with its own
#' regression of the outcome on the predictors, and the model estimates the
#' regressions, the class sizes and every row's posterior class membership at
#' once. This is clusterwise, latent class or mixture regression (DeSarbo and
#' Cron 1988; Wedel and DeSarbo 1995).
#'
#' @section Nesting:
#' `class_level` and `n_group_classes` choose one of three likelihoods.
#' \describe{
#'   \item{Single level (`id = NULL`)}{Every row is its own unit and has its own
#'     class.}
#'   \item{`class_level = "group"`}{Every row of a group (`id`) belongs to the
#'     same class, so a class is a kind of *group*: the regression mixture for
#'     repeated measures, `flexmix(y ~ x | id)`. Each group's evidence is the
#'     product of its rows' densities.}
#'   \item{`class_level = "observation"` with `id` and `n_group_classes >= 2`}{
#'     Rows keep their own class, and a second-level group class shifts how
#'     probable each regression class is within its groups (Vermunt 2003).
#'     With `n_group_classes = 1` the groups only enter the `"robust"` standard
#'     errors, which are then clustered on `id`.}
#' }
#'
#' @section Families:
#' `"gaussian"` (identity link, a residual standard deviation per class, or
#' one shared under `variance = "equal"`), `"binomial"` (logit link; the
#' outcome is 0/1, logical, a two-level factor whose second level is the
#' success, or `cbind(successes, failures)`) and `"poisson"` (log link;
#' `offset()` terms in the formula are honoured).
#'
#' A mixture of Bernoulli regressions with one trial per unit is not
#' identified -- any mixture of Bernoullis is again a Bernoulli -- so a binary
#' outcome requires `class_level = "group"` (several rows per class
#' assignment) or binomial counts with more than one trial per row. The
#' request is refused otherwise.
#'
#' @param formula A model formula `outcome ~ predictors`. Factors, interactions
#'   and `offset()` are supported.
#' @param data A data frame.
#' @param n_classes Number of regression classes.
#' @param family `"gaussian"`, `"binomial"` or `"poisson"`.
#' @param id `NULL`, or the name of the column identifying groups (persons,
#'   schools, ...).
#' @param class_level `"observation"` (each row has its own class) or
#'   `"group"` (all rows of a group share a class; needs `id`).
#' @param n_group_classes Number of second-level group classes, for
#'   `class_level = "observation"` with `id`. `1` fits no second level.
#' @param common `NULL`, or a one-sided formula naming predictor terms whose
#'   coefficients are shared by every class, such as `~ age`. Each term must
#'   appear in `formula`.
#' @param membership One-sided formula of covariates predicting class
#'   membership (a multinomial logit, the first class as the reference).
#'   Row-level for `"observation"`; constant within groups for `"group"`.
#' @param group_membership One-sided formula of group-level covariates
#'   predicting the group class, for the two-level model. Must be constant
#'   within groups.
#' @param variance `"varying"` (a residual standard deviation per class) or
#'   `"equal"`. Gaussian family only.
#' @param n_starts Number of random starts, besides one start built from the
#'   residuals of a pooled regression. Every start runs 50 EM iterations; the
#'   better half (at least two) continue to convergence.
#' @param max_iter Maximum EM iterations per start.
#' @param tol Relative convergence tolerance on the log likelihood.
#' @param min_variance Floor for a Gaussian residual variance, relative to the
#'   outcome's variance. A class that reaches it is fitting too few rows
#'   exactly; such starts are set aside.
#' @param seed `NULL` or an integer seed. The caller's random state is
#'   restored afterwards.
#' @param missing `"error"` refuses rows with missing values in any variable
#'   the model uses; `"omit"` drops them with a `latents_rows_dropped`
#'   warning stating how many.
#' @param select_start `"likelihood"` keeps the start with the highest
#'   likelihood; `"converged"` prefers the best start that converged.
#' @param vcov_type Covariance of the estimates stored with the fit:
#'   `"observed"` (inverse observed information), `"robust"` (sandwich,
#'   clustered on the top-level unit, or on `id` for a single-level fit given
#'   `id`) or `"opg"` (outer product of the scores). `"none"` skips inference.
#' @param weights `NULL`, or the name of a numeric column of `data` with a
#'   sampling weight per independent unit: per row without `id`, per `id`
#'   group otherwise (constant within it). Pseudo maximum likelihood with the
#'   weights scaled to sum to the number of units; `vcov_type` defaults to
#'   `"robust"` and refuses `"observed"` and `"opg"`.
#'
#' @return An object of class `latents_mixture_regression`. Read it with
#'   [as.data.frame()] (the coefficient table) or [get_results()] (every
#'   table, by name), and use [predict()], [plot()], [summary()], [coef()],
#'   [vcov()], [confint()], [logLik()] and [nobs()] on it. The tables are
#'   described on [get_results.latents_mixture_regression()].
#'
#' @section Conditions:
#'   `latents_bad_argument` for an inconsistent specification;
#'   `latents_bad_data` for an outcome outside the family's support or a
#'   covariate that varies within a group where it may not;
#'   `latents_not_identified` for a binary outcome with one trial per class
#'   assignment; `latents_missing_data` for missing values under
#'   `missing = "error"`; `latents_rows_dropped` (warning) under
#'   `missing = "omit"`; `latents_no_valid_start` when every start
#'   degenerated; `latents_unconverged` (warning) when the selected start did
#'   not converge; `latents_degenerate_start` (warning) when some starts
#'   degenerated; `latents_separation` (warning) when a coefficient diverged.
#'
#' @references
#' DeSarbo, W. S., & Cron, W. L. (1988). A maximum likelihood methodology for
#' clusterwise linear regression. *Journal of Classification*, 5, 249--282.
#'
#' Wedel, M., & DeSarbo, W. S. (1995). A mixture likelihood approach for
#' generalized linear models. *Journal of Classification*, 12, 21--55.
#'
#' Follmann, D. A., & Lambert, D. (1991). Identifiability of finite mixtures
#' of logistic regression models. *Journal of Statistical Planning and
#' Inference*, 27, 375--381.
#'
#' Vermunt, J. K. (2003). Multilevel latent class models. *Sociological
#' Methodology*, 33, 213--239.
#'
#' Leisch, F. (2004). FlexMix: A general framework for finite mixture models
#' and latent class regression in R. *Journal of Statistical Software*, 11(8).
#'
#' @seealso [enumerate_regressions()] to compare numbers of classes.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' fit
#' as.data.frame(fit)
#' get_results(fit, "classes")
#'
#' # One class per student, several rows per student:
#' by_student <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                                  id = "student", class_level = "group",
#'                                  n_starts = 3, seed = 1)
#' get_results(by_student, "fit")
#' @export
mixture_regression <- function(formula, data, n_classes,
                   family = c("gaussian", "binomial", "poisson"),
                   id = NULL, class_level = c("observation", "group"),
                   n_group_classes = 1L, common = NULL, membership = ~1,
                   group_membership = ~1, variance = c("varying", "equal"),
                   n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                   min_variance = 1e-6, seed = NULL,
                   missing = c("error", "omit"),
                   select_start = c("likelihood", "converged"),
                   vcov_type = c("observed", "robust", "opg", "none"),
                   weights = NULL) {
  vcov_defaulted <- missing(vcov_type)
  family <- match.arg(family)
  class_level <- match.arg(class_level)
  variance <- match.arg(variance)
  missing <- match.arg(missing)
  select_start <- match.arg(select_start)
  vcov_type <- match.arg(vcov_type)
  stopifnot(
    "`formula` must be a two-sided formula" =
      inherits(formula, "formula") && length(formula) == 3L,
    "`data` must be a data frame" = is.data.frame(data),
    "`n_classes` must be a single positive integer" =
      .mixture_is_count(n_classes),
    "`n_group_classes` must be a single positive integer" =
      .mixture_is_count(n_group_classes),
    "`id` must be NULL or a single column name" =
      is.null(id) || (is.character(id) && length(id) == 1L && !is.na(id)),
    "`common` must be NULL or a one-sided formula" =
      is.null(common) || (inherits(common, "formula") && length(common) == 2L),
    "`membership` must be a one-sided formula" =
      inherits(membership, "formula") && length(membership) == 2L,
    "`group_membership` must be a one-sided formula" =
      inherits(group_membership, "formula") && length(group_membership) == 2L,
    "`n_starts` must be a single non-negative integer" =
      length(n_starts) == 1L && is.finite(n_starts) && n_starts >= 0,
    "`max_iter` must be a single non-negative integer" =
      length(max_iter) == 1L && is.finite(max_iter) && max_iter >= 0,
    "`tol` must be a single positive number" =
      is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0,
    "`min_variance` must be a single positive number" =
      is.numeric(min_variance) && length(min_variance) == 1L &&
      is.finite(min_variance) && min_variance > 0,
    "`seed` must be NULL or a single number" =
      is.null(seed) || (is.numeric(seed) && length(seed) == 1L)
  )
  n_classes <- as.integer(n_classes)
  n_group_classes <- as.integer(n_group_classes)
  spec <- .mixture_spec(formula, data, n_classes, family, id, class_level,
                       n_group_classes, common, membership, group_membership,
                       variance, min_variance, missing)
  if (!is.null(weights)) {
    spec <- .mixture_weight_spec(spec, data, weights)
    if (!identical(vcov_type, "none")) {
      vcov_type <- .latents_weighted_vcov(TRUE, vcov_type, vcov_defaulted)
    }
  }
  results <- .mixture_with_seed(seed, .mixture_run_starts(
    spec, as.integer(n_starts), as.integer(max_iter), tol))
  fit <- .mixture_select(spec, results, select_start)
  fit$call <- match.call()
  fit$settings <- list(n_starts = as.integer(n_starts),
                       max_iter = as.integer(max_iter), tol = tol,
                       min_variance = min_variance, seed = seed,
                       select_start = select_start)
  fit$weights <- weights
  fit$sampling_weights <- spec$sampling_weights
  if (!identical(vcov_type, "none")) {
    fit$inference <- .mixture_inference(fit, vcov_type)
  }
  fit
}

#' Attach sampling weights to a mixture-regression specification
#'
#' The independent units are the rows, or the `id` groups when `id` is given
#' (they cluster the rows even when classes sit on the rows). Weights are read
#' on the rows the complete-case pass kept.
#'
#' @param spec The specification.
#' @param data The data frame passed to the fit.
#' @param weights Name of the weight column.
#' @return `spec` with `sampling_weights` (per unit), `row_weights` and
#'   `weights`.
#' @noRd
.mixture_weight_spec <- function(spec, data, weights) {
  if (is.character(weights) && length(weights) == 1L && weights %in% names(data) &&
      length(spec$kept_rows) < nrow(data) && anyNA(data[[weights]])) {
    .latents_refuse_weights("missing weights; drop those rows first")
  }
  frame <- data[spec$kept_rows, , drop = FALSE]
  unit_index <- if (is.null(spec$group_index)) seq_len(spec$n) else spec$group_index
  n_units <- if (is.null(spec$group_index)) spec$n else spec$n_groups
  unit_weights <- .latents_sampling_weights(frame, weights, unit_index, n_units)
  spec$sampling_weights <- unit_weights
  spec$row_weights <- unit_weights[unit_index]
  spec$weights <- weights
  spec
}

#' Is a value a single positive whole number?
#' @noRd
.mixture_is_count <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 &&
    abs(x - round(x)) < sqrt(.Machine$double.eps)
}

#' Evaluate code under a seed, restoring the caller's random state
#' @param seed `NULL` or a number.
#' @param code Expression, evaluated lazily.
#' @return The value of `code`.
#' @noRd
.mixture_with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())  # nolint: object_name_linter. R's name for the RNG state.
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE, after = FALSE)
  set.seed(seed)
  code
}

#' Build the model specification from a formula and data
#'
#' Resolves the outcome, the class-specific and shared designs, the offset,
#' the membership designs and the group index, after one complete-case pass
#' over every variable the model uses.
#'
#' @return A list consumed by the engine; see the fields assigned below.
#' @noRd
.mixture_spec <- function(formula, data, n_classes, family, id, class_level,
                         n_group_classes, common, membership,
                         group_membership, variance, min_variance, missing) {
  bad_argument <- function(message) {
    stop(errorCondition(message, class = "latents_bad_argument", call = NULL))
  }
  nesting <- if (identical(class_level, "group")) "group" else if (
    !is.null(id) && n_group_classes > 1L) "two-level" else "observation"
  if (identical(class_level, "group") && is.null(id)) {
    bad_argument(paste("`class_level = \"group\"` assigns every row of a group",
                       "to one class, so it needs `id` naming the groups."))
  }
  if (n_group_classes > 1L && is.null(id)) {
    bad_argument(sprintf(paste(
      "`n_group_classes = %d` needs groups to classify; pass `id`."),
      n_group_classes))
  }
  if (n_group_classes > 1L && identical(class_level, "group")) {
    bad_argument(paste(
      "Group classes shift the class shares of rows within a group; with",
      "`class_level = \"group\"` a group already has a single class. Use",
      "`class_level = \"observation\"` for the two-level model."))
  }
  if (!identical(variance, "varying") && !identical(family, "gaussian")) {
    bad_argument("`variance` applies to the gaussian family only.")
  }
  if (!is.null(id) && !id %in% names(data)) {
    bad_argument(sprintf("`id` names column `%s`, which `data` does not have.",
                         id))
  }
  if (!identical(nesting, "two-level") &&
      length(all.vars(group_membership)) > 0L) {
    bad_argument(paste("`group_membership` applies to the two-level model",
                       "only (`id` with `n_group_classes >= 2`)."))
  }
  if (!isTRUE(attr(stats::terms(membership), "intercept") == 1L)) {
    bad_argument("`membership` must keep its intercept.")
  }

  used <- unique(c(all.vars(formula), all.vars(common %||% ~1),
                   all.vars(membership), all.vars(group_membership), id))
  absent <- setdiff(used, names(data))
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf("`data` has no column%s %s.",
                                if (length(absent) > 1L) "s" else "",
                                paste(sprintf("`%s`", absent), collapse = ", ")),
                        class = "latents_bad_data", call = NULL))
  }
  complete <- stats::complete.cases(data[used])
  if (!all(complete)) {
    if (identical(missing, "error")) {
      stop(errorCondition(sprintf(paste(
        "%d of %d rows have a missing value in a variable the model uses.",
        "Describe the missingness first; pass `missing = \"omit\"` to drop",
        "those rows."), sum(!complete), length(complete)),
        class = "latents_missing_data", call = NULL))
    }
    warning(warningCondition(sprintf(paste(
      "Dropped %d of %d rows with a missing value in a model variable",
      "(complete-case analysis)."), sum(!complete), length(complete)),
      class = "latents_rows_dropped", call = NULL))
  }
  kept_rows <- which(complete)
  data <- data[kept_rows, , drop = FALSE]
  n <- nrow(data)

  frame <- stats::model.frame(formula, data = data, na.action = stats::na.fail,
                              drop.unused.levels = TRUE)
  model_terms <- stats::terms(frame)
  full_design <- stats::model.matrix(model_terms, frame)
  offset <- stats::model.offset(frame) %||% rep(0, n)
  response <- .mixture_response(stats::model.response(frame), family)

  common_names <- character()
  if (!is.null(common)) {
    common_design <- stats::model.matrix(
      stats::update(common, ~ . + 0), data = data)
    common_names <- colnames(common_design)
    unknown <- setdiff(common_names, colnames(full_design))
    if (length(unknown) > 0L) {
      bad_argument(sprintf(paste(
        "`common` names %s, which `formula` does not produce; every shared",
        "term must also be in `formula`."),
        paste(sprintf("`%s`", unknown), collapse = ", ")))
    }
    if ("(Intercept)" %in% common_names) {
      bad_argument("The intercept cannot be shared through `common`.")
    }
  }
  varying_names <- setdiff(colnames(full_design), common_names)
  if (length(varying_names) == 0L) {
    bad_argument("At least one term must vary across classes.")
  }
  x <- full_design[, varying_names, drop = FALSE]
  z <- full_design[, common_names, drop = FALSE]

  group_values <- if (!is.null(id)) data[[id]] else NULL
  group_levels <- if (!is.null(id)) sort(unique(group_values)) else NULL
  group_index <- if (!is.null(id)) match(group_values, group_levels) else NULL
  n_groups <- length(group_levels)

  membership_design <- .mixture_membership_design(membership, data)
  group_membership_design <- .mixture_membership_design(group_membership, data)
  row_membership <- membership_design$matrix
  w <- row_membership
  v <- matrix(1, max(n_groups, 1L), 1L, dimnames = list(NULL, "(Intercept)"))
  if (identical(nesting, "group")) {
    w <- .mixture_group_design(row_membership, group_index, "membership")
  }
  if (identical(nesting, "two-level")) {
    w <- row_membership[, colnames(row_membership) != "(Intercept)",
                        drop = FALSE]
    v <- .mixture_group_design(group_membership_design$matrix,
                              group_index, "group_membership")
  }

  if (identical(family, "binomial") && all(response$trials == 1) &&
      !identical(nesting, "group")) {
    stop(errorCondition(paste(
      "A mixture of logistic regressions with one binary outcome per class",
      "assignment is not identified: a mixture of Bernoulli distributions is",
      "itself a Bernoulli distribution (Follmann and Lambert 1991). Give each",
      "class assignment several outcomes with `id` and",
      "`class_level = \"group\"`, or supply binomial counts as",
      "`cbind(successes, failures)`."),
      class = "latents_not_identified", call = NULL))
  }
  units <- switch(nesting, observation = n, n_groups)
  if (units <= n_classes) {
    stop(errorCondition(sprintf(paste(
      "%d %s cannot support %d classes."), units,
      if (identical(nesting, "observation")) "rows" else "groups", n_classes),
      class = "latents_bad_data", call = NULL))
  }
  if (identical(nesting, "two-level") && n_groups <= n_group_classes) {
    stop(errorCondition(sprintf("%d groups cannot support %d group classes.",
                                n_groups, n_group_classes),
                        class = "latents_bad_data", call = NULL))
  }

  y_variance <- if (identical(family, "gaussian")) stats::var(response$y) else 1
  list(
    family = family, nesting = nesting, variance = variance,
    n = n, n_classes = n_classes,
    n_group_classes = if (identical(nesting, "two-level")) n_group_classes else 1L,
    y = response$y, trials = response$trials,
    log_normalizer = switch(family,
      gaussian = rep(0, n),
      binomial = lchoose(response$trials, response$y),
      poisson = -lgamma(response$y + 1)),
    x = x, z = z, offset = as.numeric(offset),
    w = w, v = v,
    intercept_only = ncol(w) == 1L && identical(colnames(w), "(Intercept)"),
    group_intercept_only = ncol(v) == 1L &&
      identical(colnames(v), "(Intercept)"),
    id = id, group_index = group_index, group_levels = group_levels,
    n_groups = n_groups,
    min_variance = min_variance * max(y_variance, .Machine$double.eps),
    kept_rows = kept_rows,
    terms = model_terms, xlevels = stats::.getXlevels(model_terms, frame),
    contrasts = attr(full_design, "contrasts"), model_data = data[used],
    formula = formula, common = common, membership = membership,
    membership_design = membership_design,
    group_membership_design = group_membership_design,
    group_membership = group_membership, response_name =
      deparse1(formula[[2L]]))
}

#' Retain the membership design's coding for prediction
#' @noRd
.mixture_membership_design <- function(formula, data) {
  frame <- stats::model.frame(formula, data, na.action = stats::na.fail)
  terms <- stats::terms(frame)
  design <- stats::model.matrix(terms, frame)
  list(matrix = design, terms = terms,
       xlevels = stats::.getXlevels(terms, frame),
       contrasts = attr(design, "contrasts"))
}

#' Resolve the outcome to counts and trials
#' @param response The model response (vector, factor or two-column matrix).
#' @param family Family name.
#' @return A list of numeric `y` and `trials`.
#' @noRd
.mixture_response <- function(response, family) {
  bad_data <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  whole <- function(v) all(abs(v - round(v)) < sqrt(.Machine$double.eps))
  if (identical(family, "binomial")) {
    if (is.matrix(response)) {
      if (ncol(response) != 2L) {
        bad_data("A binomial matrix outcome must be cbind(successes, failures).")
      }
      y <- as.numeric(response[, 1L])
      trials <- y + as.numeric(response[, 2L])
    } else if (is.factor(response)) {
      if (nlevels(response) != 2L) {
        bad_data("A factor outcome for the binomial family needs two levels.")
      }
      y <- as.numeric(response == levels(response)[2L])
      trials <- rep(1, length(y))
    } else {
      y <- as.numeric(response)
      trials <- rep(1, length(y))
    }
    if (any(y < 0) || any(y > trials) || !whole(y) || !whole(trials)) {
      bad_data(paste("A binomial outcome must be 0/1, logical, a two-level",
                     "factor, or whole-number counts with successes not",
                     "exceeding trials."))
    }
    return(list(y = y, trials = trials))
  }
  if (is.matrix(response) || is.factor(response) || !is.numeric(response)) {
    bad_data(sprintf("The %s family needs a numeric outcome.", family))
  }
  y <- as.numeric(response)
  if (identical(family, "poisson") && (any(y < 0) || !whole(y))) {
    bad_data("A poisson outcome must be non-negative whole numbers.")
  }
  list(y = y, trials = rep(1, length(y)))
}

#' Collapse a row-level design to one row per group, checking constancy
#' @param design Row-level design matrix.
#' @param group_index Integer group of each row.
#' @param argument Argument name, for the message.
#' @return A `G x r` design matrix.
#' @noRd
.mixture_group_design <- function(design, group_index, argument) {
  first <- match(seq_len(max(group_index)), group_index)
  collapsed <- design[first, , drop = FALSE]
  spread <- abs(design - collapsed[group_index, , drop = FALSE])
  if (any(spread > sqrt(.Machine$double.eps) * (1 + abs(design)))) {
    stop(errorCondition(sprintf(paste(
      "`%s` covariates must be constant within each group, because they",
      "predict a group-level class."), argument),
      class = "latents_bad_data", call = NULL))
  }
  rownames(collapsed) <- NULL
  collapsed
}
