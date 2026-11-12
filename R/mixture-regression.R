#' Estimate regression relationships that vary across latent classes
#'
#' @description
#' `mixture_regression()` asks whether the relationship between an outcome
#' and its predictors differs across latent classes. The conditional outcome
#' distribution is a mixture of class-specific regressions, so classes can
#' differ in outcome level, in predictor associations or in trajectory over
#' time. Regression parameters and class probabilities are estimated jointly
#' by maximum likelihood.
#'
#' Predictors in `formula` describe variation within classes; covariates in
#' `membership` describe who belongs to which class. A coefficient may be
#' shared by every class (`common`). Class membership can belong to each row
#' or to every row of a unit (`class_level`); with a time predictor and
#' shared membership the classes are trajectories, and `random` adds
#' within-class variation around them (a growth mixture model). Outcomes can
#' be Gaussian, binomial, Poisson, negative binomial or ordinal.
#'
#' @param formula A model formula `outcome ~ predictors`. Factors,
#'   interactions and `offset()` are supported.
#' @param data A data frame.
#' @param n_classes Number of regression classes.
#' @param family `"gaussian"` (default), `"binomial"`, `"poisson"`,
#'   `"negative_binomial"` or `"ordinal"`. See *Families*.
#' @param id `NULL` (default), or the column identifying the units (persons,
#'   schools) the rows belong to.
#' @param class_level `"observation"` (default): each row has its own class.
#'   `"group"`: all rows of a unit share a class (needs `id`). See *Nesting*.
#' @param n_group_classes Number of second-level group classes (default 1,
#'   none). For `class_level = "observation"` with `id`, classes of the
#'   units: a unit's class sets the class probabilities of its rows. For a
#'   growth mixture with `cluster`, classes of the clusters.
#' @param common Predictor terms whose coefficients are shared by every
#'   class, as names (`"age"`) or a one-sided formula (`~ age`). Each term
#'   must appear in `formula`; the intercept cannot be shared.
#' @param membership Covariates of class membership, a multinomial logit
#'   with the first class as reference: names or a one-sided formula (`~ 1`,
#'   the default, for none). Row-level under `"observation"`, constant within
#'   units under `"group"`.
#' @param group_membership Covariates of the group class, as names or a
#'   one-sided formula; constant within the units that carry group classes
#'   (the `id` units, or the clusters under `cluster`).
#' @param variance `"varying"` (default): a residual standard deviation (or
#'   negative-binomial dispersion) per class; `"equal"`: one shared. Gaussian
#'   and negative-binomial families only.
#' @param n_starts Number of random starts (default 10), besides one start
#'   built from the residuals of a pooled regression. Each start runs up to
#'   50 EM iterations; the better half (at least two) continue to
#'   convergence.
#' @param max_iter Maximum EM iterations per start (default 1000).
#' @param tol Convergence tolerance (default `1e-8`), per row: EM stops when
#'   a step changes the log likelihood by no more than `tol` times the number
#'   of rows of positive weight.
#' @param min_variance Floor on a Gaussian residual variance, as a fraction
#'   of the outcome's variance (default `1e-6`). A start in which a class's
#'   residual variance reaches the floor is excluded from selection. In a
#'   growth mixture the same floor bounds each random effect's conditional
#'   variance; a random-effect variance at the floor is kept and warned
#'   about.
#' @param seed Optional integer seed; when given, the caller's random-number
#'   state is restored on exit.
#' @param missing `"error"` (default) rejects rows with a missing value in
#'   any variable the model uses; `"omit"` drops them with a warning saying
#'   how many.
#' @param select_start `"likelihood"` (default): the start with the highest
#'   likelihood; `"converged"`: the best start that converged.
#' @param vcov_type Covariance of the estimates stored with the fit:
#'   `"observed"` (default, inverse observed information), `"robust"`
#'   (sandwich, clustered on the independent units), `"opg"` (outer product
#'   of the scores) or `"none"`. A weighted fit takes `"robust"` or
#'   `"none"` only.
#' @param weights Name of a numeric column holding a sampling weight per
#'   independent unit (per row without `id`, per unit with it): pseudo
#'   maximum likelihood with robust standard errors. A unit of weight zero
#'   is absent from the fit and every criterion but keeps its row and
#'   posterior.
#' @param random `NULL` (default), or the random effects within classes:
#'   names such as `"time"` (a random intercept and a random slope on time)
#'   or `"intercept"`, or a one-sided formula such as `~ 1 + time` (`~ 0 +
#'   time` for a slope alone). See *Growth mixture models*.
#' @param random_covariance How the random-effect covariance differs across
#'   classes: `"varying"` (default, one per class), `"equal"` (one shared
#'   matrix, `lcmm::hlme(nwg = FALSE)`) or `"proportional"` (a shared matrix
#'   times a class-specific scale, the last class's being one, `nwg = TRUE`).
#' @param random_diagonal `TRUE` for uncorrelated random effects (`idiag =
#'   TRUE` in lcmm); default `FALSE`.
#' @param cluster Column grouping the units of `id` into clusters, such as
#'   schools, for the multilevel growth mixture model (with `random`): each
#'   cluster belongs to one of `n_group_classes` group classes, which sets
#'   the trajectory-class probabilities of its persons. The clusters are the
#'   independent units of the likelihood, the standard errors and the BIC.
#'   See *Growth mixture models*.
#' @param residual_correlation `"none"` (default): independent residuals
#'   within a unit. With `random`, `"ar"` adds a continuous-time
#'   autoregressive process (covariance `sd^2 * exp(-rate * |t - s|)`) and
#'   `"brownian"` Brownian motion (`sd^2 * min(t, s)`) to the residuals,
#'   shared by the classes: `lcmm::hlme(cor = AR(time))` and `cor =
#'   BM(time)`. A process the data cannot separate from the residual
#'   variance or the random effects stays on its bound, warns, and has no
#'   standard error. These fits refine one class placement from the model
#'   without the process; `n_best_replicated` is `NA`, and `"starts"` reports
#'   the correlated fit rather than the initialization's replication count.
#' @param residual_time Numeric time column of the residual process; default
#'   the first numeric variable in `random`, otherwise in the predictors.
#'   Brownian motion needs non-negative times.
#'
#' @return A `latents_mixture_regression` object, or `latents_growth_mixture`
#'   when `random` is given. `as.data.frame()` returns the coefficient table
#'   and [get_results()] every table by name; [predict()], `plot()` (and
#'   [plot_trajectories()], [plot_coefficients()]), `summary()`, `coef()`,
#'   `vcov()`, `confint()`, `logLik()` and `nobs()` apply. The within-class
#'   coefficients are conditional associations under the specified model,
#'   not causal effects; check the starts table and the classification
#'   diagnostics before reading substance into the classes.
#'
#' @section Nesting:
#' \describe{
#'   \item{`id = NULL`}{Each row is an independent unit.}
#'   \item{`class_level = "group"`}{All rows of a unit share a class, the
#'     regression mixture for repeated measures (`flexmix(y ~ x | id)`); a
#'     unit's likelihood is the product of its rows' densities.}
#'   \item{`class_level = "observation"` with `id` and
#'     `n_group_classes >= 2`}{Each row has its own class, and a group class
#'     of its unit sets the class probabilities. With
#'     `n_group_classes = 1`, `id` clusters the robust standard errors and
#'     is the unit of `weights`.}
#' }
#'
#' @section Growth mixture models:
#' With repeated measurements (`id`, `class_level = "group"`) and time in
#' `formula`, the classes are trajectories: without random effects a latent
#' class growth model, with `random` a growth mixture model. For unit `i`
#' in class `k`,
#'
#' \deqn{y_i = X_i \beta_k + Z_i b_i + e_i, \quad b_i \sim N(0, G_k), \quad
#'   e_i \sim N(0, \sigma_k^2 I),}{y_i = X_i beta_k + Z_i b_i + e_i,
#'   b_i ~ N(0, G_k), e_i ~ N(0, sigma_k^2 I),}
#'
#' with `Z` the design of `random`. Random effects need the Gaussian family,
#' `id`, `class_level = "group"` and, for a random intercept, repeated rows
#' in at least some units. With `residual_correlation`, `e_i` also carries the
#' autoregressive or Brownian process over `residual_time`.
#'
#' With `cluster`, persons are nested in clusters and each cluster belongs
#' to a group class `h` that sets the trajectory-class logits of its
#' persons (a class-specific intercept with shared covariate slopes):
#'
#' \deqn{L = \prod_j \sum_h \omega_h \prod_{i \in j} \sum_k \pi_{k|h}
#'   f_k(y_i)}{L = prod_j sum_h omega_h prod_(i in j) sum_k pi_(k|h) f_k(y_i)}
#'
#' `group_membership` names
#' cluster-level covariates of the group class; the group classes are in
#' `get_results(fit, "group_classes")` and `"clusters"`. With
#' `n_group_classes = 1` the fit is comparable by BIC with fits that add
#' group classes.
#'
#' @section Families:
#' `"gaussian"`: identity link, a residual standard deviation per class or
#' one shared. `"binomial"`: logit link; the outcome is 0/1, logical, a
#' two-level factor (second level the success) or `cbind(successes,
#' failures)`. `"poisson"`: log link, `offset()` honoured.
#' `"negative_binomial"`: log link, NB2 (variance `mu + alpha * mu^2`) with a
#' dispersion per class or one shared; a dispersion estimated at zero is the
#' Poisson limit, reported with a warning and no standard error.
#'
#' `"ordinal"`: a cumulative-logit (proportional-odds) regression within each
#' class, \eqn{P(Y \le c \mid x, k) = \mathrm{plogis}(t_{kc}
#' - x'\beta_k)}{P(Y <= c | x, k) = plogis(t_kc - x' beta_k)}, with ordered thresholds per class in place of the intercept
#' and the sign convention of `MASS::polr()`: a positive slope raises the
#' probability of a higher category. The outcome is an ordered factor, a
#' factor (levels in order) or whole numbers; the categories are those
#' observed, at least two. The formula keeps its intercept, which the
#' thresholds replace. Thresholds appear in the
#' coefficient table as `"threshold:<lower>|<upper>"`, and `exp_estimate` of a
#' slope is a cumulative odds ratio. The class mean of an ordinal outcome
#' (fitted values, trajectories) is the expected category score with
#' categories scored 1 to C; `predict(type = "probabilities")` gives the
#' category probabilities.
#'
#' Identification. Call the thing that carries a class the unit: a row, or
#' an `id` group under `class_level = "group"`. An ordinal mixture whose
#' units have one outcome each and no predictors is not identified (a
#' mixture of categorical distributions is itself one) and is refused. A
#' binary mixture with one trial per unit, and an ordinal one with one
#' outcome per unit and two categories, are refused too: whether such a
#' model is identified depends on the covariate design, which the package
#' does not assess. Give each unit several
#' outcomes (`class_level = "group"`) or, for binomial, counts with more than
#' one trial per row. A single-level ordinal mixture is also refused when it
#' has more free parameters than the distinct covariate patterns can
#' identify.
#' Repeated binary or ordinal outcomes without predictors must also have
#' enough distinct response-count patterns to identify the class parameters.
#' A growth model is refused when its random-effects and residual covariance
#' components cannot be separated by the observed designs. With class-specific
#' residual variances, two common waves do not identify an unrestricted random
#' intercept and slope covariance. A diagonal covariance can be identified
#' there, as can proportional class covariances with a shared residual variance.
#'
#' Unused factor levels in membership and random-effects formulas are dropped
#' before fitting. In output tables, an identifier whose name conflicts with
#' a result column is prefixed with `id_` or `cluster_` (for example, an
#' `id = "class"` column is returned as `id_class`, beside the latent `class`).
#' A numeric suffix is added if that name is also occupied.
#'
#' Every refusal and warning is a classed condition; [latents-conditions]
#' lists them.
#'
#' @seealso [enumerate_regressions()] to compare numbers of classes;
#'   [plot_trajectories()], [plot_coefficients()] and, for fits that classify
#'   persons, [plot_classification()].
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
#'
#' \donttest{
#' # Students nested in schools, each school in one of two group classes:
#' schools <- mixture_regression(score ~ wave, growth_schools, n_classes = 2,
#'                               id = "student", class_level = "group",
#'                               random = "wave", random_covariance = "equal",
#'                               cluster = "school", n_group_classes = 2,
#'                               n_starts = 2, seed = 1)
#' get_results(schools, "group_classes")
#' }
#' @export
mixture_regression <- function(formula, data, n_classes,
                   family = c("gaussian", "binomial", "poisson", "negative_binomial",
                              "ordinal"),
                   id = NULL, class_level = c("observation", "group"),
                   n_group_classes = 1L, common = NULL, membership = ~1,
                   group_membership = ~1, variance = c("varying", "equal"),
                   n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                   min_variance = 1e-6, seed = NULL,
                   missing = c("error", "omit"),
                   select_start = c("likelihood", "converged"),
                   vcov_type = c("observed", "robust", "opg", "none"),
                   weights = NULL, random = NULL,
                   random_covariance = c("varying", "equal", "proportional"),
                   random_diagonal = FALSE, cluster = NULL,
                   residual_correlation = c("none", "ar", "brownian"),
                   residual_time = NULL) {
  vcov_defaulted <- missing(vcov_type)
  random_covariance <- .latents_match_formal(random_covariance)
  residual_correlation <- .latents_match_formal(residual_correlation)
  # Variable names may stand in for a one-sided formula, so no `~` is needed.
  common <- .mixture_as_formula(common, "common")
  membership <- .mixture_as_formula(membership, "membership")
  group_membership <- .mixture_as_formula(group_membership, "group_membership")
  random <- .mixture_as_formula(random, "random")
  family <- .latents_match_formal(family)
  class_level <- .latents_match_formal(class_level)
  variance <- .latents_match_formal(variance)
  missing <- .latents_match_formal(missing)
  select_start <- .latents_match_formal(select_start)
  vcov_type <- .latents_match_formal(vcov_type)
  .latents_check_data_frame(data, substitute(data))
  .latents_require(
    "`formula` must be a two-sided formula" =
      inherits(formula, "formula") && length(formula) == 3L,
    "`data` must be a data frame" = is.data.frame(data),
    "`n_classes` must be a single positive integer" =
      .mixture_is_count(n_classes),
    "`n_group_classes` must be a single positive integer" =
      .mixture_is_count(n_group_classes),
    "`id` must be NULL or a single column name" =
      is.null(id) || (is.character(id) && length(id) == 1L && !is.na(id)),
    "`cluster` must be NULL or a single column name" =
      is.null(cluster) || (is.character(cluster) && length(cluster) == 1L &&
                             !is.na(cluster)),
    "`common` must be NULL, variable names or a one-sided formula" =
      is.null(common) || (inherits(common, "formula") && length(common) == 2L),
    "`membership` must be variable names or a one-sided formula" =
      inherits(membership, "formula") && length(membership) == 2L,
    "`group_membership` must be variable names or a one-sided formula" =
      inherits(group_membership, "formula") && length(group_membership) == 2L,
    "`n_starts` must be a single non-negative integer" =
      is.numeric(n_starts) && length(n_starts) == 1L && is.finite(n_starts) &&
      n_starts >= 0 && n_starts == floor(n_starts) && n_starts <= .Machine$integer.max,
    "`max_iter` must be a single non-negative integer" =
      is.numeric(max_iter) && length(max_iter) == 1L && is.finite(max_iter) &&
      max_iter >= 0 && max_iter == floor(max_iter) && max_iter <= .Machine$integer.max,
    "`tol` must be a single positive number" =
      is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0,
    "`min_variance` must be a single positive number" =
      is.numeric(min_variance) && length(min_variance) == 1L &&
      is.finite(min_variance) && min_variance > 0,
    "`seed` must be NULL or a single number" =
      is.null(seed) || (is.numeric(seed) && length(seed) == 1L)
  )
  .multilpa_check_seed(seed)
  n_classes <- as.integer(n_classes)
  n_group_classes <- as.integer(n_group_classes)
  if (!is.null(random) || !is.null(cluster)) {
    .growth_check_arguments(random, family, id, class_level, n_group_classes,
                            random_diagonal, cluster)
  }
  if (!identical(residual_correlation, "none")) {
    residual_time <- .growth_residual_time(residual_correlation, residual_time, random,
                                          formula, data)
  } else if (!is.null(residual_time)) {
    stop(errorCondition(paste(
      "`residual_time` is the time scale of correlated residuals: give",
      "`residual_correlation = \"ar\"` or `\"brownian\"` as well."),
      class = "latents_bad_argument", call = NULL))
  }
  if (!is.null(cluster) && !is.null(weights)) {
    .latents_refuse_weights("a multilevel growth mixture (`cluster`)")
  }
  spec <- .mixture_spec(formula, data, n_classes, family, id, class_level,
                       n_group_classes, common, membership, group_membership,
                       variance, min_variance, missing, random, cluster,
                       support = .mixture_row_support(data, weights))
  if (!is.null(weights)) {
    spec <- .mixture_weight_spec(spec, data, weights, min_variance)
    if (!identical(vcov_type, "none")) {
      # Only a growth fit (`random`) has a bootstrap to name.
      vcov_type <- .latents_weighted_vcov(TRUE, vcov_type, vcov_defaulted,
                                          bootstrap = !is.null(random))
    }
  }
  .mixture_check_group_capacity(spec)
  if (!is.null(random)) {
    carrying <- (spec$row_weights %||% rep(1, spec$n)) > 0
    if ("(Intercept)" %in% colnames(spec$random_design) &&
        all(tabulate(spec$group_index[carrying], spec$n_groups) <= 1L)) {
      stop(errorCondition(paste(
        "Every person of positive weight has only one observation. A random",
        "intercept variance and a residual variance are then indistinguishable;",
        "repeated observations are needed for a growth mixture."),
        class = "latents_not_identified", call = NULL))
    }
    spec$random_covariance <- random_covariance
    spec$random_diagonal <- random_diagonal
    # Fitted on a standardized outcome and mapped back exactly, so the search
    # and its stopping rules do not depend on the outcome's units (G2).
    standard <- .growth_standardize(spec)
    .growth_check_covariance_design(standard$spec)
    fit <- .mixture_with_seed(seed, .growth_fit(
      standard$spec, as.integer(n_starts), as.integer(max_iter), tol, select_start))
    if (!identical(residual_correlation, "none")) {
      # The model without the process places the classes; the process is then
      # estimated with everything else (growth-mixture-correlated.R).
      times <- data[[residual_time]][spec$kept_rows]
      if (anyNA(times) || any(!is.finite(times))) {
        stop(errorCondition("`residual_time` must be finite on every fitted row.",
                            class = "latents_bad_data", call = NULL))
      }
      process <- list(residual_correlation = residual_correlation,
                      residual_time = residual_time, residual_times = times)
      spec[names(process)] <- process
      spec$model_data[[residual_time]] <- times
      standard$spec[names(process)] <- process
      standard$spec$model_data[[residual_time]] <- times
      fit <- .growth_correlated_refit(standard$spec, fit, as.integer(max_iter), tol)
    }
    fit <- .growth_unstandardize_fit(fit, spec, standard)
    fit$call <- match.call()
    fit$settings <- list(n_starts = as.integer(n_starts),
                         max_iter = as.integer(max_iter), tol = tol,
                         min_variance = min_variance, seed = seed,
                         select_start = select_start)
    fit$weights <- weights
    fit$sampling_weights <- spec$sampling_weights
    if (!identical(vcov_type, "none")) {
      fit$inference <- .growth_inference(fit, vcov_type)
    }
    return(fit)
  }
  # A unit of weight zero is absent: every start, and the selection among
  # them, sees the units of positive weight only (as the growth mixture does).
  # Its posterior is then restored from the selected estimates.
  carrying <- .mixture_carrying_spec(spec)
  results <- .mixture_with_seed(seed, .mixture_run_starts(
    carrying, as.integer(n_starts), as.integer(max_iter), tol))
  fit <- .mixture_select(carrying, results, select_start,
                        tie = tol * max(1, .mixture_observed(carrying)))
  if (carrying$n < spec$n) {
    fit$spec <- spec
    fit$expectation <- .mixture_expectation(spec, fit$params, weighted = FALSE)
  }
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

#' Read variable names as a one-sided formula
#'
#' `"time"` means `~ 1 + time` and `c("age", "sex")` means `~ 1 + age + sex`;
#' `"intercept"` alone means `~ 1`. Formulas pass through unchanged, so
#' anything names cannot say (no intercept, interactions, transformations) is
#' still written as a formula.
#'
#' A formula built from names, or one without variables (the default
#' `~ 1`), gets the base environment: its variables all come from the data,
#' and the environment it was made in -- the fitting function's own frame,
#' holding the data and every start's result -- would otherwise be stored
#' (and serialized) with the fit.
#'
#' @param value `NULL`, a character vector or a formula.
#' @param argument The argument's name, for the message.
#' @return `NULL` or a one-sided formula; raises `latents_bad_argument`.
#' @noRd
.mixture_as_formula <- function(value, argument) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "formula")) {
    if (length(all.vars(value)) == 0L) environment(value) <- baseenv()
    return(value)
  }
  if (!is.character(value) || length(value) == 0L || anyNA(value) ||
      !all(nzchar(value))) {
    stop(errorCondition(sprintf(paste(
      "`%s` must be variable names such as \"time\" (or \"intercept\"), or a",
      "one-sided formula such as `~ 1 + time`."), argument),
      class = "latents_bad_argument", call = NULL))
  }
  terms <- setdiff(value, c("intercept", "1"))
  if (length(terms) == 0L) return(stats::as.formula("~ 1", env = baseenv()))
  stats::reformulate(sprintf("`%s`", terms), env = baseenv())
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
#' @param min_variance Optional relative Gaussian variance floor, recomputed
#'   on the rows of positive weight rather than on absent units.
#' @return `spec` with `sampling_weights` (per unit), `row_weights` and
#'   `weights`.
#' @noRd
.mixture_weight_spec <- function(spec, data, weights, min_variance = NULL) {
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
  if (identical(spec$family, "gaussian") && !is.null(min_variance)) {
    carrying <- spec$row_weights > 0
    reference <- .mixture_gaussian_reference(
      spec$y[carrying], spec$offset[carrying],
      cbind(spec$x, spec$z)[carrying, , drop = FALSE])
    spec$min_variance <- min_variance * reference
  }
  spec
}

#' Reference variance and constant-response check for Gaussian regression
#'
#' A constant raw response need not have zero residual variance: an offset or
#' a design without an intercept can leave nonzero residuals. Refuse only a
#' constant adjusted response the design can fit identically. For a constant
#' raw response the adjusted spread (or its squared level without an intercept)
#' supplies a positive, outcome-unit-aware variance reference.
#' @param y,offset Response and fixed offset on rows of positive weight.
#' @param design The full regression design on those rows.
#' @return A positive variance reference in outcome-squared units.
#' @noRd
.mixture_gaussian_reference <- function(y, offset, design) {
  stopifnot(is.numeric(y), length(y) > 0L, length(offset) == length(y),
            is.matrix(design), nrow(design) == length(y))
  adjusted <- y - offset
  constant <- max(adjusted) == min(adjusted)
  if (constant && ("(Intercept)" %in% colnames(design) || adjusted[[1L]] == 0)) {
    stop(errorCondition(paste(
      "The outcome after its offset takes the same value in every row of positive weight,",
      "and the regression fits that constant exactly, leaving no residual variation."),
      class = "latents_bad_data", call = NULL))
  }
  reference <- stats::var(y)
  if (!is.finite(reference) || reference <= 0) {
    reference <- if (constant) mean(adjusted^2) else stats::var(adjusted)
  }
  max(reference, .Machine$double.xmin)
}

#' Is a value a single positive whole number?
#' @noRd
.mixture_is_count <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 &&
    x == floor(x) && x <= .Machine$integer.max
}

#' Evaluate code under a seed, restoring the caller's exposed RNG seed
#'
#' The expression form of `.latents_local_seed()`, for callers that scope a
#' seed to one expression rather than to their whole body.
#' @param seed `NULL` or a number.
#' @param code Expression, evaluated lazily.
#' @return The value of `code`.
#' @noRd
.mixture_with_seed <- function(seed, code) {
  .latents_local_seed(seed, after = FALSE)
  code
}

#' Build the model specification from a formula and data
#'
#' Resolves the outcome, the class-specific and shared designs, the offset,
#' the membership designs and the group index, after one complete-case pass
#' over every variable the model uses.
#'
#' @param support `NULL`, or one logical per row of `data`: `FALSE` for a row
#'   of sampling weight zero. Factor levels and ordinal categories are then
#'   those of the rows of positive weight, and a zero-weight row carrying a
#'   level none of them has is left out (`latents_rows_dropped`): no weighted
#'   observation informs its parameter, which would otherwise be counted in
#'   the BIC and make the robust covariance rank deficient.
#' @return A list consumed by the engine; see the fields assigned below.
#' @noRd
.mixture_spec <- function(formula, data, n_classes, family, id, class_level,
                         n_group_classes, common, membership,
                         group_membership, variance, min_variance, missing,
                         random = NULL, cluster = NULL, support = NULL) {
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
  if (n_group_classes > 1L && identical(class_level, "group") && is.null(cluster)) {
    bad_argument(paste(
      "Group classes shift the class shares of rows within a group; with",
      "`class_level = \"group\"` a group already has a single class. Use",
      "`class_level = \"observation\"` for the two-level model."))
  }
  if (!identical(variance, "varying") &&
      !family %in% c("gaussian", "negative_binomial")) {
    bad_argument(paste("`variance` applies to the gaussian and negative_binomial",
                       "families only."))
  }
  if (!is.null(id) && !id %in% names(data)) {
    bad_argument(sprintf("`id` names column `%s`, which `data` does not have.",
                         id))
  }
  if (!identical(nesting, "two-level") && is.null(cluster) &&
      length(all.vars(group_membership)) > 0L) {
    bad_argument(paste("`group_membership` applies to the two-level model",
                       "only (`id` with `n_group_classes >= 2`)."))
  }
  if (!isTRUE(attr(stats::terms(membership), "intercept") == 1L)) {
    bad_argument("`membership` must keep its intercept.")
  }

  used <- unique(c(all.vars(formula), all.vars(common %||% ~1),
                   all.vars(membership), all.vars(group_membership),
                   all.vars(random %||% ~1), id, cluster))
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
  formulas <- list(formula, membership, group_membership, random)
  outside <- if (is.null(support)) rep(FALSE, length(kept_rows)) else
    .mixture_outside_support(data[kept_rows, , drop = FALSE],
                             support[kept_rows], formulas,
                             ordinal = identical(family, "ordinal"))
  kept_rows <- kept_rows[!outside]
  data <- data[kept_rows, , drop = FALSE]
  n <- nrow(data)
  # Refuse a one-level factor before announcing the dropped rows, so a
  # subgroup selected by zero weights is told the real reason only.
  .mixture_require_levels(data, formulas, weighted = !is.null(support))
  if (any(outside)) {
    warning(warningCondition(sprintf(paste(
      "Left out %d row%s of weight zero carrying a factor level or outcome",
      "category that no row of positive weight has: the model has no",
      "parameter for it."), sum(outside), if (sum(outside) > 1L) "s" else ""),
      class = "latents_rows_dropped", call = NULL))
  }

  frame <- stats::model.frame(formula, data = data, na.action = stats::na.fail,
                              drop.unused.levels = TRUE)
  model_terms <- stats::terms(frame)
  full_design <- stats::model.matrix(model_terms, frame)
  offset <- stats::model.offset(frame) %||% rep(0, n)
  if (any(!is.finite(full_design)) || any(!is.finite(offset))) {
    stop(errorCondition("The evaluated regression design and offsets must be finite.",
                        class = "latents_bad_data", call = NULL))
  }
  response <- .mixture_response(stats::model.response(frame), family)
  ordinal <- identical(family, "ordinal")
  if (ordinal && !isTRUE(attr(model_terms, "intercept") == 1L)) {
    bad_argument(paste("An ordinal regression's thresholds take the intercept's",
                       "place, so `formula` must keep its intercept (no `- 1`",
                       "or `0 +`)."))
  }

  common_names <- if (is.null(common)) character() else
    .mixture_common_columns(common, model_terms, full_design)
  # The thresholds stand in for an ordinal regression's intercept, and they
  # always vary across classes.
  varying_names <- setdiff(colnames(full_design),
                           c(common_names, if (ordinal) "(Intercept)"))
  if (length(varying_names) == 0L && !ordinal) {
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

  if (n_classes > 1L && identical(family, "binomial") && all(response$trials == 1) &&
      !identical(nesting, "group")) {
    stop(errorCondition(paste(
      "A mixture of logistic regressions with one binary outcome per class",
      "assignment is not supported. Identification depends on the covariate",
      "design and model constraints. Give each",
      "class assignment several outcomes with `id` and",
      "`class_level = \"group\"`, or supply binomial counts as",
      "`cbind(successes, failures)`."),
      class = "latents_not_identified", call = NULL))
  }
  if (ordinal) {
    .ordinal_check_identified(length(response$levels), n_classes, nesting,
                              length(varying_names) + length(common_names))
  }
  # Classes sit on rows except under `class_level = "group"`; group classes
  # (two-level) sit on groups and are checked below.
  units <- if (identical(nesting, "group")) n_groups else n
  if (units <= n_classes) {
    stop(errorCondition(sprintf(paste(
      "%d %s cannot support %d classes."), units,
      if (identical(nesting, "group")) "groups" else "rows", n_classes),
      class = "latents_bad_data", call = NULL))
  }
  if (identical(nesting, "two-level") && n_groups <= n_group_classes) {
    stop(errorCondition(sprintf("%d groups cannot support %d group classes.",
                                n_groups, n_group_classes),
                        class = "latents_bad_data", call = NULL))
  }

  y_variance <- if (identical(family, "gaussian"))
    .mixture_gaussian_reference(response$y, offset, full_design) else 1
  clusters <- if (is.null(cluster)) NULL else
    .growth_cluster_spec(data[[cluster]], group_index, row_membership,
                         group_membership_design$matrix, n_group_classes)
  random_frame <- if (is.null(random)) NULL else stats::model.frame(
    random, data = data, na.action = stats::na.fail, drop.unused.levels = TRUE)
  random_design <- if (is.null(random)) NULL else
    stats::model.matrix(stats::terms(random_frame), random_frame)
  if (any(!is.finite(random_design))) {
    stop(errorCondition("The evaluated random-effects design must be finite.",
                        class = "latents_bad_data", call = NULL))
  }
  spec <- c(list(
    family = family, nesting = nesting, variance = variance,
    n = n, n_classes = n_classes,
    n_group_classes = if (identical(nesting, "two-level")) n_group_classes else 1L,
    y = response$y, trials = response$trials,
    log_normalizer = switch(family,
      gaussian = rep(0, n),
      binomial = lchoose(response$trials, response$y),
      poisson = ,
      negative_binomial = -lgamma(response$y + 1),
      ordinal = rep(0, n)),
    x = x, z = z, offset = as.numeric(offset),
    w = w, v = v,
    intercept_only = ncol(w) == 1L && identical(colnames(w), "(Intercept)"),
    group_intercept_only = ncol(v) == 1L &&
      identical(colnames(v), "(Intercept)"),
    id = id, group_index = group_index, group_levels = group_levels,
    n_groups = n_groups,
    min_variance = min_variance * max(y_variance, .Machine$double.xmin),
    kept_rows = kept_rows,
    terms = model_terms, xlevels = stats::.getXlevels(model_terms, frame),
    contrasts = attr(full_design, "contrasts"), model_data = data[used],
    formula = formula, common = common, membership = membership,
    membership_design = membership_design,
    group_membership_design = group_membership_design,
    group_membership = group_membership, response_name =
      deparse1(formula[[2L]]),
    random = random,
    random_design = random_design,
    random_terms = if (is.null(random)) NULL else stats::terms(random_frame),
    random_xlevels = if (is.null(random)) NULL else
      stats::.getXlevels(stats::terms(random_frame), random_frame),
    cluster = cluster), clusters,
    if (ordinal) list(n_categories = length(response$levels),
                      category_levels = response$levels,
                      category_values = response$values,
                      category_kind = response$kind))
  if (ordinal) .ordinal_check_patterns(spec)
  spec
}

#' Refuse repeated categorical mixtures with too few response patterns
#'
#' With no predictors and constant offsets, repeated observations of a person
#' are exchangeable within class. Their count vector determines the person's
#' likelihood, so it bounds the number of identifiable parameters.
#' @param spec The regression specification.
#' @return `NULL`, invisibly; raises `latents_not_identified`.
#' @noRd
.mixture_check_group_capacity <- function(spec) {
  stopifnot(is.list(spec))
  if (!identical(spec$nesting, "group") || spec$n_classes < 2L ||
      !spec$family %in% c("binomial", "ordinal") ||
      !is.null(spec$random) || !is.null(spec$cluster) ||
      !spec$intercept_only || ncol(spec$z) > 0L ||
      any(colnames(spec$x) != "(Intercept)") ||
      length(unique(spec$offset)) != 1L) return(invisible(NULL))
  groups <- spec$group_index[(spec$row_weights %||% rep(1, spec$n)) > 0]
  if (identical(spec$family, "binomial")) {
    trials <- spec$trials[(spec$row_weights %||% rep(1, spec$n)) > 0]
    capacity <- max(rowsum(trials, groups, reorder = FALSE))
  } else {
    largest <- max(tabulate(groups, spec$n_groups))
    capacity <- choose(largest + spec$n_categories - 1L,
                       spec$n_categories - 1L) - 1
  }
  parameters <- .mixture_count_parameters(spec)
  if (capacity < parameters) {
    stop(errorCondition(sprintf(paste(
      "The repeated %s outcomes without predictors identify at most %g",
      "response-count parameters, but this model has %d. Use fewer classes",
      "or a design with more repeated outcomes."), spec$family, capacity,
      parameters), class = "latents_not_identified", call = NULL))
  }
  invisible(NULL)
}

#' Which rows of a weighted fit carry positive weight
#'
#' Read before the weights are validated, so that only an exact zero marks a
#' row as absent; a missing or negative weight is left for
#' `.latents_sampling_weights()` to refuse.
#' @param data The data frame passed to the fit.
#' @param weights `NULL`, or the weight argument as given.
#' @return `NULL` when no row has weight zero, otherwise one logical per row
#'   of `data`, `FALSE` for a row of weight zero.
#' @noRd
.mixture_row_support <- function(data, weights) {
  if (!is.character(weights) || length(weights) != 1L || is.na(weights) ||
      !weights %in% names(data) || !is.numeric(data[[weights]])) {
    return(NULL)
  }
  value <- data[[weights]]
  supported <- is.na(value) | value != 0
  if (all(supported)) NULL else supported
}

#' Zero-weight rows whose categorical values the rows of positive weight lack
#'
#' Checks every factor, character or logical variable of the model frames
#' (the predictors, the membership covariates and the random effects), and an
#' ordinal outcome.
#' @param data The complete rows.
#' @param carrying One logical per row, `FALSE` for weight zero.
#' @param formulas The model's formulas; `NULL` entries are skipped.
#' @param ordinal Whether the outcome is ordinal (its categories then count).
#' @return One logical per row, `TRUE` for a row to leave out.
#' @noRd
.mixture_outside_support <- function(data, carrying, formulas, ordinal) {
  none <- rep(FALSE, nrow(data))
  # With no row of positive weight there is no support to restrict to: the
  # weights check refuses the fit (`latents_bad_weights`) once the design is
  # built, as it did before zero-weight rows were dropped.
  if (all(carrying) || !any(carrying)) return(none)
  columns <- .mixture_categorical_columns(data, formulas, ordinal = ordinal,
                                          logical = TRUE)
  Reduce(`|`, lapply(columns, function(column) {
    !carrying & !(column %in% column[carrying])
  }), none)
}

#' The categorical variables of the model's formulas
#'
#' @param data The rows to evaluate the formulas on.
#' @param formulas The model's formulas; `NULL` entries are skipped.
#' @param ordinal Whether an ordinal outcome counts (otherwise the outcome is
#'   left out).
#' @param logical Whether logical variables count.
#' @return A named list of the factor and character (and, with `logical`,
#'   logical) columns of the model frames, named by the variable.
#' @noRd
.mixture_categorical_columns <- function(data, formulas, ordinal = FALSE,
                                         logical = FALSE) {
  unlist(lapply(Filter(Negate(is.null), formulas), function(formula) {
    frame <- stats::model.frame(formula, data, na.action = stats::na.pass)
    categorical <- vapply(frame, function(column) {
      is.null(dim(column)) && (is.factor(column) || is.character(column) ||
                                 (logical && is.logical(column)))
    }, logical(1))
    response <- attr(attr(frame, "terms"), "response")
    if (isTRUE(response > 0L)) categorical[response] <- ordinal
    as.list(frame[categorical])
  }), recursive = FALSE)
}

#' Refuse a factor predictor with a single observed level
#'
#' A factor or character predictor (in the regression, membership or
#' random-effects formula) needs two distinct values among the rows fitted:
#' with one there is nothing to contrast, and base R's model matrix stops
#' with an unclassed error. With sampling weights the fitted rows are those
#' of positive weight (and any zero-weight row sharing their levels), so a
#' subgroup selected by zero weights whose indicator is in the formula lands
#' here.
#' @param data The fitted rows.
#' @param formulas The model's formulas; `NULL` entries are skipped.
#' @param weighted Whether rows of weight zero were left out first (only
#'   changes the message).
#' @return `NULL`, invisibly; raises `latents_bad_data`.
#' @noRd
.mixture_require_levels <- function(data, formulas, weighted = FALSE) {
  columns <- .mixture_categorical_columns(data, formulas)
  single <- vapply(columns, function(column) {
    length(unique(as.character(column))) < 2L
  }, logical(1))
  if (any(single)) {
    offending <- unique(names(columns)[single])
    several <- length(offending) > 1L
    stop(errorCondition(sprintf(paste(
      "%s %s a single value among the %s rows, so %s nothing to contrast.",
      "Drop %s from the model or fit the subgroup without it."),
      paste(sprintf("`%s`", offending), collapse = ", "),
      if (several) "take" else "takes",
      if (weighted) "fitted (positive-weight)" else "fitted",
      if (several) "they have" else "it has",
      if (several) "them" else "it"),
      class = "latents_bad_data", call = NULL))
  }
  invisible(NULL)
}

#' The design columns of the terms shared by every class
#'
#' A `common` term is matched to a term of `formula` by the variables it
#' involves (so `~ b:a` matches `a:b`), and contributes every design column
#' that term produces in `formula`'s own coding: a factor shares its contrast
#' columns, not the full dummy set a formula without intercept would build.
#' Columns follow the order of the terms in `common`.
#'
#' @param common A one-sided formula.
#' @param model_terms The terms of the fitted formula.
#' @param full_design Its model matrix (with the `assign` attribute).
#' @return A character vector of column names of `full_design`; raises
#'   `latents_bad_argument` for a term `formula` does not have, or for a
#'   `common` with no term at all (a shared intercept).
#' @noRd
.mixture_common_columns <- function(common, model_terms, full_design) {
  term_keys <- function(terms) {
    factors <- attr(terms, "factors")
    if (length(factors) == 0L) return(character())
    vapply(seq_len(ncol(factors)), function(j) {
      paste(sort(rownames(factors)[factors[, j] > 0]), collapse = "\r")
    }, character(1))
  }
  common_terms <- stats::terms(common)
  labels <- attr(common_terms, "term.labels")
  if (length(labels) == 0L) {
    # `"intercept"` and `~ 1` used to be accepted and to change nothing: the
    # intercept stayed class-specific without a word (REG-8).
    stop(errorCondition(paste(
      "`common` names no predictor term. An intercept shared by every class",
      "is not supported (`\"intercept\"` and `~ 1` ask for one); give the",
      "predictor terms whose coefficients the classes share, or leave",
      "`common = NULL`."), class = "latents_bad_argument", call = NULL))
  }
  position <- match(term_keys(common_terms), term_keys(model_terms))
  if (anyNA(position)) {
    stop(errorCondition(sprintf(paste(
      "`common` names %s, which `formula` does not produce; every shared",
      "term must also be in `formula`."),
      paste(sprintf("`%s`", labels[is.na(position)]), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  assign <- attr(full_design, "assign")
  columns <- unlist(lapply(position, function(term) which(assign == term)))
  colnames(full_design)[as.integer(columns)]
}

#' Retain the membership design's coding for prediction
#' @noRd
.mixture_membership_design <- function(formula, data) {
  frame <- stats::model.frame(formula, data, na.action = stats::na.fail,
                              drop.unused.levels = TRUE)
  terms <- stats::terms(frame)
  design <- stats::model.matrix(terms, frame)
  if (any(!is.finite(design))) {
    stop(errorCondition("The evaluated membership design must be finite.",
                        class = "latents_bad_data", call = NULL))
  }
  list(matrix = design, terms = terms,
       xlevels = stats::.getXlevels(terms, frame),
       contrasts = attr(design, "contrasts"))
}

#' Resolve the outcome to counts and trials
#' @param response The model response (vector, factor or two-column matrix).
#' @param family Family name.
#' @param levels For the ordinal family on new data: the fitted categories.
#' @return A list of numeric `y` and `trials` (and, for the ordinal family,
#'   the categories; see `.ordinal_response()`).
#' @noRd
.mixture_response <- function(response, family, levels = NULL) {
  if (identical(family, "ordinal")) return(.ordinal_response(response, levels))
  bad_data <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  whole <- function(v) all(abs(v - round(v)) < sqrt(.Machine$double.eps))
  if (identical(family, "binomial")) {
    if (!is.numeric(response) && !is.logical(response) && !is.factor(response)) {
      bad_data("A binomial outcome must be numeric, logical or a two-level factor.")
    }
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
    if (any(!is.finite(y)) || any(!is.finite(trials)) ||
        any(y < 0) || any(y > trials) || !whole(y) || !whole(trials)) {
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
  if (any(!is.finite(y))) bad_data("The outcome must contain only finite values.")
  if (family %in% c("poisson", "negative_binomial") && (any(y < 0) || !whole(y))) {
    bad_data(sprintf("A %s outcome must be non-negative whole numbers.", family))
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
