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
#' @section Growth mixture models:
#' With repeated measures of persons (`id`, `class_level = "group"`), a time
#' predictor in `formula` gives every class its own trajectory: latent class
#' growth analysis (Nagin 2005). `random` adds random effects *within*
#' classes, so persons scatter around their class trajectory: the growth
#' mixture model (Verbeke and Lesaffre 1996; Muthen and Shedden 1999). Within
#' class k, person i's outcomes are
#' \deqn{y_i = X_i \beta_k + Z_i b_i + e_i, \quad b_i \sim N(0, G_k),
#' \quad e_i \sim N(0, \sigma_k^2 I),}{y_i = X_i beta_k + Z_i b_i + e_i, b_i ~ N(0, G_k), e_i ~ N(0, sigma_k^2 I),}
#' with `Z` built from `random`. The likelihood is exact (Gaussian, no
#' numerical integration); estimation is EM with the random effects as missing
#' data, finished by quasi-Newton with analytic scores. Such a fit has class
#' `latents_growth_mixture`; its tables, plots and inference are described on
#' [get_results.latents_growth_mixture()] and
#' [plot.latents_growth_mixture()]. Random effects need the gaussian family,
#' `id` and `class_level = "group"`.
#'
#' With `cluster`, persons are nested in clusters (students in schools) and
#' the multilevel growth mixture model is fitted: each cluster belongs to a
#' group class h, and the trajectory-class logits of its persons have a
#' group-class-specific intercept and shared covariate slopes,
#' \deqn{L = \prod_j \sum_h \omega_h \prod_{i \in j} \sum_k \pi_{k|h} f_k(y_i)}{L = prod_j sum_h omega_h prod_(i in j) sum_k pi_(k|h) f_k(y_i)}
#' (Asparouhov and Muthen 2008; Vermunt 2003). `group_membership` names
#' cluster-level covariates of the group class. The group classes are read
#' with `get_results(fit, "group_classes")` and `"clusters"`.
#'
#' @section Families:
#' `"gaussian"` (identity link, a residual standard deviation per class, or
#' one shared under `variance = "equal"`), `"binomial"` (logit link; the
#' outcome is 0/1, logical, a two-level factor whose second level is the
#' success, or `cbind(successes, failures)`), `"poisson"` (log link;
#' `offset()` terms in the formula are honoured) and `"negative_binomial"`
#' (log link, NB2: variance mu + alpha mu^2 with a dispersion alpha per class,
#' or one shared under `variance = "equal"`, for overdispersed counts; a
#' dispersion estimated at zero is the Poisson limit and raises
#' `latents_boundary`).
#'
#' `"ordinal"` fits a cumulative-logit (proportional-odds) regression within
#' each class (McCullagh 1980),
#' \deqn{P(Y \le c \mid x, k) = F(t_{kc} - x'\beta_k), \quad c = 1, \ldots, C - 1,}{P(Y <= c | x, k) = plogis(t_kc - x' beta_k), c = 1, ..., C - 1,}
#' with ordered thresholds \eqn{t_{k1} < \ldots < t_{k,C-1}}{t_k1 < ... < t_k,C-1}
#' per class in place of the intercept, and the parametrization of
#' `MASS::polr()`: a positive slope moves the class towards higher
#' categories. The outcome is an ordered factor, a factor (its levels taken
#' in order) or whole-number categories; the categories are those observed,
#' at least two. `common` terms share their slopes across classes, and
#' `offset()` is honoured; the formula's intercept is replaced by the
#' thresholds, so it must not be removed. The thresholds appear in the
#' coefficient table as terms `"threshold:<lower>|<upper>"`, and
#' `exp_estimate` of a slope is the cumulative odds ratio of a higher
#' category. The class "mean" of an ordinal outcome (fitted values,
#' `predict()`, trajectories) is the expected category score
#' \eqn{\sum_c c P(Y = c)}{sum_c c P(Y = c)}, the categories scored 1 to C in
#' order; `predict(type = "probabilities")` gives the category probabilities.
#' With one outcome per class assignment the classes are identified by how
#' continuous predictors shift the category probabilities: a mixture with no
#' predictors, or with two categories (then the logistic mixture of
#' `"binomial"`), is refused with `latents_not_identified` unless
#' `class_level = "group"` gives each class assignment several outcomes.
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
#' @param family `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"negative_binomial"` or `"ordinal"` (proportional-odds cumulative logit;
#'   see *Families*).
#' @param id `NULL`, or the name of the column identifying groups (persons,
#'   schools, ...).
#' @param class_level `"observation"` (each row has its own class) or
#'   `"group"` (all rows of a group share a class; needs `id`).
#' @param n_group_classes Number of second-level group classes, for
#'   `class_level = "observation"` with `id`. `1` fits no second level.
#' @param common `NULL`, or the predictor terms whose coefficients are
#'   shared by every class: variable names such as `"age"`, or a one-sided
#'   formula such as `~ age`. Each term must appear in `formula`.
#' @param membership Covariates predicting class membership (a multinomial
#'   logit, the first class as the reference): variable names such as
#'   `"age"`, or a one-sided formula (`~ 1`, the default, for none).
#'   Row-level for `"observation"`; constant within groups for `"group"`.
#' @param group_membership Group-level covariates predicting the group class,
#'   for the two-level model, as names or a one-sided formula. Must be
#'   constant within groups.
#' @param variance `"varying"` (a residual standard deviation, or a
#'   negative-binomial dispersion, per class) or `"equal"` (one shared).
#'   Gaussian and negative-binomial families only.
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
#' @param random `NULL` (no random effects), or the random effects within
#'   classes: variable names such as `"time"` (a random intercept and a random
#'   slope on `time`) or `"intercept"` (random intercepts only), or a
#'   one-sided formula such as `~ 1 + time` (needed for anything names cannot
#'   say, such as `~ 0 + time`, a slope without a random intercept). See the
#'   growth mixture section.
#' @param random_covariance How the random-effect covariance differs across
#'   classes: `"varying"` (one per class), `"equal"` (one shared matrix, as
#'   `lcmm::hlme()` with `nwg = FALSE`) or `"proportional"` (a shared matrix
#'   times a class-specific scale, the last class's being one; `nwg = TRUE`).
#' @param random_diagonal `TRUE` for uncorrelated random effects (a diagonal
#'   covariance, as `idiag = TRUE` in lcmm).
#' @param cluster `NULL`, or the name of a column grouping the persons (`id`)
#'   into clusters, such as schools, for the multilevel growth mixture model
#'   (with `random`): each cluster belongs to one of `n_group_classes` group
#'   classes, which shifts how probable each trajectory class is for its
#'   persons. The clusters are the independent units of the likelihood, the
#'   standard errors and the BIC; with `n_group_classes = 1` the model is the
#'   single-level growth mixture with clusters as those units, comparable by
#'   BIC with more group classes. See the growth mixture section.
#' @return An object of class `latents_mixture_regression`, or
#'   `latents_growth_mixture` when `random` is given. Read it with
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
#'   assignment, or an ordinal mixture with one outcome per class assignment
#'   and two categories or no predictors; `latents_missing_data` for missing values under
#'   `missing = "error"`; `latents_rows_dropped` (warning) under
#'   `missing = "omit"`; `latents_no_valid_start` when every start
#'   degenerated; `latents_unconverged` (warning) when the selected start did
#'   not converge; `latents_degenerate_start` (warning) when some starts
#'   degenerated; `latents_separation` (warning) when a coefficient diverged.
#'
#' @references
#' Asparouhov, T., & Muthen, B. (2008). Multilevel mixture models. In G. R.
#' Hancock & K. M. Samuelsen (Eds.), *Advances in latent variable mixture
#' models* (pp. 27--51). Information Age.
#'
#' McCullagh, P. (1980). Regression models for ordinal data. *Journal of the
#' Royal Statistical Society B*, 42, 109--142.
#'
#' Muthen, B., & Shedden, K. (1999). Finite mixture modeling with mixture
#' outcomes using the EM algorithm. *Biometrics*, 55, 463--469.
#'
#' Nagin, D. S. (2005). *Group-based modeling of development*. Harvard
#' University Press.
#'
#' Proust-Lima, C., Philipps, V., & Liquet, B. (2017). Estimation of extended
#' mixed models using latent classes and latent processes: the R package
#' lcmm. *Journal of Statistical Software*, 78(2), 1--56.
#'
#' Verbeke, G., & Lesaffre, E. (1996). A linear mixed-effects model with
#' heterogeneity in the random-effects population. *Journal of the American
#' Statistical Association*, 91, 217--221.
#'
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
                   random_diagonal = FALSE, cluster = NULL) {
  vcov_defaulted <- missing(vcov_type)
  random_covariance <- match.arg(random_covariance)
  # Variable names may stand in for a one-sided formula, so no `~` is needed.
  caller <- parent.frame()
  common <- .mixture_as_formula(common, "common", caller)
  membership <- .mixture_as_formula(membership, "membership", caller)
  group_membership <- .mixture_as_formula(group_membership, "group_membership", caller)
  random <- .mixture_as_formula(random, "random", caller)
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
  if (!is.null(cluster) && !is.null(weights)) {
    .latents_refuse_weights("a multilevel growth mixture (`cluster`)")
  }
  spec <- .mixture_spec(formula, data, n_classes, family, id, class_level,
                       n_group_classes, common, membership, group_membership,
                       variance, min_variance, missing, random, cluster)
  if (!is.null(weights)) {
    spec <- .mixture_weight_spec(spec, data, weights)
    if (!identical(vcov_type, "none")) {
      vcov_type <- .latents_weighted_vcov(TRUE, vcov_type, vcov_defaulted)
    }
  }
  if (!is.null(random)) {
    spec$random_covariance <- random_covariance
    spec$random_diagonal <- random_diagonal
    fit <- .mixture_with_seed(seed, .growth_fit(
      spec, as.integer(n_starts), as.integer(max_iter), tol, select_start))
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

#' Read variable names as a one-sided formula
#'
#' `"time"` means `~ 1 + time` and `c("age", "sex")` means `~ 1 + age + sex`;
#' `"intercept"` alone means `~ 1`. Formulas pass through unchanged, so
#' anything names cannot say (no intercept, interactions, transformations) is
#' still written as a formula.
#'
#' @param value `NULL`, a character vector or a formula.
#' @param argument The argument's name, for the message.
#' @param env The environment the formula is evaluated in.
#' @return `NULL` or a one-sided formula; raises `latents_bad_argument`.
#' @noRd
.mixture_as_formula <- function(value, argument, env = parent.frame()) {
  if (is.null(value) || inherits(value, "formula")) return(value)
  if (!is.character(value) || length(value) == 0L || anyNA(value) ||
      !all(nzchar(value))) {
    stop(errorCondition(sprintf(paste(
      "`%s` must be variable names such as \"time\" (or \"intercept\"), or a",
      "one-sided formula such as `~ 1 + time`."), argument),
      class = "latents_bad_argument", call = NULL))
  }
  terms <- setdiff(value, c("intercept", "1"))
  if (length(terms) == 0L) return(stats::as.formula("~ 1", env = env))
  stats::reformulate(sprintf("`%s`", terms), env = env)
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
    x == floor(x) && x <= .Machine$integer.max
}

#' Evaluate code under a seed, restoring the caller's random state
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
#' @return A list consumed by the engine; see the fields assigned below.
#' @noRd
.mixture_spec <- function(formula, data, n_classes, family, id, class_level,
                         n_group_classes, common, membership,
                         group_membership, variance, min_variance, missing,
                         random = NULL, cluster = NULL) {
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
  data <- data[kept_rows, , drop = FALSE]
  n <- nrow(data)

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
      "assignment is not identified: a mixture of Bernoulli distributions is",
      "itself a Bernoulli distribution (Follmann and Lambert 1991). Give each",
      "class assignment several outcomes with `id` and",
      "`class_level = \"group\"`, or supply binomial counts as",
      "`cbind(successes, failures)`."),
      class = "latents_not_identified", call = NULL))
  }
  if (ordinal) {
    .ordinal_check_identified(length(response$levels), n_classes, nesting,
                              length(varying_names) + length(common_names))
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
  clusters <- if (is.null(cluster)) NULL else
    .growth_cluster_spec(data[[cluster]], group_index, row_membership,
                         group_membership_design$matrix, n_group_classes)
  random_design <- if (is.null(random)) NULL else stats::model.matrix(random, data = data)
  if (any(!is.finite(random_design))) {
    stop(errorCondition("The evaluated random-effects design must be finite.",
                        class = "latents_bad_data", call = NULL))
  }
  c(list(
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
    min_variance = min_variance * max(y_variance, .Machine$double.eps),
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
    cluster = cluster), clusters,
    if (ordinal) list(n_categories = length(response$levels),
                      category_levels = response$levels,
                      category_values = response$values,
                      category_kind = response$kind))
}

#' Retain the membership design's coding for prediction
#' @noRd
.mixture_membership_design <- function(formula, data) {
  frame <- stats::model.frame(formula, data, na.action = stats::na.fail)
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
