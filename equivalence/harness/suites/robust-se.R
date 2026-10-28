# Robust (CR1 sandwich) and outer-product standard errors of Gaussian mixtures,
# and every column of the table parameter_inference() returns (EQ-13).
#
# Before this suite the robust path was compared only for two-level categorical
# LCA (Mplus MLR) and mixture regression (flexmix), and `vcov_type = "opg"`
# only for categorical models (poLCA, glca). The Gaussian path, the package's
# default family, had no robust or opg reference, and parameter_inference()'s
# Wald columns (statistic, p_value, p_adjusted, conf_low, conf_high) are not
# in get_results(fit, "all"), so no ledger ever counted them.
#
# Reference (robust-se-reference.R; nothing there calls the package):
#
#   - The likelihood is written from its definition in base R: Cholesky
#     Gaussian log densities, profile probabilities within each group class
#     (or a multinomial logit in the row covariates, last profile the
#     baseline), group-class probabilities (or a multinomial logit in the group
#     covariates, last group class the baseline), summed over group classes
#     per group (Vermunt 2003). It is parameterised by the free natural
#     parameters (means, variances or covariance elements, the first K - 1
#     probabilities of each block, or the logits), not the package's
#     log-variance / log-Cholesky / logit chart. At a maximum the observed,
#     sandwich and outer-product covariances all transform with the Jacobian
#     of the reparameterisation, so the natural-scale standard errors are the
#     same in either chart.
#   - Per-group scores are numDeriv Jacobians of the group log-likelihoods and
#     the information is minus the numDeriv Hessian of their sum (Richardson
#     extrapolation). observed = A^-1; robust = A^-1 B A^-1 G / (G - 1), the
#     CR1 factor of sandwich::vcovCL() that ?parameter_inference documents;
#     opg = B^-1 (BHHH); the MLR scaling correction tr(A^-1 B) / q.
#   - The arithmetic is calibrated first on a one-profile Gaussian with 40
#     clusters, where the robust errors of the mean and of the variance are
#     exactly sandwich::vcovCL(type = "HC0", cadjust = TRUE) of lm(y ~ 1) and
#     of lm(e^2 ~ 1): the information is block diagonal at the maximum, so
#     each sandwich entry is the cluster sum of squared score sums over the
#     squared information, which is what vcovCL computes for an intercept.
#   - The reference point is the reference likelihood's own maximum, reached
#     by Newton steps (score below 1e-6 standard errors). The faithful case
#     starts from mclust's EM continued to its fixed point (the mclust suite's
#     refinement), so its estimates are another program's; the simulated cases
#     start from the package's estimates, which only selects the mode and the
#     class labels: at least one Newton step is always taken, so the point
#     compared is where the reference's own score vanishes, never the
#     package's estimates themselves, and the `estimate` column checks the
#     package lands on it. A row per fit reports that score, in standard
#     errors, at the point used.
#
# What is compared:
#
#   - Every column of every parameter_inference() table, through the verb
#     ledger: `estimate` and `standard_error` against the reference above;
#     `statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high` against
#     the definitions in ?parameter_inference (Wald z, two-sided normal p,
#     stats::p.adjust() over the tests the table carries, logit intervals for
#     probabilities, log intervals for variances and covariance diagonals,
#     symmetric otherwise) applied to the package's own estimate and standard
#     error, so that arithmetic is checked to rounding while the standard
#     errors are checked against their independent reference.
#   - The `level`, `adjust` and `vcov_type` attributes against the arguments,
#     and `scaling_correction` (where the table carries it) against the
#     reference's tr(A^-1 B) / q, as extra rows.
#   - Every output of every fit, through compare_fit(), with the reference's
#     posteriors, group posteriors, group log-likelihoods, row entropies
#     within group classes and observed-information standard errors.
#
# Cases: the one-profile calibration (robust); mclust's faithful data under
# VVV with two components, one row per unit (observed, robust and opg); a
# simulated two-level Gaussian with two profiles and two group classes, 40
# groups of 12 (robust and opg); and a simulated membership regression with a
# row covariate, at two levels with a group covariate and at one group class
# with the 40 groups as clusters (robust). Levels and p-value adjustments vary
# across the calls so that each documented option is exercised.

source(file.path("equivalence", "harness", "robust-se-reference.R"), local = environment())

#' Tolerances for parameter_inference() tables
#'
#' Estimates: two maxima of the same likelihood, both carried to a score far
#' below their standard errors, the default 1e-6. Standard errors: numDeriv's
#' Richardson Hessian is good to about 1e-7 relative, and the package's
#' information is a central difference of its analytic score with step 1e-4,
#' 1e-4 as the harness uses for every numerical standard error. The Wald
#' columns are computed from the package's own estimate and standard error,
#' the same arithmetic on the same inputs, so they agree to rounding: 1e-10.
.robust_se_tolerances <- c(
  default = 1e-6, `parameter_inference$standard_error` = 1e-4,
  `parameter_inference$statistic` = 1e-10, `parameter_inference$p_value` = 1e-10,
  `parameter_inference$p_adjusted` = 1e-10, `parameter_inference$conf_low` = 1e-10,
  `parameter_inference$conf_high` = 1e-10)

suite_robust_se <- function() {
  require_suite_packages(
    "numDeriv", "sandwich", "mclust",
    reason = paste("the reference scores and information (numDeriv), the cluster-robust",
                   "calibration (sandwich::vcovCL) and the faithful reference estimates (mclust)"))
  # Mclust() resolves its model-search function by name, so the namespace has
  # to be attached rather than only loaded.
  if (!"package:mclust" %in% search()) {
    suppressMessages(library(mclust))
    on.exit(detach("package:mclust", unload = FALSE), add = TRUE)
  }
  two_level <- .robust_se_two_level_data()
  covariates <- .robust_se_covariate_data()
  stack_cases(c(
    .robust_se_calibration(two_level),
    .robust_se_faithful(),
    .robust_se_case(
      "two-level Gaussian, 2 profiles x 2 group classes, 40 groups",
      fit = multilpa(two_level, c("y1", "y2"), "school", n_profiles = 2L,
                     n_group_classes = 2L, n_starts = 5L, seed = 1L, tol = 1e-15,
                     max_iter = 20000L),
      model = rse_model(c("y1", "y2"), 2L, 2L),
      data = rse_data(two_level, c("y1", "y2"), "school"),
      frame = two_level[c("y1", "y2")], groups = two_level$school,
      calls = list(list(vcov_type = "robust", level = 0.95, adjust = "bonferroni"),
                   list(vcov_type = "opg", level = 0.95, adjust = "hochberg"))),
    .robust_se_case(
      "covariates, row x and group z, 2 profiles x 2 group classes, 40 groups",
      fit = multilpa(covariates, c("y1", "y2"), "school", n_profiles = 2L,
                     n_group_classes = 2L, profile_covariates = "x",
                     group_covariates = "z", n_starts = 5L, seed = 1L, tol = 1e-15,
                     max_iter = 20000L),
      model = rse_model(c("y1", "y2"), 2L, 2L, profile_covariates = "x",
                        group_covariates = "z"),
      data = rse_data(covariates, c("y1", "y2"), "school", "x", "z"),
      frame = covariates[c("y1", "y2")], groups = covariates$school,
      calls = list(list(vcov_type = "robust", level = 0.9, adjust = "BY"))),
    .robust_se_case(
      "covariates, row x, one group class, 40 clusters",
      fit = multilpa(covariates, c("y1", "y2"), "school", n_profiles = 2L,
                     n_group_classes = 1L, profile_covariates = "x", n_starts = 5L,
                     seed = 1L, tol = 1e-15, max_iter = 20000L),
      model = rse_model(c("y1", "y2"), 2L, 1L, profile_covariates = "x"),
      data = rse_data(covariates, c("y1", "y2"), "school", "x"),
      frame = covariates[c("y1", "y2")], groups = covariates$school,
      calls = list(list(vcov_type = "robust", level = 0.95, adjust = "hommel")))))
}

#' Run code under a fixed seed and restore the caller's stream
#' @param seed A whole number.
#' @param code Code to evaluate.
#' @return The value of `code`.
.robust_se_with_seed <- function(seed, code) {
  saved <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(saved)) {
    rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", saved, envir = globalenv()), add = TRUE, after = FALSE)
  set.seed(seed)
  code
}

#' Two Gaussian profiles mixed differently in two kinds of group
#'
#' 40 groups of 12; group classes alternate, with profile shares 0.8 / 0.2
#' and 0.2 / 0.8, so the group classes are genuinely there. Variances differ by
#' profile and indicator.
#' @return A `data.frame` with `school`, `y1` and `y2`.
.robust_se_two_level_data <- function() {
  .robust_se_with_seed(13L, {
    school <- rep(seq_len(40L), each = 12L)
    group_class <- rep(rep(1:2, length.out = 40L), each = 12L)
    first <- stats::runif(480L) < ifelse(group_class == 1L, 0.8, 0.2)
    data.frame(school = school,
               y1 = stats::rnorm(480L, ifelse(first, -1.5, 1.5), ifelse(first, 1, 0.7)),
               y2 = stats::rnorm(480L, ifelse(first, -1, 1), ifelse(first, 0.8, 1.2)))
  })
}

#' Membership driven by a row covariate and a group's class
#'
#' 40 groups of 12. A group's class follows a logit in the group covariate
#' `z` (0.5 + z); a row's profile follows -1 + 2 [group class 2] + 0.8 x.
#' @return A `data.frame` with `school`, `x`, `z`, `y1` and `y2`.
.robust_se_covariate_data <- function() {
  .robust_se_with_seed(29L, {
    school <- rep(seq_len(40L), each = 12L)
    z_group <- stats::rnorm(40L)
    second_class <- stats::runif(40L) >= stats::plogis(0.5 + z_group)
    x <- stats::rnorm(480L)
    second <- stats::runif(480L) <
      stats::plogis(-1 + 2 * rep(second_class, each = 12L) + 0.8 * x)
    data.frame(school = school, x = x, z = rep(z_group, each = 12L),
               y1 = stats::rnorm(480L, ifelse(second, 2, -2), 0.8),
               y2 = stats::rnorm(480L, ifelse(second, 1.5, -1.5), 0.8))
  })
}

#' Compare one fit and its parameter_inference() tables with the reference
#'
#' @param label Label of the fit.
#' @param fit The fitted model.
#' @param model The matching `rse_model`.
#' @param data The `rse_data()` list.
#' @param frame The indicator columns, in the fit's row order.
#' @param groups Each row's group label, or `NULL` for one row per unit.
#' @param calls A list of `parameter_inference()` settings: `vcov_type`,
#'   `level` and `adjust`.
#' @param start Starting coordinates for the reference's Newton steps; by
#'   default the package's estimates, which select the mode and the labels.
#' @param extra Further comparison rows for the fit.
#' @param dataset Dataset name for the report.
#' @param source What the reference is, for the report.
#' @return A list of `harness_case`s: the fit, then one per call.
.robust_se_case <- function(label, fit, model, data, frame, groups, calls,
                            start = NULL, extra = NULL, dataset = "simulated",
                            source = .robust_se_source()) {
  group_log_likelihood <- rse_group_log_likelihood(model, data)
  start <- start %||% rse_theta(model, parameter_inference(fit))
  maximum <- rse_maximise(group_log_likelihood, start)
  covariances <- rse_covariances(group_log_likelihood, maximum$theta)
  fit_case <- compare_fit(
    label, fit,
    rse_fit_reference(model, maximum$theta, data, frame, covariances$observed, groups),
    extra = rbind(extra, compare_values(
      sprintf("%s: reference score at its maximum, in standard errors", label),
      0, maximum$displacement, tolerance = 1e-6)))
  verb_cases <- lapply(calls, function(call) {
    obtained <- parameter_inference(fit, vcov_type = call$vcov_type, level = call$level,
                                    adjust = call$adjust)
    verb_label <- sprintf("%s: parameter_inference(vcov_type = \"%s\", level = %s, adjust = \"%s\")",
                          label, call$vcov_type, call$level, call$adjust)
    compare_verb(
      verb_label,
      list(parameter_inference = rse_keyed_inference(obtained)),
      list(parameter_inference = rse_expected_inference(
        model, maximum$theta, covariances[[call$vcov_type]], obtained, call$level,
        call$adjust)),
      tolerance = .robust_se_tolerances,
      extra = .robust_se_attributes(verb_label, obtained, call, covariances$scaling))
  })
  lapply(c(list(fit_case), verb_cases), function(case) {
    case$comparisons <- transform(case$comparisons, source = source, dataset = dataset,
                                  precision = "machine")
    case
  })
}

#' The reference's description for the report
#' @param origin Where the reference estimates came from, if not the
#'   reference likelihood's own Newton steps alone.
#' @return A single string.
.robust_se_source <- function(origin = NULL) {
  paste(c(origin, sprintf("base-R likelihood, numDeriv %s scores and information",
                          utils::packageVersion("numDeriv"))), collapse = "; ")
}

#' The table's attributes against the arguments and the reference
#'
#' `level`, `adjust` and `vcov_type` must record the call. `scaling_correction`
#' is the MLR factor tr(A^-1 B) / q under `"robust"` and `"opg"`, a ratio of
#' numDeriv quantities compared at 1e-4 as a standard error is, and `NA` under
#' `"observed"` (`?parameter_inference`). A covariate fit's table carries no
#' such attribute (its method documents none), so no row is made for it.
#' @param label Label for the rows.
#' @param obtained The verb's result.
#' @param call The settings it was called with.
#' @param scaling The reference's tr(A^-1 B) / q.
#' @return A comparison `data.frame`.
.robust_se_attributes <- function(label, obtained, call, scaling) {
  recorded <- function(name) {
    as.numeric(identical(attr(obtained, name), call[[name]]))
  }
  rows <- rbind(
    compare_values(sprintf("%s: attribute level", label), call$level,
                   as.numeric(attr(obtained, "level") %||% NA_real_),
                   tolerance = .equivalence_tolerances[["machine"]]),
    compare_values(sprintf("%s: attribute %s records the call", label,
                           c("adjust", "vcov_type")),
                   c(1, 1), c(recorded("adjust"), recorded("vcov_type")), tolerance = 0.5))
  correction <- attr(obtained, "scaling_correction")
  if (is.null(correction)) return(rows)
  if (identical(call$vcov_type, "observed")) {
    return(rbind(rows, compare_values(
      sprintf("%s: attribute scaling_correction is NA", label), 1,
      as.numeric(is.na(correction)), tolerance = 0.5)))
  }
  rbind(rows, compare_values(sprintf("%s: attribute scaling_correction, tr(A^-1 B) / q", label),
                             scaling, correction, tolerance = 1e-4, scale = "relative"))
}

#' Calibrate the reference on a one-profile Gaussian with clusters
#'
#' With one profile the maximum is the sample mean and the variance with
#' divisor n, and the information is block diagonal there, so the CR1
#' sandwich of the mean is vcovCL(lm(y ~ 1), type = "HC0", cadjust = TRUE)
#' and that of the variance is the same for lm(e^2 ~ 1), e the residuals.
#' Those two values fix the reference's arithmetic exactly. The reference's
#' numDeriv derivatives of a function quadratic in the mean and smooth in the
#' variance are good to about 1e-10 relative, so 1e-8 (relative) is pinned.
#' The package's errors are compared with the same values directly at 1e-6:
#' its information is a central difference of the analytic score with step
#' 1e-4, whose truncation error is of order 1e-8.
#'
#' @param frame The two-level fixture.
#' @return A list of `harness_case`s.
.robust_se_calibration <- function(frame) {
  label <- "calibration: one profile, 40 clusters"
  y <- frame$y1
  residual <- y - mean(y)
  cluster_se <- function(response) {
    model <- stats::lm(response ~ 1)
    sqrt(diag(sandwich::vcovCL(model, cluster = frame$school, type = "HC0",
                               cadjust = TRUE)))[[1L]]
  }
  closed_form <- c(mean = cluster_se(y), variance = cluster_se(residual^2))
  model <- rse_model("y1", 1L, 1L)
  data <- rse_data(frame, "y1", "school")
  theta <- c(mean(y), mean(residual^2))
  reference_se <- rse_standard_errors(
    model, theta, rse_covariances(rse_group_log_likelihood(model, data), theta)$robust)
  fit <- multilpa(frame, "y1", "school", n_profiles = 1L, n_group_classes = 1L,
                  n_starts = 1L, seed = 1L, tol = 1e-15, max_iter = 20000L)
  package <- parameter_inference(fit, vcov_type = "robust")
  package_se <- package$standard_error[match(c("mean", "variance"), package$parameter)]
  calibration <- rbind(
    compare_values(sprintf("%s: reference robust SE of the %s against sandwich::vcovCL",
                           label, names(closed_form)),
                   unname(closed_form), reference_se[model$labels$parameter %in%
                                                       c("mean", "variance")],
                   tolerance = 1e-8, scale = "relative"),
    compare_values(sprintf("%s: package robust SE of the %s against sandwich::vcovCL",
                           label, names(closed_form)),
                   unname(closed_form), package_se, tolerance = 1e-6, scale = "relative"))
  .robust_se_case(label, fit, model, data, frame["y1"], groups = frame$school,
                  calls = list(list(vcov_type = "robust", level = 0.95, adjust = "none")),
                  start = theta, extra = calibration,
                  source = .robust_se_source(sprintf("sandwich %s vcovCL",
                                                     utils::packageVersion("sandwich"))))
}

#' mclust's faithful data, VVV with two components, one row per unit
#'
#' The reference point starts from mclust's EM continued to its fixed point,
#' as the mclust suite refines it, aligned to the fit's profiles on the means.
#' @return A list of `harness_case`s.
.robust_se_faithful <- function() {
  frame <- as.data.frame(datasets::faithful)[c("eruptions", "waiting")]
  frame$unit <- seq_len(nrow(frame))
  indicators <- c("eruptions", "waiting")
  values <- as.matrix(frame[indicators])
  reference <- mclust::Mclust(values, G = 2L, modelNames = "VVV", verbose = FALSE,
                              control = mclust::emControl(tol = c(1e-15, 1e-15),
                                                          itmax = c(1e5, 1e5)))
  refined <- mclust::em(data = values, modelName = "VVV", parameters = reference$parameters,
                        control = mclust::emControl(tol = 0, itmax = 20000L))
  stopifnot("mclust's refining EM must finish without an error code" =
              identical(as.numeric(attr(refined, "returnCode")), 0))
  fit <- multilpa(frame, vars = indicators, id = "unit", n_profiles = 2L,
                  n_group_classes = 1L, n_starts = 10L, seed = 1L, tol = 1e-15,
                  max_iter = 20000L, covariance_model = "full", variance_model = "varying")
  means <- t(refined$parameters$mean)
  colnames(means) <- indicators
  order <- align_classes(fit$means[, indicators], means)
  covariances <- refined$parameters$variance$sigma[, , order, drop = FALSE]
  dimnames(covariances) <- list(indicators, indicators, NULL)
  model <- rse_model(indicators, 2L, 1L, covariance = "full")
  start <- rse_theta_from_parameters(model, means[order, , drop = FALSE], covariances,
                                     matrix(refined$parameters$pro[order], 1L))
  .robust_se_case("faithful, VVV, 2 components, one row per unit", fit, model,
                  rse_data(frame, indicators, "unit"), frame[indicators], groups = NULL,
                  calls = list(list(vcov_type = "observed", level = 0.95, adjust = "none"),
                               list(vcov_type = "robust", level = 0.9, adjust = "holm"),
                               list(vcov_type = "opg", level = 0.99, adjust = "BH")),
                  start = start, dataset = "faithful",
                  source = .robust_se_source(sprintf("mclust %s estimates",
                                                     utils::packageVersion("mclust"))))
}
