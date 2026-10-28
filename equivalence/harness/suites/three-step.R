# Three-step distal outcomes and covariates against sandwich, tidySEM, depmixS4
# and likelihoods written from their definitions.
#
# three_step() and r3step() are not fits: they read a fit's posteriors, assign
# every unit to its modal class, estimate the classification error matrix D
# and relate an outcome or covariates to the classes with D held fixed
# (Vermunt 2010; Bakk, Tekle & Vermunt 2013; Di Mari, Oberski & Vermunt 2016).
# The step-one posteriors are therefore the input of every reference here,
# read from `get_results(fit, "posteriors")` and `"group_posteriors"`; their
# own equivalence is the business of the fit suites (mclust, glca,
# mplus-twolevel-lca, depmixs4). From the posteriors on, nothing is taken from
# the package: the modal classes, D, the weights and every step-three estimate
# and standard error are recomputed in `three-step-reference.R`. Every value
# column of every table the two verbs return is compared, through
# `compare_verb()` and its ledger.
#
# References and the correspondence of parameterisations:
#
#   method = "proportional" | "modal"
#     lm(z ~ 0 + class, weights = w) on the outcome stacked once per class,
#     with the unit's posterior (or modal indicator) as weight: the
#     coefficients are the weighted class means. sandwich::vcovCL(cluster =
#     ~unit cluster, type = "HC0", cadjust = TRUE) is the CR1 sandwich the
#     package documents (G / (G - 1) for G clusters); vcov_type =
#     "independent" makes every unit its own cluster.
#   method = "bch"
#     lm() refuses the negative BCH weights, so the class means and their
#     sandwich are written as an estimating equation (score w_ic (z_i - m_c),
#     bread sum_i w_ic, CR1). Those lines are checked against vcovCL on the
#     proportional and modal weights of every case before they are used. The
#     weights are D^-1[W_i, c] with D from tidySEM's
#     `classification_probs_mostlikely()` (an independent construction of the
#     same matrix; it is checked against the D written here).
#   method = "ml"
#     depmixS4::mix(list(W ~ 1, z ~ 1), family = list(multinomial("identity"),
#     gaussian())): the assigned class W is an indicator of the true class
#     with response probabilities fixed at D, and the class sizes are fixed
#     at the average posteriors, as the package documents. variance = "equal"
#     is fit(equal = ). depmixS4's maximum is carried to the fixed point by
#     Newton steps on the same likelihood written in base R with its score and
#     Hessian in closed form (checked against numDeriv), after the two
#     likelihoods are shown to agree at depmixS4's estimates. The standard
#     errors are the CR1 sandwich of that likelihood in (means, log sd); the
#     means' block does not depend on how the sd is charted.
#   r3step(covariates = )
#     depmixS4::mix(list(W ~ 1), prior = ~ covariates) with W's response
#     probabilities fixed at D. depmixS4's prior logits are against the first
#     class; the package's against the last: b_pkg[, k] = b_dep[, k] -
#     b_dep[, K]. "observed" is the inverse closed-form Hessian, "robust" the
#     sandwich of the cluster-summed closed-form unit scores times G / (G -
#     1), both in a chart whose covariates are centred and then mapped back
#     (intercept - sum_j b_j c_j). One case has two covariates (the package
#     fits in an orthonormal basis of them and maps back), adjust = "holm"
#     and ci_level = 0.9. depmixS4::standardError() is not used:
#     the review of 2026-10-27 (probe p11) found its fdHess Hessian about 11%
#     off on this model; that figure is the review's, not re-measured here.
#   r3step(by_group_class = TRUE)
#     No program fits it, so the likelihood sum_g log sum_h eta_h prod_i
#     sum_k pi_hk(x_i) D[k, W_i] is written here and maximised from the
#     step-one solution (group-class profile mixes as intercepts, group-class
#     shares), which keeps the group classes in step one's order.
#   r3step(initial_covariates =, transition_covariates =)
#     Step two classifies each occasion from its own indicators with the lta
#     fit's (diagonal Gaussian) measurement model and the profile shares
#     pooled over occasions. Step three is depmixS4::depmix(list(W ~ 1), prior
#     = ~ initial, transition = ~ transition) with W fixed at D. depmixS4
#     reads a transition covariate at the origin row, the package at the
#     destination occasion, so depmixS4 is given the next row's value. Its
#     logits (against the first profile of each row) become the package's
#     (initial against the last profile, transitions against staying) by
#     differencing, and are carried to the fixed point with the recursion's
#     score from a forward-backward pass (Fisher's identity). Two and three
#     profiles; the r3step result's `classification_errors` attribute (step
#     two's D) is compared with the reference's D.
#
# Contrasts (contrast = "pairs") are m_b - m_a for a < b, with standard error
# sqrt(V_aa + V_bb - 2 V_ab) from the same covariance, Wald p-values and
# stats::p.adjust(method = "BH") across the pairs (one srl case uses
# "bonferroni" and ci_level = 0.9); r3step's p_adjusted covers the covariate
# terms only, as documented.
#
# Data: srl (300 rows, three profiles, one unit per row, distal test_anxiety,
# covariate intrinsic_value, which is not an indicator); course_engagement
# (1,422 rows of 106 students, two profiles x two group classes, distal and
# covariate previous_grade, and its student mean at the group level, where an
# outcome must be constant within a group); study_hours (150 students x 6
# weeks, lta on hours and score with two and three profiles, initial
# covariate motivation, transition covariate sleep). Step-one fits use
# n_starts = 3 and seed = 1.

source(file.path("equivalence", "harness", "three-step-reference.R"), local = environment())

suite_three_step <- function() {
  require_suite_packages(
    "sandwich", "tidySEM", "depmixS4", "Rsolnp", "numDeriv",
    reason = paste("the cluster-robust sandwich of the weighted means, an independent",
                   "classification error matrix, the ML step-three, R3STEP and latent",
                   "Markov maxima with D fixed, and the reference numerical derivatives"))
  stack_cases(c(.three_step_srl(), .three_step_engagement(), .three_step_chain()))
}

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

#' Tolerances of the closed-form methods
#'
#' bch, proportional and modal are weighted means and a sandwich: the package
#' and the reference do the same arithmetic in another order (lm's QR, a
#' cluster sum), so they agree to rounding. 1e-10 leaves room for the QR's
#' cancellation and nothing else.
.three_step_closed <- c(default = 1e-10)

#' Tolerances of the standard errors and what is formed from them
#'
#' Standard errors, intervals, Wald statistics and p-values. The reference's
#' derivatives are exact (closed form for ML and R3STEP, numDeriv's Jacobian
#' of a closed-form score for the chain, numDeriv's Hessian for the
#' group-class model); the package's Hessians are central differences of its
#' analytic score (step 1e-5, error about 1e-10 relative). What is formed
#' from an estimate inherits its convergence gap, divided by a standard error,
#' so each family takes its estimates' bound or more: 1e-5 by default.
#' @param table The verb's table name.
#' @param tolerance The bound.
#' @return Named tolerances.
.three_step_error_tolerances <- function(table, tolerance = 1e-5) {
  columns <- c("standard_error", "conf_low", "conf_high", "statistic", "p_value",
               "p_adjusted")
  stats::setNames(rep(tolerance, length(columns)), paste0(table, "$", columns))
}

#' Tolerances of the ML step-three estimates
#'
#' Means, variances and expected class sizes: 1e-8, both sides at the fixed
#' point. The package polishes its EM by Newton steps until its scaled score
#' is below 1e-9 (1 + |log-likelihood|), which admits a gap of order 1e-8 to
#' 1e-7 standard errors; the reference carries depmixS4's maximum on by Newton
#' steps on the closed-form score, accepted on the score's norm, to a score
#' below 1e-13. Measured 2026-10-27: the two agree to 4.3e-13 (the package's
#' last Newton step lands quadratically, well inside its rule).
#'
#' Standard errors and what is formed from them: 1e-6. A gap of 1e-8 in a
#' mean whose standard error is 0.05 moves its Wald statistic by 2e-7, and
#' the package's central-difference Hessian adds about 1e-10. Measured
#' 2026-10-27: 3.7e-11.
.three_step_ml <- c(default = 1e-8, .three_step_error_tolerances("three_step", 1e-6))

#' Tolerances of the pooled and group-class R3STEP coefficients
#'
#' Coefficients: 1e-5, a convergence difference in standard-error units. The
#' package stops its pooled and group-class R3STEP searches with optim(BFGS,
#' reltol = 1e-12) and does not polish them, so on these likelihoods (-37 to
#' -916) it stops once a step gains less than about 4e-11 to 9e-10, and a
#' deficit delta leaves an estimate up to sqrt(2 delta), about 4e-5, standard
#' errors away. Measured 2026-10-27 with the closed-form score: the
#' reference's score is below 1e-13; the largest gaps are 3.4e-8 (srl),
#' 1.1e-7 (course_engagement, groups) and 1.2e-6 (by_group_class). The same
#' maximum is checked directly by the likelihood rows (1e-8). Standard errors
#' and what is formed from them inherit the gap: 1e-5.
.three_step_r3 <- c(default = 1e-8, `r3step$estimate` = 1e-5,
                    .three_step_error_tolerances("r3step"))

#' Tolerances of the three-step latent Markov coefficients
#'
#' Coefficients, and the standard errors and what is formed from them: 1e-6.
#' The package's chain search is optim(BFGS, reltol = 1e-14), unpolished.
#' Measured 2026-10-27 with the reference's forward-backward score: the
#' reference's score is 5e-15; at the package's estimates it is 2.0e-7 (two
#' profiles) and 8.3e-7 (three), the Newton step from there predicts the
#' observed gap to 1e-15, and the gaps are 1.0e-8 and 1.9e-7 (1.9e-7
#' standard errors; the standard errors themselves agree to 9e-9). 1e-6
#' sits five times above the larger gap. The 1e-8 the review of 2026-10-27
#' proposed was measured against a reference that had itself stalled 2e-8
#' short, and is below what the package reaches. D, step two's classification
#' error matrix, is a closed form of the occasion posteriors in both: 1e-12.
.three_step_chain_tolerance <- c(default = 1e-8, `r3step$estimate` = 1e-6,
                                 .three_step_error_tolerances("r3step", 1e-6))

#' Stationarity a reference must reach before it is compared
#'
#' Largest absolute score at the reference's estimate. A closed-form score is
#' exact to rounding (measured 4e-15 to 7e-14 on these cases), so 1e-8 is a
#' floor well clear of noise that still refuses a stalled refinement (the
#' review of 2026-10-27 found one at 9.85e-6). numDeriv's gradient of a
#' likelihood of size 1,000 carries errors near 1e-7, so the one reference
#' differentiated numerically (group-class R3STEP, measured 2.3e-8) is held
#' to 1e-6. A closed-form score must match numDeriv's gradient to 1e-4
#' relative (its Richardson error, measured up to 9.5e-6; a wrong formula
#' shows differences of order one) and a closed-form Hessian numDeriv's
#' Jacobian of that score to 1e-8 (measured 1.5e-11 to 1.9e-10).
.three_step_stationary <- c(`closed form` = 1e-8, numerical = 1e-6)

# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

#' srl: three profiles, one unit per row
#' @return A list of `harness_case`s.
.three_step_srl <- function() {
  data <- srl
  data$unit <- seq_len(nrow(data))
  fit <- multilpa(data, vars = c("cognitive_strategies", "self_efficacy", "self_regulation"),
                  id = "unit", n_profiles = 3, n_group_classes = 1, n_starts = 3, seed = 1)
  posteriors <- .three_step_posteriors(fit)
  z <- data$test_anxiety
  label <- "srl, 3 profiles"
  distal <- function(method, contrast = "none", variance = NULL, adjust = "BH",
                     ci_level = 0.95) {
    arguments <- list(fit, data, "test_anxiety", method = method, contrast = contrast)
    if (!is.null(variance)) arguments$variance <- variance
    default_level <- isTRUE(all.equal(ci_level, 0.95))
    if (!identical(adjust, "BH")) arguments$adjust <- adjust
    if (!default_level) arguments$ci_level <- ci_level
    obtained <- do.call(three_step, arguments)
    .three_step_distal_case(sprintf("%s: three_step(method = \"%s\"%s%s%s%s)", label, method,
                                    if (is.null(variance)) "" else
                                      sprintf(", variance = \"%s\"", variance),
                                    if (identical(contrast, "pairs")) ", contrast = \"pairs\"" else "",
                                    if (identical(adjust, "BH")) "" else
                                      sprintf(", adjust = \"%s\"", adjust),
                                    if (default_level) "" else
                                      sprintf(", ci_level = %s", ci_level)),
                            obtained, posteriors, z, data$unit, "individuals", method,
                            contrast, variance %||% "varying", "srl", adjust, ci_level)
  }
  # Two covariates exercise the package's orthonormal basis and its
  # back-transform, and the term order; holm and 90% intervals the
  # correction and the interval arguments.
  covariates <- c("intrinsic_value", "test_anxiety")
  c(lapply(c("bch", "proportional", "modal", "ml"), distal),
    list(distal("ml", variance = "equal"), distal("bch", contrast = "pairs"),
         distal("ml", contrast = "pairs"),
         distal("bch", contrast = "pairs", adjust = "bonferroni", ci_level = 0.9),
         distal("ml", ci_level = 0.9)),
    lapply(c("observed", "robust"), function(vcov_type) {
      .three_step_r3_case(sprintf("%s: r3step(vcov_type = \"%s\")", label, vcov_type),
                          r3step(fit, data, "intrinsic_value", vcov_type = vcov_type),
                          as.matrix(data["intrinsic_value"]), posteriors, data$unit,
                          "individuals", vcov_type, "srl")
    }),
    list(.three_step_r3_case(
      sprintf(paste("%s: r3step(c(\"intrinsic_value\", \"test_anxiety\"), vcov_type =",
                    "\"robust\", adjust = \"holm\", ci_level = 0.9)"), label),
      r3step(fit, data, covariates, vcov_type = "robust", adjust = "holm", ci_level = 0.9),
      as.matrix(data[covariates]), posteriors, data$unit, "individuals", "robust", "srl",
      adjust = "holm", ci_level = 0.9)))
}

#' course_engagement: two profiles x two group classes, rows nested in students
#' @return A list of `harness_case`s.
.three_step_engagement <- function() {
  data <- course_engagement
  # A group-level outcome and covariate must be constant within a group.
  data$student_grade <- stats::ave(data$previous_grade, data$student)
  fit <- multilpa(data, vars = c("browse", "lectures", "forum_read", "forum_post",
                                 "attendance"),
                  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1)
  posteriors <- .three_step_posteriors(fit)
  groups <- unique(data$student)
  group_posteriors <- .three_step_group_posteriors(fit, groups)
  group_grade <- data$student_grade[!duplicated(data$student)]
  label <- "course_engagement, 2 profiles x 2 group classes"
  individuals <- function(method, vcov_type = "cluster", contrast = "none",
                          variance = NULL) {
    arguments <- list(fit, data, "previous_grade", method = method, vcov_type = vcov_type,
                      contrast = contrast)
    if (!is.null(variance)) arguments$variance <- variance
    cluster <- if (identical(vcov_type, "cluster")) data$student else seq_len(nrow(data))
    .three_step_distal_case(
      sprintf("%s: three_step(method = \"%s\", vcov_type = \"%s\"%s%s)", label, method,
              vcov_type, if (is.null(variance)) "" else sprintf(", variance = \"%s\"", variance),
              if (identical(contrast, "pairs")) ", contrast = \"pairs\"" else ""),
      do.call(three_step, arguments), posteriors, data$previous_grade, cluster,
      "individuals", method, contrast, variance %||% "varying", "course_engagement")
  }
  at_groups <- function(method) {
    .three_step_distal_case(
      sprintf("%s: three_step(level = \"groups\", method = \"%s\")", label, method),
      three_step(fit, data, "student_grade", level = "groups", method = method),
      group_posteriors, group_grade, seq_along(groups), "groups", method, "none",
      "varying", "course_engagement")
  }
  pooled <- function(level, vcov_type) {
    covariate <- if (identical(level, "groups")) "student_grade" else "previous_grade"
    .three_step_r3_case(
      sprintf("%s: r3step(level = \"%s\", vcov_type = \"%s\")", label, level, vcov_type),
      r3step(fit, data, covariate, level = level, vcov_type = vcov_type),
      if (identical(level, "groups")) cbind(student_grade = group_grade) else
        as.matrix(data["previous_grade"]),
      if (identical(level, "groups")) group_posteriors else posteriors,
      if (identical(level, "groups")) seq_along(groups) else data$student,
      level, vcov_type, "course_engagement")
  }
  c(lapply(c("bch", "proportional", "modal", "ml"), individuals),
    lapply(c("bch", "proportional", "modal", "ml"), individuals, vcov_type = "independent"),
    list(individuals("ml", vcov_type = "independent", variance = "equal"),
         individuals("proportional", contrast = "pairs")),
    lapply(c("bch", "proportional", "modal", "ml"), at_groups),
    lapply(c("observed", "robust"), pooled, level = "individuals"),
    lapply(c("observed", "robust"), pooled, level = "groups"),
    lapply(c("observed", "robust"), function(vcov_type) {
      .three_step_group_case(
        sprintf("%s: r3step(by_group_class = TRUE, vcov_type = \"%s\")", label, vcov_type),
        r3step(fit, data, "previous_grade", by_group_class = TRUE, vcov_type = vcov_type),
        fit, as.matrix(data["previous_grade"]), posteriors, data$student, vcov_type)
    }))
}

#' study_hours: the three-step latent Markov model, two and three profiles
#' @return A list of `harness_case`s.
.three_step_chain <- function() {
  data <- study_hours[order(study_hours$student, study_hours$week), ]
  vars <- c("hours", "score")
  unlist(lapply(2:3, function(n_profiles) {
    fit <- lta(data, vars = vars, id = "student", n_profiles = n_profiles, time = "week",
               n_starts = 3, seed = 1)
    profiles <- get_results(fit, "profiles")
    covariances <- get_results(fit, "covariances")
    stopifnot("the reference classifies occasions under a diagonal Gaussian model" =
                all(covariances$covariance[covariances$indicator !=
                                             covariances$indicator_2] == 0))
    means <- unclass(stats::xtabs(mean ~ profile + indicator, profiles))[, vars]
    variances <- unclass(stats::xtabs(variance ~ profile + indicator, profiles))[, vars]
    shares <- colMeans(.three_step_posteriors(fit))
    reference <- three_step_chain_reference(data, "student", vars, means, variances, shares,
                                            "motivation", "sleep")
    .three_step_check_depmix(reference, sprintf("study_hours chain, %d profiles", n_profiles))
    source_label <- sprintf("depmixS4 %s; numDeriv %s", utils::packageVersion("depmixS4"),
                            utils::packageVersion("numDeriv"))
    lapply(c("observed", "robust"), function(vcov_type) {
      obtained <- r3step(fit, data, initial_covariates = "motivation",
                         transition_covariates = "sleep", vcov_type = vcov_type)
      label <- sprintf(paste("study_hours, %d profiles: r3step(initial_covariates =",
                             "\"motivation\", transition_covariates = \"sleep\",",
                             "vcov_type = \"%s\")"), n_profiles, vcov_type)
      expected <- three_step_coefficient_table(reference$labels, reference$estimates,
                                               reference[[vcov_type]], reference$tested)
      case <- compare_verb(
        label, list(r3step = keyed_table(.three_step_chain_keys(obtained),
                                         c("block", "from", "to", "term"))),
        list(r3step = .three_step_chain_keys(expected)),
        tolerance = .three_step_chain_tolerance,
        extra = rbind(
          compare_values(sprintf("%s: maximised log-likelihood", label),
                         reference$log_likelihood, attr(obtained, "log_likelihood"),
                         tolerance = 1e-8, scale = "relative"),
          # The initial rows have no origin; the key stands in for it only
          # where the package reports none.
          compare_values(sprintf("%s: `from` is missing exactly on the initial rows", label),
                         1, as.numeric(identical(is.na(obtained$from),
                                                 obtained$block == "initial")),
                         tolerance = 0.5),
          compare_block(sprintf("%s: attribute classification_errors (D, rows true)", label),
                        as.vector(reference$errors),
                        as.vector(unclass(attr(obtained, "classification_errors"))),
                        tolerance = 1e-12)))
      .three_step_tag(case, source_label, "study_hours")
    })
  }), recursive = FALSE)
}

#' Key the chain's rows by block, origin, destination and term
#'
#' The initial rows have no origin (`from` is missing), and a join drops a
#' missing key, so the origin is written `"(initial)"` there on both sides.
#' @param table An r3step() transition table or its expected counterpart.
#' @return The table with `from` filled on the initial rows.
.three_step_chain_keys <- function(table) {
  table$from[table$block == "initial" & is.na(table$from)] <- "(initial)"
  table
}

# ---------------------------------------------------------------------------
# One verb call against its reference
# ---------------------------------------------------------------------------

#' Compare one three_step() table
#' @param label Label of the call.
#' @param obtained The table three_step() returned.
#' @param posteriors Units-by-classes posteriors at the call's level.
#' @param z The outcome, one value per unit.
#' @param cluster The independent unit of each unit.
#' @param level,method,contrast,variance,adjust,ci_level As in the call.
#' @param dataset Dataset name for the report.
#' @return A `harness_case`.
.three_step_distal_case <- function(label, obtained, posteriors, z, cluster, level, method,
                                    contrast, variance, dataset, adjust = "BH",
                                    ci_level = 0.95) {
  keys <- c("level", "method", "class", if (identical(contrast, "pairs")) "reference_class")
  modal <- three_step_modal(posteriors)
  errors <- three_step_error_matrix(posteriors, modal)
  .three_step_check_sandwich(z, posteriors, modal, cluster, label)
  if (identical(method, "ml")) {
    reference <- three_step_ml_reference(z, posteriors, cluster, variance)
    .three_step_check_depmix(reference, label)
    estimates <- reference$means
    covariance <- reference$covariance
    tolerance <- .three_step_ml
    source_label <- sprintf("depmixS4 %s; numDeriv %s", utils::packageVersion("depmixS4"),
                            utils::packageVersion("numDeriv"))
  } else {
    weights <- switch(method,
      proportional = posteriors,
      modal = vapply(seq_len(ncol(posteriors)), function(class) as.numeric(modal == class),
                     numeric(length(modal))),
      bch = {
        # tidySEM labels its result from the posteriors' column names.
        named <- posteriors
        colnames(named) <- paste0("class_", seq_len(ncol(posteriors)))
        independent <- unname(as.matrix(tidySEM:::classification_probs_mostlikely(named)))
        stopifnot("tidySEM's error matrix must equal the one written here" =
                    max(abs(independent - unname(errors))) < 1e-12)
        solve(independent)[modal, , drop = FALSE]
      })
    reference <- if (identical(method, "bch")) {
      three_step_estimating_sandwich(z, weights, cluster)
    } else three_step_lm_sandwich(z, weights, cluster)
    estimates <- reference$estimates
    covariance <- reference$covariance
    tolerance <- .three_step_closed
    source_label <- if (identical(method, "bch")) {
      sprintf("tidySEM %s; estimating-equation sandwich", utils::packageVersion("tidySEM"))
    } else sprintf("sandwich %s (vcovCL)", utils::packageVersion("sandwich"))
  }
  expected <- if (identical(contrast, "pairs")) {
    three_step_pairs_table(level, method, estimates, covariance, adjust = adjust,
                           ci_level = ci_level)
  } else {
    three_step_mean_table(level, method, estimates, covariance,
                          effective_n = if (identical(method, "ml")) reference$effective_n else
                            colSums(weights)^2 / colSums(weights^2),
                          variance = if (identical(method, "ml")) reference$variances,
                          ci_level = ci_level)
  }
  case <- compare_verb(label, list(three_step = keyed_table(obtained, keys)),
                       list(three_step = expected), tolerance = tolerance)
  .three_step_tag(case, source_label, dataset)
}

#' Compare one pooled r3step() table
#' @param label Label of the call.
#' @param obtained The table r3step() returned.
#' @param covariates Units-by-covariates matrix at the call's level.
#' @param posteriors,cluster As for [.three_step_distal_case()].
#' @param level,vcov_type,adjust,ci_level As in the call.
#' @param dataset Dataset name for the report.
#' @return A `harness_case`.
.three_step_r3_case <- function(label, obtained, covariates, posteriors, cluster, level,
                                vcov_type, dataset, adjust = "BH", ci_level = 0.95) {
  reference <- three_step_r3_reference(covariates, posteriors, cluster)
  .three_step_check_depmix(reference, label)
  n_free <- ncol(posteriors) - 1L
  prefix <- if (identical(level, "groups")) "group_class_" else "profile_"
  labels <- data.frame(level = level,
                       outcome = rep(paste0(prefix, seq_len(n_free)), each = length(reference$terms)),
                       term = rep(reference$terms, n_free), stringsAsFactors = FALSE)
  expected <- three_step_coefficient_table(labels, reference$estimates, reference[[vcov_type]],
                                           tested = labels$term != "(Intercept)",
                                           adjust = adjust, ci_level = ci_level)
  case <- compare_verb(label, list(r3step = keyed_table(obtained, c("level", "outcome", "term"))),
                       list(r3step = expected), tolerance = .three_step_r3)
  .three_step_tag(case, sprintf("depmixS4 %s; numDeriv %s", utils::packageVersion("depmixS4"),
                                utils::packageVersion("numDeriv")), dataset)
}

#' Compare one r3step(by_group_class = TRUE) table
#' @param label Label of the call.
#' @param obtained The table r3step() returned.
#' @param fit The step-one fit, for its group-class profile mixes and shares.
#' @param covariates Rows-by-covariates matrix.
#' @param posteriors Rows-by-profiles marginal posteriors.
#' @param group Group label of each row.
#' @param vcov_type As in the call.
#' @return A `harness_case`.
.three_step_group_case <- function(label, obtained, fit, covariates, posteriors, group,
                                   vcov_type) {
  mixes <- get_results(fit, "profile_probabilities")
  profile_mix <- unclass(stats::xtabs(probability ~ group_class + profile, mixes))
  group_shares <- mixes$group_class_probability[!duplicated(mixes$group_class)][
    order(unique(mixes$group_class))]
  reference <- three_step_group_reference(covariates, posteriors, group, profile_mix,
                                          group_shares)
  stopifnot("the reference must reach a stationary maximum" =
              max(abs(reference$newton$gradient)) <
                .three_step_stationary[[reference$newton$kind]] &&
              reference$newton$concave)
  n_free <- ncol(posteriors) - 1L
  n_group_classes <- nrow(profile_mix)
  labels <- rbind(
    data.frame(level = "individuals",
               outcome = rep(paste0("profile_", seq_len(n_free)), each = n_group_classes),
               term = rep(paste0("group_class_", seq_len(n_group_classes)), n_free)),
    data.frame(level = "individuals",
               outcome = rep(paste0("profile_", seq_len(n_free)), each = ncol(covariates)),
               term = rep(colnames(covariates), n_free)),
    data.frame(level = rep("groups", n_group_classes - 1L),
               outcome = paste0("group_class_", seq_len(n_group_classes - 1L)),
               term = rep("(Intercept)", n_group_classes - 1L)))
  expected <- three_step_coefficient_table(
    labels, reference$estimates, reference[[vcov_type]],
    tested = labels$level == "individuals" & labels$term %in% colnames(covariates))
  case <- compare_verb(
    label, list(r3step = keyed_table(obtained, c("level", "outcome", "term"))),
    list(r3step = expected), tolerance = .three_step_r3,
    extra = compare_values(sprintf("%s: maximised log-likelihood", label),
                           reference$log_likelihood, attr(obtained, "log_likelihood"),
                           tolerance = 1e-8, scale = "relative"))
  .three_step_tag(case, sprintf("likelihood written here; numDeriv %s",
                                utils::packageVersion("numDeriv")), "course_engagement")
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' A fit's profile posteriors as a rows-by-profiles matrix in data order
#' @param fit A fitted model.
#' @return A numeric matrix.
.three_step_posteriors <- function(fit) {
  wide <- get_results(fit, "posteriors", format = "wide")
  stopifnot("the posteriors must come in the fit's row order" =
              identical(as.integer(wide$row), seq_len(nrow(wide))))
  unname(as.matrix(wide[grep("^posterior_profile_", names(wide))]))
}

#' A fit's group posteriors as a groups-by-group-classes matrix
#' @param fit A two-level fit.
#' @param groups Group labels in the order wanted (first appearance).
#' @return A numeric matrix.
.three_step_group_posteriors <- function(fit, groups) {
  long <- get_results(fit, "group_posteriors")
  wide <- unclass(stats::xtabs(posterior ~ group + group_class, long))
  unname(wide[as.character(groups), , drop = FALSE])
}

#' Check the estimating-equation sandwich against sandwich::vcovCL
#'
#' The BCH reference is the estimating equation alone; it is trusted only
#' because, on the proportional and modal weights of the same case, it equals
#' vcovCL (and lm's coefficients) to rounding.
#' @return `NULL`, invisibly, or an error naming the case.
.three_step_check_sandwich <- function(z, posteriors, modal, cluster, label) {
  indicator <- vapply(seq_len(ncol(posteriors)), function(class) as.numeric(modal == class),
                      numeric(length(modal)))
  agree <- vapply(list(posteriors, indicator), function(weights) {
    by_hand <- three_step_estimating_sandwich(z, weights, cluster)
    by_lm <- three_step_lm_sandwich(z, weights, cluster)
    max(abs(by_hand$estimates - by_lm$estimates)) < 1e-10 &&
      max(abs(by_hand$covariance - by_lm$covariance)) < 1e-12
  }, logical(1))
  if (!all(agree)) {
    stop(sprintf("%s: the estimating-equation sandwich does not reproduce vcovCL", label),
         call. = FALSE)
  }
  invisible(NULL)
}

#' Check that a depmixS4 maximum and the likelihood written here agree
#'
#' At depmixS4's estimates the two likelihoods must be the same function
#' (1e-9 relative); the Newton refinement must not lose likelihood and must
#' stay in depmixS4's basin (gain below 1e-6); it must end at a maximum (a
#' negative definite Hessian) whose score is within
#' [.three_step_stationary]; and a closed-form score and Hessian must match
#' numDeriv's (1e-4 and 1e-8 relative, see [.three_step_stationary]).
#' @return `NULL`, invisibly, or an error naming the case.
.three_step_check_depmix <- function(reference, label) {
  same_function <- abs(reference$own_at_depmix - reference$depmix_log_likelihood) <=
    1e-9 * abs(reference$depmix_log_likelihood)
  gain <- reference$log_likelihood - reference$depmix_log_likelihood
  score <- max(abs(reference$newton$gradient))
  stationary <- score < .three_step_stationary[[reference$newton$kind]] &&
    isTRUE(reference$newton$concave)
  derivatives <- reference$derivative_check %||% c(score = 0, hessian = 0)
  exact <- derivatives[["score"]] < 1e-4 && derivatives[["hessian"]] < 1e-8
  if (!same_function || gain < -1e-9 || gain > 1e-6 || !stationary || !exact) {
    stop(sprintf(paste("%s: reference check failed (likelihoods %.12g vs %.12g, Newton",
                       "gain %.3g, largest score %.3g, concave %s, derivative check",
                       "%.3g / %.3g)"), label, reference$own_at_depmix,
                 reference$depmix_log_likelihood, gain, score, reference$newton$concave,
                 derivatives[["score"]], derivatives[["hessian"]]),
         call. = FALSE)
  }
  invisible(NULL)
}

#' Label a case's rows with their source and dataset
#' @param case A `harness_case`.
#' @param source,dataset Labels.
#' @return The case.
.three_step_tag <- function(case, source, dataset) {
  case$comparisons <- transform(case$comparisons, source = source, dataset = dataset,
                                precision = "machine")
  case
}
