# Two-level latent class analysis against multilevLCA (Lyrvall, Di Mari, Bakk,
# Oser & Kuha 2025, Multivariate Behavioral Research 60(4):731-747).
#
# multilevLCA fits the same nonparametric two-level mixture and is the only
# other CRAN package to publish peer-reviewed worked numbers on bundled public
# data. Its parameter blocks map onto multilpa's as
#
#   vOmega <-> group_probabilities      P(group class)
#   mPi    <-> profile_probabilities    P(profile | group class), transposed
#   mPhi   <-> response_probabilities   P(item = 1 | profile), binary items
#
# and its two information criteria correspond to multilpa's two sample-size
# conventions: BIClow to bic_individual, BIChigh to bic_groups.
#
# `fixedpars = 0` is required. multilevLCA defaults to `fixedpars = 1`, a
# two-step estimator whose measurement model is held fixed from step one, which
# is deliberately not the one-step maximum likelihood multilpa computes.
#
# Every output of every fit is compared through `compare_fit()`. The reference
# is multilevLCA's solution with its extended output (`extout = TRUE`): its
# estimates and likelihood, and its posteriors evaluated at those estimates by
# the base-R two-level likelihood in two-level-lca-reference.R, after that
# evaluation has been shown to reproduce multilevLCA's own likelihood, its
# marginal profile posteriors (`mPMsumX`), its per-group log-likelihoods
# (`vLLK`) and its group posteriors (`mPW`). Two of multilevLCA's outputs are
# not used as values, for reasons established on dataTOY:
#
#   - `mPW` is floored at 1e-5 and renormalised for display in most rows (it
#     reports 9.9999e-6 where the posterior is 1e-57), while its EM uses the unfloored
#     values (`vOmega` is exactly the mean of the unfloored posteriors, 2/3 on
#     dataTOY 3 x 2, against 0.6666633 for the floored ones). It is checked
#     through the same display step. For the same reason `ICL_BIChigh` and
#     `R2entr_high`, computed from the floored values, are not compared
#     (ICL_BIChigh is 0.0075 above the unfloored value on dataTOY 2 x 2), and
#     `R2entr_low` is an
#     entropy R-squared against the class-proportion entropy (0.822 where the
#     package's relative entropy against log K is 0.826), a different quantity.
#   - Its standard errors (`SEs`, `Varmat`) are `solve(crossprod(mScore))`, an
#     outer product of one score row per individual, which treats the rows of
#     a group as independent; and the group-class column of `mScore` does not
#     sum to zero at its maximum (666.7 on dataTOY 3 x 2), so it is not the
#     gradient of the likelihood. The package has no such estimator, so
#     nothing is compared with them. Instead, as in the glca suite,
#     `responses$probability_standard_error` and the class-probability
#     standard errors of `parameter_inference()` are compared with the inverse
#     of a numDeriv Hessian of the base-R likelihood at multilevLCA's
#     estimates, and `vcov_type = "opg"` with the per-group BHHH estimate where
#     the package documents it as available (more groups than parameters).

source(file.path("equivalence", "harness", "two-level-lca-reference.R"), local = environment())

suite_multilevlca <- function() {
  require_suite_packages(
    "multilevLCA", "numDeriv",
    reason = "the third independent two-level implementation, the dataTOY and dataIEA data, and the reference observed information")
  stack_cases(list(
    .multilevlca_case(3L, 2L), .multilevlca_case(2L, 2L), .multilevlca_case(3L, 3L, tolerance = c(default = 1e-5)),
    .multilevlca_published()))
}

#' Fit multilevLCA, run on to its own fixed point when asked
#'
#' multilevLCA stops when the relative change in the log-likelihood falls below
#' `tol` (1e-8 by default), which on dataTOY 3 x 2 leaves its reported
#' posteriors 1.6e-6 from the values at its reported estimates. Even a `tol`
#' of 1e-15 is not enough on a flat likelihood: on dataTOY 3 x 3 the
#' log-likelihood stops changing in double precision at iteration 170 (its
#' `eps` is 0) while P(group class 1) is still moving along a ridge,
#' 0.5103328 there, 0.51033365 at 300 iterations and the same to 12 digits at
#' 600 and 1000. Given `iterations`, `tol = -1`, which no change satisfies,
#' makes its own EM run that many steps from the same start (the same seed),
#' as the mclust suite continues mclust's EM.
#'
#' @param data,items,profiles,id,group_classes As for `multilevLCA::multiLCA()`.
#' @param tol Its stopping tolerance, when `iterations` is `NULL`.
#' @param iterations A fixed number of EM steps, or `NULL`.
#' @return The `multiLCA` result with extended output, without the data
#'   columns (which would turn the posterior matrices into character ones for a
#'   character group identifier).
.multilevlca_fit <- function(data, items, profiles, id, group_classes, tol = 1e-8,
                             iterations = NULL) {
  set.seed(2023L)
  multilevLCA::multiLCA(data, Y = items, iT = profiles, id_high = id, iM = group_classes,
                        fixedpars = 0, verbose = FALSE, extout = TRUE, dataout = FALSE,
                        tol = if (is.null(iterations)) tol else -1,
                        maxIter = iterations %||% 100000L)
}

#' Compare one dataTOY configuration against multilevLCA
#' @param profiles Number of lower-level latent classes.
#' @param group_classes Number of higher-level latent classes.
#' @param tolerance Named tolerances for [compare_fit()]. The 3 x 3 model sits on
#'   the flat ridge described at [.multilevlca_fit()]: both programs stop on it
#'   where the log-likelihood no longer changes in double precision, and on
#'   x86 Linux the package stopped 1.6e-6 to 3.5e-6 along it in the group
#'   posteriors and the quantities built on them (CI, 2026-10-28; arm64 macOS
#'   agreed to 1e-6). That case is held to 1e-5.
#' @return A `harness_case`.
.multilevlca_case <- function(profiles, group_classes, tolerance = numeric()) {
  observed <- get(utils::data(list = "dataTOY", package = "multilevLCA",
                              envir = environment()), envir = environment())
  items <- paste0("Y_", seq_len(10))
  reference <- .multilevlca_fit(observed, items, profiles, "id_high", group_classes,
                                iterations = 2000L)
  # The package's tightest stopping rule, as in the mclust and glca suites.
  fit <- multilpa(observed, vars = items, id = "id_high",
                  n_profiles = profiles, n_group_classes = group_classes,
                  categorical = items, n_starts = 30L, seed = 1L,
                  tol = 1e-15, max_iter = 20000L)
  label <- sprintf("dataTOY, %d profiles x %d group classes", profiles, group_classes)
  log_likelihood <- utils::tail(as.numeric(reference$LLKSeries), 1L)
  shortfall <- log_likelihood - as.numeric(logLik(fit))
  if (abs(shortfall) > 1e-4) {
    # Different modes: nothing but the shortfall is comparable, as in the
    # mclust suite.
    rows <- compare_values(sprintf("%s: log-likelihood shortfall against multilevLCA", label),
                           0, max(shortfall, 0), tolerance = 1e-5)
    inventory <- fit_outputs(fit)
    case <- harness_case(rows, output_ledger(label, inventory, rows, exempt = stats::setNames(
      rep("multilevLCA stopped at a different mode, so its estimates describe another solution",
          length(unique(inventory$table))), unique(inventory$table))))
  } else {
    built <- .multilevlca_reference(reference, observed, items, "id_high", fit)
    alignment <- two_level_alignment(built$reference$responses,
                                     built$reference$profile_probabilities, fit)
    opg <- if (fit$n_groups > fit$n_parameters) {
      two_level_error_rows(label, fit, built$errors$opg, alignment, "opg")
    } else {
      refused <- tryCatch(parameter_inference(fit, vcov_type = "opg"),
                          latents_too_few_groups = function(condition) condition)
      compare_values(sprintf("%s: OPG refused with %d groups for %d parameters", label,
                             fit$n_groups, fit$n_parameters),
                     1, as.numeric(inherits(refused, "latents_too_few_groups")),
                     tolerance = 0.5)
    }
    case <- compare_fit(
      label, fit, built$reference, tolerance = tolerance,
      extra = rbind(
        compare_values(sprintf("%s: log-likelihood shortfall against multilevLCA", label),
                       0, max(shortfall, 0), tolerance = 1e-5),
        # multilevLCA's own criteria, as it reports them.
        compare_values(sprintf("%s: AIC as multilevLCA reports it", label),
                       reference$AIC, fit$aic, tolerance = 1e-4,
                       output = "information_criteria$aic"),
        compare_values(sprintf("%s: BIClow (individual sample size) as multilevLCA reports it",
                               label),
                       reference$BIClow, fit$bic_individual, tolerance = 1e-4,
                       output = "information_criteria$bic_individual"),
        compare_values(sprintf("%s: BIChigh (group sample size) as multilevLCA reports it",
                               label),
                       reference$BIChigh, fit$bic_groups, tolerance = 1e-4,
                       output = "information_criteria$bic_groups"),
        .with_output_columns(two_level_error_rows(label, fit, built$errors$observed,
                                                  alignment, "observed",
                                                  levels = c("profile", "group"))),
        .with_output_columns(opg)))
  }
  case$comparisons <- transform(case$comparisons,
                                source = sprintf("multilevLCA %s",
                                                 utils::packageVersion("multilevLCA")),
                                dataset = "dataTOY", precision = "machine")
  case
}

#' multilevLCA's solution in the harness's reference shape
#'
#' @param reference A `multiLCA` result with extended output.
#' @param analysis The rows multilevLCA analysed (it deletes incomplete rows),
#'   in the order the package fit was given them.
#' @param items The binary indicators.
#' @param id The group column.
#' @param fit The package fit, for its category labels and probability bound
#'   only; no estimate is read from it.
#' @param standard_errors Whether to compute the reference standard errors.
#' @return A list: `reference`, a `harness_reference`, and `errors`, from
#'   [two_level_standard_errors()] or `NULL`.
.multilevlca_reference <- function(reference, analysis, items, id, fit,
                                   standard_errors = TRUE) {
  levels <- fit$categorical_levels[items]
  stopifnot("multilevLCA's mPhi is P(item = 1), the second of the fit's categories" =
              all(vapply(levels, identical, logical(1), c("0", "1"))))
  profiles <- ncol(reference$mPhi)
  rho <- stats::setNames(lapply(items, function(item) {
    positive <- as.numeric(reference$mPhi[sprintf("P(%s|C)", item), ])
    matrix(c(1 - positive, positive), profiles, dimnames = list(NULL, levels[[item]]))
  }), items)
  delta <- as.numeric(reference$vOmega)
  gamma <- t(unname(reference$mPi))
  codes <- .two_level_codes(analysis, levels)
  groups <- as.character(analysis[[id]])

  at_reference <- two_level_evaluate(delta, gamma, rho, codes, groups)
  log_likelihood <- utils::tail(as.numeric(reference$LLKSeries), 1L)
  .two_level_check("multilevLCA's log-likelihood", at_reference$log_likelihood,
                   log_likelihood, 1e-9 * abs(log_likelihood))
  # multilevLCA reports the rows sorted by group (a stable sort of the group
  # factor) and the groups in sorted order; the tolerances below allow the one
  # EM step by which its reported posteriors trail its estimates.
  sorted <- order(as.numeric(factor(analysis[[id]])))
  .two_level_check("multilevLCA's marginal profile posteriors",
                   at_reference$posteriors[sorted, ],
                   reference$mPMsumX[, sprintf("C%d", seq_len(profiles))], 1e-6)
  group_ids <- rownames(reference$vLLK)
  .two_level_check("multilevLCA's per-group log-likelihoods",
                   at_reference$group_log_likelihood[group_ids] / abs(reference$vLLK[, 1L]),
                   reference$vLLK[, 1L] / abs(reference$vLLK[, 1L]), 1e-8)
  # mPW is floored at 1e-5 and renormalised before it is returned (a group
  # posterior of 1e-57 is reported as 1e-5 / (1 + 1e-5) with two group
  # classes), but not in every row: on dataIEA a posterior of 8.3e-19 (CHL) is
  # reported as it is while smaller ones are floored. Each group is checked
  # against ours either as it is or put through that display step.
  ours <- unname(at_reference$group_posteriors[group_ids, , drop = FALSE])
  theirs <- unname(reference$mPW)
  displayed <- pmax(ours, 1e-5)
  displayed <- displayed / rowSums(displayed)
  row_gap <- pmin(apply(abs(ours - theirs), 1L, max), apply(abs(displayed - theirs), 1L, max))
  .two_level_check("multilevLCA's group posteriors, as is or after its display floor",
                   row_gap, rep(0, length(row_gap)), 1e-6)

  # The package bounds response probabilities at `min_probability`; the
  # reference is projected onto that bound, as in the glca suite, and every
  # posterior is evaluated there. On these data no probability is below it,
  # so the projection is the identity.
  bounded <- two_level_floor(rho, fit$min_probability)
  at_bound <- two_level_evaluate(delta, gamma, bounded, codes, groups)
  errors <- if (standard_errors) two_level_standard_errors(delta, gamma, bounded, codes, groups)
  list(reference = harness_reference(
         log_likelihood = log_likelihood, n_parameters = nrow(reference$parvec),
         posteriors = at_bound$posteriors, groups = groups, responses = bounded,
         standard_errors = if (standard_errors) list(probability = errors$observed$probability),
         profile_probabilities = gamma, group_probabilities = delta,
         group_posteriors = at_bound$group_posteriors,
         group_log_likelihood = at_bound$group_log_likelihood,
         row_conditional_entropy = at_bound$row_conditional_entropy,
         data = analysis[items]),
       errors = errors)
}

#' Published two-level results for the ICCS 2016 citizenship data
#'
#' Lyrvall et al. (2025) report a five-class, two-group-class model for twelve
#' binary citizenship-norm items measured on students within countries.
#' multilevLCA deletes incomplete rows, leaving 83,060 of 90,221.
#'
#' multilpa is scored at multilevLCA's estimates rather than re-optimized. That
#' is both faster and a stronger claim: it tests whether the two packages write
#' down the same likelihood and count the same free parameters, instead of
#' testing whether two different multi-start schemes reach the same mode on a
#' five-by-two landscape.
#'
#' Two multilevLCA runs, because two claims are made, as in the mclust suite.
#' The published values come from its default tolerance (1e-8), and are checked
#' against multilpa scored at that run's estimates. Every output is compared
#' with multilpa scored at a refined run's estimates, where multilevLCA's
#' reported posteriors, which trail its estimates by one EM step, are close
#' enough to the values at its estimates to check this file's evaluation
#' against them.
#'
#' @return A `harness_case`.
.multilevlca_published <- function() {
  observed <- get(utils::data(list = "dataIEA", package = "multilevLCA",
                              envir = environment()), envir = environment())
  items <- c("obey", "rights", "local", "work", "envir", "vote",
             "history", "respect", "news", "protest", "discuss", "party")
  complete <- observed[stats::complete.cases(observed[, c(items, "COUNTRY")]), , drop = FALSE]
  as_published <- .multilevlca_fit(observed, items, 5L, "COUNTRY", 2L, tol = 1e-8)
  evaluated <- .multilevlca_scored(as_published, complete, items)
  label <- "dataIEA 5 profiles x 2 group classes"
  against_package <- rbind(
    compare_values(sprintf("%s: analysis sample size", label),
                   83060, evaluated$n_informative, tolerance = 0.5),
    compare_values(sprintf("%s: log-likelihood at multilevLCA's estimates", label),
                   utils::tail(as.numeric(as_published$LLKSeries), 1L),
                   as.numeric(logLik(evaluated)), tolerance = 1e-4)
  ) |> transform(source = sprintf("multilevLCA %s", utils::packageVersion("multilevLCA")),
                 dataset = "dataIEA", precision = "machine")
  against_paper <- compare_printed(
    sprintf("%s: published %s", label,
            c("BIC, individual sample size", "BIC, group sample size", "AIC")),
    c(869122.92, 868554.62, 868479.34),
    c(evaluated$bic_individual, evaluated$bic_groups, evaluated$aic),
    digits = 2) |>
    transform(source = "Lyrvall et al. (2025) Multivar Behav Res 60(4):731-747",
              dataset = "dataIEA", precision = "printed")

  # Run on with `tol = 1e-15`, where it stops by itself after 353 steps with
  # its posteriors 3.2e-8 behind its estimates; every reference value is
  # evaluated at those estimates, which are what multilpa is scored at, so
  # this is enough here. (A fixed 2000 steps, as on dataTOY, would take
  # three minutes on 83,060 rows.)
  refined <- .multilevlca_fit(observed, items, 5L, "COUNTRY", 2L, tol = 1e-15)
  scored <- .multilevlca_scored(refined, complete, items)
  built <- .multilevlca_reference(refined, complete, items, "COUNTRY", scored,
                                  standard_errors = FALSE)
  every_output <- compare_fit(
    sprintf("%s, scored at multilevLCA's refined estimates", label), scored, built$reference,
    exempt = c(`responses$probability_standard_error` = paste(
      "scored at multilevLCA's estimates without iterating (max_iter = 0), so the fit is",
      "not converged and the package reports no standard errors by design",
      "(?parameter_inference: latents_no_converge); they are compared on dataTOY")))
  every_output$comparisons <- transform(
    every_output$comparisons,
    source = sprintf("multilevLCA %s", utils::packageVersion("multilevLCA")),
    dataset = "dataIEA", precision = "machine")
  stack_cases(list(against_package, against_paper, every_output))
}

#' Score multilpa at a multilevLCA solution without iterating
#' @param reference A `multiLCA` result.
#' @param complete The rows multilevLCA analysed.
#' @param items The binary indicators.
#' @return A `multilpa` fit whose estimates are multilevLCA's.
.multilevlca_scored <- function(reference, complete, items) {
  # mPhi holds P(item = 1 | class); multilpa carries both categories in sorted
  # level order, so category one is the zero response.
  responses <- lapply(items, function(item) {
    positive <- as.numeric(reference$mPhi[sprintf("P(%s|C)", item), ])
    # `cbind()` names the second column after the expression it was given, and
    # multilpa refuses a start whose response labels do not match the data's
    # own category levels. The matrix carries no labels.
    unname(cbind(1 - positive, positive))
  })
  multilpa(complete, vars = items, id = "COUNTRY",
           n_profiles = ncol(reference$mPhi), n_group_classes = nrow(reference$vOmega),
           categorical = items, n_starts = 1L, max_iter = 0L,
           start = starting_values(list(
             profile_probabilities = t(unname(reference$mPi)),
             group_probabilities = as.numeric(reference$vOmega),
             response_probabilities = responses)))
}
