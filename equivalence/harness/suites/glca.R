# Nonparametric two-level latent class analysis against glca (Kim, Jeon, Chang
# & Chung 2022, Applied Psychological Measurement 46(5):439-441).
#
# glca is an independent implementation of the same model: its documentation
# states that the "level-2 (group-level) latent class is categorized by the
# prevalence of level-1 latent class for group variable", and it cites Vermunt
# (2003) directly. Its `measure.inv = TRUE` default is multilpa's shared
# measurement model. The parameter blocks map one to one:
#
#   glca delta  <-> multilpa group_probabilities     P(group class)
#   glca gamma  <-> multilpa profile_probabilities   P(profile | group class)
#   glca rho    <-> multilpa response_probabilities  P(item response | profile)
#
# Both datasets carry missing item responses, so this also checks that the two
# observed-data likelihoods agree, not only the complete-data ones. gss08 also
# has three rows with no observed item, which glca drops and multilpa keeps
# (unclassified, out of the individual sample size), and solutions with
# response probabilities on their bound.
#
# Every output of every fit is compared through `compare_fit()`. The reference
# is glca's refined solution: its estimates and likelihood as glca reports
# them, and its posteriors evaluated at those estimates by the base-R two-level
# likelihood in two-level-lca-reference.R, after that evaluation has been shown
# to reproduce glca's own likelihood, profile posteriors (`posterior$class`,
# the marginal posteriors P(c | y_g)) and group posteriors
# (`posterior$cluster`) to 1e-9. The same evaluation supplies what glca does not
# report: the per-group log-likelihood, the posteriors of the rows glca drops
# (their group's prior), and the standard errors.
#
# Standard errors. glca reports outer-product-of-gradients (BHHH) standard
# errors: `glca_score()` inverts `t(S) %*% S` of its per-group scores
# `GetUDScore()`. They are not used as a reference, because glca's response-
# probability scores are wrong on groups with missing responses: on nyts18
# (3 x 2) they equal minus the numerical per-group gradient exactly in the
# groups with complete responses and differ in every group with a missing one
# (largest gap 0.97), and they sum to 2.69 rather than zero at glca's own
# maximum (`glca$convergence$score`). Its delta standard error, 0.0912, comes
# from that matrix; the BHHH value from correct scores is 0.0878. On gss08
# there are 9 groups and 17 or 23 parameters, so the cross-product is singular
# and glca's generalized inverse returns zeros. Instead:
#
#   - `responses$probability_standard_error` (observed information) is compared
#     with the inverse of a numDeriv Hessian of the base-R likelihood at glca's
#     estimates;
#   - the class-probability standard errors of `parameter_inference()` are
#     compared with the same, as extra rows;
#   - `parameter_inference(vcov_type = "opg")` is compared with the BHHH
#     estimate from numDeriv per-group scores of the same likelihood, where the
#     package documents it as available (more groups than parameters), and its
#     refusal (`latents_too_few_groups`) is checked where it is not.

source(file.path("equivalence", "harness", "two-level-lca-reference.R"), local = environment())

suite_glca <- function() {
  require_suite_packages(
    "glca", "numDeriv",
    reason = "the independent Vermunt (2003) two-level implementation, the nyts18 and gss08 data, and the reference observed information")
  cases <- list(
    list(dataset = "nyts18", id = "SCH_ID",
         items = c("ECIGT", "ECIGAR", "ESLT", "EELCIGT", "EHOOKAH"),
         profiles = 3L, group_classes = 2L),
    list(dataset = "nyts18", id = "SCH_ID",
         items = c("ECIGT", "ECIGAR", "ESLT", "EELCIGT", "EHOOKAH"),
         profiles = 2L, group_classes = 2L),
    list(dataset = "gss08", id = "REGION",
         items = c("DEFECT", "HLTH", "RAPE", "POOR", "SINGLE", "NOMORE"),
         profiles = 3L, group_classes = 2L),
    list(dataset = "gss08", id = "REGION",
         items = c("DEFECT", "HLTH", "RAPE", "POOR", "SINGLE", "NOMORE"),
         profiles = 2L, group_classes = 3L)
  )
  stack_cases(lapply(cases, function(case) do.call(.glca_case, case)))
}

#' Run one glca comparison
#' @param dataset Name of a glca dataset.
#' @param id Name of the level-2 identifier column.
#' @param items Character vector of binary item names.
#' @param profiles Number of level-1 latent classes.
#' @param group_classes Number of level-2 latent classes.
#' @return A `harness_case`.
.glca_case <- function(dataset, id, items, profiles, group_classes) {
  observed <- get(utils::data(list = dataset, package = "glca",
                              envir = environment()), envir = environment())
  formula <- stats::as.formula(sprintf("glca::item(%s) ~ 1", paste(items, collapse = ", ")))
  reference <- glca::glca(formula, group = observed[[id]], data = observed,
                          nclass = profiles, ncluster = group_classes,
                          n.init = 20L, seed = 1L, maxiter = 20000L,
                          eps = 1e-10, verbose = FALSE)
  # glca stops when no parameter moves by more than `eps` in one EM step; where
  # the likelihood is flat that is still short of the fixed point (on nyts18,
  # 3 x 2, continuing moves a parameter by 1.4e-9). Its own EM is continued from
  # its own solution, as the mclust and poLCA suites refine theirs, so the
  # reference is glca's maximum and not where its stopping rule fired. The
  # refining `eps` is 1e-13 because rounding keeps gss08 (2 x 3) moving by
  # about 1e-14 per step indefinitely; at 1e-13 every case stops, and the
  # parameters are then the same as at 1e-15 to 4e-14 wherever 1e-15 stops.
  reference <- glca::glca(formula, group = observed[[id]], data = observed,
                          nclass = profiles, ncluster = group_classes,
                          init.param = reference$param, seed = 1L, maxiter = 20000L,
                          eps = 1e-13, verbose = FALSE)
  stopifnot("glca's refining EM must converge" = isTRUE(reference$convergence$converged))
  # The package's tightest stopping rule, as in the mclust suite. At 1e-12 its
  # posteriors stop 5e-7 from glca's refined ones on nyts18 (3 x 2), and a
  # bivariate-residual chi-square, which squares posterior-weighted
  # differences, then differs by 1.5e-6 relative; at 1e-15 the gaps are 1e-7
  # and 2e-7.
  fit <- multilpa(observed, vars = items, id = id,
                  n_profiles = profiles, n_group_classes = group_classes,
                  categorical = items, n_starts = 30L, seed = 1L,
                  tol = 1e-15, max_iter = 20000L, missing = "fiml")

  label <- sprintf("%s, %d profiles x %d group classes", dataset, profiles, group_classes)
  shortfall <- reference$gof$loglik - as.numeric(logLik(fit))
  # Two multi-start searches can stop at different modes; then nothing but the
  # shortfall is comparable, and the ledger says why, as in the mclust suite.
  if (abs(shortfall) > 1e-4) {
    rows <- compare_values(sprintf("%s: log-likelihood shortfall against glca", label),
                           0, max(shortfall, 0), tolerance = 1e-5)
    inventory <- fit_outputs(fit)
    case <- harness_case(rows, output_ledger(label, inventory, rows, exempt = stats::setNames(
      rep("glca stopped at a different mode, so its estimates describe another solution",
          length(unique(inventory$table))), unique(inventory$table))))
  } else {
    built <- .glca_reference(reference, observed, id, items, fit)
    alignment <- two_level_alignment(built$reference$responses,
                                     built$reference$profile_probabilities, fit)
    # BHHH needs more groups than parameters; the package documents a refusal
    # otherwise, which is checked instead.
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
      label, fit, built$reference,
      extra = rbind(
        compare_values(sprintf("%s: log-likelihood shortfall against glca", label),
                       0, max(shortfall, 0), tolerance = 1e-5),
        # glca's own criteria and entropy, computed by glca over the rows it
        # keeps, which are the package's informative rows.
        compare_values(sprintf("%s: AIC as glca reports it", label),
                       reference$gof$AIC, fit$aic, tolerance = 1e-4,
                       output = "information_criteria$aic"),
        compare_values(sprintf("%s: BIC as glca reports it (n = %d rows with a response)",
                               label, reference$model$N),
                       reference$gof$BIC, fit$bic_individual, tolerance = 1e-4,
                       output = "information_criteria$bic_individual"),
        compare_values(sprintf("%s: relative entropy as glca reports it", label),
                       reference$gof$entropy,
                       .glca_individual_relative_entropy(fit), tolerance = 1e-6,
                       output = "entropy$relative_entropy"),
        .with_output_columns(two_level_error_rows(label, fit, built$errors$observed,
                                                  alignment, "observed",
                                                  levels = c("profile", "group"))),
        .with_output_columns(opg)))
  }
  case$comparisons <- transform(case$comparisons,
                                source = sprintf("glca %s", utils::packageVersion("glca")),
                                dataset = dataset, precision = "machine")
  case
}

#' The individual-level relative entropy a fit reports
#' @param fit A two-level fit.
#' @return One number.
.glca_individual_relative_entropy <- function(fit) {
  entropy <- get_results(fit, "entropy")
  entropy$relative_entropy[entropy$level == "individuals"]
}

#' glca's solution in the harness's reference shape
#'
#' @param reference A refined `glca` fit.
#' @param observed The data it was fitted to.
#' @param id The group column.
#' @param items The indicators.
#' @param fit The package fit, for its category labels and probability bound
#'   only; no estimate is read from it.
#' @return A list: `reference`, a `harness_reference`, and `errors`, the
#'   reference standard errors from [two_level_standard_errors()].
.glca_reference <- function(reference, observed, id, items, fit) {
  levels <- fit$categorical_levels[items]
  stopifnot("glca codes categories in factor-level order, which is the fit's order" =
              identical(unname(levels), lapply(items, function(item) levels(observed[[item]]))))
  rho <- stats::setNames(lapply(items, function(item) {
    block <- unname(reference$param$rho[[item]])
    colnames(block) <- levels[[item]]
    block
  }), items)
  delta <- as.numeric(reference$param$delta)
  gamma <- unname(reference$param$gamma)
  codes <- .two_level_codes(observed, levels)
  groups <- as.character(observed[[id]])

  # Everything glca reports, reproduced at glca's estimates, before any of this
  # code's own values are used. glca's profile posteriors are the marginal
  # posteriors P(c | y_g): its M-step for rho (`UpRhoML`) weights by them.
  at_glca <- two_level_evaluate(delta, gamma, rho, codes, groups)
  .two_level_check("glca's log-likelihood", at_glca$log_likelihood, reference$gof$loglik,
                   1e-9 * abs(reference$gof$loglik))
  # glca returns each group's rows under their data row names, in group order,
  # and leaves out rows with no observed item.
  reported <- match(unlist(lapply(reference$datalist$y, rownames)), row.names(observed))
  stopifnot("glca's rows must be distinct rows of the data" =
              !anyNA(reported) && !anyDuplicated(reported))
  reported_posteriors <- do.call(rbind, lapply(reference$posterior$class, as.matrix))
  .two_level_check("glca's profile posteriors", at_glca$posteriors[reported, ],
                   reported_posteriors, 1e-9)
  group_posteriors <- as.matrix(reference$posterior$cluster)
  .two_level_check("glca's group posteriors",
                   at_glca$group_posteriors[rownames(group_posteriors), ], group_posteriors, 1e-9)
  stopifnot("glca drops only rows with no observed item" =
              setequal(setdiff(seq_len(nrow(observed)), reported),
                       which(rowSums(!is.na(codes)) == 0L)))

  # The package bounds response probabilities at `min_probability`; glca's
  # maximum puts some at zero (gss08). The reference is glca's solution
  # projected onto that bound, and every posterior and per-group likelihood is
  # evaluated there; without a probability below the bound (nyts18) that is
  # glca's own solution, and its posteriors are glca's to the 1e-9 checked
  # above. Evaluating at the projection matters for the bivariate residuals:
  # Cramer's V is a square root of a chi-square, so a table cell holding the
  # 1e-10 posterior mass the bound leaves gives V of order 1e-5 where glca's
  # zero gives V of 1e-65 (gss08, 3 x 2).
  bounded <- two_level_floor(rho, fit$min_probability)
  at_bound <- two_level_evaluate(delta, gamma, bounded, codes, groups)
  errors <- two_level_standard_errors(delta, gamma, bounded, codes, groups)
  list(reference = harness_reference(
         log_likelihood = reference$gof$loglik, n_parameters = reference$model$npar,
         posteriors = at_bound$posteriors, groups = groups, responses = bounded,
         standard_errors = list(probability = errors$observed$probability),
         profile_probabilities = gamma, group_probabilities = delta,
         group_posteriors = at_bound$group_posteriors,
         group_log_likelihood = at_bound$group_log_likelihood,
         row_conditional_entropy = at_bound$row_conditional_entropy,
         data = observed[items]),
       errors = errors)
}
