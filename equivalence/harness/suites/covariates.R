# One-step membership-covariate models (`multilpa(profile_covariates =,
# group_covariates =, profile_slopes =)`) against four independent programs.
#
# Until this suite the covariate model entered no ledger: reference_tables()
# built `profile_probabilities` and never `coefficients`, so a wrong sign
# convention, reference class or standard error of a covariate effect -- the
# result a user reports -- would have gone unseen. Six fits on three simulated
# fixtures, each compared on every output through `compare_fit()`, and the
# table `parameter_inference()` returns through `compare_verb()`:
#
#   fit                         reference                          package call
#   single level, x             poLCA 1.6 (Linzer & Lewis 2011)    profile_covariates = "x"
#   two level, x and z          glca 1.4 (Kim et al. 2022),        profile_covariates = c("x", "z"),
#                                 coeff.inv = TRUE                   profile_slopes = "shared"
#   two level, x                glca, coeff.inv = FALSE            profile_slopes = "group_class"
#   two level, x; z on groups   multilevLCA 2.1 one-step           profile_slopes = "group_class",
#                                 (Lyrvall et al. 2025)              group_covariates = "z"
#   two level, x; z on groups   multilevLCA one-step,              profile_slopes = "shared",
#                                 fixedslopes = TRUE                 group_covariates = "z"
#   Poisson, Poisson, Gaussian  depmixS4 1.5 mix(prior = ~ x)      count = c("y1", "y2"),
#                                 (Visser & Speekenbrink 2010)       profile_covariates = "x"
#
# Correspondence of the parameterisations (each read into a full logit column
# per class, against any baseline; reference_tables() re-expresses them
# against the package's last class):
#
#   poLCA      `coeff`, terms x (R - 1): log P(class r) / P(class 1).
#   glca       `beta$Level1$Cluster w`, terms x (C - 1): log P(class c | w) /
#              P(class C), the intercept and level-1 slopes of cluster w
#              (`coeff.inv = TRUE` holds the slopes equal across clusters).
#              `beta$Level2`: the slopes of the level-2 covariates, which glca
#              puts in the CLASS logit, shared across clusters (`glca_em()`:
#              `exp(x %*% beta1 + z %*% beta2)`), not in the cluster logit; its
#              cluster probabilities `delta` take no covariate. A glca level-2
#              covariate is therefore a package profile covariate constant
#              within groups, not a group covariate. The package's group
#              covariates are compared against multilevLCA instead.
#   multilevLCA `cGamma[, , m]` (fixedslopes = FALSE), terms x (T - 1):
#              log P(C_t | G_m) / P(C_1); `mGamma_fixslope` (fixedslopes =
#              TRUE) holds the intercept of each G_m and one slope; `mAlpha`,
#              terms x (M - 1): log P(G_m) / P(G_1) on the higher-level
#              covariates. `mPhi`: P(item = 1 | C).
#   depmixS4   `prior = ~ x`: multinomial logits against the first state; the
#              Poisson response's intercept is log(mean), the Gaussian's
#              intercept and sd are the mean and the square root of the
#              variance.
#
# Each program's solution is evaluated by the base-R likelihood in
# covariates-reference.R, which must reproduce the program's own
# log-likelihood, posteriors and (where reported) group posteriors and group
# log-likelihoods to 1e-9 before any of its values is used. That code then
# supplies what the programs do not report: the within-class row entropies
# the criteria of a two-level model charge, and the standard errors (inverse
# numDeriv Hessian in a centred, unbounded chart).
#
# References are refined to the fixed point by the program's own EM where it
# can be continued (poLCA's E step, M step and Newton step for the logits,
# glca restarted from its own solution, depmixS4's EM), and by Newton steps
# on the base-R likelihood for multilevLCA, whose compiled EM cannot be
# restarted from its estimates; the displacement is reported. The package fits
# with `tol = 1e-16`, the "stop when the likelihood stops rising" rule the
# polca suite uses for the same reason: at 1e-14 a two-level fit stopped
# 1.7e-6 from the maximum along the flat group-logit direction (standard error
# 0.41), outside the 1e-6 the coefficients are compared at.
#
# Linzer, D. A., & Lewis, J. B. (2011). poLCA: An R package for polytomous
#   variable latent class analysis. Journal of Statistical Software, 42(10).
# Kim, Y., Jeon, S., Chang, C., & Chung, H. (2022). glca: An R package for
#   multiple-group latent class analysis. Applied Psychological Measurement,
#   46(5), 439-441.
# Lyrvall, J., Di Mari, R., Bakk, Z., Oser, J., & Kuha, J. (2025). Multilevel
#   latent class analysis: State-of-the-art methodologies and their
#   implementation in the R package multilevLCA. Multivariate Behavioral
#   Research.
# Visser, I., & Speekenbrink, M. (2010). depmixS4: An R package for hidden
#   Markov models. Journal of Statistical Software, 36(7).

source(file.path("equivalence", "harness", "covariates-reference.R"), local = environment())

suite_covariates <- function() {
  require_suite_packages(
    "poLCA", "glca", "multilevLCA", "depmixS4", "numDeriv", "MASS",
    reason = paste("the four independent covariate-model references (single level,",
                   "two level with class covariates, two level with group covariates,",
                   "count and Gaussian indicators) and the reference observed information"))
  two_level <- .covariates_two_level_data()
  stack_cases(list(
    .covariates_polca_case(.covariates_single_level_data()),
    .covariates_glca_case(two_level, slopes = "shared", profile_covariates = c("x", "z")),
    .covariates_glca_case(two_level, slopes = "group_class", profile_covariates = "x"),
    .covariates_multilevlca_case(two_level, slopes = "group_class"),
    .covariates_multilevlca_case(two_level, slopes = "shared"),
    .covariates_depmixs4_case(.covariates_count_data())))
}

# ---------------------------------------------------------------------------
# Fixtures (base R, fixed seeds)
# ---------------------------------------------------------------------------

#' Single level: 600 rows, five binary items, two classes, slope 1.2 on x
#' @return A data frame: `y1`..`y5` (factors, levels "1", "2"), `x`, `unit`.
.covariates_single_level_data <- function() {
  set.seed(11L)
  n <- 600L
  x <- stats::rnorm(n)
  class <- ifelse(stats::runif(n) < stats::plogis(-0.3 + 1.2 * x), 2L, 1L)
  .covariates_items(class, x, unit = sprintf("unit_%03d", seq_len(n)))
}

#' Five binary items from the classes, coded 1/2, beside the covariates
#' @param class The true class of each row.
#' @param ... Further columns.
#' @return A data frame.
.covariates_items <- function(class, ...) {
  probability <- rbind(c(0.8, 0.7, 0.8, 0.2, 0.3), c(0.2, 0.3, 0.25, 0.8, 0.7))
  items <- vapply(seq_len(5L), function(j) {
    ifelse(stats::runif(length(class)) < probability[class, j], 1L, 2L)
  }, integer(length(class)))
  colnames(items) <- sprintf("y%d", seq_len(5L))
  out <- data.frame(items, ...)
  out[colnames(items)] <- lapply(out[colnames(items)], factor, levels = c("1", "2"))
  out
}

#' Two level: 40 groups of 15 with a genuine group-class effect
#'
#' Group class 1 or 2 with probability `plogis(z_g)` (delta = (0.5, 0.5) at
#' z = 0), class prevalence 0.25 against 0.75 by group class, and a slope on
#' x of 1.6 in group class 1 and 0.8 in group class 2. Without the
#' group-class effect the group classes collapse, as they did in the probe.
#' Groups are labelled `school_01`.. so that a group label is never a row
#' number.
#' @return A data frame: `y1`..`y5`, `x`, `z` (constant within group), `g`.
.covariates_two_level_data <- function() {
  set.seed(12L)
  n_groups <- 40L
  group <- rep(seq_len(n_groups), each = 15L)
  z <- stats::rnorm(n_groups)
  group_class <- ifelse(stats::runif(n_groups) < stats::plogis(z), 1L, 2L)
  x <- stats::rnorm(length(group))
  logit <- stats::qlogis(c(0.25, 0.75))[group_class[group]] +
    c(1.6, 0.8)[group_class[group]] * x
  class <- ifelse(stats::runif(length(group)) < stats::plogis(logit), 2L, 1L)
  .covariates_items(class, x = x, z = z[group], g = sprintf("school_%02d", group))
}

#' Counts and a Gaussian: 400 rows, two Poisson indicators and one Gaussian
#' @return A data frame: `y1`, `y2` (counts), `w`, `x`, `unit`.
.covariates_count_data <- function() {
  set.seed(3L)
  n <- 400L
  x <- stats::rnorm(n)
  class <- ifelse(stats::runif(n) < stats::plogis(0.2 + x), 2L, 1L)
  data.frame(y1 = stats::rpois(n, c(2, 7)[class]), y2 = stats::rpois(n, c(5, 1)[class]),
             w = stats::rnorm(n, c(0, 2)[class]), x = x,
             unit = sprintf("unit_%03d", seq_len(n)))
}

#' The package fit every case makes
#'
#' `tol = 1e-16` makes the stopping rule "the likelihood stopped rising" (see
#' the header); five starts, as the references use five.
#' @return A `multilpa_covariates` fit.
.covariates_fit <- function(data, vars, id, n_group_classes, profile_covariates,
                            group_covariates = character(), slopes = "shared",
                            categorical = character(), count = character()) {
  fit <- multilpa(data, vars = vars, id = id, n_profiles = 2L,
                  n_group_classes = n_group_classes, categorical = categorical,
                  count = count, profile_covariates = profile_covariates,
                  group_covariates = group_covariates, profile_slopes = slopes,
                  n_starts = 5L, seed = 1L, tol = 1e-16, max_iter = 20000L)
  stopifnot("the package fit must converge" = isTRUE(fit$converged))
  fit
}

# ---------------------------------------------------------------------------
# The comparison every case shares
# ---------------------------------------------------------------------------

#' Compare one fit with one reference model on every output
#'
#' @param label Label for the rows.
#' @param fit The package fit.
#' @param model The reference, as a covariate model (covariates-reference.R)
#'   in the reference's class order, at its refined maximum.
#' @param context From [covariate_context()].
#' @param slopes `"shared"` or `"group_class"`.
#' @param n_parameters The reference program's own parameter count.
#' @param data The fitted data, for the residual tables.
#' @param vars The indicators.
#' @param grouped Whether the fit has a group level (several group classes).
#' @param extra Further comparison rows (outside the fit's tables).
#' @return A `harness_case`.
.covariates_compare <- function(label, fit, model, context, slopes, n_parameters, data,
                                vars, grouped, extra = NULL) {
  model <- .covariates_align(model, fit, context)
  at <- covariate_evaluate(model, context)
  errors <- covariate_standard_errors(model, context, slopes)
  se <- errors$standard_errors
  by_profile <- function(block, parameter, term_of) {
    if (is.null(block)) return(NULL)
    keys <- sprintf("measurement|profile_%d|%s|%s", rep(seq_len(nrow(block)), ncol(block)),
                    rep(term_of(colnames(block)), each = nrow(block)), parameter)
    matrix(unname(se[keys]), nrow(block), dimnames = dimnames(block))
  }
  standard_errors <- Filter(Negate(is.null), list(
    probability = if (length(model$responses) > 0L) {
      lapply(stats::setNames(nm = names(model$responses)), function(item) {
        by_profile(model$responses[[item]], "response", function(c) paste0(item, ":", c))
      })
    },
    mean = by_profile(model$means, "mean", identity),
    variance = by_profile(model$variances, "variance", identity),
    count_mean = by_profile(model$count_means, "count_mean", identity)))
  reference <- harness_reference(
    log_likelihood = at$log_likelihood, n_parameters = n_parameters,
    posteriors = at$posteriors,
    groups = if (grouped) context$group_ids[context$group_index],
    responses = if (length(model$responses) > 0L) model$responses,
    means = model$means, variances = model$variances, count_means = model$count_means,
    standard_errors = standard_errors,
    coefficients = Filter(Negate(is.null), covariate_coefficients(model, slopes)),
    group_posteriors = at$group_posteriors,
    group_log_likelihood = at$group_log_likelihood,
    row_conditional_entropy = if (grouped) at$row_conditional_entropy,
    data = data[vars])
  extra <- rbind(
    extra,
    # The default for the model row is 1e-6 relative to |log L|, a gap of 2e-3
    # here. Two EMs at one fixed point agree to second order in their
    # parameter gap (1e-7), so the likelihood is claimed to 1e-8 absolute.
    compare_values(sprintf("%s: maximised log-likelihood", label), at$log_likelihood,
                   as.numeric(logLik(fit)), tolerance = 1e-8),
    compare_values(sprintf("%s: free parameters in the reference's chart", label),
                   errors$n_parameters, fit$n_parameters, tolerance = 0.5),
    # The reference must be the maximum of the base-R likelihood to well
    # inside the 1e-6 its estimates are compared at: one Newton step from it
    # (inverse information times score) moves no chart coordinate by more
    # than 1e-7. A raw score is not the measure, since its size depends on
    # the coordinate's curvature (a Poisson log mean's information is 1,400,
    # a group logit's 1.2).
    compare_values(sprintf("%s: Newton distance from the reference to the base-R maximum", label),
                   0, errors$newton_step, tolerance = 1e-7))
  case <- compare_fit(
    label, fit, reference,
    # Standard errors come from numerical second derivatives in both
    # programs; two such estimates agree to about 1e-6, and 1e-4 is the
    # harness's standing tolerance for every standard-error column.
    tolerance = c(`count_means$standard_error` = 1e-4),
    extra = extra)
  stack_cases(list(case, .covariates_inference(label, fit, model, se, slopes)))
}

#' Put a reference model's classes in the fit's order
#'
#' Profiles are matched on the measurement parameters (response
#' probabilities, Gaussian means, log Poisson means), group classes on the
#' group posteriors, both by exhaustive search; the fit is read for this
#' matching only.
#'
#' @return The permuted model.
.covariates_align <- function(model, fit, context) {
  signature <- function(responses, means, counts) {
    blocks <- Filter(Negate(is.null), c(responses, list(means, if (!is.null(counts)) log(counts))))
    do.call(cbind, lapply(blocks, function(block) unname(as.matrix(block))))
  }
  ours <- signature(
    lapply(stats::setNames(nm = names(model$responses)), function(item) {
      fit$response_probabilities[[item]][, colnames(model$responses[[item]]), drop = FALSE]
    }),
    if (!is.null(model$means)) fit$means[, colnames(model$means), drop = FALSE],
    if (!is.null(model$count_means)) {
      fit$count_means[, colnames(model$count_means), drop = FALSE]
    })
  theirs <- signature(model$responses, model$means, model$count_means)
  model <- covariate_permute(model, align_classes(ours, theirs))
  if (length(model$profile) > 1L) {
    at <- covariate_evaluate(model, context)
    # Both in the order the groups first appear in the data.
    model <- covariate_permute(model, seq_len(ncol(model$profile[[1L]])),
                               align_classes(t(unname(fit$group_posteriors)),
                                             t(unname(at$group_posteriors))))
  }
  model
}

#' Compare `parameter_inference()` with the reference, every column
#'
#' Rows are the natural parameters of [covariate_natural()]: estimates from
#' the reference, standard errors from its observed information. The other
#' columns follow `?parameter_inference` from those two: `statistic` is
#' `estimate / standard_error`, `NA` for a parameter whose null value is on
#' the boundary of its space (a probability, a variance, a Poisson mean);
#' `p_value` is two-sided normal; `p_adjusted` equals it under the default
#' `adjust = "none"`; intervals are 95% Wald, formed on the logit scale for a
#' probability and the log scale for a variance and a Poisson mean (see the
#' note in the suite's report: the help page names variances only), and
#' `estimate +/- z * standard_error` otherwise.
#'
#' @param label Label of the fit.
#' @param fit The package fit.
#' @param model The aligned reference model.
#' @param se Its standard errors, named like [covariate_natural()].
#' @param slopes `"shared"` or `"group_class"`.
#' @return A `harness_case`.
.covariates_inference <- function(label, fit, model, se, slopes) {
  natural <- covariate_natural(model, slopes)
  keys <- c("level", "outcome", "term", "parameter")
  parts <- do.call(rbind, strsplit(names(natural), "|", fixed = TRUE))
  estimate <- unname(natural)
  standard_error <- unname(se[names(natural)])
  parameter <- parts[, 4L]
  z <- stats::qnorm(0.975)
  bounded <- parameter %in% c("response", "variance", "count_mean")
  statistic <- ifelse(bounded, NA_real_, estimate / standard_error)
  probability <- parameter == "response"
  positive <- parameter %in% c("variance", "count_mean")
  interval <- function(sign) {
    out <- estimate + sign * z * standard_error
    out[probability] <- stats::plogis(stats::qlogis(estimate[probability]) + sign * z *
      standard_error[probability] / (estimate[probability] * (1 - estimate[probability])))
    out[positive] <- exp(log(estimate[positive]) + sign * z *
                           standard_error[positive] / estimate[positive])
    out
  }
  p_value <- 2 * stats::pnorm(-abs(statistic))
  expected <- data.frame(level = parts[, 1L], outcome = parts[, 2L], term = parts[, 3L],
                         parameter = parameter, estimate = estimate,
                         standard_error = standard_error, statistic = statistic,
                         p_value = p_value, p_adjusted = p_value,
                         conf_low = interval(-1), conf_high = interval(1),
                         stringsAsFactors = FALSE)
  obtained <- keyed_table(as.data.frame(parameter_inference(fit)), keys)
  # Every column after the estimate is a function of the standard error, and
  # inherits its 1e-4 (a relative gap d in the error moves a statistic of 8
  # by 8 d, which is 8 d / 8 = d on the harness's relative scale).
  derived <- c("standard_error", "statistic", "p_value", "p_adjusted", "conf_low", "conf_high")
  # The harness compares a cell on |r - o| / max(|r|, 1), which is absolute
  # for a standard error below one: a response probability's error of 0.03
  # could be 0.3% off and pass. The design claims standard errors to 1e-4
  # relative, so that is stated once more, outside the ledger, on every
  # standard error of the table.
  matched <- obtained$standard_error[match(do.call(paste, c(expected[keys], sep = "|")),
                                           do.call(paste, c(obtained[keys], sep = "|")))]
  relative <- compare_values(
    sprintf("%s: largest relative gap of any parameter_inference() standard error (%d)", label,
            nrow(expected)),
    0, max(abs(matched / expected$standard_error - 1)), tolerance = 1e-4)
  compare_verb(sprintf("%s: parameter_inference()", label),
               list(parameter_inference = obtained),
               list(parameter_inference = expected),
               tolerance = stats::setNames(rep(1e-4, length(derived)),
                                           paste0("parameter_inference$", derived)),
               extra = relative)
}

#' Name a case's rows by their source
#' @param case A `harness_case`.
#' @param source,dataset Labels.
#' @return The case.
.covariates_label <- function(case, source, dataset) {
  case$comparisons <- transform(case$comparisons, source = source, dataset = dataset,
                                precision = "machine")
  case
}

# ---------------------------------------------------------------------------
# poLCA: single level
# ---------------------------------------------------------------------------

#' Single level against poLCA
#' @param data From [.covariates_single_level_data()].
#' @return A `harness_case`.
.covariates_polca_case <- function(data) {
  items <- sprintf("y%d", seq_len(5L))
  set.seed(20261027L)
  searched <- poLCA::poLCA(cbind(y1, y2, y3, y4, y5) ~ x, data = data, nclass = 2L,
                           nrep = 5L, maxiter = 20000L, tol = 1e-12, verbose = FALSE,
                           calc.se = FALSE)
  design <- stats::model.matrix(~ x, data)
  refined <- .covariates_polca_refine(searched, design)
  fit <- .covariates_fit(data, items, "unit", 1L, "x", categorical = items)
  context <- covariate_context(data, levels = fit$categorical_levels[items],
                               profile_covariates = "x")
  levels <- lapply(data[items], levels)
  model <- list(
    responses = stats::setNames(lapply(seq_along(items), function(j) {
      block <- unname(refined$probs[[j]])
      colnames(block) <- levels[[j]]
      block[, fit$categorical_levels[[items[[j]]]], drop = FALSE]
    }), items),
    # poLCA's logits are against its first class.
    profile = list(cbind(0, refined$coefficients)))
  dimnames(model$profile[[1L]]) <- list(colnames(design), NULL)
  at <- covariate_evaluate(model, context)
  covariate_check("poLCA's log-likelihood", at$log_likelihood, refined$llik,
                  1e-9 * abs(refined$llik))
  covariate_check("poLCA's posteriors", at$posteriors, refined$posterior, 1e-9)
  label <- "single level, profile covariate x (poLCA)"
  case <- .covariates_compare(
    label, fit, model, context, "shared", n_parameters = searched$npar, data = data,
    vars = items, grouped = FALSE,
    extra = rbind(
      compare_values(sprintf("%s: poLCA's log-likelihood at its fixed point", label),
                     refined$llik, as.numeric(logLik(fit)), tolerance = 1e-8),
      compare_values(sprintf("%s: poLCA's parameter count", label), searched$npar,
                     fit$n_parameters, tolerance = 0.5),
      .covariates_polca_opg(label, fit, refined, design)))
  .covariates_label(case, sprintf("poLCA %s", utils::packageVersion("poLCA")),
                    "single-level-covariate")
}

#' Continue poLCA's covariate EM from its own solution to the fixed point
#'
#' poLCA's loop with covariates (`poLCA()`, `S > 1`): its E step
#' (`poLCA.postClass.C`), M step for the probabilities (`poLCA.probHat.C`),
#' and one Newton step for the logits (`poLCA.dLL2dBeta.C`, `MASS::ginv`,
#' `poLCA.updatePrior`), until no parameter moves by more than 1e-14.
#'
#' @param searched A `poLCA` fit with covariates.
#' @param design Its design matrix (intercept and covariates).
#' @param max_steps Bound on the number of steps; reaching it is an error.
#' @param tolerance Largest parameter change at the fixed point.
#' @return A list: `probs`, `coefficients` (terms by classes 2..R, against
#'   class 1), `llik` (poLCA's expression), `posterior` and `prior` at the
#'   fixed point, `codes`, and `steps`.
.covariates_polca_refine <- function(searched, design, max_steps = 200000L,
                                     tolerance = 1e-14) {
  from <- function(name) utils::getFromNamespace(name, "poLCA")
  posterior_of <- from("poLCA.postClass.C")
  m_step <- from("poLCA.probHat.C")
  logit_derivatives <- from("poLCA.dLL2dBeta.C")
  update_prior <- from("poLCA.updatePrior")
  likelihoods <- from("poLCA.ylik.C")
  n <- nrow(searched$y)
  codes <- matrix(vapply(searched$y, as.integer, integer(n)), n,
                  dimnames = list(NULL, names(searched$y)))
  k <- length(searched$P)
  vp <- from("poLCA.vectorize")(searched$probs)
  b <- as.vector(searched$coeff)
  # EM is a fixed-point iteration: each step starts from the last one. Reduce()
  # rather than a loop with a break; once settled each element returns the
  # state untouched.
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    prior <- update_prior(state$b, design, k)
    posterior <- posterior_of(prior, state$vp, codes)
    vp <- state$vp
    vp$vecprobs <- m_step(posterior, codes, state$vp)
    derivatives <- logit_derivatives(posterior, prior, design)
    b <- as.vector(state$b + MASS::ginv(-derivatives$hess) %*% derivatives$grad)
    moved <- max(abs(c(vp$vecprobs - state$vp$vecprobs, b - state$b)))
    list(vp = vp, b = b, settled = moved <= tolerance, steps = index)
  }, seq_len(max_steps), init = list(vp = vp, b = b, settled = FALSE, steps = 0L))
  stopifnot("poLCA's EM must reach its fixed point within the step bound" = state$settled)
  prior <- update_prior(state$b, design, k)
  llik <- sum(log(rowSums(prior * likelihoods(state$vp, codes))) - log(.Machine$double.xmax))
  stopifnot("refining must not lose likelihood" = llik >= searched$llik - 1e-9)
  list(probs = from("poLCA.unvectorize")(state$vp),
       coefficients = matrix(state$b, ncol(design), k - 1L),
       llik = llik, posterior = posterior_of(prior, state$vp, codes), prior = prior,
       vp = state$vp, codes = codes, steps = state$steps)
}

#' poLCA's own coefficient standard errors against the package's OPG
#'
#' `poLCA.se()` inverts the cross-product of the per-observation scores (the
#' outer-product, BHHH, estimator) and reports the logits' errors from it, so
#' it is compared with `parameter_inference(vcov_type = "opg")`, as rows
#' outside the ledger. With two classes the package's logit (profile 1
#' against profile 2) is plus or minus poLCA's (class 2 against class 1),
#' whatever the class order, and has the same error.
#'
#' @return Comparison rows with `output = NA`.
.covariates_polca_opg <- function(label, fit, refined, design) {
  stopifnot("the sign argument above needs two classes" = length(refined$probs[[1L]][, 1L]) == 2L)
  errors <- utils::getFromNamespace("poLCA.se", "poLCA")(
    y = refined$codes, x = design, probs = refined$probs, prior = refined$prior,
    rgivy = refined$posterior)
  opg <- as.data.frame(parameter_inference(fit, vcov_type = "opg"))
  ours <- opg$standard_error[match(paste0("profile|profile_1|", colnames(design)),
                                   paste(opg$level, opg$outcome, opg$term, sep = "|"))]
  compare_block(sprintf("%s: poLCA's standard errors (OPG) of the logits, against the package's OPG",
                        label), as.vector(errors$b), ours, tolerance = 1e-4)
}

# ---------------------------------------------------------------------------
# glca: two level, covariates of class membership
# ---------------------------------------------------------------------------

#' Two level with class covariates against glca
#' @param data From [.covariates_two_level_data()].
#' @param slopes `"shared"` (glca `coeff.inv = TRUE`) or `"group_class"`
#'   (`coeff.inv = FALSE`).
#' @param profile_covariates The covariates of class membership; a level-2
#'   one (constant within groups) enters glca's class logit too (see header).
#' @return A `harness_case`.
.covariates_glca_case <- function(data, slopes, profile_covariates) {
  items <- sprintf("y%d", seq_len(5L))
  formula <- stats::as.formula(sprintf("glca::item(%s) ~ %s", paste(items, collapse = ", "),
                                       paste(profile_covariates, collapse = " + ")))
  invariant <- identical(slopes, "shared")
  searched <- glca::glca(formula, group = data$g, data = data, nclass = 2L, ncluster = 2L,
                         n.init = 5L, seed = 1L, maxiter = 20000L, eps = 1e-10,
                         verbose = FALSE, coeff.inv = invariant)
  # glca stops when no parameter moves by more than `eps`; its own EM is
  # continued from its own solution, as the glca suite refines. 2e-13 rather
  # than the glca suite's 1e-13: with coeff.inv = FALSE the Newton step for the
  # logits (through MASS::ginv) keeps a coefficient moving by about 1e-13 per
  # step indefinitely, and 50,000 steps at 1e-13 never stopped; at 2e-13 both
  # cases stop within 130 steps, and the coefficients then differ from those
  # at 5e-13 by 2e-13 at most (probe of 2026-10-27).
  refined <- glca::glca(formula, group = data$g, data = data, nclass = 2L, ncluster = 2L,
                        init.param = searched$param, seed = 1L, maxiter = 5000L,
                        eps = 2e-13, verbose = FALSE, coeff.inv = invariant)
  stopifnot("glca's refining EM must converge" = isTRUE(refined$convergence$converged))
  fit <- .covariates_fit(data, items, "g", 2L, profile_covariates, slopes = slopes,
                         categorical = items)
  context <- covariate_context(data, groups = data$g, levels = fit$categorical_levels[items],
                               profile_covariates = profile_covariates)
  param <- refined$param
  stopifnot("glca codes categories in factor-level order, the fit's order" =
              identical(unname(fit$categorical_levels[items]),
                        lapply(items, function(item) levels(data[[item]]))))
  level_two <- param$beta$Level2
  model <- list(
    responses = stats::setNames(lapply(items, function(item) {
      block <- unname(param$rho[[item]])
      colnames(block) <- levels(data[[item]])
      block
    }), items),
    # Class c against the last class within each cluster, the level-2 slopes
    # shared by every cluster.
    profile = lapply(seq_along(param$beta$Level1), function(w) {
      block <- cbind(rbind(param$beta$Level1[[w]], if (length(level_two) > 0L) level_two), 0)
      stopifnot("glca's terms must be the design's" =
                  identical(rownames(block), colnames(context$profile_design)))
      block
    }),
    group = matrix(log(as.numeric(param$delta)), 1L,
                   dimnames = list("(Intercept)", NULL)))
  at <- covariate_evaluate(model, context)
  covariate_check("glca's log-likelihood", at$log_likelihood, refined$gof$loglik,
                  1e-9 * abs(refined$gof$loglik))
  clusters <- as.matrix(refined$posterior$cluster)
  covariate_check("glca's group posteriors", at$group_posteriors[rownames(clusters), ],
                  clusters, 1e-9)
  reported <- match(unlist(lapply(refined$datalist$y, rownames)), row.names(data))
  stopifnot("glca reports every row once" =
              !anyNA(reported) && setequal(reported, seq_len(nrow(data))))
  covariate_check("glca's profile posteriors (marginal over clusters)", at$posteriors[reported, ],
                  do.call(rbind, lapply(refined$posterior$class, as.matrix)), 1e-9)
  label <- sprintf("two level, profile covariates %s, %s slopes (glca)",
                   paste(profile_covariates, collapse = " and "), slopes)
  case <- .covariates_compare(
    label, fit, model, context, slopes, n_parameters = refined$model$npar, data = data,
    vars = items, grouped = TRUE,
    extra = rbind(
      compare_values(sprintf("%s: glca's log-likelihood at its fixed point", label),
                     refined$gof$loglik, as.numeric(logLik(fit)), tolerance = 1e-8),
      compare_values(sprintf("%s: glca's parameter count", label), refined$model$npar,
                     fit$n_parameters, tolerance = 0.5)))
  .covariates_label(case, sprintf("glca %s", utils::packageVersion("glca")),
                    sprintf("two-level-covariate, %s slopes", slopes))
}

# ---------------------------------------------------------------------------
# multilevLCA: two level, covariates at both levels
# ---------------------------------------------------------------------------

#' Two level with class and group covariates against multilevLCA's one-step
#' estimator
#' @param data From [.covariates_two_level_data()].
#' @param slopes `"group_class"` (`fixedslopes = FALSE`) or `"shared"`
#'   (`fixedslopes = TRUE`).
#' @return A `harness_case`.
.covariates_multilevlca_case <- function(data, slopes) {
  items <- sprintf("y%d", seq_len(5L))
  coded <- data
  # multilevLCA takes items coded from 0: category "1" is 0, "2" is 1.
  coded[items] <- lapply(data[items], function(item) as.integer(item) - 1L)
  theirs <- multilevLCA::multiLCA(coded, Y = items, iT = 2L, id_high = "g", iM = 2L,
                                  Z = "x", Zh = "z", fixedpars = 0,
                                  fixedslopes = identical(slopes, "shared"), extout = TRUE,
                                  maxIter = 20000L, tol = 1e-14, NRtol = 1e-12,
                                  verbose = FALSE)
  fit <- .covariates_fit(data, items, "g", 2L, "x", group_covariates = "z", slopes = slopes,
                         categorical = items)
  context <- covariate_context(data, groups = data$g, levels = fit$categorical_levels[items],
                               profile_covariates = "x", group_covariates = "z")
  stopifnot("the items' categories must be 1 and 2" =
              all(vapply(fit$categorical_levels[items], identical, logical(1), c("1", "2"))))
  phi <- theirs$mPhi
  terms <- colnames(context$profile_design)
  # Logits of class 2 against class 1, one block per higher-level class.
  logits <- if (identical(slopes, "shared")) {
    # `mGamma_fixslope` labels its second row
    # "gamma(Intercept_C)+gamma(Intercept|G2)", but the value is the G2
    # deviation alone: G2's intercept is the first row plus the second. Read
    # as labelled, the likelihood at multilevLCA's estimates is 10.05 below the
    # one it reports; read so, it is the reported one to 2e-13, and the
    # package's own intercepts are the same two numbers (probe of 2026-10-27).
    fixed <- theirs$mGamma_fixslope
    list(matrix(c(fixed[1L, 1L], fixed[3L, 1L]), ncol = 1L),
         matrix(c(fixed[1L, 1L] + fixed[2L, 1L], fixed[3L, 1L]), ncol = 1L))
  } else {
    lapply(seq_len(2L), function(m) matrix(theirs$cGamma[, , m], ncol = 1L))
  }
  model <- list(
    responses = stats::setNames(lapply(items, function(item) {
      p <- phi[sprintf("P(%s|C)", item), ]
      cbind(`1` = 1 - p, `2` = p)
    }), items),
    profile = lapply(logits, function(block) {
      matrix(cbind(0, block), length(terms), dimnames = list(terms, NULL))
    }),
    group = matrix(cbind(0, theirs$mAlpha), 2L, dimnames = list(c("(Intercept)", "z"), NULL)))
  log_likelihood <- as.numeric(utils::tail(theirs$LLKSeries, 1L))
  at <- covariate_evaluate(model, context)
  covariate_check("multilevLCA's log-likelihood", at$log_likelihood, log_likelihood,
                  1e-9 * abs(log_likelihood))
  # Its group posteriors and group log-likelihoods are those of its last E
  # step, one update before the estimates it reports (their sum is the
  # reported likelihood to 2e-13, while single groups sit up to 1.2e-7 from
  # the evaluation at the estimates, with fixedslopes = TRUE). They are
  # checked at 1e-6, which that lag fits inside and a misread block does not
  # (a swapped higher-level class moves a group posterior by up to 1).
  groups_reported <- theirs$mPW[, c("G1", "G2"), drop = FALSE]
  covariate_check("multilevLCA's group posteriors (its last E step)",
                  at$group_posteriors[rownames(groups_reported), ], groups_reported, 1e-6)
  covariate_check("multilevLCA's group log-likelihoods (its last E step)",
                  at$group_log_likelihood[rownames(theirs$vLLK)], theirs$vLLK[, 1L], 1e-6)
  # The parameter vector multilevLCA reports, where it reports one
  # (`fixedslopes = FALSE`); otherwise the count its AIC implies,
  # AIC = -2 log L + 2 q. The two disagree where both exist: the AIC of the
  # 16-parameter group_class model implies 15 (probe of 2026-10-27), so the
  # AIC count is used only where it is the only one.
  n_parameters <- if (!is.null(theirs$parvec)) length(theirs$parvec) else
    round((theirs$AIC + 2 * log_likelihood) / 2)
  refined <- covariate_newton(model, context, slopes)
  label <- sprintf("two level, profile covariate x with %s slopes, group covariate z (multilevLCA)",
                   slopes)
  case <- .covariates_compare(
    label, fit, refined$model, context, slopes, n_parameters = n_parameters,
    data = data, vars = items, grouped = TRUE,
    extra = rbind(
      compare_values(sprintf("%s: multilevLCA's log-likelihood", label), log_likelihood,
                     as.numeric(logLik(fit)), tolerance = 1e-8),
      compare_values(sprintf("%s: multilevLCA's parameter count", label),
                     n_parameters, fit$n_parameters, tolerance = 0.5),
      # At tol = 1e-14 multilevLCA stops when its likelihood stops changing,
      # so it is short of the maximum by rounding (4e-13 at |log L| = 1848).
      # Along the flattest direction of the chart (information eigenvalue 1.2)
      # that admits a gap of sqrt(2 * 4e-13 / 1.2) = 8e-7, and the map to the
      # natural parameters stretches a chart step by at most the standard
      # deviation of x or z (about 1). 1e-5 is that with a margin of ten; a
      # move beyond it would mean the Newton steps reached another point.
      compare_values(sprintf("%s: Newton refinement's largest move from multilevLCA's estimates",
                             label), 0, refined$displacement, tolerance = 1e-5)))
  .covariates_label(case, sprintf("multilevLCA %s", utils::packageVersion("multilevLCA")),
                    sprintf("two-level-group-covariate, %s slopes", slopes))
}

# ---------------------------------------------------------------------------
# depmixS4: counts and a Gaussian, prior on x
# ---------------------------------------------------------------------------

#' Poisson and Gaussian indicators against depmixS4's mixture
#' @param data From [.covariates_count_data()].
#' @return A `harness_case`.
.covariates_depmixs4_case <- function(data) {
  vars <- c("y1", "y2", "w")
  specified <- depmixS4::mix(list(y1 ~ 1, y2 ~ 1, w ~ 1), data = data, nstates = 2L,
                             family = list(stats::poisson(), stats::poisson(),
                                           stats::gaussian()),
                             prior = ~ x, initdata = data)
  # Five random starts, as the package uses five.
  searched <- lapply(seq_len(5L), function(seed) {
    set.seed(seed)
    depmixS4::fit(specified, verbose = FALSE, emcontrol = depmixS4::em.control(
      tol = 1e-12, maxit = 5000L, random.start = TRUE))
  })
  best <- searched[[which.max(vapply(searched, function(one) {
    as.numeric(depmixS4::logLik(one))
  }, numeric(1)))]]
  refined <- .covariates_depmixs4_refine(best)
  fit <- .covariates_fit(data, vars, "unit", 1L, "x", count = c("y1", "y2"))
  context <- covariate_context(data, counts = c("y1", "y2"), gaussian = "w",
                               profile_covariates = "x")
  pars <- depmixS4::getpars(refined)
  # The prior's intercepts for every state, then its slopes; then one block
  # per state: log mean of y1, log mean of y2, mean and sd of w.
  prior <- matrix(pars[seq_len(4L)], 2L, 2L, byrow = TRUE,
                  dimnames = list(colnames(context$profile_design), NULL))
  blocks <- matrix(pars[-seq_len(4L)], 2L, 4L, byrow = TRUE)
  model <- list(
    count_means = matrix(exp(blocks[, c(1L, 2L)]), 2L, dimnames = list(NULL, c("y1", "y2"))),
    means = matrix(blocks[, 3L], 2L, dimnames = list(NULL, "w")),
    variances = matrix(blocks[, 4L]^2, 2L, dimnames = list(NULL, "w")),
    profile = list(prior))
  log_likelihood <- as.numeric(depmixS4::logLik(refined))
  at <- covariate_evaluate(model, context)
  covariate_check("depmixS4's log-likelihood", at$log_likelihood, log_likelihood,
                  1e-9 * abs(log_likelihood))
  covariate_check("depmixS4's posteriors", at$posteriors,
                  as.matrix(depmixS4::posterior(refined, type = "smoothing")), 1e-9)
  n_parameters <- attr(depmixS4::logLik(refined), "df")
  label <- "single level, Poisson and Gaussian indicators, profile covariate x (depmixS4)"
  case <- .covariates_compare(
    label, fit, model, context, "shared", n_parameters = n_parameters, data = data,
    vars = vars, grouped = FALSE,
    extra = rbind(
      compare_values(sprintf("%s: depmixS4's log-likelihood at its fixed point", label),
                     log_likelihood, as.numeric(logLik(fit)), tolerance = 1e-8),
      compare_values(sprintf("%s: depmixS4's parameter count", label), n_parameters,
                     fit$n_parameters, tolerance = 0.5)))
  .covariates_label(case, sprintf("depmixS4 %s", utils::packageVersion("depmixS4")),
                    "count-gaussian-covariate")
}

#' Continue depmixS4's EM from its own solution until the parameters settle
#'
#' depmixS4 stops on a relative likelihood change; as the depmixs4 suite
#' does, its own EM is run on, one step per `fit()` from the current
#' parameters, until no parameter moves by more than 1e-12.
#'
#' @param theirs A fitted depmixS4 mixture.
#' @param max_steps Bound on the number of steps; reaching it is an error.
#' @return The refined fit.
.covariates_depmixs4_refine <- function(theirs, max_steps = 20000L) {
  # At its fixed point depmixS4 warns "Log likelihood decreased on iteration
  # 1 from L to L'" whenever a step's likelihood differs from the last by
  # rounding (the two print identically to 15 digits). That one warning, with
  # a decrease below 1e-9 of |log L|, is muffled; any other, or a real
  # decrease, reaches the harness.
  rounding_decrease <- function(condition) {
    message_text <- conditionMessage(condition)
    if (!startsWith(message_text, "Log likelihood decreased")) return(FALSE)
    values <- as.numeric(regmatches(message_text, gregexpr("-?[0-9]+[.][0-9]+(e-?[0-9]+)?",
                                                           message_text))[[1L]])
    length(values) == 2L && values[[1L]] - values[[2L]] <= 1e-9 * abs(values[[1L]])
  }
  step <- function(object) {
    withCallingHandlers(
      depmixS4::fit(object, verbose = FALSE, emcontrol = depmixS4::em.control(
        maxit = 1L, tol = 0, crit = "absolute", random.start = FALSE)),
      warning = function(condition) {
        if (rounding_decrease(condition)) invokeRestart("muffleWarning")
      })
  }
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    moved <- step(state$model)
    change <- max(abs(depmixS4::getpars(moved) - depmixS4::getpars(state$model)))
    list(model = moved, settled = change < 1e-12, steps = index)
  }, seq_len(max_steps), init = list(model = theirs, settled = FALSE, steps = 0L))
  stopifnot("depmixS4's EM must settle within the step bound" = state$settled,
            "refining depmixS4's solution must not lose likelihood" =
              as.numeric(depmixS4::logLik(state$model)) >=
                as.numeric(depmixS4::logLik(theirs)) - 1e-9)
  state$model
}
