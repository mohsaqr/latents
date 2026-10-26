# Growth mixture models (mixture_regression() with `random =`) against
# lcmm::hlme (Proust-Lima, Philipps and Liquet 2017, Journal of Statistical
# Software 78(2):1-56), and the one-class limit against nlme::lme (Pinheiro and
# Bates 2000) and lme4::lmer (Bates, Maechler, Bolker and Walker 2015), both by
# maximum likelihood.
#
# Until this suite growth mixtures were compared on chosen quantities
# (equivalence/test-growth-*-lcmm.R). Here every output column of a growth fit
# (`get_results(fit, "all")`) is compared or exempted with its reason, through
# the ledger of outputs.R.
#
# hlme and the package's growth mixture are the same model:
#   y_i | class k ~ N(X_i beta_k, Z_i G_k Z_i' + sigma^2 I),
#   P(class k | w_i) = softmax(w_i' gamma)_k,
# hlme's `mixture` <-> the package's class-specific fixed effects, `classmb`
# <-> `membership`, `random` <-> `random` (hlme always adds the intercept), and
# one residual variance in both (`variance = "equal"`). hlme's covariance is
# shared by the classes (`nwg = FALSE`, the package's "equal") or proportional
# (`nwg = TRUE`, G_k = w_k^2 G with w_K = 1). Two package settings have no hlme
# counterpart and are not in this suite: class-specific residual variances
# (`variance = "varying"`; hlme has one `stderr`) and class-specific
# covariances of more than one random effect (`random_covariance = "varying"`
# with an intercept and a slope; hlme's classes differ only by a scalar). For a
# random intercept alone a proportional 1 x 1 covariance is any pair of
# class variances, so the package's "varying" is compared with `nwg = TRUE`
# there.
#
# Parameterisations: hlme's membership logits are against its last class and
# the package's against its first (gamma_k = l_k - l_1); hlme reports the
# residual standard deviation (`stderr`) and the random-effect covariance
# (`varcov`, upper triangle by column), and estimates the covariance through
# its Cholesky factor, so its `V` covers `cholesky`, not `varcov`.
#
# Every case makes two claims, so that a disagreement says whether the package
# computes wrongly or searches badly:
#   (A) one mode. hlme is started from the package's estimates (`B =`, mapped
#       by `.lcmm_start()`), read there with `maxiter = 0` (the same likelihood
#       at the same point), and iterated by its own Marquardt algorithm to its
#       fixed point with convB = convL = convG = 1e-12 (1e-10 on case 6, see
#       there). The free parameters (fixed effects, membership logits,
#       random-effect variances, SDs and correlations, shares, residual SD)
#       are compared with that fixed point in units of their standard errors,
#       to 1e-3; every other output with hlme read at the package's own
#       estimates, the same point: hlme's `loglik`, its posteriors (`pprob`),
#       its class-specific marginal and subject-specific predictions (`pred`),
#       its trajectories (`predictY`). What hlme does not report is written
#       here from its definition at those estimates: the class-specific random
#       effect predictions (checked against hlme's `predRE`, their posterior
#       mean), the criteria, the classification summaries and every standard
#       error. Standard errors come from the observed information of a
#       from-scratch base-R likelihood (`.lcmm_likelihood()`), whose Hessian
#       (at the package's point; at hlme's mode for hlme's own errors) is taken
#       by Richardson extrapolation (numDeriv) and carried to every reported
#       function by the delta method.
#   (B) search. The package's own starts against hlme's best over
#       `gridsearch()` (10 short runs from random draws around the one-class
#       fit, the best continued to convergence), continued here to hlme's own
#       fixed point. A shortfall is a finding about the package's search; a
#       higher package maximum is recorded as a row that is not an output of
#       the fit.
#
# The standard errors hlme reports (`V`, from its Marquardt algorithm) are
# compared as rows outside the ledger: see `.lcmm_own_error_rows()` for what
# they were found to be.
#
# Findings (2026-10-26, lcmm 2.2.1, nlme 3.1.168, lme4 1.1.38, package 0.97.1):
#   * Computation: every function of the parameters agrees with hlme at one
#     mode -- log likelihood at the package's estimates to 6e-15 relative,
#     posteriors, predictions and criteria to 1e-6 or better, the package's
#     standard errors with the reference observed information to 7.1e-7 --,
#     and the one-class fit with nlme and lme4 (log likelihood to 1.5e-15
#     relative, estimates and BLUPs to 7.9e-8). The fit in units of 1e4 is the
#     unit fit rescaled (estimates, posteriors and predictions to 6.5e-8,
#     standard errors to 2.7e-9).
#   * Stopping (an observation, not a defect): on membership-covariate the
#     package's quasi-Newton finish stops one Newton step short of hlme's
#     fixed point. Its largest score is 3.8e-5, its log likelihood 3.7e-11
#     below the fixed point (the reference information predicts 3.6e-11 for
#     that step), and one membership logit sits 1.55e-6 from it -- 7.5e-6 of
#     its standard error, so no estimate, interval or test can differ. `tol`
#     does not move it: 1e-12, 1e-14 and 1e-15 give the same point. This is
#     why the free parameters are compared with the fixed point in standard
#     errors (`.lcmm_parameter_distance()`), and everything computed from
#     them at the package's own point.
#   * hlme: its `loglik` uses pi in single precision (see
#     `.lcmm_log_likelihood()`), 1.391e-8 per row below the exact likelihood;
#     with convB = convL = convG = 1e-12 on an outcome in units of 1e4 it
#     diverges (conv = 4, NaN estimates) from a point it accepts at 1e-10;
#     its `V` is the inverse observed information by numerical second
#     differences, within 2.2e-3 of the Richardson-extrapolated one.

suite_lcmm <- function() {
  require_suite_packages(
    "lcmm", "numDeriv", "nlme", "lme4",
    reason = paste("the growth mixture reference (hlme), the reference observed",
                   "information, and the one-class mixed-model references"))
  cases <- .lcmm_cases()
  done <- lapply(cases, .lcmm_case)
  names(done) <- vapply(cases, `[[`, character(1), "name")
  # The outcome in units of 1e4 against the same data in units of one, both
  # fitted by the package: estimates must move by the units and nothing else.
  scaled <- done[["scale-1e4"]]
  scaled$case$comparisons <- rbind(
    scaled$case$comparisons,
    .lcmm_annotate(.lcmm_equivariance(done[["shared-covariance"]]$fit, scaled$fit, 1e4),
                   "latents (unit-scale fit)", "scale-1e4"))
  done[["scale-1e4"]] <- scaled
  stack_cases(lapply(done, `[[`, "case"))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' The cases: a generating model, what is fitted, and the hlme counterpart
#' @return A list of case specifications.
#' @noRd
.lcmm_cases <- function() {
  slope_effects <- matrix(c(0.8, 0.05, 0.05, 0.06), 2L)
  base <- list(times = 0:5, drop = 0, scale = 1, membership = ~ 1,
               membership_slopes = 0, random = ~ 1 + time, lcmm_random = ~ time,
               covariance = "equal", nwg = FALSE, effects = slope_effects,
               effect_scale = NULL, residual_sd = 0.6, convergence = 1e-12)
  case <- function(...) utils::modifyList(base, list(...))
  list(
    # 1. One class: a linear mixed model, also against nlme and lme4.
    case(name = "one-class", seed = 401L, n = 150L, k = 1L, intercepts = 2,
         slopes = 0.4, shares = 1),
    # 2. Random intercept, a class-specific variance (SD 0.5 against 1.2).
    case(name = "varying-intercept", seed = 402L, n = 200L, k = 2L,
         intercepts = c(1, 4), slopes = c(0.5, -0.3), shares = c(0.55, 0.45),
         random = ~ 1, lcmm_random = ~ 1, covariance = "varying", nwg = TRUE,
         effects = matrix(1), effect_scale = c(0.5, 1.2)),
    # 3. Random intercept and slope, one covariance shared by the classes.
    case(name = "shared-covariance", seed = 403L, n = 200L, k = 2L,
         intercepts = c(1.5, 4), slopes = c(0.6, -0.4), shares = c(0.6, 0.4)),
    # 4. Three classes whose shares depend on a person covariate x.
    case(name = "membership-covariate", seed = 404L, n = 240L, k = 3L,
         intercepts = c(0, 3, 6), slopes = c(0.6, 0, -0.5), shares = c(1, 1, 1) / 3,
         membership = ~ x, membership_slopes = c(0, 1, -1)),
    # 5. An unbalanced panel: occasions missing at random (40%), every person
    #    keeping at least one, so panels run from one to six occasions.
    case(name = "unbalanced", seed = 405L, n = 220L, k = 2L,
         intercepts = c(1.5, 4), slopes = c(0.6, -0.4), shares = c(0.6, 0.4),
         drop = 0.4),
    # 6. Case 3's data in units of 1e4. hlme's convergence criteria are held
    #    at 1e-10, not 1e-12: started from the package's solution with
    #    convB = convL = convG = 1e-12 hlme stops with conv = 4 ("Infinite
    #    parameters"), every estimate NaN and loglik -1e9, while at 1e-10 it
    #    converges (conv = 1) in 8 iterations (measured 2026-10-26, lcmm
    #    2.2.1). That is hlme's Marquardt step on variances of order 1e8,
    #    not a property of the model.
    case(name = "scale-1e4", seed = 403L, n = 200L, k = 2L,
         intercepts = c(1.5, 4), slopes = c(0.6, -0.4), shares = c(0.6, 0.4),
         scale = 1e4, convergence = 1e-10))
}

#' Evaluate an expression under a seed, leaving the caller's stream untouched
#' @noRd
.lcmm_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
    rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}

#' Simulate one case
#' @return A data frame sorted by person and time: `idn` (numeric, for hlme),
#'   `id` (zero-padded text, for the package, so both sort alike), `time`,
#'   `x` (a person covariate) and `y`.
#' @noRd
.lcmm_simulate <- function(case) {
  .lcmm_with_seed(case$seed, {
    n <- case$n
    k <- case$k
    x <- stats::rnorm(n)
    logits <- outer(x, rep_len(case$membership_slopes, k)) +
      matrix(log(case$shares), n, k, byrow = TRUE)
    prior <- exp(logits) / rowSums(exp(logits))
    class <- vapply(seq_len(n), function(i) sample.int(k, 1L, prob = prior[i, ]),
                    integer(1))
    q <- ncol(case$effects)
    effects <- matrix(stats::rnorm(n * q), n) %*% chol(case$effects)
    effects <- effects * (case$effect_scale %||% rep(1, k))[class]
    times <- case$times
    frame <- data.frame(idn = rep(seq_len(n), each = length(times)),
                        time = rep(times, n))
    person <- frame$idn
    slope_effect <- if (q > 1L) effects[person, 2L] else 0
    frame$y <- case$intercepts[class[person]] + case$slopes[class[person]] * frame$time +
      effects[person, 1L] + slope_effect * frame$time +
      stats::rnorm(nrow(frame), 0, case$residual_sd)
    keep <- stats::runif(nrow(frame)) >= case$drop | !duplicated(person)
    frame <- frame[keep, ]
    frame$y <- frame$y * case$scale
    frame$id <- sprintf("p%03d", frame$idn)
    frame$x <- x[frame$idn]
    row.names(frame) <- NULL
    stopifnot("every person keeps an occasion" = length(unique(frame$idn)) == n)
    frame
  })
}

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

#' Fit one case with both programs and compare
#' @return A list with `case` (a `harness_case`) and `fit`.
#' @noRd
.lcmm_case <- function(case) {
  frame <- .lcmm_simulate(case)
  label <- sprintf("%s, %d class%s", case$name, case$k, if (case$k > 1L) "es" else "")
  # Warnings are not suppressed: the harness counts them and records the first.
  fit <- mixture_regression(y ~ time, frame, n_classes = case$k, id = "id",
                            class_level = "group", membership = case$membership,
                            random = case$random, random_covariance = case$covariance,
                            variance = "equal", n_starts = 10L, seed = case$seed,
                            tol = 1e-12, max_iter = 5000L)
  setup <- .lcmm_setup(frame, case)
  search <- .lcmm_search(label, case, frame, fit)
  same <- .lcmm_same_mode(label, case, frame, setup, fit)
  rows <- rbind(.with_output_columns(same$comparisons), search,
                if (case$k == 1L) .lcmm_mixed_model_rows(label, frame, fit))
  source <- sprintf("lcmm %s", utils::packageVersion("lcmm"))
  list(case = harness_case(.lcmm_annotate(rows, source, case$name), same$ledger), fit = fit)
}

#' Add the columns the harness reads off every comparison row
#' @noRd
.lcmm_annotate <- function(rows, source, dataset) {
  rows <- .with_output_columns(rows)
  rows$source <- source
  rows$dataset <- dataset
  rows$precision <- "machine"
  rows
}

#' hlme's log likelihood on the exact Gaussian constant
#'
#' hlme's `loglik` sits below the likelihood of its own estimates by the same
#' amount per observed row, 1.391e-8: exactly what a Gaussian constant
#' -0.5 log(2 pi) taken with pi in single precision (3.14159274..., the float
#' nearest 3.14159265) costs, 0.5 log(fl(pi) / pi). The from-scratch
#' likelihood here and the package both use R's double `pi`; every case checks
#' the correction against the from-scratch likelihood at hlme's mode, to 1e-10
#' relative (see `.lcmm_same_mode()`), so it is a measured constant and not an
#' assumption. On 1,200 rows it is 1.67e-5, which without the correction reads
#' as the package's maximum lying 1.67e-5 above hlme's.
#'
#' @param model An hlme fit.
#' @param n_rows Observed rows it was fitted to.
#' @noRd
.lcmm_log_likelihood <- function(model, n_rows) {
  single_pi <- 3.1415927410125732
  model$loglik + 0.5 * n_rows * log(single_pi / pi)
}

#' Call hlme for a case
#'
#' @param start hlme's starting vector, or an hlme fit.
#' @param maxiter Iterations; 0 evaluates the likelihood at `start`.
#' @param tight Use convB = convL = convG = `case$convergence` (1e-12 unless
#'   the case says otherwise) instead of hlme's 1e-4.
#' @noRd
.lcmm_hlme <- function(case, frame, start, maxiter = 500L, tight = TRUE) {
  arguments <- list(fixed = y ~ time, random = case$lcmm_random, subject = "idn",
                    ng = case$k, data = frame, B = start, maxiter = maxiter,
                    verbose = FALSE)
  if (case$k > 1L) {
    arguments <- c(arguments, list(mixture = ~ time, classmb = case$membership,
                                   nwg = case$nwg))
  }
  if (tight) {
    arguments <- c(arguments, list(convB = case$convergence, convL = case$convergence,
                                   convG = case$convergence))
  }
  do.call(lcmm::hlme, arguments)
}

# ---------------------------------------------------------------------------
# Claim B: search
# ---------------------------------------------------------------------------

#' The package's own search against hlme's best over a grid search
#'
#' `gridsearch()` runs 10 short hlme fits (30 iterations) from random draws
#' around the one-class fit and continues the best at hlme's default
#' convergence (1e-4); that is then continued here from its own estimates with
#' convB = convL = convG = 1e-12 (1e-10 on case 6), so the reference is hlme's fixed point and
#' not where its default stopping rule fired.
#'
#' @return Comparison rows, none an output of the fit.
#' @noRd
.lcmm_search <- function(label, case, frame, fit) {
  hlme <- lcmm::hlme
  minit <- eval(bquote(hlme(y ~ time, random = .(case$lcmm_random), subject = "idn",
                            ng = 1, data = frame, verbose = FALSE)))
  best <- if (case$k == 1L) minit else .lcmm_with_seed(case$seed, {
    # gridsearch() reads its model unevaluated and calls it by name in this
    # frame, so the call is built around the local `hlme`.
    model <- bquote(hlme(y ~ time, mixture = ~ time, random = .(case$lcmm_random),
                         subject = "idn", ng = .(case$k), nwg = .(case$nwg),
                         classmb = .(case$membership), data = frame, verbose = FALSE))
    eval(bquote(lcmm::gridsearch(.(model), rep = 10, maxiter = 30, minit = minit)))
  })
  refined <- .lcmm_hlme(case, frame, best$best)
  shortfall <- .lcmm_log_likelihood(refined, nrow(frame)) - fit$log_likelihood
  rows <- list(
    compare_values(sprintf("%s: hlme's best search converged by its own criteria (conv = 1), %s",
                           label, c("at its default tolerance",
                                    sprintf("when continued to %g", case$convergence))),
                   c(1, 1), c(best$conv, refined$conv), tolerance = 0.5),
    compare_values(sprintf("%s: shortfall of the package's own search against hlme's best",
                           label), 0, max(shortfall, 0), tolerance = 1e-6))
  if (-shortfall > 1e-6) {
    rows <- c(rows, list(compare_values(
      sprintf("%s: package maximum above hlme's best (by %.3g; hlme's search)", label,
              -shortfall), -shortfall, -shortfall, tolerance = 1e-6)))
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Claim A: one mode
# ---------------------------------------------------------------------------

#' hlme's starting vector from the package's parameters
#'
#' Positional, in hlme's order: membership logits against the last class
#' (term by term, classes 1 to K-1), class-specific fixed effects (term by
#' term), the covariance's upper triangle by column (or, for a proportional
#' covariance, the last class's variance and each other class's scale), and the
#' residual standard deviation.
#' @noRd
.lcmm_start <- function(fit, case) {
  params <- fit$params
  k <- case$k
  membership <- if (k > 1L) {
    logits <- params$gamma - params$gamma[, k]
    as.vector(t(logits[, -k, drop = FALSE]))
  }
  covariances <- latents:::.growth_covariances(fit$spec, params)
  last <- covariances[[k]]
  random <- last[upper.tri(last, diag = TRUE)]
  scales <- if (isTRUE(case$nwg)) vapply(covariances[-k], function(g) {
    sqrt(g[1L, 1L] / last[1L, 1L])
  }, numeric(1))
  c(membership, as.vector(t(params$beta)), random, scales, sqrt(params$sigma2[[1L]]))
}

#' The natural parameters of an hlme estimate vector
#' @return A list with `beta` (terms by classes), `gamma` (membership terms by
#'   classes, column 1 zero), `covariances` (one matrix per class) and
#'   `sigma2`.
#' @noRd
.lcmm_natural <- function(best, model, case, setup) {
  counts <- model$N
  k <- case$k
  q <- ncol(setup$z_names)
  p_w <- length(setup$w_names)
  p_x <- length(setup$x_names)
  block <- function(i) best[sum(counts[seq_len(i - 1L)]) + seq_len(counts[[i]])]
  logits <- if (k > 1L) cbind(matrix(block(1L), p_w, k - 1L, byrow = TRUE), 0) else
    matrix(0, p_w, 1L)
  covariance <- matrix(0, q, q)
  covariance[upper.tri(covariance, diag = TRUE)] <- block(3L)
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  scales <- if (counts[[4L]] > 0L) c(block(4L), 1) else rep(1, k)
  list(beta = matrix(block(2L), p_x, k, byrow = TRUE,
                     dimnames = list(setup$x_names, NULL)),
       gamma = logits - logits[, 1L],
       covariances = lapply(scales, function(s) covariance * s^2),
       sigma2 = best[["stderr"]]^2)
}

#' Claim A: every output against hlme at one mode
#' @return A list with `comparisons` and `ledger`.
#' @noRd
.lcmm_same_mode <- function(label, case, frame, setup, fit) {
  start <- .lcmm_start(fit, case)
  at_package <- .lcmm_hlme(case, frame, start, maxiter = 0L, tight = FALSE)
  mode <- .lcmm_hlme(case, frame, start)
  stopifnot("hlme must take the starting vector in the order built here" =
              length(start) == length(mode$best))
  # Two references, each for what it isolates (the rule of
  # `.hard_depmixs4_same_mode()`): the free parameters against hlme's fixed
  # point, in standard errors (`.lcmm_parameter_distance()`); every function
  # of them -- posteriors, predictions, criteria, standard errors, intervals
  # -- against hlme and the from-scratch likelihood evaluated at the
  # package's own estimates (`at_package`, hlme with `maxiter = 0`), the same
  # point, so a stopping rule cannot pass for a computing error.
  reference <- .lcmm_reference(mode, case, setup)
  at_point <- .lcmm_reference(at_package, case, setup)
  expected <- .lcmm_expected(at_point, at_package, case, setup, fit)
  at_mode_tables <- .lcmm_expected(reference, mode, case, setup, fit)
  # At the package's point the parameter columns equal its estimates by
  # construction, so they are compared with the fixed point below instead.
  expected[names(.lcmm_parameters)] <- lapply(names(.lcmm_parameters), function(table) {
    out <- expected[[table]]
    if (!is.null(out)) out[names(.lcmm_parameters[[table]])] <- NULL
    out
  })
  expected <- Filter(Negate(is.null), expected)
  exempt <- .lcmm_exemptions(fit)
  obtained <- fit_tables(fit)
  # Tolerances: estimates 1e-6; standard errors, and the statistics,
  # p-values and intervals that inherit them, 1e-4 (a numerical Hessian here
  # against the package's numerical Jacobian of its analytic score).
  inference <- c("std_error", "statistic", "p_value", "conf_low", "conf_high")
  tolerance <- c(default = 1e-6,
                 stats::setNames(rep(1e-4, 2L * length(inference)),
                                 c(paste0("coefficients$", inference),
                                   paste0("membership$", inference))),
                 stats::setNames(rep(1e-4, 3L), paste0("random$", inference[c(1L, 4L, 5L)])),
                 stats::setNames(rep(1e-4, 3L), paste0("trajectories$", inference[c(1L, 4L, 5L)])),
                 `classes$share_std_error` = 1e-4, `classes$share_conf_low` = 1e-4,
                 `classes$share_conf_high` = 1e-4, `classes$residual_sd_std_error` = 1e-4)
  comparisons <- compare_tables(label, obtained, .drop_exempt(expected, exempt), tolerance)
  comparisons <- rbind(comparisons,
                       .lcmm_parameter_distance(label, obtained, at_mode_tables))
  modal <- if (case$k > 1L) at_package$pprob$class else rep(1L, setup$n)
  gain <- mode$loglik - at_package$loglik
  claims <- rbind(
    compare_values(sprintf("%s: hlme's log likelihood at the package's estimates", label),
                   .lcmm_log_likelihood(at_package, setup$n_rows), fit$log_likelihood,
                   tolerance = 1e-10, scale = "relative"),
    compare_values(sprintf("%s: log likelihood hlme gains from the package's solution",
                           label), 0, gain, tolerance = 1e-6),
    compare_values(sprintf("%s: hlme converged from the package's solution (conv = 1)",
                           label), 1, mode$conv, tolerance = 0.5),
    compare_values(sprintf(paste("%s: the from-scratch likelihood at hlme's mode equals",
                                 "hlme's loglik with the exact Gaussian constant"), label),
                   .lcmm_log_likelihood(mode, setup$n_rows), reference$log_likelihood,
                   tolerance = 1e-10, scale = "relative"),
    compare_values(sprintf(paste("%s: hlme's loglik shortfall per observed row against the",
                                 "from-scratch likelihood, 0.5 log(fl(pi) / pi)"), label),
                   0.5 * log(3.1415927410125732 / pi),
                   (reference$log_likelihood - mode$loglik) / setup$n_rows,
                   tolerance = 1e-11),
    compare_values(sprintf("%s: free parameters, against hlme's", label),
                   length(mode$best), fit$n_parameters, tolerance = 0.5),
    compare_values(sprintf(paste("%s: the from-scratch likelihood at the package's",
                                 "estimates equals the package's"), label),
                   at_point$log_likelihood, fit$log_likelihood,
                   tolerance = 1e-10, scale = "relative"),
    compare_block(sprintf("%s: modal class of every person, against hlme's pprob", label),
                  modal,
                  max.col(fit$expectation$posterior, ties.method = "first"),
                  tolerance = 0.5),
    compare_block(sprintf(paste("%s: posterior mean of the class-specific random-effect",
                                "predictions written here, against hlme's predRE"), label),
                  as.vector(as.matrix(at_package$predRE[-1L])),
                  as.vector(at_point$posterior_effects), tolerance = 1e-8),
    compare_values(sprintf(paste("%s: hlme's own %s (persons) at the package's estimates,",
                                 "with the exact Gaussian constant, against the package's"),
                           label, c("AIC", "BIC")),
                   c(at_package$AIC, at_package$BIC) -
                     2 * (.lcmm_log_likelihood(at_package, setup$n_rows) - at_package$loglik),
                   c(obtained$fit$aic, obtained$fit$bic),
                   tolerance = 1e-6, scale = "relative"),
    .lcmm_own_error_rows(label, mode, case, setup, reference))
  comparisons <- rbind(comparisons, .with_output_columns(claims))
  list(comparisons = comparisons,
       ledger = output_ledger(label, fit_outputs(obtained), comparisons, exempt))
}

#' Outputs hlme cannot speak to, with the reason
#' @noRd
.lcmm_exemptions <- function(fit) {
  flag <- paste("this package's warning flag for a degenerate random-effect covariance,",
                "with its own threshold (unit-tested); hlme has no such flag")
  out <- c(`fit$converged` = "this optimiser's own stopping flag",
           `fit$n_best_replicated` = "how many of this program's random starts reached the best value",
           `fit$random_boundary` = flag, `random$boundary` = flag)
  outputs <- fit_outputs(fit)$output
  out[names(out) %in% c(outputs, unique(sub("[$].*", "", outputs)))]
}

# ---------------------------------------------------------------------------
# The reference likelihood, written from the model's definition
# ---------------------------------------------------------------------------

#' Designs of a case, persons grouped by their pattern of occasions
#'
#' The fixed and random designs depend on time alone, so persons observed at
#' the same occasions share one marginal covariance per class: the likelihood
#' is evaluated one pattern at a time.
#' @noRd
.lcmm_setup <- function(frame, case) {
  stopifnot("the frame must be sorted by person and time" =
              !is.unsorted(frame$idn) && all(tapply(frame$time, frame$idn,
                                                    function(t) !is.unsorted(t))))
  persons <- unique(frame$id)
  index <- match(frame$id, persons)
  fixed <- stats::model.matrix(~ time, frame)
  random <- stats::model.matrix(case$random, frame)
  w <- stats::model.matrix(case$membership, frame[!duplicated(index), , drop = FALSE])
  key <- vapply(split(frame$time, index), paste, character(1), collapse = ",")
  groups <- split(seq_along(persons), factor(key, unique(key)))
  patterns <- lapply(groups, function(members) {
    rows <- which(index == members[[1L]])
    list(persons = members, x = fixed[rows, , drop = FALSE], z = random[rows, , drop = FALSE],
         y = matrix(frame$y[index %in% members], length(members), length(rows), byrow = TRUE))
  })
  list(patterns = patterns, persons = persons, n = length(persons), w = w,
       x_names = colnames(fixed), z_names = random[0L, , drop = FALSE],
       w_names = colnames(w), n_rows = nrow(frame), frame = frame)
}

#' Each person's log density, and predicted random effects, in one class
#' @return A list with `log_density` (one per person) and `effects` (persons
#'   by random effects: G Z' V^-1 (y - X beta)).
#' @noRd
.lcmm_class_density <- function(setup, beta, covariance, sigma2) {
  per_pattern <- lapply(setup$patterns, function(pattern) {
    m <- nrow(pattern$x)
    root <- chol(pattern$z %*% covariance %*% t(pattern$z) + diag(sigma2, m))
    residual <- pattern$y - matrix(as.vector(pattern$x %*% beta), nrow(pattern$y), m,
                                   byrow = TRUE)
    standardised <- backsolve(root, t(residual), transpose = TRUE)
    list(log_density = -0.5 * (m * log(2 * pi) + 2 * sum(log(diag(root))) +
                                 colSums(standardised^2)),
         effects = t(covariance %*% t(pattern$z) %*% backsolve(root, standardised)))
  })
  order <- unlist(lapply(setup$patterns, `[[`, "persons"))
  log_density <- numeric(setup$n)
  log_density[order] <- unlist(lapply(per_pattern, `[[`, "log_density"))
  effects <- matrix(0, setup$n, ncol(setup$z_names))
  effects[order, ] <- do.call(rbind, lapply(per_pattern, `[[`, "effects"))
  list(log_density = log_density, effects = effects)
}

#' The model's likelihood in the natural parameters, from scratch
#'
#' theta = (beta by class, membership logits of classes 2..K against class 1,
#' the covariance's lower triangle by column for each covariance block -- one
#' block shared, or one per class for a class-specific covariance --, the
#' residual variance). At a maximum the inverse observed information
#' transforms with the Jacobian of any reparameterisation, so its errors are
#' those of the natural parameters whatever scale either program optimises on.
#' @noRd
.lcmm_likelihood <- function(case, setup) {
  k <- case$k
  p_x <- length(setup$x_names)
  p_w <- length(setup$w_names)
  q <- ncol(setup$z_names)
  cells <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  blocks <- if (identical(case$covariance, "varying")) k else 1L
  sizes <- c(beta = p_x * k, gamma = p_w * (k - 1L), covariance = nrow(cells) * blocks,
             sigma2 = 1L)
  ends <- cumsum(sizes)
  at <- function(name) ends[[name]] - sizes[[name]] + seq_len(sizes[[name]])
  unpack <- function(theta) {
    values <- matrix(theta[at("covariance")], nrow(cells), blocks)
    matrices <- lapply(seq_len(blocks), function(b) {
      m <- matrix(0, q, q)
      m[cells] <- values[, b]
      m[cells[, 2:1, drop = FALSE]] <- values[, b]
      m
    })
    list(beta = matrix(theta[at("beta")], p_x, k),
         gamma = cbind(0, matrix(theta[at("gamma")], p_w, k - 1L)),
         covariances = matrices[if (blocks == 1L) rep(1L, k) else seq_len(k)],
         sigma2 = theta[at("sigma2")])
  }
  pack <- function(natural) {
    c(as.vector(natural$beta), as.vector(natural$gamma[, -1L]),
      as.vector(vapply(natural$covariances[seq_len(blocks)], function(m) m[cells],
                       numeric(nrow(cells)))),
      natural$sigma2)
  }
  evaluate <- function(theta) {
    p <- unpack(theta)
    linear <- setup$w %*% p$gamma
    log_prior <- linear - apply(linear, 1L, max)
    log_prior <- log_prior - log(rowSums(exp(log_prior)))
    classes <- lapply(seq_len(k), function(j) {
      .lcmm_class_density(setup, p$beta[, j], p$covariances[[j]], p$sigma2)
    })
    joint <- log_prior + vapply(classes, `[[`, numeric(setup$n), "log_density")
    joint <- matrix(joint, setup$n, k)
    peak <- apply(joint, 1L, max)
    total <- peak + log(rowSums(exp(joint - peak)))
    list(log_likelihood = sum(total), posterior = exp(joint - total),
         prior = exp(log_prior), effects = lapply(classes, `[[`, "effects"))
  }
  list(unpack = unpack, pack = pack, evaluate = evaluate,
       log_likelihood = function(theta) evaluate(theta)$log_likelihood, at = at)
}

#' The reference at hlme's mode: estimates, posteriors, effects, covariance
#' @noRd
.lcmm_reference <- function(mode, case, setup) {
  likelihood <- .lcmm_likelihood(case, setup)
  natural <- .lcmm_natural(mode$best, mode, case, setup)
  theta <- likelihood$pack(natural)
  state <- likelihood$evaluate(theta)
  # Richardson extrapolation from a first step of 1% of each coordinate (the
  # default 10% can leave a covariance block positive definite only barely).
  hessian <- numDeriv::hessian(likelihood$log_likelihood, theta,
                               method.args = list(d = 0.01, r = 4))
  vcov <- chol2inv(chol(-hessian))
  delta <- function(f) {
    jacobian <- numDeriv::jacobian(f, theta, method.args = list(d = 0.01, r = 4))
    list(estimate = f(theta),
         std_error = sqrt(diag(jacobian %*% vcov %*% t(jacobian))))
  }
  posterior_effects <- Reduce(`+`, lapply(seq_len(case$k), function(j) {
    state$effects[[j]] * state$posterior[, j]
  }))
  list(theta = theta, natural = natural, likelihood = likelihood, vcov = vcov,
       delta = delta, log_likelihood = state$log_likelihood,
       posterior = state$posterior, prior = state$prior, effects = state$effects,
       posterior_effects = posterior_effects)
}

# ---------------------------------------------------------------------------
# The expected tables, in the package's shapes
# ---------------------------------------------------------------------------

#' Wald columns from estimates and standard errors (95%)
#' @noRd
.lcmm_wald <- function(estimate, std_error) {
  z <- stats::qnorm(0.975)
  statistic <- estimate / std_error
  data.frame(estimate = estimate, std_error = std_error, statistic = statistic,
             p_value = 2 * stats::pnorm(-abs(statistic)),
             conf_low = estimate - z * std_error, conf_high = estimate + z * std_error)
}

#' A delta-method interval formed on a link scale and carried back (95%)
#' @noRd
.lcmm_linked <- function(estimate, std_error, scale) {
  z <- stats::qnorm(0.975)
  link <- switch(scale, log = log, atanh = atanh, logit = stats::qlogis)
  inverse <- switch(scale, log = exp, atanh = tanh, logit = stats::plogis)
  slope <- switch(scale, log = 1 / estimate, atanh = 1 / (1 - estimate^2),
                  logit = 1 / (estimate * (1 - estimate)))
  linked <- abs(slope) * std_error
  data.frame(estimate = estimate, std_error = std_error,
             conf_low = inverse(link(estimate) - z * linked),
             conf_high = inverse(link(estimate) + z * linked))
}

#' Every table hlme's mode supplies, built from its outputs and definitions
#'
#' `fit` is read for labels only (term and class names, the identifier's
#' name); no estimate is taken from it.
#' @noRd
.lcmm_expected <- function(reference, mode, case, setup, fit) {
  k <- case$k
  classes <- sprintf("class_%d", seq_len(k))
  natural <- reference$natural
  likelihood <- reference$likelihood
  delta <- reference$delta
  unpack <- likelihood$unpack
  tables <- list()

  beta <- delta(function(v) as.vector(unpack(v)$beta))
  tables$coefficients <- data.frame(
    class = rep(classes, each = length(setup$x_names)),
    term = rep(setup$x_names, k), .lcmm_wald(beta$estimate, beta$std_error))

  if (k > 1L) {
    gamma <- delta(function(v) as.vector(unpack(v)$gamma[, -1L]))
    tables$membership <- data.frame(
      class = rep(classes[-1L], each = length(setup$w_names)),
      term = rep(setup$w_names, k - 1L), .lcmm_wald(gamma$estimate, gamma$std_error))
  }

  # Random effects: one block when shared ("all"), else one per class.
  terms <- colnames(setup$z_names)
  q <- length(terms)
  shared <- identical(case$covariance, "equal")
  blocks <- if (shared) 1L else k
  pairs <- which(upper.tri(diag(q)), arr.ind = TRUE)
  tables$random <- do.call(rbind, lapply(seq_len(blocks), function(b) {
    covariance <- function(v) unpack(v)$covariances[[b]]
    variance <- delta(function(v) diag(covariance(v)))
    deviation <- delta(function(v) sqrt(diag(covariance(v))))
    correlation <- if (nrow(pairs) > 0L) delta(function(v) {
      g <- covariance(v)
      g[pairs] / sqrt(diag(g)[pairs[, 1L]] * diag(g)[pairs[, 2L]])
    })
    rbind(
      data.frame(class = if (shared) "all" else classes[b], parameter = "variance",
                 term = terms, .lcmm_linked(variance$estimate, variance$std_error, "log")),
      data.frame(class = if (shared) "all" else classes[b], parameter = "sd",
                 term = terms, .lcmm_linked(deviation$estimate, deviation$std_error, "log")),
      if (!is.null(correlation)) data.frame(
        class = if (shared) "all" else classes[b], parameter = "correlation",
        term = paste(terms[pairs[, 1L]], terms[pairs[, 2L]], sep = ":"),
        .lcmm_linked(correlation$estimate, correlation$std_error, "atanh")))
  }))

  # hlme reports no posteriors for one class: every person is in it.
  posterior <- if (k > 1L) as.matrix(mode$pprob[sprintf("prob%d", seq_len(k))]) else
    matrix(1, setup$n, 1L)
  modal <- if (k > 1L) mode$pprob$class else rep(1L, setup$n)
  share <- delta(function(v) {
    linear <- setup$w %*% unpack(v)$gamma
    colMeans(exp(linear) / rowSums(exp(linear)))
  })
  shares <- .lcmm_linked(share$estimate, share$std_error, "logit")
  residual <- delta(function(v) sqrt(unpack(v)$sigma2))
  tables$classes <- data.frame(
    class = classes, share = shares$estimate, share_std_error = shares$std_error,
    share_conf_low = shares$conf_low, share_conf_high = shares$conf_high,
    expected_persons = colSums(posterior), assigned_persons = tabulate(modal, k),
    average_posterior = vapply(seq_len(k), function(j) {
      if (any(modal == j)) mean(posterior[modal == j, j]) else NA_real_
    }, numeric(1)),
    residual_sd = abs(mode$best[["stderr"]]),
    residual_sd_std_error = residual$std_error)

  n_persons <- setup$n
  q_count <- length(mode$best)
  log_likelihood <- .lcmm_log_likelihood(mode, setup$n_rows)
  deviance <- -2 * log_likelihood
  entropy_sum <- -sum(ifelse(posterior > 0, posterior * log(posterior), 0))
  tables$fit <- data.frame(
    n_persons = n_persons, n_observations = setup$n_rows,
    random_diagonal = isTRUE(mode$idiag == 1),
    n_parameters = q_count, log_likelihood = log_likelihood,
    aic = deviance + 2 * q_count, bic = deviance + log(n_persons) * q_count,
    bic_rows = deviance + log(setup$n_rows) * q_count,
    sabic = deviance + log((n_persons + 2) / 24) * q_count,
    icl = deviance + log(n_persons) * q_count + 2 * entropy_sum,
    entropy = if (k > 1L) 1 - entropy_sum / (n_persons * log(k)) else NA_real_,
    smallest_share = min(colSums(posterior) / n_persons))

  id <- fit$spec$id
  assignments <- data.frame(.row = seq_len(n_persons), probability = apply(posterior, 1L, max))
  assignments[sprintf("posterior_class_%d", seq_len(k))] <- as.data.frame(posterior)
  tables$assignments <- assignments

  effects <- t(vapply(seq_len(n_persons), function(i) {
    reference$effects[[modal[[i]]]][i, ]
  }, numeric(q))) |> matrix(n_persons, q)
  tables$random_effects <- stats::setNames(
    data.frame(setup$persons, classes[modal], effects, stringsAsFactors = FALSE),
    c(id, "class", terms))

  times <- setup$frame$time
  grid <- seq(min(times), max(times), length.out = 60L)
  # predictY() asks for every covariate of the model, the membership one
  # included; the class-specific mean trajectories do not depend on it.
  newdata <- data.frame(time = grid, x = 0)
  marginal <- lcmm::predictY(mode, newdata, var.time = "time")$pred
  design <- cbind(1, grid)
  tables$trajectories <- do.call(rbind, lapply(seq_len(k), function(j) {
    covariance <- reference$vcov[likelihood$at("beta"), likelihood$at("beta")]
    block <- (j - 1L) * length(setup$x_names) + seq_along(setup$x_names)
    std_error <- sqrt(rowSums((design %*% covariance[block, block]) * design))
    out <- data.frame(class = classes[j], time = grid,
                      .lcmm_wald(as.vector(marginal[, j]), std_error)[
                        c("estimate", "std_error", "conf_low", "conf_high")])
    out
  }))

  predictions <- mode$pred
  person <- match(predictions$idn, unique(setup$frame$idn))
  row_class <- modal[person]
  # With one class hlme's predictions carry no class suffix.
  pick <- function(prefix) {
    columns <- as.matrix(predictions[if (k > 1L) sprintf("%s%d", prefix, seq_len(k)) else prefix])
    columns[cbind(seq_len(nrow(columns)), row_class)]
  }
  tables$individual <- stats::setNames(
    data.frame(setup$persons[person], setup$frame$time, classes[row_class],
               predictions$obs, pick("pred_m"),
               pick("pred_ss"), stringsAsFactors = FALSE),
    c(id, "time", "class", "observed", "class_mean", "predicted"))
  tables
}

#' The output columns that are free parameters of the model, with their errors
#'
#' Everything else a growth fit reports is a function of these and the data.
#' @noRd
.lcmm_parameters <- list(
  coefficients = c(estimate = "std_error"), membership = c(estimate = "std_error"),
  random = c(estimate = "std_error"),
  classes = c(share = "share_std_error", residual_sd = "residual_sd_std_error"))

#' The free parameters against hlme's fixed point, in standard errors
#'
#' Each estimate's distance from hlme's fixed point divided by the estimate's
#' own standard error (the fit's `std_error` column, observed information),
#' held to 1e-3: within a thousandth of a standard error no estimate,
#' interval or test computed from the fit can differ, and the package stops
#' on its own rule (a quasi-Newton finish), as every reference program stops
#' on its own. Machine-precision agreement of the computation is established
#' separately, by comparing every function of the parameters at the package's
#' own point. A cell with no standard error (none exists, or it is zero: a
#' one-class share) keeps an absolute comparison with the fixed point at the
#' harness default, 1e-6, in its own block under the same output key, so the
#' ledger counts every value once.
#'
#' @param obtained The fit's tables, from `fit_tables()`.
#' @param at_mode_tables Tables built from hlme's fixed point.
#' @return Comparison rows with `output` set.
#' @noRd
.lcmm_parameter_distance <- function(label, obtained, at_mode_tables) {
  rows <- lapply(names(.lcmm_parameters), function(table) {
    have <- obtained[[table]]
    want <- at_mode_tables[[table]]
    if (is.null(want) || nrow(have) == 0L) return(NULL)
    keys <- .output_keys(want, table)
    as_key <- function(frame) do.call(paste, c(lapply(frame[keys], as.character), sep = "\r"))
    at <- match(as_key(have), as_key(want))
    stopifnot("every parameter cell must be in the fixed-point reference" = !anyNA(at))
    columns <- .lcmm_parameters[[table]]
    lapply(names(columns), function(column) {
      estimate <- have[[column]]
      reference <- want[[column]][at]
      error <- have[[columns[[column]]]]
      with_error <- is.finite(error) & error > 0
      output <- paste0(table, "$", column)
      rbind(
        if (any(with_error)) compare_block(
          sprintf("%s: %s, distance from hlme's fixed point in standard errors",
                  label, output),
          rep(0, sum(with_error)),
          abs(estimate[with_error] - reference[with_error]) / error[with_error],
          tolerance = 1e-3, output = output),
        if (any(!with_error)) compare_block(
          sprintf("%s: %s, no standard error, against hlme's fixed point", label, output),
          reference[!with_error], estimate[!with_error],
          tolerance = .output_default_tolerances[["default"]], output = output))
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

# ---------------------------------------------------------------------------
# hlme's own standard errors
# ---------------------------------------------------------------------------

#' hlme's `V` against the reference observed information
#'
#' hlme's `V` is a covariance of its estimates in hlme's own coordinates: the
#' Cholesky factor of the random-effect covariance (its `cholesky`) in place
#' of `varcov`, logits against the last class, the residual SD. The errors of
#' the package's natural parameters follow by the delta method (a Jacobian
#' taken here numerically), and are compared, as ratios, with the reference
#' errors, outside the ledger.
#'
#' What hlme reports, measured on the six cases (2026-10-26, lcmm 2.2.1): its
#' errors are within 2.2e-3 of the inverse observed information (the
#' reference's Richardson-extrapolated Hessian) on every parameter, and within
#' 4.5e-4 on all but two (membership-covariate: 2.2e-3 on class 1's slope,
#' 1.4e-3 on class 3's membership intercept); they are 0.12 to 0.17 (largest
#' relative gap per case)
#' from the outer-product-of-scores (OPG) errors, whose own gap to the observed
#' ones is 0.13 to 0.21. hlme's `V` is therefore the inverse observed
#' information, computed by numerical second differences less accurate than
#' Richardson extrapolation: the tolerance, 3e-3 relative, is the largest gap
#' measured with a margin, and is forty times below what separates the
#' observed information from the OPG, so it still identifies the estimator.
#' The package's own errors (`vcov_type = "observed"`) are compared with the
#' reference at 1e-4 in the ledger.
#' @noRd
.lcmm_own_error_rows <- function(label, mode, case, setup, reference) {
  n_best <- length(mode$best)
  covariance <- matrix(0, n_best, n_best)
  covariance[upper.tri(covariance, diag = TRUE)] <- mode$V
  covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
  counts <- mode$N
  random_at <- sum(counts[seq_len(2L)]) + seq_len(counts[[3L]])
  q <- ncol(setup$z_names)
  internal <- mode$best
  internal[random_at] <- mode$cholesky
  # hlme's internal coordinates to the natural ones.
  to_natural <- function(values) {
    factor <- matrix(0, q, q)
    factor[upper.tri(factor, diag = TRUE)] <- values[random_at]
    covariance_values <- crossprod(factor)
    values[random_at] <- covariance_values[upper.tri(covariance_values, diag = TRUE)]
    .lcmm_natural(values, mode, case, setup)
  }
  wanted <- function(values) {
    natural <- to_natural(values)
    c(as.vector(natural$beta), as.vector(natural$gamma[, -1L]),
      sqrt(natural$sigma2),
      unlist(lapply(natural$covariances[if (identical(case$covariance, "equal")) 1L else
        seq_len(case$k)], diag)))
  }
  jacobian <- numDeriv::jacobian(wanted, internal, method.args = list(d = 0.01, r = 4))
  theirs <- sqrt(diag(jacobian %*% covariance %*% t(jacobian)))
  unpack <- reference$likelihood$unpack
  ours <- reference$delta(function(v) {
    p <- unpack(v)
    c(as.vector(p$beta), as.vector(p$gamma[, -1L]), sqrt(p$sigma2),
      unlist(lapply(p$covariances[if (identical(case$covariance, "equal")) 1L else
        seq_len(case$k)], diag)))
  })$std_error
  groups <- rep(c("fixed effects", "membership logits", "residual SD",
                  "random-effect variances"),
                c(length(reference$natural$beta), length(reference$natural$gamma[, -1L]),
                  1L, length(theirs) - length(reference$natural$beta) -
                    length(reference$natural$gamma[, -1L]) - 1L))
  do.call(rbind, lapply(unique(groups), function(group) {
    compare_block(sprintf(paste("%s: hlme's own V-based standard errors, %s, as a ratio to",
                                "the reference observed-information ones"), label, group),
                  rep(1, sum(groups == group)), theirs[groups == group] / ours[groups == group],
                  tolerance = 3e-3)
  }))
}

# ---------------------------------------------------------------------------
# One class: nlme and lme4 by maximum likelihood
# ---------------------------------------------------------------------------

#' The one-class fit against nlme::lme and lme4::lmer (ML)
#' @return Comparison rows, outside the ledger (the ledger is hlme's).
#' @noRd
.lcmm_mixed_model_rows <- function(label, frame, fit) {
  tables <- fit_tables(fit)
  coefficients <- tables$coefficients$estimate
  random <- tables$random
  variances <- random$estimate[random$parameter == "variance"]
  correlation <- random$estimate[random$parameter == "correlation"]
  effects <- as.matrix(tables$random_effects[c("(Intercept)", "time")])
  lme <- nlme::lme(y ~ time, random = ~ time | id, data = frame, method = "ML",
                   control = nlme::lmeControl(maxIter = 500, msMaxIter = 500,
                                              tolerance = 1e-12, msTol = 1e-14,
                                              niterEM = 100))
  lme_cov <- nlme::getVarCov(lme)
  lmer <- lme4::lmer(y ~ time + (time | id), data = frame, REML = FALSE,
                     control = lme4::lmerControl(optimizer = "nloptwrap",
                                                 optCtrl = list(xtol_abs = 1e-14,
                                                                ftol_abs = 1e-14,
                                                                xtol_rel = 0, ftol_rel = 0,
                                                                maxeval = 1e5)))
  lmer_cov <- lme4::VarCorr(lmer)$id
  rows <- function(program, log_likelihood, fixed, covariance, sigma, blups) {
    rbind(
      compare_values(sprintf("%s: %s log likelihood (ML)", label, program),
                     log_likelihood, fit$log_likelihood, tolerance = 1e-8, scale = "relative"),
      compare_block(sprintf("%s: %s fixed effects", label, program), unname(fixed),
                    coefficients, tolerance = 1e-6),
      compare_block(sprintf("%s: %s random-effect variances and correlation", label, program),
                    c(diag(covariance), stats::cov2cor(covariance)[1L, 2L]),
                    c(variances, correlation), tolerance = 1e-6),
      compare_values(sprintf("%s: %s residual SD", label, program), sigma,
                     tables$classes$residual_sd, tolerance = 1e-6),
      compare_block(sprintf("%s: %s predicted random effects (BLUPs)", label, program),
                    as.vector(blups), as.vector(effects), tolerance = 1e-6))
  }
  rbind(
    rows("nlme::lme", as.numeric(stats::logLik(lme)), nlme::fixef(lme),
         unname(matrix(lme_cov, 2L)), lme$sigma,
         as.matrix(nlme::ranef(lme)[tables$random_effects$id, ])),
    rows("lme4::lmer", as.numeric(stats::logLik(lmer)), lme4::fixef(lmer),
         unname(matrix(lmer_cov, 2L)), stats::sigma(lmer),
         as.matrix(lme4::ranef(lmer)$id[tables$random_effects$id, ])))
}

# ---------------------------------------------------------------------------
# Equivariance
# ---------------------------------------------------------------------------

#' The fit of an outcome in units of `scale` against the unit-scale fit
#'
#' Coefficients and their errors scale by `scale`, variances by its square,
#' standard deviations by it; correlations, shares, logits and posteriors do
#' not move; the log likelihood falls by (rows) log(scale).
#' @return Comparison rows, outside the ledger.
#' @noRd
.lcmm_equivariance <- function(unit, scaled, scale) {
  a <- fit_tables(unit)
  b <- fit_tables(scaled)
  label <- sprintf("scale-1e4 against the unit-scale fit (outcome x %g)", scale)
  power <- c(variance = 2, sd = 1, correlation = 0)[a$random$parameter]
  rbind(
    compare_values(sprintf("%s: log likelihood + rows x log(scale)", label),
                   unit$log_likelihood,
                   scaled$log_likelihood + a$fit$n_observations * log(scale),
                   tolerance = 1e-8, scale = "relative"),
    compare_block(sprintf("%s: coefficients / scale", label),
                  a$coefficients$estimate, b$coefficients$estimate / scale, tolerance = 1e-6),
    compare_block(sprintf("%s: coefficient standard errors / scale", label),
                  a$coefficients$std_error, b$coefficients$std_error / scale,
                  tolerance = 1e-4),
    compare_block(sprintf("%s: random-effect estimates / scale^power", label),
                  a$random$estimate, b$random$estimate / scale^power, tolerance = 1e-6),
    compare_block(sprintf("%s: random-effect standard errors / scale^power", label),
                  a$random$std_error, b$random$std_error / scale^power, tolerance = 1e-4),
    compare_block(sprintf("%s: residual SD / scale", label),
                  a$classes$residual_sd, b$classes$residual_sd / scale, tolerance = 1e-6),
    compare_block(sprintf("%s: class shares", label), a$classes$share, b$classes$share,
                  tolerance = 1e-6),
    compare_block(sprintf("%s: membership logits", label), a$membership$estimate,
                  b$membership$estimate, tolerance = 1e-6),
    compare_block(sprintf("%s: posteriors", label),
                  unlist(a$assignments[grep("^posterior_", names(a$assignments))]),
                  unlist(b$assignments[grep("^posterior_", names(b$assignments))]),
                  tolerance = 1e-6),
    compare_block(sprintf("%s: predicted random effects / scale", label),
                  unlist(a$random_effects[c("(Intercept)", "time")]),
                  unlist(b$random_effects[c("(Intercept)", "time")]) / scale,
                  tolerance = 1e-6))
}
