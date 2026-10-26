# mixture_regression() against flexmix (Leisch 2004, Journal of Statistical
# Software 11(8); Gruen & Leisch 2008, Journal of Statistical Software 28(4)),
# on every output column of the fit.
#
# `equivalence/test-mixture-regression-flexmix.R` compares chosen quantities
# (the likelihood at flexmix's estimates, the coefficients, the parameter
# count, flexmix's BIC). This suite compares every value column of
# `get_results(fit, "all")`, or exempts it with its reason, for seven
# simulated fixtures: Gaussian with two and with three classes and a
# concomitant model (`membership`), binomial counts, Poisson with a
# concomitant model, Gaussian with a coefficient shared by the classes
# (`common`), row-level classes with cluster-robust standard errors (`id`,
# `vcov_type = "robust"`), and one class per group (`class_level = "group"`,
# flexmix's grouped formula `y ~ x | id`) with a group-level concomitant.
#
# Families flexmix cannot fit as the package does are left out: flexmix
# (2.3.20) has no cumulative-logit (ordinal) driver and no negative-binomial
# one -- FLXMRglm's families are gaussian, binomial, poisson and Gamma, and
# FLXMRmultinom is the unordered multinomial logit. The package's own tests
# compare those families with MASS::polr() and MASS::glm.nb().
#
# Every case makes two claims, as the hard depmixS4 suite does:
#   (A) computation at one mode. mixture_regression() takes no starting
#       values, so flexmix is started from the package's solution: its EM is
#       run from the package's posterior matrix (flexmix's `cluster =`) to its
#       fixed point ([.flexmix_fixed_point()]); the log-likelihood that gains
#       says whether the package stopped at the mode. The free parameters are
#       compared with that fixed point, and every function of them
#       (posteriors, assignments, fitted values, counts, classification,
#       criteria, standard errors, Wald columns) with flexmix evaluated at the
#       package's own estimates ([.flexmix_at_estimates()]), so that neither
#       program's stopping rule is charged to the other's computation
#       ([.flexmix_same_mode()]). Estimates and their functions to 1e-6,
#       standard errors and the columns computed from them to 1e-4 (below).
#   (B) search. The package's own random starts against the best of ten
#       flexmix searches from flexmix's own random initial partitions, the
#       best refined to its fixed point. A shortfall of the package's maximum
#       is a finding about its search; a higher package maximum is recorded as
#       a row that is not an output of the fit.
#
# Findings (2026-10-26, flexmix 2.3.20, package 0.97.1): no disagreement.
#   * Free parameters against flexmix's fixed point: within 3.8e-7 (a
#     concomitant logit on gaussian-3-membership), class coefficients within
#     1.9e-7, residual SDs 1.4e-7; flexmix's EM gains at most 7.5e-12 in
#     log-likelihood from the package's solution.
#   * At the package's estimates flexmix's log-likelihood, posteriors, fitted
#     values, counts, classification and criteria agree to rounding (8.9e-16
#     in a posterior, 3.6e-15 in a residual); standard errors to 8.5e-9, and
#     flexmix's refit(method = "optim") errors to the reference to 1.0e-7.
#   * Search: on all seven fixtures the ten flexmix searches reach one
#     maximum and the package's own search reaches it (shortfall 8.2e-12 or
#     less). These are well-separated fixtures; they test computation, not
#     search on hard surfaces.
#
# Two flexmix conventions are not maximum likelihood, and the reference
# replaces them by the maximum-likelihood step inside flexmix's own EM:
#   * Gaussian variance. flexmix's FLXMRglm M-step sets
#     sigma = sqrt(sum(w r^2) / mean(w) / (nrow(x) - rank)), a degrees-of-freedom
#     divisor, so its fixed point is not the maximum of the likelihood it
#     reports (the existing flexmix test records this). The reference driver
#     keeps flexmix's weighted least squares and sets sigma^2 = sum(w r^2) /
#     sum(w), the ML step ([.flexmix_gaussian_ml()]). For the shared-coefficient
#     case the same holds for FLXMRglmfix(varFix = TRUE).
#   * Concomitant M-step precision. FLXPmultinom solves the weighted
#     multinomial logit with nnet's BFGS from zero at nnet's defaults
#     (reltol 1e-8, maxit 100). Its EM, run to the end from the package's
#     solution, then stops 5.1e-5 from the tightly solved fixed point in a
#     concomitant logit and 4.3e-6 in a class coefficient (gaussian-2- and
#     -3-membership; measured 2026-10-26). The reference passes reltol 1e-15
#     and maxit 10000 to the same nnet call ([.flexmix_concomitant()]). The
#     estimator is unchanged; only its solver tolerance is.
#
# The shared-coefficient case uses equal variances. With class-specific
# variances flexmix's FLXMRglmfix stacks the classes into one regression
# weighted by the posteriors alone, not by the posteriors over each class's
# variance, so its M-step for the shared coefficient is not the ML step and
# its fixed point is not the package's.
#
# Standard errors. flexmix offers two estimators through `refit()`:
#   * `method = "mstep"`: the weighted GLM standard errors of one M-step,
#     which treat the posteriors as known (the complete-data information).
#     On these fixtures they are 5 to 38 percent below the observed-information
#     errors of the class coefficients (0.9 percent below to 0.2 above on
#     gaussian-group-classes, whose posteriors are nearly 0 or 1).
#   * `method = "optim"` (the default): BFGS from flexmix's estimates and the
#     inverse of optim's finite-difference Hessian of flexmix's analytic
#     gradient -- the inverse observed information, in flexmix's coordinates
#     (log sigma, concomitant logits). It agrees with the reference below to
#     1.2e-6 relative or better on every fixture, FLXMRglmfix included
#     (measured 2026-10-26).
# The package reports the inverse observed information (`vcov_type =
# "observed"`) or its cluster sandwich with the CR1 factor G / (G - 1)
# (`"robust"`), for which flexmix has no counterpart. Like is compared with
# like through a reference written here: the log-likelihood from scratch in
# base R in the natural parameters (class coefficients, shared coefficients,
# residual SDs, concomitant logits), its Hessian by Richardson extrapolation
# (numDeriv), and for the robust case the unit scores by numDeriv's Jacobian,
# summed within clusters. At a maximum the inverse information transforms with
# the Jacobian of any reparameterisation, so these are the standard errors of
# the reported quantities whatever scale the package works on; the class
# priors' errors come from the logits by the softmax Jacobian. flexmix's
# `refit(method = "optim")` errors of the class and concomitant coefficients
# are checked against the reference's observed-information errors in every
# case, outside the ledger, so the reference is shown to be flexmix's own
# estimator. The package differentiates its analytic score numerically
# (central differences, step a thousandth of a standard error), and the
# reference is a second difference of the likelihood; both carry truncation
# error, so standard errors and every column computed from them (statistic,
# p value, interval, adjusted p value) are held to 1e-4, the harness's
# standard-error tolerance.

suite_flexmix <- function() {
  require_suite_packages(
    "flexmix", "nnet", "numDeriv",
    reason = paste("the regression-mixture reference, its concomitant M-step solver,",
                   "and the reference observed information"))
  stack_cases(lapply(.flexmix_fixtures(), .flexmix_case))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' Run code under a seed and restore the caller's random stream
#' @param seed Integer seed.
#' @param code Code to evaluate.
#' @return The value of `code`.
#' @noRd
.flexmix_with_seed <- function(seed, code) {
  saved <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(saved)) {
    if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", saved, globalenv()), add = TRUE)
  set.seed(seed)
  code
}

#' The simulated fixtures
#'
#' Each fixture names the package call (`formula`, `membership`, `common`,
#' `id`, `class_level`, `variance`, `vcov_type`), the flexmix formula, the
#' design of the from-scratch likelihood (`terms`, the class-specific
#' predictors including the intercept), and a simulator.
#'
#' @return A list of fixtures.
#' @noRd
.flexmix_fixtures <- function() {
  base <- list(family = "gaussian", membership = NULL, common = NULL, id = NULL,
               class_level = "observation", variance = "varying",
               vcov_type = "observed", terms = ~ x)
  fixture <- function(...) utils::modifyList(base, list(...))
  # Class sizes and lines well apart, so every class holds a fifth of the
  # rows or more and flexmix's searches find the same maximum as the package.
  rows <- function(n, seed, classes) {
    .flexmix_with_seed(seed, {
      frame <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n), r = stats::rnorm(n))
      frame$truth <- classes(frame)
      frame
    })
  }
  two_by_z <- function(frame) 1L + (stats::runif(nrow(frame)) < stats::plogis(0.3 + 0.9 * frame$z))
  list(
    fixture(name = "gaussian-2-membership", seed = 501L, n_classes = 2L,
            formula = y ~ x, flexmix_formula = y ~ x, membership = ~ z,
            simulate = function() {
              frame <- rows(400L, 501L, two_by_z)
              .flexmix_with_seed(5011L, {
                frame$y <- ifelse(frame$truth == 1L, 1 + 2 * frame$x, -1 - frame$x) +
                  stats::rnorm(nrow(frame), sd = ifelse(frame$truth == 1L, 0.7, 1))
              })
              frame
            }),
    fixture(name = "gaussian-3-membership", seed = 502L, n_classes = 3L,
            formula = y ~ x, flexmix_formula = y ~ x, membership = ~ z,
            simulate = function() {
              frame <- rows(600L, 502L, function(frame) {
                logits <- cbind(0, 0.2 + 0.8 * frame$z, -0.2 - 0.8 * frame$z)
                p <- exp(logits) / rowSums(exp(logits))
                vapply(seq_len(nrow(p)), function(i) sample.int(3L, 1L, prob = p[i, ]),
                       integer(1))
              })
              .flexmix_with_seed(5021L, {
                intercept <- c(3, 0, -3)[frame$truth]
                slope <- c(1.5, -1, 0.5)[frame$truth]
                frame$y <- intercept + slope * frame$x +
                  stats::rnorm(nrow(frame), sd = c(0.8, 0.6, 1)[frame$truth])
              })
              frame
            }),
    # Binomial counts with 4 to 10 trials a row: one Bernoulli trial per class
    # assignment is refused by the package (not identified), several are not.
    fixture(name = "binomial-counts", seed = 503L, n_classes = 2L, family = "binomial",
            formula = cbind(s, f) ~ x, flexmix_formula = cbind(s, f) ~ x,
            simulate = function() {
              frame <- rows(400L, 503L, function(frame) 1L + (stats::runif(nrow(frame)) < 0.4))
              .flexmix_with_seed(5031L, {
                frame$trials <- sample(4:10, nrow(frame), replace = TRUE)
                frame$s <- stats::rbinom(nrow(frame), frame$trials, stats::plogis(
                  ifelse(frame$truth == 1L, -1 + 1.2 * frame$x, 1 - frame$x)))
              })
              frame$f <- frame$trials - frame$s
              frame
            }),
    fixture(name = "poisson-membership", seed = 504L, n_classes = 2L, family = "poisson",
            formula = y ~ x, flexmix_formula = y ~ x, membership = ~ z,
            simulate = function() {
              frame <- rows(400L, 504L, two_by_z)
              .flexmix_with_seed(5041L, {
                frame$y <- stats::rpois(nrow(frame), exp(ifelse(frame$truth == 1L,
                                                                0.5 + 0.8 * frame$x,
                                                                2 - 0.4 * frame$x)))
              })
              frame
            }),
    # `r` has one slope in both classes (flexmix's FLXMRglmfix(fixed = ~ r)).
    fixture(name = "gaussian-common-equal", seed = 505L, n_classes = 2L,
            formula = y ~ x + r, flexmix_formula = y ~ x, common = ~ r,
            variance = "equal",
            simulate = function() {
              frame <- rows(400L, 505L, function(frame) 1L + (stats::runif(nrow(frame)) < 0.4))
              .flexmix_with_seed(5051L, {
                frame$y <- ifelse(frame$truth == 1L, 1 + 2 * frame$x, -1.5 - frame$x) +
                  0.6 * frame$r + stats::rnorm(nrow(frame), sd = 0.8)
              })
              frame
            }),
    # Rows in 80 clusters of five, each cluster sharing a random intercept, so
    # the cluster sandwich differs from the inverse information.
    fixture(name = "gaussian-cluster-robust", seed = 506L, n_classes = 2L,
            formula = y ~ x, flexmix_formula = y ~ x, id = "id",
            vcov_type = "robust",
            simulate = function() {
              frame <- rows(400L, 506L, function(frame) 1L + (stats::runif(nrow(frame)) < 0.45))
              frame$id <- sprintf("c%03d", rep(seq_len(80L), each = 5L))
              .flexmix_with_seed(5061L, {
                shared <- stats::rnorm(80L, sd = 0.6)[rep(seq_len(80L), each = 5L)]
                frame$y <- ifelse(frame$truth == 1L, 1 + 2 * frame$x, -1 - frame$x) +
                  shared + stats::rnorm(nrow(frame), sd = 0.8)
              })
              frame
            }),
    # One class per group of six rows; the class depends on a group-level `w`.
    fixture(name = "gaussian-group-classes", seed = 507L, n_classes = 2L,
            formula = y ~ x, flexmix_formula = y ~ x | id, membership = ~ w,
            id = "id", class_level = "group",
            simulate = function() {
              .flexmix_with_seed(507L, {
                groups <- 90L
                size <- 6L
                w <- stats::rnorm(groups)
                truth <- 1L + (stats::runif(groups) < stats::plogis(-0.2 + w))
                index <- rep(seq_len(groups), each = size)
                frame <- data.frame(id = sprintf("g%03d", index), w = w[index],
                                    x = stats::rnorm(groups * size), truth = truth[index])
                frame$y <- ifelse(frame$truth == 1L, 1 + frame$x, -0.5 - 1.5 * frame$x) +
                  stats::rnorm(nrow(frame))
                frame
              })
            })
  )
}

# ---------------------------------------------------------------------------
# flexmix, with the maximum-likelihood steps
# ---------------------------------------------------------------------------

#' flexmix's Gaussian driver with the maximum-likelihood variance step
#'
#' flexmix's own weighted least squares (`lm.wfit`) for the coefficients;
#' only the divisor of the residual variance changes, from flexmix's
#' `mean(w) (nrow(x) - rank)` to `sum(w)`.
#'
#' @param common `NULL`, or the one-sided formula of the shared terms, which
#'   selects FLXMRglmfix with one variance shared by the classes.
#' @return A flexmix model driver.
#' @noRd
.flexmix_gaussian_ml <- function(common) {
  if (is.null(common)) {
    model <- flexmix::FLXMRglm(family = "gaussian")
    model@fit <- function(x, y, w, component) {
      fit <- stats::lm.wfit(x, y, w = w)
      model@defineComponent(para = list(
        coef = stats::coef(fit), df = ncol(x) + 1,
        sigma = sqrt(sum(fit$weights * fit$residuals^2) / sum(fit$weights))))
    }
    return(model)
  }
  model <- flexmix::FLXMRglmfix(family = "gaussian", fixed = common, varFix = TRUE)
  # flexmix stacks the classes into one regression whose weights are the
  # posteriors; with one shared variance that is the ML step for every
  # coefficient, and sum(w) over the stack is the number of rows.
  model@fit <- function(x, y, w, incidence, variance, ...) {
    stopifnot("the ML driver covers one variance shared by every class" =
                length(variance) == 1L && variance == nrow(incidence))
    fit <- stats::lm.wfit(x, y, w = w)
    coefs <- stats::coef(fit)
    names(coefs) <- colnames(incidence)
    sigma <- sqrt(sum(fit$weights * fit$residuals^2) / sum(fit$weights))
    # flexmix's own degrees-of-freedom bookkeeping, unchanged.
    df <- rowSums(incidence / rep(colSums(incidence), each = nrow(incidence))) +
      rep(1 / variance, variance)
    lapply(seq_len(nrow(incidence)), function(class) {
      model@defineComponent(list(coef = coefs[as.logical(incidence[class, ])],
                                 sigma = sigma, df = df[[class]]))
    })
  }
  model
}

#' The flexmix model driver for a fixture
#' @param fixture One fixture.
#' @return A flexmix model driver.
#' @noRd
.flexmix_driver <- function(fixture) {
  if (identical(fixture$family, "gaussian")) return(.flexmix_gaussian_ml(fixture$common))
  stopifnot("shared coefficients are covered for the Gaussian family only" =
              is.null(fixture$common))
  flexmix::FLXMRglm(family = fixture$family)
}

#' flexmix's multinomial concomitant model, solved to precision
#'
#' The M-step (`@fit`) and the final coefficient fit (`@refit`) are
#' FLXPmultinom's own nnet calls, with nnet's relative tolerance at 1e-15
#' (from its default 1e-8), its iteration bound at 10000 (from 100) and its
#' absolute tolerance at 0 (from 1e-4; the criterion is a cross-entropy and
#' never reaches it). Without a concomitant formula flexmix's constant prior
#' (FLXPconstant) is used, whose M-step is the column means of the posteriors.
#'
#' @param membership `NULL` or a one-sided formula.
#' @return A flexmix concomitant model.
#' @noRd
.flexmix_concomitant <- function(membership) {
  if (is.null(membership)) return(flexmix::FLXPconstant())
  concomitant <- flexmix::FLXPmultinom(membership)
  multinom_fit <- environment(concomitant@fit)$multinom.fit
  concomitant@fit <- function(x, y, w, ...) {
    multinom_fit(x, y, w, maxit = 10000L, reltol = 1e-15, abstol = 0)$fitted.values
  }
  concomitant@refit <- function(x, y, w, ...) {
    if (missing(w) || is.null(w)) w <- rep(1, nrow(y))
    rownames(y) <- rownames(x) <- NULL
    fit <- nnet::multinom(y ~ 0 + x, weights = w, data = list(y = y, x = x),
                          trace = FALSE, maxit = 10000L, reltol = 1e-15, abstol = 0)
    fit$coefnames <- colnames(x)
    fit$vcoefnames <- fit$coefnames[seq_along(fit$coefnames)]
    fit
  }
  concomitant
}

#' Run flexmix's EM from a posterior matrix
#'
#' `tolerance = 0` makes flexmix run exactly `iterations` EM steps;
#' `minprior = 0` keeps every class (flexmix's default drops one whose prior
#' falls under 0.05).
#'
#' @param fixture One fixture.
#' @param data The data.
#' @param posterior Rows-by-classes starting posteriors, or `NULL` for
#'   flexmix's random initial partition.
#' @param iterations Number of EM steps.
#' @param tolerance flexmix's relative log-likelihood tolerance.
#' @return A `flexmix` object.
#' @noRd
.flexmix_run <- function(fixture, data, posterior, iterations, tolerance = 0) {
  flexmix::flexmix(fixture$flexmix_formula, data = data, k = fixture$n_classes,
                   cluster = posterior, model = .flexmix_driver(fixture),
                   concomitant = .flexmix_concomitant(fixture$membership),
                   control = list(iter.max = iterations, tolerance = tolerance,
                                  minprior = 0, verbose = 0))
}

#' Every free parameter of a flexmix fit, as one vector
#' @param fitted A `flexmix` object.
#' @return A numeric vector.
#' @noRd
.flexmix_parameter_vector <- function(fitted) {
  concomitant <- flexmix::parameters(fitted, which = "concomitant")
  c(as.vector(flexmix::parameters(fitted)), as.vector(as.matrix(concomitant)),
    flexmix::prior(fitted))
}

#' flexmix's fixed point from a starting posterior
#'
#' flexmix's EM is run in blocks of 200 steps until one further step moves no
#' parameter by more than 1e-9, or for at most five blocks. Restarting from the
#' posterior continues the EM exactly: flexmix's state between steps is the
#' posterior matrix, from which each step begins with its M-step.
#'
#' With a concomitant model the step never falls to 1e-9: nnet's BFGS stops on
#' a change in its criterion, which locates the minimiser only to about
#' sqrt(eps |f| / curvature), some 4e-8 for a cross-entropy of a few hundred
#' over a few hundred rows. On gaussian-2-membership successive blocks gave
#' steps of 1.7e-9, 6.5e-8, 2.1e-8, 1.5e-8, 4.8e-9 and 5.4e-8, not decreasing
#' (measured 2026-10-26): the noise floor of the concomitant M-step, not a
#' slow approach. On gaussian-3-membership the last step was 2.2e-7, on the
#' other concomitant fixtures 9e-9 to 3.5e-8, and 4e-16 where there is no
#' concomitant model and flexmix's M-steps are exact (2.8e-9 on
#' binomial-counts, whose EM had not finished contracting after 1000 steps;
#' measured 2026-10-26). The last step is returned and held in claim A to the
#' estimate tolerance, 1e-6: the reference is then at a fixed point to within
#' what the estimates are compared to.
#'
#' @param fixture One fixture.
#' @param data The data.
#' @param posterior Starting posteriors.
#' @return A list with the `fitted` flexmix object at the fixed point and
#'   `step`, the largest change one more EM step makes.
#' @noRd
.flexmix_fixed_point <- function(fixture, data, posterior) {
  state <- Reduce(function(state, block) {
    if (state$settled) return(state)
    fitted <- .flexmix_run(fixture, data, state$posterior, 200L)
    following <- .flexmix_run(fixture, data, flexmix::posterior(fitted), 1L)
    step <- max(abs(.flexmix_parameter_vector(following) -
                      .flexmix_parameter_vector(fitted)))
    list(fitted = fitted, posterior = flexmix::posterior(fitted), step = step,
         settled = step <= 1e-9)
  }, seq_len(5L), init = list(fitted = NULL, posterior = posterior, step = Inf,
                              settled = FALSE))
  list(fitted = state$fitted, step = state$step)
}

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

#' Fit the package to a fixture
#'
#' `tol = 1e-14`: the package stops on a relative log-likelihood change, and at
#' 1e-12 it stopped 1.3e-6 from the maximum in a concomitant logit on
#' gaussian-2-membership (log-likelihood 6.6e-11 below it), the size of the
#' estimate tolerance; at 1e-15 it was 4.9e-8 from it (Newton's maximum of the
#' from-scratch likelihood; measured 2026-10-26).
#'
#' @param fixture One fixture.
#' @param data The data.
#' @return A `latents_mixture_regression` fit.
#' @noRd
.flexmix_package_fit <- function(fixture, data) {
  arguments <- list(formula = fixture$formula, data = data, n_classes = fixture$n_classes,
                    family = fixture$family, id = fixture$id,
                    class_level = fixture$class_level, variance = fixture$variance,
                    vcov_type = fixture$vcov_type, n_starts = 20L, seed = fixture$seed,
                    tol = 1e-14, max_iter = 20000L)
  if (!is.null(fixture$membership)) arguments$membership <- fixture$membership
  if (!is.null(fixture$common)) arguments$common <- fixture$common
  do.call(mixture_regression, arguments)
}

#' Compare one fixture against flexmix
#' @param fixture One fixture.
#' @return A `harness_case`.
#' @noRd
.flexmix_case <- function(fixture) {
  data <- fixture$simulate()
  label <- sprintf("%s, %d classes", fixture$name, fixture$n_classes)
  fit <- .flexmix_package_fit(fixture, data)
  same_mode <- .flexmix_same_mode(label, fixture, data, fit)
  search <- .flexmix_search(label, fixture, data, fit)
  case <- harness_case(rbind(same_mode$comparisons, .with_output_columns(search)),
                       same_mode$ledger)
  case$comparisons <- transform(
    case$comparisons, source = sprintf("flexmix %s", utils::packageVersion("flexmix")),
    dataset = fixture$name, precision = "machine")
  case
}

#' Outputs that are not estimates of the model, with their reasons
#' @noRd
.flexmix_exemptions <- c(
  `fit$converged` = "this optimiser's own stopping flag",
  `fit$iterations` = "this optimiser's iteration count",
  `fit$n_starts` = "how many random starts this program ran: a setting of its search",
  `fit$n_best_replicated` = "how many of this program's random starts reached the best value")

#' The tolerances of this suite
#'
#' Estimates and everything computed from estimates and posteriors alone at
#' the harness default, 1e-6. Standard errors, and the columns computed from
#' them, at 1e-4: both programs take numerical derivatives (see the header).
#' @noRd
.flexmix_tolerances <- c(
  default = 1e-6,
  `coefficients$std_error` = 1e-4, `coefficients$statistic` = 1e-4,
  `coefficients$p_value` = 1e-4, `coefficients$conf_low` = 1e-4,
  `coefficients$conf_high` = 1e-4, `coefficients$p_adjusted` = 1e-4,
  `classes$sigma_std_error` = 1e-4, `classes$prior_std_error` = 1e-4,
  `membership$std_error` = 1e-4, `membership$statistic` = 1e-4,
  `membership$p_value` = 1e-4, `membership$conf_low` = 1e-4,
  `membership$conf_high` = 1e-4)

#' The output columns that are free parameters (or the class priors, a
#' function of the intercept logits alone)
#'
#' These are compared with flexmix's fixed point; every other column with
#' flexmix evaluated at the package's own estimates.
#' @noRd
.flexmix_parameter_columns <- list(
  coefficients = "estimate", classes = c("sigma", "prior"), membership = "estimate")

#' Claim A: every output against flexmix at the same mode
#'
#' Two references are built from flexmix, and each output is compared with the
#' one that isolates what it can be wrong about, as in the hard depmixS4
#' suite:
#'
#'   * the free parameters against flexmix's fixed point, reached by
#'     flexmix's own EM from the package's posteriors: agreement says both
#'     programs report the same mode;
#'   * every function of them -- posteriors, assignments, fitted values,
#'     counts, classification, criteria, standard errors and the Wald columns
#'     -- against flexmix evaluated at the package's own estimates
#'     ([.flexmix_at_estimates()]): the same point, so agreement says the
#'     package computes each of them correctly, whatever its stopping rule
#'     and whatever the noise floor of flexmix's concomitant M-step.
#'
#' At the fixed point instead, a gap of 3.8e-7 in the posteriors on
#' gaussian-2-membership (flexmix's fixed point 8.2e-8 from Newton's maximum
#' in the parameters, the package's 4.6e-8) became 1.8e-6 in the fitted
#' values, which multiply the posteriors by class means 4.5 apart (measured
#' 2026-10-26).
#'
#' @return A list with `comparisons` and `ledger`.
#' @noRd
.flexmix_same_mode <- function(label, fixture, data, fit) {
  obtained <- fit_tables(fit)
  design <- .flexmix_design(fixture, data)
  posterior <- as.matrix(obtained$assignments[sprintf("probability_class_%d",
                                                      seq_len(fixture$n_classes))])
  refined <- .flexmix_fixed_point(fixture, data, posterior)
  ours <- .flexmix_package_estimates(fixture, design, obtained)
  at_package <- .flexmix_at_estimates(refined$fitted, fixture, data, design, ours)
  reference <- .flexmix_reference(fixture, design, at_package, ours)
  expected <- .flexmix_expected(fixture, data, design, reference)
  # At the package's own point the parameter columns are the package's
  # estimates by construction; they are compared with the fixed point instead.
  expected[names(.flexmix_parameter_columns)] <- lapply(
    names(.flexmix_parameter_columns), function(table) {
      out <- expected[[table]]
      if (!is.null(out)) out[.flexmix_parameter_columns[[table]]] <- NULL
      out
    })
  expected <- Filter(Negate(is.null), expected)
  mode <- .flexmix_mode_tables(fixture, .flexmix_estimates(refined$fitted, fixture, design))
  claims <- rbind(
    compare_values(sprintf("%s: flexmix log likelihood at the package's estimates", label),
                   at_package@logLik, fit$log_likelihood, tolerance = 1e-10,
                   scale = "relative"),
    compare_values(sprintf("%s: log likelihood flexmix's EM gains from the package's solution",
                           label), 0, refined$fitted@logLik - fit$log_likelihood,
                   tolerance = 1e-8),
    # Held to the estimate tolerance, not to 1e-9: the concomitant M-step's
    # noise floor reaches 2.2e-7, see .flexmix_fixed_point().
    compare_values(sprintf("%s: largest move of one more flexmix EM step at its fixed point",
                           label), 0, refined$step, tolerance = 1e-6),
    # The from-scratch likelihood whose Hessian gives the reference standard
    # errors must reproduce flexmix's at the same point before its curvature
    # is trusted.
    compare_values(sprintf("%s: from-scratch log likelihood at the package's estimates", label),
                   at_package@logLik, reference$scratch_log_likelihood,
                   tolerance = 1e-10, scale = "relative"))
  comparisons <- rbind(
    compare_tables(label, obtained, expected, .flexmix_tolerances),
    compare_tables(sprintf("%s, against flexmix's fixed point", label), obtained, mode,
                   .flexmix_tolerances),
    .with_output_columns(claims),
    .flexmix_refit_check(label, fixture, at_package, reference))
  list(comparisons = comparisons,
       ledger = output_ledger(label, fit_outputs(fit), comparisons, .flexmix_exemptions))
}

#' Claim B: the package's own search against flexmix's best
#'
#' Ten flexmix searches, each from flexmix's own random initial partition
#' under its own seed, at flexmix's default `minprior` (a class falling under
#' 0.05 is dropped, and such a search, no longer a fit with this many classes,
#' is set aside and counted). The best is refined to its fixed point by
#' [.flexmix_fixed_point()], so the reference is a maximum flexmix's EM
#' reaches, not where its stopping rule fired.
#'
#' @return Comparison rows, none of them an output of the fit.
#' @noRd
.flexmix_search <- function(label, fixture, data, fit) {
  searches <- lapply(seq_len(10L), function(search) {
    .flexmix_with_seed(fixture$seed * 100L + search, flexmix::flexmix(
      fixture$flexmix_formula, data = data, k = fixture$n_classes,
      model = .flexmix_driver(fixture),
      concomitant = .flexmix_concomitant(fixture$membership),
      control = list(iter.max = 2000L, tolerance = 1e-10, verbose = 0)))
  })
  kept <- Filter(function(searched) searched@k == fixture$n_classes, searches)
  stopifnot("at least one flexmix search must keep every class" = length(kept) > 0L)
  log_likelihoods <- vapply(kept, function(searched) searched@logLik, numeric(1))
  best <- kept[[which.max(log_likelihoods)]]
  refined <- .flexmix_fixed_point(fixture, data, flexmix::posterior(best))$fitted@logLik
  rows <- list(compare_values(
    sprintf(paste("%s: shortfall of the package's own search against flexmix (best of %d",
                  "flexmix searches keeping %d classes, of 10; %d within 1e-6 of the best",
                  "before refining)"),
            label, length(kept), fixture$n_classes,
            sum(log_likelihoods >= max(log_likelihoods) - 1e-6)),
    0, max(refined - fit$log_likelihood, 0), tolerance = 1e-6))
  excess <- fit$log_likelihood - refined
  if (excess > 1e-6) {
    rows <- c(rows, list(compare_values(
      sprintf("%s: package maximum above flexmix's best (by %.3g; flexmix's search)",
              label, excess), excess, excess, tolerance = 1e-6)))
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Estimates, and flexmix at the package's estimates
# ---------------------------------------------------------------------------

#' The design of the from-scratch likelihood, built from the data
#' @param fixture One fixture.
#' @param data The data.
#' @return A list of the response, design matrices and unit indices.
#' @noRd
.flexmix_design <- function(fixture, data) {
  binomial <- identical(fixture$family, "binomial")
  unit <- if (identical(fixture$class_level, "group")) match(data[[fixture$id]],
                                                             unique(data[[fixture$id]])) else
    seq_len(nrow(data))
  first <- !duplicated(unit)
  membership <- stats::model.matrix(fixture$membership %||% ~ 1, data)
  list(family = fixture$family, n_classes = fixture$n_classes,
       y = if (binomial) data$s else data$y,
       trials = if (binomial) data$s + data$f else NULL,
       x = stats::model.matrix(fixture$terms, data),
       z = if (is.null(fixture$common)) matrix(0, nrow(data), 0L) else
         stats::model.matrix(fixture$common, data)[, -1L, drop = FALSE],
       w = membership[first, , drop = FALSE],
       unit = unit, first = first, n_units = sum(first),
       cluster = if (is.null(fixture$id)) NULL else
         match(data[[fixture$id]], unique(data[[fixture$id]])),
       equal_variance = identical(fixture$variance, "equal"))
}

#' Estimates in one shape
#'
#' `beta` terms by classes, `common` named by term, `sigma` one per class
#' (Gaussian), `gamma` the concomitant logits, membership terms by classes 2 to
#' K (the intercept logits of the class priors when there are no covariates),
#' and `prior` the class priors when there are no covariates.
#' @noRd
.flexmix_shape <- function(beta, common, sigma, gamma, prior) {
  list(beta = beta, common = common, sigma = sigma, gamma = gamma, prior = prior)
}

#' The estimates of a flexmix fit
#' @param fitted A `flexmix` object, components in the package's class order.
#' @param fixture One fixture.
#' @param design From [.flexmix_design()].
#' @return Estimates, as [.flexmix_shape()].
#' @noRd
.flexmix_estimates <- function(fitted, fixture, design) {
  k <- fixture$n_classes
  stopifnot("flexmix must keep every class" = fitted@k == k)
  estimates <- flexmix::parameters(fitted)
  terms <- colnames(design$x)
  shared_terms <- colnames(design$z)
  # paste0() of an empty vector is "coef.", not empty.
  common <- if (length(shared_terms) == 0L) numeric() else
    stats::setNames(estimates[paste0("coef.", shared_terms), 1L], shared_terms)
  prior <- flexmix::prior(fitted)
  gamma <- if (is.null(fixture$membership)) {
    matrix(log(prior[-1L]) - log(prior[[1L]]), 1L, dimnames = list("(Intercept)", NULL))
  } else {
    concomitant <- as.matrix(flexmix::parameters(fitted, which = "concomitant"))
    concomitant[colnames(design$w), -1L, drop = FALSE]
  }
  .flexmix_shape(
    beta = matrix(estimates[paste0("coef.", terms), ], length(terms), k,
                  dimnames = list(terms, NULL)),
    common = common,
    sigma = if (identical(fixture$family, "gaussian")) unname(estimates["sigma", ]) else
      numeric(),
    gamma = gamma, prior = if (is.null(fixture$membership)) prior)
}

#' The package's estimates, read from its reported tables
#' @param fixture One fixture.
#' @param design From [.flexmix_design()].
#' @param obtained The fit's tables, from `fit_tables()`.
#' @return Estimates, as [.flexmix_shape()].
#' @noRd
.flexmix_package_estimates <- function(fixture, design, obtained) {
  k <- fixture$n_classes
  classes <- sprintf("class_%d", seq_len(k))
  coefficients <- obtained$coefficients
  lookup <- stats::setNames(coefficients$estimate,
                            paste(coefficients$class, coefficients$term, sep = "|"))
  terms <- colnames(design$x)
  shared_terms <- colnames(design$z)
  beta <- matrix(lookup[paste(rep(classes, each = length(terms)), terms, sep = "|")],
                 length(terms), k, dimnames = list(terms, NULL))
  # paste() with an empty vector gives "common|", not an empty vector.
  common <- if (length(shared_terms) == 0L) numeric() else
    stats::setNames(lookup[paste("common", shared_terms, sep = "|")], shared_terms)
  prior <- obtained$classes$prior[match(classes, obtained$classes$class)]
  gamma <- if (is.null(fixture$membership)) {
    matrix(log(prior[-1L]) - log(prior[[1L]]), 1L, dimnames = list("(Intercept)", NULL))
  } else {
    membership <- obtained$membership
    logits <- stats::setNames(membership$estimate,
                              paste(membership$class, membership$term, sep = "|"))
    matrix(logits[paste(rep(classes[-1L], each = ncol(design$w)), colnames(design$w),
                        sep = "|")], ncol(design$w), k - 1L,
           dimnames = list(colnames(design$w), NULL))
  }
  out <- .flexmix_shape(
    beta = beta, common = common,
    sigma = if (identical(fixture$family, "gaussian"))
      obtained$classes$sigma[match(classes, obtained$classes$class)] else numeric(),
    gamma = gamma, prior = if (is.null(fixture$membership)) prior)
  stopifnot("every package estimate must be found in its tables" =
              !anyNA(unlist(out)))
  out
}

#' flexmix evaluated at given estimates
#'
#' The flexmix object at its fixed point with its components rebuilt from the
#' estimates by the driver's own `defineComponent()` (as flexmix's FLXboot
#' rebuilds them from bootstrap estimates), its concomitant coefficients or
#' prior replaced, and its posteriors and log-likelihood recomputed by
#' flexmix's `posterior(newdata =)`, which runs flexmix's E-step: component
#' log densities, flexmix's grouping of rows, and `prior(newdata =)`.
#'
#' @param fitted flexmix's fixed point.
#' @param fixture One fixture.
#' @param data The data.
#' @param design From [.flexmix_design()].
#' @param estimates As [.flexmix_shape()].
#' @return A `flexmix` object whose components, prior, posteriors and
#'   log-likelihood are those at `estimates`.
#' @noRd
.flexmix_at_estimates <- function(fitted, fixture, data, design, estimates) {
  at <- fitted
  driver <- fitted@model[[1L]]
  at@components <- lapply(seq_len(fixture$n_classes), function(class) {
    component <- fitted@components[[class]][[1L]]
    parameters <- component@parameters
    values <- c(estimates$beta[, class], estimates$common)
    stopifnot("every flexmix coefficient must have a package estimate" =
                all(names(parameters$coef) %in% names(values)))
    parameters$coef[] <- values[names(parameters$coef)]
    if (identical(fixture$family, "gaussian")) parameters$sigma <- estimates$sigma[[class]]
    list(driver@defineComponent(c(parameters, list(df = component@df))))
  })
  if (is.null(fixture$membership)) {
    at@prior <- estimates$prior
  } else {
    at@concomitant@coef[rownames(estimates$gamma), -1L] <- estimates$gamma
  }
  scaled <- flexmix::posterior(at, newdata = data)
  unscaled <- flexmix::posterior(at, newdata = data, unscaled = TRUE)
  at@posterior <- list(scaled = scaled, unscaled = unscaled)
  # flexmix's grouped posteriors repeat each group's value on its rows; the
  # likelihood counts each group once, as flexmix's EM does (`groupfirst`).
  at@logLik <- sum(log(rowSums(unscaled[design$first, , drop = FALSE])))
  at
}

# ---------------------------------------------------------------------------
# The from-scratch likelihood and the reference standard errors
# ---------------------------------------------------------------------------

#' Natural parameters as one vector
#' @param design From [.flexmix_design()].
#' @param estimates As [.flexmix_shape()].
#' @return Class coefficients (terms within class), shared coefficients,
#'   residual SDs (Gaussian; one when shared), concomitant logits (terms
#'   within class, classes 2 to K).
#' @noRd
.flexmix_pack <- function(design, estimates) {
  c(as.vector(estimates$beta), unname(estimates$common),
    if (identical(design$family, "gaussian"))
      if (design$equal_variance) estimates$sigma[[1L]] else estimates$sigma,
    as.vector(estimates$gamma))
}

#' Split a natural-parameter vector into its blocks
#' @param design From [.flexmix_design()].
#' @param theta As [.flexmix_pack()] returns.
#' @return A list of `beta`, `common`, `sigma` (one per class) and `gamma`.
#' @noRd
.flexmix_unpack <- function(design, theta) {
  k <- design$n_classes
  p <- ncol(design$x)
  q <- ncol(design$z)
  s <- if (!identical(design$family, "gaussian")) 0L else if (design$equal_variance) 1L else k
  r <- ncol(design$w)
  at <- cumsum(c(0L, p * k, q, s))
  list(beta = matrix(theta[seq_len(p * k)], p, k),
       common = theta[at[[2L]] + seq_len(q)],
       sigma = if (s == 0L) numeric() else rep_len(theta[at[[3L]] + seq_len(s)], k),
       gamma = matrix(theta[at[[4L]] + seq_len(r * (k - 1L))], r, k - 1L))
}

#' Each unit's log-likelihood contribution, from scratch
#'
#' Within class k a row's density is the family's at the linear predictor
#' `x beta_k + z common`; a unit's classes are a multinomial logit of its
#' concomitant covariates with the first class as reference; a unit is a row,
#' or a group whose rows share their class. Log-sum-exp over the classes.
#'
#' @param design From [.flexmix_design()].
#' @param theta Natural parameters, as [.flexmix_pack()] returns.
#' @return The vector of unit contributions.
#' @noRd
.flexmix_unit_log_likelihood <- function(design, theta) {
  parameters <- .flexmix_unpack(design, theta)
  k <- design$n_classes
  shared <- if (ncol(design$z) > 0L) drop(design$z %*% parameters$common) else 0
  eta <- design$x %*% parameters$beta + shared
  density <- vapply(seq_len(k), function(class) {
    switch(design$family,
      gaussian = stats::dnorm(design$y, eta[, class], parameters$sigma[[class]], log = TRUE),
      binomial = stats::dbinom(design$y, design$trials, stats::plogis(eta[, class]),
                               log = TRUE),
      poisson = stats::dpois(design$y, exp(eta[, class]), log = TRUE))
  }, numeric(nrow(eta)))
  by_unit <- rowsum(matrix(density, nrow(eta)), design$unit, reorder = TRUE)
  logits <- cbind(0, design$w %*% parameters$gamma)
  log_prior <- logits - .flexmix_log_sum_exp(logits)
  .flexmix_log_sum_exp(log_prior + by_unit)
}

#' Row-wise log-sum-exp
#' @param m A matrix.
#' @return A vector.
#' @noRd
.flexmix_log_sum_exp <- function(m) {
  peak <- apply(m, 1L, max)
  peak + log(rowSums(exp(m - peak)))
}

#' Reference values at one point: flexmix's functions and the from-scratch
#' standard errors
#'
#' Posteriors, class means, log-likelihood, parameter count, AIC, row BIC and
#' relative entropy are flexmix's (`posterior()`, `fitted()`, `logLik`, `df`,
#' `EIC()`, and stats::AIC()/BIC() of flexmix's logLik object, whose `df` and
#' `nobs` are flexmix's) at the flexmix object given. Standard errors are the
#' from-scratch observed information (or its cluster sandwich) at `estimates`.
#'
#' @param fixture One fixture.
#' @param design From [.flexmix_design()].
#' @param at A `flexmix` object at `estimates` ([.flexmix_at_estimates()]).
#' @param estimates As [.flexmix_shape()].
#' @return A list of reference values.
#' @noRd
.flexmix_reference <- function(fixture, design, at, estimates) {
  theta <- .flexmix_pack(design, estimates)
  unit_log_likelihood <- function(t) .flexmix_unit_log_likelihood(design, t)
  information <- -numDeriv::hessian(function(t) sum(unit_log_likelihood(t)), theta)
  bread <- solve(information)
  observed_errors <- .flexmix_unpack(design, sqrt(diag(bread)))
  covariance <- if (identical(fixture$vcov_type, "robust")) {
    # Unit scores (rows here: classes are row-level) summed within clusters,
    # with the CR1 factor G / (G - 1).
    scores <- rowsum(numDeriv::jacobian(unit_log_likelihood, theta), design$cluster)
    clusters <- nrow(scores)
    bread %*% crossprod(scores) %*% bread * clusters / (clusters - 1)
  } else bread
  gamma_at <- length(theta) - length(estimates$gamma) + seq_along(estimates$gamma)
  flexmix_log_likelihood <- flexmix::logLik(at)
  posterior <- flexmix::posterior(at)
  class_means <- flexmix::fitted(at)
  if (identical(fixture$family, "binomial")) class_means <- class_means * design$trials
  list(estimates = estimates, errors = .flexmix_unpack(design, sqrt(diag(covariance))),
       observed_errors = observed_errors,
       gamma_covariance = covariance[gamma_at, gamma_at, drop = FALSE],
       posterior = posterior, unit_posterior = posterior[design$first, , drop = FALSE],
       class_means = class_means, log_likelihood = at@logLik, n_parameters = at@df,
       # flexmix's logLik object carries its df and its nobs (rows); the
       # log-likelihood is recomputed here because `at` was not fitted.
       aic = -2 * at@logLik + 2 * attr(flexmix_log_likelihood, "df"),
       bic_rows = -2 * at@logLik +
         attr(flexmix_log_likelihood, "df") * log(attr(flexmix_log_likelihood, "nobs")),
       n_observations = attr(flexmix_log_likelihood, "nobs"),
       entropy = flexmix::EIC(at),
       scratch_log_likelihood = sum(unit_log_likelihood(theta)))
}

#' flexmix's own observed-information errors against the reference's
#'
#' `refit(method = "optim")` at the package's estimates: the class (and shared)
#' coefficients' errors, and the concomitant coefficients' errors where there
#' is a concomitant model, against the reference's inverse observed
#' information (the bread, also in the robust case, whose reported errors are
#' the sandwich). A check that the reference is flexmix's estimator, outside
#' the ledger. flexmix's Hessian is optim's finite difference of its analytic
#' gradient (`ndeps` 1e-3), so the two are held to the standard-error
#' tolerance, 1e-4; they agree to 1.2e-6 relative or better (measured
#' 2026-10-26).
#'
#' @param label Label of the case.
#' @param fixture One fixture.
#' @param at flexmix at the package's estimates.
#' @param reference From [.flexmix_reference()].
#' @return Comparison rows with `output = NA`.
#' @noRd
.flexmix_refit_check <- function(label, fixture, at, reference) {
  # The refitted object holds one coefficient matrix (estimate, standard
  # error, z, p) per class and per concomitant logit block; summary() would
  # only print them.
  refitted <- flexmix::refit(at, method = "optim")
  k <- fixture$n_classes
  terms <- rownames(reference$estimates$beta)
  shared_terms <- names(reference$estimates$common)
  # One coefficient table per class; a shared coefficient appears in every
  # class's table, and is read from the first.
  components <- refitted@components[[1L]]
  refit_errors <- c(
    unlist(lapply(seq_len(k), function(class) {
      unname(components[[class]][terms, "Std. Error"])
    })),
    unname(components[[1L]][shared_terms, "Std. Error"]))
  observed <- reference$observed_errors
  rows <- compare_block(
    sprintf("%s: flexmix refit(method = \"optim\") coefficient errors against the reference",
            label),
    c(as.vector(observed$beta), observed$common), refit_errors, tolerance = 1e-4)
  if (is.null(fixture$membership)) return(rows)
  membership_terms <- rownames(reference$estimates$gamma)
  concomitant_errors <- unlist(lapply(refitted@concomitant, function(block) {
    unname(block[membership_terms, "Std. Error"])
  }))
  rbind(rows, compare_block(
    sprintf("%s: flexmix refit(method = \"optim\") concomitant errors against the reference",
            label),
    as.vector(observed$gamma), concomitant_errors, tolerance = 1e-4))
}

# ---------------------------------------------------------------------------
# Expected tables
# ---------------------------------------------------------------------------

#' Wald columns, by definition
#' @param estimate,std_error Vectors.
#' @return A `data.frame` with `estimate`, `std_error`, `statistic`,
#'   `p_value`, `conf_low`, `conf_high`.
#' @noRd
.flexmix_wald <- function(estimate, std_error) {
  critical <- stats::qnorm(0.975)
  statistic <- estimate / std_error
  data.frame(estimate = estimate, std_error = std_error, statistic = statistic,
             p_value = 2 * stats::pnorm(-abs(statistic)),
             conf_low = estimate - critical * std_error,
             conf_high = estimate + critical * std_error)
}

#' The parameter columns at flexmix's fixed point, in the fit's shapes
#' @param fixture One fixture.
#' @param estimates flexmix's fixed point, as [.flexmix_shape()].
#' @return A named list of tables.
#' @noRd
.flexmix_mode_tables <- function(fixture, estimates) {
  k <- fixture$n_classes
  classes <- sprintf("class_%d", seq_len(k))
  terms <- rownames(estimates$beta)
  shared_terms <- names(estimates$common)
  class_rows <- data.frame(class = classes)
  if (identical(fixture$family, "gaussian")) class_rows$sigma <- estimates$sigma
  if (is.null(fixture$membership)) class_rows$prior <- estimates$prior
  Filter(Negate(is.null), list(
    coefficients = data.frame(
      class = c(rep(classes, each = length(terms)), rep("common", length(shared_terms))),
      term = c(rep(terms, k), shared_terms),
      estimate = c(as.vector(estimates$beta), unname(estimates$common))),
    classes = if (ncol(class_rows) > 1L) class_rows,
    membership = data.frame(
      model = "class", class = rep(classes[-1L], each = nrow(estimates$gamma)),
      term = rep(rownames(estimates$gamma), k - 1L), reference = "class_1",
      estimate = as.vector(estimates$gamma))))
}

#' Every expected table, in the fit's shapes, from one reference point
#'
#' Definitions as in `?get_results` and `.reference_level_tables()` of the
#' harness, whose average-posterior, count and entropy definitions are those
#' of the regression tables: modal class by first maximum, average posteriors
#' of the units assigned to each class (`NA` for an empty class), shares as
#' column means of the unit posteriors.
#'
#' @param fixture One fixture.
#' @param data The data.
#' @param design From [.flexmix_design()].
#' @param reference From [.flexmix_reference()].
#' @return A named list of tables.
#' @noRd
.flexmix_expected <- function(fixture, data, design, reference) {
  k <- fixture$n_classes
  classes <- sprintf("class_%d", seq_len(k))
  estimates <- reference$estimates
  errors <- reference$errors
  terms <- rownames(estimates$beta)
  shared_terms <- names(estimates$common)
  non_gaussian <- !identical(fixture$family, "gaussian")

  coefficient_rows <- data.frame(
    class = c(rep(classes, each = length(terms)), rep("common", length(shared_terms))),
    term = c(rep(terms, k), shared_terms),
    .flexmix_wald(c(as.vector(estimates$beta), unname(estimates$common)),
                  c(as.vector(errors$beta), errors$common)))
  slopes <- coefficient_rows$term != "(Intercept)"
  coefficient_rows$p_adjusted <- NA_real_
  coefficient_rows$p_adjusted[slopes] <- stats::p.adjust(coefficient_rows$p_value[slopes],
                                                         method = "BH")
  if (non_gaussian) coefficient_rows$exp_estimate <- exp(coefficient_rows$estimate)

  unit_posterior <- reference$unit_posterior
  n_units <- nrow(unit_posterior)
  units <- .reference_level_tables(unit_posterior, "units", seq_len(n_units))
  quality <- units$classification
  class_rows <- data.frame(class = classes, share = quality$estimated_proportion,
                           count = quality$estimated_n, n_assigned = quality$n_modal,
                           mean_posterior = quality$average_posterior)
  if (!non_gaussian) {
    class_rows$sigma <- estimates$sigma
    class_rows$sigma_std_error <- errors$sigma
  }
  if (is.null(fixture$membership)) {
    p <- estimates$prior
    jacobian <- (diag(p, k) - tcrossprod(p))[, -1L, drop = FALSE]
    class_rows$prior <- p
    class_rows$prior_std_error <- sqrt(diag(jacobian %*% reference$gamma_covariance %*%
                                              t(jacobian)))
  }

  # Without covariates the membership table holds the intercept logits of the
  # class priors.
  membership_terms <- rownames(estimates$gamma)
  membership_rows <- data.frame(
    model = "class", class = rep(classes[-1L], each = length(membership_terms)),
    term = rep(membership_terms, k - 1L), reference = "class_1",
    .flexmix_wald(as.vector(estimates$gamma), as.vector(errors$gamma)))
  membership_rows$odds_ratio <- exp(membership_rows$estimate)

  posterior <- reference$posterior
  rows <- seq_len(nrow(posterior))
  modal <- max.col(posterior, ties.method = "first")
  log_likelihood <- reference$log_likelihood
  q <- reference$n_parameters
  bic <- -2 * log_likelihood + q * log(n_units)
  fit_row <- data.frame(
    n_group_classes = 1L, n_observations = reference$n_observations,
    n_groups = if (is.null(fixture$id)) NA_integer_ else length(unique(data[[fixture$id]])),
    log_likelihood = log_likelihood, n_parameters = q, aic = reference$aic, bic = bic,
    bic_rows = reference$bic_rows,
    sabic = -2 * log_likelihood + q * log((n_units + 2) / 24),
    icl = bic + 2 * units$entropy$entropy_sum, entropy = reference$entropy,
    # No second level, so no group-class entropy.
    group_entropy = NA_real_, smallest_share = min(quality$estimated_proportion))

  assignments <- data.frame(.row = rows, row = rows, class = modal,
                            posterior = posterior[cbind(rows, modal)])
  assignments[sprintf("probability_class_%d", seq_len(k))] <- as.data.frame(unname(posterior))

  means <- reference$class_means
  fitted_rows <- data.frame(row = rows, observed = design$y,
                            fitted = rowSums(posterior * means),
                            fitted_modal = means[cbind(rows, modal)], class = modal)
  fitted_rows$residual <- fitted_rows$observed - fitted_rows$fitted

  averages <- units$average_posteriors
  classification <- data.frame(assigned = classes[averages$assigned_class],
                               class = classes[averages$class],
                               mean_posterior = averages$average_posterior,
                               n_assigned = averages$n_assigned)

  groups <- if (identical(fixture$class_level, "group")) {
    unit_modal <- max.col(unit_posterior, ties.method = "first")
    out <- data.frame(group = data[[fixture$id]][design$first],
                      n_rows = tabulate(design$unit, n_units), class = unit_modal,
                      posterior = unit_posterior[cbind(seq_len(n_units), unit_modal)])
    names(out)[[1L]] <- fixture$id
    out[sprintf("probability_class_%d", seq_len(k))] <- as.data.frame(unname(unit_posterior))
    out
  }

  Filter(Negate(is.null), list(
    coefficients = coefficient_rows, classes = class_rows, membership = membership_rows,
    fit = fit_row, assignments = assignments, fitted = fitted_rows,
    classification = classification, groups = groups))
}
