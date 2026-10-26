# Adversarial latent transition data against depmixS4 (Visser & Speekenbrink
# 2010, Journal of Statistical Software 36(7):1-21).
#
# The depmixS4 suite compares lta() on four well-behaved simulations: states
# 0.75 persistent, means two standard deviations apart, six or seven
# occasions. This suite simulates the data on which hidden Markov software goes
# wrong: transitions within 0.005 of the boundary, a state holding 5% of the
# occasions, sequences of two occasions, sequences of 2 to 20 occasions,
# states whose responses overlap heavily, a category one state never uses, a
# transition that is structurally zero, an indicator measured in units of
# 1e-4, states that differ only in their variances, and one state more than
# the data were generated with.
#
# Every case makes two separate claims, so that a disagreement says whether
# the package computes wrongly or searches badly:
#   (A) computation at one mode. lta() takes no starting values (multilpa()
#       has `start =`, lta() does not), so the two programs are put at the same
#       mode from the package's side: depmixS4 is set to the package's
#       estimates, its log-likelihood there is read (the same likelihood
#       function), and its own EM map is iterated from there to its fixed point
#       (`.hard_depmixs4_fixed_point()`). The log-likelihood that gains says
#       whether the package stopped at a fixed point. The free parameters are
#       compared with depmixS4's fixed point in units of their standard errors
#       (to 1e-3; absolutely, to 1e-6, where Wald inference is refused), and
#       every function of them
#       (posteriors, counts, classification tables, criteria, residuals,
#       diagnostics, thresholds, standard errors) with depmixS4 evaluated at
#       the package's own estimates, so a stopping rule cannot pass for a
#       computing error (see `.hard_depmixs4_same_mode()`). The tables come from
#       the depmixS4 suite's builders, sourced below, not copied.
#   (B) search. The package's own random starts (30, tol 1e-12) against
#       depmixS4's best over five of its own random starts, each taken to a
#       fixed point. A shortfall of the package's maximum is a finding about
#       its optimiser; a maximum above depmixS4's is recorded as a row that is
#       not an output of the fit, and is a finding about depmixS4's search.
#
# Findings (2026-10-26, depmixS4 1.5.1, package 0.97.1):
#   * No computing disagreement on any fixture: every function of the
#     parameters agrees with depmixS4 at the same point (posteriors, counts
#     and log-likelihoods to 4e-15; standard errors with the boundary-free
#     reference to 2e-9, to 2.4e-7 on variance-only and 2.3e-5 on
#     weak-separation, whose information has condition number 2.1e4).
#   * weak-separation (an observation, not a defect): the package stops 5.7e-5
#     from the maximum in the means (4.1e-5 in the transitions), 5.5e-9 below
#     it in log-likelihood -- 8.5e-5 of a standard error in the means and
#     9.4e-5 at most. The surface is flat (condition number 2.1e4 in the
#     reference's coordinates, EM contraction 0.99929 per step), and the
#     log-likelihood stopping rule (`tol = 1e-12`, an absolute change of
#     4.3e-9) fires there; the package, like Mplus, Latent GOLD, depmixS4,
#     mclust and poLCA, stops on a log-likelihood change. `tol = 1e-14` still
#     stops 3.7e-5 away; plain EM (`acceleration = "none"`) at 1e-12 stops
#     1.3e-3 away after 12,763 iterations. No estimate, interval or test can
#     differ at that distance, so the parameters are compared in standard
#     errors (see `.hard_depmixs4_parameter_distance()`).
#   * depmixS4 itself: on weak-separation all five of its searches reach
#     `maxit = 5000` without converging (reported only in the fitted object's
#     `message` slot) 1e-3 below the maximum, and are continued here; on
#     over-extracted its best search is 1.83 below the package's maximum
#     (0.50 below over 20 seeds), which the same-mode claim shows to be a
#     fixed point of depmixS4's own EM.
#   * The package floors every probability at `min_probability` (1e-10) and
#     refuses Wald standard errors when one sits there; depmixS4 lets it fall
#     to zero. The log-likelihood that costs is within the bound stated in
#     `.hard_depmixs4_same_mode()`.
#
# Standard errors. The depmixS4 suite's reference Hessian is taken in the
# natural coordinates with the last entry of every distribution as the
# complement; numDeriv's first Richardson step is 10% of each value, so on a
# distribution with an entry under a tenth of another the step leaves the
# simplex. On the rare-state fixture that returned NaN for every standard
# error (measured 2026-10-26). The reference here is the same forward
# recursion written in coordinates with no boundary -- baseline-category
# logits for every distribution, standardised means, log variances -- and
# mapped back by the delta method. At a maximum the inverse information
# transforms with the Jacobian of the reparameterisation, so the means' and
# variances' errors are the same estimand, and the two references are checked
# against each other wherever the natural one is defined.
#
# depmixS4's own `standardError()` (nlme::fdHess) is not used: the depmixS4
# suite documents that it is 1 to 15 percent off on regular data, and its
# replacement there (`.depmixs4_own_errors()`, numDeriv in depmixS4's own
# natural parameters) has the boundary-crossing step above.

sys.source(file.path("equivalence", "harness", "suites", "depmixs4.R"), envir = environment())

# The depmixS4 suite's natural-coordinate reference, kept under its own name
# for the cross-check; `.depmixs4_reference()` resolves
# `.depmixs4_standard_errors` in this environment and so uses the
# boundary-free version defined at the end of this file.
.depmixs4_natural_errors <- .depmixs4_standard_errors

suite_hard_depmixs4 <- function() {
  require_suite_packages(
    "depmixS4", "numDeriv",
    reason = paste("the independent hidden-Markov reference these adversarial",
                   "transition fixtures are compared against, and the reference",
                   "observed information"))
  stack_cases(lapply(.hard_depmixs4_fixtures(), .hard_depmixs4_case))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' A transition matrix with one persistence on the diagonal
#' @param k Number of states.
#' @param stay Probability of staying.
#' @return A `k` by `k` matrix, rows summing to one.
#' @noRd
.hard_depmixs4_stay <- function(k, stay) {
  transition <- matrix((1 - stay) / (k - 1), k, k)
  diag(transition) <- stay
  transition
}

#' The adversarial fixtures
#'
#' Each fixture is a generating model and the number of states fitted. Two
#' Gaussian indicators `y1`, `y2` (and a three-category `c1` where given), as
#' the depmixS4 suite's builders expect.
#'
#' @return A list of fixture specifications.
#' @noRd
.hard_depmixs4_fixtures <- function() {
  two <- function(m1, m2) cbind(c(-m1, m1), c(-m2, m2))
  sds <- function(k, value = 1) matrix(value, k, 2L)
  base <- list(categorical = NULL, scale = c(y1 = 1, y2 = 1), min_variance = 1e-6)
  fixture <- function(...) utils::modifyList(base, list(...))
  # Sequence lengths drawn once, with their own seed, so the fixture is the
  # same whatever the simulation seed draws.
  unequal <- local({
    set.seed(2040L)
    c(2L, 20L, sample(2:20, 98L, replace = TRUE))
  })
  list(
    # 1. Near-absorbing: stay 0.995, measurement weak enough (means 1 SD
    #    apart on y1) that persistence carries much of the information. Over
    #    150 x 12 occasions about 8 moves are expected.
    fixture(name = "near-absorbing", seed = 301L, profiles = 2L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.995),
            means = two(0.5, 0.4), sds = sds(2L), lengths = rep(12L, 150L)),
    # 2. A rare state: initial share 0.05 and stationary share 0.048
    #    (pi_3 = 0.02 (pi_1 + pi_2) / 0.4).
    fixture(name = "rare-state", seed = 302L, profiles = 3L,
            initial = c(0.5, 0.45, 0.05),
            transition = rbind(c(0.85, 0.13, 0.02), c(0.14, 0.84, 0.02),
                               c(0.20, 0.20, 0.60)),
            means = cbind(c(-1.5, 0, 2), c(-1, 0, 1.5)), sds = sds(3L, 0.8),
            lengths = rep(8L, 150L)),
    # 3. Two occasions per sequence: one transition each, 500 sequences.
    fixture(name = "two-occasions", seed = 303L, profiles = 2L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.8),
            means = two(0.8, 0.6), sds = sds(2L), lengths = rep(2L, 500L)),
    # 4. Sequence lengths from 2 to 20 occasions, both extremes present.
    fixture(name = "unequal-2-to-20", seed = 304L, profiles = 3L,
            initial = c(0.4, 0.3, 0.3), transition = .hard_depmixs4_stay(3L, 0.8),
            means = cbind(c(-2, 0, 2), c(-1.5, 0, 1.5)), sds = sds(3L, 0.9),
            lengths = unequal),
    # 5. Weak separation: means a quarter of a standard deviation from zero,
    #    Mahalanobis distance between the states about 0.6.
    fixture(name = "weak-separation", seed = 305L, profiles = 2L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.8),
            means = two(0.25, 0.2), sds = sds(2L), lengths = rep(6L, 250L)),
    # 6. State 1 never uses category "c": its probability's maximum is on the
    #    floor.
    fixture(name = "category-on-floor", seed = 306L, profiles = 2L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.8),
            means = two(1, 0.8), sds = sds(2L), lengths = rep(6L, 120L),
            categorical = rbind(c(0.6, 0.4, 0), c(0.2, 0.3, 0.5))),
    # 7. Three states, the move 1 -> 3 structurally zero. The outer states are
    #    7 SD apart on y1, so no misclassified pair puts expected mass on the
    #    move and its maximum is on the bound: the package holds it at its
    #    floor under seeds 1 to 4 and 307 (with means -2, 0, 2 and SD 0.8 it
    #    was interior, 0.0033, under seed 307; measured 2026-10-26).
    fixture(name = "structural-zero", seed = 307L, profiles = 3L,
            initial = c(0.4, 0.3, 0.3),
            transition = rbind(c(0.8, 0.2, 0), c(0.1, 0.8, 0.1), c(0.1, 0.1, 0.8)),
            means = cbind(c(-2.5, 0, 2.5), c(-2, 0, 2)), sds = sds(3L, 0.7),
            lengths = rep(8L, 120L)),
    # 8. y1 in units of 1e-4: its variances are about 1e-8. The variance floor
    #    is absolute (`min_variance`, 1e-6 by default), so it is carried into
    #    the new units with the square of the scale, as the invariance suite
    #    does; left at 1e-6 it sits a hundred times above every y1 variance and
    #    the fit is bound-active by construction, not by the data.
    fixture(name = "scale-1e-4", seed = 308L, profiles = 2L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.8),
            means = two(1, 0.8), sds = sds(2L), lengths = rep(6L, 100L),
            scale = c(y1 = 1e-4, y2 = 1), min_variance = 1e-6 * 1e-4^2),
    # 9. States with equal means that differ only in their spread (SD 0.5
    #    against 1.5): no location separates them. The reference builders
    #    match states by their estimated means, which is safe only because
    #    depmixS4 is started from the package's solution, so its states are
    #    in the package's order and its means within 1e-8 of the package's,
    #    far closer than the estimated means of the two states are to each
    #    other.
    fixture(name = "variance-only", seed = 309L, profiles = 2L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.85),
            means = two(0, 0), sds = cbind(c(0.5, 1.5), c(0.5, 1.5)),
            lengths = rep(8L, 150L)),
    # 10. One state too many: three states fitted to two-state data, so the
    #     likelihood has ridges and several modes. The package's 30 starts end
    #     at five different maxima and reach the best (-2915.286, with the
    #     move 2 -> 3 on its bound) three times; depmixS4's searches stop at
    #     lower ones (measured 2026-10-26).
    fixture(name = "over-extracted", seed = 310L, profiles = 3L,
            initial = c(0.5, 0.5), transition = .hard_depmixs4_stay(2L, 0.8),
            means = two(1, 0.8), sds = sds(2L), lengths = rep(6L, 150L))
  )
}

#' Simulate one fixture
#' @param fixture One element of [.hard_depmixs4_fixtures()].
#' @return A data frame with `g`, `t`, the generating `state`, `y1`, `y2` and
#'   optionally `c1`, ordered by group and occasion.
#' @noRd
.hard_depmixs4_simulate <- function(fixture) {
  set.seed(fixture$seed)
  k <- length(fixture$initial)
  frame <- do.call(rbind, lapply(seq_along(fixture$lengths), function(group) {
    path <- unlist(Reduce(function(previous, step) {
      sample.int(k, 1L, prob = fixture$transition[previous, ])
    }, seq_len(fixture$lengths[[group]] - 1L),
    init = sample.int(k, 1L, prob = fixture$initial), accumulate = TRUE))
    data.frame(g = group, t = seq_along(path), state = path)
  }))
  draw <- function(column) {
    fixture$scale[[column]] * stats::rnorm(nrow(frame), fixture$means[frame$state, column],
                                           fixture$sds[frame$state, column])
  }
  frame$y1 <- draw(1L)
  frame$y2 <- draw(2L)
  if (!is.null(fixture$categorical)) {
    frame$c1 <- factor(vapply(frame$state, function(state) {
      sample(c("a", "b", "c"), 1L, prob = fixture$categorical[state, ])
    }, character(1)), levels = c("a", "b", "c"))
  }
  stopifnot("the simulated sequences must have the lengths asked for" =
              identical(as.vector(table(frame$g)), as.integer(fixture$lengths)))
  frame
}

# ---------------------------------------------------------------------------
# One case: claim A (same mode) and claim B (search)
# ---------------------------------------------------------------------------

#' Compare one adversarial fixture against depmixS4
#' @param fixture One element of [.hard_depmixs4_fixtures()].
#' @return A `harness_case`.
#' @noRd
.hard_depmixs4_case <- function(fixture) {
  frame <- .hard_depmixs4_simulate(fixture)
  categorical <- !is.null(fixture$categorical)
  indicators <- c("y1", "y2", if (categorical) "c1")
  k <- fixture$profiles
  label <- sprintf("%s, %d states", fixture$name, k)
  # Warnings are not suppressed: the harness counts them and records the first.
  fit <- lta(frame, indicators, "g", n_profiles = k, time = "t",
             categorical = if (categorical) "c1" else character(),
             n_starts = 30, seed = fixture$seed, max_iter = 5000, tol = 1e-12,
             min_variance = fixture$min_variance)
  formulas <- lapply(indicators, function(name) stats::as.formula(paste(name, "~ 1")))
  families <- c(list(stats::gaussian(), stats::gaussian()),
                if (categorical) list(depmixS4::multinomial("identity")))
  ntimes <- as.integer(table(frame$g)[unique(as.character(frame$g))])
  model <- depmixS4::depmix(formulas, data = frame, nstates = k, family = families,
                            ntimes = ntimes)

  search <- .hard_depmixs4_search(label, model, fit, fixture$seed,
                                  if (categorical) length(fit$categorical_levels$c1) else 0L)
  same_mode <- .hard_depmixs4_same_mode(label, model, fit, frame, indicators, ntimes,
                                        categorical)
  case <- harness_case(rbind(.with_output_columns(same_mode$comparisons), search),
                       same_mode$ledger)
  case$comparisons <- transform(
    case$comparisons, source = sprintf("depmixS4 %s", utils::packageVersion("depmixS4")),
    dataset = fixture$name, precision = "machine")
  case
}

#' Claim B: the package's own search against depmixS4's best
#'
#' depmixS4's `fit()` draws a random start by default (`em.control(random.start
#' = TRUE)`); five seeds give five independent searches, and the best of them
#' is depmixS4's maximum.
#'
#' A search that reaches `maxit` is not at a maximum. depmixS4 says so only in
#' the fitted object's `message` slot -- no warning, no flag -- and on
#' weak-separation all five searches stopped there, 1e-3 below the maximum
#' (EM contracts by 0.99929 per step on that fixture; measured 2026-10-26).
#' Such a search is continued to its fixed point by depmixS4's own EM map
#' ([.hard_depmixs4_fixed_point()]), so the reference is a maximum depmixS4's
#' algorithm reaches from its own start, and the label counts how many were
#' continued.
#'
#' @return Comparison rows, none of them an output of the fit.
#' @noRd
.hard_depmixs4_search <- function(label, model, fit, seed, n_categories) {
  searches <- lapply(seq_len(5L), function(start) {
    set.seed(seed * 100L + start)
    searched <- depmixS4::fit(model, verbose = FALSE, emcontrol =
      depmixS4::em.control(maxit = 5000, tol = 1e-12))
    # The text depmixS4's em() stores when its iteration bound is reached.
    stopped <- grepl("without convergence", searched@message, fixed = TRUE)
    if (stopped) searched <- .hard_depmixs4_fixed_point(searched, fit$n_profiles, n_categories)
    list(log_likelihood = as.numeric(depmixS4::logLik(searched)), stopped = stopped)
  })
  depmix_ll <- vapply(searches, `[[`, numeric(1), "log_likelihood")
  stopped <- sum(vapply(searches, `[[`, logical(1), "stopped"))
  best <- max(depmix_ll)
  # The label says how many of depmixS4's searches reached its best, so a
  # multimodal fixture is visible in the report.
  rows <- list(compare_values(
    sprintf(paste("%s: shortfall of the package's own search against depmixS4 (best of 5",
                  "depmixS4 searches, reached by %d; %d stopped at maxit and were continued",
                  "to their fixed point)"),
            label, sum(depmix_ll >= best - 1e-6), stopped),
    0, max(best - fit$log_likelihood, 0), tolerance = 1e-6))
  excess <- fit$log_likelihood - best
  if (excess > 1e-6) {
    # The package found a higher maximum than any of depmixS4's searches: a
    # finding about depmixS4's search, not a package defect.
    rows <- c(rows, list(compare_values(
      sprintf("%s: package maximum above depmixS4's best (by %.3g; depmixS4's search)",
              label, excess), excess, excess, tolerance = 1e-6)))
  }
  do.call(rbind, rows)
}

#' The output columns that are free parameters of the model
#'
#' Everything else a fit reports is a function of these and of the data.
#' @noRd
.hard_depmixs4_parameters <- list(
  profiles = c("mean", "variance"), responses = "probability",
  transitions = "probability", initial = "probability")

#' Claim A: every output of the fit against depmixS4 at the same mode
#'
#' Two references are built from depmixS4, and each output is compared with the
#' one that isolates what it can be wrong about:
#'
#'   * the free parameters (means, variances, initial, transition and response
#'     probabilities) against depmixS4's fixed point, reached by depmixS4's own
#'     EM from the package's solution, as a distance in standard errors (see
#'     [.hard_depmixs4_parameter_distance()]): agreement says both programs
#'     report the same mode for every inferential purpose;
#'   * every function of them -- posteriors, counts, sequence log-likelihoods,
#'     classification tables, criteria, residuals, diagnostics, standard
#'     deviations, thresholds, standard errors -- against depmixS4 evaluated at
#'     the package's own estimates: the same point, so agreement says the
#'     package computes each of them correctly, whatever its stopping rule.
#'
#' Comparing the functions at depmixS4's fixed point instead would charge them
#' with the distance between the two stopping points: on the near-absorbing
#' and rare-state fixtures the package stops about 2e-6 short of the mode (one
#' Newton step on the reference recursion, in standard-deviation and logit
#' units; log-likelihood 3e-12 and 5e-11 below it), which is within its
#' `tol = 1e-12` and moved the posteriors by 1.4e-6 (measured 2026-10-26).
#'
#' @return A list with `comparisons` and `ledger`.
#' @noRd
.hard_depmixs4_same_mode <- function(label, model, fit, frame, indicators, ntimes,
                                     categorical) {
  started <- depmixS4::setpars(model, .depmixs4_pars(fit, categorical))
  at_ours <- as.numeric(depmixS4::logLik(started))
  fixed <- .hard_depmixs4_fixed_point(started, fit$n_profiles,
                                      if (categorical) length(fit$categorical_levels$c1) else 0L)
  at_package <- .depmixs4_reference(started, fit, frame, indicators, ntimes, categorical)
  at_mode <- .depmixs4_reference(fixed, fit, frame, indicators, ntimes, categorical)
  expected <- .depmixs4_expected(at_package, fit, frame, ntimes)
  at_mode_tables <- .depmixs4_expected(at_mode, fit, frame, ntimes)
  # The parameter columns are compared with the fixed point below, in standard
  # errors; at the package's own point they would equal its estimates by
  # construction, so they are taken out of the same-point comparison.
  tables <- intersect(names(.hard_depmixs4_parameters),
                      names(Filter(Negate(is.null), expected)))
  expected[tables] <- lapply(tables, function(table) {
    out <- expected[[table]]
    out[.hard_depmixs4_parameters[[table]]] <- NULL
    out
  })
  bound <- .hard_depmixs4_bound(at_mode$estimates)
  # Refining from a fixed point gains nothing, so a gain says the package
  # stopped short of its mode. A probability the package holds at its floor
  # (min_probability) is the one exception: depmixS4 has no floor, and moving
  # such a probability from the floor to zero rescales the other entries of
  # its distribution by 1 / (1 - min_probability), which gains at most
  # min_probability on each occasion it governs, so at most
  # min_probability * n_observations per bound probability (7.2e-8 on
  # category-on-floor, where 1.9e-8 is gained).
  gain_tolerance <- 1e-8 + length(bound) * fit$min_probability * fit$n_observations
  claims <- rbind(
    compare_values(sprintf("%s: log likelihood at the package's estimates", label),
                   at_ours, fit$log_likelihood, tolerance = 1e-8, scale = "relative"),
    compare_values(sprintf("%s: log likelihood depmixS4's EM gains from the package's solution",
                           label), 0, as.numeric(depmixS4::logLik(fixed)) - at_ours,
                   tolerance = gain_tolerance))
  exempt <- .hard_depmixs4_exemptions(bound, fit)
  expected <- .drop_exempt(expected, exempt)
  obtained <- fit_tables(fit)
  obtained$model$.row <- 1L
  comparisons <- compare_tables(label, obtained, expected, .output_default_tolerances)
  comparisons <- rbind(comparisons,
                       .hard_depmixs4_parameter_distance(label, fit, obtained, at_mode_tables),
                       .with_output_columns(claims),
                       .hard_depmixs4_reference_check(label, at_package, frame, ntimes))
  list(comparisons = comparisons,
       ledger = output_ledger(label, fit_outputs(fit), comparisons, exempt))
}

#' The free parameters against depmixS4's fixed point, in standard errors
#'
#' Each estimate's distance from depmixS4's fixed point is divided by the
#' estimate's own standard error (`parameter_inference(fit, boundary = "fix")`,
#' observed information), and held to 1e-3: agreement to a thousandth of a
#' standard error means no estimate, interval or test computed from the fit can
#' differ. That two programs stopping on a log-likelihood change (the rule of
#' this package, Mplus, Latent GOLD, depmixS4, mclust and poLCA) stop within
#' 1e-6 of each other in absolute terms is not a property of either: on
#' weak-separation the package stops 5.7e-5 from the maximum in a mean, which
#' is 8.5e-5 of its standard error (9.4e-5 at most, for a transition
#' probability; measured 2026-10-26). Machine-precision agreement of the computation is
#' established separately, by comparing every function of the parameters with
#' depmixS4 at the package's own point.
#'
#' A cell with no standard error -- every cell of a fit where Wald inference is
#' refused because a probability is on its bound (`latents_boundary_fit`) --
#' keeps an absolute comparison with the fixed point at the harness default
#' (1e-6), in its own block under the same output key, so the ledger counts
#' every value once.
#'
#' The inference rows are matched to the output cells by name: `mean` and
#' `variance` have `outcome` the profile and `term` the indicator; `response`
#' has `term` `indicator:category`; `initial_probability` has `term` the group
#' class; `transition_probability` has `outcome` the profile moved to and
#' `term` `group_class:profile moved from` (checked against the transitions
#' table: the rows of one term sum to one).
#'
#' @param obtained The fit's tables, from `fit_tables()`.
#' @param at_mode_tables Tables built from depmixS4's fixed point.
#' @return Comparison rows with `output` set.
#' @noRd
.hard_depmixs4_parameter_distance <- function(label, fit, obtained, at_mode_tables) {
  inference <- tryCatch(parameter_inference(fit, boundary = "fix"),
                        latents_boundary_fit = function(condition) NULL)
  lookup <- if (is.null(inference)) character() else stats::setNames(
    inference$standard_error,
    paste(inference$parameter, inference$outcome, inference$term, sep = "|"))
  key_of <- function(table, frame) {
    profile <- function(p) sprintf("profile_%d", as.integer(p))
    switch(table,
      profiles = paste(NA, profile(frame$profile), frame$indicator, sep = "|"),
      responses = paste("response", profile(frame$profile),
                        paste(frame$indicator, frame$category, sep = ":"), sep = "|"),
      initial = paste("initial_probability", profile(frame$profile),
                      sprintf("group_class_%d", as.integer(frame$group_class)), sep = "|"),
      transitions = paste("transition_probability", profile(frame$to),
                          sprintf("group_class_%d:%s", as.integer(frame$group_class),
                                  profile(frame$from)), sep = "|"))
  }
  if (!is.null(inference)) {
    transitions <- subset(inference, parameter == "transition_probability")
    stopifnot("transition inference rows must be grouped by the profile moved from" =
                all(abs(tapply(transitions$estimate, transitions$term, sum) - 1) < 1e-8))
  }
  tables <- intersect(names(.hard_depmixs4_parameters), names(at_mode_tables))
  rows <- lapply(tables, function(table) {
    have <- obtained[[table]]
    want <- at_mode_tables[[table]]
    keys <- .output_keys(want, table)
    as_key <- function(frame) do.call(paste, c(lapply(frame[keys], as.character), sep = "\r"))
    at <- match(as_key(have), as_key(want))
    stopifnot("every parameter cell must be in the fixed-point reference" = !anyNA(at))
    lapply(.hard_depmixs4_parameters[[table]], function(column) {
      cell_key <- key_of(table, have)
      if (identical(table, "profiles")) cell_key <- sub("^NA", column, cell_key)
      error <- unname(lookup[cell_key])
      error[!is.finite(error) | error <= 0] <- NA_real_
      estimate <- have[[column]]
      reference <- want[[column]][at]
      with_error <- !is.na(error)
      output <- paste0(table, "$", column)
      rbind(
        if (any(with_error)) compare_block(
          sprintf("%s: %s, distance from depmixS4's fixed point in standard errors",
                  label, output),
          rep(0, sum(with_error)),
          abs(estimate[with_error] - reference[with_error]) / error[with_error],
          tolerance = 1e-3, output = output),
        if (any(!with_error)) compare_block(
          sprintf("%s: %s, no standard error (Wald inference refused), against depmixS4's fixed point",
                  label, output),
          reference[!with_error], estimate[!with_error],
          tolerance = .output_default_tolerances[["default"]], output = output))
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

#' depmixS4's fixed point, from a starting solution
#'
#' depmixS4's EM map is iterated from `started` until one step moves no
#' parameter by more than 1e-12 (`.depmixs4_refine()`'s rule). On a flat
#' likelihood the plain iteration is too slow for the harness: on
#' weak-separation its contraction rate is 0.99929 per step and it needs
#' 14,627 steps (145 s, measured 2026-10-26). The map is therefore iterated in
#' SQUAREM cycles (Varadhan & Roland 2008, scheme S3): two depmixS4 EM steps,
#' an extrapolation along them, and a depmixS4 EM step from the extrapolated
#' point, kept only when every probability and standard deviation stays
#' positive and the log-likelihood does not fall below the second step's. The
#' extrapolation is an affine combination of depmixS4 parameter vectors, so
#' every distribution still sums to one. Acceleration changes how fast the
#' fixed point is reached, not which point it is: the result is accepted only
#' after `.depmixs4_refine()` finds that a plain depmixS4 EM step from it moves
#' nothing by more than 1e-12.
#'
#' @param started A depmixS4 model with its parameters set.
#' @param k Number of states.
#' @param n_categories Categories of the one categorical response, or 0.
#' @param max_cycles Bound on the SQUAREM cycles.
#' @return The depmixS4 model at its fixed point.
#' @noRd
.hard_depmixs4_fixed_point <- function(started, k, n_categories, max_cycles = 300L) {
  em <- function(object) {
    depmixS4::fit(object, verbose = FALSE, emcontrol = depmixS4::em.control(
      maxit = 0, tol = 0, crit = "absolute", random.start = FALSE))
  }
  log_likelihood <- function(object) as.numeric(depmixS4::logLik(object))
  # depmixS4's layout: prior, transition rows, then per state the Gaussian
  # intercept and SD of y1 and y2 and the category probabilities.
  positive <- c(rep(TRUE, k + k^2),
                rep(c(FALSE, TRUE, FALSE, TRUE, rep(TRUE, n_categories)), k))
  stopifnot("the positivity mask must cover depmixS4's parameters" =
              length(positive) == length(depmixS4::getpars(started)))
  state <- Reduce(function(state, cycle) {
    if (state$settled) return(state)
    p0 <- depmixS4::getpars(state$model)
    one <- em(state$model)
    p1 <- depmixS4::getpars(one)
    if (max(abs(p1 - p0)) < 1e-12) return(list(model = state$model, settled = TRUE))
    two <- em(one)
    p2 <- depmixS4::getpars(two)
    r <- p1 - p0
    v <- p2 - 2 * p1 + p0
    alpha <- min(-sqrt(sum(r^2) / sum(v^2)), -1)
    jumped <- p0 - 2 * alpha * r + alpha^2 * v
    following <- two
    if (is.finite(alpha) && all(jumped[positive] > 0)) {
      candidate <- em(depmixS4::setpars(state$model, jumped))
      if (is.finite(log_likelihood(candidate)) &&
          log_likelihood(candidate) >= log_likelihood(two)) following <- candidate
    }
    list(model = following, settled = FALSE)
  }, seq_len(max_cycles), init = list(model = started, settled = FALSE))
  .depmixs4_refine(state$model, max_steps = 2000L)
}

#' Outputs a reference at a boundary cannot speak to, with the reason
#'
#' A probability is on its bound when depmixS4's refined estimate is below
#' 1e-8: EM moves a probability whose maximum is zero geometrically towards
#' zero and never reaches it, while the package holds it at `min_probability`
#' (1e-10). There the log-likelihood has no curvature in the logit, the
#' observed information is singular, and the package refuses Wald inference
#' (`latents_boundary_fit`, documented); the reference has no standard error
#' to compare either. (Thresholds next to such a probability, `qlogis(1 - p)`,
#' are compared at the package's own estimates, where both programs hold the
#' same `p`; at depmixS4's fixed point `p` keeps falling towards zero and the
#' threshold has no limit.)
#'
#' @param bound From [.hard_depmixs4_bound()] at depmixS4's fixed point.
#' @param fit The fit.
#' @return A named character vector of reasons.
#' @noRd
.hard_depmixs4_exemptions <- function(bound, fit) {
  if (length(bound) == 0L) return(character())
  where <- sprintf("%s (depmixS4's refined value %s, the package's %s)",
                   names(bound), format(bound, digits = 3),
                   format(.hard_depmixs4_bound_package(fit, names(bound)), digits = 3))
  reason <- sprintf(paste(
    "no standard error exists at this mode: %s is on its bound, where the",
    "information is singular; the package refuses Wald inference",
    "(latents_boundary_fit) and the reference has no curvature to invert"),
    paste(where, collapse = "; "))
  outputs <- fit_outputs(fit)$output
  errors <- intersect(c("profiles$mean_standard_error", "profiles$variance_standard_error",
                        "responses$probability_standard_error"), outputs)
  stats::setNames(rep(reason, length(errors)), errors)
}

#' Probabilities of a reference solution that sit on their lower bound
#' @param estimates Aligned estimates from [.depmixs4_reference()].
#' @return Named numeric vector of the bound values (empty when interior).
#' @noRd
.hard_depmixs4_bound <- function(estimates) {
  k <- length(estimates$initial)
  states <- seq_len(k)
  values <- c(
    stats::setNames(estimates$initial, sprintf("initial[%d]", states)),
    stats::setNames(as.vector(t(estimates$transition)),
                    sprintf("transition[%d->%d]", rep(states, each = k), rep(states, k))),
    if (!is.null(estimates$responses)) {
      block <- estimates$responses$c1
      stats::setNames(as.vector(block),
                      sprintf("c1[state %d, %s]", rep(states, ncol(block)),
                              rep(colnames(block), each = k)))
    })
  values[values < 1e-8]
}

#' The package's value of each bound probability, for the stated reason
#' @param fit The fit.
#' @param names Names from [.hard_depmixs4_bound()].
#' @return A numeric vector.
#' @noRd
.hard_depmixs4_bound_package <- function(fit, names) {
  k <- fit$n_profiles
  states <- seq_len(k)
  values <- c(
    stats::setNames(fit$initial_probabilities[1L, ], sprintf("initial[%d]", states)),
    stats::setNames(as.vector(t(fit$transition_probabilities[, , 1L])),
                    sprintf("transition[%d->%d]", rep(states, each = k), rep(states, k))),
    if (!is.null(fit$response_probabilities$c1)) {
      block <- fit$response_probabilities$c1
      stats::setNames(as.vector(block),
                      sprintf("c1[state %d, %s]", rep(states, ncol(block)),
                              rep(colnames(block), each = k)))
    })
  unname(values[names])
}

#' The two standard-error references against each other
#'
#' Where the depmixS4 suite's natural-coordinate Hessian stays inside the
#' parameter space its errors and this suite's logit-coordinate errors are the
#' same estimand computed two ways; their agreement is a check on the
#' reference, outside the ledger. Both are numerical second derivatives, so
#' they are held to the harness's standard-error tolerance (1e-4, see
#' `.output_default_tolerances`); they agree to 1e-8 where the information is
#' well conditioned, to 1.4e-7 on variance-only and to 1.4e-5 on
#' weak-separation, whose information has condition number 2.1e4 in the logit
#' coordinates (measured 2026-10-26).
#' Where a natural-coordinate step leaves the parameter space the natural
#' reference is undefined and no row is written.
#'
#' @return Comparison rows with `output = NA`, or `NULL`.
#' @noRd
.hard_depmixs4_reference_check <- function(label, reference, frame, ntimes) {
  logit <- reference$standard_errors
  if (is.null(logit)) return(NULL)
  estimates <- reference$estimates
  # numDeriv's first Richardson step on a coordinate x is 0.1 |x|, or 1e-4
  # when |x| is below its zero.tol, sqrt(.Machine$double.eps / 7e-7) = 1.8e-5.
  # Probabilities: the complement (the last entry) absorbs the steps of the
  # free ones, so it must exceed a tenth of their sum. Variances: one below
  # zero.tol takes the absolute step of 1e-4 and goes negative (y1 in units of
  # 1e-4 has variances near 1e-8, and the natural reference returns NaN).
  rows <- c(list(estimates$initial), lapply(seq_len(nrow(estimates$transition)),
                                            function(i) estimates$transition[i, ]),
            if (!is.null(estimates$responses))
              lapply(seq_len(nrow(estimates$responses$c1)),
                     function(i) estimates$responses$c1[i, ]))
  inside <- all(vapply(rows, function(p) {
    p[[length(p)]] > 0.1 * sum(p[-length(p)]) && all(p > 0)
  }, logical(1))) && all(estimates$variances > sqrt(.Machine$double.eps / 7e-7))
  if (!inside) return(NULL)
  natural <- .depmixs4_natural_errors(estimates, frame, ntimes, reference$log_likelihood)
  out <- rbind(
    compare_block(sprintf("%s: logit- against natural-coordinate reference SEs, means", label),
                  natural$mean, logit$mean, tolerance = 1e-4),
    compare_block(sprintf("%s: logit- against natural-coordinate reference SEs, variances", label),
                  natural$variance, logit$variance, tolerance = 1e-4))
  if (!is.null(natural$probability)) {
    out <- rbind(out, compare_block(
      sprintf("%s: logit- against natural-coordinate reference SEs, category probabilities",
              label), natural$probability$c1, logit$probability$c1, tolerance = 1e-4))
  }
  out
}

# ---------------------------------------------------------------------------
# Reference standard errors in coordinates with no boundary
# ---------------------------------------------------------------------------

#' Observed-information standard errors of a latent transition model, from scratch
#'
#' The same scaled forward recursion as the depmixS4 suite's reference, over
#' every sequence at once, parameterised without a boundary: baseline-category
#' logits (last entry the baseline) for the initial distribution, every
#' transition row and every state's category probabilities; each mean as
#' `(mean - centre) / spread` of its column; each variance as
#' `log(variance / spread^2)`. Its Hessian at depmixS4's estimates is taken by
#' Richardson extrapolation (numDeriv) and inverted; the means and variances
#' map back by their (diagonal) Jacobian and the category probabilities by the
#' softmax Jacobian `diag(p) - p p'`, which also gives the last category's
#' error.
#'
#' A solution with a probability on its bound has no standard error (see
#' [.hard_depmixs4_exemptions()]), and `NULL` is returned.
#'
#' @param estimates Aligned estimates from [.depmixs4_reference()].
#' @param frame The data.
#' @param ntimes Sequence lengths, in data order.
#' @param log_likelihood depmixS4's maximum, which the recursion must
#'   reproduce at its estimates before its Hessian is trusted.
#' @return `NULL`, or a list with `mean`, `variance` and, with a categorical
#'   indicator, `probability`, as `harness_reference()` takes them.
#' @noRd
.depmixs4_standard_errors <- function(estimates, frame, ntimes, log_likelihood) {
  if (length(.hard_depmixs4_bound(estimates)) > 0L) return(NULL)
  k <- length(estimates$initial)
  continuous <- colnames(estimates$means)
  j <- length(continuous)
  y <- as.matrix(frame[continuous])
  centre <- colMeans(y)
  spread <- apply(y, 2L, stats::sd)
  codes <- if (!is.null(estimates$responses)) as.integer(frame$c1)
  n_categories <- if (is.null(codes)) 0L else ncol(estimates$responses$c1)
  last <- cumsum(ntimes)
  slot <- matrix(t(vapply(seq_along(ntimes), function(i) {
    c(last[[i]] - ntimes[[i]] + seq_len(ntimes[[i]]),
      rep(NA_integer_, max(ntimes) - ntimes[[i]]))
  }, integer(max(ntimes)))), length(ntimes))

  sizes <- c(initial = k - 1L, transition = k * (k - 1L), mean = k * j,
             variance = k * j, probability = k * max(n_categories - 1L, 0L))
  ends <- cumsum(sizes)
  at <- function(name) ends[[name]] - sizes[[name]] + seq_len(sizes[[name]])
  # Rows of baseline-category logits to rows of probabilities.
  softmax_rows <- function(logits, width) {
    full <- cbind(matrix(logits, ncol = width - 1L), 0)
    full <- full - apply(full, 1L, max)
    exp(full) / rowSums(exp(full))
  }
  logits_of <- function(p) log(p[, -ncol(p), drop = FALSE]) - log(p[, ncol(p)])
  unpack <- function(theta) {
    list(initial = as.vector(softmax_rows(matrix(theta[at("initial")], 1L), k)),
         transition = softmax_rows(matrix(theta[at("transition")], k, byrow = TRUE), k),
         means = sweep(sweep(matrix(theta[at("mean")], k, j), 2L, spread, "*"),
                       2L, centre, "+"),
         variances = sweep(exp(matrix(theta[at("variance")], k, j)), 2L, spread^2, "*"),
         probability = if (n_categories > 0L)
           softmax_rows(matrix(theta[at("probability")], k), n_categories))
  }
  log_likelihood_at <- function(theta) {
    p <- unpack(theta)
    log_emission <- matrix(vapply(seq_len(k), function(state) {
      gaussian <- rowSums(vapply(seq_len(j), function(column) {
        stats::dnorm(y[, column], p$means[state, column], sqrt(p$variances[state, column]),
                     log = TRUE)
      }, numeric(nrow(y))))
      if (is.null(codes)) gaussian else gaussian + log(p$probability[state, codes])
    }, numeric(nrow(y))), nrow(y))
    emit <- function(rows) {
      block <- log_emission[rows, , drop = FALSE]
      peak <- apply(block, 1L, max)
      list(density = exp(block - peak), peak = peak)
    }
    start <- emit(slot[, 1L])
    alpha <- sweep(start$density, 2L, p$initial, "*")
    scale <- rowSums(alpha)
    state <- list(alpha = alpha / scale, total = start$peak + log(scale))
    # One matrix product per occasion, over the sequences that reach it.
    state <- Reduce(function(state, occasion) {
      rows <- slot[, occasion]
      active <- !is.na(rows)
      step <- emit(rows[active])
      alpha <- (state$alpha[active, , drop = FALSE] %*% p$transition) * step$density
      scale <- rowSums(alpha)
      state$alpha[active, ] <- alpha / scale
      state$total[active] <- state$total[active] + step$peak + log(scale)
      state
    }, seq_len(ncol(slot))[-1L], init = state)
    sum(state$total)
  }
  theta <- c(logits_of(matrix(estimates$initial, 1L)),
             as.vector(t(logits_of(estimates$transition))),
             as.vector(sweep(sweep(estimates$means, 2L, centre, "-"), 2L, spread, "/")),
             as.vector(log(sweep(estimates$variances, 2L, spread^2, "/"))),
             if (n_categories > 0L) as.vector(logits_of(estimates$responses$c1)))
  stopifnot("the logit-coordinate forward recursion must reproduce depmixS4's maximum" =
              abs(log_likelihood_at(theta) - log_likelihood) < 1e-8 * abs(log_likelihood))
  covariance <- solve(-numDeriv::hessian(log_likelihood_at, theta))
  errors <- sqrt(diag(covariance))
  out <- list(
    mean = sweep(matrix(errors[at("mean")], k, j), 2L, spread, "*"),
    variance = estimates$variances * matrix(errors[at("variance")], k, j))
  dimnames(out$mean) <- dimnames(out$variance) <- list(NULL, continuous)
  if (n_categories > 0L) {
    free <- matrix(at("probability"), k)
    out$probability <- list(c1 = t(vapply(seq_len(k), function(state) {
      p <- estimates$responses$c1[state, ]
      jacobian <- (diag(p) - tcrossprod(p))[, -n_categories, drop = FALSE]
      block <- covariance[free[state, ], free[state, ], drop = FALSE]
      sqrt(diag(jacobian %*% block %*% t(jacobian)))
    }, numeric(n_categories))))
    dimnames(out$probability$c1) <- list(NULL, colnames(estimates$responses$c1))
  }
  out
}
