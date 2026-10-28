# Mixture latent transition analysis against LMest (Bartolucci, Pandolfi &
# Pennoni 2017, Journal of Statistical Software 81(4):1-38).
#
# `lta(n_group_classes = K1)` with categorical indicators fits the mixed latent
# Markov model of Van de Pol & Langeheine (1990): each sequence belongs to one
# of K1 group classes, each class has its own initial distribution and
# homogeneous transition matrix over K2 states, and the state-conditional
# response probabilities are shared by every class and occasion. LMest's
# `lmestMixed(k1 = K1, k2 = K2)` fits the same model by EM (`lmmixed()`). The
# parameter blocks map one to one:
#
#   LMest la[u]            <-> group_probabilities[u]
#   LMest Piv[v, u]        <-> initial_probabilities[u, v]
#   LMest Pi[v0, v1, u]    <-> transition_probabilities[v0, v1, u] (rows: from)
#   LMest Psi[y + 1, v, j] <-> response_probabilities[[j]][v, y]
#   LMest np, lk           <-> n_parameters, log_likelihood
#
# LMest's count, (K1 - 1) + K1 (K2^2 - 1) + K2 r b, is the package's: 13 for
# two classes, two states and three binary items, 19 for three-category ones,
# 35 for three states and three-category items.
#
# The depmixS4 suites never have more than one group class, so until this suite
# the transitions per group class, the sequences' group-class posteriors and
# the group rows of the classification tables of a transition fit were not
# compared with anything. LMest reports its estimates and likelihood only, so
# everything else is computed by a forward-backward recursion over (group
# class x state) written from the model's definition in base R
# (lta-mixed-reference.R), after it has reproduced LMest's own log-likelihood
# and its group-class posteriors per sequence (`W`) at LMest's estimates.
#
# Every case makes three claims, as hard-depmixs4 does, so that a disagreement
# says whether the package computes wrongly, stops early or searches badly:
#
#   (A) computation. lta() takes no starting values and lmestMixed() takes none
#       either (`start` is 0, deterministic, or 1, random), so the two
#       programs cannot be put at one point through their APIs. Every function
#       of the parameters -- posteriors at both levels, sequence
#       log-likelihoods, expected transition counts, prevalences, the
#       classification tables of both levels, BCH weights, entropy, every
#       information criterion, residuals, thresholds, the sequence tables --
#       is therefore compared with the base-R recursion evaluated at the
#       package's own estimates, to 1e-10 on the max(|x|, 1) scale (rounding
#       is 1e-15), so a stopping rule cannot pass for a computing error. Standard errors are compared
#       there too, with the inverse of a numDeriv Hessian of the same
#       recursion in baseline-category logits, to 1e-4 relative.
#   (B) the mode. The free parameters (initial, transition, response and
#       group-class probabilities) are compared with LMest's solution carried
#       to its EM fixed point, in units of their standard errors (to 1e-3).
#       LMest stops when one EM step raises the log-likelihood by less than
#       `tol` times its size (1e-10 here, about 2e-7), which leaves it short
#       of the fixed point; it accepts no starting values, so its solution is
#       continued by the same EM map written in base R
#       (`.lta_mixed_fixed_point()`).
#   (C) search. The package's own maximum against the best of three LMest
#       searches (its deterministic start and two random ones), refined as in
#       (B).
#
# `parameter_inference()` holds the standard errors of the transition,
# initial and group-class probabilities, which no fit table reports; its
# table is compared on every column through the verb ledger, against rows
# generated from the model's dimensions so that a dropped row is caught.
#
# LMest prints its iteration log with cat(); it is captured with
# capture.output(), which discards console text and leaves warnings and
# errors to reach the harness.
#
# Van de Pol, F. & Langeheine, R. (1990). Mixed Markov latent class models.
# Sociological Methodology, 20, 213-247.

source(file.path("equivalence", "harness", "lta-mixed-reference.R"), local = environment())

suite_lta_mixed <- function() {
  require_suite_packages(
    "LMest", "numDeriv",
    reason = paste("the independent mixed latent Markov implementation this suite",
                   "compares group-class transitions against, and the reference",
                   "observed information"))
  stack_cases(lapply(.lta_mixed_fixtures(), .lta_mixed_case))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' The fixtures: two generating mixtures of chains, three measurement models
#'
#' The first two: 200 sequences of four occasions, two group classes (shares
#' 0.6 and 0.4) and two states, whose chains differ in direction: the first
#' persistent (0.8, 0.75 on the diagonal), the second switching (0.7, 0.6
#' off it). Three items with two
#' categories in one fixture and three in the other; LMest needs every item to
#' have the same number of categories. Each seed was chosen from six
#' (401-406, 411-416) as the one whose maximum is interior (every probability
#' above 0.03), reached by all ten starts of a tuning fit, and has the
#' smallest largest standard error (0.09 and 0.13; measured 2026-10-27). Near
#' a boundary no standard error exists; hard-depmixs4 covers that case.
#'
#' The third fixture has two group classes and three states, so that no
#' array of the model is square in its group-class and state axes: a
#' transposed initial table, a transition array read with its class axis in
#' the wrong place, or a term parsed into the wrong axis changes a shape or a
#' value instead of passing unseen, as it could when K1 = K2. Its group
#' classes are a persistent chain and a cyclic one (0.7 on the step 1 -> 2 ->
#' 3 -> 1); 300 sequences of three-category items. Its seed was chosen from
#' twelve (421-432) by the same rule; it is the only one of the twelve with
#' every probability above 0.03 (0.032) and all ten starts at the maximum,
#' largest standard error 0.08 (measured 2026-10-27).
#' @noRd
.lta_mixed_fixtures <- function() {
  chain <- list(la = c(0.6, 0.4), initial = rbind(c(0.65, 0.35), c(0.35, 0.65)),
                transition = list(rbind(c(0.8, 0.2), c(0.25, 0.75)),
                                  rbind(c(0.3, 0.7), c(0.6, 0.4))))
  list(
    c(list(name = "binary-items", seed = 405L,
           responses = list(rbind(c(0.9, 0.1), c(0.1, 0.9)),
                            rbind(c(0.85, 0.15), c(0.15, 0.85)),
                            rbind(c(0.9, 0.1), c(0.15, 0.85)))), chain),
    c(list(name = "three-category-items", seed = 412L,
           responses = list(rbind(c(0.7, 0.2, 0.1), c(0.1, 0.25, 0.65)),
                            rbind(c(0.65, 0.25, 0.1), c(0.15, 0.3, 0.55)),
                            rbind(c(0.75, 0.15, 0.1), c(0.1, 0.2, 0.7)))), chain),
    list(name = "three-states", seed = 424L, n = 300L, la = c(0.55, 0.45),
         initial = rbind(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6)),
         transition = list(rbind(c(0.8, 0.1, 0.1), c(0.1, 0.8, 0.1), c(0.1, 0.1, 0.8)),
                           rbind(c(0.2, 0.7, 0.1), c(0.1, 0.2, 0.7), c(0.7, 0.1, 0.2))),
         responses = list(rbind(c(0.75, 0.15, 0.1), c(0.15, 0.7, 0.15), c(0.1, 0.15, 0.75)),
                          rbind(c(0.7, 0.2, 0.1), c(0.2, 0.65, 0.15), c(0.1, 0.2, 0.7)),
                          rbind(c(0.8, 0.1, 0.1), c(0.15, 0.7, 0.15), c(0.1, 0.1, 0.8)))))
}

#' Simulate one fixture
#' @param fixture One element of [.lta_mixed_fixtures()].
#' @param n Number of sequences.
#' @param occasions Occasions per sequence.
#' @return A data frame with `id`, `occ` and integer items `y1`, `y2`, `y3`
#'   coded from 0 (LMest's coding), ordered by sequence and occasion.
#' @noRd
.lta_mixed_simulate <- function(fixture, n = fixture[["n"]] %||% 200L, occasions = 4L) {
  set.seed(fixture$seed)
  k2 <- ncol(fixture$initial)
  group_class <- sample.int(length(fixture$la), n, replace = TRUE, prob = fixture$la)
  paths <- vapply(seq_len(n), function(i) {
    u <- group_class[[i]]
    unlist(Reduce(function(previous, step) {
      sample.int(k2, 1L, prob = fixture$transition[[u]][previous, ])
    }, seq_len(occasions - 1L), init = sample.int(k2, 1L, prob = fixture$initial[u, ]),
    accumulate = TRUE))
  }, integer(occasions))
  state <- as.vector(paths)
  frame <- data.frame(id = rep(seq_len(n), each = occasions),
                      occ = rep(seq_len(occasions), times = n))
  items <- lapply(fixture$responses, function(block) {
    vapply(state, function(s) sample.int(ncol(block), 1L, prob = block[s, ]) - 1L, integer(1))
  })
  names(items) <- sprintf("y%d", seq_along(items))
  frame <- data.frame(frame, items)
  stopifnot("every item must use every category, as LMest's coding assumes" =
              all(vapply(seq_along(items), function(j) {
                identical(sort(unique(items[[j]])), seq_len(ncol(fixture$responses[[j]])) - 1L)
              }, logical(1))))
  frame
}

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

#' Compare one fixture against LMest
#' @param fixture One element of [.lta_mixed_fixtures()].
#' @return A `harness_case`.
#' @noRd
.lta_mixed_case <- function(fixture) {
  occasions <- 4L
  frame <- .lta_mixed_simulate(fixture, occasions = occasions)
  items <- grep("^y[0-9]+$", names(frame), value = TRUE)
  # The package takes the items as factors, so each is categorical with the
  # categories LMest numbers from 0.
  categorical_frame <- frame
  categorical_frame[items] <- lapply(frame[items], factor)
  k1 <- length(fixture$la)
  k2 <- ncol(fixture$initial)
  label <- sprintf("%s, %d group classes x %d states", fixture$name, k1, k2)
  # Warnings are not suppressed: the harness counts them and records the first.
  fit <- lta(categorical_frame, items, "id", n_profiles = k2, n_group_classes = k1,
             time = "occ", categorical = items, n_starts = 10, seed = fixture$seed,
             max_iter = 5000, tol = 1e-12)
  stopifnot("the package's categories must be LMest's codes 0, 1, ..." =
              all(vapply(items, function(item) {
                identical(fit$categorical_levels[[item]],
                          as.character(seq_along(fit$categorical_levels[[item]]) - 1L))
              }, logical(1))))
  layout <- .lta_mixed_codes(frame, fit$categorical_levels[items], occasions)

  searches <- .lta_mixed_lmest(frame, items, k1, k2, fit$categorical_levels)
  lmest_ll <- vapply(searches, `[[`, numeric(1), "log_likelihood")
  best <- searches[[which.max(lmest_ll)]]
  # The recursion must reproduce LMest's own log-likelihood at LMest's
  # estimates before any of its other output is trusted. Both sum the same
  # logs over the same sequences, so they agree to rounding (measured at most
  # 1.1e-15 relative on the three fixtures, 2026-10-27); 1e-9 relative is far
  # below anything a wrong recursion could produce.
  at_lmest_evaluated <- .lta_mixed_evaluate(best$parameters, layout)
  at_lmest <- at_lmest_evaluated$log_likelihood
  evaluator_check <- rbind(
    compare_values(
      sprintf("%s: base-R recursion against LMest's log likelihood, both at LMest's estimates (reference self-check)",
              label), best$log_likelihood, at_lmest, tolerance = 1e-9, scale = "relative"),
    # LMest's W is computed by the E-step that follows its last M-step, so it
    # is at the parameters LMest returns. It ties the recursion's
    # group-class posteriors to LMest directly, not only through the
    # likelihood. Same quantity at the same point: they differ by rounding
    # (measured at most 1.3e-15, 2026-10-27); 1e-9 is far below any error in
    # the recursion and far above rounding.
    compare_block(
      sprintf("%s: base-R recursion's group-class posteriors against LMest's W, both at LMest's estimates (reference self-check)",
              label), best$group_posteriors, at_lmest_evaluated$group_posteriors,
      tolerance = 1e-9))

  # (B) LMest's solution continued to its EM fixed point, in the fit's order.
  refined <- .lta_mixed_fixed_point(best$parameters, layout)
  at_mode <- .lta_mixed_align(refined$parameters, fit, items)
  mode_ll <- .lta_mixed_log_likelihood(at_mode, layout)
  stopifnot("continuing LMest's EM must not lose likelihood" = mode_ll >= at_lmest - 1e-9)

  # (A) every function of the parameters at the package's own estimates.
  at_package <- .lta_mixed_fit_parameters(fit, items)
  evaluated <- .lta_mixed_evaluate(at_package, layout)
  errors <- .lta_mixed_standard_errors(at_package, layout)
  expected <- .lta_mixed_expected(evaluated, at_package, best$n_parameters, fit, frame,
                                  categorical_frame, items)
  obtained <- fit_tables(fit)
  obtained$model$.row <- 1L

  claims <- rbind(
    compare_values(sprintf("%s: log likelihood at the package's estimates (base-R recursion)", label),
                   evaluated$log_likelihood, fit$log_likelihood, tolerance = 1e-8,
                   scale = "relative"),
    .lta_mixed_search(label, lmest_ll, mode_ll, refined$cycles, fit$log_likelihood),
    # LMest's own stopping point, as it reports it. Its rule (a relative gain
    # under 1e-10, about 2e-7 per step here) leaves it below the fixed point by
    # the gain the refinement measures (7.8e-7, 3.7e-6 and 2.5e-6,
    # 2026-10-27); 1e-5 bounds that shortfall, as the plan for this suite
    # states (on the p09 probe it was 1.4e-5 at relative 8.8e-9).
    compare_values(sprintf("%s: LMest's reported log likelihood at its own stop (refining it gains %.2g)",
                           label, mode_ll - best$log_likelihood),
                   best$log_likelihood, fit$log_likelihood, tolerance = 1e-5))
  # Every table here is the recursion and the package computing the same
  # function at the same parameters, so they differ by rounding only: the
  # largest gap measured on any column was 6e-15 on the max(|x|, 1) scale
  # (2026-10-27). The harness default (1e-6) would let an absolute error of
  # 1e-4 through on a criterion near 1700, so the default here is 1e-10,
  # four orders above rounding; the columns the harness tolerates by name
  # (standard errors, thresholds, OCC) keep their tolerances.
  at_point <- .output_default_tolerances
  at_point[["default"]] <- 1e-10
  comparisons <- rbind(
    compare_tables(label, obtained, expected, at_point),
    .lta_mixed_parameter_distance(label, obtained, at_mode, errors),
    .lta_mixed_response_errors(label, obtained, errors),
    .with_output_columns(rbind(evaluator_check, claims)))
  fit_case <- harness_case(comparisons, output_ledger(label, fit_outputs(fit), comparisons))
  categories <- lapply(stats::setNames(nm = items), function(item) {
    as.character(sort(unique(frame[[item]])))
  })
  inference_case <- .lta_mixed_inference(label, fit, at_package, at_mode, errors, k1, k2,
                                         categories)
  case <- stack_cases(list(fit_case, inference_case))
  case$comparisons <- transform(
    case$comparisons, source = sprintf("LMest %s", utils::packageVersion("LMest")),
    dataset = fixture$name, precision = "machine")
  case
}

#' LMest's searches
#'
#' `lmestMixed()` with its deterministic start (`start = 0`) and two random
#' ones (`start = 1`, seeds 1 and 2), each stopped at a relative gain of 1e-10
#' or 5000 iterations. Its estimates are put in this suite's shapes, still in
#' LMest's own state and group-class order.
#'
#' Every sequence is given weight 1 (`weights`). Unweighted, LMest pools
#' sequences with the same responses into one row with a frequency, so its
#' group-class posteriors `W` come one per distinct sequence, in an order of
#' its own; with unit weights there is one row per sequence, in `id` order,
#' and the likelihood, the EM map and the estimates are the same (LMest's
#' `lk` is `sum(yv * log(fm))` either way).
#'
#' @param frame The simulated data, items coded from 0.
#' @param items Item names.
#' @param k1,k2 Group classes and states.
#' @param levels The fit's category labels, which name the columns.
#' @return A list of searches, each with `log_likelihood`, `n_parameters`,
#'   `parameters` and `group_posteriors` (LMest's `W`, sequences by group
#'   classes).
#' @noRd
.lta_mixed_lmest <- function(frame, items, k1, k2, levels) {
  formula <- stats::as.formula(paste(paste(items, collapse = " + "), "~ NULL"))
  starts <- list(c(start = 0, seed = 1), c(start = 1, seed = 1), c(start = 1, seed = 2))
  unit_weights <- rep(1, length(unique(frame$id)))
  lapply(starts, function(start) {
    # capture.output() keeps LMest's iteration table out of the harness log;
    # it is console text, not a condition.
    utils::capture.output(fitted <- LMest::lmestMixed(
      responsesFormula = formula, data = frame, index = c("id", "occ"), k1 = k1, k2 = k2,
      start = start[["start"]], weights = unit_weights, seed = start[["seed"]], tol = 1e-10,
      maxit = 5000))
    stopifnot("LMest must report one group-class posterior row per sequence" =
                identical(dim(fitted$W), c(length(unit_weights), as.integer(k1))))
    responses <- lapply(stats::setNames(seq_along(items), items), function(j) {
      block <- t(matrix(fitted$Psi[, , j], dim(fitted$Psi)[1L], k2))
      colnames(block) <- levels[[items[[j]]]]
      block
    })
    stopifnot("LMest's arrays must be states by group classes" =
                identical(dim(fitted$Piv), c(as.integer(k2), as.integer(k1))) &&
                  identical(dim(fitted$Pi), as.integer(c(k2, k2, k1))))
    list(log_likelihood = as.numeric(fitted$lk), n_parameters = fitted$np,
         parameters = list(la = as.vector(fitted$la),
                           initial = t(matrix(fitted$Piv, k2, k1)),
                           transition = array(as.vector(fitted$Pi), c(k2, k2, k1)),
                           responses = responses),
         group_posteriors = unname(as.matrix(fitted$W)))
  })
}

#' The package's estimates, in the reference's shapes
#'
#' Read from the fit's public fields to evaluate the reference recursion at
#' the package's point; no package function is called.
#' @noRd
.lta_mixed_fit_parameters <- function(fit, items) {
  k2 <- fit$n_profiles
  k1 <- fit$n_group_classes
  # The documented shapes. On the K1 != K2 fixture a field stored with its
  # axes swapped stops the case here, before it could be read into the wrong
  # cells; when K1 = K2 the likelihood at these estimates catches it.
  stopifnot("the fit's initial probabilities must be group classes by states" =
              identical(dim(fit$initial_probabilities), as.integer(c(k1, k2))),
            "the fit's transition probabilities must be from x to x group class" =
              identical(dim(fit$transition_probabilities), as.integer(c(k2, k2, k1))))
  responses <- lapply(stats::setNames(nm = items), function(item) {
    block <- unname(fit$response_probabilities[[item]])
    colnames(block) <- fit$categorical_levels[[item]]
    block
  })
  list(la = unname(fit$group_probabilities),
       initial = matrix(unname(fit$initial_probabilities), k1, k2),
       transition = array(as.vector(fit$transition_probabilities), c(k2, k2, k1)),
       responses = responses)
}

#' Put a reference solution's states and group classes in the fit's order
#'
#' States are matched on the response probabilities, then group classes on
#' their initial distributions and transition matrices, each by exhaustive
#' search; the fit's estimates are read for this matching only.
#' @noRd
.lta_mixed_align <- function(parameters, fit, items) {
  ours <- .lta_mixed_fit_parameters(fit, items)
  profile_signature <- function(p) do.call(cbind, unname(p$responses[items]))
  states <- align_classes(profile_signature(ours), profile_signature(parameters))
  parameters$responses <- lapply(parameters$responses, function(block) block[states, , drop = FALSE])
  parameters$initial <- parameters$initial[, states, drop = FALSE]
  parameters$transition <- parameters$transition[states, states, , drop = FALSE]
  class_signature <- function(p) {
    cbind(p$initial, t(vapply(seq_along(p$la), function(u) as.vector(p$transition[, , u]),
                              numeric(length(p$transition[, , 1L])))))
  }
  classes <- align_classes(class_signature(ours), class_signature(parameters))
  parameters$la <- parameters$la[classes]
  parameters$initial <- parameters$initial[classes, , drop = FALSE]
  parameters$transition <- parameters$transition[, , classes, drop = FALSE]
  parameters
}

# ---------------------------------------------------------------------------
# Expected tables
# ---------------------------------------------------------------------------

#' Every table the fit reports, from the recursion at the package's point
#'
#' The generic builder (`reference_tables()`) writes the tables a transition
#' fit shares with an ordinary mixture: as in the depmixS4 suite, the
#' per-occasion posteriors are the marginal smoothed state probabilities, each
#' sequence is a group and the sequence log-likelihood is the group
#' log-likelihood; here the group posteriors are the sequences' group-class
#' posteriors, and a transition fit charges each level its own entropy. The
#' transition tables are built here. The free-parameter columns are compared
#' with the fixed point instead ([.lta_mixed_parameter_distance()]) and the
#' response standard errors on the relative scale
#' ([.lta_mixed_response_errors()]), so they are left out here.
#'
#' @param evaluated `.lta_mixed_evaluate()` at the package's estimates.
#' @param parameters Those estimates.
#' @param n_parameters LMest's parameter count.
#' @return A named list of expected tables.
#' @noRd
.lta_mixed_expected <- function(evaluated, parameters, n_parameters, fit, frame,
                                categorical_frame, items) {
  k1 <- length(parameters$la)
  k2 <- ncol(parameters$initial)
  ids <- unique(frame$id)
  z <- evaluated$group_posteriors
  dimnames(z) <- list(as.character(ids), NULL)
  shared <- harness_reference(
    log_likelihood = evaluated$log_likelihood, n_parameters = n_parameters,
    posteriors = evaluated$posteriors, groups = frame$id,
    responses = parameters$responses,
    # Passed so that the model row counts the group classes; the
    # ordinary-mixture table it also builds is dropped below.
    profile_probabilities = parameters$initial, group_probabilities = parameters$la,
    group_posteriors = z,
    group_log_likelihood = stats::setNames(evaluated$sequence_log_likelihood,
                                           as.character(ids)),
    data = categorical_frame[items])
  tables <- reference_tables(shared, fit)
  tables$profile_probabilities <- NULL
  tables$assignments$occ <- frame$occ
  tables$responses$probability <- NULL
  # A row of transitions out of a is sum_b E[a -> b]: the recursion's two
  # routes to it (pair and single-occasion posteriors) must agree.
  stopifnot("expected moves out of an origin must be its occupancy before the last occasion" =
              isTRUE(all.equal(evaluated$origin_occupancy,
                               apply(evaluated$counts, c(1L, 3L), sum), tolerance = 1e-10)))
  modal_profile <- max.col(evaluated$posteriors, ties.method = "first")
  modal_class <- max.col(z, ties.method = "first")
  tables$transitions <- data.frame(
    group_class = rep(seq_len(k1), each = k2 * k2),
    from = rep(rep(seq_len(k2), each = k2), times = k1),
    to = rep(seq_len(k2), times = k2 * k1),
    expected_count = as.vector(aperm(evaluated$counts, c(2L, 1L, 3L))),
    stable = rep(as.vector(t(diag(TRUE, k2))), times = k1),
    # `?lta`: a row is estimated unless no sequence of its group class is
    # expected to occupy its origin before the final occasion, which leaves
    # the row without information. Every origin here is occupied (every
    # fixture is interior), so only the TRUE branch is exercised; an empty
    # origin lies on the boundary, where this suite's standard errors do not
    # exist.
    estimated = rep(as.vector(evaluated$origin_occupancy > 0), each = k2))
  tables$initial <- data.frame(group_class = rep(seq_len(k1), each = k2),
                               profile = rep(seq_len(k2), times = k1),
                               prevalence = as.vector(t(evaluated$prevalence)))
  # Every key: the rows-found check compares the modal state of each
  # occasion and the modal group class of each sequence.
  tables$sequences <- data.frame(group = frame$id,
                                 group_class = modal_class[match(frame$id, ids)],
                                 time = frame$occ, profile = modal_profile)
  # Balanced sequences: every one observed at every occasion, none with gaps.
  lengths <- as.vector(table(factor(frame$id, ids)))
  occasions <- length(unique(frame$occ))
  per_class <- lapply(seq_len(k1), function(u) {
    in_class <- lengths[modal_class == u]
    data.frame(group_class = u, groups = length(in_class), observations = sum(in_class),
               mean_length = mean(in_class), median_length = stats::median(in_class),
               shortest = min(in_class), longest = max(in_class),
               complete = sum(in_class == occasions), gaps = 0L)
  })
  tables$sequence_summary <- do.call(rbind, per_class)
  tables$sequence_lengths <- data.frame(group = ids, group_class = modal_class,
                                        occasions = lengths, observations = lengths,
                                        complete = lengths == occasions)
  tables
}

# ---------------------------------------------------------------------------
# The parameters, the standard errors and the search
# ---------------------------------------------------------------------------

#' A parameter's value and standard error for each row of a fit table
#'
#' @param table `"responses"`, `"transitions"` or `"initial"`.
#' @param frame The fit's table.
#' @param column The column compared.
#' @param parameters,errors Parameter sets in the reference's shapes.
#' @return A list with `value` and `error`, one per row of `frame`.
#' @noRd
.lta_mixed_cells <- function(table, frame, column, parameters, errors) {
  pick <- function(p) {
    if (identical(column, "group_class_probability")) return(p$la[frame$group_class])
    switch(table,
      responses = vapply(seq_len(nrow(frame)), function(row) {
        p$responses[[frame$indicator[[row]]]][frame$profile[[row]], frame$category[[row]]]
      }, numeric(1)),
      transitions = p$transition[cbind(frame$from, frame$to, frame$group_class)],
      initial = p$initial[cbind(frame$group_class, frame$profile)])
  }
  list(value = pick(parameters), error = pick(errors))
}

#' The free parameters against LMest's fixed point, in standard errors
#'
#' Each estimate's distance from the fixed point of LMest's solution, divided
#' by the reference standard error at the package's point, is held to 1e-3:
#' at a thousandth of a standard error no estimate, interval or test computed
#' from the fit can differ. Machine-precision agreement of the computation is
#' established separately, at the package's own point. (The package stops
#' within 6e-12 of the refined maximum in log-likelihood on all three
#' fixtures, and its estimates are within 2.7e-6 standard errors of it,
#' measured 2026-10-27.)
#'
#' @return Comparison rows with `output` set.
#' @noRd
.lta_mixed_parameter_distance <- function(label, obtained, at_mode, errors) {
  parameters <- list(responses = "probability",
                     transitions = c("probability", "group_class_probability"),
                     initial = c("probability", "group_class_probability"))
  rows <- lapply(names(parameters), function(table) {
    have <- obtained[[table]]
    lapply(parameters[[table]], function(column) {
      cells <- .lta_mixed_cells(table, have, column, at_mode, errors)
      output <- paste0(table, "$", column)
      compare_block(
        sprintf("%s: %s, distance from LMest's fixed point in standard errors", label, output),
        rep(0, nrow(have)), abs(have[[column]] - cells$value) / cells$error,
        tolerance = 1e-3, output = output)
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

#' The response standard errors against the reference, relative
#'
#' Both are the inverse observed information: the package's from its analytic
#' derivatives, the reference from a Richardson-extrapolated numerical Hessian
#' of the recursion. They agree to 1e-7 relative here (measured 2026-10-27);
#' 1e-4 relative is the harness's standard-error tolerance, stated relative
#' because these errors are near 0.02, where the absolute default would allow
#' 0.5 percent.
#' @noRd
.lta_mixed_response_errors <- function(label, obtained, errors) {
  have <- obtained$responses
  reference <- vapply(seq_len(nrow(have)), function(row) {
    errors$responses[[have$indicator[[row]]]][have$profile[[row]], have$category[[row]]]
  }, numeric(1))
  compare_block(
    sprintf("%s: responses$probability_standard_error, relative to the numDeriv reference", label),
    rep(0, nrow(have)), abs(have$probability_standard_error / reference - 1),
    tolerance = 1e-4, output = "responses$probability_standard_error")
}

#' The package's search against LMest's best
#'
#' @param lmest_ll The log-likelihoods LMest's searches report.
#' @param mode_ll The best of them, continued to its fixed point.
#' @param cycles SQUAREM cycles that took.
#' @param fit_ll The package's maximum.
#' @return Comparison rows, none an output of the fit.
#' @noRd
.lta_mixed_search <- function(label, lmest_ll, mode_ll, cycles, fit_ll) {
  reached <- sum(lmest_ll >= max(lmest_ll) - 1e-4)
  rows <- list(compare_values(
    sprintf(paste("%s: shortfall of the package's own search against LMest's best",
                  "(best of %d LMest searches, reached by %d; refined to its fixed point",
                  "in %d SQUAREM cycles)"), label, length(lmest_ll), reached, cycles),
    0, max(mode_ll - fit_ll, 0), tolerance = 1e-6))
  excess <- fit_ll - mode_ll
  if (excess > 1e-6) {
    # A higher maximum than any of LMest's searches: a finding about LMest's
    # search, not a package defect.
    rows <- c(rows, list(compare_values(
      sprintf("%s: package maximum above LMest's best (by %.3g; LMest's search)",
              label, excess), excess, excess, tolerance = 1e-6)))
  }
  do.call(rbind, rows)
}

#' The rows `parameter_inference()` must report, from the model's dimensions
#'
#' Written from the model, not read from the package's table, so a row the
#' package drops is a row the comparison misses: K1 group-class
#' probabilities, K1 K2 initial probabilities, K1 K2^2 transition
#' probabilities and K2 sum_j C_j response probabilities (26 rows for two
#' group classes, two states and three binary items; 32 for three-category
#' items; 53 for three states). The key columns spell each row as
#' `?parameter_inference` documents it: the outcome is the state (or the
#' group class, for its probability), the term the group class, the origin
#' (`group_class_u:profile_a`) or the item and category (`item:category`).
#' The group-class probabilities have no term; see [.lta_mixed_inference()]
#' for its spelling. The integer columns locate each row's parameter.
#'
#' @param k1,k2 Group classes and states, from the fixture.
#' @param categories Named list: each item's category codes, from the data.
#' @return A data frame with `level`, `outcome`, `term`, `parameter` and the
#'   locating columns `item`, `category`, `group_class`, `from` and `state`.
#' @noRd
.lta_mixed_inference_grid <- function(k1, k2, categories) {
  cells <- function(level, outcome, term, parameter, item = NA_character_,
                    category = NA_character_, group_class = NA_integer_, from = NA_integer_,
                    state = NA_integer_) {
    data.frame(level = level, outcome = outcome, term = term, parameter = parameter,
               item = item, category = category, group_class = group_class, from = from,
               state = state, stringsAsFactors = FALSE)
  }
  profile <- function(v) sprintf("profile_%d", v)
  group_class <- function(u) sprintf("group_class_%d", u)
  responses <- lapply(names(categories), function(item) {
    grid <- expand.grid(category = categories[[item]], state = seq_len(k2),
                        stringsAsFactors = FALSE)
    cells("measurement", profile(grid$state), paste0(item, ":", grid$category), "response",
          item = item, category = grid$category, state = grid$state)
  })
  initial <- expand.grid(state = seq_len(k2), group_class = seq_len(k1))
  transition <- expand.grid(state = seq_len(k2), from = seq_len(k2), group_class = seq_len(k1))
  grid <- do.call(rbind, c(responses, list(
    cells("group", group_class(seq_len(k1)), "(none)", "probability",
          group_class = seq_len(k1)),
    cells("profile", profile(initial$state), group_class(initial$group_class),
          "initial_probability", group_class = initial$group_class, state = initial$state),
    cells("profile", profile(transition$state),
          paste0(group_class(transition$group_class), ":", profile(transition$from)),
          "transition_probability", group_class = transition$group_class,
          from = transition$from, state = transition$state))))
  stopifnot("the grid must hold one row per probability of the model" =
              nrow(grid) == k1 + k1 * k2 + k1 * k2^2 + k2 * sum(lengths(categories)))
  grid
}

#' One parameter set's value at each row of the grid
#' @param grid From [.lta_mixed_inference_grid()].
#' @param p Parameters (or standard errors) in the reference's shapes.
#' @return A numeric vector, one per row of `grid`.
#' @noRd
.lta_mixed_grid_values <- function(grid, p) {
  vapply(seq_len(nrow(grid)), function(row) {
    at <- grid[row, ]
    switch(at$parameter,
      response = p$responses[[at$item]][at$state, at$category],
      probability = p$la[[at$group_class]],
      initial_probability = p$initial[at$group_class, at$state],
      transition_probability = p$transition[at$from, at$state, at$group_class])
  }, numeric(1))
}

#' `parameter_inference()` on every column, through the verb ledger
#'
#' One row per probability of the model, the rows generated from its
#' dimensions ([.lta_mixed_inference_grid()]), so a missing row fails the
#' rows-found check and every column comparison: `estimate` against the
#' fixed point in standard errors (as the fit's tables are),
#' `standard_error` relative to the numDeriv reference, and the rest at the
#' package's point. The documented convention (`?parameter_inference`): a
#' probability's null value lies on the boundary of its space, so
#' `statistic`, `p_value` and `p_adjusted` are `NA`, and its Wald interval is
#' formed on the logit scale, `plogis(qlogis(p) -/+ z se / (p (1 - p)))` with
#' `z = qnorm(0.975)`; the reference interval takes the reference standard
#' error and the package's estimate read from the fit's fields.
#'
#' The group-class probabilities have no `term`. The rows-found check of
#' `compare_tables()` counts a row with a missing key as absent, so `term` is
#' spelled `"(none)"` on both sides; it is a key, never a compared value.
#'
#' @param k1,k2,categories As for [.lta_mixed_inference_grid()].
#' @return A `harness_case`.
#' @noRd
.lta_mixed_inference <- function(label, fit, at_package, at_mode, errors, k1, k2, categories) {
  table <- parameter_inference(fit)
  table$term[is.na(table$term)] <- "(none)"
  keys <- c("level", "outcome", "term", "parameter")
  obtained <- list(parameter_inference = keyed_table(table, keys))
  have <- obtained$parameter_inference
  grid <- .lta_mixed_inference_grid(k1, k2, categories)
  # The package's row for each grid row, NA where it reports none; a missing
  # row then compares as a disagreement in every block below.
  row_key <- function(frame) do.call(paste, c(frame[keys], sep = "\r"))
  matched <- match(row_key(grid), row_key(have))
  estimate <- have$estimate[matched]
  standard_error <- have$standard_error[matched]
  at_estimates <- .lta_mixed_grid_values(grid, at_package)
  reference_error <- .lta_mixed_grid_values(grid, errors)
  z <- stats::qnorm(0.975)
  half_width <- z * reference_error / (at_estimates * (1 - at_estimates))
  expected <- list(parameter_inference = data.frame(
    grid[keys], statistic = NA_real_, p_value = NA_real_, p_adjusted = NA_real_,
    conf_low = stats::plogis(stats::qlogis(at_estimates) - half_width),
    conf_high = stats::plogis(stats::qlogis(at_estimates) + half_width)))
  extra <- rbind(
    # Every row of the model is reported, and nothing else.
    compare_values(sprintf("%s: parameter_inference() rows, against the model's count", label),
                   nrow(grid), nrow(have), tolerance = 0.5),
    # The table must report the fit's own estimates in the right cells: the
    # same numbers, so 1e-10 as at the package's point. This is what catches
    # a value put in the wrong row (a group-class axis read as a state axis).
    compare_block(sprintf("%s: parameter_inference$estimate, against the fit's own estimates",
                          label),
                  at_estimates, estimate, tolerance = 1e-10,
                  output = "parameter_inference$estimate"),
    compare_block(sprintf("%s: parameter_inference$estimate, distance from LMest's fixed point in standard errors",
                          label),
                  rep(0, nrow(grid)),
                  abs(estimate - .lta_mixed_grid_values(grid, at_mode)) / reference_error,
                  tolerance = 1e-3, output = "parameter_inference$estimate"),
    compare_block(sprintf("%s: parameter_inference$standard_error, relative to the numDeriv reference",
                          label),
                  rep(0, nrow(grid)), abs(standard_error / reference_error - 1),
                  tolerance = 1e-4, output = "parameter_inference$standard_error"))
  compare_verb(sprintf("%s: parameter_inference()", label), obtained, expected, extra = extra)
}
