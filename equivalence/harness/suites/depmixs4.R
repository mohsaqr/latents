# Latent transition analysis against depmixS4 (Visser & Speekenbrink 2010,
# Journal of Statistical Software 36(7):1-21).
#
# depmixS4 fits hidden Markov and latent Markov models by expectation
# maximization over the same likelihood: an initial state distribution, a
# homogeneous first-order transition matrix, and state-conditional response
# distributions shared across occasions. Its `ntimes` argument carries the
# sequence lengths that multilpa reads from `time =` and the group column. The
# parameter blocks map one to one:
#
#   depmixS4 prior      <-> multilpa initial_probabilities
#   depmixS4 transition <-> multilpa transition_probabilities  (rows: from)
#   depmixS4 response   <-> multilpa means/variances or response_probabilities
#
# Two different claims are checked for every case, because they can fail
# separately. Setting depmixS4's parameters to multilpa's estimates and reading
# its log likelihood compares the two likelihood *functions* and nothing else.
# Letting depmixS4 search from its own start compares the two *optimizers*, and
# only then does agreement say that both found the same maximum.
#
# Beyond those two claims every table the fit reports is rebuilt from
# depmixS4's solution and compared cell by cell (see outputs.R):
#
#   * Tables an ordinary mixture also reports (profiles, covariances,
#     responses, posteriors, group posteriors, assignments, counts, the
#     classification tables, BCH weights, entropy, information criteria, the
#     model row, residuals, Gaussian diagnostics) come from the generic
#     `reference_tables()`. A transition fit defines them on the same objects
#     -- the per-occasion posteriors are the smoothed marginal state
#     probabilities (`subject_posteriors`, R/transitions.R), each sequence is a
#     group with one group class, and the sequence log-likelihood is the group
#     log-likelihood -- so the builder is fed depmixS4's smoothed posteriors
#     (`forwardbackward()$gamma`) and its per-sequence log-likelihoods (the
#     scaling constants of its forward pass) and is correct as written. Only
#     `profile_probabilities`, a table a transition fit does not have, is
#     dropped.
#   * The transition tables are built here: `transitions` from depmixS4's
#     transition matrix and its two-slice posteriors (`xi`), `initial` from its
#     prior and the mean smoothed posterior, `sequences` from the modal
#     smoothed state of every occasion, and the three sequence descriptions
#     from the `ntimes` and occasions depmixS4 was given.
#   * Standard errors: the reference the ledger counts is the
#     Richardson-extrapolated Hessian (numDeriv) of a forward-algorithm
#     likelihood written here in base R, at depmixS4's estimates, as the mclust
#     suite does for a Gaussian mixture. depmixS4's own estimator,
#     `standardError()`, targets the same quantity (observed information under
#     its sum-to-one constraints) but differences with `nlme::fdHess()`, which
#     is 1 to 15 percent off on gaussian-3-profile and 4-fold off for a mean
#     near zero (see `.depmixs4_own_errors()`); it is compared with fdHess
#     replaced by numDeriv and everything else depmixS4's, as rows outside the
#     ledger (`output = NA`).
#   * depmixS4's solution is first refined with its own EM to the fixed point
#     (see `.depmixs4_refine()`); the second claim is checked on the unrefined
#     fit, which is what its own search found.

suite_depmixs4 <- function() {
  require_suite_packages(
    "depmixS4", "numDeriv",
    reason = paste("the independent hidden-Markov reference this suite compares",
                   "transitions against, and the reference observed information"))
  cases <- list(
    list(dataset = "gaussian-2-profile", profiles = 2L, categorical = FALSE,
         groups = 80L, occasions = 6L, seed = 101L),
    list(dataset = "gaussian-3-profile", profiles = 3L, categorical = FALSE,
         groups = 90L, occasions = 7L, seed = 102L),
    list(dataset = "gaussian-unbalanced", profiles = 2L, categorical = FALSE,
         groups = 70L, occasions = 6L, seed = 103L, drop = 0.2),
    list(dataset = "mixed-2-profile", profiles = 2L, categorical = TRUE,
         groups = 80L, occasions = 6L, seed = 104L)
  )
  stack_cases(lapply(cases, function(case) do.call(.depmixs4_case, case)))
}

#' Simulate one latent transition dataset
#' @param profiles Number of states.
#' @param groups Number of sequences.
#' @param occasions Sequence length before any thinning.
#' @param seed Simulation seed.
#' @param categorical Whether to add a three-category indicator.
#' @param drop Proportion of interior observations to delete, making the
#'   sequences unequal in length.
#' @return A data frame with `g`, `t`, two Gaussian indicators and optionally a
#'   categorical one, ordered by group and occasion.
#' @noRd
.depmixs4_data <- function(profiles, groups, occasions, seed, categorical,
                           drop = 0) {
  set.seed(seed)
  stay <- 0.75
  transition <- matrix((1 - stay) / (profiles - 1), profiles, profiles)
  diag(transition) <- stay
  frame <- do.call(rbind, lapply(seq_len(groups), function(group) {
    path <- unlist(Reduce(function(previous, step) {
      sample.int(profiles, 1L, prob = transition[previous, ])
    }, seq_len(occasions - 1L), init = sample.int(profiles, 1L),
    accumulate = TRUE))
    data.frame(g = group, t = seq_len(occasions), state = path)
  }))
  centres <- seq(-2, 2, length.out = profiles)
  frame$y1 <- stats::rnorm(nrow(frame), centres[frame$state], 0.8)
  frame$y2 <- stats::rnorm(nrow(frame), 0.75 * centres[frame$state], 0.9)
  if (categorical) {
    shares <- c(0.7, 0.2, 0.1)
    frame$c1 <- factor(vapply(frame$state, function(state) {
      sample(c("a", "b", "c"), 1L,
             prob = shares[1L + (seq_len(3L) + state - 2L) %% 3L])
    }, character(1)), levels = c("a", "b", "c"))
  }
  if (drop > 0) {
    interior <- which(frame$t > 1L & frame$t < occasions)
    frame <- frame[-sample(interior, floor(drop * length(interior))), ]
  }
  frame[order(frame$g, frame$t), ]
}

#' Assemble depmixS4's parameter vector from a multilpa transition fit
#' @param fit A fitted `multilpa_transitions` model with one group class.
#' @param categorical Whether a categorical indicator is present.
#' @return A numeric vector in depmixS4's `getpars()` order.
#' @noRd
.depmixs4_pars <- function(fit, categorical) {
  responses <- unlist(lapply(seq_len(fit$n_profiles), function(profile) {
    gaussian <- as.vector(rbind(fit$means[profile, ],
                                sqrt(fit$variances[profile, ])))
    if (!categorical) return(gaussian)
    c(gaussian, fit$response_probabilities[["c1"]][profile, ])
  }), use.names = FALSE)
  c(fit$initial_probabilities[1L, ],
    as.vector(t(matrix(fit$transition_probabilities[, , 1L],
                       fit$n_profiles, fit$n_profiles))),
    responses)
}

#' Compare one dataset against depmixS4
#' @param dataset Label carried into the report.
#' @param profiles Number of states.
#' @param categorical Whether to include a categorical indicator.
#' @param groups Number of sequences.
#' @param occasions Sequence length.
#' @param seed Simulation seed.
#' @param drop Proportion of interior observations deleted.
#' @return A `harness_case`.
#' @noRd
.depmixs4_case <- function(dataset, profiles, categorical, groups, occasions,
                           seed, drop = 0) {
  frame <- .depmixs4_data(profiles, groups, occasions, seed, categorical, drop)
  indicators <- c("y1", "y2", if (categorical) "c1")
  # Warnings are not suppressed here. The harness counts every warning a suite
  # raises and records the first one in suites.csv, so a start that hit its
  # iteration cap is visible in the evidence rather than hidden by the call
  # that produced it.
  fit <- lta(
    frame, indicators, "g", n_profiles = profiles, time = "t",
    categorical = if (categorical) "c1" else character(),
    n_starts = 8, seed = seed, max_iter = 2000, tol = 1e-12)
  formulas <- lapply(indicators, function(name) stats::as.formula(paste(name, "~ 1")))
  families <- c(list(stats::gaussian(), stats::gaussian()),
                if (categorical) list(depmixS4::multinomial("identity")))
  ntimes <- as.integer(table(frame$g)[unique(as.character(frame$g))])
  model <- depmixS4::depmix(formulas, data = frame, nstates = profiles,
                            family = families, ntimes = ntimes)
  label <- sprintf("%s, %d states", dataset, profiles)

  # Claim one: the same likelihood function. depmixS4 evaluated exactly at
  # multilpa's estimates must return multilpa's log likelihood.
  at_ours <- as.numeric(depmixS4::logLik(
    depmixS4::setpars(model, .depmixs4_pars(fit, categorical))))
  same_function <- compare_values(
    quantity = sprintf("%s: log likelihood at multilpa's estimates", label),
    reference = at_ours, obtained = fit$log_likelihood,
    tolerance = 1e-8, scale = "relative")

  # Claim two: the same maximum, found independently. depmixS4 searches from
  # its own start; agreement then says both optimizers reached the same mode.
  set.seed(seed)
  theirs <- depmixS4::fit(model, verbose = FALSE, emcontrol =
    depmixS4::em.control(maxit = 2000, tol = 1e-12))
  same_maximum <- compare_values(
    quantity = sprintf("%s: maximized log likelihood", label),
    reference = as.numeric(depmixS4::logLik(theirs)),
    obtained = fit$log_likelihood, tolerance = 1e-6, scale = "relative")
  claims <- rbind(same_function, same_maximum)

  shortfall <- as.numeric(depmixS4::logLik(theirs)) - fit$log_likelihood
  case <- if (abs(shortfall) > 1e-4) {
    # Different modes: depmixS4's estimates describe another solution, so only
    # the two claims are comparable, and the ledger says why nothing else is.
    inventory <- fit_outputs(fit)
    harness_case(claims, output_ledger(label, inventory, claims, exempt = stats::setNames(
      rep("depmixS4 stopped at a different mode, so its estimates describe another solution",
          length(unique(inventory$table))), unique(inventory$table))))
  } else {
    # depmixS4 stops on a relative change of 1e-12 in the log-likelihood, about
    # 2e-9 in absolute terms here. Where the likelihood is flat that leaves its
    # estimates short of EM's fixed point: on gaussian-3-profile the middle
    # state's y1 mean stopped at -0.02302896, while 33 more steps of depmixS4's
    # own EM move it to -0.02303057, beside multilpa's -0.02303045, for a gain
    # of 7e-10 in the log-likelihood. The solution is refined with
    # depmixS4's own EM until its parameters stop moving, as the mclust suite
    # refines mclust's, so the reference is depmixS4's estimate and not where
    # its stopping rule fired.
    refined <- .depmixs4_refine(theirs)
    reference <- .depmixs4_reference(refined, fit, frame, indicators, ntimes, categorical)
    obtained <- fit_tables(fit)
    obtained$model$.row <- 1L
    comparisons <- compare_tables(label, obtained,
                                  .depmixs4_expected(reference, fit, frame, ntimes),
                                  .output_default_tolerances)
    comparisons <- rbind(comparisons, claims,
                         .depmixs4_own_errors(label, refined, reference, obtained))
    harness_case(comparisons, output_ledger(label, fit_outputs(fit), comparisons))
  }
  case$comparisons <- transform(
    case$comparisons, source = sprintf("depmixS4 %s", utils::packageVersion("depmixS4")),
    dataset = dataset, precision = "machine")
  case
}

#' Continue depmixS4's EM from its own solution to the fixed point
#'
#' `fit()` with `maxit = 0` and no random start performs exactly one EM step
#' from the model's current parameters (one M-step on the current posteriors,
#' then the E-step), and with no likelihood comparison left to make it raises
#' no warning. Steps are taken until no parameter moves by more than 1e-12;
#' the stopping rule is on the parameters, because on a flat likelihood the
#' log-likelihood stops changing (at 1e-12) well before they do.
#'
#' @param theirs A fitted depmixS4 model.
#' @param max_steps Bound on the EM steps.
#' @return The refined fitted model. A classed warning is raised if the bound
#'   is reached, and refining must not lose likelihood.
#' @noRd
.depmixs4_refine <- function(theirs, max_steps = 5000L) {
  step <- function(object) {
    depmixS4::fit(object, verbose = FALSE, emcontrol = depmixS4::em.control(
      maxit = 0, tol = 0, crit = "absolute", random.start = FALSE))
  }
  # Reduce() rather than a loop with a break: once the parameters have stopped
  # moving each remaining element returns the state untouched.
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    moved <- step(state$model)
    change <- max(abs(depmixS4::getpars(moved) - depmixS4::getpars(state$model)))
    list(model = moved, settled = change < 1e-12, steps = index)
  }, seq_len(max_steps), init = list(model = theirs, settled = FALSE, steps = 0L))
  if (!state$settled) {
    warning(warningCondition(sprintf(
      "depmixS4's EM was still moving after %d refining steps", max_steps),
      class = "multilpa_harness_no_converge"))
  }
  stopifnot("refining depmixS4's solution must not lose likelihood" =
              as.numeric(depmixS4::logLik(state$model)) >=
                as.numeric(depmixS4::logLik(theirs)) - 1e-9)
  state$model
}

#' depmixS4's solution, in the fit's state order
#'
#' depmixS4 lays its parameters out as the prior, then the transition matrix
#' by row (from), then one block per state holding, for each response in
#' formula order, the Gaussian intercept and SD or the category
#' probabilities. Its states are put in the fit's order by matching the
#' Gaussian means (labels are arbitrary in both programs).
#'
#' `forwardbackward()` returns `xi[t, to, from]` = P(state at t is `from`,
#' state at t + 1 is `to`) (its C code fills `rep(alpha, each = ns) * (B *
#' beta * A)`), with `NA` at the last occasion of every sequence; the two
#' marginal identities are asserted before the counts are used, so a change
#' of that layout fails here rather than producing wrong counts.
#'
#' @return A list of aligned estimates, posteriors, counts and standard errors.
#' @noRd
.depmixs4_reference <- function(theirs, fit, frame, indicators, ntimes, categorical) {
  k <- fit$n_profiles
  continuous <- c("y1", "y2")
  pars <- depmixS4::getpars(theirs)
  width <- (length(pars) - k - k^2) / k
  stopifnot("depmixS4 parameter blocks must be equally sized" = width == floor(width))
  blocks <- matrix(pars[k + k^2 + seq_len(k * width)], k, width, byrow = TRUE)
  gaussian_columns <- seq_len(2L * length(continuous))
  means <- blocks[, gaussian_columns[c(TRUE, FALSE)], drop = FALSE]
  sds <- blocks[, gaussian_columns[c(FALSE, TRUE)], drop = FALSE]
  dimnames(means) <- dimnames(sds) <- list(NULL, continuous)
  order_states <- align_classes(fit$means[, continuous, drop = FALSE], means)

  responses <- if (categorical) {
    levels_c1 <- fit$categorical_levels[["c1"]]
    probabilities <- blocks[, length(gaussian_columns) + seq_along(levels_c1), drop = FALSE]
    dimnames(probabilities) <- list(NULL, levels_c1)
    list(c1 = probabilities[order_states, , drop = FALSE])
  }
  initial <- pars[seq_len(k)][order_states]
  transition <- matrix(pars[k + seq_len(k^2)], k, k, byrow = TRUE)[order_states, order_states]

  passes <- depmixS4::forwardbackward(theirs)
  last <- cumsum(ntimes)
  first <- last - ntimes + 1L
  inner <- setdiff(seq_len(nrow(frame)), last)
  pairs <- aperm(passes$xi[inner, , , drop = FALSE], c(1L, 3L, 2L))  # [t, from, to]
  stopifnot(
    "depmixS4's xi must sum over `to` to the smoothed posterior at t" =
      max(abs(apply(pairs, c(1L, 2L), sum) - passes$gamma[inner, , drop = FALSE])) < 1e-10,
    "depmixS4's xi must sum over `from` to the smoothed posterior at t + 1" =
      max(abs(apply(pairs, c(1L, 3L), sum) - passes$gamma[inner + 1L, , drop = FALSE])) < 1e-10)
  counts <- apply(pairs, c(2L, 3L), sum)[order_states, order_states]
  # The forward pass rescales alpha at every occasion by `sca`, so a
  # sequence's log-likelihood is minus the sum of the logs of its scalings.
  sequence_log_likelihood <- vapply(seq_along(ntimes), function(i) {
    -sum(log(passes$sca[first[[i]]:last[[i]]]))
  }, numeric(1))
  stopifnot("per-sequence log-likelihoods must add to depmixS4's total" =
              abs(sum(sequence_log_likelihood) -
                    as.numeric(depmixS4::logLik(theirs))) < 1e-8)

  estimates <- list(initial = initial, transition = transition,
                    means = means[order_states, , drop = FALSE],
                    variances = sds[order_states, , drop = FALSE]^2,
                    responses = responses)
  list(estimates = estimates,
       posteriors = passes$gamma[, order_states, drop = FALSE],
       counts = counts, sequence_log_likelihood = sequence_log_likelihood,
       log_likelihood = as.numeric(depmixS4::logLik(theirs)),
       n_parameters = depmixS4::freepars(theirs),
       order_states = order_states,
       standard_errors = .depmixs4_standard_errors(estimates, frame, ntimes,
                                                   as.numeric(depmixS4::logLik(theirs))))
}

#' Every table the fit reports, rebuilt from depmixS4's solution
#'
#' @param reference From [.depmixs4_reference()].
#' @param fit The fitted model, for its labels and dimensions only.
#' @param frame The data both programs were fitted to.
#' @param ntimes The sequence lengths depmixS4 was given.
#' @return A named list of expected tables.
#' @noRd
.depmixs4_expected <- function(reference, fit, frame, ntimes) {
  k <- fit$n_profiles
  estimates <- reference$estimates
  groups <- unique(frame$g)
  shared <- harness_reference(
    log_likelihood = reference$log_likelihood, n_parameters = reference$n_parameters,
    posteriors = reference$posteriors, groups = frame$g,
    means = estimates$means, variances = estimates$variances,
    responses = estimates$responses, standard_errors = reference$standard_errors,
    group_probabilities = 1,
    group_posteriors = matrix(1, length(groups), 1L,
                              dimnames = list(as.character(groups), NULL)),
    group_log_likelihood = stats::setNames(reference$sequence_log_likelihood,
                                           as.character(groups)),
    data = frame[fit$vars])
  tables <- reference_tables(shared, fit)
  # The ordinary-mixture profile-share table; a transition fit reports the
  # initial distribution and the prevalence instead, compared below.
  tables$profile_probabilities <- NULL
  # The occasion of each row is part of the assignments a transition fit
  # reports, and must be the occasion depmixS4 placed that row at.
  tables$assignments$t <- frame$t

  # Modal smoothed state of every occasion (first state on a tie), the
  # package's `sequences` profile; it is a key of that table, so it is
  # compared through the rows-found check.
  modal <- max.col(reference$posteriors, ties.method = "first")
  stable <- as.vector(t(diag(TRUE, k)))
  tables$transitions <- data.frame(
    group_class = 1L, from = rep(seq_len(k), each = k), to = rep(seq_len(k), times = k),
    probability = as.vector(t(estimates$transition)),
    expected_count = as.vector(t(reference$counts)),
    stable = stable,
    # A transition row is estimable when its origin is expected to be left at
    # least once in some sequence.
    estimated = rep(rowSums(reference$counts) > 0, each = k),
    group_class_probability = 1)
  tables$initial <- data.frame(
    group_class = 1L, profile = seq_len(k), probability = estimates$initial,
    # Prevalence: the share of all occasions expected in each state.
    prevalence = colMeans(reference$posteriors), group_class_probability = 1)
  tables$sequences <- data.frame(group = frame$g, group_class = 1L, time = frame$t,
                                 profile = modal)

  # The sequence descriptions follow from the data depmixS4 was given: its
  # `ntimes`, and the occasions observed anywhere. On the observed-occasion
  # layout a sequence occupies as many occasions as it has observations, and
  # is complete when that is the longest; it has a gap when its occasions
  # skip a position seen elsewhere in the data.
  positions <- sort(unique(frame$t))
  position <- match(frame$t, positions)
  span <- tapply(position, frame$g, function(p) max(p) - min(p) + 1L)[as.character(groups)]
  # One row, which compare_tables() matches by position, so the group-class
  # key is left out: it would otherwise be compared as a value.
  tables$sequence_summary <- data.frame(
    groups = length(ntimes), observations = sum(ntimes),
    mean_length = mean(ntimes), median_length = stats::median(ntimes),
    shortest = min(ntimes), longest = max(ntimes),
    complete = sum(ntimes == length(positions)), gaps = sum(ntimes < span))
  tables$sequence_lengths <- data.frame(group = groups, group_class = 1L,
                                        occasions = ntimes, observations = ntimes,
                                        complete = ntimes == max(ntimes))
  tables
}

#' Observed-information standard errors of a latent transition model, from scratch
#'
#' The log-likelihood is a scaled forward recursion written here in base R,
#' over every sequence at once, parameterised by the free natural parameters:
#' the first `K - 1` initial probabilities, the first `K - 1` entries of each
#' transition row, the Gaussian means and variances, and the first `C - 1`
#' category probabilities of each state. Its Hessian at depmixS4's estimates is
#' taken by Richardson extrapolation (numDeriv). At a maximum the inverse
#' observed information transforms with the Jacobian of any
#' reparameterisation, so these are the standard errors of the natural
#' parameters whatever scale either program works on; the last category's
#' error is that of one minus the sum of the free ones.
#'
#' @param estimates Aligned estimates from [.depmixs4_reference()].
#' @param frame The data.
#' @param ntimes Sequence lengths, in data order.
#' @param log_likelihood depmixS4's maximum, which the recursion must
#'   reproduce at its estimates before its Hessian is trusted.
#' @return A list with `mean`, `variance` and, with a categorical indicator,
#'   `probability`, as `harness_reference()` takes them.
#' @noRd
.depmixs4_standard_errors <- function(estimates, frame, ntimes, log_likelihood) {
  k <- length(estimates$initial)
  continuous <- colnames(estimates$means)
  j <- length(continuous)
  y <- as.matrix(frame[continuous])
  codes <- if (!is.null(estimates$responses)) as.integer(frame$c1)
  n_categories <- if (is.null(codes)) 0L else ncol(estimates$responses$c1)
  # Groups by occasions, holding the row of each observation; rows are already
  # in group and occasion order, so the occasions are each sequence's own.
  last <- cumsum(ntimes)
  slot <- t(vapply(seq_along(ntimes), function(i) {
    c(last[[i]] - ntimes[[i]] + seq_len(ntimes[[i]]),
      rep(NA_integer_, max(ntimes) - ntimes[[i]]))
  }, integer(max(ntimes))))
  slot <- matrix(slot, length(ntimes))

  sizes <- c(initial = k - 1L, transition = k * (k - 1L), mean = k * j,
             variance = k * j, probability = k * max(n_categories - 1L, 0L))
  ends <- cumsum(sizes)
  part <- function(theta, name) theta[ends[[name]] - sizes[[name]] + seq_len(sizes[[name]])]
  complete <- function(free, width) {
    free <- matrix(free, ncol = width - 1L)
    cbind(free, 1 - rowSums(free))
  }
  unpack <- function(theta) {
    list(initial = as.vector(complete(part(theta, "initial"), k)),
         transition = complete(matrix(part(theta, "transition"), k, byrow = TRUE), k),
         means = matrix(part(theta, "mean"), k, j),
         variances = matrix(part(theta, "variance"), k, j),
         probability = if (n_categories > 0L)
           complete(matrix(part(theta, "probability"), k), n_categories))
  }
  log_likelihood_at <- function(theta) {
    p <- unpack(theta)
    log_emission <- vapply(seq_len(k), function(state) {
      gaussian <- rowSums(vapply(seq_len(j), function(column) {
        stats::dnorm(y[, column], p$means[state, column], sqrt(p$variances[state, column]),
                     log = TRUE)
      }, numeric(nrow(y))))
      if (is.null(codes)) gaussian else gaussian + log(p$probability[state, codes])
    }, numeric(nrow(y)))
    log_emission <- matrix(log_emission, nrow(y))
    emit <- function(rows) {
      block <- log_emission[rows, , drop = FALSE]
      peak <- apply(block, 1L, max)
      list(density = exp(block - peak), peak = peak)
    }
    start <- emit(slot[, 1L])
    alpha <- sweep(start$density, 2L, p$initial, "*")
    scale <- rowSums(alpha)
    state <- list(alpha = alpha / scale, total = start$peak + log(scale))
    # Each step updates only the sequences that reach this occasion; the
    # recursion runs over occasions, the short dimension, with all sequences
    # as one matrix product.
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
  theta <- c(estimates$initial[-k], as.vector(t(estimates$transition[, -k, drop = FALSE])),
             as.vector(estimates$means), as.vector(estimates$variances),
             if (n_categories > 0L) as.vector(estimates$responses$c1[, -n_categories]))
  stopifnot("the base-R forward recursion must reproduce depmixS4's maximum" =
              abs(log_likelihood_at(theta) - log_likelihood) < 1e-8 * abs(log_likelihood))
  covariance <- solve(-numDeriv::hessian(log_likelihood_at, theta))
  at <- function(name) ends[[name]] - sizes[[name]] + seq_len(sizes[[name]])
  errors <- sqrt(diag(covariance))
  out <- list(mean = matrix(errors[at("mean")], k, j, dimnames = list(NULL, continuous)),
              variance = matrix(errors[at("variance")], k, j,
                                dimnames = list(NULL, continuous)))
  if (n_categories > 0L) {
    free <- matrix(at("probability"), k)
    out$probability <- list(c1 = t(vapply(seq_len(k), function(state) {
      block <- covariance[free[state, ], free[state, ], drop = FALSE]
      c(sqrt(diag(block)), sqrt(sum(block)))
    }, numeric(n_categories))))
    dimnames(out$probability$c1) <- list(NULL, colnames(estimates$responses$c1))
  }
  out
}

#' depmixS4's own standard-error estimator against the package's
#'
#' depmixS4's `standardError()` is `vcov()`: `hessian2vcov()` (the bordered
#' inverse that imposes its linear sum-to-one constraints, `conMat`) applied to
#' `hessian()`, which is `nlme::fdHess()` of minus its log-likelihood in its
#' own parameters (probabilities, Gaussian intercepts and SDs). fdHess takes a
#' step proportional to each parameter (relative step `eps^(1/3)`,
#' `minAbsPar = 0`), so a parameter near zero is differenced over a step of
#' about 1e-7 and the second difference is dominated by rounding: on
#' gaussian-3-profile its SE of the middle state's y1 mean (estimate -0.023)
#' is 0.287 where the observed information gives 0.0720, and its other SEs
#' there are 1 to 15 percent off. That is a limitation of the finite
#' difference, not a different estimand, so the comparison here keeps
#' depmixS4's likelihood (`logLik()` after `setpars()`), its parameterisation
#' and its constrained inverse, and replaces only fdHess with a
#' Richardson-extrapolated Hessian (numDeriv). It is independent of the base-R
#' recursion behind the ledger's reference, and carries `output = NA` so the
#' ledger does not count the same values twice. An SD's error maps to the
#' variance's by the delta method, `2 * sd * se(sd)`.
#'
#' @return A comparison `data.frame`.
#' @noRd
.depmixs4_own_errors <- function(label, theirs, reference, obtained) {
  k <- length(reference$estimates$initial)
  pars <- depmixS4::getpars(theirs)
  # The rows of depmixS4's constraint matrix that vcov() keeps: equality
  # constraints that involve at least one parameter.
  equality <- theirs@lin.upper - theirs@lin.lower == 0 &
    apply(theirs@conMat != 0, 1L, any)
  minus_log_likelihood <- function(values) {
    -as.numeric(depmixS4::logLik(depmixS4::setpars(theirs, values)))
  }
  covariance <- depmixS4:::hessian2vcov(numDeriv::hessian(minus_log_likelihood, pars),
                                        theirs@conMat[equality, , drop = FALSE])
  se <- sqrt(diag(covariance))
  width <- (length(pars) - k - k^2) / k
  block <- function(values) matrix(values[k + k^2 + seq_len(k * width)], k, width,
                                   byrow = TRUE)[reference$order_states, , drop = FALSE]
  value_blocks <- block(pars)
  error_blocks <- block(se)
  profiles <- obtained$profiles
  column_of <- match(profiles$indicator, c("y1", "y2"))
  at <- cbind(profiles$profile, 2L * column_of - 1L)
  at_sd <- cbind(profiles$profile, 2L * column_of)
  source_label <- "depmixS4 vcov() with a Richardson Hessian"
  rows <- list(
    compare_block(sprintf("%s: %s, means", label, source_label),
                  error_blocks[at], profiles$mean_standard_error, tolerance = 1e-4),
    compare_block(sprintf("%s: %s, variances (delta method)", label, source_label),
                  2 * value_blocks[at_sd] * error_blocks[at_sd],
                  profiles$variance_standard_error, tolerance = 1e-4))
  responses <- obtained$responses
  if (!is.null(responses) && nrow(responses) > 0L) {
    category <- match(responses$category, colnames(reference$estimates$responses$c1))
    rows <- c(rows, list(compare_block(
      sprintf("%s: %s, category probabilities", label, source_label),
      error_blocks[cbind(responses$profile, 4L + category)],
      responses$probability_standard_error, tolerance = 1e-4)))
  }
  do.call(rbind, rows)
}
