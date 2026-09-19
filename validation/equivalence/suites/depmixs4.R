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

suite_depmixs4 <- function() {
  if (!requireNamespace("depmixS4", quietly = TRUE)) {
    stop("suite 'depmixs4' needs the depmixS4 package")
  }
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
  do.call(rbind, lapply(cases, function(case) do.call(.depmixs4_case, case)))
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
#' @return A `data.frame` of compared quantities.
#' @noRd
.depmixs4_case <- function(dataset, profiles, categorical, groups, occasions,
                           seed, drop = 0) {
  frame <- .depmixs4_data(profiles, groups, occasions, seed, categorical, drop)
  indicators <- c("y1", "y2", if (categorical) "c1")
  fit <- suppressWarnings(fit_transitions(
    frame, indicators, "g", n_profiles = profiles, time = "t",
    categorical = if (categorical) "c1" else character(),
    n_starts = 8, seed = seed, max_iter = 2000, tol = 1e-12))
  formulas <- lapply(indicators, function(name) stats::as.formula(paste(name, "~ 1")))
  families <- c(list(stats::gaussian(), stats::gaussian()),
                if (categorical) list(depmixS4::multinomial("identity")))
  model <- depmixS4::depmix(formulas, data = frame, nstates = profiles,
                            family = families,
                            ntimes = as.integer(table(frame$g)[unique(as.character(frame$g))]))

  # Claim one: the same likelihood function. depmixS4 evaluated exactly at
  # multilpa's estimates must return multilpa's log likelihood.
  at_ours <- as.numeric(depmixS4::logLik(
    depmixS4::setpars(model, .depmixs4_pars(fit, categorical))))
  same_function <- compare_values(
    quantity = "log likelihood at multilpa's estimates",
    reference = at_ours, obtained = fit$log_likelihood,
    tolerance = 1e-8, scale = "relative")

  # Claim two: the same maximum, found independently. depmixS4 searches from
  # its own start; agreement then says both optimizers reached the same mode.
  set.seed(seed)
  theirs <- suppressWarnings(depmixS4::fit(model, verbose = FALSE, emcontrol =
    depmixS4::em.control(maxit = 2000, tol = 1e-12)))
  their_pars <- depmixS4::getpars(theirs)
  same_maximum <- compare_values(
    quantity = "maximized log likelihood",
    reference = as.numeric(depmixS4::logLik(theirs)),
    obtained = fit$log_likelihood, tolerance = 1e-6, scale = "relative")

  # State labels are arbitrary in both packages, so the estimates are compared
  # only after each is put in its own order of the first indicator's mean.
  # depmixS4 lays its parameters out as the prior, then the transition matrix
  # by row, then one equally sized block per state. The block width is derived
  # rather than assumed, because it depends on how many response families the
  # model carries; y1's intercept is the first entry of each block.
  block <- (length(their_pars) - profiles - profiles^2) / profiles
  stopifnot("depmixS4 parameter blocks must be equally sized" =
              block == floor(block))
  their_means <- their_pars[profiles + profiles^2 +
                              (seq_len(profiles) - 1L) * block + 1L]
  their_order <- order(their_means)
  our_order <- order(fit$means[, "y1"])
  their_transition <- matrix(their_pars[profiles + seq_len(profiles^2)],
                             profiles, profiles, byrow = TRUE)[their_order, their_order]
  our_transition <- matrix(fit$transition_probabilities[, , 1L], profiles,
                           profiles)[our_order, our_order]
  same_transitions <- compare_values(
    quantity = sprintf("transition %d to %d",
                       rep(seq_len(profiles), each = profiles),
                       rep(seq_len(profiles), times = profiles)),
    reference = as.vector(t(their_transition)),
    obtained = as.vector(t(our_transition)),
    tolerance = 1e-5)
  same_means <- compare_values(
    quantity = sprintf("state %d mean of y1", seq_len(profiles)),
    reference = sort(their_means), obtained = fit$means[our_order, "y1"],
    tolerance = 1e-5)

  rows <- rbind(same_function, same_maximum, same_transitions, same_means)
  data.frame(source = "depmixS4", dataset = dataset, precision = "machine",
             rows, stringsAsFactors = FALSE)
}
