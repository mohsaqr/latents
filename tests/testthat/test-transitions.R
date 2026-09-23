# Latent transition analysis.
#
# The likelihood is checked against explicit enumeration of every state path
# rather than against itself, so a bug in the forward-backward recursion cannot
# be hidden by the recursion being used on both sides of the comparison.

.transition_fixture <- function(n_groups = 9L, occasions = 4L, seed = 99L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(n_groups), each = occasions),
                      t = rep(seq_len(occasions), times = n_groups))
  frame$y1 <- stats::rnorm(nrow(frame))
  frame$y2 <- stats::rnorm(nrow(frame))
  frame$c1 <- sample(c("lo", "hi"), nrow(frame), replace = TRUE)
  frame
}

.transition_recovery_data <- function(seed, n_groups = 120L, occasions = 8L) {
  set.seed(seed)
  transition <- matrix(c(0.85, 0.15, 0.25, 0.75), 2L, 2L, byrow = TRUE)
  frame <- do.call(rbind, lapply(seq_len(n_groups), function(group) {
    path <- unlist(Reduce(function(previous, step) {
      sample.int(2L, 1L, prob = transition[previous, ])
    }, seq_len(occasions - 1L), init = sample.int(2L, 1L, prob = c(0.7, 0.3)),
    accumulate = TRUE))
    data.frame(g = group, t = seq_len(occasions), state = path)
  }))
  frame$y1 <- stats::rnorm(nrow(frame), c(-2, 2)[frame$state], 0.8)
  frame$y2 <- stats::rnorm(nrow(frame), c(-1.5, 1.5)[frame$state], 0.8)
  frame
}

# Sums the complete-data likelihood over every state path of every group and
# every group class. Shares only the model definition with the recursion.
.transition_enumerate <- function(continuous, codes, group_index, time, fit) {
  n_profiles <- fit$n_profiles
  n_rows <- length(group_index)
  density <- vapply(seq_len(n_profiles), function(profile) {
    gaussian <- if (ncol(continuous) == 0L) numeric(n_rows) else
      rowSums(vapply(seq_len(ncol(continuous)), function(indicator) {
        value <- stats::dnorm(continuous[, indicator], fit$means[profile, indicator],
                              sqrt(fit$variances[profile, indicator]), log = TRUE)
        ifelse(is.na(value), 0, value)
      }, numeric(n_rows)))
    categorical <- if (is.null(codes)) numeric(n_rows) else
      rowSums(vapply(seq_along(fit$response_probabilities), function(indicator) {
        value <- log(fit$response_probabilities[[indicator]][profile, codes[, indicator]])
        ifelse(is.na(value), 0, value)
      }, numeric(n_rows)))
    gaussian + categorical
  }, numeric(n_rows))
  sum(vapply(split(seq_along(group_index), group_index), function(rows) {
    rows <- rows[order(time[rows])]
    paths <- as.matrix(expand.grid(rep(list(seq_len(n_profiles)), length(rows))))
    log(sum(vapply(seq_len(fit$n_group_classes), function(type) {
      initial <- fit$initial_probabilities[type, ]
      transition <- matrix(fit$transition_probabilities[, , type],
                           n_profiles, n_profiles)
      fit$group_probabilities[type] * sum(apply(paths, 1L, function(path) {
        moves <- if (length(path) < 2L) 1 else
          prod(transition[cbind(path[-length(path)], path[-1L])])
        initial[path[1L]] * moves * exp(sum(density[cbind(rows, path)]))
      }))
    }, numeric(1))))
  }, numeric(1)))
}

test_that("the recursion reproduces the likelihood enumeration gives", {
  data <- .transition_fixture()
  continuous <- as.matrix(data[, c("y1", "y2")])
  cases <- list(
    list(profiles = 2L, types = 1L), list(profiles = 3L, types = 1L),
    list(profiles = 2L, types = 2L), list(profiles = 3L, types = 2L))
  invisible(lapply(cases, function(case) {
    fit <- quietly(lta(
      data, c("y1", "y2"), "g", n_profiles = case$profiles, time = "t",
      n_group_classes = case$types, n_starts = 2, seed = 4, max_iter = 25))
    expect_equal(fit$log_likelihood,
                 .transition_enumerate(continuous, NULL, fit$group_index,
                                       data$t, fit),
                 tolerance = 1e-10,
                 label = sprintf("%d profiles, %d group classes",
                                 case$profiles, case$types))
  }))
})

test_that("enumeration agrees for categorical, mixed and incomplete indicators", {
  data <- .transition_fixture()
  mixed <- quietly(lta(
    data, c("y1", "y2", "c1"), "g", n_profiles = 2L, time = "t",
    categorical = "c1", n_starts = 2, seed = 8, max_iter = 25))
  expect_equal(mixed$log_likelihood,
               .transition_enumerate(as.matrix(data[, c("y1", "y2")]),
                                     mixed$categorical_data, mixed$group_index,
                                     data$t, mixed), tolerance = 1e-10)

  pure <- quietly(lta(
    data, "c1", "g", n_profiles = 2L, time = "t", categorical = "c1",
    n_starts = 2, seed = 8, max_iter = 25))
  expect_equal(pure$log_likelihood,
               .transition_enumerate(matrix(numeric(0), nrow(data), 0L),
                                     pure$categorical_data, pure$group_index,
                                     data$t, pure), tolerance = 1e-10)
  expect_identical(pure$measurement_model, "categorical")

  incomplete <- data
  incomplete$y1[c(2L, 9L, 20L)] <- NA
  incomplete$y2[c(5L, 9L)] <- NA
  fiml <- quietly(lta(
    incomplete, c("y1", "y2", "c1"), "g", n_profiles = 2L, time = "t",
    categorical = "c1", missing = "fiml", n_starts = 2, seed = 8, max_iter = 25))
  expect_equal(fiml$log_likelihood,
               .transition_enumerate(as.matrix(incomplete[, c("y1", "y2")]),
                                     fiml$categorical_data, fiml$group_index,
                                     incomplete$t, fiml), tolerance = 1e-10)
})

test_that("known transition probabilities are recovered across seeds", {
  recovered <- do.call(rbind, lapply(1:3, function(seed) {
    data <- .transition_recovery_data(seed)
    fit <- quietly(lta(
      data, c("y1", "y2"), "g", n_profiles = 2L, time = "t", n_starts = 4,
      seed = seed))
    # Profile labels are arbitrary, so align them by the first indicator mean
    # before the estimates are compared with the values that generated them.
    aligned <- order(fit$means[, "y1"])
    matrix(fit$transition_probabilities[aligned, aligned, 1L], 2L, 2L)
  }))
  staying <- cbind(recovered[c(TRUE, FALSE), 1L], recovered[c(FALSE, TRUE), 2L])
  expect_equal(mean(staying[, 1L]), 0.85, tolerance = 0.05)
  expect_equal(mean(staying[, 2L]), 0.75, tolerance = 0.05)
})

test_that("the ordering comes from time, so input row order changes nothing", {
  data <- .transition_fixture()
  set.seed(2)
  shuffled <- data[sample.int(nrow(data)), ]
  # Evaluated at one fixed parameter set, so that the comparison isolates the
  # sequence layout from the start values the two orderings would otherwise
  # produce through k-means.
  parameters <- list(
    means = matrix(c(-1, 1, -0.5, 0.5), 2L, 2L), variances = matrix(1, 2L, 2L),
    initial_probabilities = matrix(c(0.6, 0.4), 1L, 2L),
    transition_probabilities = array(c(0.8, 0.3, 0.2, 0.7), c(2L, 2L, 1L)),
    group_probabilities = 1)
  evaluate <- function(frame) {
    groups <- .multilpa_prepare_groups(frame$g)
    layout <- .multilpa_sequence_layout(groups$index, frame$t, groups$n, "observed")
    .multilpa_transition_expectation(as.matrix(frame[, c("y1", "y2")]),
                                     layout, parameters)
  }
  straight <- evaluate(data)
  scrambled <- evaluate(shuffled)
  expect_equal(straight$log_likelihood, scrambled$log_likelihood)
  expect_equal(straight$subject_posteriors,
               scrambled$subject_posteriors[match(rownames(data),
                                                  rownames(shuffled)), ])
})

test_that("enough starts reach the same maximum from either row order", {
  data <- .transition_fixture()
  set.seed(2)
  shuffled <- data[sample.int(nrow(data)), ]
  # The transition likelihood is more multimodal than the cross-sectional one,
  # and k-means starts depend on row order, so agreement is a claim about the
  # restart policy rather than about the algebra. Two starts is not enough.
  maximize <- function(frame) {
    quietly(lta(frame, c("y1", "y2"), "g", n_profiles = 2L,
      time = "t", n_starts = 10, seed = 3, max_iter = 400))$log_likelihood
  }
  expect_equal(maximize(data), maximize(shuffled), tolerance = 1e-6)
})

test_that("a grid of occasions differs from consecutive ones only when a wave is skipped", {
  data <- .transition_fixture(n_groups = 8L, occasions = 5L, seed = 31L)
  fit <- function(frame, occasions) {
    quietly(lta(frame, c("y1", "y2"), "g", n_profiles = 2L,
      time = "t", n_starts = 2, seed = 6, max_iter = 25, occasions = occasions))
  }
  balanced_observed <- fit(data, "observed")
  balanced_grid <- fit(data, "grid")
  expect_equal(balanced_observed$log_likelihood, balanced_grid$log_likelihood)
  expect_true(balanced_observed$balanced)

  gapped <- data[-c(3L, 13L, 14L, 23L), ]
  gapped_observed <- fit(gapped, "observed")
  gapped_grid <- fit(gapped, "grid")
  expect_false(gapped_grid$balanced)
  # Under the grid a skipped position still spans a transition, so a gapped
  # group's span exceeds the rows it contributes; consecutive occasions never do.
  grid_lengths <- get_results(gapped_grid, "sequence_lengths")
  observed_lengths <- get_results(gapped_observed, "sequence_lengths")
  expect_true(any(grid_lengths$occasions > grid_lengths$observations))
  expect_identical(observed_lengths$occasions, observed_lengths$observations)
  expect_identical(sum(grid_lengths$observations), nrow(gapped))
  # A skipped wave costs a second transition under the grid, so the two
  # likelihoods must not coincide; a no-op option would be a silent defect.
  expect_false(isTRUE(all.equal(gapped_observed$log_likelihood,
                                gapped_grid$log_likelihood)))
})

test_that("the fitted quantities satisfy their own definitions", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2"), "g", n_profiles = 3L, time = "t", n_group_classes = 2L,
    n_starts = 2, seed = 5, max_iter = 60))
  expect_equal(rowSums(fit$initial_probabilities), rep(1, 2),
               ignore_attr = TRUE)
  expect_equal(unname(apply(fit$transition_probabilities, c(1L, 3L), sum)),
               matrix(1, 3L, 2L))
  expect_equal(unname(rowSums(fit$subject_posteriors)), rep(1, nrow(data)))
  expect_equal(unname(rowSums(fit$group_posteriors)), rep(1, 9L))
  expect_equal(sum(fit$group_probabilities), 1)
  expect_true(all(fit$transition_probabilities > 0))
  # EM cannot decrease the likelihood it maximizes.
  expect_false(is.unsorted(fit$log_likelihood_history))
  expect_equal(fit$log_likelihood, sum(fit$group_log_likelihood))
  # Every expected transition count is one transition of one group.
  expect_equal(sum(fit$transition_counts), sum(fit$sequence_lengths - 1L))
})

test_that("the free parameter count is the one the model actually has", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2"), "g", n_profiles = 3L, time = "t", n_group_classes = 2L,
    n_starts = 2, seed = 5, max_iter = 20))
  # 3 profiles x 2 indicators means and variances, one group-class split,
  # two initial distributions and two 3x3 transition matrices.
  expect_identical(fit$n_parameters, 12L + 1L + 2L * 2L + 2L * 3L * 2L)
  expect_identical(attr(logLik(fit), "df"), fit$n_parameters)
  expect_identical(nobs(fit), 9L)
  expect_identical(attr(logLik(fit), "nobs"), 9L)
  expect_equal(as.numeric(logLik(fit)), fit$log_likelihood)
  expect_equal(fit$aic, -2 * fit$log_likelihood + 2 * fit$n_parameters)
})

test_that("broken contracts raise their own condition classes", {
  data <- .transition_fixture()
  expect_error(lta(data, c("y1", "y2"), "g", n_profiles = 1L,
                               time = "t", n_starts = 1, seed = 1),
               class = "latents_bad_transition")
  one_occasion <- data[data$t == 1L, ]
  expect_error(lta(one_occasion, c("y1", "y2"), "g", n_profiles = 2L,
                               time = "t", n_starts = 1, seed = 1),
               class = "latents_bad_transition")
  repeated <- data
  repeated$t[2L] <- 1L
  expect_error(lta(repeated, c("y1", "y2"), "g", n_profiles = 2L,
                               time = "t", n_starts = 1, seed = 1),
               class = "latents_bad_time")
  absent <- data
  absent$t[3L] <- NA
  expect_error(lta(absent, c("y1", "y2"), "g", n_profiles = 2L,
                               time = "t", n_starts = 1, seed = 1),
               class = "latents_bad_time")
  expect_error(lta(data, c("y1", "y2"), "g", n_profiles = 2L,
                               time = NULL, n_starts = 1, seed = 1))
  expect_error(get_results(data, "transitions"))
})

test_that("the accessors return the tidy tables they promise", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2", "c1"), "g", n_profiles = 2L, time = "t",
    categorical = "c1", n_group_classes = 2L, n_starts = 2, seed = 5,
    max_iter = 30))

  moves <- get_results(fit, "transitions")
  expect_s3_class(moves, "data.frame")
  expect_named(moves, c("group_class", "from", "to", "probability",
                        "expected_count", "stable", "estimated",
                        "group_class_probability"))
  expect_identical(nrow(moves), 8L)
  expect_identical(moves$stable, moves$from == moves$to)
  expect_true(is.logical(moves$stable) && !anyNA(moves$stable))
  expect_true(is.logical(moves$estimated) && !anyNA(moves$estimated))
  expect_equal(moves$probability,
               as.vector(aperm(fit$transition_probabilities, c(2L, 1L, 3L))))
  # Coercion gives the measurement model for every fitted family now, so the
  # transition matrix is asked for by name.
  expect_identical(get_results(fit, "transitions"), moves)
  expect_identical(as.data.frame(fit), get_results(fit, "profiles"))
  # group_class_probability is a group-class attribute repeated down its rows:
  # it must be constant within a class and sum to one over the classes.
  shares <- split(moves$group_class_probability, moves$group_class)
  expect_true(all(vapply(shares, function(share) length(unique(share)),
                         integer(1)) == 1L))
  expect_equal(sum(vapply(shares, `[[`, numeric(1), 1L)), 1)

  initial <- get_results(fit, "initial")
  expect_named(initial, c("group_class", "profile", "probability", "prevalence",
                          "group_class_probability"))
  expect_equal(as.vector(tapply(initial$probability, initial$group_class, sum)),
               rep(1, 2))

  lengths <- get_results(fit, "sequence_lengths")
  expect_named(lengths, c("group", "group_class", "occasions", "observations",
                          "complete"))
  expect_identical(nrow(lengths), 9L)
  expect_true(all(lengths$complete))
  expect_identical(lengths$occasions, lengths$observations)
  expect_identical(sum(lengths$observations), nrow(data))

  expect_identical(nrow(get_results(fit, "profiles")), 4L)
  expect_identical(nrow(get_results(fit, "responses")), 4L)
  expect_identical(nrow(get_results(fit, "posteriors")),
                   nrow(data) * fit$n_profiles)
  expect_identical(nrow(get_results(fit, "posteriors", format = "wide")),
                   nrow(data))
  expect_identical(nrow(get_results(fit, "group_posteriors")),
                   9L * fit$n_group_classes)
  expect_s3_class(get_results(fit, "information_criteria"), "data.frame")
  expect_s3_class(get_results(fit, "entropy"), "data.frame")
  expect_s3_class(get_results(fit, "classification"), "data.frame")
  expect_output(print(fit), "Latent transition model")
})

test_that("transitions() restricts by argument instead of by bracket", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t",
    n_group_classes = 2L, n_starts = 2, seed = 5, max_iter = 30))

  every <- get_results(fit, "transitions")
  moves <- get_results(fit, "transitions", stable = FALSE)
  stays <- get_results(fit, "transitions", stable = TRUE)

  expect_named(moves, names(every))
  expect_false(any(moves$stable))
  expect_true(all(stays$stable))
  expect_identical(nrow(moves) + nrow(stays), nrow(every))
  # the kept rows keep their order and their values, and are renumbered
  expect_identical(row.names(stays), as.character(seq_len(nrow(stays))))
  expect_equal(stays$probability, every$probability[every$stable])
  # a staying probability per group class and profile
  expect_identical(nrow(stays), fit$n_group_classes * fit$n_profiles)

  # every fitted row is estimated here, and the arguments combine
  expect_true(all(get_results(fit, "transitions", estimated = TRUE)$estimated))
  expect_identical(nrow(get_results(fit, "transitions", estimated = TRUE)), nrow(every))
  expect_identical(nrow(get_results(fit, "transitions", estimated = FALSE)), 0L)
  expect_identical(get_results(fit, "transitions", estimated = TRUE, stable = FALSE), moves)
  # and a restriction belongs to the table that defines it, so asking another
  # table for it is refused by name rather than dropped
  expect_error(get_results(fit, "profiles", stable = FALSE),
               class = "latents_bad_argument")

  expect_error(get_results(fit, "transitions", stable = NA),
               "`stable` must be NULL, TRUE or FALSE")
  expect_error(get_results(fit, "transitions", estimated = c(TRUE, FALSE)),
               "`estimated` must be NULL, TRUE or FALSE")
})

test_that("the summary keeps named fields only, and reports its own tables", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t",
    n_group_classes = 2L, n_starts = 2, seed = 5, max_iter = 30))
  digest <- summary(fit)

  # A Gaussian, diagonal-covariance fit has neither `covariances` nor
  # `response_probabilities`; selecting them by name must not leave NA-named
  # NULL holes in the summary.
  expect_false(anyNA(names(digest)))
  expect_false(any(vapply(digest, is.null, logical(1))))
  expect_false("response_probabilities" %in% names(digest))

  expect_identical(as.data.frame(digest), get_results(fit, "profiles"))
  expect_identical(get_results(digest, "transitions"),
                   get_results(fit, "transitions"))
  # A summary serves its stored tables, so a restriction that would rebuild
  # one is refused rather than silently ignored.
  expect_error(get_results(digest, "transitions", stable = TRUE),
               class = "latents_bad_argument")
  expect_identical(get_results(digest, "initial"),
                   get_results(fit, "initial"))
  expect_identical(get_results(digest, "starts"), fit$starts)
  expect_error(get_results(digest, "nonsense"), class = "latents_bad_argument")

  # a fit missing a mandatory field is refused, not summarized with a hole
  broken <- fit
  broken$aic <- NULL
  expect_error(summary(broken), class = "latents_incomplete_fit")
})

test_that("the shared sequence and diagnostic verbs accept a transition fit", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t", n_group_classes = 2L,
    n_starts = 2, seed = 5, max_iter = 30))
  long <- get_results(fit, "sequences")
  expect_named(long, c("group", "group_class", "time", "profile"))
  expect_identical(nrow(long), nrow(data))
  expect_identical(nrow(get_results(fit, "sequences", format = "wide")), 9L)
  expect_identical(nrow(get_results(fit, "sequence_summary")), 2L)
  expect_s3_class(get_results(fit, "information_criteria"), "data.frame")
  expect_s3_class(get_results(fit, "entropy"), "data.frame")
  expect_s3_class(get_results(fit, "classification", level = "individuals"), "data.frame")
})

test_that("groups of one occasion contribute an initial state and no transition", {
  data <- .transition_fixture()
  short <- rbind(data, data.frame(g = 10L, t = 1L, y1 = 0.4, y2 = -0.2,
                                  c1 = "lo"))
  fit <- quietly(lta(
    short, c("y1", "y2"), "g", n_profiles = 2L, time = "t", n_starts = 2,
    seed = 7, max_iter = 25))
  expect_identical(unname(fit$sequence_lengths[10L]), 1L)
  expect_equal(sum(fit$transition_counts), sum(fit$sequence_lengths - 1L))
  expect_equal(fit$log_likelihood,
               .transition_enumerate(as.matrix(short[, c("y1", "y2")]), NULL,
                                     fit$group_index, short$t, fit),
               tolerance = 1e-10)
})

test_that("the per-occasion offset keeps extreme densities from erasing the priors", {
  set.seed(4)
  data <- data.frame(g = rep(seq_len(6), each = 5L),
                     t = rep(seq_len(5), times = 6L))
  spread <- 1e7
  x <- cbind(stats::rnorm(nrow(data), 0, spread),
             stats::rnorm(nrow(data), 0, spread))
  parameters <- list(
    means = matrix(c(-spread, spread, -spread, spread), 2L, 2L),
    variances = matrix(1e-6, 2L, 2L),
    initial_probabilities = matrix(c(0.9, 0.1), 1L, 2L),
    transition_probabilities = array(c(0.99, 0.01, 0.01, 0.99), c(2L, 2L, 1L)),
    group_probabilities = 1)
  groups <- .multilpa_prepare_groups(data$g)
  layout <- .multilpa_sequence_layout(groups$index, data$t, groups$n, "observed")
  expectation <- .multilpa_transition_expectation(x, layout, parameters)
  # The recursion accumulates log densities across occasions, so its working
  # scale here is about -2e21. Without the offset the log transition
  # probabilities fall below the precision of that sum and the posteriors stop
  # being probabilities; with it they remain normalized.
  expect_lt(min(expectation$group_log_likelihood), -1e20)
  expect_true(is.finite(expectation$log_likelihood))
  expect_equal(rowSums(expectation$subject_posteriors), rep(1, nrow(data)))
  expect_true(all(expectation$subject_posteriors <= 1))
})

test_that("the emission layout removes an offset it can add back exactly", {
  data <- .transition_fixture(n_groups = 5L, occasions = 3L, seed = 12L)
  groups <- .multilpa_prepare_groups(data$g)
  layout <- .multilpa_sequence_layout(groups$index, data$t, groups$n, "observed")
  log_density <- matrix(stats::rnorm(nrow(data) * 2L, -50, 20), nrow(data), 2L)
  emission <- .multilpa_sequence_emission(log_density, layout)
  # Every occasion's largest log density is taken to zero for every group.
  expect_equal(lapply(emission$log_density, function(block) apply(block, 1L, max)),
               rep(list(rep(0, groups$n)), layout$n_occasions))
  # What was taken out is the sum, over occasions, of each group's largest
  # log density at that occasion, computed here without the function under test.
  reference <- Reduce(`+`, lapply(seq_len(layout$n_occasions), function(occasion) {
    rows <- layout$slot[, occasion]
    vapply(seq_along(rows), function(group) {
      if (is.na(rows[group])) 0 else max(log_density[rows[group], ])
    }, numeric(1))
  }))
  expect_equal(emission$offset, reference)
  # Adding each occasion's offset back recovers the densities it was given.
  round_trip <- do.call(rbind, lapply(seq_len(layout$n_occasions), function(occasion) {
    rows <- layout$slot[, occasion]
    observed <- !is.na(rows)
    original <- log_density[rows[observed], , drop = FALSE]
    emission$log_density[[occasion]][observed, , drop = FALSE] +
      apply(original, 1L, max) - original
  }))
  expect_equal(round_trip, matrix(0, nrow(data), 2L))
})

test_that("the generics either answer or refuse, and never answer emptily", {
  data <- .transition_fixture()
  fit <- quietly(lta(
    data, c("y1", "y2", "c1"), "g", n_profiles = 2L, time = "t",
    categorical = "c1", n_starts = 2, seed = 5, max_iter = 30))

  estimates <- coef(fit)
  expect_true(is.numeric(estimates) && !is.null(names(estimates)))
  expect_false(anyNA(estimates))
  # Every free quantity is present: means, variances, responses, initial
  # probabilities, transitions and the group-class share.
  # Names follow the package-wide `level.parameter.outcome.term` grammar, shared
  # with coef.multilpa() and coef.multilpa_covariates(), so one regex works on
  # every fit class.
  expect_identical(sum(grepl("^measurement[.]mean[.]", names(estimates))), 4L)
  expect_identical(sum(grepl("^measurement[.]variance[.]", names(estimates))), 4L)
  expect_identical(sum(grepl("^measurement[.]response[.]", names(estimates))), 4L)
  expect_identical(
    sum(grepl("^profile[.]initial_probability[.]", names(estimates))), 2L)
  expect_identical(
    sum(grepl("^profile[.]transition_probability[.]", names(estimates))), 4L)
  expect_equal(
    unname(estimates[grepl("^profile[.]transition_probability[.]", names(estimates))]),
    get_results(fit, "transitions")$probability)
  expect_identical(sum(grepl("^group[.]probability[.]", names(estimates))), 1L)
  # Every name parses back into the four tidy fields: level, parameter, outcome
  # and a term that may itself contain dots or colons.
  expect_true(all(lengths(regmatches(names(estimates),
                                     gregexpr("[.]", names(estimates)))) >= 2L))

  digest <- summary(fit)
  expect_s3_class(digest, "summary_multilpa_transitions")
  expect_identical(digest$n_profiles, 2L)
  expect_output(print(digest), "-- transitions", fixed = TRUE)

  # Refusals are classed, so that a caller can catch them, and are raised
  # instead of a default method returning nothing useful.
  expect_error(vcov(fit), class = "latents_no_inference")
  expect_error(confint(fit), class = "latents_no_inference")
  expect_error(parameter_inference(fit, data), class = "latents_no_inference")
  expect_s3_class(draw(plot(fit, what = "transitions")), "multilpa_transitions")
  expect_s3_class(draw(plot(fit, what = "profiles")), "multilpa_transitions")
})
