# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# lta() fits, their accessors and transition grids against depmixS4 (Visser and
# Speekenbrink, 2010), a hidden Markov program written separately: a latent
# transition model with Gaussian indicators is a hidden Markov model with one
# Gaussian response per indicator. Two ways:
# - depmixS4's forward-backward pass at latents' own estimates (setpars(), no
#   fitting): the log likelihood, the state posteriors and the expected
#   transition counts (xi) of one E-step, which agree to rounding;
# - depmixS4's EM continued from latents' estimates to a relative change of
#   1e-12 (and, for the grids, from random starts), for the maximum.
# latents stops at a relative change of 1e-12 in the likelihood; its
# estimates then sit within about 1e-6 of the maximum (the square root),
# hence 1e-5 on parameters that depmixS4 reaches by a second EM.
# depmixS4 lets a covariate at occasion t govern the move from t to t + 1;
# latents uses the covariate of the destination occasion, so depmixS4 is
# given the covariate one occasion ahead (`z_ahead`).

transition_panel <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  n <- 150L
  occasions <- 4L
  kind <- 1L + (stats::runif(n) < 0.5)
  z <- matrix(stats::rnorm(n * occasions), n)
  states <- do.call(cbind, Reduce(function(previous, t) {
    move <- stats::runif(n) < stats::plogis(c(-1.8, -0.2)[kind] + 0.8 * z[, t])
    ifelse(move, 3L - previous, previous)
  }, seq_len(occasions)[-1L], init = 1L + (stats::runif(n) < 0.4), accumulate = TRUE))
  s <- as.vector(t(states))
  data.frame(id = rep(seq_len(n), each = occasions), time = rep(seq_len(occasions), n),
             z = as.vector(t(z)), w = rep(stats::rnorm(n), each = occasions),
             y1 = stats::rnorm(length(s), c(-1, 1.2)[s], c(0.8, 1)[s]),
             y2 = stats::rnorm(length(s), c(0, 1.5)[s], 0.9))
})
# The panel is sorted by id and time: the covariate of the next occasion of
# the same sequence (0 on the last occasion, where no move follows).
transition_panel$z_ahead <- stats::ave(transition_panel$z, transition_panel$id, FUN = \(v) c(v[-1L], 0))
panel_n <- 150L
panel_states <- c("profile_1", "profile_2")

homogeneous_fit <- lta(transition_panel, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                       n_starts = 2, seed = 1, tol = 1e-12)
covariate_fit <- lta(transition_panel, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                     transition_covariates = "z", n_starts = 2, seed = 1, tol = 1e-12)
grouped_fit <- quietly(lta(transition_panel, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                           n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-12))

# depmixS4 prints a line when its EM stops; the line is not needed.
depmix_fit <- function(model, random_start, seed = 1L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  printed <- utils::capture.output(fitted <- depmixS4::fit(
    model, emcontrol = depmixS4::em.control(maxit = 100000L, tol = 1e-12, crit = "relative",
                                             random.start = random_start),
    verbose = FALSE))
  fitted
}

# The Gaussian responses of depmixS4 per state: y1 intercept and sd, y2 intercept and sd.
response_parameters <- function(estimates, k, profile_of = seq_len(k)) {
  unlist(lapply(profile_of, \(p) {
    c(estimates[[sprintf("measurement.mean.profile_%d.y1", p)]],
      sqrt(estimates[[sprintf("measurement.variance.profile_%d.y1", p)]]),
      estimates[[sprintf("measurement.mean.profile_%d.y2", p)]],
      sqrt(estimates[[sprintf("measurement.variance.profile_%d.y2", p)]]))
  }))
}

# A two-state depmixS4 model of the panel; `transition` a formula for the moves.
panel_model <- function(transition = ~ 1, nstates = 2L) {
  depmixS4::depmix(list(y1 ~ 1, y2 ~ 1), data = transition_panel, nstates = nstates,
                   family = list(stats::gaussian(), stats::gaussian()), ntimes = rep(4L, panel_n),
                   transition = transition)
}

# The parameters of the homogeneous fit in depmixS4's order: initial
# probabilities, transition rows (from profile 1, from profile 2), responses.
homogeneous_parameters <- function(estimates) {
  at <- \(name, ...) estimates[[sprintf(name, ...)]]
  initial <- vapply(1:2, \(k) at("profile.initial_probability.profile_%d.group_class_1", k), numeric(1))
  moves <- unlist(lapply(1:2, \(from) vapply(1:2, \(to) {
    at("profile.transition_probability.profile_%d.group_class_1:profile_%d", to, from)
  }, numeric(1))))
  c(initial, moves, response_parameters(estimates, 2L))
}

# Expected transition counts from depmixS4's xi, normalised by row. xi[t, j, i]
# is the posterior probability of profile i at row t and profile j at row
# t + 1 (destination second, origin third: its sums over the second index are
# the posteriors of row t); the last occasion of a sequence carries no move.
xi_network <- function(xi, states = 1:2) {
  counts <- outer(states, states, Vectorize(\(i, j) sum(xi[, j, i])))
  matrix(counts / rowSums(counts), length(states), dimnames = list(panel_states, panel_states))
}
moves_only <- function(xi) xi[rep(c(TRUE, TRUE, TRUE, FALSE), panel_n), , , drop = FALSE]

posterior_table <- function(gamma) {
  data.frame(row = rep(seq_len(nrow(gamma)), ncol(gamma)), profile = rep(seq_len(ncol(gamma)), each = nrow(gamma)),
             posterior = as.vector(gamma))
}

test_that("homogeneous lta(): depmixS4's forward-backward pass at latents' estimates", {
  skip_without_reference("depmixS4")
  at_latents <- depmixS4::setpars(panel_model(), homogeneous_parameters(coef(homogeneous_fit)))
  passed <- depmixS4::forwardbackward(at_latents)
  # One E-step evaluated by two programs on the same parameters: rounding only.
  expect_reference(as.numeric(logLik(homogeneous_fit)), passed$logLike, reference = "depmixS4",
                   covers = c("lta()", "logLik.multilpa_transitions()"),
                   quantity = "log likelihood at latents' estimates, homogeneous", tolerance = 1e-10,
                   measure = "relative")
  # Probabilities, many near 0: absolute.
  expect_reference(get_results(homogeneous_fit, "posteriors"), posterior_table(passed$gamma), reference = "depmixS4",
                   covers = c("lta()", "get_results.multilpa_transitions()"),
                   quantity = "row posteriors, homogeneous", tolerance = 1e-10, measure = "absolute",
                   key = c("row", "profile"))
  # get_tna() normalises the expected transition counts of the fit's E-step.
  network <- get_tna(homogeneous_fit)
  expect_reference(network$weights, xi_network(moves_only(passed$xi)), reference = "depmixS4",
                   covers = c("get_tna()", "get_tna.multilpa_transitions()"),
                   quantity = "transition network weights (expected counts from xi), homogeneous",
                   tolerance = 1e-10, measure = "relative")
  # Its initial distribution is the fitted one; at the maximum it equals the
  # mean posterior of the first occasion (the EM fixed point), to the fit's
  # convergence: 1e-6.
  first <- passed$gamma[rep(c(TRUE, FALSE, FALSE, FALSE), panel_n), , drop = FALSE]
  expect_reference(network$inits, stats::setNames(colMeans(first), panel_states), reference = "depmixS4",
                   covers = c("get_tna()", "get_tna.multilpa_transitions()"),
                   quantity = "initial distribution of the network, homogeneous", tolerance = 1e-6,
                   measure = "relative")
})

test_that("homogeneous lta(): the estimates are the maximum depmixS4's EM reaches", {
  skip_without_reference("depmixS4")
  continued <- depmix_fit(depmixS4::setpars(panel_model(), homogeneous_parameters(coef(homogeneous_fit))),
                          random_start = FALSE)
  # From random starts too: the maximum is the one depmixS4 finds by itself (best of three).
  random <- lapply(1:3, \(seed) depmix_fit(panel_model(), random_start = TRUE, seed = seed))
  best <- max(vapply(random, \(f) as.numeric(stats4::logLik(f)), numeric(1)))
  expect_reference(rep(as.numeric(logLik(homogeneous_fit)), 2L), c(continued = as.numeric(stats4::logLik(continued)), random = best),
                   reference = "depmixS4", covers = c("lta()", "logLik.multilpa_transitions()"),
                   quantity = "maximised log likelihood, homogeneous (continued; best of three random starts)",
                   tolerance = 1e-8, measure = "relative")
  estimates <- depmixS4::getpars(continued)
  names(estimates) <- NULL
  # depmixS4's order: initial (2), transition rows (4), responses (8).
  natural <- c(sprintf("profile.initial_probability.profile_%d.group_class_1", 1:2),
               sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d",
                       c(1L, 2L, 1L, 2L), c(1L, 1L, 2L, 2L)))
  # 1e-5: both EMs stop at a relative change of 1e-12 (see the head of the file).
  expect_reference(coef(homogeneous_fit)[natural], stats::setNames(estimates[1:6], natural),
                   reference = "depmixS4", covers = c("lta()", "coef.multilpa_transitions()"),
                   quantity = "initial and transition probabilities, homogeneous", tolerance = 1e-5,
                   measure = "relative")
  sds <- estimates[c(8L, 10L, 12L, 14L)]
  profiles <- data.frame(profile = rep(1:2, each = 2), indicator = rep(c("y1", "y2"), 2),
                         mean = estimates[c(7L, 9L, 11L, 13L)], variance = sds^2, standard_deviation = sds)
  expect_reference(as.data.frame(homogeneous_fit), profiles, reference = "depmixS4",
                   covers = c("lta()", "as.data.frame.multilpa_transitions()"),
                   quantity = "profile means, variances and sds, homogeneous", tolerance = 1e-5,
                   measure = "relative", key = c("profile", "indicator"))
  expect_reference(get_results(summary(homogeneous_fit)), profiles, reference = "depmixS4",
                   covers = c("summary.multilpa_transitions()", "get_results.summary_multilpa_transitions()"),
                   quantity = "profile means, variances and sds, get_results(summary()), homogeneous",
                   tolerance = 1e-5, measure = "relative", key = c("profile", "indicator"))
  expect_reference(as.data.frame(summary(homogeneous_fit)), profiles, reference = "depmixS4",
                   covers = c("summary.multilpa_transitions()", "as.data.frame.summary_multilpa_transitions()"),
                   quantity = "profile means, variances and sds, summary(), homogeneous",
                   tolerance = 1e-5, measure = "relative", key = c("profile", "indicator"))
  expect_reference(get_results(homogeneous_fit, "profiles"), profiles[c("profile", "indicator", "mean", "variance")],
                   reference = "depmixS4", covers = c("lta()", "get_results.multilpa_transitions()"),
                   quantity = "profile table, get_results(), homogeneous", tolerance = 1e-5,
                   measure = "relative", key = c("profile", "indicator"))
})

# The covariate fit's parameters in depmixS4's order: initial probabilities;
# for each origin the multinomial logits against profile 1 (intercepts of
# profiles 1 and 2, then slopes on z_ahead); responses. From profile 2 the
# logit of profile 2 against profile 1 is minus latents' logit of leaving.
covariate_parameters <- function(estimates) {
  at <- \(name) estimates[[name]]
  initial <- stats::plogis(at("initial.group_class_1.profile_1.(Intercept)"))
  leave_1 <- c(at("transition.group_class_1.profile_1->profile_2.(Intercept)"),
               at("transition.group_class_1.profile_1->profile_2.z"))
  leave_2 <- c(at("transition.group_class_1.profile_2->profile_1.(Intercept)"),
               at("transition.group_class_1.profile_2->profile_1.z"))
  responses <- unlist(lapply(1:2, \(p) {
    c(at(sprintf("mean.all.profile_%d.y1", p)), exp(at(sprintf("log_variance.all.profile_%d.y1", p)) / 2),
      at(sprintf("mean.all.profile_%d.y2", p)), exp(at(sprintf("log_variance.all.profile_%d.y2", p)) / 2))
  }))
  c(initial, 1 - initial, 0, leave_1[[1L]], 0, leave_1[[2L]], 0, -leave_2[[1L]], 0, -leave_2[[2L]], responses)
}

test_that("covariate lta(): depmixS4's forward-backward pass at latents' estimates", {
  skip_without_reference("depmixS4")
  at_latents <- depmixS4::setpars(panel_model(~ z_ahead), covariate_parameters(coef(covariate_fit)))
  passed <- depmixS4::forwardbackward(at_latents)
  expect_reference(as.numeric(logLik(covariate_fit)), passed$logLike, reference = "depmixS4",
                   covers = c("lta()", "logLik.multilpa_transitions()"),
                   quantity = "log likelihood at latents' estimates, covariate transitions", tolerance = 1e-10,
                   measure = "relative")
  expect_reference(get_results(covariate_fit, "posteriors"), posterior_table(passed$gamma), reference = "depmixS4",
                   covers = c("lta()", "get_results.multilpa_lta()"),
                   quantity = "row posteriors, covariate transitions", tolerance = 1e-10, measure = "absolute",
                   key = c("row", "profile"))
})

test_that("covariate lta(): the estimates and transition table are the maximum depmixS4's EM reaches", {
  skip_without_reference("depmixS4")
  continued <- depmix_fit(depmixS4::setpars(panel_model(~ z_ahead), covariate_parameters(coef(covariate_fit))),
                          random_start = FALSE)
  expect_reference(as.numeric(logLik(covariate_fit)), as.numeric(stats4::logLik(continued)), reference = "depmixS4",
                   covers = c("lta()", "logLik.multilpa_transitions()"),
                   quantity = "maximised log likelihood, covariate transitions", tolerance = 1e-8,
                   measure = "relative")
  p <- depmixS4::getpars(continued)
  names(p) <- NULL
  expected <- c("initial.group_class_1.profile_1.(Intercept)" = stats::qlogis(p[[1L]]),
                "transition.group_class_1.profile_1->profile_2.(Intercept)" = p[[4L]],
                "transition.group_class_1.profile_1->profile_2.z" = p[[6L]],
                "transition.group_class_1.profile_2->profile_1.(Intercept)" = -p[[8L]],
                "transition.group_class_1.profile_2->profile_1.z" = -p[[10L]],
                "mean.all.profile_1.y1" = p[[11L]], "log_variance.all.profile_1.y1" = 2 * log(p[[12L]]),
                "mean.all.profile_1.y2" = p[[13L]], "log_variance.all.profile_1.y2" = 2 * log(p[[14L]]),
                "mean.all.profile_2.y1" = p[[15L]], "log_variance.all.profile_2.y1" = 2 * log(p[[16L]]),
                "mean.all.profile_2.y2" = p[[17L]], "log_variance.all.profile_2.y2" = 2 * log(p[[18L]]))
  # Unconstrained coordinates, some near zero (log variances, a mean of 0.13): scaled, 1e-5.
  expect_reference(coef(covariate_fit)[names(expected)], expected, reference = "depmixS4",
                   covers = c("lta()", "coef.multilpa_lta()"),
                   quantity = "unconstrained estimates, covariate transitions", tolerance = 1e-5,
                   measure = "scaled")
  # depmixS4's trDens holds each row's move probabilities, trDens[t, to, from]
  # for the move from row t to row t + 1 (its forward pass multiplies it into
  # the column of the origin); averaged over the sequences at each destination
  # occasion.
  moves <- expand.grid(to = 1:2, from = 1:2, occasion = 2:4)
  origin <- transition_panel$time
  moves$probability <- mapply(\(to, from, occasion) mean(continued@trDens[origin == occasion - 1L, to, from]),
                              moves$to, moves$from, moves$occasion)
  expected_moves <- data.frame(group_class = "group_class_1", occasion = moves$occasion,
                               from = panel_states[moves$from], to = panel_states[moves$to],
                               probability = moves$probability)
  expect_reference(as.data.frame(covariate_fit), expected_moves, reference = "depmixS4",
                   covers = c("lta()", "as.data.frame.multilpa_lta()"),
                   quantity = "transition probabilities by occasion, covariate transitions", tolerance = 1e-5,
                   measure = "relative", key = c("group_class", "occasion", "from", "to"))
  expect_reference(as.data.frame(summary(covariate_fit)), expected_moves, reference = "depmixS4",
                   covers = c("summary.multilpa_lta()", "as.data.frame.summary_multilpa_lta()"),
                   quantity = "transition probabilities by occasion, summary(), covariate transitions",
                   tolerance = 1e-5, measure = "relative", key = c("group_class", "occasion", "from", "to"))
})

# Two group classes of a homogeneous transition model are one hidden Markov
# model of four states (class h, profile k) = 2 (h - 1) + k whose transitions
# never cross classes and whose responses depend on the profile alone.
grouped_parameters <- function(estimates) {
  at <- \(name, ...) estimates[[sprintf(name, ...)]]
  cells <- expand.grid(k = 1:2, h = 1:2)
  initial <- mapply(\(k, h) at("group.probability.group_class_%d", h) *
                      at("profile.initial_probability.profile_%d.group_class_%d", k, h), cells$k, cells$h)
  moves <- unlist(lapply(seq_len(nrow(cells)), \(row) {
    from <- cells$k[row]
    h <- cells$h[row]
    vapply(seq_len(nrow(cells)), \(column) {
      if (cells$h[column] != h) 0
      else at("profile.transition_probability.profile_%d.group_class_%d:profile_%d", cells$k[column], h, from)
    }, numeric(1))
  }))
  c(initial, moves, response_parameters(estimates, 2L, profile_of = cells$k))
}

test_that("lta() with two group classes: depmixS4's forward-backward pass of the four-state model", {
  skip_without_reference("depmixS4")
  at_latents <- depmixS4::setpars(panel_model(nstates = 4L), grouped_parameters(coef(grouped_fit)))
  passed <- depmixS4::forwardbackward(at_latents)
  expect_reference(as.numeric(logLik(grouped_fit)), passed$logLike, reference = "depmixS4",
                   covers = c("lta()", "logLik.multilpa_transitions()"),
                   quantity = "log likelihood at latents' estimates, two group classes", tolerance = 1e-10,
                   measure = "relative")
  xi <- moves_only(passed$xi)
  # Pooled: the expected moves of profile k to l summed over both classes.
  pooled <- xi[, 1:2, 1:2, drop = FALSE] + xi[, 3:4, 3:4, drop = FALSE]
  expect_reference(get_tna(grouped_fit)$weights, xi_network(pooled), reference = "depmixS4",
                   covers = c("get_tna()", "get_tna.multilpa_transitions()"),
                   quantity = "pooled transition network weights, two group classes", tolerance = 1e-10,
                   measure = "relative")
  # Per class: that class's expected moves; the network is the class's
  # estimated matrix, equal to them at the maximum to the fit's convergence: 1e-6.
  groups <- get_group_tna(grouped_fit)
  table_of <- \(h, weights) data.frame(class = h, from = rep(panel_states, 2), to = rep(panel_states, each = 2),
                                       weight = as.vector(weights))
  expect_reference(do.call(rbind, lapply(1:2, \(h) table_of(h, groups[[h]]$weights))),
                   do.call(rbind, lapply(1:2, \(h) {
                     block <- c(2L * h - 1L, 2L * h)
                     table_of(h, xi_network(xi[, block, block, drop = FALSE]))
                   })),
                   reference = "depmixS4", covers = c("get_group_tna()", "get_group_tna.multilpa_transitions()"),
                   quantity = "transition network weights per group class, two group classes", tolerance = 1e-6,
                   measure = "relative", key = c("class", "from", "to"))
})

# Sixty sequences of three occasions, two states that persist with probability 0.8
# (the transition grid of test-accessors-selection.R). Three profiles overfit
# these data: their highest maximum (-591.04, found by depmixS4 in 1 of 20
# random starts and by lta() with 300 starts) puts 1.9 units in a profile with
# a variance of 0.008, and lower maxima at -592.58 and -593.89 are common. So
# each candidate is compared with depmixS4's EM continued from its estimates
# (the maximum it reports), and the two-profile candidate also with depmixS4's
# best of ten random starts.
transition_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  n <- 60L
  first <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  states <- Reduce(\(previous, occasion) ifelse(stats::runif(n) < 0.8, previous, 3L - previous),
                   2:3, accumulate = TRUE, init = first)
  state <- as.vector(t(do.call(cbind, states)))
  data.frame(id = rep(seq_len(n), each = 3L), time = rep(1:3, n),
             y1 = stats::rnorm(3L * n, c(0, 2.5)[state]), y2 = stats::rnorm(3L * n, c(0, 2)[state]))
})
transition_grid <- quietly(enumerate_classes(transition_data, c("y1", "y2"), "id", time = "time",
                                             n_profiles = 2:3, n_group_classes = 1,
                                             n_starts = 3, seed = 1))

# A homogeneous candidate's parameters in depmixS4's order (k profiles).
candidate_parameters <- function(estimates, k) {
  at <- \(name, ...) estimates[[sprintf(name, ...)]]
  initial <- vapply(seq_len(k), \(j) at("profile.initial_probability.profile_%d.group_class_1", j), numeric(1))
  moves <- unlist(lapply(seq_len(k), \(from) vapply(seq_len(k), \(to) {
    at("profile.transition_probability.profile_%d.group_class_1:profile_%d", to, from)
  }, numeric(1))))
  c(initial, moves, response_parameters(estimates, k))
}

grid_row <- function(fitted, k) {
  log_likelihood <- stats4::logLik(fitted)
  # depmixS4's BIC takes the number of rows (occasions), latents' bic_individual.
  data.frame(candidate = sprintf("profiles%d_groups1", k), log_likelihood = as.numeric(log_likelihood),
             n_parameters = attr(log_likelihood, "df"), aic = stats4::AIC(fitted),
             bic_individual = stats4::BIC(fitted))
}

test_that("the transition grid's maxima, parameter counts and criteria are depmixS4's", {
  skip_without_reference("depmixS4")
  model_of <- \(k) depmixS4::depmix(list(y1 ~ 1, y2 ~ 1), data = transition_data, nstates = k,
                                    family = list(stats::gaussian(), stats::gaussian()), ntimes = rep(3L, 60L))
  expected <- do.call(rbind, lapply(2:3, \(k) {
    estimates <- coef(get_candidate(transition_grid, n_profiles = k))
    grid_row(depmix_fit(depmixS4::setpars(model_of(k), candidate_parameters(estimates, k)), random_start = FALSE), k)
  }))
  # Both EMs stop at a relative change of 1e-12 or below: 1e-8 on the likelihood and criteria.
  expect_reference(as.data.frame(transition_grid), expected, reference = "depmixS4",
                   covers = c("enumerate_classes()", "as.data.frame.latents_transition_enumeration()"),
                   quantity = "maximum, parameters and criteria per candidate (depmixS4 continued)", tolerance = 1e-8,
                   measure = "relative", key = "candidate")
  expect_reference(summary(transition_grid), expected, reference = "depmixS4",
                   covers = c("enumerate_classes()", "summary.latents_transition_enumeration()"),
                   quantity = "maximum, parameters and criteria per candidate, summary() (depmixS4 continued)",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  expect_reference(get_results(transition_grid), expected, reference = "depmixS4",
                   covers = c("enumerate_classes()", "get_results.latents_transition_enumeration()"),
                   quantity = "maximum, parameters and criteria per candidate, get_results() (depmixS4 continued)",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  # grid$<candidate> is the candidate's fit.
  expect_reference(as.numeric(logLik(transition_grid$profiles2_groups1)), expected$log_likelihood[[1L]],
                   reference = "depmixS4", covers = c("enumerate_classes()", "$.latents_transition_enumeration()"),
                   quantity = "log likelihood of grid$profiles2_groups1 (depmixS4 continued)", tolerance = 1e-8,
                   measure = "relative")
  fits <- lapply(1:10, \(seed) depmix_fit(model_of(2L), random_start = TRUE, seed = seed))
  best <- fits[[which.max(vapply(fits, \(f) as.numeric(stats4::logLik(f)), numeric(1)))]]
  expect_reference(subset(as.data.frame(transition_grid), candidate == "profiles2_groups1"), grid_row(best, 2L),
                   reference = "depmixS4",
                   covers = c("enumerate_classes()", "as.data.frame.latents_transition_enumeration()"),
                   quantity = "two-profile maximum and criteria (depmixS4, best of ten random starts)",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
})

test_that("predict() of lta() fits is not implemented (a declared API gap)", {
  # The reference a result would be compared with: depmixS4's forward-backward
  # posteriors of the new sequences, as above for the fitted ones.
  expect_defect_error(predict(homogeneous_fit, newdata = transition_panel),
                      defect = "predict() is not implemented for latent transition models: it always raises latents_unsupported_prediction; FINDINGS.md (latents, API gaps)",
                      covers = "predict.multilpa_transitions()",
                      quantity = "profile posteriors of new sequences, homogeneous")
  expect_defect_error(predict(covariate_fit, newdata = transition_panel),
                      defect = "predict() is not implemented for latent transition models with covariates: it always raises latents_unsupported_prediction; FINDINGS.md (latents, API gaps)",
                      covers = "predict.multilpa_lta()",
                      quantity = "profile posteriors of new sequences, covariate transitions")
})
