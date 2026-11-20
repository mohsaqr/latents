# The accessors of lta() fits (classes multilpa_transitions and multilpa_lta)
# and the transition networks get_tna() / get_group_tna() hand to tna, on a
# simulated panel: 150 sequences of four occasions, two Gaussian indicators,
# two profiles. References, all written here in base R:
# - the likelihood summed over every profile path (16 per sequence), as a
#   function of unconstrained coordinates (means, log variances, initial and
#   transition logits against the reference profile / against staying);
# - its maximum, reached by Newton steps from latents' estimates with
#   numDeriv's gradient and Hessian, for the estimates;
# - numDeriv's Hessian at that maximum, inverted, for the covariance; natural-
#   scale standard errors by the delta method (numDeriv's Jacobian); Wald
#   intervals and z tests from those standard errors;
# - posteriors and expected transition counts from the same path enumeration
#   (Bayes' rule over the paths);
# - for multilpa() fits made with `time`, row posteriors from Bayes' rule over
#   the group classes, their modal sequences, and the consecutive-pair counts
#   tna's relative model is defined by.

transition_panel <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  n <- 150L
  occasions <- 4L
  # Two kinds of sequences, one stickier than the other, so that a model with
  # two group classes has something to find.
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

# The panel is balanced and sorted by id and time: one row per sequence, one
# column per occasion.
panel_n <- 150L
panel_y1 <- matrix(transition_panel$y1, panel_n, 4L, byrow = TRUE)
panel_y2 <- matrix(transition_panel$y2, panel_n, 4L, byrow = TRUE)
panel_z <- matrix(transition_panel$z, panel_n, 4L, byrow = TRUE)
panel_w <- matrix(transition_panel$w, panel_n, 4L, byrow = TRUE)[, 1L]

homogeneous_fit <- lta(transition_panel, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                       n_starts = 2, seed = 1, tol = 1e-12)
covariate_fit <- lta(transition_panel, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                     transition_covariates = "z", n_starts = 2, seed = 1, tol = 1e-12)

# ---- Reference code -------------------------------------------------------

# Every profile path of four occasions, one row per path.
panel_paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
# The moves each path makes: a 2 x 2 count (from, to) per path.
panel_path_moves <- lapply(seq_len(nrow(panel_paths)), \(p) {
  path <- panel_paths[p, ]
  outer(1:2, 1:2, Vectorize(\(k, l) sum(path[-4L] == k & path[-1L] == l)))
})

rows_log_sum_exp <- function(values) {
  top <- apply(values, 1L, max)
  top + log(rowSums(exp(values - top)))
}

# Log density of each sequence's rows under each profile (sequences x
# occasions), diagonal normal: `means` and `variances` are profile x indicator.
panel_log_density <- function(means, variances) {
  lapply(1:2, \(k) {
    stats::dnorm(panel_y1, means[k, 1L], sqrt(variances[k, 1L]), log = TRUE) +
      stats::dnorm(panel_y2, means[k, 2L], sqrt(variances[k, 2L]), log = TRUE)
  })
}

# Log probability of every path of every sequence (sequences x paths).
# `log_initial`: sequences x profiles; `log_move(from, t)`: sequences x
# destination profiles, the log probabilities of the move into occasion t.
panel_path_log <- function(log_density, log_initial, log_move) {
  vapply(seq_len(nrow(panel_paths)), \(p) {
    s <- panel_paths[p, ]
    log_initial[, s[1L]] + log_density[[s[1L]]][, 1L] +
      Reduce(`+`, lapply(2:4, \(t) log_move(s[t - 1L], t)[, s[t]] + log_density[[s[t]]][, t]))
  }, numeric(panel_n))
}

# Path log probabilities in unconstrained coordinates, in the order of
# coef() of a multilpa_lta fit: four means (profile 1 y1, y2, profile 2 y1,
# y2), four log variances, the initial logit of profile 1 against profile 2,
# then for each origin the logit of leaving against staying (intercept, and
# slope on z when `covariate`).
panel_unconstrained_paths <- function(phi, covariate) {
  means <- matrix(phi[1:4], 2L, 2L, byrow = TRUE)
  variances <- matrix(exp(phi[5:8]), 2L, 2L, byrow = TRUE)
  intercepts <- if (covariate) phi[c(10L, 12L)] else phi[10:11]
  slopes <- if (covariate) phi[c(11L, 13L)] else c(0, 0)
  initial <- c(phi[[9L]], 0)
  initial <- initial - log(sum(exp(initial)))
  log_move <- function(from, t) {
    eta <- intercepts[from] + slopes[from] * panel_z[, t]
    leave <- stats::plogis(eta, log.p = TRUE)
    stay <- stats::plogis(-eta, log.p = TRUE)
    if (from == 1L) cbind(stay, leave) else cbind(leave, stay)
  }
  panel_path_log(panel_log_density(means, variances),
                 matrix(initial, panel_n, 2L, byrow = TRUE), log_move)
}

panel_log_likelihood <- function(phi, covariate) {
  sum(rows_log_sum_exp(panel_unconstrained_paths(phi, covariate)))
}

# The maximum of that likelihood: two Newton steps from `start` with
# numDeriv's gradient and Hessian (quadratic convergence from a point already
# within about 1e-6 of it).
panel_maximum <- function(start, covariate) {
  value <- \(phi) panel_log_likelihood(phi, covariate)
  Reduce(function(phi, step) {
    phi - solve(numDeriv::hessian(value, phi), numDeriv::grad(value, phi))
  }, 1:2, init = unname(start))
}

# Posterior probability of each profile at each row (Bayes' rule over the
# paths), as get_results(, "posteriors") lays it out: row, profile, posterior.
panel_posteriors <- function(path_log) {
  weights <- exp(path_log - rows_log_sum_exp(path_log))
  cells <- expand.grid(t = 1:4, k = 1:2)
  do.call(rbind, lapply(seq_len(nrow(cells)), \(i) {
    t <- cells$t[i]
    k <- cells$k[i]
    data.frame(row = (seq_len(panel_n) - 1L) * 4L + t, profile = k,
               posterior = rowSums(weights[, panel_paths[, t] == k, drop = FALSE]))
  }))
}

# The homogeneous fit's natural parameters (coef() of a multilpa_transitions
# fit) as a function of the unconstrained coordinates, and back.
homogeneous_names <- c(
  sprintf("measurement.mean.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2")),
  sprintf("measurement.variance.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2")),
  "group.probability.group_class_1",
  sprintf("profile.initial_probability.profile_%d.group_class_1", 1:2),
  sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d",
          c(1L, 2L, 1L, 2L), c(1L, 1L, 2L, 2L)))
homogeneous_natural <- function(phi) {
  stats::setNames(c(phi[1:4], exp(phi[5:8]), 1, stats::plogis(phi[[9L]]), stats::plogis(-phi[[9L]]),
                    stats::plogis(-phi[[10L]]), stats::plogis(phi[[10L]]),
                    stats::plogis(phi[[11L]]), stats::plogis(-phi[[11L]])),
                  homogeneous_names)
}
homogeneous_unconstrained <- function(estimates) {
  value <- \(name) estimates[[name]]
  p <- \(to, from) value(sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d", to, from))
  c(unname(estimates[homogeneous_names[1:4]]), log(unname(estimates[homogeneous_names[5:8]])),
    log(value(homogeneous_names[10L]) / value(homogeneous_names[11L])),
    log(p(2, 1) / p(1, 1)), log(p(1, 2) / p(2, 2)))
}

covariate_names <- c(
  sprintf("mean.all.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2")),
  sprintf("log_variance.all.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2")),
  "initial.group_class_1.profile_1.(Intercept)",
  "transition.group_class_1.profile_1->profile_2.(Intercept)",
  "transition.group_class_1.profile_1->profile_2.z",
  "transition.group_class_1.profile_2->profile_1.(Intercept)",
  "transition.group_class_1.profile_2->profile_1.z")

# Profile means and variances as the profile tables lay them out.
profile_table <- function(means, variances, labels) {
  data.frame(profile = labels[rep(1:2, each = 2)], indicator = rep(c("y1", "y2"), 2),
             mean = as.vector(t(means)), variance = as.vector(t(variances)))
}

# Wald columns from an estimate and its standard error.
wald_columns <- function(estimate, standard_error) {
  statistic <- estimate / standard_error
  data.frame(standard_error = standard_error, statistic = statistic,
             conf_low = estimate - stats::qnorm(0.975) * standard_error,
             conf_high = estimate + stats::qnorm(0.975) * standard_error)
}

# A Wald table against its reference, in two comparisons: the standard
# errors relative; the z statistics and interval bounds scaled by
# max(1, |reference|), because a bound (or a statistic) can sit near zero,
# where a relative measure would judge the rounding of a number close to 0.
expect_wald <- function(actual, expected, keys, quantity, covers, tolerance) {
  expect_reference(actual, expected[c(keys, "standard_error")], reference = "numDeriv", covers = covers,
                   quantity = sprintf("%s: standard errors", quantity), tolerance = tolerance,
                   measure = "relative", key = keys)
  expect_reference(actual, expected[c(keys, "statistic", "conf_low", "conf_high")], reference = "numDeriv",
                   covers = covers, quantity = sprintf("%s: statistics and bounds", quantity),
                   tolerance = tolerance, measure = "scaled", key = keys)
}

# ---- multilpa_transitions ------------------------------------------------

test_that("homogeneous lta(): estimates, standard errors and intervals equal the path likelihood's", {
  skip_without_reference("numDeriv")
  estimates <- coef(homogeneous_fit)
  maximum <- panel_maximum(homogeneous_unconstrained(estimates), covariate = FALSE)
  natural <- homogeneous_natural(maximum)
  # latents stops at a relative change of 1e-12 in the likelihood; the
  # estimates then sit within about 1e-6 of the maximum Newton's method finds.
  expect_reference(estimates, natural, reference = "independent: maximum of the path likelihood, base R",
                   covers = c("lta()", "coef.multilpa_transitions()"),
                   quantity = "natural estimates, homogeneous", tolerance = 1e-6, measure = "relative",
                   note = "Newton steps with numDeriv's gradient and Hessian, started at latents' estimates.")
  expected_profiles <- profile_table(matrix(natural[1:4], 2L, byrow = TRUE),
                                     matrix(natural[5:8], 2L, byrow = TRUE), 1:2)
  expected_profiles$standard_deviation <- sqrt(expected_profiles$variance)
  expect_reference(as.data.frame(homogeneous_fit), expected_profiles,
                   reference = "independent: maximum of the path likelihood, base R",
                   covers = c("lta()", "as.data.frame.multilpa_transitions()"),
                   quantity = "profile means and variances, homogeneous", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))
  summary_fit <- summary(homogeneous_fit)
  expect_reference(as.data.frame(summary_fit), expected_profiles,
                   reference = "independent: maximum of the path likelihood, base R",
                   covers = c("summary.multilpa_transitions()", "as.data.frame.summary_multilpa_transitions()"),
                   quantity = "profile means and variances, summary(), homogeneous", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))
  expect_reference(get_results(summary_fit), expected_profiles,
                   reference = "independent: maximum of the path likelihood, base R",
                   covers = c("summary.multilpa_transitions()", "get_results.summary_multilpa_transitions()"),
                   quantity = "profile means and variances, get_results(summary()), homogeneous",
                   tolerance = 1e-6, measure = "relative", key = c("profile", "indicator"))

  # Natural-scale covariance by the delta method from the inverse of
  # numDeriv's Hessian in the unconstrained coordinates (the observed
  # information is invariant to the coordinates at the maximum).
  value <- \(phi) panel_log_likelihood(phi, covariate = FALSE)
  covariance <- solve(-numDeriv::hessian(value, maximum))
  jacobian <- numDeriv::jacobian(homogeneous_natural, maximum)
  errors <- stats::setNames(sqrt(pmax(diag(jacobian %*% covariance %*% t(jacobian)), 0)),
                            homogeneous_names)
  # 1e-5 relative: the package differentiates its analytic scores with a
  # step of 1e-4, and numDeriv's second derivatives are good to about 1e-7.
  expect_reference(as.data.frame(summary_fit),
                   data.frame(profile = expected_profiles$profile, indicator = expected_profiles$indicator,
                              mean_standard_error = unname(errors[homogeneous_names[1:4]]),
                              variance_standard_error = unname(errors[homogeneous_names[5:8]])),
                   reference = "numDeriv",
                   covers = c("summary.multilpa_transitions()", "as.data.frame.summary_multilpa_transitions()"),
                   quantity = "standard errors of the profile means and variances, homogeneous",
                   tolerance = 1e-5, measure = "relative", key = c("profile", "indicator"))

  # Wald inference on the natural scale: means symmetric; variances on the
  # log scale; probabilities on the logit scale (the delta-method error of
  # the transformed estimate, mapped back).
  inference <- subset(parameter_inference(homogeneous_fit), level != "group")
  names_used <- sprintf("%s.%s.%s.%s", inference$level, inference$parameter, inference$outcome, inference$term)
  estimate <- unname(natural[names_used])
  error <- unname(errors[names_used])
  z <- stats::qnorm(0.975)
  is_mean <- inference$parameter == "mean"
  is_variance <- inference$parameter == "variance"
  # One bound per row and side (-1 low, +1 high), on the parameter's scale.
  bound <- function(side) {
    vapply(seq_along(estimate), \(i) {
      e <- estimate[i]
      s <- error[i]
      if (is_mean[i]) e + side * z * s
      else if (is_variance[i]) e * exp(side * z * s / e)
      else stats::plogis(stats::qlogis(e) + side * z * s / (e * (1 - e)))
    }, numeric(1))
  }
  conf_low <- bound(-1)
  conf_high <- bound(1)
  statistic <- ifelse(is_mean, estimate / error, NA_real_)
  expected_inference <- data.frame(level = inference$level, outcome = inference$outcome,
                                   term = inference$term, parameter = inference$parameter,
                                   standard_error = error, statistic = statistic,
                                   conf_low = conf_low, conf_high = conf_high)
  # 1e-5, as for the summary's standard errors.
  expect_wald(inference, expected_inference, c("level", "outcome", "term", "parameter"),
              "Wald inference, natural scale, homogeneous",
              c("parameter_inference()", "parameter_inference.multilpa_transitions()"), 1e-5)
  # p-values of the means run down to 1e-130, where a 1e-6 difference in a z
  # of 24 is a large relative one: the difference is absolute.
  means_only <- subset(inference, parameter == "mean", select = c(level, outcome, term, parameter, p_value, p_adjusted))
  expect_reference(means_only,
                   data.frame(level = means_only$level, outcome = means_only$outcome, term = means_only$term,
                              parameter = means_only$parameter,
                              p_value = 2 * stats::pnorm(-abs(statistic[is_mean])),
                              p_adjusted = 2 * stats::pnorm(-abs(statistic[is_mean]))),
                   reference = "numDeriv",
                   covers = c("parameter_inference()", "parameter_inference.multilpa_transitions()"),
                   quantity = "z-test p-values of the means, homogeneous (no adjustment)",
                   tolerance = 1e-5, measure = "absolute", key = c("level", "outcome", "term", "parameter"))
  expect_reference(confint(homogeneous_fit)[names_used, ],
                   matrix(c(conf_low, conf_high), ncol = 2L,
                          dimnames = list(names_used, c("2.5 %", "97.5 %"))),
                   reference = "numDeriv", covers = c("confint.multilpa_transitions()", "lta()"),
                   quantity = "95% intervals, natural scale, homogeneous", tolerance = 1e-5,
                   measure = "scaled")
})

test_that("homogeneous lta(): nobs() is the number of sequences", {
  expect_reference(nobs(homogeneous_fit), length(unique(transition_panel$id)),
                   reference = "independent: number of sequences (distinct ids) in the data",
                   covers = c("lta()", "nobs.multilpa_transitions()"),
                   quantity = "number of observations, homogeneous", tolerance = 0)
})

# ---- multilpa_lta --------------------------------------------------------

test_that("covariate lta(): estimates, covariance and Wald tables equal the path likelihood's", {
  skip_without_reference("numDeriv")
  estimates <- coef(covariate_fit)
  maximum <- stats::setNames(panel_maximum(estimates[covariate_names], covariate = TRUE), covariate_names)
  # As above: latents' convergence at 1e-12 leaves the estimates within about
  # 1e-6 of the maximum.
  expect_reference(estimates, maximum, reference = "independent: maximum of the path likelihood, base R",
                   covers = c("lta()", "coef.multilpa_lta()"),
                   quantity = "unconstrained estimates, covariate transitions", tolerance = 1e-6,
                   measure = "relative",
                   note = "Newton steps with numDeriv's gradient and Hessian, started at latents' estimates.")
  expect_reference(get_results(covariate_fit, "profiles"),
                   profile_table(matrix(maximum[1:4], 2L, byrow = TRUE),
                                 matrix(exp(maximum[5:8]), 2L, byrow = TRUE), c("profile_1", "profile_2")),
                   reference = "independent: maximum of the path likelihood, base R",
                   covers = c("lta()", "get_results.multilpa_lta()"),
                   quantity = "profile means and variances, covariate transitions", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))

  # The transition table: the model-implied probabilities of each move,
  # averaged over the sequences at each occasion (z of the destination
  # occasion).
  intercepts <- maximum[c(10L, 12L)]
  slopes <- maximum[c(11L, 13L)]
  moves <- expand.grid(to = 1:2, from = 1:2, occasion = 2:4)
  moves$probability <- mapply(\(to, from, t) {
    leave <- stats::plogis(intercepts[[from]] + slopes[[from]] * panel_z[, t])
    mean(if (to == from) 1 - leave else leave)
  }, moves$to, moves$from, moves$occasion)
  expected_moves <- data.frame(group_class = "group_class_1", occasion = moves$occasion,
                               from = sprintf("profile_%d", moves$from), to = sprintf("profile_%d", moves$to),
                               probability = moves$probability)
  expect_reference(as.data.frame(covariate_fit), expected_moves,
                   reference = "independent: maximum of the path likelihood, base R",
                   covers = c("lta()", "as.data.frame.multilpa_lta()"),
                   quantity = "transition probabilities by occasion, covariate transitions",
                   tolerance = 1e-6, measure = "relative", key = c("group_class", "occasion", "from", "to"))
  summary_fit <- summary(covariate_fit)
  expect_reference(as.data.frame(summary_fit), expected_moves,
                   reference = "independent: maximum of the path likelihood, base R",
                   covers = c("summary.multilpa_lta()", "as.data.frame.summary_multilpa_lta()"),
                   quantity = "transition probabilities by occasion, summary(), covariate transitions",
                   tolerance = 1e-6, measure = "relative", key = c("group_class", "occasion", "from", "to"))

  value <- \(phi) panel_log_likelihood(phi, covariate = TRUE)
  covariance <- solve(-numDeriv::hessian(value, unname(maximum)))
  dimnames(covariance) <- list(covariate_names, covariate_names)
  # 1e-5 relative on the standard errors: numDeriv's second derivatives are
  # good to about 1e-7, the package's information is analytic. The
  # correlations rest on near-zero off-diagonal entries: absolute.
  expect_reference(sqrt(diag(vcov(covariate_fit))), sqrt(diag(covariance)), reference = "numDeriv",
                   covers = c("lta()", "vcov.multilpa_lta()"),
                   quantity = "standard errors, covariate transitions", tolerance = 1e-5, measure = "relative")
  expect_reference(stats::cov2cor(vcov(covariate_fit)), stats::cov2cor(covariance), reference = "numDeriv",
                   covers = c("lta()", "vcov.multilpa_lta()"),
                   quantity = "correlations of the estimates, covariate transitions", tolerance = 1e-5,
                   measure = "absolute")

  errors <- sqrt(diag(covariance))
  # Wald columns of a coefficient table, from the coefficient's name.
  expected_wald <- function(table, names_used, keys) {
    cbind(table[keys], wald_columns(unname(maximum[names_used]), unname(errors[names_used])))
  }
  coefficients <- get_results(covariate_fit, "transition_coefficients")
  coefficient_names <- sprintf("transition.%s.%s->%s.%s", coefficients$group_class, coefficients$from,
                               coefficients$to, coefficients$term)
  coefficient_keys <- c("group_class", "from", "to", "term")
  expect_wald(coefficients, expected_wald(coefficients, coefficient_names, coefficient_keys),
              coefficient_keys, "transition coefficients", c("lta()", "get_results.multilpa_lta()"), 1e-5)
  expect_wald(get_results(summary_fit, "transition_coefficients"),
              expected_wald(coefficients, coefficient_names, coefficient_keys), coefficient_keys,
              "transition coefficients, get_results(summary())",
              c("summary.multilpa_lta()", "get_results.summary_multilpa_lta()"), 1e-5)
  initial <- get_results(covariate_fit, "initial_coefficients")
  initial_names <- sprintf("initial.%s.%s.%s", initial$group_class, initial$profile, initial$term)
  expect_wald(initial, expected_wald(initial, initial_names, c("group_class", "profile", "term")),
              c("group_class", "profile", "term"), "initial coefficients",
              c("lta()", "get_results.multilpa_lta()"), 1e-5)
  inference <- parameter_inference(covariate_fit)
  inference_names <- ifelse(inference$block == "transition",
                            sprintf("transition.%s.%s->%s.%s", inference$group_class, inference$from,
                                    inference$to, inference$term),
                            sprintf("initial.%s.%s.%s", inference$group_class, inference$to, inference$term))
  inference_keys <- c("block", "group_class", "from", "to", "term")
  expected_inference <- expected_wald(inference, inference_names, inference_keys)
  expect_wald(inference, expected_inference, inference_keys, "Wald inference, covariate transitions",
              c("parameter_inference()", "parameter_inference.multilpa_lta()"), 1e-5)
  # p-values run down to 1e-12: absolute (no adjustment by default).
  p_values <- 2 * stats::pnorm(-abs(expected_inference$statistic))
  expect_reference(subset(inference, select = c(block, group_class, from, to, term, p_value, p_adjusted)),
                   cbind(inference[inference_keys], p_value = p_values, p_adjusted = p_values),
                   reference = "numDeriv",
                   covers = c("parameter_inference()", "parameter_inference.multilpa_lta()"),
                   quantity = "z-test p-values, covariate transitions (no adjustment)", tolerance = 1e-5,
                   measure = "absolute", key = inference_keys)
})

test_that("covariate lta(): posteriors are Bayes' rule over the paths, and nobs() counts the sequences", {
  # The posteriors at latents' own estimates: one evaluation, rounding only.
  # Probabilities, many near 0: absolute.
  expect_reference(get_results(covariate_fit, "posteriors"),
                   panel_posteriors(panel_unconstrained_paths(coef(covariate_fit)[covariate_names], covariate = TRUE)),
                   reference = "independent: posteriors by Bayes' rule over every profile path, base R",
                   covers = c("lta()", "get_results.multilpa_lta()"),
                   quantity = "row posteriors, covariate transitions", tolerance = 1e-10,
                   measure = "absolute", key = c("row", "profile"))
  expect_reference(nobs(covariate_fit), length(unique(transition_panel$id)),
                   reference = "independent: number of sequences (distinct ids) in the data",
                   covers = c("lta()", "nobs.multilpa_lta()"),
                   quantity = "number of observations, covariate transitions", tolerance = 0)
})

# ---- Transition networks -------------------------------------------------

test_that("get_tna() and get_group_tna() of an lta() fit carry its enumerated expected transitions", {
  skip_without_reference("tna")
  fit <- quietly(lta(transition_panel, c("y1", "y2"), "id", n_profiles = 2, time = "time",
                     n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-12))
  estimates <- coef(fit)
  value <- \(name) estimates[[name]]
  means <- matrix(vapply(sprintf("measurement.mean.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2")),
                         value, numeric(1)), 2L, byrow = TRUE)
  variances <- matrix(vapply(sprintf("measurement.variance.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2")),
                             value, numeric(1)), 2L, byrow = TRUE)
  log_density <- panel_log_density(means, variances)
  # Joint log probability of class h and every path, per sequence.
  joint <- lapply(1:2, \(h) {
    initial <- vapply(1:2, \(k) value(sprintf("profile.initial_probability.profile_%d.group_class_%d", k, h)),
                      numeric(1))
    move <- outer(1:2, 1:2, Vectorize(\(k, l) {
      value(sprintf("profile.transition_probability.profile_%d.group_class_%d:profile_%d", l, h, k))
    }))
    log(value(sprintf("group.probability.group_class_%d", h))) +
      panel_path_log(log_density, matrix(log(initial), panel_n, 2L, byrow = TRUE),
                     \(from, t) matrix(log(move[from, ]), panel_n, 2L, byrow = TRUE))
  })
  total <- rows_log_sum_exp(do.call(cbind, joint))
  weights <- lapply(joint, \(j) exp(j - total))
  counts <- lapply(weights, \(w) {
    Reduce(`+`, Map(\(p, moves) sum(w[, p]) * moves, seq_len(nrow(panel_paths)), panel_path_moves))
  })
  first <- lapply(weights, \(w) vapply(1:2, \(k) sum(w[, panel_paths[, 1L] == k]), numeric(1)))
  states <- c("profile_1", "profile_2")
  as_matrix <- \(m) matrix(m / rowSums(m), 2L, dimnames = list(states, states))

  network <- get_tna(fit)
  # Expected counts summed over the classes and normalised by row, at the
  # fit's own estimates: the same sums in another order.
  expect_reference(network$weights, as_matrix(Reduce(`+`, counts)),
                   reference = "independent: expected transition counts enumerated over every path, base R",
                   covers = c("get_tna()", "get_tna.multilpa_transitions()"),
                   quantity = "pooled transition network weights", tolerance = 1e-10, measure = "relative")
  # The initial distribution is the class-weighted average of the class
  # initial probabilities; at the maximum it equals the average posterior of
  # the first occasion (the EM fixed point), to the fit's convergence: 1e-6.
  expect_reference(network$inits, stats::setNames(Reduce(`+`, first) / panel_n, states),
                   reference = "independent: expected transition counts enumerated over every path, base R",
                   covers = c("get_tna()", "get_tna.multilpa_transitions()"),
                   quantity = "pooled initial distribution", tolerance = 1e-6, measure = "relative")
  groups <- get_group_tna(fit)
  # Each class's network is its estimated matrix; at the maximum that is the
  # class's expected counts normalised by row (the M step), to the fit's
  # convergence: 1e-6. The class initial distributions run down to 1e-8 (a
  # class that starts in one profile): absolute for those.
  actual <- do.call(rbind, lapply(1:2, \(h) {
    data.frame(class = h, from = rep(states, 2), to = rep(states, each = 2),
               weight = as.vector(groups[[h]]$weights))
  }))
  expected <- do.call(rbind, lapply(1:2, \(h) {
    data.frame(class = h, from = rep(states, 2), to = rep(states, each = 2),
               weight = as.vector(as_matrix(counts[[h]])))
  }))
  expect_reference(actual, expected,
                   reference = "independent: expected transition counts enumerated over every path, base R",
                   covers = c("get_group_tna()", "get_group_tna.multilpa_transitions()"),
                   quantity = "transition network weights per group class", tolerance = 1e-6,
                   measure = "relative", key = c("class", "from", "to"))
  expect_reference(
    do.call(rbind, lapply(1:2, \(h) data.frame(class = h, profile = states, initial = as.vector(groups[[h]]$inits)))),
    do.call(rbind, lapply(1:2, \(h) data.frame(class = h, profile = states, initial = first[[h]] / sum(first[[h]])))),
    reference = "independent: expected transition counts enumerated over every path, base R",
    covers = c("get_group_tna()", "get_group_tna.multilpa_transitions()"),
    quantity = "initial distribution per group class", tolerance = 1e-6, measure = "absolute",
    key = c("class", "profile"))
})

# Modal profile sequences of a two-level profile model, written from Bayes'
# rule: `log_share[[h]]` the log profile probabilities in group class h,
# `log_class` the log class probabilities (sequences x classes). Returns the
# modal profile of every row (sequences x occasions) and the modal class of
# every sequence.
modal_sequences <- function(log_density, log_share, log_class) {
  per_class <- lapply(1:2, \(h) {
    terms <- lapply(1:2, \(k) log_share[[h]][k] + log_density[[k]])
    row_total <- log(exp(terms[[1L]]) + exp(terms[[2L]]))
    list(terms = terms, row_total = row_total, group = log_class[, h] + rowSums(row_total))
  })
  group_total <- rows_log_sum_exp(cbind(per_class[[1L]]$group, per_class[[2L]]$group))
  class_posterior <- lapply(per_class, \(part) exp(part$group - group_total))
  profile_posterior <- lapply(1:2, \(k) {
    Reduce(`+`, lapply(1:2, \(h) class_posterior[[h]] * exp(per_class[[h]]$terms[[k]] - per_class[[h]]$row_total)))
  })
  list(states = ifelse(profile_posterior[[1L]] >= profile_posterior[[2L]], 1L, 2L),
       class = ifelse(class_posterior[[1L]] >= class_posterior[[2L]], 1L, 2L))
}

# tna's relative model of sequences: consecutive pairs counted and
# normalised by row; the initial distribution is the share of first states.
sequence_network <- function(states) {
  counts <- outer(1:2, 1:2, Vectorize(\(k, l) sum(states[, -4L] == k & states[, -1L] == l)))
  labels <- c("profile_1", "profile_2")
  list(weights = matrix(counts / rowSums(counts), 2L, dimnames = list(labels, labels)),
       inits = stats::setNames(tabulate(states[, 1L], 2L) / nrow(states), labels))
}

# The networks of a multilpa() fit, as one table: class (0 for the pooled
# network), from, to, weight; and the initial distributions.
network_table <- function(network, class) {
  data.frame(class = class, from = rep(rownames(network$weights), 2),
             to = rep(colnames(network$weights), each = 2), weight = as.vector(network$weights))
}

check_profile_networks <- function(fit, log_share, log_class, covers) {
  estimates <- coef(fit)
  means <- matrix(estimates[sprintf("measurement.mean.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2"))],
                  2L, byrow = TRUE)
  variances <- matrix(estimates[sprintf("measurement.variance.profile_%d.%s", rep(1:2, each = 2), c("y1", "y2"))],
                      2L, byrow = TRUE)
  modal <- modal_sequences(panel_log_density(means, variances), log_share, log_class)
  setting <- class(fit)[1L]
  # Ratios of integer counts: rounding only.
  expect_reference(network_table(get_tna(fit), 0L), network_table(sequence_network(modal$states), 0L),
                   reference = "independent: consecutive pairs of modal sequences from Bayes' rule, base R",
                   covers = c("get_tna()", covers[["pooled"]]),
                   quantity = sprintf("pooled network weights, %s", setting), tolerance = 1e-12,
                   measure = "relative", key = c("class", "from", "to"))
  expect_reference(get_tna(fit)$inits, sequence_network(modal$states)$inits,
                   reference = "independent: consecutive pairs of modal sequences from Bayes' rule, base R",
                   covers = c("get_tna()", covers[["pooled"]]),
                   quantity = sprintf("pooled initial distribution, %s", setting), tolerance = 1e-12,
                   measure = "relative")
  groups <- get_group_tna(fit)
  classes <- sort(unique(modal$class))
  expect_reference(
    do.call(rbind, lapply(seq_along(classes), \(i) network_table(groups[[i]], classes[i]))),
    do.call(rbind, lapply(classes, \(h) {
      network_table(sequence_network(modal$states[modal$class == h, , drop = FALSE]), h)
    })),
    reference = "independent: consecutive pairs of modal sequences from Bayes' rule, base R",
    covers = c("get_group_tna()", covers[["grouped"]]),
    quantity = sprintf("network weights per modal group class, %s", setting), tolerance = 1e-12,
    measure = "relative", key = c("class", "from", "to"))
}

test_that("get_tna() and get_group_tna() of a multilpa() fit count its modal sequences", {
  skip_without_reference("tna")
  fit <- quietly(multilpa(transition_panel, c("y1", "y2"), "id", n_profiles = 2, n_group_classes = 2,
                          time = "time", n_starts = 2, seed = 1, tol = 1e-10))
  estimates <- coef(fit)
  log_share <- lapply(1:2, \(h) {
    log(estimates[sprintf("profile.probability.profile_%d.group_class_%d", 1:2, h)])
  })
  log_class <- matrix(log(estimates[sprintf("group.probability.group_class_%d", 1:2)]),
                      panel_n, 2L, byrow = TRUE)
  check_profile_networks(fit, log_share, log_class,
                         c(pooled = "get_tna.multilpa()", grouped = "get_group_tna.multilpa()"))
})

test_that("get_tna() and get_group_tna() of a multilpa() fit with a group covariate count its modal sequences", {
  skip_without_reference("tna")
  fit <- quietly(multilpa(transition_panel, c("y1", "y2"), "id", n_profiles = 2, n_group_classes = 2,
                          time = "time", group_covariates = "w", n_starts = 2, seed = 1, tol = 1e-10))
  estimates <- coef(fit)
  # Profile 1 against profile 2 within each class; class 1 against class 2 on w.
  log_share <- lapply(1:2, \(h) {
    logits <- c(estimates[[sprintf("profile.coefficient.profile_1.group_class_%d", h)]], 0)
    logits - log(sum(exp(logits)))
  })
  eta <- estimates[["group.coefficient.group_class_1.(Intercept)"]] +
    estimates[["group.coefficient.group_class_1.w"]] * panel_w
  log_class <- cbind(stats::plogis(eta, log.p = TRUE), stats::plogis(-eta, log.p = TRUE))
  check_profile_networks(fit, log_share, log_class,
                         c(pooled = "get_tna.multilpa_covariates()",
                           grouped = "get_group_tna.multilpa_covariates()"))
})
