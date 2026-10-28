# Ordinal and count indicators, single level and in a latent transition
# model, against depmixS4 (Visser & Speekenbrink 2010, Journal of Statistical
# Software 36(7):1-21), MASS::glm.nb (Venables & Ripley 2002) and
# VGAM::vglm(family = acat()) (Yee 2015, Vector Generalized Linear and
# Additive Models, Springer).
#
# No suite compared the `ordinal` table (probability, location,
# location_standard_error) or the `count_means` table (mean, dispersion and
# their standard errors) before, and both paths have had real defects
# (17-digit ordinal labels, review CAT-F3; the NB dispersion at its floor,
# MR-F4). The model (`?multilpa`, R/extra-indicators.R) and its references:
#
#   Poisson counts with a Gaussian indicator. depmixS4 fits the same mixture
#   (`mix()`, one state per profile) and the same latent Markov model
#   (`depmix()`) by EM with `poisson()` and `gaussian()` responses, missing
#   responses integrated out. Its Poisson intercept is the log mean, its
#   Gaussian response an intercept and an SD:
#
#     depmixS4 pr / prior        <-> profile_probabilities / initial
#     depmixS4 transition (row)  <-> transitions (from)
#     depmixS4 poisson intercept <-> log(count_means$mean)
#     depmixS4 gaussian sd       <-> sqrt(profiles$variance)
#
#   Negative-binomial counts (NB2: variance mu + alpha mu^2). No R package
#   fits an NB2 mixture (flexmix 2.3.20 has no FLXMRnegbin), so the reference
#   is an EM whose M-step is MASS::glm.nb on posterior weights -- one
#   intercept-only fit per class for a varying dispersion, `y ~ 0 + class` on
#   the data stacked over classes for a shared one; glm.nb's theta is 1 /
#   alpha -- and whose E-step is stats::dnbinom (ordinal-count-reference.R).
#
#   Ordinal indicators (adjacent-category logit with shared intercepts and one
#   location per profile, the last profile's zero). This is exactly
#   VGAM::vglm(y ~ class, family = acat(parallel = TRUE)) with the last class
#   as the baseline: VGAM's eta_j = log P(j + 1) / P(j) = (a_(j+1) - a_j) +
#   eta_c, so its class coefficient is the location. The reference is an EM
#   with that M-step on posterior weights and the closed-form E-step.
#
# The two EM references start from the classes that generated the data, not
# from anything the package produced, and run to their fixed point (no
# parameter moving by 1e-10), so they are found independently: agreement of
# the maximised likelihood says both optimisers reached the same mode, and
# agreement of every output that the package computes the same tables from
# the same solution.
#
# For every case the likelihood function is checked on its own first: the
# reference likelihood (depmixS4's, or the base-R one) evaluated at the
# package's estimates must return the package's log-likelihood.
#
# Every output of every fit is then compared through compare_fit() (single
# level: profiles, covariances, responses, ordinal, count_means, class
# shares, posteriors, assignments, counts, the classification tables, BCH
# weights, entropy, information criteria, the model row, residuals of the
# `categorical` and `count` kinds, Gaussian diagnostics) or, for the latent
# transition fit, which lta() returns as the general `multilpa_lta` family,
# compare_fit_tables() with tables built in ordinal-count-reference.R.
# Standard errors are the inverse of a Richardson-extrapolated Hessian
# (numDeriv) of a likelihood written from scratch in base R, at the
# reference's estimates (`ordinal_count_standard_errors()`,
# `ordinal_count_transition_errors()`).

source(file.path("equivalence", "harness", "ordinal-count-reference.R"), local = environment())

suite_ordinal_count <- function() {
  require_suite_packages(
    "depmixS4", "MASS", "VGAM", "numDeriv",
    reason = paste("the Poisson mixture and latent Markov references, the NB2 and",
                   "adjacent-category M-steps, and the reference observed information"))
  # Standard errors: a Richardson Hessian of the base-R likelihood at the
  # reference's estimates against the package's own observed information at
  # its estimates; two numerical second derivatives agree to about 1e-5.
  error_tolerance <- c(`profiles$mean_standard_error` = 1e-4,
                       `profiles$variance_standard_error` = 1e-4,
                       `responses$probability_standard_error` = 1e-4,
                       `ordinal$location_standard_error` = 1e-4,
                       `count_means$standard_error` = 1e-4,
                       `count_means$dispersion_standard_error` = 1e-4)
  poisson <- .ordinal_count_poisson_data(n = 400L, seed = 501L)
  missing_cells <- .ordinal_count_poisson_data(n = 400L, seed = 501L, mcar = 0.1)
  negative_binomial <- .ordinal_count_nb_data(n = 400L, seed = 502L)
  ordinal <- .ordinal_count_ordinal_data(n = 500L, seed = 503L)
  stack_cases(list(
    .ordinal_count_depmix_case("poisson-gaussian", poisson, error_tolerance),
    .ordinal_count_depmix_case("poisson-gaussian-mcar10", missing_cells, error_tolerance),
    .ordinal_count_em_case("nb2-varying", negative_binomial, error_tolerance,
                           count = c("n1", "n2"), dispersion = "varying"),
    .ordinal_count_em_case("nb2-equal", negative_binomial, error_tolerance,
                           count = c("n1", "n2"), dispersion = "equal"),
    .ordinal_count_em_case("ordinal-binary", ordinal, error_tolerance,
                           ordinal = c("o1", "o2"), categorical = "b"),
    .ordinal_count_lta_case("lta-poisson-gaussian", n_groups = 120L, n_occasions = 4L,
                            seed = 504L)))
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' Two Poisson indicators and one Gaussian, two classes
#' @param n Rows.
#' @param seed Simulation seed.
#' @param mcar Share of cells of each indicator deleted completely at random.
#' @return A list with `data` (an `id` column and the indicators) and `class`,
#'   the generating class of each row.
.ordinal_count_poisson_data <- function(n, seed, mcar = 0) {
  set.seed(seed)
  class <- sample.int(2L, n, replace = TRUE, prob = c(0.6, 0.4))
  data <- data.frame(id = seq_len(n), c1 = stats::rpois(n, c(1, 4)[class]),
                     c2 = stats::rpois(n, c(3, 0.8)[class]),
                     y = stats::rnorm(n, c(0, 1.5)[class]))
  if (mcar > 0) {
    data[c("c1", "c2", "y")] <- lapply(data[c("c1", "c2", "y")], function(column) {
      column[sample.int(n, round(mcar * n))] <- NA
      column
    })
  }
  list(data = data, class = class)
}

#' Two NB2 indicators with dispersions 0.5 and 0.2, two classes
#' @inheritParams .ordinal_count_poisson_data
#' @return As for [.ordinal_count_poisson_data()].
.ordinal_count_nb_data <- function(n, seed) {
  set.seed(seed)
  class <- sample.int(2L, n, replace = TRUE, prob = c(0.6, 0.4))
  data <- data.frame(id = seq_len(n),
                     n1 = stats::rnbinom(n, mu = c(2, 6)[class], size = 1 / 0.5),
                     n2 = stats::rnbinom(n, mu = c(5, 1.5)[class], size = 1 / 0.2))
  list(data = data, class = class)
}

#' Two four-category ordinal indicators and one binary, two classes
#'
#' Categories are drawn from the adjacent-category model with intercepts
#' (0, 0.5, 0.2, -0.8) and locations of opposite sign in the two classes.
#'
#' @inheritParams .ordinal_count_poisson_data
#' @return As for [.ordinal_count_poisson_data()].
.ordinal_count_ordinal_data <- function(n, seed) {
  set.seed(seed)
  class <- sample.int(2L, n, replace = TRUE, prob = c(0.55, 0.45))
  draw <- function(location) {
    p <- exp(c(0, 0.5, 0.2, -0.8) + location * (0:3))
    sample.int(4L, 1L, prob = p / sum(p))
  }
  data <- data.frame(id = seq_len(n),
                     o1 = vapply(c(-0.6, 0.8)[class], draw, integer(1)),
                     o2 = vapply(c(0.7, -0.5)[class], draw, integer(1)),
                     b = factor(ifelse(stats::runif(n) < c(0.3, 0.75)[class], "yes", "no")))
  list(data = data, class = class)
}

#' Latent Markov sequences with two Poisson indicators and one Gaussian
#' @param n_groups Sequences.
#' @param n_occasions Occasions per sequence.
#' @param seed Simulation seed.
#' @return A data frame with `g`, `t`, `c1`, `c2` and `y`, in sequence and
#'   occasion order.
.ordinal_count_lta_data <- function(n_groups, n_occasions, seed) {
  set.seed(seed)
  transition <- matrix(c(0.8, 0.2, 0.25, 0.75), 2L, byrow = TRUE)
  frame <- do.call(rbind, lapply(seq_len(n_groups), function(group) {
    path <- unlist(Reduce(function(previous, step) {
      sample.int(2L, 1L, prob = transition[previous, ])
    }, seq_len(n_occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE))
    data.frame(g = group, t = seq_len(n_occasions), state = path)
  }))
  frame$c1 <- stats::rpois(nrow(frame), c(1, 4)[frame$state])
  frame$c2 <- stats::rpois(nrow(frame), c(3, 0.8)[frame$state])
  frame$y <- stats::rnorm(nrow(frame), c(0, 1.5)[frame$state])
  frame[c("g", "t", "c1", "c2", "y")]
}

# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

#' Poisson and Gaussian mixture against depmixS4's mix()
#' @param dataset Label.
#' @param fixture From [.ordinal_count_poisson_data()].
#' @param error_tolerance Tolerances of the standard-error columns.
#' @return A `harness_case`.
.ordinal_count_depmix_case <- function(dataset, fixture, error_tolerance) {
  data <- fixture$data
  vars <- c("c1", "c2", "y")
  layout <- ordinal_count_layout(data[vars], gaussian = "y", count = c("c1", "c2"))
  incomplete <- anyNA(data[vars])
  # The package's fits run to tol = 1e-15, as the mclust suite's do. At
  # 1e-12 (a relative change in the log-likelihood) its NB2 fit stopped with
  # a largest likelihood gradient of 3e-5 and posteriors 1e-6 from the fixed
  # point, which a default tolerance of 1e-6 then reports; at 1e-15 the
  # gradient is 7e-7. That is the package's stopping rule, not its solution.
  fit <- multilpa(data, vars = vars, id = "id", n_profiles = 2L, n_group_classes = 1L,
                  count = c("c1", "c2"), missing = if (incomplete) "fiml" else "error",
                  n_starts = 5L, seed = 1L, tol = 1e-15, max_iter = 20000L)
  model <- depmixS4::mix(list(c1 ~ 1, c2 ~ 1, y ~ 1), data = data, nstates = 2L,
                         family = list(stats::poisson(), stats::poisson(), stats::gaussian()))
  label <- sprintf("%s, 2 profiles", dataset)

  # Claim one: the same likelihood function, depmixS4's evaluated at the
  # package's estimates.
  ours <- c(fit$profile_probabilities[1L, ],
            unlist(lapply(seq_len(2L), function(profile) {
              c(log(fit$count_means[profile, c("c1", "c2")]), fit$means[profile, "y"],
                sqrt(fit$variances[profile, "y"]))
            }), use.names = FALSE))
  at_ours <- as.numeric(.oc_depmix_quiet(depmixS4::logLik(depmixS4::setpars(model, ours))))
  # Claim two: the same maximum, depmixS4 searching from its own start.
  set.seed(1L)
  theirs <- .oc_depmix_quiet(depmixS4::fit(model, verbose = FALSE, emcontrol =
    depmixS4::em.control(maxit = 5000L, tol = 1e-12)))
  maximum <- as.numeric(.oc_depmix_quiet(depmixS4::logLik(theirs)))
  claims <- rbind(
    compare_values(sprintf("%s: depmixS4 log likelihood at the package's estimates", label),
                   at_ours, fit$log_likelihood, tolerance = 1e-10, scale = "relative"),
    compare_values(sprintf("%s: maximised log likelihood, depmixS4's own search", label),
                   maximum, fit$log_likelihood, tolerance = 1e-6, scale = "relative"))

  case <- if (abs(maximum - fit$log_likelihood) > 1e-4) {
    .ordinal_count_other_mode(label, fit, claims, "depmixS4")
  } else {
    # depmixS4 stops on a relative change of 1e-12 in the log-likelihood,
    # which leaves its estimates up to 1e-6 short of EM's fixed point here
    # (on a fixture like this one, one more step moved a Poisson intercept by
    # 1.1e-6), so its own EM is
    # continued to the fixed point first.
    refined <- ordinal_count_depmix_refine(theirs)
    pars <- depmixS4::getpars(refined)
    blocks <- matrix(pars[-(1:2)], 2L, byrow = TRUE)
    estimates <- c(list(pro = pars[1:2]), ordinal_count_depmix_responses(blocks, layout))
    log_likelihood <- as.numeric(.oc_depmix_quiet(depmixS4::logLik(refined)))
    posteriors <- as.matrix(depmixS4::posterior(refined, type = "smoothing"))
    e <- ordinal_count_e_step(estimates, layout)
    stopifnot(
      "the base-R likelihood must reproduce depmixS4's" =
        abs(e$log_likelihood - log_likelihood) < 1e-9 * abs(log_likelihood),
      "the base-R posteriors must reproduce depmixS4's" =
        max(abs(e$posteriors - posteriors)) < 1e-10)
    reference <- harness_reference(
      log_likelihood = log_likelihood, n_parameters = depmixS4::freepars(refined),
      posteriors = posteriors, means = estimates$means, variances = estimates$variances,
      standard_errors = ordinal_count_standard_errors(estimates, layout, log_likelihood),
      profile_probabilities = matrix(estimates$pro, 1L), group_probabilities = 1,
      data = data[vars], count_means = estimates$count_means)
    compare_fit(label, fit, reference, tolerance = c(default = 1e-6, error_tolerance),
                extra = claims)
  }
  case$comparisons <- transform(
    case$comparisons, source = sprintf("depmixS4 %s", utils::packageVersion("depmixS4")),
    dataset = dataset, precision = "machine")
  case
}

#' NB2 or ordinal mixture against the EM with other programs' M-steps
#' @param dataset Label.
#' @param fixture A fixture list.
#' @param error_tolerance Tolerances of the standard-error columns.
#' @param categorical,ordinal,count Indicator names of each kind.
#' @param dispersion `"varying"` or `"equal"` for NB2 counts.
#' @return A `harness_case`.
.ordinal_count_em_case <- function(dataset, fixture, error_tolerance,
                                   categorical = character(), ordinal = character(),
                                   count = character(), dispersion = "varying") {
  data <- fixture$data
  # Categorical indicators are listed before ordinal ones. The package pairs
  # discrete indicators for its residuals as categorical then ordinal
  # (R/bivariate-residuals.R, `.multilpa_discrete_items()`), while
  # reference_tables() orders the pairs as the fit's `vars`; the two
  # orientations name the same pair differently, and only agree when `vars`
  # lists the categorical ones first. The statistics are symmetric in the pair.
  vars <- c(categorical, ordinal, count)
  negative_binomial <- length(count) > 0L
  layout <- ordinal_count_layout(data[vars], categorical = categorical, ordinal = ordinal,
                                 count = count,
                                 dispersion = if (negative_binomial) dispersion else "none")
  # tol = 1e-15: see .ordinal_count_depmix_case().
  fit <- multilpa(data, vars = vars, id = "id", n_profiles = 2L, n_group_classes = 1L,
                  categorical = categorical, ordinal = ordinal, count = count,
                  count_model = if (negative_binomial) "negative_binomial" else "poisson",
                  count_dispersion = dispersion, n_starts = 5L, seed = 1L, tol = 1e-15,
                  max_iter = 20000L)
  label <- sprintf("%s, 2 profiles", dataset)

  # Claim one: the same likelihood function, the base-R one evaluated at the
  # package's estimates.
  ours <- list(pro = fit$profile_probabilities[1L, ],
               responses = fit$response_probabilities[categorical],
               ordinal = stats::setNames(lapply(ordinal, function(item) {
                 .oc_ordinal_probabilities(fit$ordinal_intercepts[[item]],
                                           fit$ordinal_locations[, item])
               }), ordinal),
               count_means = fit$count_means,
               # A shared dispersion may be stored once; the likelihood
               # takes one row per profile.
               count_dispersion = if (!is.null(fit$count_dispersion))
                 fit$count_dispersion[rep_len(seq_len(nrow(fit$count_dispersion)), 2L), ,
                                      drop = FALSE])
  at_ours <- ordinal_count_e_step(ours, layout)$log_likelihood

  # Claim two: the same maximum, found by the reference EM started from the
  # generating classes (softened, so that no M-step sees a zero weight).
  start <- 0.1 + 0.8 * diag(2L)[fixture$class, ]
  em <- ordinal_count_em(start, layout)
  m_steps <- if (negative_binomial) sprintf("MASS %s glm.nb", utils::packageVersion("MASS")) else
    sprintf("VGAM %s vglm(acat)", utils::packageVersion("VGAM"))
  claims <- rbind(
    compare_values(sprintf("%s: base-R log likelihood at the package's estimates", label),
                   at_ours, fit$log_likelihood, tolerance = 1e-10, scale = "relative"),
    compare_values(sprintf("%s: maximised log likelihood, reference EM from the generating classes (%d steps)",
                           label, em$steps),
                   em$log_likelihood, fit$log_likelihood, tolerance = 1e-6, scale = "relative"))

  case <- if (abs(em$log_likelihood - fit$log_likelihood) > 1e-4) {
    .ordinal_count_other_mode(label, fit, claims, "the reference EM")
  } else {
    # A location's standard error depends on which profile is the baseline
    # (the last, whose location is zero), so the reference is put in the
    # fit's profile order, matched on the posteriors, before its errors are
    # taken; compare_fit()'s own alignment is then the identity.
    aligned <- ordinal_count_permute(em$estimates, em$posteriors, align_classes(
      t(fit$subject_posteriors), t(em$posteriors)))
    estimates <- aligned$estimates
    reference <- harness_reference(
      log_likelihood = em$log_likelihood,
      n_parameters = ordinal_count_n_parameters(estimates, layout),
      posteriors = aligned$posteriors, responses = estimates$responses,
      standard_errors = ordinal_count_standard_errors(estimates, layout, em$log_likelihood),
      profile_probabilities = matrix(estimates$pro, 1L), group_probabilities = 1,
      data = data[vars], ordinal = estimates$ordinal, count_means = estimates$count_means,
      count_dispersion = estimates$count_dispersion)
    compare_fit(label, fit, reference, tolerance = c(default = 1e-6, error_tolerance),
                extra = claims)
  }
  case$comparisons <- transform(
    case$comparisons, source = sprintf("EM with %s M-step", m_steps), dataset = dataset,
    precision = "machine")
  case
}

#' Different modes: only the claims are comparable
#' @param label Case label.
#' @param fit The package's fit.
#' @param claims The two claim rows.
#' @param reference_name Who found the other mode.
#' @return A `harness_case` whose ledger exempts every table with the reason.
.ordinal_count_other_mode <- function(label, fit, claims, reference_name) {
  inventory <- fit_outputs(fit)
  harness_case(claims, output_ledger(label, inventory, claims, exempt = stats::setNames(
    rep(sprintf("%s stopped at a different mode, so its estimates describe another solution",
                reference_name), length(unique(inventory$table))), unique(inventory$table))))
}

#' Latent transitions with Poisson and Gaussian indicators against depmixS4
#' @param dataset Label.
#' @param n_groups,n_occasions,seed Passed to [.ordinal_count_lta_data()].
#' @return A `harness_case`.
.ordinal_count_lta_case <- function(dataset, n_groups, n_occasions, seed) {
  frame <- .ordinal_count_lta_data(n_groups, n_occasions, seed)
  vars <- c("c1", "c2", "y")
  k <- 2L
  layout <- ordinal_count_layout(frame[vars], gaussian = "y", count = c("c1", "c2"))
  # tol = 1e-15: see .ordinal_count_depmix_case().
  fit <- lta(frame, vars, "g", n_profiles = k, time = "t", count = c("c1", "c2"),
             n_starts = 5L, seed = seed, tol = 1e-15, max_iter = 20000L)
  ntimes <- as.integer(table(frame$g)[as.character(unique(frame$g))])
  model <- depmixS4::depmix(list(c1 ~ 1, c2 ~ 1, y ~ 1), data = frame, nstates = k,
                            family = list(stats::poisson(), stats::poisson(),
                                          stats::gaussian()), ntimes = ntimes)
  label <- sprintf("%s, %d states", dataset, k)

  # The package's estimates, read from its tables, in depmixS4's layout.
  tables <- get_results(fit, "all")
  profile_of <- function(labels) as.integer(sub("^profile_", "", labels))
  transition <- matrix(0, k, k)
  transition[cbind(profile_of(tables$transitions$from), profile_of(tables$transitions$to))] <-
    tables$transitions$probability
  initial <- tables$initial$probability[order(profile_of(tables$initial$profile))]
  count_means <- matrix(0, k, 2L, dimnames = list(NULL, c("c1", "c2")))
  count_means[cbind(tables$count_means$profile,
                    match(tables$count_means$indicator, c("c1", "c2")))] <- tables$count_means$mean
  gaussian <- tables$profiles[order(profile_of(tables$profiles$profile)), , drop = FALSE]
  ours <- c(initial, as.vector(t(transition)),
            as.vector(rbind(log(count_means[, "c1"]), log(count_means[, "c2"]),
                            gaussian$mean, sqrt(gaussian$variance))))
  log_likelihood_ours <- tables$fit$log_likelihood

  # Claim one: the same likelihood function. Claim two: the same maximum.
  at_ours <- as.numeric(depmixS4::logLik(depmixS4::setpars(model, ours)))
  set.seed(seed)
  theirs <- depmixS4::fit(model, verbose = FALSE, emcontrol =
    depmixS4::em.control(maxit = 5000L, tol = 1e-12))
  maximum <- as.numeric(depmixS4::logLik(theirs))
  claims <- rbind(
    compare_values(sprintf("%s: depmixS4 log likelihood at the package's estimates", label),
                   at_ours, log_likelihood_ours, tolerance = 1e-10, scale = "relative"),
    compare_values(sprintf("%s: maximised log likelihood, depmixS4's own search", label),
                   maximum, log_likelihood_ours, tolerance = 1e-6, scale = "relative"))

  case <- if (abs(maximum - log_likelihood_ours) > 1e-4) {
    inventory <- fit_outputs(fit)
    harness_case(claims, output_ledger(label, inventory, claims, exempt = stats::setNames(
      rep("depmixS4 stopped at a different mode, so its estimates describe another solution",
          length(unique(inventory$table))), unique(inventory$table))))
  } else {
    refined <- ordinal_count_depmix_refine(theirs)
    pars <- depmixS4::getpars(refined)
    blocks <- matrix(pars[-seq_len(k + k^2)], k, byrow = TRUE)
    emission <- ordinal_count_depmix_responses(blocks, layout)
    # States are arbitrary in both programs: depmixS4's are put in the fit's
    # order by matching the Gaussian means and log count means.
    order_states <- align_classes(cbind(gaussian$mean, log(count_means)),
                                  cbind(emission$means, log(emission$count_means)))
    emission <- lapply(emission, function(block) block[order_states, , drop = FALSE])
    log_likelihood <- as.numeric(depmixS4::logLik(refined))
    reference <- list(
      initial = pars[seq_len(k)][order_states],
      transition = matrix(pars[k + seq_len(k^2)], k, k, byrow = TRUE)[order_states, order_states],
      emission = emission, log_likelihood = log_likelihood,
      n_parameters = depmixS4::freepars(refined),
      posteriors = depmixS4::forwardbackward(refined)$gamma[, order_states, drop = FALSE])
    reference$errors <- ordinal_count_transition_errors(
      reference$initial, reference$transition, emission, layout, ntimes, log_likelihood)
    # The variance floor the package applies by default: a millionth of the
    # smallest indicator variance (`?lta`, `min_variance`). It only sets the
    # Gaussian diagnostics' boundary flag.
    expected <- ordinal_count_lta_tables(reference, layout, frame, ntimes,
                                         min_variance = 1e-6 * stats::var(frame$y))
    # The Wald columns are the logit over its standard error, so they inherit
    # the standard error's tolerance (z up to 7 here: 1e-4 relative).
    wald <- c("standard_error", "statistic", "p_value", "conf_low", "conf_high")
    compare_fit_tables(label, fit, expected,
                       tolerance = c(default = 1e-6, stats::setNames(
                         rep(1e-4, 2L * length(wald)),
                         c(paste0("transition_coefficients$", wald),
                           paste0("initial_coefficients$", wald)))),
                       extra = claims)
  }
  case$comparisons <- transform(
    case$comparisons, source = sprintf("depmixS4 %s", utils::packageVersion("depmixS4")),
    dataset = dataset, precision = "machine")
  case
}
