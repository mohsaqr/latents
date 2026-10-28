# The general latent transition engine (class `multilpa_lta`) against
# depmixS4 (Visser & Speekenbrink 2010, Journal of Statistical Software
# 36(7):1-21).
#
# The depmixs4 and hard-depmixs4 suites compare the homogeneous engine
# (`multilpa_transitions`). `transition_covariates`, `initial_covariates` and
# `transitions = "occasion"` move a fit to the general engine, whose tables
# (transition and initial logit coefficients with Wald columns, model-implied
# transition probabilities per occasion) no other suite compares. depmixS4
# fits the same likelihood with multinomial-logit transition and prior models,
# so its parameters map one to one:
#
#   depmixS4 prior = ~ w (mlogit, base: state 1)
#       <-> initial_coefficients (logit of each profile against the last)
#       b_k = B[, k] - B[, K]
#   depmixS4 transition = ~ x_lead (mlogit per origin state i, base: state 1)
#       <-> transition_coefficients (logit of moving to j against staying in i)
#       a_{i -> j} = C_i[, j] - C_i[, i]
#   depmixS4 gaussian (Intercept), sd <-> profiles mean, sqrt(variance)
#   depmixS4 multinomial mlogit on a binary factor (base: category "0")
#       <-> responses probability, P("1") = plogis(coefficient)
#
# Which occasion's covariate governs a move. The package reads the covariate
# at the destination occasion (`?lta`, `transition_covariates`). depmixS4
# stores the transition of the move t -> t + 1 on row t (`trDens[t, to,
# from]`), so its covariate is the destination's value written on the origin
# row (`x_lead`); the last row of each sequence governs no move and holds 0.
# The unshifted mapping gives a log likelihood 11.9 lower (probe p07 of the
# 2026-10-27 review). `transitions = "occasion"` estimates one logit per move
# and destination occasion (terms `occasion_2`, ...), which is depmixS4's
# `transition = ~ 0 + o2 + o3 + o4` with `o_t` the indicator that the move
# written on that row ends at occasion t: the same cell-mean coding, so no
# column is unidentified (a factor of the unshifted occasion would leave the
# last occasion's level on unused rows only).
#
# Every case makes three claims, so that a disagreement says what is wrong:
#   (A) computation at one point. depmixS4 is set to the package's estimates
#       (`setpars()`), and every function of the parameters is read from it
#       there: the log likelihood and criteria, the smoothed posteriors
#       (`forwardbackward()$gamma`), the model-implied transition
#       probabilities averaged over the groups at each occasion (`trDens`),
#       the initial distribution averaged over groups (`init`), the Gaussian
#       diagnostics, and the observed-information standard errors -- a
#       Richardson Hessian (numDeriv) of depmixS4's own `logLik()` over its
#       free parameters, mapped to the package's logits by the delta method
#       (the map is linear in the logits, so the delta method is exact for
#       them). depmixS4's `standardError()` is not used: its `nlme::fdHess`
#       is 1-15 percent off (depmixs4 suite) and gave 9.72 against 10.90 in
#       probe p11. The Wald columns (statistic, p value, interval) and the odds
#       ratios are written here from those estimates and errors.
#   (B) the same maximum. The free parameters are compared with the maximum
#       of depmixS4's likelihood, reached from the package's solution, in
#       units of their standard errors (1e-3). depmixS4's own EM is not used
#       to get there: its multinomial-logit M-step is an inner `nnet`
#       optimisation, and on a draft of the occasion case 300 steps of it
#       (`.depmixs4_refine()`) from its own converged solution were still
#       moving parameters by more than 1e-12 (6.3 s, measured 2026-10-27),
#       where the package's solution is already at the maximum to 1e-6 of a
#       standard error. Newton steps on
#       depmixS4's `logLik()` (numDeriv gradient, the Hessian of claim A) are
#       taken instead until a step moves no parameter by more than 1e-5 of
#       its standard error, as the Mplus suite carries Mplus's estimates to
#       the maximum by Newton steps (see `.lta_general_newton()` for why the
#       rule is in standard errors).
#       The log likelihood that gains is a claim too.
#   (C) search. The package's own starts (5, tol 1e-12) against the best of
#       five depmixS4 random-start searches (`fit()`, `random.start = TRUE`).
#
# Two datasets, both simulated here in base R with fixed seeds:
#   covariates  N = 150 sequences of 5 occasions, two Gaussian indicators,
#               K = 2, a transition covariate x varying over occasions and an
#               initial covariate w (the generator of probes p07/p08, whose
#               depmixS4 maximum was -2324.06073573 with 14 parameters);
#   occasion    N = 200 sequences of 4 occasions, three binary items, K = 2,
#               `transitions = "occasion"`, no covariates.

suite_lta_general <- function() {
  require_suite_packages(
    "depmixS4", "numDeriv",
    reason = paste("the independent hidden-Markov reference with multinomial-logit",
                   "transition and prior models, and the reference observed information"))
  stack_cases(list(.lta_general_case(.lta_general_covariates()),
                   .lta_general_case(.lta_general_occasion())))
}

# ---------------------------------------------------------------------------
# Datasets and their depmixS4 models
# ---------------------------------------------------------------------------

#' The covariate case: probes p07/p08's generator
#'
#' Draws in the order the probe's loop drew them (initial state, then one
#' uniform per move, sequence by sequence), so the data are the probe's.
#' Transition: staying in profile 1 has logit 1.5 + 0.8 x, in profile 2
#' 1.5 - 0.5 x, at the destination occasion's x; initial: profile 1 with
#' probability plogis(0.3 + w).
#' @return A case specification: data, the lta() settings, the depmixS4
#'   model, and the term correspondence.
#' @noRd
.lta_general_covariates <- function() {
  set.seed(5)
  n <- 150L
  n_occasions <- 5L
  id <- rep(seq_len(n), each = n_occasions)
  occ <- rep(seq_len(n_occasions), n)
  x <- stats::rnorm(n * n_occasions)
  w <- rep(stats::rnorm(n), each = n_occasions)
  states <- unlist(lapply(seq_len(n), function(i) {
    offset <- (i - 1L) * n_occasions
    first <- if (stats::runif(1) < stats::plogis(0.3 + w[offset + 1L])) 1L else 2L
    Reduce(function(previous, step) {
      stay <- stats::plogis(1.5 + 0.8 * x[offset + step] * (previous == 1L) -
                              0.5 * x[offset + step] * (previous == 2L))
      if (stats::runif(1) < stay) previous else 3L - previous
    }, seq_len(n_occasions)[-1L], init = first, accumulate = TRUE)
  }))
  y1 <- stats::rnorm(n * n_occasions, c(0, 2)[states])
  y2 <- stats::rnorm(n * n_occasions, c(1, -1)[states], c(1, 0.7)[states])
  frame <- data.frame(id, occ, x, w, y1, y2)
  # The destination's covariate on the origin row; the last row of a sequence
  # governs no move.
  frame$x_lead <- stats::ave(frame$x, frame$id, FUN = function(v) c(v[-1L], 0))
  list(dataset = "covariates", frame = frame, k = 2L, ntimes = rep(n_occasions, n),
       continuous = c("y1", "y2"), categorical = character(),
       transitions = "homogeneous", transition_covariates = "x", initial_covariates = "w",
       model = function(frame, ntimes) depmixS4::depmix(
         list(y1 ~ 1, y2 ~ 1), data = frame, nstates = 2L,
         family = list(stats::gaussian(), stats::gaussian()),
         transition = ~ x_lead, prior = ~ w,
         initdata = frame[frame$occ == 1L, , drop = FALSE], ntimes = ntimes),
       # package term -> depmixS4 term, in depmixS4's column order
       initial_terms = c(`(Intercept)` = "(Intercept)", w = "w"),
       transition_terms = c(`(Intercept)` = "(Intercept)", x = "x_lead"))
}

#' The occasion case: occasion-specific transitions, binary items
#'
#' Staying has probability 0.85, 0.75 and 0.6 on the moves into occasions 2,
#' 3 and 4; profile 1 starts with probability 0.6; the items are endorsed
#' with probabilities (0.8, 0.7, 0.75) in profile 1 and (0.2, 0.3, 0.25) in
#' profile 2.
#' @return A case specification, as [.lta_general_covariates()].
#' @noRd
.lta_general_occasion <- function() {
  set.seed(6)
  n <- 200L
  n_occasions <- 4L
  frame <- data.frame(id = rep(seq_len(n), each = n_occasions),
                      occ = rep(seq_len(n_occasions), n))
  stay <- c(0.85, 0.75, 0.6)
  states <- unlist(lapply(seq_len(n), function(i) {
    first <- if (stats::runif(1) < 0.6) 1L else 2L
    Reduce(function(previous, step) {
      if (stats::runif(1) < stay[[step - 1L]]) previous else 3L - previous
    }, seq_len(n_occasions)[-1L], init = first, accumulate = TRUE)
  }))
  endorse <- rbind(c(0.8, 0.7, 0.75), c(0.2, 0.3, 0.25))
  items <- c("u1", "u2", "u3")
  frame[items] <- lapply(seq_along(items), function(j) {
    factor(stats::rbinom(nrow(frame), 1L, endorse[states, j]), levels = c(0L, 1L))
  })
  # o_t: the move written on this row ends at occasion t (cell-mean coding of
  # the destination occasion; the last row of a sequence governs no move).
  lead <- stats::ave(frame$occ, frame$id, FUN = function(v) c(v[-1L], NA))
  frame[sprintf("o%d", seq_len(n_occasions)[-1L])] <- lapply(
    seq_len(n_occasions)[-1L], function(t) as.numeric(!is.na(lead) & lead == t))
  # A column of ones makes the prior a multinomial logit with an intercept, as
  # the package parameterises it; `prior = ~ 1` would make it a probability
  # vector under a sum constraint.
  frame$one <- 1
  list(dataset = "occasion", frame = frame, k = 2L, ntimes = rep(n_occasions, n),
       continuous = character(), categorical = items,
       transitions = "occasion", transition_covariates = character(),
       initial_covariates = character(),
       model = function(frame, ntimes) depmixS4::depmix(
         list(u1 ~ 1, u2 ~ 1, u3 ~ 1), data = frame, nstates = 2L,
         family = rep(list(depmixS4::multinomial("mlogit")), 3L),
         transition = ~ 0 + o2 + o3 + o4, prior = ~ 0 + one,
         initdata = frame[frame$occ == 1L, , drop = FALSE], ntimes = ntimes),
       initial_terms = c(`(Intercept)` = "one"),
       transition_terms = c(occasion_2 = "o2", occasion_3 = "o3", occasion_4 = "o4"))
}

# ---------------------------------------------------------------------------
# The parameter correspondence
# ---------------------------------------------------------------------------

#' The package's parameters, as one depmixS4 parameter vector
#'
#' depmixS4's `getpars()` order: the prior's terms-by-states coefficients by
#' term, each origin state's transition coefficients likewise, then for each
#' state every response in formula order (Gaussian intercept and SD; the
#' category coefficients of a multinomial logit). The base state (1) and the
#' base category hold zero.
#' @param tables The fit's tables, from `get_results(fit, "all")`.
#' @param case A case specification.
#' @param levels Category labels of each categorical indicator, in the order
#'   depmixS4 codes them.
#' @return A numeric vector of `getpars()`'s length.
#' @noRd
.lta_general_depmix_pars <- function(tables, case, levels) {
  k <- case$k
  profiles <- sprintf("profile_%d", seq_len(k))
  initial <- tables$initial_coefficients
  logits <- matrix(0, length(case$initial_terms), k,
                   dimnames = list(names(case$initial_terms), profiles))
  logits[cbind(initial$term, initial$profile)] <- initial$estimate
  prior <- logits - logits[, 1L]
  moves <- tables$transition_coefficients
  transition <- lapply(profiles, function(origin) {
    out <- matrix(0, length(case$transition_terms), k,
                  dimnames = list(names(case$transition_terms), profiles))
    rows <- moves[moves$from == origin, , drop = FALSE]
    out[cbind(rows$term, rows$to)] <- rows$estimate
    out - out[, 1L]
  })
  responses <- lapply(profiles, function(profile) {
    gaussian <- unlist(lapply(case$continuous, function(indicator) {
      row <- tables$profiles[tables$profiles$profile == profile &
                               tables$profiles$indicator == indicator, , drop = FALSE]
      stopifnot("one profile row per indicator" = nrow(row) == 1L)
      c(row$mean, sqrt(row$variance))
    }))
    categorical <- unlist(lapply(case$categorical, function(indicator) {
      rows <- tables$responses[tables$responses$profile == profile &
                                 tables$responses$indicator == indicator, , drop = FALSE]
      probability <- rows$probability[match(levels[[indicator]], rows$category)]
      stopifnot("one response row per category" = !anyNA(probability))
      log(probability / probability[[1L]])
    }))
    c(gaussian, categorical)
  })
  c(as.vector(t(prior)), unlist(lapply(transition, function(block) as.vector(t(block)))),
    unlist(responses))
}

#' depmixS4's parameters, as the package reports them
#'
#' The inverse of [.lta_general_depmix_pars()]: logits re-expressed against
#' the package's references (the last profile for the initial distribution,
#' staying for a transition), SDs squared, category coefficients turned into
#' probabilities. Keys: `initial|profile|term`, `transition|from|to|term`,
#' `mean|profile|indicator`, `variance|profile|indicator`,
#' `probability|profile|indicator|category`.
#' @param pars A depmixS4 parameter vector.
#' @param case A case specification.
#' @param levels Category labels of each categorical indicator.
#' @return A named numeric vector.
#' @noRd
.lta_general_package_values <- function(pars, case, levels) {
  k <- case$k
  profiles <- sprintf("profile_%d", seq_len(k))
  n_initial <- length(case$initial_terms)
  n_transition <- length(case$transition_terms)
  block <- function(values, terms) {
    matrix(values, length(terms), k, byrow = TRUE, dimnames = list(names(terms), profiles))
  }
  prior <- block(pars[seq_len(n_initial * k)], case$initial_terms)
  used <- n_initial * k
  transition <- lapply(seq_len(k), function(i) {
    block(pars[used + (i - 1L) * n_transition * k + seq_len(n_transition * k)],
          case$transition_terms)
  })
  used <- used + k * n_transition * k
  initial <- unlist(lapply(profiles[-k], function(profile) {
    stats::setNames(prior[, profile] - prior[, k],
                    sprintf("initial|%s|%s", profile, rownames(prior)))
  }))
  moves <- unlist(lapply(seq_len(k), function(i) {
    unlist(lapply(setdiff(seq_len(k), i), function(j) {
      stats::setNames(transition[[i]][, j] - transition[[i]][, i],
                      sprintf("transition|%s|%s|%s", profiles[[i]], profiles[[j]],
                              rownames(transition[[i]])))
    }))
  }))
  widths <- c(rep(2L, length(case$continuous)),
              vapply(case$categorical, function(item) length(levels[[item]]), integer(1)))
  per_state <- sum(widths)
  measurement <- unlist(lapply(seq_len(k), function(state) {
    values <- pars[used + (state - 1L) * per_state + seq_len(per_state)]
    ends <- cumsum(widths)
    pieces <- lapply(seq_along(widths), function(r) {
      values[ends[[r]] - widths[[r]] + seq_len(widths[[r]])]
    })
    gaussian <- unlist(lapply(seq_along(case$continuous), function(r) {
      indicator <- case$continuous[[r]]
      stats::setNames(c(pieces[[r]][[1L]], pieces[[r]][[2L]]^2),
                      sprintf(c("mean|%s|%s", "variance|%s|%s"), profiles[[state]], indicator))
    }))
    categorical <- unlist(lapply(seq_along(case$categorical), function(c_index) {
      indicator <- case$categorical[[c_index]]
      coefficients <- pieces[[length(case$continuous) + c_index]]
      probability <- exp(coefficients - max(coefficients))
      stats::setNames(probability / sum(probability),
                      sprintf("probability|%s|%s|%s", profiles[[state]], indicator,
                              levels[[indicator]]))
    }))
    c(gaussian, categorical)
  }))
  c(initial, moves, measurement)
}

#' The key of each row of an output table, in [.lta_general_package_values()]'s
#' naming
#' @noRd
.lta_general_keys <- function(table, name, column = NULL) {
  switch(name,
    initial_coefficients = sprintf("initial|%s|%s", table$profile, table$term),
    transition_coefficients = sprintf("transition|%s|%s|%s", table$from, table$to, table$term),
    profiles = sprintf("%s|%s|%s", column, table$profile, table$indicator),
    responses = sprintf("probability|%s|%s|%s", table$profile, table$indicator,
                        table$category))
}

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

#' Compare one dataset's general transition fit with depmixS4
#' @param case A case specification.
#' @return A `harness_case`.
#' @noRd
.lta_general_case <- function(case) {
  frame <- case$frame
  vars <- c(case$continuous, case$categorical)
  label <- sprintf("%s, %d states", case$dataset, case$k)
  # Warnings are not suppressed: the harness counts them and records the first.
  fit <- lta(frame, vars = vars, id = "id", n_profiles = case$k, time = "occ",
             categorical = case$categorical, transitions = case$transitions,
             transition_covariates = case$transition_covariates,
             initial_covariates = case$initial_covariates,
             n_starts = 5L, seed = 1L, tol = 1e-12)
  tables <- get_results(fit, "all")
  levels <- lapply(stats::setNames(nm = case$categorical), function(item) levels(frame[[item]]))
  model <- case$model(frame, case$ntimes)
  pars_template <- depmixS4::getpars(model)
  free <- !depmixS4::getpars(model, "fixed")
  # The layout assumed above, checked against depmixS4's own names for the
  # prior and transition blocks.
  expected_names <- c(
    as.vector(t(outer(case$initial_terms, sprintf(".St%d", seq_len(case$k)), paste0))),
    rep(as.vector(t(outer(case$transition_terms, sprintf(".St%d", seq_len(case$k)), paste0))),
        case$k))
  stopifnot("depmixS4's prior and transition blocks must be laid out as assumed" =
              identical(names(pars_template)[seq_along(expected_names)], expected_names))

  at_package <- .lta_general_depmix_pars(tables, case, levels)
  stopifnot("the package's parameters must fill depmixS4's vector" =
              length(at_package) == length(pars_template),
            "depmixS4's fixed (base) coefficients must be zero" =
              all(abs(at_package[!free]) < 1e-14))
  fill <- function(theta) {
    out <- at_package
    out[free] <- theta
    out
  }
  log_likelihood <- function(theta) {
    as.numeric(depmixS4::logLik(depmixS4::setpars(model, fill(theta))))
  }
  values <- function(theta) .lta_general_package_values(fill(theta), case, levels)
  theta_package <- at_package[free]

  # (A) The shared point. The Hessian is depmixS4's likelihood's, over its
  # free parameters; at a maximum the inverse information transforms with the
  # Jacobian, here of the map to the package's reported parameters.
  shared <- depmixS4::setpars(model, at_package)
  information <- -numDeriv::hessian(log_likelihood, theta_package)
  covariance <- solve(information)
  jacobian <- numDeriv::jacobian(values, theta_package)
  errors <- stats::setNames(sqrt(diag(jacobian %*% covariance %*% t(jacobian))),
                            names(values(theta_package)))
  identified <- all(eigen(information, symmetric = TRUE, only.values = TRUE)$values > 0)

  # (B) The maximum of depmixS4's likelihood, from the package's solution.
  maximum <- .lta_general_newton(log_likelihood, theta_package, covariance)
  at_maximum <- values(maximum$theta)

  expected <- .lta_general_expected(shared, fit, case, values(theta_package), errors,
                                    identified)
  distances <- .lta_general_distances(label, tables, at_maximum, errors)
  claims <- rbind(
    compare_values(sprintf(paste("%s: log likelihood depmixS4's likelihood gains from the",
                                 "package's solution to its maximum (%d Newton steps)"),
                           label, maximum$steps),
                   0, maximum$log_likelihood - log_likelihood(theta_package),
                   tolerance = 1e-8),
    .lta_general_search(label, model, fit))
  case_result <- compare_fit_tables(
    label, fit, expected,
    tolerance = c(
      # The likelihood and the criteria built from it, evaluated by both
      # programs at the same parameters: rounding only.
      `fit$log_likelihood` = 1e-8, `fit$aic` = 1e-8, `fit$bic` = 1e-8,
      `fit$bic_individual` = 1e-8,
      # Observed-information standard errors and the Wald columns built from
      # them: the numerical-Hessian budget of the depmixs4 suite.
      `transition_coefficients$standard_error` = 1e-4,
      `transition_coefficients$statistic` = 1e-4,
      `transition_coefficients$p_value` = 1e-4,
      `transition_coefficients$conf_low` = 1e-4,
      `transition_coefficients$conf_high` = 1e-4,
      `initial_coefficients$standard_error` = 1e-4,
      `initial_coefficients$statistic` = 1e-4,
      `initial_coefficients$p_value` = 1e-4,
      `initial_coefficients$conf_low` = 1e-4,
      `initial_coefficients$conf_high` = 1e-4),
    extra = rbind(distances, .with_output_columns(claims)))
  case_result$comparisons <- transform(
    case_result$comparisons,
    source = sprintf("depmixS4 %s", utils::packageVersion("depmixS4")),
    dataset = case$dataset, precision = "machine")
  case_result
}

#' Newton steps on depmixS4's log likelihood to its maximum
#'
#' The information is held at its value at the start (the package's
#' solution), so each step is `covariance %*% gradient`; started within a
#' small fraction of a standard error, that converges geometrically. Steps are
#' taken until none moves a parameter by more than 1e-5 of its standard
#' error, two orders inside the 1e-3 the parameters are held to. An absolute
#' rule cannot be met: a Richardson gradient of a log likelihood of order 1e3
#' is resolved only to about 1e-5 (on the covariates case two step settings
#' differ by 2.5e-6 at the same point, measured 2026-10-27), so successive
#' steps from the maximum wander by 3e-8 to 2e-7 in the parameters, 5e-7 to
#' 1.1e-6 of a standard error, without settling.
#' @param log_likelihood Function of the free parameters.
#' @param theta Starting free parameters.
#' @param covariance Inverse observed information at `theta`.
#' @param max_steps Bound on the steps.
#' @return A list with `theta`, `log_likelihood` and `steps`; a classed warning
#'   when the bound is reached.
#' @noRd
.lta_general_newton <- function(log_likelihood, theta, covariance, max_steps = 10L) {
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    step <- as.vector(covariance %*% numDeriv::grad(log_likelihood, state$theta))
    list(theta = state$theta + step,
         settled = max(abs(step) / sqrt(diag(covariance))) < 1e-5, steps = index)
  }, seq_len(max_steps), init = list(theta = theta, settled = FALSE, steps = 0L))
  if (!state$settled) {
    warning(warningCondition(sprintf(
      "Newton steps on depmixS4's likelihood were still moving after %d steps", max_steps),
      class = "multilpa_harness_no_converge"))
  }
  list(theta = state$theta, log_likelihood = log_likelihood(state$theta),
       steps = state$steps)
}

#' Every table the fit reports, from depmixS4 at the package's estimates
#'
#' The free-parameter columns (coefficient estimates, means, variances,
#' response probabilities) are left out: they are compared with depmixS4's
#' maximum in standard errors ([.lta_general_distances()]). Their tables keep
#' their keys, so the rows-found check still runs.
#' @param shared depmixS4's model set to the package's estimates.
#' @param fit The fit, for its settings (`min_variance`) only.
#' @param case A case specification.
#' @param at_package [.lta_general_package_values()] at the shared point.
#' @param errors Reference standard errors, keyed the same way.
#' @param identified Whether depmixS4's information is positive definite.
#' @return A named list of expected tables.
#' @noRd
.lta_general_expected <- function(shared, fit, case, at_package, errors, identified) {
  frame <- case$frame
  k <- case$k
  profiles <- sprintf("profile_%d", seq_len(k))
  passes <- depmixS4::forwardbackward(shared)
  posteriors <- passes$gamma
  log_likelihood <- as.numeric(depmixS4::logLik(shared))
  n_parameters <- depmixS4::freepars(shared)
  z <- stats::qnorm(0.975)
  coefficient_table <- function(prefix, columns) {
    keys <- grep(sprintf("^%s[|]", prefix), names(at_package), value = TRUE)
    pieces <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
    estimate <- unname(at_package[keys])
    error <- unname(errors[keys])
    statistic <- estimate / error
    data.frame(group_class = "group_class_1",
               stats::setNames(as.data.frame(pieces[, -1L, drop = FALSE],
                                             stringsAsFactors = FALSE), columns),
               standard_error = error, statistic = statistic,
               p_value = 2 * stats::pnorm(-abs(statistic)),
               conf_low = estimate - z * error, conf_high = estimate + z * error,
               odds_ratio = exp(estimate), estimated = identified,
               stringsAsFactors = FALSE)
  }
  # trDens[row, to, from] is the move written on `row`: from that row's
  # occasion into the next. The package's table averages each occasion's moves
  # over the groups that make them.
  lead <- stats::ave(frame$occ, frame$id, FUN = function(v) c(v[-1L], NA))
  destinations <- sort(unique(stats::na.omit(lead)))
  cells <- expand.grid(from = seq_len(k), to = seq_len(k))
  transitions <- do.call(rbind, lapply(destinations, function(t) {
    rows <- which(!is.na(lead) & lead == t)
    data.frame(group_class = "group_class_1", occasion = t,
               from = profiles[cells$from], to = profiles[cells$to],
               probability = vapply(seq_len(nrow(cells)), function(cell) {
                 mean(shared@trDens[rows, cells$to[[cell]], cells$from[[cell]]])
               }, numeric(1)),
               estimated = identified, stringsAsFactors = FALSE)
  }))
  n_groups <- length(case$ntimes)
  n_observations <- sum(case$ntimes)
  tables <- list(
    transitions = transitions,
    transition_coefficients = coefficient_table("transition", c("from", "to", "term")),
    initial = data.frame(group_class = "group_class_1", profile = profiles,
                         probability = colMeans(shared@init), stringsAsFactors = FALSE),
    initial_coefficients = coefficient_table("initial", c("profile", "term")),
    # One chain: every sequence in it.
    group_classes = data.frame(group_class = "group_class_1", weight = 1,
                               count = n_groups, stringsAsFactors = FALSE),
    assignments = data.frame(
      .row = seq_len(nrow(frame)), occ = frame$occ,
      profile = max.col(posteriors, ties.method = "first"), group_class = 1L,
      uncertainty = 1 - apply(posteriors, 1L, max),
      stats::setNames(as.data.frame(posteriors), sprintf("posterior_profile_%d", seq_len(k)))),
    fit = data.frame(n_profiles = k, n_group_classes = 1L,
                     log_likelihood = log_likelihood, n_parameters = n_parameters,
                     aic = -2 * log_likelihood + 2 * n_parameters,
                     bic = -2 * log_likelihood + n_parameters * log(n_groups),
                     bic_individual = -2 * log_likelihood + n_parameters * log(n_observations),
                     n_groups = n_groups, n_observations = n_observations))
  if (length(case$continuous) > 0L) {
    variances <- matrix(at_package[sprintf("variance|%s|%s", rep(profiles, length(case$continuous)),
                                           rep(case$continuous, each = k))],
                        k, length(case$continuous), dimnames = list(NULL, case$continuous))
    tables$profiles <- data.frame(occasion = NA_integer_,
                                  profile = rep(profiles, length(case$continuous)),
                                  indicator = rep(case$continuous, each = k),
                                  stringsAsFactors = FALSE)
    tables$gaussian_diagnostics <- .reference_gaussian_diagnostics(
      list(data = frame[case$continuous], posteriors = posteriors, variances = variances),
      case$continuous, fit$min_variance)
  }
  if (length(case$categorical) > 0L) {
    keys <- grep("^probability[|]", names(at_package), value = TRUE)
    pieces <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
    tables$responses <- data.frame(occasion = NA_integer_, profile = pieces[, 2L],
                                   indicator = pieces[, 3L], category = pieces[, 4L],
                                   stringsAsFactors = FALSE)
  }
  tables
}

#' The free parameters against depmixS4's maximum, in standard errors
#'
#' Each reported parameter's distance from the maximum of depmixS4's
#' likelihood, divided by its reference standard error, is held to 1e-3: at a
#' thousandth of a standard error no estimate, interval or test computed from
#' the fit can differ (the hard-depmixs4 suite's rule).
#' @param tables The fit's tables.
#' @param at_maximum [.lta_general_package_values()] at depmixS4's maximum.
#' @param errors Reference standard errors.
#' @return Comparison rows, one per output column, with `output` set.
#' @noRd
.lta_general_distances <- function(label, tables, at_maximum, errors) {
  columns <- list(transition_coefficients = "estimate", initial_coefficients = "estimate",
                  profiles = c("mean", "variance"), responses = "probability")
  columns <- columns[intersect(names(columns), names(Filter(function(table) {
    is.data.frame(table) && nrow(table) > 0L
  }, tables)))]
  rows <- lapply(names(columns), function(name) {
    table <- tables[[name]]
    lapply(columns[[name]], function(column) {
      keys <- .lta_general_keys(table, name, column)
      stopifnot("every parameter cell must have a depmixS4 counterpart" =
                  all(keys %in% names(at_maximum)))
      compare_block(
        sprintf("%s: %s$%s, distance from depmixS4's maximum in standard errors",
                label, name, column),
        rep(0, length(keys)), abs(table[[column]] - at_maximum[keys]) / errors[keys],
        tolerance = 1e-3, output = paste0(name, "$", column))
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

#' Claim C: the package's own search against depmixS4's best of five
#' @return Comparison rows, none of them an output of the fit.
#' @noRd
.lta_general_search <- function(label, model, fit) {
  searched <- vapply(seq_len(5L), function(start) {
    set.seed(700L + start)
    result <- depmixS4::fit(model, verbose = FALSE, emcontrol = depmixS4::em.control(
      maxit = 5000, tol = 1e-12, random.start = TRUE))
    as.numeric(depmixS4::logLik(result))
  }, numeric(1))
  best <- max(searched)
  rows <- list(compare_values(
    sprintf(paste("%s: shortfall of the package's own search against depmixS4 (best of 5",
                  "random-start depmixS4 searches, reached by %d)"),
            label, sum(searched >= best - 1e-6)),
    0, max(best - fit$log_likelihood, 0), tolerance = 1e-6))
  excess <- fit$log_likelihood - best
  if (excess > 1e-6) {
    rows <- c(rows, list(compare_values(
      sprintf("%s: package maximum above depmixS4's best (by %.3g; depmixS4's search)",
              label, excess), excess, excess, tolerance = 1e-6)))
  }
  do.call(rbind, rows)
}
