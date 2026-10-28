# Ordinal and count indicators, written from their definitions in base R.
#
# Sourced by suites/ordinal-count.R. Three things are here, none of which
# calls the package:
#
#   1. A single-level mixture likelihood over Gaussian, categorical, ordinal
#      and count indicators, and a first-order latent Markov (forward)
#      likelihood over the same indicators. Each indicator is conditionally
#      independent given the class, and a missing value contributes a factor
#      of one:
#
#        Gaussian     dnorm(y, mu_c, sqrt(var_c))
#        categorical  rho_c(y)
#        ordinal      adjacent-category logit (`?multilpa`, `ordinal`):
#                     log P(m | c) = a_m + eta_c (m - 1) - log sum_l exp(a_l +
#                     eta_c (l - 1)), a_1 = 0 and the last class's eta = 0
#        Poisson      dpois(y, mu_c)
#        NB2          dnbinom(y, size = 1 / alpha, mu = mu_c): mean mu, variance
#                     mu + alpha mu^2, alpha per class or shared
#
#      Its numerical Hessian (numDeriv, Richardson extrapolation) at a
#      reference's estimates gives the observed-information standard errors
#      in the natural parameters.
#
#   2. An EM whose M-step is done by other programs' maximum-likelihood
#      fitters, given posterior weights: MASS::glm.nb for the NB2 means and
#      dispersions (one intercept-only fit per class when the dispersion
#      varies, one `y ~ 0 + class` fit on the data stacked over classes when it
#      is shared), VGAM::vglm(family = acat(parallel = TRUE)) for an ordinal
#      indicator (`y ~ class` on the posterior-weighted category table, the
#      last class the baseline: VGAM's eta_j = log P(j + 1) / P(j) is the
#      package's a_(j+1) - a_j + eta_c), and the closed forms for the class
#      shares, Gaussian moments and category proportions. Its E-step is the
#      base-R likelihood above. Run to its fixed point it is a reference
#      maximum that owes nothing to the package except where it starts.
#
#   3. The tables of a general latent transition fit (`multilpa_lta`), built
#      from a depmixS4 solution and the forward likelihood.

# ---------------------------------------------------------------------------
# Layout and densities
# ---------------------------------------------------------------------------

#' The indicators of a fixture, coded for the likelihood
#'
#' Category labels are read from the data: a factor's levels, or an ordinal
#' indicator's sorted distinct values, which is how the package labels them.
#'
#' @param data The indicator data frame.
#' @param gaussian,categorical,ordinal,count Indicator names of each kind.
#' @param dispersion `"none"` (Poisson counts), `"varying"` or `"equal"` (NB2).
#' @param k Number of classes.
#' @return A list describing the indicators.
ordinal_count_layout <- function(data, gaussian = character(), categorical = character(),
                                 ordinal = character(), count = character(),
                                 dispersion = c("none", "varying", "equal"), k = 2L) {
  dispersion <- match.arg(dispersion)
  levels_of <- function(name) {
    value <- data[[name]]
    if (is.factor(value)) levels(value) else as.character(sort(unique(value[!is.na(value)])))
  }
  code <- function(names) {
    stats::setNames(lapply(names, function(name) {
      match(as.character(data[[name]]), levels_of(name))
    }), names)
  }
  list(data = data, k = k, gaussian = gaussian, categorical = categorical,
       ordinal = ordinal, count = count, dispersion = dispersion,
       levels = stats::setNames(lapply(c(categorical, ordinal), levels_of),
                                c(categorical, ordinal)),
       codes = code(c(categorical, ordinal)))
}

#' log f(y_i | class) summed over indicators, rows by classes
#'
#' @param estimates A list with any of `means`, `variances` (classes by
#'   Gaussian indicators), `responses`, `ordinal` (named lists of classes by
#'   categories probability matrices), `count_means` and `count_dispersion`
#'   (classes by count indicators).
#' @param layout From [ordinal_count_layout()].
#' @return A numeric matrix, one row per data row and one column per class.
ordinal_count_log_density <- function(estimates, layout) {
  data <- layout$data
  k <- layout$k
  per_class <- function(values, density) {
    observed <- !is.na(values)
    out <- matrix(0, length(values), k)
    out[observed, ] <- vapply(seq_len(k), function(class) density(values[observed], class),
                              numeric(sum(observed)))
    out
  }
  gaussian <- lapply(layout$gaussian, function(name) {
    per_class(data[[name]], function(y, class) {
      stats::dnorm(y, estimates$means[class, name], sqrt(estimates$variances[class, name]),
                   log = TRUE)
    })
  })
  categories <- lapply(c(layout$categorical, layout$ordinal), function(name) {
    block <- c(estimates$responses, estimates$ordinal)[[name]]
    per_class(layout$codes[[name]], function(code, class) log(block[class, code]))
  })
  counts <- lapply(layout$count, function(name) {
    per_class(data[[name]], function(y, class) {
      mu <- estimates$count_means[class, name]
      if (identical(layout$dispersion, "none")) stats::dpois(y, mu, log = TRUE) else
        stats::dnbinom(y, size = 1 / estimates$count_dispersion[class, name], mu = mu,
                       log = TRUE)
    })
  })
  Reduce(`+`, c(gaussian, categories, counts), matrix(0, nrow(data), k))
}

#' Log-sum-exp over the columns of each row
#' @param m A numeric matrix.
#' @return One value per row.
.oc_log_sum_exp <- function(m) {
  peak <- apply(m, 1L, max)
  peak + log(rowSums(exp(m - peak)))
}

#' The E-step of the single-level mixture
#' @param estimates As for [ordinal_count_log_density()], plus `pro`.
#' @param layout From [ordinal_count_layout()].
#' @return A list with `posteriors` (rows by classes), `row_log_likelihood`
#'   and `log_likelihood`.
ordinal_count_e_step <- function(estimates, layout) {
  joint <- sweep(ordinal_count_log_density(estimates, layout), 2L, log(estimates$pro), "+")
  rows <- .oc_log_sum_exp(joint)
  list(posteriors = exp(joint - rows), row_log_likelihood = rows,
       log_likelihood = sum(rows))
}

# ---------------------------------------------------------------------------
# Ordinal parameters
# ---------------------------------------------------------------------------

#' Category probabilities from intercepts and locations
#' @param intercepts `a_2, ..., a_M`.
#' @param locations One location per class, the last zero.
#' @return Classes by categories.
.oc_ordinal_probabilities <- function(intercepts, locations) {
  scores <- seq_len(length(intercepts) + 1L) - 1
  logits <- outer(locations, scores) +
    matrix(c(0, intercepts), length(locations), length(scores), byrow = TRUE)
  exp(logits - .oc_log_sum_exp(logits))
}

#' Intercepts and locations from category probabilities
#'
#' Under the model `log P(m | c) - log P(m | K) = eta_c (m - 1) + const_c`, so
#' the location is the slope of that difference across categories and the
#' intercepts are the last class's log odds against the first category.
#'
#' @param block Classes by categories, the last class the baseline.
#' @return A list with `intercepts` (`M - 1`) and `locations` (`K - 1`).
.oc_ordinal_parameters <- function(block) {
  k <- nrow(block)
  m <- ncol(block)
  against_last <- log(block) - matrix(log(block[k, ]), k, m, byrow = TRUE)
  locations <- (against_last[, m] - against_last[, 1L]) / (m - 1L)
  list(intercepts = unname(log(block[k, -1L]) - log(block[k, 1L])),
       locations = unname(locations[-k]))
}

# ---------------------------------------------------------------------------
# Free natural parameters
# ---------------------------------------------------------------------------

#' The free natural parameters of a set of estimates, as a nested list
#'
#' Class shares and category probabilities drop their last element (one
#' minus the sum of the rest); an ordinal indicator is its intercepts and its
#' first `K - 1` locations; a shared dispersion appears once per indicator.
#' `utils::relist()` inverts [unlist()] of this list.
#'
#' @param estimates As for [ordinal_count_e_step()]; `pro` may be absent (a
#'   transition model's emission parameters).
#' @param layout From [ordinal_count_layout()].
#' @return A named list of numeric vectors and matrices.
ordinal_count_free <- function(estimates, layout) {
  k <- layout$k
  drop_last <- function(block) block[, -ncol(block), drop = FALSE]
  free <- list(
    pro = if (!is.null(estimates$pro)) estimates$pro[-k],
    means = if (length(layout$gaussian) > 0L) unname(estimates$means[, layout$gaussian, drop = FALSE]),
    variances = if (length(layout$gaussian) > 0L)
      unname(estimates$variances[, layout$gaussian, drop = FALSE]),
    responses = if (length(layout$categorical) > 0L)
      lapply(estimates$responses[layout$categorical], function(block) unname(drop_last(block))),
    ordinal = if (length(layout$ordinal) > 0L)
      lapply(estimates$ordinal[layout$ordinal], .oc_ordinal_parameters),
    count_means = if (length(layout$count) > 0L)
      unname(estimates$count_means[, layout$count, drop = FALSE]),
    count_dispersion = switch(layout$dispersion,
      none = NULL,
      varying = unname(estimates$count_dispersion[, layout$count, drop = FALSE]),
      equal = unname(estimates$count_dispersion[1L, layout$count])))
  Filter(Negate(is.null), free)
}

#' Estimates from free natural parameters
#' @param free A list shaped like [ordinal_count_free()]'s.
#' @param layout From [ordinal_count_layout()].
#' @return Estimates as [ordinal_count_e_step()] takes them.
ordinal_count_estimates <- function(free, layout) {
  k <- layout$k
  complete <- function(block, labels) {
    block <- matrix(block, k)
    out <- cbind(block, 1 - rowSums(block))
    colnames(out) <- labels
    out
  }
  named <- function(block, names) matrix(block, k, length(names), dimnames = list(NULL, names))
  list(
    pro = if (!is.null(free$pro)) c(free$pro, 1 - sum(free$pro)),
    means = if (!is.null(free$means)) named(free$means, layout$gaussian),
    variances = if (!is.null(free$variances)) named(free$variances, layout$gaussian),
    responses = if (!is.null(free$responses))
      stats::setNames(lapply(layout$categorical, function(name) {
        complete(free$responses[[name]], layout$levels[[name]])
      }), layout$categorical),
    ordinal = if (!is.null(free$ordinal))
      stats::setNames(lapply(layout$ordinal, function(name) {
        block <- .oc_ordinal_probabilities(free$ordinal[[name]]$intercepts,
                                           c(free$ordinal[[name]]$locations, 0))
        colnames(block) <- layout$levels[[name]]
        block
      }), layout$ordinal),
    count_means = if (!is.null(free$count_means)) named(free$count_means, layout$count),
    count_dispersion = switch(layout$dispersion,
      none = NULL,
      varying = named(free$count_dispersion, layout$count),
      equal = named(rep(free$count_dispersion, each = k), layout$count)))
}

#' Number of free parameters of a single-level mixture
#' @param estimates,layout As for [ordinal_count_free()].
#' @return A whole number.
ordinal_count_n_parameters <- function(estimates, layout) {
  length(unlist(ordinal_count_free(estimates, layout)))
}

#' Inverse observed information of a single-level mixture, from scratch
#'
#' The likelihood is [ordinal_count_e_step()] in the free natural parameters
#' of [ordinal_count_free()]; its Hessian is taken by Richardson extrapolation
#' at the reference's estimates, after the likelihood has been shown to
#' reproduce the reference's maximum. At a maximum the inverse observed
#' information transforms with the Jacobian of any reparameterisation, so the
#' errors are those of the natural parameters on whatever scale the package
#' works.
#'
#' @param estimates,layout As for [ordinal_count_free()].
#' @param log_likelihood The reference's maximum.
#' @return A list as `harness_reference(standard_errors = )` takes it.
ordinal_count_standard_errors <- function(estimates, layout, log_likelihood) {
  skeleton <- ordinal_count_free(estimates, layout)
  theta <- unlist(skeleton)
  at <- function(values) {
    ordinal_count_e_step(ordinal_count_estimates(utils::relist(values, skeleton), layout),
                         layout)$log_likelihood
  }
  stopifnot("the base-R likelihood must reproduce the reference's maximum" =
              abs(at(theta) - log_likelihood) < 1e-8 * abs(log_likelihood))
  covariance <- solve(-numDeriv::hessian(at, theta))
  .oc_errors(covariance, skeleton, layout)
}

#' Standard errors of every natural parameter, from a covariance of the free ones
#' @param covariance Covariance of `unlist(skeleton)`.
#' @param skeleton From [ordinal_count_free()].
#' @param layout From [ordinal_count_layout()].
#' @return A list as `harness_reference(standard_errors = )` takes it.
.oc_errors <- function(covariance, skeleton, layout) {
  k <- layout$k
  # Position of every free parameter in unlist(skeleton), shaped like it.
  index <- utils::relist(seq_along(unlist(skeleton)), skeleton)
  diagonal <- sqrt(diag(covariance))
  error_of <- function(positions, names) {
    matrix(diagonal[positions], k, length(names), dimnames = list(NULL, names))
  }
  out <- list()
  if (!is.null(index$means)) {
    out$mean <- error_of(index$means, layout$gaussian)
    out$variance <- error_of(index$variances, layout$gaussian)
  }
  if (!is.null(index$responses)) {
    # The last category is one minus the others: its variance is the sum of
    # the free block's covariance.
    out$probability <- stats::setNames(lapply(layout$categorical, function(name) {
      free <- matrix(index$responses[[name]], k)
      block <- t(vapply(seq_len(k), function(class) {
        cells <- covariance[free[class, ], free[class, ], drop = FALSE]
        c(sqrt(diag(cells)), sqrt(sum(cells)))
      }, numeric(ncol(free) + 1L)))
      matrix(block, k, dimnames = list(NULL, layout$levels[[name]]))
    }), layout$categorical)
  }
  if (!is.null(index$ordinal)) {
    # The last class's location is fixed at zero and has no error.
    out$location <- matrix(vapply(layout$ordinal, function(name) {
      c(diagonal[index$ordinal[[name]]$locations], NA_real_)
    }, numeric(k)), k, dimnames = list(NULL, layout$ordinal))
  }
  if (!is.null(index$count_means)) out$count_mean <- error_of(index$count_means, layout$count)
  if (!is.null(index$count_dispersion)) {
    positions <- if (identical(layout$dispersion, "equal"))
      rep(index$count_dispersion, each = k) else index$count_dispersion
    out$count_dispersion <- error_of(positions, layout$count)
  }
  out
}

# ---------------------------------------------------------------------------
# The EM whose M-step is other programs' weighted fits
# ---------------------------------------------------------------------------

#' One M-step, given posterior weights
#'
#' Rows with the same observed value and class are merged, their weights
#' summed: a weighted likelihood is a sum over rows, so this changes no
#' estimate and keeps each fit to a few dozen rows.
#'
#' @param posteriors Rows by classes.
#' @param layout From [ordinal_count_layout()].
#' @return Estimates as [ordinal_count_e_step()] takes them.
ordinal_count_m_step <- function(posteriors, layout) {
  data <- layout$data
  k <- layout$k
  classes <- seq_len(k)
  # One row per (value, class) with its summed posterior weight; `class` is a
  # factor whose first level is the last class, VGAM's and glm's baseline.
  weighted_table <- function(values) {
    observed <- !is.na(values)
    stacked <- data.frame(value = rep(values[observed], k),
                          class = rep(classes, each = sum(observed)),
                          weight = as.vector(posteriors[observed, , drop = FALSE]))
    merged <- stats::aggregate(weight ~ value + class, data = stacked, FUN = sum)
    merged$class <- factor(merged$class, levels = c(k, classes[-k]))
    merged
  }
  gaussian <- lapply(stats::setNames(nm = layout$gaussian), function(name) {
    y <- data[[name]]
    observed <- !is.na(y)
    w <- posteriors[observed, , drop = FALSE]
    mean <- colSums(w * y[observed]) / colSums(w)
    variance <- colSums(w * outer(y[observed], mean, "-")^2) / colSums(w)
    list(mean = mean, variance = variance)
  })
  responses <- lapply(stats::setNames(nm = layout$categorical), function(name) {
    code <- layout$codes[[name]]
    observed <- !is.na(code)
    totals <- rowsum(posteriors[observed, , drop = FALSE],
                     factor(code[observed], seq_along(layout$levels[[name]])), reorder = TRUE)
    block <- t(totals) / colSums(totals)
    matrix(block, k, dimnames = list(NULL, layout$levels[[name]]))
  })
  ordinal <- lapply(stats::setNames(nm = layout$ordinal), function(name) {
    table <- weighted_table(layout$codes[[name]])
    table$value <- factor(table$value, seq_along(layout$levels[[name]]), ordered = TRUE)
    # Convergence is judged on the coefficients: on the deviance VGAM stops
    # at a half-step and warns, two digits short of the maximum.
    fitted <- VGAM::vglm(value ~ class, family = VGAM::acat(parallel = TRUE),
                         data = table, weights = weight,
                         control = VGAM::vglm.control(epsilon = 1e-11, maxit = 200L,
                                                      criterion = "coefficients"))
    block <- VGAM::predictvglm(fitted, type = "response",
                           newdata = data.frame(class = factor(classes, levels(table$class))))
    matrix(unname(as.matrix(block)), k, dimnames = list(NULL, layout$levels[[name]]))
  })
  counts <- lapply(stats::setNames(nm = layout$count), function(name) {
    table <- weighted_table(data[[name]])
    if (identical(layout$dispersion, "equal")) {
      fitted <- .oc_glm_nb(value ~ 0 + class, table)
      coefficients <- stats::coef(fitted$fit)[sprintf("class%d", classes)]
      return(list(mean = unname(exp(coefficients)), dispersion = rep(1 / fitted$theta, k)))
    }
    by_class <- lapply(classes, function(class) {
      rows <- table[as.integer(as.character(table$class)) == class, , drop = FALSE]
      if (identical(layout$dispersion, "varying")) {
        fitted <- .oc_glm_nb(value ~ 1, rows)
        c(exp(unname(stats::coef(fitted$fit))), 1 / fitted$theta)
      } else {
        c(sum(rows$weight * rows$value) / sum(rows$weight), NA_real_)
      }
    })
    list(mean = vapply(by_class, `[[`, numeric(1), 1L),
         dispersion = vapply(by_class, `[[`, numeric(1), 2L))
  })
  named <- function(list_of, field, names) {
    if (length(names) == 0L) return(NULL)
    matrix(vapply(list_of, `[[`, numeric(k), field), k, dimnames = list(NULL, names))
  }
  list(pro = colMeans(posteriors),
       means = named(gaussian, "mean", layout$gaussian),
       variances = named(gaussian, "variance", layout$gaussian),
       responses = if (length(responses) > 0L) responses,
       ordinal = if (length(ordinal) > 0L) ordinal,
       count_means = named(counts, "mean", layout$count),
       count_dispersion = if (!identical(layout$dispersion, "none"))
         named(counts, "dispersion", layout$count))
}

#' A weighted NB2 fit by MASS, its theta taken to MASS's own maximum
#'
#' glm.nb() alternates an IRLS fit of the means with `theta.ml()`, which it
#' calls with its default step tolerance `.Machine$double.eps^0.25`; that
#' leaves theta about 1e-8 short of its maximum (measured against a
#' `theta.ml()` run to a Newton step of 1e-15 on a fixture like the suite's), more than the
#' reference EM's fixed-point tolerance. Its convergence `epsilon` cannot
#' close the gap: it also bounds the change in a log-likelihood of a few
#' hundred, which floating point rounding keeps above 1e-13 ("alternation
#' limit reached"). So glm.nb runs at epsilon 1e-10, and theta is then taken
#' to its maximum by MASS's own `theta.ml()` at glm.nb's fitted means, until
#' its Newton step is below 1e-12 (smaller steps are lost to rounding, and
#' theta.ml then reports "iteration limit reached"). The
#' means need nothing more: with one mean per class they solve the weighted
#' score equation for any theta.
#'
#' @param formula `value ~ 1` or `value ~ 0 + class`.
#' @param table Rows with `value`, `class` and `weight`.
#' @return A list with the glm.nb `fit` and the polished `theta`.
.oc_glm_nb <- function(formula, table) {
  fitted <- MASS::glm.nb(formula, data = table, weights = weight,
                         control = stats::glm.control(epsilon = 1e-10, maxit = 200L))
  theta <- MASS::theta.ml(table$value, stats::fitted(fitted), sum(table$weight),
                          table$weight, limit = 200L, eps = 1e-12)
  list(fit = fitted, theta = as.vector(theta))
}

#' Put a reference's classes in another order
#' @param estimates As for [ordinal_count_e_step()].
#' @param posteriors Rows by classes.
#' @param order A permutation from `align_classes()`.
#' @return A list with the permuted `estimates` and `posteriors`.
ordinal_count_permute <- function(estimates, posteriors, order) {
  rows <- function(block) if (is.null(block)) NULL else block[order, , drop = FALSE]
  estimates$pro <- estimates$pro[order]
  estimates[c("means", "variances", "count_means", "count_dispersion")] <-
    lapply(estimates[c("means", "variances", "count_means", "count_dispersion")], rows)
  # An absent block stays NULL, not an empty list.
  if (!is.null(estimates$responses)) estimates$responses <- lapply(estimates$responses, rows)
  if (!is.null(estimates$ordinal)) estimates$ordinal <- lapply(estimates$ordinal, rows)
  list(estimates = estimates, posteriors = posteriors[, order, drop = FALSE])
}

#' Run the reference EM from given posteriors to its fixed point
#'
#' Stops when no free natural parameter moves by more than `tolerance` in an
#' EM step; the bound on steps raises a classed warning, and the EM must not
#' lose likelihood from its first M-step on.
#'
#' @param posteriors Starting posteriors, rows by classes.
#' @param layout From [ordinal_count_layout()].
#' @param tolerance Largest parameter change at the fixed point.
#' @param max_steps Bound on the EM steps.
#' @return A list with `estimates`, `posteriors`, `log_likelihood` and `steps`.
ordinal_count_em <- function(posteriors, layout, tolerance = 1e-10, max_steps = 3000L) {
  first <- ordinal_count_m_step(posteriors, layout)
  start <- list(estimates = first, e = ordinal_count_e_step(first, layout),
                settled = FALSE, steps = 1L)
  # Reduce() rather than a loop with a break: once the parameters have stopped
  # moving each remaining element returns the state untouched.
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    estimates <- ordinal_count_m_step(state$e$posteriors, layout)
    e <- ordinal_count_e_step(estimates, layout)
    stopifnot("the reference EM must not lose likelihood" =
                e$log_likelihood >= state$e$log_likelihood - 1e-9)
    change <- max(abs(unlist(ordinal_count_free(estimates, layout)) -
                        unlist(ordinal_count_free(state$estimates, layout))))
    list(estimates = estimates, e = e, settled = change < tolerance, steps = index + 1L)
  }, seq_len(max_steps), init = start)
  if (!state$settled) {
    warning(warningCondition(sprintf(
      "the reference EM was still moving after %d steps", max_steps),
      class = "multilpa_harness_no_converge"))
  }
  list(estimates = state$estimates, posteriors = state$e$posteriors,
       log_likelihood = state$e$log_likelihood, steps = state$steps)
}

# ---------------------------------------------------------------------------
# depmixS4
# ---------------------------------------------------------------------------

#' Run a depmixS4 call, muffling only its complete-case `nobs()` warning
#'
#' With a missing response depmixS4's `nobs()` warns that it counts complete
#' cases. `nobs()` only feeds its own AIC and BIC; its likelihood and EM
#' integrate the missing value out (checked: its maximum is the base-R
#' likelihood's). Every other warning passes through to the harness.
#'
#' @param expression The call.
#' @return Its value.
.oc_depmix_quiet <- function(expression) {
  withCallingHandlers(expression, warning = function(condition) {
    if (identical(conditionMessage(condition),
                  "missing values detected; nobs is number of cases with complete data")) {
      invokeRestart("muffleWarning")
    }
  })
}

#' Continue depmixS4's EM from its own solution to the fixed point
#'
#' `fit()` with `maxit = 0` and no random start performs one EM step from the
#' model's current parameters. Steps are taken until no parameter moves by
#' more than 1e-12 (as the depmixs4 suite refines).
#'
#' @param theirs A fitted depmixS4 model (`mix` or `depmix`).
#' @param max_steps Bound on the EM steps.
#' @return The refined model.
ordinal_count_depmix_refine <- function(theirs, max_steps = 5000L) {
  step <- function(object) {
    .oc_depmix_quiet(depmixS4::fit(object, verbose = FALSE, emcontrol = depmixS4::em.control(
      maxit = 0, tol = 0, crit = "absolute", random.start = FALSE)))
  }
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    moved <- step(state$model)
    change <- max(abs(depmixS4::getpars(moved) - depmixS4::getpars(state$model)))
    list(model = moved, settled = change < 1e-12)
  }, seq_len(max_steps), init = list(model = theirs, settled = FALSE))
  if (!state$settled) {
    warning(warningCondition(sprintf(
      "depmixS4's EM was still moving after %d refining steps", max_steps),
      class = "multilpa_harness_no_converge"))
  }
  stopifnot("refining depmixS4's solution must not lose likelihood" =
              as.numeric(.oc_depmix_quiet(depmixS4::logLik(state$model))) >=
                as.numeric(.oc_depmix_quiet(depmixS4::logLik(theirs))) - 1e-9)
  state$model
}

#' depmixS4's response parameters as estimates
#'
#' depmixS4 lays out one block per state holding, for each response in
#' formula order, a Poisson log-mean intercept, or a Gaussian intercept and
#' SD.
#'
#' @param blocks States by response parameters, in depmixS4's order.
#' @param layout From [ordinal_count_layout()]; responses are the counts then
#'   the Gaussian indicators, as the suite's formulas list them.
#' @return Estimates without `pro`.
ordinal_count_depmix_responses <- function(blocks, layout) {
  k <- layout$k
  n_count <- length(layout$count)
  gaussian <- blocks[, n_count + seq_len(2L * length(layout$gaussian)), drop = FALSE]
  list(count_means = matrix(exp(blocks[, seq_len(n_count)]), k,
                            dimnames = list(NULL, layout$count)),
       means = matrix(gaussian[, c(TRUE, FALSE)], k, dimnames = list(NULL, layout$gaussian)),
       variances = matrix(gaussian[, c(FALSE, TRUE)]^2, k,
                          dimnames = list(NULL, layout$gaussian)))
}

# ---------------------------------------------------------------------------
# Latent transitions
# ---------------------------------------------------------------------------

#' The forward log-likelihood of a first-order latent Markov model
#'
#' Rows are in sequence and occasion order. The recursion is scaled at every
#' occasion and runs over occasions, the short dimension, with all sequences
#' as one matrix product.
#'
#' @param log_emission Rows by states, log f(y_it | state).
#' @param initial Initial state probabilities.
#' @param transition States by states, rows `from`.
#' @param ntimes Sequence lengths, in data order.
#' @return Each sequence's log-likelihood.
ordinal_count_forward <- function(log_emission, initial, transition, ntimes) {
  last <- cumsum(ntimes)
  slot <- matrix(t(vapply(seq_along(ntimes), function(i) {
    c(last[[i]] - ntimes[[i]] + seq_len(ntimes[[i]]),
      rep(NA_integer_, max(ntimes) - ntimes[[i]]))
  }, integer(max(ntimes)))), length(ntimes))
  emit <- function(rows) {
    block <- log_emission[rows, , drop = FALSE]
    peak <- apply(block, 1L, max)
    list(density = exp(block - peak), peak = peak)
  }
  start <- emit(slot[, 1L])
  alpha <- sweep(start$density, 2L, initial, "*")
  scale <- rowSums(alpha)
  state <- list(alpha = alpha / scale, total = start$peak + log(scale))
  state <- Reduce(function(state, occasion) {
    rows <- slot[, occasion]
    active <- !is.na(rows)
    step <- emit(rows[active])
    alpha <- (state$alpha[active, , drop = FALSE] %*% transition) * step$density
    scale <- rowSums(alpha)
    state$alpha[active, ] <- alpha / scale
    state$total[active] <- state$total[active] + step$peak + log(scale)
    state
  }, seq_len(ncol(slot))[-1L], init = state)
  state$total
}

#' Initial and transition probabilities from multinomial logits
#'
#' The initial logits are against the last state; each transition row's
#' against staying in its own state, as the package's coefficient tables
#' report them.
#'
#' @param initial_logits `K - 1` logits.
#' @param transition_logits States by `K - 1` logits, the `to` states other
#'   than the row's own in order.
#' @return A list with `initial` and `transition`.
.oc_transition_probabilities <- function(initial_logits, transition_logits) {
  k <- length(initial_logits) + 1L
  softmax <- function(v) exp(v - max(v)) / sum(exp(v - max(v)))
  transition <- t(vapply(seq_len(k), function(from) {
    logits <- numeric(k)
    logits[-from] <- transition_logits[from, ]
    softmax(logits)
  }, numeric(k)))
  list(initial = softmax(c(initial_logits, 0)), transition = transition)
}

#' Standard errors of the initial and transition logits, from scratch
#'
#' The forward likelihood over every sequence, in the logits of
#' [.oc_transition_probabilities()] and the free natural emission parameters
#' of [ordinal_count_free()]; its Hessian is taken by Richardson
#' extrapolation at the reference's estimates, after the recursion has been
#' shown to reproduce the reference's maximum.
#'
#' @param initial,transition The reference's probabilities.
#' @param emission Its emission estimates (no `pro`).
#' @param layout From [ordinal_count_layout()].
#' @param ntimes Sequence lengths.
#' @param log_likelihood The reference's maximum.
#' @return A list with `initial` (`K - 1` errors) and `transition` (states by
#'   `K - 1` errors, laid out as the logits).
ordinal_count_transition_errors <- function(initial, transition, emission, layout, ntimes,
                                            log_likelihood) {
  k <- length(initial)
  initial_logits <- log(initial[-k]) - log(initial[k])
  transition_logits <- t(vapply(seq_len(k), function(from) {
    log(transition[from, -from]) - log(transition[from, from])
  }, numeric(k - 1L)))
  transition_logits <- matrix(transition_logits, k)
  skeleton <- list(initial = initial_logits, transition = transition_logits,
                   emission = ordinal_count_free(emission, layout))
  at <- function(values) {
    parameters <- utils::relist(values, skeleton)
    chain <- .oc_transition_probabilities(parameters$initial,
                                          matrix(parameters$transition, k))
    log_emission <- ordinal_count_log_density(
      ordinal_count_estimates(parameters$emission, layout), layout)
    sum(ordinal_count_forward(log_emission, chain$initial, chain$transition, ntimes))
  }
  theta <- unlist(skeleton)
  stopifnot("the base-R forward recursion must reproduce the reference's maximum" =
              abs(at(theta) - log_likelihood) < 1e-8 * abs(log_likelihood))
  errors <- sqrt(diag(solve(-numDeriv::hessian(at, theta))))
  index <- utils::relist(seq_along(theta), skeleton)
  list(initial = errors[index$initial],
       transition = matrix(errors[index$transition], k))
}

#' A Wald coefficient row, as the package's coefficient tables report it
#' @param estimate,standard_error Logit and its error.
#' @param level Confidence level.
#' @return A `data.frame` of the derived columns.
.oc_wald <- function(estimate, standard_error, level = 0.95) {
  statistic <- estimate / standard_error
  half <- stats::qnorm(1 - (1 - level) / 2) * standard_error
  data.frame(estimate = estimate, standard_error = standard_error, statistic = statistic,
             p_value = 2 * stats::pnorm(-abs(statistic)), conf_low = estimate - half,
             conf_high = estimate + half, odds_ratio = exp(estimate), estimated = TRUE)
}

#' Every table of a general transition fit, from a depmixS4 solution
#'
#' @param reference A list with `initial`, `transition`, `emission`,
#'   `posteriors`, `log_likelihood`, `n_parameters` and `errors` (from
#'   [ordinal_count_transition_errors()]), states in the fit's order.
#' @param layout From [ordinal_count_layout()].
#' @param frame The data, in sequence and occasion order, with `t`.
#' @param ntimes Sequence lengths.
#' @param min_variance The variance floor, for the Gaussian diagnostics'
#'   boundary flag.
#' @return A named list of expected tables, shaped as `fit_tables()` reports a
#'   `multilpa_lta` fit.
ordinal_count_lta_tables <- function(reference, layout, frame, ntimes, min_variance) {
  k <- layout$k
  label <- function(i) sprintf("profile_%d", i)
  posteriors <- reference$posteriors
  inner <- setdiff(seq_len(nrow(frame)), cumsum(ntimes))
  # A transition row is estimable when its origin is occupied, in expectation,
  # at an occasion that has a successor.
  occupied <- colSums(posteriors[inner, , drop = FALSE]) > 0
  pairs <- expand.grid(to = seq_len(k), from = seq_len(k))
  moves <- pairs[pairs$from != pairs$to, , drop = FALSE]
  # The logits of row `from` are laid out over the other states in order.
  move_column <- vapply(seq_len(nrow(moves)), function(i) {
    match(moves$to[[i]], seq_len(k)[-moves$from[[i]]])
  }, integer(1))
  logit <- function(i) {
    log(reference$transition[moves$from[[i]], moves$to[[i]]]) -
      log(reference$transition[moves$from[[i]], moves$from[[i]]])
  }
  transition_estimates <- vapply(seq_len(nrow(moves)), logit, numeric(1))
  transition_errors <- reference$errors$transition[cbind(moves$from, move_column)]
  initial_estimates <- log(reference$initial[-k]) - log(reference$initial[k])
  log_likelihood <- reference$log_likelihood
  q <- reference$n_parameters
  emission <- reference$emission
  n_count <- length(layout$count)
  shared <- harness_reference(log_likelihood = log_likelihood, n_parameters = q,
                              posteriors = posteriors, means = emission$means,
                              variances = emission$variances, data = frame[layout$gaussian])
  list(
    transitions = data.frame(group_class = "group_class_1", occasion = NA_integer_,
                             from = label(pairs$from), to = label(pairs$to),
                             probability = reference$transition[cbind(pairs$from, pairs$to)],
                             estimated = occupied[pairs$from], stringsAsFactors = FALSE),
    transition_coefficients = data.frame(
      group_class = "group_class_1", from = label(moves$from), to = label(moves$to),
      term = "(Intercept)", .oc_wald(transition_estimates, transition_errors),
      stringsAsFactors = FALSE),
    initial = data.frame(group_class = "group_class_1", profile = label(seq_len(k)),
                         probability = reference$initial, stringsAsFactors = FALSE),
    initial_coefficients = data.frame(
      group_class = "group_class_1", profile = label(seq_len(k - 1L)), term = "(Intercept)",
      .oc_wald(initial_estimates, reference$errors$initial), stringsAsFactors = FALSE),
    profiles = data.frame(occasion = NA_integer_,
                          profile = label(rep(seq_len(k), length(layout$gaussian))),
                          indicator = rep(layout$gaussian, each = k),
                          mean = as.vector(emission$means), variance = as.vector(emission$variances),
                          stringsAsFactors = FALSE),
    count_means = data.frame(occasion = NA_integer_, profile = rep(seq_len(k), n_count),
                             indicator = rep(layout$count, each = k),
                             mean = as.vector(emission$count_means), stringsAsFactors = FALSE),
    group_classes = data.frame(group_class = "group_class_1", weight = 1,
                               count = length(ntimes), stringsAsFactors = FALSE),
    assignments = keyed_table(data.frame(
      .row = seq_len(nrow(frame)), t = frame$t,
      profile = max.col(posteriors, ties.method = "first"), group_class = 1L,
      uncertainty = 1 - apply(posteriors, 1L, max),
      stats::setNames(as.data.frame(posteriors), sprintf("posterior_profile_%d", seq_len(k)))),
      ".row"),
    fit = data.frame(n_profiles = k, n_group_classes = 1L, log_likelihood = log_likelihood,
                     n_parameters = q, aic = -2 * log_likelihood + 2 * q,
                     bic = -2 * log_likelihood + q * log(length(ntimes)),
                     bic_individual = -2 * log_likelihood + q * log(nrow(frame)),
                     n_groups = length(ntimes), n_observations = nrow(frame)),
    gaussian_diagnostics = .reference_gaussian_diagnostics(shared, layout$gaussian,
                                                           min_variance))
}
