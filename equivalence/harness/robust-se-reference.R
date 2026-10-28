# From-scratch reference for the observed, robust (CR1 sandwich) and outer-
# product (BHHH) covariances of Gaussian mixtures, single and two level, with
# and without membership covariates, and for the Wald columns that
# `parameter_inference()` derives from them. Used by suites/robust-se.R.
#
# Nothing here calls the package. The likelihood is written in base R from the
# model's definition (Vermunt 2003 for two levels):
#
#   log L_g = log sum_d delta_d(z_g) prod_{i in g} sum_k gamma_k(x_i | d) f_k(y_i)
#
# with f_k the multivariate normal density (Cholesky log density), gamma the
# profile probabilities within group class d (a multinomial logit in x_i, last
# profile the baseline, for a covariate model) and delta the group-class
# probabilities (a multinomial logit in z_g, last group class the baseline).
# One group class and one row per group is the single-level mixture.
#
# Every model is parameterised by its free natural parameters: means,
# variances or covariance elements, the first K - 1 profile probabilities in
# each group class and the first D - 1 group-class probabilities, or the
# membership logits of a covariate model. The reported natural vector is a
# linear function of these (the last probability is one minus the others).
#
# Covariances, with s_g the per-group score (numDeriv Jacobian of the group
# log-likelihoods, Richardson extrapolation, with the step rule of
# `.rse_derivative()`) and A minus the numDeriv Hessian of the total:
#
#   observed  A^-1
#   robust    A^-1 B A^-1 x G / (G - 1),  B = sum_g s_g s_g'  (CR1, the
#             default of sandwich::vcovCL(), as ?parameter_inference states)
#   opg       B^-1
#   scaling   tr(A^-1 B) / q  (MLR scaling correction, B without CR1)
#
# At a maximum (zero total score) all three transform with the Jacobian of any
# reparameterisation -- the Hessian's second-derivative term is multiplied by
# the zero score -- so the natural chart used here gives the same standard
# errors as the package's log-variance / logit chart. That is why the
# reference point must be the reference likelihood's own maximum:
# `rse_maximise()` takes at least one Newton step, and continues until the
# score is below 1e-6 standard errors.

#' Describe a Gaussian mixture's free and reported parameters
#'
#' @param indicators Names of the continuous indicators.
#' @param n_profiles,n_group_classes Numbers of profiles and group classes.
#' @param covariance `"diagonal"` (one variance per profile and indicator) or
#'   `"full"` (an unrestricted covariance matrix per profile).
#' @param profile_covariates,group_covariates Covariate names of a membership
#'   regression; empty for a model with free class probabilities.
#' @return A list of class `rse_model`: the dimensions and `labels`, a
#'   `data.frame` with one row per reported natural parameter: `level`,
#'   `outcome`, `term` and `parameter` as `?parameter_inference` documents
#'   them, `kind` (`"probability"`, `"positive"` or `"real"`, the scale a Wald
#'   interval is formed on) and `free` (whether it is a coordinate of the
#'   parameter vector `theta`).
rse_model <- function(indicators, n_profiles, n_group_classes = 1L,
                      covariance = c("diagonal", "full"),
                      profile_covariates = character(), group_covariates = character()) {
  covariance <- match.arg(covariance)
  stopifnot(
    "`indicators` must be distinct names" = is.character(indicators) &&
      length(indicators) > 0L && !anyDuplicated(indicators),
    "`n_profiles` and `n_group_classes` must be positive whole numbers" =
      n_profiles >= 1L && n_group_classes >= 1L,
    "group covariates need several group classes" =
      length(group_covariates) == 0L || n_group_classes > 1L)
  k <- as.integer(n_profiles)
  d <- as.integer(n_group_classes)
  j <- length(indicators)
  row <- function(level, outcome, term, parameter, kind, free = TRUE) {
    data.frame(level = level, outcome = outcome, term = term, parameter = parameter,
               kind = kind, free = free, stringsAsFactors = FALSE)
  }
  profiles <- sprintf("profile_%d", seq_len(k))
  group_classes <- sprintf("group_class_%d", seq_len(d))
  means <- row("measurement", rep(profiles, each = j), rep(indicators, times = k),
               "mean", "real")
  cells <- if (covariance == "diagonal") cbind(seq_len(j), seq_len(j)) else
    which(lower.tri(diag(j), diag = TRUE), arr.ind = TRUE)
  on_diagonal <- cells[, 1L] == cells[, 2L]
  spread <- row("measurement", rep(profiles, each = nrow(cells)),
                rep(if (covariance == "diagonal") indicators else
                  paste(indicators[cells[, 1L]], indicators[cells[, 2L]], sep = ":"),
                  times = k),
                if (covariance == "diagonal") "variance" else "covariance",
                rep(ifelse(on_diagonal, "positive", "real"), times = k))
  covariates <- length(profile_covariates) > 0L || length(group_covariates) > 0L
  membership <- if (!covariates) {
    rbind(row("profile", rep(profiles, times = d), rep(group_classes, each = k),
              "probability", "probability", free = rep(seq_len(k) < k, times = d)),
          row("group", group_classes, NA_character_, "probability", "probability",
              free = seq_len(d) < d))
  } else {
    intercepts <- if (d == 1L) "(Intercept)" else group_classes
    profile_terms <- c(intercepts, profile_covariates)
    group_terms <- c("(Intercept)", group_covariates)
    rbind(row("profile", rep(profiles[-k], each = length(profile_terms)),
              rep(profile_terms, times = k - 1L), "coefficient", "real"),
          if (d > 1L) row("group", rep(group_classes[-d], each = length(group_terms)),
                          rep(group_terms, times = d - 1L), "coefficient", "real"))
  }
  labels <- rbind(means, spread, membership)
  row.names(labels) <- NULL
  structure(list(indicators = indicators, n_profiles = k, n_group_classes = d,
                 covariance = covariance, cells = cells,
                 profile_covariates = profile_covariates,
                 group_covariates = group_covariates, covariates = covariates,
                 labels = labels),
            class = "rse_model")
}

#' The reported natural parameters from the free coordinates
#'
#' A dropped profile probability is one minus the others of its group class,
#' and a dropped group-class probability one minus the other group classes;
#' everything else is a coordinate.
#'
#' @param model An `rse_model`.
#' @param theta The free coordinates, in the order of `model$labels[free, ]`.
#' @return The natural vector, in the order of `model$labels`.
rse_natural <- function(model, theta) {
  labels <- model$labels
  natural <- numeric(nrow(labels))
  natural[labels$free] <- theta
  dependent <- which(!labels$free)
  natural[dependent] <- vapply(dependent, function(position) {
    siblings <- labels$level == labels$level[[position]] &
      labels$parameter == "probability" & labels$free &
      (if (identical(labels$level[[position]], "profile"))
        labels$term %in% labels$term[[position]] else TRUE)
    1 - sum(natural[siblings])
  }, numeric(1))
  natural
}

#' The free coordinates of another program's estimates
#'
#' @param model An `rse_model` without membership covariates.
#' @param means Profiles-by-indicators matrix with indicator names.
#' @param covariances Indicators-by-indicators-by-profiles array.
#' @param gamma Group-classes-by-profiles matrix of profile probabilities.
#' @param delta Group-class probabilities.
#' @return `theta`, in the order of `model$labels[free, ]`.
rse_theta_from_parameters <- function(model, means, covariances, gamma, delta = 1) {
  stopifnot("membership covariates have logits, not probabilities" = !model$covariates)
  labels <- model$labels
  class_number <- function(label, prefix) as.integer(sub(prefix, "", label, fixed = TRUE))
  natural <- vapply(seq_len(nrow(labels)), function(i) {
    profile <- if (identical(labels$level[[i]], "group")) NA_integer_ else
      class_number(labels$outcome[[i]], "profile_")
    switch(labels$parameter[[i]],
           mean = means[profile, labels$term[[i]]],
           variance = covariances[labels$term[[i]], labels$term[[i]], profile],
           covariance = {
             pair <- strsplit(labels$term[[i]], ":", fixed = TRUE)[[1L]]
             covariances[pair[[1L]], pair[[2L]], profile]
           },
           probability = if (is.na(profile)) {
             delta[[class_number(labels$outcome[[i]], "group_class_")]]
           } else gamma[class_number(labels$term[[i]], "group_class_"), profile])
  }, numeric(1))
  natural[labels$free]
}

#' The free coordinates of a reported natural vector
#' @param model An `rse_model`.
#' @param table A `data.frame` with `level`, `outcome`, `term`, `parameter`
#'   and `estimate`, such as an inference table.
#' @return `theta`, in the order of `model$labels[free, ]`.
rse_theta <- function(model, table) {
  free <- model$labels[model$labels$free, , drop = FALSE]
  at <- match(.rse_key(free), .rse_key(table))
  stopifnot("every free parameter must be in the table" = !anyNA(at))
  table$estimate[at]
}

#' One text key per parameter, for joining on the four label columns
#' @param table A `data.frame` with the four label columns.
#' @return A character vector.
.rse_key <- function(table) {
  paste(table$level, table$outcome, ifelse(is.na(table$term), "<none>", table$term),
        table$parameter, sep = "\r")
}

#' The model's parameters as arrays
#' @param model An `rse_model`.
#' @param natural The natural vector from [rse_natural()].
#' @return A list with `means` (profiles x indicators), `covariances`
#'   (indicators x indicators x profiles) and either `gamma` (group classes x
#'   profiles) and `delta`, or the logit blocks `profile_logits` (terms x
#'   (K - 1)) and `group_logits` (terms x (D - 1)).
.rse_unpack <- function(model, natural) {
  labels <- model$labels
  k <- model$n_profiles
  d <- model$n_group_classes
  j <- length(model$indicators)
  means <- matrix(natural[labels$parameter == "mean"], k, j, byrow = TRUE,
                  dimnames = list(NULL, model$indicators))
  spread <- matrix(natural[labels$parameter %in% c("variance", "covariance")], k,
                   nrow(model$cells), byrow = TRUE)
  covariances <- vapply(seq_len(k), function(profile) {
    m <- matrix(0, j, j)
    m[model$cells] <- spread[profile, ]
    m[model$cells[, 2:1, drop = FALSE]] <- spread[profile, ]
    m
  }, matrix(0, j, j))
  covariances <- array(covariances, c(j, j, k),
                       dimnames = list(model$indicators, model$indicators, NULL))
  out <- list(means = means, covariances = covariances)
  if (!model$covariates) {
    # Rows are group classes, each holding its K profile probabilities.
    out$gamma <- matrix(natural[labels$level == "profile"], d, k, byrow = TRUE)
    out$delta <- natural[labels$level == "group"]
    return(out)
  }
  profile_rows <- labels$level == "profile"
  out$profile_logits <- matrix(natural[profile_rows], ncol = k - 1L,
                               dimnames = list(unique(labels$term[profile_rows]), NULL))
  group_rows <- labels$level == "group"
  out$group_logits <- if (d > 1L) matrix(natural[group_rows], ncol = d - 1L,
                                         dimnames = list(unique(labels$term[group_rows]), NULL))
  out
}

#' Row-wise log-sum-exp
#' @param m A numeric matrix.
#' @return One value per row.
.rse_log_sum_exp <- function(m) {
  peak <- apply(m, 1L, max)
  peak + log(rowSums(exp(m - peak)))
}

#' Everything the likelihood is built from, at one parameter value
#'
#' @param model An `rse_model`.
#' @param theta The free coordinates.
#' @param data A list from [rse_data()].
#' @return A list with `group_log_likelihood` (one value per group, named by
#'   group), `group_posteriors` (groups x group classes, rows named by group),
#'   `posteriors` (rows x profiles, the marginal posterior
#'   `sum_d P(d | y_g) P(k | d, y_i)`), `row_conditional_entropy`
#'   (`sum_d P(d | y_g) H(k | d, y_i)`) and the unpacked `parameters`.
rse_evaluate <- function(model, theta, data) {
  parameters <- .rse_unpack(model, rse_natural(model, theta))
  k <- model$n_profiles
  d <- model$n_group_classes
  y <- data$y
  # log f_k(y_i), from the Cholesky factor of each profile's covariance.
  log_density <- vapply(seq_len(k), function(profile) {
    root <- chol(matrix(parameters$covariances[, , profile], ncol(y)))
    standardised <- backsolve(root, t(y) - parameters$means[profile, ], transpose = TRUE)
    -0.5 * colSums(standardised^2) - sum(log(diag(root))) - 0.5 * ncol(y) * log(2 * base::pi)
  }, numeric(nrow(y)))
  log_density <- matrix(log_density, nrow(y), k)
  # log gamma_k(x_i | d), one rows x profiles matrix per group class.
  log_profile <- lapply(seq_len(d), function(group_class) {
    if (!model$covariates) {
      return(matrix(log(parameters$gamma[group_class, ]), nrow(y), k, byrow = TRUE))
    }
    logits <- parameters$profile_logits
    intercept <- if (d == 1L) logits["(Intercept)", ] else
      logits[sprintf("group_class_%d", group_class), ]
    slopes <- logits[model$profile_covariates, , drop = FALSE]
    eta <- cbind(matrix(intercept, nrow(y), k - 1L, byrow = TRUE) + data$x %*% slopes, 0)
    eta - .rse_log_sum_exp(eta)
  })
  # log delta_d(z_g), groups x group classes.
  n_groups <- length(data$group_labels)
  log_group <- if (!model$covariates) {
    matrix(log(parameters$delta), n_groups, d, byrow = TRUE)
  } else if (d == 1L) {
    matrix(0, n_groups, 1L)
  } else {
    logits <- parameters$group_logits
    eta <- cbind(matrix(logits["(Intercept)", ], n_groups, d - 1L, byrow = TRUE) +
                   data$z %*% logits[model$group_covariates, , drop = FALSE], 0)
    eta - .rse_log_sum_exp(eta)
  }
  joint <- lapply(log_profile, function(lp) lp + log_density)
  row_given_class <- vapply(joint, .rse_log_sum_exp, numeric(nrow(y)))
  row_given_class <- matrix(row_given_class, nrow(y), d)
  by_group <- log_group + rowsum(row_given_class, data$group, reorder = TRUE)
  group_log_likelihood <- .rse_log_sum_exp(by_group)
  group_posteriors <- exp(by_group - group_log_likelihood)
  # P(k | d, y_i) for each group class, weighted by P(d | y_g).
  within <- lapply(seq_len(d), function(group_class) {
    exp(joint[[group_class]] - row_given_class[, group_class])
  })
  weight <- group_posteriors[data$group, , drop = FALSE]
  posteriors <- Reduce(`+`, Map(function(p, w) p * w, within, split(weight, col(weight))))
  conditional_entropy <- Reduce(`+`, Map(function(p, w) {
    # 0 log 0 = 0 for a profile a row cannot belong to.
    terms <- p * log(p)
    terms[p == 0] <- 0
    -w * rowSums(terms)
  }, within, split(weight, col(weight))))
  dimnames(group_posteriors) <- list(data$group_labels, NULL)
  list(group_log_likelihood = stats::setNames(group_log_likelihood, data$group_labels),
       group_posteriors = group_posteriors, posteriors = posteriors,
       row_conditional_entropy = conditional_entropy, parameters = parameters)
}

#' A data set in the shape the reference likelihood reads
#' @param frame The data frame the fit was made from.
#' @param indicators Indicator columns.
#' @param id Group identifier column; each row is its own group when it is
#'   unique.
#' @param profile_covariates,group_covariates Covariate columns (a group
#'   covariate is read from each group's first row).
#' @return A list with `y`, `group` (each row's group number, groups numbered
#'   in order of first appearance), `group_labels`, `x` and `z`.
rse_data <- function(frame, indicators, id, profile_covariates = character(),
                     group_covariates = character()) {
  labels <- unique(as.character(frame[[id]]))
  group <- match(as.character(frame[[id]]), labels)
  first <- match(seq_along(labels), group)
  list(y = as.matrix(frame[indicators]), group = group, group_labels = labels,
       x = as.matrix(frame[profile_covariates]),
       z = as.matrix(frame[first, group_covariates, drop = FALSE]))
}

#' The group log-likelihoods as a function of the free coordinates
#' @inheritParams rse_evaluate
#' @return A function of `theta` returning one value per group, with a
#'   `bounded` attribute flagging the coordinates that are variances or
#'   probabilities (see [.rse_derivative()]).
rse_group_log_likelihood <- function(model, data) {
  structure(function(theta) unname(rse_evaluate(model, theta, data)$group_log_likelihood),
            bounded = model$labels$kind[model$labels$free] != "real")
}

#' numDeriv derivatives with a step that does not collapse near zero
#'
#' numDeriv's step is relative, `d |theta_j|`, so a mean near zero (a mixture
#' of profiles at -1.5 and 1.5) is differenced over steps of 1e-6 and less,
#' where rounding in a log-likelihood of order 1e3 dominates. Each unbounded
#' coordinate (a mean, a covariance off the diagonal, a logit) is therefore
#' differenced in the shifted chart `u_j = 1 + (theta_j - theta0_j) /
#' max(|theta0_j|, 1)`, evaluated at `u = 1`: its step is relative for a large
#' coordinate and absolute for a small one. A variance or a probability keeps
#' numDeriv's relative step (scale `|theta0_j|`), which never leaves its
#' range: the Hessian's first step is a tenth of the coordinate. The
#' derivatives are carried back by the diagonal chain rule.
#'
#' @param f A function of `theta`, scalar (`gradient`, `hessian`) or vector
#'   (`jacobian`) valued.
#' @param theta The point.
#' @param what Which derivative.
#' @param bounded Which coordinates are variances or probabilities.
#' @return The gradient vector, the Hessian matrix or the Jacobian matrix
#'   (rows the values of `f`, columns the coordinates).
.rse_derivative <- function(f, theta, what = c("gradient", "hessian", "jacobian"),
                            bounded = rep(FALSE, length(theta))) {
  what <- match.arg(what)
  stopifnot("`bounded` must flag every coordinate" = length(bounded) == length(theta),
            "a bounded coordinate must be positive" = all(theta[bounded] > 0))
  scale <- ifelse(bounded, abs(theta), pmax(abs(theta), 1))
  shifted <- function(u) f(theta + (u - 1) * scale)
  at <- rep(1, length(theta))
  switch(what,
         gradient = numDeriv::grad(shifted, at) / scale,
         hessian = numDeriv::hessian(shifted, at) / tcrossprod(scale),
         jacobian = sweep(numDeriv::jacobian(shifted, at), 2L, scale, `/`))
}

#' Carry a point to the reference likelihood's maximum by Newton steps
#'
#' Each step solves the numDeriv Hessian against the numDeriv gradient. The
#' stopping rule is the score measured in standard errors, `max |g_j| *
#' sqrt([A^-1]_jj)`, the displacement of the estimate from the stationary
#' point in its own standard errors. The Richardson gradient of a total
#' log-likelihood of order 1e3 carries rounding noise of about 1e-7, which
#' is 3e-8 to 7e-8 standard errors on the suite's fixtures (measured by
#' iterating past convergence), so the rule cannot be set below that floor.
#' 1e-6 standard errors sits above it, and moves an estimate by less than
#' 1e-6 of its standard error, far inside the 1e-6 tolerance on estimates and
#' with no visible effect on a standard error. At least one step is taken.
#'
#' @param group_log_likelihood A function from [rse_group_log_likelihood()].
#' @param theta The starting coordinates.
#' @param max_steps Bound on the number of Newton steps.
#' @return A list with `theta`, `steps` and `displacement`; an error of class
#'   `rse_no_maximum` when the bound is reached or the Hessian is not
#'   negative definite.
rse_maximise <- function(group_log_likelihood, theta, max_steps = 25L) {
  total <- function(t) sum(group_log_likelihood(t))
  bounded <- attr(group_log_likelihood, "bounded")
  step_from <- function(theta, steps) {
    gradient <- .rse_derivative(total, theta, "gradient", bounded)
    information <- -.rse_derivative(total, theta, "hessian", bounded)
    root <- tryCatch(chol(information), error = function(condition) {
      stop(errorCondition(sprintf("the reference Hessian is not negative definite: %s",
                                  conditionMessage(condition)),
                          class = "rse_no_maximum", call = NULL))
    })
    covariance <- chol2inv(root)
    displacement <- max(abs(gradient) * sqrt(diag(covariance)))
    # At least one step is always taken, so a start that already satisfies
    # the rule (the package's own estimates) is still replaced by the point
    # the reference's own Newton step gives.
    if (steps >= 1L && displacement < 1e-6) {
      return(list(theta = theta, steps = steps, displacement = displacement))
    }
    if (steps >= max_steps) {
      stop(errorCondition(sprintf(
        "Newton refinement did not reach the maximum in %d steps (score %.2e SE)",
        max_steps, displacement), class = "rse_no_maximum", call = NULL))
    }
    step_from(theta + as.vector(covariance %*% gradient), steps + 1L)
  }
  step_from(theta, 0L)
}

#' The three covariances of the free coordinates
#'
#' @param group_log_likelihood A function from [rse_group_log_likelihood()].
#' @param theta The reference's maximum.
#' @return A list with `observed`, `robust`, `opg` (q x q) and `scaling`.
rse_covariances <- function(group_log_likelihood, theta) {
  bounded <- attr(group_log_likelihood, "bounded")
  scores <- .rse_derivative(group_log_likelihood, theta, "jacobian", bounded)
  information <- -.rse_derivative(function(t) sum(group_log_likelihood(t)), theta, "hessian",
                                  bounded)
  bread <- solve(information)
  meat <- crossprod(scores)
  n_groups <- nrow(scores)
  list(observed = bread,
       robust = bread %*% meat %*% bread * n_groups / (n_groups - 1),
       opg = solve(meat),
       scaling = sum(diag(bread %*% meat)) / length(theta))
}

#' Natural-scale standard errors by the delta method
#'
#' The natural vector is linear in `theta`, so its Jacobian is exact.
#'
#' @param model An `rse_model`.
#' @param theta The free coordinates.
#' @param covariance A covariance of `theta`.
#' @return One standard error per row of `model$labels`; a probability that is
#'   fixed by the model (the group class of a one-group-class model) gets 0.
rse_standard_errors <- function(model, theta, covariance) {
  jacobian <- numDeriv::jacobian(function(t) rse_natural(model, t), theta)
  sqrt(pmax(diag(jacobian %*% covariance %*% t(jacobian)), 0))
}

#' Wald columns from the documented definitions
#'
#' `?parameter_inference`: statistic `estimate / standard_error`, two-sided
#' normal p-value, both `NA` for a parameter whose null lies on the boundary
#' (a variance, a probability, a covariance diagonal); `p_adjusted` is
#' [stats::p.adjust()] with the stated method over the tests the table
#' carries. Intervals for a probability strictly inside (0, 1) are formed on
#' the logit scale, for a positive parameter on the log scale (standard error
#' carried by the delta method, `se / (p (1 - p))` and `se / v`), and
#' otherwise as `estimate +/- z se`.
#'
#' @param estimate,standard_error Numeric vectors.
#' @param kind `"probability"`, `"positive"` or `"real"`, per row.
#' @param level Confidence level.
#' @param adjust A [stats::p.adjust()] method.
#' @return A `data.frame` with `statistic`, `p_value`, `p_adjusted`,
#'   `conf_low` and `conf_high`.
rse_wald <- function(estimate, standard_error, kind, level, adjust) {
  bounded <- kind != "real"
  statistic <- ifelse(bounded, NA_real_, estimate / standard_error)
  p_value <- 2 * stats::pnorm(-abs(statistic))
  tested <- !is.na(p_value)
  p_adjusted <- rep(NA_real_, length(p_value))
  p_adjusted[tested] <- stats::p.adjust(p_value[tested], method = adjust)
  z <- stats::qnorm((1 + level) / 2)
  # Each scale is applied only to the rows it is defined for, so that no log
  # of a negative estimate is ever taken.
  conf_low <- estimate - z * standard_error
  conf_high <- estimate + z * standard_error
  logit <- kind == "probability" & estimate > 0 & estimate < 1
  centre <- stats::qlogis(estimate[logit])
  spread <- z * standard_error[logit] / (estimate[logit] * (1 - estimate[logit]))
  conf_low[logit] <- stats::plogis(centre - spread)
  conf_high[logit] <- stats::plogis(centre + spread)
  logarithm <- kind == "positive" & estimate > 0
  spread <- z * standard_error[logarithm] / estimate[logarithm]
  conf_low[logarithm] <- estimate[logarithm] * exp(-spread)
  conf_high[logarithm] <- estimate[logarithm] * exp(spread)
  data.frame(statistic = statistic, p_value = p_value, p_adjusted = p_adjusted,
             conf_low = conf_low, conf_high = conf_high)
}

#' A `parameter_inference()` table with its keys declared
#'
#' The group-class probability has no term (`NA`). The harness joins keys
#' with `merge()`, which pairs `NA` with `NA`, but counts a row with a
#' missing key as absent from the fit, so the missing term is spelled
#' `"<none>"` here, on the obtained side exactly as on the expected one. No
#' value column is touched.
#'
#' @param table The verb's result.
#' @return A [keyed_table()].
rse_keyed_inference <- function(table) {
  plain <- as.data.frame(table)
  plain$term <- ifelse(is.na(plain$term), "<none>", plain$term)
  keyed_table(plain, c("level", "outcome", "term", "parameter"))
}

#' The expected `parameter_inference()` table
#'
#' `estimate` and `standard_error` come from the reference alone. The Wald
#' columns are [rse_wald()] applied to the obtained `estimate` and
#' `standard_error`: they check the documented arithmetic given the package's
#' standard errors, which the `standard_error` column checks independently.
#'
#' @param model An `rse_model`.
#' @param theta The reference's maximum.
#' @param covariance The reference covariance of `theta` for this `vcov_type`.
#' @param obtained The verb's result.
#' @param level,adjust The arguments the verb was called with.
#' @return A `data.frame` keyed like [rse_keyed_inference()].
rse_expected_inference <- function(model, theta, covariance, obtained, level, adjust) {
  labels <- model$labels
  expected <- data.frame(labels[c("level", "outcome", "term", "parameter")],
                         estimate = rse_natural(model, theta),
                         standard_error = rse_standard_errors(model, theta, covariance))
  expected$term <- ifelse(is.na(expected$term), "<none>", expected$term)
  at <- match(.rse_key(labels), .rse_key(obtained))
  stopifnot("the verb must report every parameter of the reference" = !anyNA(at))
  cbind(expected, rse_wald(obtained$estimate[at], obtained$standard_error[at],
                           labels$kind, level, adjust))
}

#' The reference in the harness's shape, for comparing the whole fit
#'
#' @param model An `rse_model`.
#' @param theta The reference's maximum.
#' @param data The list from [rse_data()].
#' @param frame The indicator data, in the fit's row order.
#' @param covariance The observed covariance of `theta`, for the measurement
#'   standard errors the fit reports.
#' @param groups Each row's group label, for a two-level fit.
#' @return A `harness_reference`.
rse_fit_reference <- function(model, theta, data, frame, covariance, groups = NULL) {
  evaluated <- rse_evaluate(model, theta, data)
  parameters <- evaluated$parameters
  errors <- rse_standard_errors(model, theta, covariance)
  labels <- model$labels
  k <- model$n_profiles
  j <- length(model$indicators)
  on_diagonal <- labels$parameter == "variance" |
    (labels$parameter == "covariance" & labels$kind == "positive")
  # The diagonal cells of a full covariance are in lower-triangle order;
  # re-order them by indicator.
  diagonal_order <- if (model$covariance == "diagonal") seq_len(j) else
    order(model$cells[model$cells[, 1L] == model$cells[, 2L], 1L])
  variance_errors <- matrix(errors[on_diagonal], k, j, byrow = TRUE)[, diagonal_order, drop = FALSE]
  standard_errors <- list(
    mean = matrix(errors[labels$parameter == "mean"], k, j, byrow = TRUE,
                  dimnames = list(NULL, model$indicators)),
    variance = matrix(variance_errors, k, j, dimnames = list(NULL, model$indicators)))
  membership <- if (!model$covariates) {
    list(profile_probabilities = parameters$gamma, group_probabilities = parameters$delta)
  } else {
    # Full columns, the baseline class's logits zero; a model with one group
    # class has no group block at all.
    list(coefficients = Filter(Negate(is.null), list(
      profile = cbind(parameters$profile_logits, 0),
      group = if (model$n_group_classes > 1L) cbind(parameters$group_logits, 0))))
  }
  do.call(harness_reference, c(list(
    log_likelihood = sum(evaluated$group_log_likelihood),
    n_parameters = sum(labels$free),
    posteriors = evaluated$posteriors, groups = groups,
    means = parameters$means, covariances = parameters$covariances,
    standard_errors = standard_errors,
    group_posteriors = if (model$n_group_classes > 1L) evaluated$group_posteriors,
    group_log_likelihood = evaluated$group_log_likelihood,
    data = frame,
    row_conditional_entropy = if (model$n_group_classes > 1L)
      evaluated$row_conditional_entropy), membership))
}
