# Membership-covariate latent class models, written from their definitions in
# base R.
#
# Sourced by suites/covariates.R. No program the covariate suite compares
# against reports every quantity the package does: poLCA, glca, multilevLCA and
# depmixS4 give no observed-information standard errors of the membership
# coefficients and the response probabilities together (poLCA's are an outer
# product, multilevLCA's are its own two-step corrected ones, glca's scores are
# documented wrong with missing data, depmixS4 differences with fdHess), and
# glca and multilevLCA give no within-class row entropies. What they do not
# report is computed here, at their estimates, from the model's definition
# (Vermunt 2003; Vermunt 2010; Lyrvall et al. 2025):
#
#   f(y_g | w_g) = sum_h P(h | w_g) prod_{i in g} sum_k P(k | h, x_i) f_k(y_i)
#   P(h | w_g)    = exp(w_g' a_h) / sum_h' exp(w_g' a_h')
#   P(k | h, x_i) = exp(x_i' b_hk) / sum_k' exp(x_i' b_hk')
#
# with f_k a product over indicators (categorical probabilities, Poisson
# means, Gaussian means and variances), a missing indicator contributing a
# factor of one. A single-level model is the special case of one group class
# with each row its own group. Every program-reported quantity this code can
# also produce (likelihood, posteriors, group posteriors, group likelihoods) is
# checked against the program's own value before it is trusted, so a slip here
# cannot pass silently.
#
# Standard errors are the inverse of minus a Richardson-extrapolated numDeriv
# Hessian of `covariate_evaluate()`'s log-likelihood, in a chart where every
# parameter is unbounded and the covariates are centred and scaled: response
# probabilities as logits against the last category, Poisson means as logs,
# Gaussian variances as logs, membership logits against the last class. The
# raw-unit chart is not used: a slope on a covariate far from zero is nearly
# collinear with its intercept, and on the probe of 2026-10-27 (p11) logits
# near 60 there gave a NaN information. At a maximum the inverse information
# transforms with the Jacobian of any reparameterisation, so `J V J'`, with `J`
# the numDeriv Jacobian of the natural parameters (probabilities, means,
# variances, logits against the last class in the covariates' own units) with
# respect to the chart, is the covariance of the natural parameters.
#
# Vermunt, J. K. (2003). Multilevel latent class models. Sociological
#   Methodology, 33(1), 213-239.
# Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#   three-step approaches. Political Analysis, 18(4), 450-469.
# Lyrvall, J., Di Mari, R., Bakk, Z., Oser, J., & Kuha, J. (2025). Multilevel
#   latent class analysis: State-of-the-art methodologies and their
#   implementation in the R package multilevLCA. Multivariate Behavioral
#   Research.

#' The data a covariate model is evaluated on
#'
#' @param data A data frame holding the indicators and covariates.
#' @param groups Group label of each row; `NULL` for a single-level model
#'   (each row its own group).
#' @param levels Named list of category labels of the categorical indicators,
#'   in the fit's category order.
#' @param counts,gaussian Names of the Poisson and Gaussian indicators.
#' @param profile_covariates,group_covariates Names of the numeric covariates
#'   of profile and of group-class membership.
#' @return A list: `codes` (rows by categorical indicators), `counts`,
#'   `gaussian` (matrices of the observed values), `profile_design` (rows by
#'   `(Intercept)` and the profile covariates), `group_design` (groups by
#'   `(Intercept)` and the group covariates), `group_index`, `group_ids` (in
#'   first-occurrence order) and the centring constants of both designs.
covariate_context <- function(data, groups = NULL, levels = list(), counts = character(),
                              gaussian = character(), profile_covariates = character(),
                              group_covariates = character()) {
  stopifnot("`data` must be a data frame" = is.data.frame(data),
            "every indicator and covariate must be a column of `data`" =
              all(c(names(levels), counts, gaussian, profile_covariates,
                    group_covariates) %in% names(data)))
  groups <- as.character(groups %||% seq_len(nrow(data)))
  group_ids <- unique(groups)
  group_index <- match(groups, group_ids)
  first_rows <- match(group_ids, groups)
  codes <- vapply(names(levels), function(item) {
    match(as.character(data[[item]]), levels[[item]])
  }, integer(nrow(data)))
  codes <- matrix(codes, nrow(data), dimnames = list(NULL, names(levels)))
  stopifnot("every observed response must be one of the given categories" =
              identical(unname(is.na(codes)), unname(is.na(as.matrix(data[names(levels)])))))
  as_matrix <- function(columns) {
    matrix(as.numeric(as.matrix(data[columns])), nrow(data),
           dimnames = list(NULL, columns))
  }
  design <- function(rows, columns) {
    cbind(`(Intercept)` = 1, as_matrix(columns)[rows, , drop = FALSE])
  }
  group_design <- design(first_rows, group_covariates)
  stopifnot("a group covariate must be constant within each group" =
              all(vapply(group_covariates, function(column) {
                all(tapply(data[[column]], group_index, function(v) length(unique(v)) == 1L))
              }, logical(1))))
  list(codes = codes, levels = levels, counts = as_matrix(counts),
       gaussian = as_matrix(gaussian),
       profile_design = design(seq_len(nrow(data)), profile_covariates),
       group_design = group_design, group_index = group_index, group_ids = group_ids,
       profile_centring = .covariate_centring(design(seq_len(nrow(data)), profile_covariates)),
       group_centring = .covariate_centring(group_design))
}

#' The linear map from a centred and scaled design to the raw one
#'
#' With `D_std = D_raw A` (each covariate column centred on its mean and
#' divided by its standard deviation, the intercept kept), the logits
#' `D_raw B_raw = D_std B_std` give `B_raw = A B_std`.
#'
#' @param design A design matrix whose first column is the intercept.
#' @return The square matrix `A`.
.covariate_centring <- function(design) {
  p <- ncol(design)
  centre <- colMeans(design)[-1L]
  spread <- apply(design, 2L, stats::sd)[-1L]
  map <- diag(p)
  if (p > 1L) {
    map[1L, -1L] <- -centre / spread
    map[-1L, -1L] <- diag(1 / spread, p - 1L)
  }
  dimnames(map) <- list(colnames(design), colnames(design))
  map
}

#' Row-wise log-sum-exp
#' @param m A numeric matrix.
#' @return One value per row.
.covariate_row_lse <- function(m) {
  top <- apply(m, 1L, max)
  top + log(rowSums(exp(m - top)))
}

#' Row-wise log softmax
#' @param m A numeric matrix of logits.
#' @return The log probabilities, same shape.
.covariate_log_softmax <- function(m) m - .covariate_row_lse(m)

#' log f_k(y_i), the measurement density of every row in every profile
#'
#' @param model A covariate model (see [covariate_evaluate()]).
#' @param context From [covariate_context()].
#' @return A rows-by-profiles matrix; a missing indicator adds zero.
.covariate_log_density <- function(model, context) {
  n <- length(context$group_index)
  k <- ncol(model$profile[[1L]])
  observed_only <- function(values, column_values) {
    out <- matrix(0, n, k)
    seen <- !is.na(column_values)
    out[seen, ] <- values[seen, , drop = FALSE]
    out
  }
  categorical <- lapply(colnames(context$codes), function(item) {
    code <- context$codes[, item]
    values <- matrix(0, n, k)
    values[!is.na(code), ] <- t(log(model$responses[[item]]))[code[!is.na(code)], , drop = FALSE]
    observed_only(values, code)
  })
  counts <- lapply(colnames(context$counts), function(item) {
    y <- context$counts[, item]
    lambda <- model$count_means[, item]
    values <- outer(ifelse(is.na(y), 0, y), log(lambda)) -
      matrix(lambda, n, k, byrow = TRUE) - lgamma(ifelse(is.na(y), 0, y) + 1)
    observed_only(values, y)
  })
  gaussian <- lapply(colnames(context$gaussian), function(item) {
    y <- context$gaussian[, item]
    centred <- outer(ifelse(is.na(y), 0, y), model$means[, item], "-")
    values <- -0.5 * (log(2 * pi) + matrix(log(model$variances[, item]), n, k, byrow = TRUE) +
                        centred^2 / matrix(model$variances[, item], n, k, byrow = TRUE))
    observed_only(values, y)
  })
  Reduce(`+`, c(categorical, counts, gaussian), matrix(0, n, k))
}

#' Evaluate a covariate model at given parameters
#'
#' @param model A list: `responses` (named list of profiles-by-categories
#'   probability matrices), `count_means` (profiles by Poisson indicators),
#'   `means` and `variances` (profiles by Gaussian indicators), `profile` (one
#'   design-columns-by-profiles logit matrix per group class, rows named as the
#'   profile design's columns) and `group` (group-design-columns by group
#'   classes, `NULL` with one group class). Logits are identified up to a
#'   constant per design column, so any baseline will do.
#' @param context From [covariate_context()].
#' @return A list: `log_likelihood`; `group_log_likelihood`, named by group;
#'   `group_posteriors`, groups by group classes with group row names;
#'   `posteriors`, the rows' marginal profile posteriors
#'   `P(k | y_g) = sum_h P(h | y_g) P(k | h, y_i, x_i)`; and
#'   `row_conditional_entropy`, `sum_h P(h | y_g) H(k | h, y_i, x_i)`.
covariate_evaluate <- function(model, context) {
  log_density <- .covariate_log_density(model, context)
  n_group_classes <- length(model$profile)
  # Given its group class, a row's profile is independent of the other rows.
  joint <- lapply(model$profile, function(logits) {
    .covariate_log_softmax(context$profile_design %*% logits) + log_density
  })
  row_mix <- vapply(joint, .covariate_row_lse, numeric(nrow(log_density)))
  row_mix <- matrix(row_mix, nrow(log_density), n_group_classes)
  per_group <- rowsum(row_mix, context$group_index, reorder = TRUE)
  log_prior <- if (n_group_classes == 1L) matrix(0, nrow(per_group), 1L) else
    .covariate_log_softmax(context$group_design %*% model$group)
  group_joint <- per_group + log_prior
  group_log_likelihood <- .covariate_row_lse(group_joint)
  group_posteriors <- exp(group_joint - group_log_likelihood)
  row_group <- group_posteriors[context$group_index, , drop = FALSE]
  within <- lapply(seq_len(n_group_classes), function(h) exp(joint[[h]] - row_mix[, h]))
  posteriors <- Reduce(`+`, lapply(seq_len(n_group_classes), function(h) {
    within[[h]] * row_group[, h]
  }))
  row_conditional_entropy <- Reduce(`+`, lapply(seq_len(n_group_classes), function(h) {
    q <- within[[h]]
    row_group[, h] * -rowSums(ifelse(q > 0, q * log(q), 0))
  }))
  dimnames(group_posteriors) <- list(context$group_ids, NULL)
  list(log_likelihood = sum(group_log_likelihood),
       group_log_likelihood = stats::setNames(as.vector(group_log_likelihood),
                                              context$group_ids),
       group_posteriors = group_posteriors, posteriors = posteriors,
       row_conditional_entropy = row_conditional_entropy)
}

#' Put a covariate model's classes in another order
#'
#' @param model A covariate model.
#' @param profiles Permutation of the profiles (`new[k] = old[profiles[k]]`).
#' @param group_classes Permutation of the group classes.
#' @return The permuted model.
covariate_permute <- function(model, profiles, group_classes = seq_along(model$profile)) {
  rows <- function(m) if (is.null(m)) NULL else m[profiles, , drop = FALSE]
  model$responses <- lapply(model$responses, rows)
  model$count_means <- rows(model$count_means)
  model$means <- rows(model$means)
  model$variances <- rows(model$variances)
  model$profile <- lapply(model$profile[group_classes], function(m) m[, profiles, drop = FALSE])
  if (!is.null(model$group)) model$group <- model$group[, group_classes, drop = FALSE]
  model
}

#' The unbounded chart the information is taken in
#'
#' @param model A covariate model, the point the chart is centred on.
#' @param context From [covariate_context()].
#' @param slopes `"shared"` (one slope per covariate and profile logit, the
#'   intercepts differing by group class) or `"group_class"` (slopes per group
#'   class too).
#' @return A list: `theta`, the model in chart coordinates, and `unpack`, the
#'   function mapping chart coordinates back to a model.
covariate_chart <- function(model, context, slopes = c("shared", "group_class")) {
  slopes <- match.arg(slopes)
  k <- ncol(model$profile[[1L]])
  n_group_classes <- length(model$profile)
  against_last <- function(m) m[, -ncol(m), drop = FALSE] - m[, ncol(m)]
  # Logits against the last class, in the centred and scaled design.
  profile_std <- lapply(model$profile, function(logits) {
    solve(context$profile_centring, against_last(logits))
  })
  if (identical(slopes, "shared") && n_group_classes > 1L) {
    spread <- max(vapply(profile_std[-1L], function(block) {
      max(abs(block[-1L, , drop = FALSE] - profile_std[[1L]][-1L, , drop = FALSE]))
    }, numeric(1)))
    stopifnot("shared slopes must be the same in every group class" = spread < 1e-8)
  }
  skeleton <- list(
    responses = lapply(model$responses, function(block) {
      log(block[, -ncol(block), drop = FALSE] / block[, ncol(block)])
    }),
    count_means = if (!is.null(model$count_means)) log(model$count_means),
    means = model$means,
    variances = if (!is.null(model$variances)) log(model$variances),
    intercepts = vapply(profile_std, function(block) block[1L, ], numeric(k - 1L)),
    slopes = if (nrow(context$profile_centring) > 1L) {
      if (identical(slopes, "shared")) profile_std[[1L]][-1L, , drop = FALSE] else
        lapply(profile_std, function(block) block[-1L, , drop = FALSE])
    },
    group = if (n_group_classes > 1L) {
      solve(context$group_centring, against_last(model$group))
    })
  skeleton <- Filter(Negate(is.null), skeleton)
  theta <- unlist(utils::as.relistable(skeleton))
  attributes(theta) <- NULL
  template <- model
  unpack <- function(theta) {
    parts <- utils::relist(theta, utils::as.relistable(skeleton))
    out <- template
    out$responses <- lapply(stats::setNames(nm = names(model$responses)), function(item) {
      logits <- cbind(parts$responses[[item]], 0)
      probabilities <- exp(logits - apply(logits, 1L, max))
      probabilities <- probabilities / rowSums(probabilities)
      dimnames(probabilities) <- dimnames(model$responses[[item]])
      probabilities
    })
    if (!is.null(model$count_means)) out$count_means <- exp(parts$count_means)
    if (!is.null(model$means)) out$means <- parts$means
    if (!is.null(model$variances)) out$variances <- exp(parts$variances)
    intercepts <- matrix(parts$intercepts, k - 1L, n_group_classes)
    out$profile <- lapply(seq_len(n_group_classes), function(h) {
      slope_block <- if (is.null(parts$slopes)) NULL else
        if (identical(slopes, "shared")) parts$slopes else parts$slopes[[h]]
      standardised <- rbind(intercepts[, h], slope_block)
      raw <- cbind(context$profile_centring %*% standardised, 0)
      dimnames(raw) <- dimnames(model$profile[[h]])
      raw
    })
    if (n_group_classes > 1L) {
      out$group <- cbind(context$group_centring %*% parts$group, 0)
      dimnames(out$group) <- dimnames(model$group)
    }
    out
  }
  list(theta = theta, unpack = unpack)
}

#' The membership coefficients under the package's term labels
#'
#' With one group class the profile terms are `(Intercept)` and the
#' covariates; with several, the intercepts are `group_class_h` and the slopes
#' the covariates (`"shared"`) or `covariate:group_class_h` (`"group_class"`).
#'
#' @param model A covariate model.
#' @param slopes As in [covariate_chart()].
#' @return A list with `profile` (terms by profiles, full columns) and `group`
#'   (group-design terms by group classes, `NULL` with one group class).
covariate_coefficients <- function(model, slopes = c("shared", "group_class")) {
  slopes <- match.arg(slopes)
  n_group_classes <- length(model$profile)
  k <- ncol(model$profile[[1L]])
  covariates <- rownames(model$profile[[1L]])[-1L]
  profile <- if (n_group_classes == 1L) model$profile[[1L]] else {
    intercepts <- matrix(t(vapply(model$profile, function(m) m[1L, ], numeric(k))),
                         n_group_classes, k)
    rownames(intercepts) <- sprintf("group_class_%d", seq_len(n_group_classes))
    slope_rows <- if (identical(slopes, "shared")) {
      model$profile[[1L]][covariates, , drop = FALSE]
    } else {
      do.call(rbind, lapply(seq_len(n_group_classes), function(h) {
        block <- model$profile[[h]][covariates, , drop = FALSE]
        rownames(block) <- sprintf("%s:group_class_%d", covariates, h)
        block
      }))
    }
    rbind(intercepts, slope_rows)
  }
  list(profile = profile, group = model$group)
}

#' Every natural parameter `parameter_inference()` reports, as a named vector
#'
#' Names are `level|outcome|term|parameter`: response probabilities (every
#' category), Gaussian means and variances, Poisson means, and the membership
#' logits against the last class.
#'
#' @param model A covariate model, in the fit's class order.
#' @param slopes As in [covariate_chart()].
#' @return A named numeric vector.
covariate_natural <- function(model, slopes = c("shared", "group_class")) {
  slopes <- match.arg(slopes)
  k <- ncol(model$profile[[1L]])
  profile_label <- sprintf("profile_%d", seq_len(k))
  measurement <- function(block, parameter, term_of) {
    if (is.null(block)) return(numeric())
    stats::setNames(as.vector(block), sprintf(
      "measurement|%s|%s|%s", rep(profile_label, ncol(block)),
      rep(term_of(colnames(block)), each = k), parameter))
  }
  responses <- unlist(lapply(names(model$responses), function(item) {
    measurement(model$responses[[item]], "response", function(categories) {
      paste0(item, ":", categories)
    })
  }))
  coefficients <- covariate_coefficients(model, slopes)
  logits <- function(block, level, prefix) {
    if (is.null(block)) return(numeric())
    against <- block[, -ncol(block), drop = FALSE] - block[, ncol(block)]
    stats::setNames(as.vector(against), sprintf(
      "%s|%s_%d|%s|coefficient", level, prefix, rep(seq_len(ncol(against)), each = nrow(against)),
      rep(rownames(block), times = ncol(against))))
  }
  c(responses, measurement(model$means, "mean", identity),
    measurement(model$variances, "variance", identity),
    measurement(model$count_means, "count_mean", identity),
    logits(coefficients$profile, "profile", "profile"),
    logits(coefficients$group, "group", "group_class"))
}

#' Observed-information standard errors of every natural parameter
#'
#' @param model A covariate model at a maximum, in the fit's class order.
#' @param context From [covariate_context()].
#' @param slopes As in [covariate_chart()].
#' @param step The first step of the Richardson sequence, absolute, in the
#'   chart (logits, log means, log variances, means of unit-scale data and
#'   coefficients per standard deviation of their covariate).
#' @return A list: `standard_errors`, named like [covariate_natural()];
#'   `n_parameters`, the chart's dimension; `newton_step`, the largest
#'   coordinate of the Newton step `V g` from `model` (the score `g` by
#'   [.covariate_score()]), its distance from the maximum in the chart.
covariate_standard_errors <- function(model, context, slopes = c("shared", "group_class"),
                                      step = 1e-1) {
  require_suite_packages("numDeriv", reason = "the reference observed information")
  slopes <- match.arg(slopes)
  chart <- covariate_chart(model, context, slopes)
  log_likelihood <- function(theta) covariate_evaluate(chart$unpack(theta), context)$log_likelihood
  # One absolute first step for every coordinate (numDeriv's own steps are
  # relative, which shrinks them to rounding noise for a coordinate near zero).
  information <- -numDeriv::hessian(log_likelihood, chart$theta,
                                    method.args = list(d = 0, eps = step, zero.tol = Inf))
  covariance <- chol2inv(chol(information))
  jacobian <- numDeriv::jacobian(function(theta) covariate_natural(chart$unpack(theta), slopes),
                                 chart$theta)
  gradient <- .covariate_score(log_likelihood, chart$theta)
  list(standard_errors = stats::setNames(sqrt(rowSums((jacobian %*% covariance) * jacobian)),
                                         names(covariate_natural(model, slopes))),
       n_parameters = length(chart$theta),
       newton_step = max(abs(covariance %*% gradient)))
}

#' The score, by Richardson extrapolation from one absolute first step of 1e-2
#'
#' numDeriv's default relative step (1e-4) divides the likelihood's rounding
#' (4e-13 at |log L| = 1848) into a score error of order 1e-8; at 1e-2 the
#' rounding contributes 4e-11 and the extrapolation removes the truncation
#' error.
#'
#' @param log_likelihood A function of the chart coordinates.
#' @param theta The point.
#' @return The gradient.
.covariate_score <- function(log_likelihood, theta) {
  numDeriv::grad(log_likelihood, theta, method.args = list(d = 0, eps = 1e-2, zero.tol = Inf))
}

#' Carry a covariate model to the maximum by Newton steps on this likelihood
#'
#' For a program whose own EM cannot be continued from its solution
#' (multilevLCA's runs in compiled code from its own starting values). The
#' observed information at the program's solution is held fixed and the score
#' recomputed at every step, which converges linearly at a rate set by how far
#' the information moves over the step; it is iterated until no chart
#' coordinate moves by more than `tolerance`.
#'
#' @param model A covariate model near a maximum.
#' @param context From [covariate_context()].
#' @param slopes As in [covariate_chart()].
#' @param max_steps Bound on the number of steps; reaching it is an error.
#' @param tolerance Largest chart step at the maximum.
#' @return A list: `model`, the refined model; `steps`; `displacement`, the
#'   largest absolute change of any natural parameter.
covariate_newton <- function(model, context, slopes = c("shared", "group_class"),
                             max_steps = 50L, tolerance = 1e-9) {
  require_suite_packages("numDeriv", reason = "the reference observed information")
  slopes <- match.arg(slopes)
  chart <- covariate_chart(model, context, slopes)
  log_likelihood <- function(theta) covariate_evaluate(chart$unpack(theta), context)$log_likelihood
  information <- -numDeriv::hessian(log_likelihood, chart$theta,
                                    method.args = list(d = 0, eps = 1e-1, zero.tol = Inf))
  factor <- chol(information)
  # With numDeriv's default step the score's rounding error held the Newton
  # step at about 1e-8; see .covariate_score().
  score <- function(theta) .covariate_score(log_likelihood, theta)
  # Reduce() rather than a loop with a break: once the step is below the
  # tolerance every remaining element returns the state untouched.
  state <- Reduce(function(state, index) {
    if (state$settled) return(state)
    move <- backsolve(factor, forwardsolve(t(factor), score(state$theta)))
    list(theta = state$theta + move, settled = max(abs(move)) < tolerance, steps = index)
  }, seq_len(max_steps), init = list(theta = chart$theta, settled = FALSE, steps = 0L))
  stopifnot("Newton refinement must settle within its step bound" = state$settled)
  refined <- chart$unpack(state$theta)
  stopifnot("refining must not lose likelihood" =
              covariate_evaluate(refined, context)$log_likelihood >=
                covariate_evaluate(model, context)$log_likelihood - 1e-9)
  list(model = refined, steps = state$steps,
       displacement = max(abs(covariate_natural(refined, slopes) -
                                covariate_natural(model, slopes))))
}

#' Check a value this code computes against the program's own
#' @param what Description for the error.
#' @param ours,theirs Numeric vectors or matrices of the same shape.
#' @param tolerance Largest absolute difference allowed.
#' @return `NULL`, invisibly; stops with the largest difference otherwise.
covariate_check <- function(what, ours, theirs, tolerance) {
  gap <- max(abs(as.vector(ours) - as.vector(theirs)))
  if (!is.finite(gap) || gap > tolerance) {
    stop(errorCondition(sprintf("%s: this code's value is %.3g from the program's (allowed %.3g)",
                                what, gap, tolerance),
                        class = "multilpa_bad_reference", call = NULL))
  }
  invisible(NULL)
}
