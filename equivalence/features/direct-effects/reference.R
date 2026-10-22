# An independent implementation of the profile model with direct effects,
# written from the formula and calling nothing from latents:
#
#   f(y_i | z_i) = sum_k pi_k prod_j N(y_ij; mu_kj + sum_c b_jc z_ic, s2_kj)
#
# (s2_kj = s2_j under EEI). ref_log_likelihood() evaluates it with dnorm()
# and log-sum-exp; ref_maximize() maximizes it by EM in which each indicator's
# class intercepts and effects are ONE weighted least-squares problem on the
# rows stacked over classes (an exact joint M-step given the variances, unlike
# latents' coordinate-wise ECM), from k-means and random-partition starts;
# ref_standard_errors() inverts numDeriv::hessian() of the log likelihood in
# the natural parameters. Sourced by equivalence/test-features-direct-effects.R.
#
# A parameter point is a list: means (k x p, at covariates 0), variances
# (k x p; equal rows under EEI), effects (a data frame indicator, covariate,
# b) and proportions (k).

#' Rows by indicators: the covariate part sum_c b_jc z_ic
ref_shift <- function(data, vars, effects) {
  vapply(vars, function(indicator) {
    rows <- effects$indicator == indicator
    if (!any(rows)) return(numeric(nrow(data)))
    as.vector(as.matrix(data[effects$covariate[rows]]) %*% effects$b[rows])
  }, numeric(nrow(data)))
}

#' Rows by classes: log pi_k + sum_j log N(.)
ref_joint <- function(point, data, vars) {
  residual <- as.matrix(data[vars]) - ref_shift(data, vars, point$effects)
  k <- length(point$proportions)
  vapply(seq_len(k), function(j) {
    log(point$proportions[[j]]) +
      rowSums(stats::dnorm(residual, matrix(point$means[j, ], nrow(data), length(vars),
                                            byrow = TRUE),
                           matrix(sqrt(point$variances[j, ]), nrow(data), length(vars),
                                  byrow = TRUE), log = TRUE))
  }, numeric(nrow(data)))
}

ref_log_sum_exp <- function(joint) {
  top <- apply(joint, 1L, max)
  top + log(rowSums(exp(joint - top)))
}

#' The log likelihood of a point
ref_log_likelihood <- function(point, data, vars) {
  sum(ref_log_sum_exp(matrix(ref_joint(point, data, vars), nrow(data))))
}

#' The posteriors of a point, rows by classes
ref_posteriors <- function(point, data, vars) {
  joint <- matrix(ref_joint(point, data, vars), nrow(data))
  exp(joint - ref_log_sum_exp(joint))
}

#' One EM M-step from posteriors, given the current variances
ref_m_step <- function(posteriors, point, data, vars, model) {
  k <- ncol(posteriors)
  n <- nrow(data)
  y <- as.matrix(data[vars])
  effects <- point$effects
  # Each indicator: weighted least squares of y on class dummies and its
  # covariates, rows stacked over classes with weights r_ik / s2_kj.
  fitted <- lapply(vars, function(indicator) {
    j <- match(indicator, vars)
    covariates <- effects$covariate[effects$indicator == indicator]
    dummies <- kronecker(diag(k), matrix(1, n, 1L))
    design <- cbind(dummies, if (length(covariates)) do.call(rbind, rep(list(
      as.matrix(data[covariates])), k)))
    weights <- as.vector(sweep(posteriors, 2L, point$variances[, j], "/"))
    response <- rep(y[, j], k)
    coefficients <- qr.solve(crossprod(design, design * weights),
                             crossprod(design, response * weights))
    list(means = coefficients[seq_len(k)], b = coefficients[-seq_len(k)])
  })
  point$means <- vapply(fitted, `[[`, numeric(k), "means")
  point$means <- matrix(point$means, k, length(vars), dimnames = list(NULL, vars))
  point$effects$b <- unlist(lapply(vars, function(indicator) {
    fitted[[match(indicator, vars)]]$b
  }))[order(unlist(lapply(vars, \(indicator) which(effects$indicator == indicator))))]
  residual <- y - ref_shift(data, vars, point$effects)
  totals <- colSums(posteriors)
  squares <- vapply(seq_len(k), function(j) {
    colSums(posteriors[, j] * sweep(residual, 2L, point$means[j, ])^2)
  }, numeric(length(vars)))
  squares <- matrix(squares, length(vars), k)
  point$variances <- if (identical(model, "EEI") || k == 1L) {
    matrix(rowSums(squares) / n, k, length(vars), byrow = TRUE)
  } else t(sweep(squares, 2L, totals, "/"))
  colnames(point$variances) <- vars
  point$proportions <- totals / n
  point
}

#' EM from one start (posteriors), to a relative change of `tol`
ref_em <- function(posteriors, data, vars, effects, model, tol = 1e-12, max_iter = 20000L) {
  k <- ncol(posteriors)
  y <- as.matrix(data[vars])
  point <- list(means = matrix(0, k, length(vars)),
                variances = matrix(apply(y, 2L, stats::var), k, length(vars), byrow = TRUE),
                effects = transform(effects, b = 0), proportions = colMeans(posteriors))
  # Two passes from the start's posteriors settle means, effects and variances.
  point <- ref_m_step(posteriors, point, data, vars, model)
  point <- ref_m_step(posteriors, point, data, vars, model)
  state <- list(point = point, log_likelihood = -Inf, iterations = 0L, converged = FALSE)
  # EM is iterative by nature: each step needs the previous point.
  repeat {
    posteriors <- ref_posteriors(state$point, data, vars)
    if (any(!is.finite(posteriors))) break
    next_point <- ref_m_step(posteriors, state$point, data, vars, model)
    if (any(next_point$variances < 1e-10 * mean(apply(y, 2L, stats::var)))) break
    value <- ref_log_likelihood(next_point, data, vars)
    change <- value - state$log_likelihood
    state <- list(point = next_point, log_likelihood = value,
                  iterations = state$iterations + 1L, converged = FALSE)
    if (abs(change) < tol * abs(value)) {
      state$converged <- TRUE
      break
    }
    if (state$iterations >= max_iter) break
  }
  state
}

#' The best EM solution over k-means and random-partition starts
ref_maximize <- function(data, vars, effects, k, model, n_starts = 12L, seed = 1L) {
  y <- as.matrix(data[vars])
  starts <- local({
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", globalenv())
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    # Three kinds of start in turn: k-means on the standardized indicators, a
    # random hard partition, random soft posteriors.
    lapply(seq_len(n_starts), function(start) {
      if (k == 1L) return(matrix(1, nrow(y), 1L))
      if (start %% 3L == 0L) {
        draws <- matrix(stats::rexp(nrow(y) * k), nrow(y), k)
        return(draws / rowSums(draws))
      }
      groups <- if (start %% 3L == 1L) {
        stats::kmeans(scale(y), k, nstart = 1L, iter.max = 50L)$cluster
      } else sample.int(k, nrow(y), replace = TRUE)
      posteriors <- matrix(0.1 / (k - 1), nrow(y), k)
      posteriors[cbind(seq_len(nrow(y)), groups)] <- 0.9
      posteriors
    })
  })
  runs <- lapply(starts, ref_em, data = data, vars = vars, effects = effects, model = model)
  runs <- Filter(\(run) is.finite(run$log_likelihood), runs)
  values <- vapply(runs, `[[`, numeric(1), "log_likelihood")
  best <- runs[[which.max(values)]]
  best$all <- values
  best
}

#' Free natural parameters of a point as a named vector: means, variances
#' (one per indicator under EEI), effects, and logits of classes 2..k
#' against class 1
ref_pack <- function(point, model) {
  k <- length(point$proportions)
  vars <- colnames(point$means) %||% paste0("v", seq_len(ncol(point$means)))
  label <- \(j, v) sprintf("mean:%s:%d", v, j)
  means <- stats::setNames(as.vector(t(point$means)),
                           as.vector(t(outer(seq_len(k), vars, label))))
  variances <- if (identical(model, "EEI") || k == 1L) {
    stats::setNames(point$variances[1L, ], sprintf("variance:%s:shared", vars))
  } else stats::setNames(as.vector(t(point$variances)),
                         as.vector(t(outer(seq_len(k), vars,
                                           \(j, v) sprintf("variance:%s:%d", v, j)))))
  effects <- stats::setNames(point$effects$b, sprintf("effect:%s:%s", point$effects$indicator,
                                                      point$effects$covariate))
  logits <- stats::setNames(log(point$proportions[-1L] / point$proportions[[1L]]),
                            sprintf("logit:%d", seq_len(k)[-1L]))
  c(means, variances, effects, logits)
}

ref_unpack <- function(values, template, model) {
  k <- length(template$proportions)
  p <- ncol(template$means)
  point <- template
  point$means[] <- matrix(values[seq_len(k * p)], k, p, byrow = TRUE)
  at <- k * p
  if (identical(model, "EEI") || k == 1L) {
    point$variances[] <- matrix(values[at + seq_len(p)], k, p, byrow = TRUE)
    at <- at + p
  } else {
    point$variances[] <- matrix(values[at + seq_len(k * p)], k, p, byrow = TRUE)
    at <- at + k * p
  }
  point$effects$b <- values[at + seq_len(nrow(template$effects))]
  at <- at + nrow(template$effects)
  logits <- c(0, values[at + seq_len(k - 1L)])
  point$proportions <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
  point
}

#' Standard errors from numDeriv's Hessian of the independent log likelihood,
#' with the class proportions' by the delta method from the logits
#'
#' numDeriv steps relative to each coordinate's value, which for a mean near
#' zero is far below the likelihood's rounding noise. The Hessian is therefore
#' taken in units of each coordinate's rough standard error (from a central
#' second difference), u = (theta - estimate) / scale, at u = 0, where numDeriv
#' steps by `eps`; the gradient in those units is the score per standard
#' error.
ref_standard_errors <- function(point, data, vars, model) {
  values <- ref_pack(point, model)
  objective <- function(v) ref_log_likelihood(ref_unpack(v, point, model), data, vars)
  centre <- objective(values)
  scale <- vapply(seq_along(values), function(i) {
    step <- 1e-3 * max(abs(values[[i]]), 0.1)
    up <- values
    up[[i]] <- up[[i]] + step
    down <- values
    down[[i]] <- down[[i]] - step
    curvature <- (objective(up) - 2 * centre + objective(down)) / step^2
    1 / sqrt(abs(curvature))
  }, numeric(1))
  standardized <- function(u) objective(values + u * scale)
  zero <- numeric(length(values))
  # Steps of 0.1 to 0.0125 standard errors: the likelihood moves by 5e-3 to
  # 8e-5, far above its rounding noise; Richardson removes the truncation.
  settings <- list(eps = 0.1, d = 0.1, r = 4, zero.tol = 1e-8)
  hessian <- numDeriv::hessian(standardized, zero, method.args = settings)
  covariance <- solve(-hessian) * outer(scale, scale)
  dimnames(covariance) <- list(names(values), names(values))
  errors <- stats::setNames(sqrt(diag(covariance)), names(values))
  is_logit <- startsWith(names(values), "logit:")
  k <- length(point$proportions)
  if (k > 1L) {
    proportions <- point$proportions
    # d pi_j / d logit_m (m = 2..k) = pi_j (1[j = m] - pi_m)
    jacobian <- vapply(seq_len(k - 1L) + 1L, \(m) proportions * ((seq_len(k) == m) -
                                                                  proportions[[m]]),
                       numeric(k))
    jacobian <- matrix(jacobian, k, k - 1L)
    proportion_errors <- sqrt(diag(jacobian %*% covariance[is_logit, is_logit, drop = FALSE] %*%
                                     t(jacobian)))
    errors <- c(errors[!is_logit], stats::setNames(proportion_errors,
                                                   sprintf("proportion:%d", seq_len(k))))
  } else errors <- errors[!is_logit]
  score <- numDeriv::grad(standardized, zero, method.args = settings)
  list(errors = errors, scaled_score = stats::setNames(score, names(values)),
       values = values)
}

#' A latents fit as a reference point (no latents function is called: the
#' fit's stored fields are read)
ref_point_from_fit <- function(fit, effects) {
  b <- fit$direct_effects[sprintf("%s:%s", effects$indicator, effects$covariate)]
  effects$b <- unname(b)
  list(means = unclass(fit$means), variances = unclass(fit$variances), effects = effects,
       proportions = as.vector(fit$profile_probabilities))
}

#' A cheap dataset for the reference battery, drawn with simulab
ref_dataset <- function(seed, n, k, model, effects, covariate_kind = "normal") {
  means <- rbind(c(0, 0, 0), c(2, 1.6, 1.8), c(-1.6, 1.8, -1.2))[seq_len(max(k, 2L)), ]
  if (k == 1L) means[2L, ] <- means[1L, ]
  sds <- if (identical(model, "EEI")) c(1, 1.2, 0.8) else
    rbind(c(1, 1, 1), c(0.8, 1.2, 0.9), c(1.2, 0.9, 1))[seq_len(max(k, 2L)), ]
  proportions <- c(0.5, 0.3, 0.2)[seq_len(max(k, 2L))]
  simulated <- simulab::simulate_lpa(n = n, means = means, sds = sds,
                                     proportions = proportions, seed = seed)
  frame <- as.data.frame(simulated)
  y <- as.matrix(frame[paste0("indicator_", 1:3)])
  colnames(y) <- paste0("y", 1:3)
  covariates <- unique(effects$covariate)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed + 5L)
  z <- vapply(covariates, \(name) switch(covariate_kind,
                                         normal = stats::rnorm(n, 0.5, 1.5),
                                         binary = as.numeric(stats::runif(n) < 0.4)),
              numeric(n))
  z <- matrix(z, n, length(covariates), dimnames = list(NULL, covariates))
  data <- data.frame(id = seq_len(n), z, y)
  data[colnames(y)] <- as.matrix(data[colnames(y)]) + ref_shift(data, colnames(y), effects)
  data
}

#' The reference battery's design: one row per dataset
ref_design <- function() {
  effect_sets <- list(
    one = data.frame(indicator = "y3", covariate = "z", b = 0.8),
    two_covariates = data.frame(indicator = "y3", covariate = c("z", "w"), b = c(0.8, -0.5)),
    two_indicators = data.frame(indicator = c("y3", "y1"), covariate = "z", b = c(0.8, -0.6)))
  grid <- expand.grid(n = c(60L, 200L, 600L, 2000L), k = 1:3, model = c("VVI", "EEI"),
                      effects = names(effect_sets), stringsAsFactors = FALSE)
  # k = 1 is the same model under VVI and EEI: keep one.
  grid <- grid[!(grid$k == 1L & grid$model == "EEI"), , drop = FALSE]
  replicate_two <- grid[grid$n %in% c(200L, 600L) & grid$k > 1L, , drop = FALSE]
  grid <- rbind(transform(grid, replicate = 1L), transform(replicate_two, replicate = 2L))
  grid$covariate_kind <- ifelse(seq_len(nrow(grid)) %% 3L == 0L, "binary", "normal")
  grid$seed <- 7000L + seq_len(nrow(grid)) * 13L
  grid$name <- sprintf("ref_%s_k%d_n%d_%s_%d", tolower(grid$model), grid$k, grid$n,
                       grid$effects, grid$replicate)
  rownames(grid) <- NULL
  list(grid = grid, effect_sets = effect_sets)
}

#' Everything the battery compares for one dataset of the design
ref_check <- function(row, effect_sets) {
  effects <- effect_sets[[row$effects]]
  data <- ref_dataset(row$seed, row$n, row$k, row$model, effects, row$covariate_kind)
  vars <- paste0("y", 1:3)
  terms <- vapply(unique(effects$indicator), \(indicator) sprintf(
    "%s ~ %s", indicator, paste(effects$covariate[effects$indicator == indicator],
                                collapse = " + ")), character(1), USE.NAMES = FALSE)
  fit <- suppressMessages(lpa(data, vars, row$k, model = row$model, n_starts = 20L, seed = 1L,
                              tol = 1e-10, max_iter = 5000L, direct_effects = terms))
  point <- ref_point_from_fit(fit, effects)
  at_latents <- ref_log_likelihood(point, data, vars)
  independent <- ref_maximize(data, vars, effects, row$k, row$model)
  summary <- data.frame(name = row$name, n = row$n, k = row$k, model = row$model,
                        effects = row$effects, covariate_kind = row$covariate_kind,
                        latents = fit$log_likelihood, reference_at_latents = at_latents,
                        reference_max = independent$log_likelihood,
                        reference_converged = independent$converged,
                        latents_converged = isTRUE(fit$converged),
                        boundary = isTRUE(fit$boundary), max_relative_se = NA_real_,
                        worst_se = NA_character_, max_scaled_score = NA_real_)
  # At an active variance bound the information is not that of an interior
  # maximum; neither program's standard errors are defined there.
  if (summary$boundary) return(summary)
  errors <- ref_standard_errors(point, data, vars, row$model)
  inference <- parameter_inference(fit)
  keyed <- function(table) {
    stats::setNames(table$standard_error, ifelse(
      table$parameter == "direct_effect", sprintf("effect:%s", table$term),
      ifelse(table$parameter == "probability", sprintf("proportion:%s", sub("profile_", "",
                                                                         table$outcome)),
             sprintf("%s:%s:%s", table$parameter, table$term,
                     ifelse(table$outcome == "shared" | (table$parameter == "variance" &
                                                          (row$model == "EEI" | row$k == 1L)),
                            "shared", sub("profile_", "", table$outcome))))))
  }
  mine <- keyed(subset(inference, level %in% c("measurement", "profile")))
  mine <- mine[!duplicated(names(mine))]
  if (row$k == 1L) mine <- mine[!startsWith(names(mine), "proportion")]
  compared <- intersect(names(errors$errors), names(mine))
  stopifnot("the reference and latents name different parameters" =
              setequal(compared, names(errors$errors)))
  relative_se <- abs(mine[compared] - errors$errors[compared]) / errors$errors[compared]
  scaled_score <- abs(errors$scaled_score)
  summary$max_relative_se <- max(relative_se)
  summary$worst_se <- names(which.max(relative_se))
  summary$max_scaled_score <- max(scaled_score)
  summary
}
