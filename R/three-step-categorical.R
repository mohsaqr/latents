# Step three for a categorical distal outcome. The class-specific probability
# of each category is estimated by the weighting methods (BCH, proportional,
# modal), as the weighted share of that category, or by the maximum-likelihood
# step-three model, whose log likelihood
#
#   L = sum_i log sum_t pi_t D[t, W_i] p_t(y_i)
#
# has the classification errors D and the class sizes pi (the average
# posterior probabilities) fixed and a distribution p_t over the categories
# in each class free (Vermunt 2010; Bakk, Tekle and Vermunt 2013), as for a
# numeric outcome and as in Mplus's DCAT and Latent GOLD's step-three model
# without an equation for the sizes. It keeps every probability in range
# where the BCH shares fall outside [0, 1]. EM reaches the maximum from any
# interior start because L is concave in the probabilities p_t(c).

#' The three-step table of a categorical distal outcome
#'
#' @param pieces The units' assignments, from `.multilpa_level_assignments()`
#'   or the trajectory equivalent.
#' @param codes The outcome as category codes, one per unit.
#' @param categories The category labels, in code order.
#' @param level,method,ci_level,contrast,adjust,vcov_type As in [three_step()].
#' @param correction Optional step-one covariance correction closure.
#' @return The table [three_step()] returns.
#' @noRd
.multilpa_three_step_categorical <- function(pieces, codes, categories, level, method,
                                             ci_level, contrast, adjust, vcov_type,
                                             correction = NULL) {
  n_classes <- pieces$n_classes
  n_categories <- length(categories)
  indicators <- outer(codes, seq_len(n_categories), "==") * 1
  if (identical(method, "ml")) {
    fit <- .multilpa_step_categorical_fit(pieces, indicators)
    probabilities <- fit$probabilities
    # A maximum with a probability of 0 or 1 sits where its logit is infinite:
    # the estimate stands, but the information is singular there.
    bound <- probabilities < 1e-8 | probabilities > 1 - 1e-8
    influence <- if (any(bound)) {
      warning(warningCondition(paste(
        "A class probability of the outcome is on its bound (0 or 1) at the",
        "maximum, so the step-three model has no standard errors, intervals or",
        "tests; the estimates are reported."), class = "latents_boundary", call = NULL))
      matrix(NA_real_, nrow(indicators), n_classes * n_categories)
    } else .multilpa_step_categorical_influence(fit, indicators)
    carrying <- fit$responsibilities
    effective_n <- colSums(fit$responsibilities)
  } else {
    weights <- .multilpa_step_weights(pieces, method)
    totals <- colSums(weights)
    if (any(!is.finite(totals)) || any(totals <= 0)) {
      stop(errorCondition("A class has no positive total outcome weight.",
                          class = "latents_inseparable_classes", call = NULL))
    }
    probabilities <- crossprod(weights, indicators) / totals
    # The influence of unit i on the share of category c in class k, laid out
    # class by class: w_ik (d_ic - p_kc) / T_k.
    influence <- do.call(cbind, lapply(seq_len(n_classes), function(class) {
      weights[, class] * sweep(indicators, 2L, probabilities[class, ], "-") / totals[class]
    }))
    carrying <- weights
    effective_n <- totals^2 / colSums(weights^2)
  }
  clustered <- .multilpa_step_clusters(influence, pieces, vcov_type)
  thin <- .multilpa_step_thin(carrying, pieces, vcov_type)
  clustered[, rep(thin, each = n_categories)] <- NA_real_
  if (!is.null(correction) && !all(is.na(clustered))) {
    if (identical(method, "ml")) {
      free <- seq_len(n_categories - 1L)
      jacobian <- do.call(rbind, lapply(seq_len(n_classes), function(class) {
        p <- probabilities[class, ]
        block <- matrix(0, n_categories, n_classes * length(free))
        columns <- (class - 1L) * length(free) + free
        block[, columns] <- diag(p)[, free, drop = FALSE] - tcrossprod(p, p[free])
        block
      }))
      transform <- jacobian %*% solve(.multilpa_step_categorical_information(fit, indicators))
      clustered <- correction(clustered, function(updated) {
        weights <- t(.multilpa_error_matrix(updated)[, updated$modal, drop = FALSE]) *
          matrix(colMeans(updated$posteriors), length(codes), n_classes, byrow = TRUE)
        joint <- weights * tcrossprod(indicators, probabilities)
        fit$responsibilities <- joint / rowSums(joint)
        colSums(.multilpa_step_categorical_scores(fit, indicators))
      }, transform)
    } else {
      clustered <- correction(clustered, function(updated) {
        weights <- .multilpa_step_weights(updated, method)
        as.vector(t(crossprod(weights, indicators) / colSums(weights)))
      })
    }
  }
  outside <- probabilities < 0 | probabilities > 1
  if (any(outside)) {
    cells <- which(outside, arr.ind = TRUE)
    warning(warningCondition(sprintf(paste(
      "The BCH weights give %s outside [0, 1] (%s): the classes are too poorly",
      "separated for the correction. They are reported without an interval or",
      "test, as is any class difference that involves them;",
      "`method = \"ml\"` keeps the probabilities in range."),
      if (nrow(cells) > 1L) "probabilities" else "a probability",
      paste(sprintf("class %d, %s", cells[, 1L], categories[cells[, 2L]]), collapse = "; ")),
      class = "latents_bch_out_of_range", call = NULL))
  }
  estimates <- as.vector(t(probabilities))
  if (identical(contrast, "overall")) {
    # On the baseline-category logits log(p_c / p_C), as Latent GOLD and Mplus
    # test a class effect on a nominal outcome; the statistic does not depend
    # on which category is the baseline. A share outside (0, 1) has no logit.
    if (any(probabilities <= 0 | probabilities >= 1)) {
      return(.multilpa_step_overall(rep(NA_real_, n_classes * (n_categories - 1L)),
                                    matrix(NA_real_, nrow(clustered),
                                           n_classes * (n_categories - 1L)),
                                    n_classes, level, method, vcov_type))
    }
    logits <- lapply(seq_len(n_classes), function(class) {
      p <- probabilities[class, ]
      columns <- (class - 1L) * n_categories + seq_len(n_categories)
      contributions <- sweep(clustered[, columns, drop = FALSE], 2L, p, "/")
      list(estimates = log(p[-n_categories] / p[n_categories]),
           clustered = contributions[, -n_categories, drop = FALSE] -
             contributions[, n_categories])
    })
    return(.multilpa_step_overall(unlist(lapply(logits, `[[`, "estimates")),
                                  do.call(cbind, lapply(logits, `[[`, "clustered")),
                                  n_classes, level, method, vcov_type))
  }
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  if (identical(contrast, "pairs")) {
    return(.multilpa_step_categorical_pairs(probabilities, clustered, categories, level,
                                            method, quantile, adjust, vcov_type,
                                            invalid = outside))
  }
  errors <- sqrt(colSums(clustered^2))
  # Intervals on the logit scale, back-transformed, so they stay in (0, 1).
  interior <- estimates > 0 & estimates < 1 & is.finite(errors)
  half <- quantile * errors / (estimates * (1 - estimates))
  logit <- stats::qlogis(pmin(pmax(estimates, 0), 1))
  result <- data.frame(
    level = level, method = method, class = rep(seq_len(n_classes), each = n_categories),
    category = rep(categories, n_classes), estimate = estimates, standard_error = errors,
    conf_low = ifelse(interior, stats::plogis(logit - half), NA_real_),
    conf_high = ifelse(interior, stats::plogis(logit + half), NA_real_),
    effective_n = rep(effective_n, each = n_categories),
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "vcov_type") <- vcov_type
  column <- .latents_column
  .latents_table(result, ci_level,
                 sprintf("Distal outcome by class (%s, %s%% CI)", method,
                         format(100 * ci_level)),
                 list(column("Class", "class", "integer"), column("Category", "category", "text"),
                      column("Probability", "estimate"),
                      column(sprintf("%s%% CI", format(100 * ci_level)),
                             c("conf_low", "conf_high"), "ci", "estimate"),
                      column("SE", "standard_error"),
                      column("Effective n", "effective_n")))
}

#' Fit the step-three model of a categorical distal outcome
#'
#' @param pieces The units' assignments: `posteriors`, `modal`, `n_classes`.
#' @param indicators Units-by-categories 0/1 matrix of the outcome.
#' @param tol Largest change in a probability at which EM stops.
#' @param max_iter Maximum number of EM iterations.
#' @return A list: `sizes` (one per class, fixed), `probabilities` (classes
#'   by categories), `responsibilities` (units by classes), `log_likelihood`,
#'   `converged`, `iterations`.
#' @noRd
.multilpa_step_categorical_fit <- function(pieces, indicators, tol = 1e-12,
                                           max_iter = 20000L) {
  # The probability that a unit of each class is assigned where this one was.
  errors <- t(.multilpa_error_matrix(pieces)[, pieces$modal, drop = FALSE])
  sizes <- colMeans(pieces$posteriors)
  e_step <- function(sizes, probabilities) {
    joint <- errors * tcrossprod(indicators, probabilities) *
      matrix(sizes, nrow(indicators), length(sizes), byrow = TRUE)
    marginal <- rowSums(joint)
    list(responsibilities = joint / marginal, log_likelihood = sum(log(marginal)))
  }
  m_step <- function(responsibilities) {
    crossprod(responsibilities, indicators) / colSums(responsibilities)
  }
  # The posterior-weighted shares, pulled towards the overall shares so that
  # no category starts at zero in a class.
  start <- m_step(pieces$posteriors)
  overall <- matrix(colMeans(indicators), nrow(start), ncol(start), byrow = TRUE)
  probabilities <- 0.9 * start + 0.1 * overall
  state <- e_step(sizes, probabilities)
  iterations <- 0L
  converged <- FALSE
  # EM is iterative by nature: each step needs the previous one's estimates,
  # so this loop cannot be vectorised.
  while (iterations < max_iter) {
    iterations <- iterations + 1L
    previous <- probabilities
    probabilities <- m_step(state$responsibilities)
    state <- e_step(sizes, probabilities)
    # EM can creep when the classes overlap: the likelihood stops changing
    # before the probabilities do, so the stop is on the parameters.
    change <- max(abs(probabilities - previous))
    if (change <= tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning(warningCondition(sprintf(paste(
      "The step-three model did not converge in %d iterations, so its class",
      "probabilities are not a maximum."), iterations),
      class = "latents_unconverged", call = NULL))
  }
  list(sizes = sizes, probabilities = m_step(state$responsibilities),
       responsibilities = state$responsibilities, log_likelihood = state$log_likelihood,
       converged = converged, iterations = iterations)
}

#' How each class's log joint density moves with the parameters
#'
#' The parameters are the category logits `theta[t, c] = log(p[t, c] / p[t, C])`
#' for `c < C`, class by class. For a unit, `d log p_t(y) / d theta[u, c] =
#' 1[t = u] (d_c - p_uc)`.
#' @param fit The list from `.multilpa_step_categorical_fit()`.
#' @param indicators Units-by-categories 0/1 matrix of the outcome.
#' @param class The class `t`.
#' @return A units-by-parameters matrix.
#' @noRd
.multilpa_step_categorical_derivative <- function(fit, indicators, class) {
  n_classes <- nrow(fit$probabilities)
  free <- seq_len(ncol(indicators) - 1L)
  theta <- matrix(0, nrow(indicators), n_classes * length(free))
  theta[, (class - 1L) * length(free) + free] <-
    sweep(indicators[, free, drop = FALSE], 2L, fit$probabilities[class, free], "-")
  theta
}

#' Per-unit scores of the categorical step-three model
#'
#' The responsibility-weighted derivatives, `s_i = sum_t r_it g_it`.
#' @param fit The list from `.multilpa_step_categorical_fit()`.
#' @param indicators Units-by-categories 0/1 matrix of the outcome.
#' @return A units-by-parameters matrix.
#' @noRd
.multilpa_step_categorical_scores <- function(fit, indicators) {
  Reduce(`+`, lapply(seq_len(nrow(fit$probabilities)), function(class) {
    fit$responsibilities[, class] * .multilpa_step_categorical_derivative(fit, indicators, class)
  }))
}

#' The observed information of the categorical step-three model
#'
#' The negative Hessian of the log likelihood, from its closed form. For a
#' mixture `l_i = log sum_t exp(eta_it)`,
#' `d2 l_i = sum_t r_it (d2 eta_it + g_it g_it') - s_i s_i'`; here `d2 eta_it`
#' is the same for every unit: for class `t`'s category logits,
#' `-(diag(p_t) - p_t p_t')`.
#' @param fit The list from `.multilpa_step_categorical_fit()`.
#' @param indicators Units-by-categories 0/1 matrix of the outcome.
#' @return A square matrix over the parameters.
#' @noRd
.multilpa_step_categorical_information <- function(fit, indicators) {
  n_classes <- nrow(fit$probabilities)
  free <- seq_len(ncol(indicators) - 1L)
  spread <- Reduce(`+`, lapply(seq_len(n_classes), function(class) {
    derivative <- .multilpa_step_categorical_derivative(fit, indicators, class)
    crossprod(derivative * fit$responsibilities[, class], derivative)
  }))
  multinomial <- function(p) diag(p, length(p)) - tcrossprod(p)
  curvature <- .latents_block_diagonal(lapply(seq_len(n_classes), function(class) {
    sum(fit$responsibilities[, class]) * multinomial(fit$probabilities[class, free])
  }))
  scores <- .multilpa_step_categorical_scores(fit, indicators)
  crossprod(scores) - spread + curvature
}

#' A block-diagonal matrix from square blocks
#' @param blocks A list of square matrices.
#' @return The block-diagonal matrix.
#' @noRd
.latents_block_diagonal <- function(blocks) {
  sizes <- vapply(blocks, nrow, integer(1))
  starts <- cumsum(sizes) - sizes
  do.call(rbind, lapply(seq_along(blocks), function(block) {
    rows <- matrix(0, sizes[[block]], sum(sizes))
    rows[, starts[[block]] + seq_len(sizes[[block]])] <- blocks[[block]]
    rows
  }))
}

#' Per-unit influence of the categorical step-three model on its probabilities
#'
#' The influence of a unit on the parameters is `I^{-1} s_i`; its category
#' logits are carried to the probabilities of every category by the softmax
#' Jacobian.
#' @param fit The list from `.multilpa_step_categorical_fit()`.
#' @param indicators Units-by-categories 0/1 matrix of the outcome.
#' @return A units-by-(classes x categories) matrix, laid out class by class.
#' @noRd
.multilpa_step_categorical_influence <- function(fit, indicators) {
  information <- .multilpa_step_categorical_information(fit, indicators)
  inverse <- tryCatch(solve(information), error = function(condition) {
    stop(errorCondition(paste(
      "The information of the step-three model is singular; the class",
      "probabilities are not identified."),
      class = "latents_singular_information", call = NULL))
  })
  influence <- .multilpa_step_categorical_scores(fit, indicators) %*% inverse
  n_classes <- nrow(fit$probabilities)
  free <- seq_len(ncol(indicators) - 1L)
  do.call(cbind, lapply(seq_len(n_classes), function(class) {
    p <- fit$probabilities[class, ]
    # d p_c / d theta_d = p_c (1[c = d] - p_d), for every category c and the
    # free logits d.
    jacobian <- diag(p)[, free, drop = FALSE] - tcrossprod(p, p[free])
    columns <- (class - 1L) * length(free) + free
    influence[, columns, drop = FALSE] %*% t(jacobian)
  }))
}

#' Differences in category probabilities between pairs of classes
#'
#' @param probabilities Classes-by-categories estimates.
#' @param clustered The scaled unit contributions, laid out class by class.
#' @param categories The category labels.
#' @param level,method,adjust,vcov_type As in [three_step()].
#' @param quantile The normal quantile for the interval.
#' @param invalid Classes-by-categories logical matrix of the probabilities
#'   reported without inference (BCH shares outside `[0, 1]`); a difference
#'   that involves one keeps its estimate and standard error but gets no
#'   statistic, p-value or interval, and does not count in the adjustment.
#' @return One row per pair of classes and category.
#' @noRd
.multilpa_step_categorical_pairs <- function(probabilities, clustered, categories, level,
                                             method, quantile, adjust, vcov_type,
                                             invalid = array(FALSE, dim(probabilities))) {
  n_classes <- nrow(probabilities)
  if (n_classes < 2L) {
    stop(errorCondition("A single class has no other class to be compared with.",
                        class = "latents_inseparable_classes", call = NULL))
  }
  n_categories <- length(categories)
  pairs <- utils::combn(n_classes, 2L)
  grid <- expand.grid(category = seq_len(n_categories), pair = seq_len(ncol(pairs)))
  class <- pairs[2L, grid$pair]
  reference <- pairs[1L, grid$pair]
  column_of <- function(classes) (classes - 1L) * n_categories + grid$category
  difference <- probabilities[cbind(class, grid$category)] -
    probabilities[cbind(reference, grid$category)]
  errors <- sqrt(colSums((clustered[, column_of(class), drop = FALSE] -
                            clustered[, column_of(reference), drop = FALSE])^2))
  involved <- invalid[cbind(class, grid$category)] | invalid[cbind(reference, grid$category)]
  valid <- is.finite(errors) & errors > 0 & !involved
  statistic <- ifelse(valid, difference / errors, NA_real_)
  raw <- ifelse(valid, 2 * stats::pnorm(-abs(statistic)), NA_real_)
  result <- data.frame(
    level = level, method = method, class = class, reference_class = reference,
    category = categories[grid$category], estimate = difference, standard_error = errors,
    statistic = statistic, p_value = raw, p_adjusted = stats::p.adjust(raw, method = adjust),
    conf_low = ifelse(involved, NA_real_, difference - quantile * errors),
    conf_high = ifelse(involved, NA_real_, difference + quantile * errors),
    row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "adjust") <- adjust
  attr(result, "vcov_type") <- vcov_type
  column <- .latents_column
  percent <- format(100 * (2 * stats::pnorm(quantile) - 1))
  .latents_table(result, NULL,
                 sprintf("Distal outcome: class differences (%s, %s%% CI)", method, percent),
                 list(column("Class", "class", "integer"),
                      column("vs", "reference_class", "integer"),
                      column("Category", "category", "text"),
                      column("Difference", "estimate"),
                      column(sprintf("%s%% CI", percent), c("conf_low", "conf_high"), "ci",
                             "estimate"),
                      column("p", "p_value", "p"),
                      column(sprintf("p (%s)", adjust), "p_adjusted", "p")))
}

#' Wald test that every class has the same outcome distribution
#'
#' The parameters are laid out class by class, `m` per class (the mean, or
#' the probabilities of all but the last category). The contrasts are each
#' class's parameters minus the first class's, and the statistic
#' `d' (L V L')^{-1} d` has `(K - 1) m` degrees of freedom, with `V` the
#' cluster-robust covariance from the summed unit contributions.
#' @param estimates The parameters, class by class.
#' @param clustered The scaled unit contributions, one column per parameter.
#' @param n_classes The number of classes.
#' @param level,method,vcov_type As in [three_step()].
#' @return A one-row table: `level`, `method`, `statistic`, `df`, `p_value`.
#' @noRd
.multilpa_step_overall <- function(estimates, clustered, n_classes, level, method,
                                   vcov_type) {
  if (n_classes < 2L) {
    stop(errorCondition("A single class has no other class to be compared with.",
                        class = "latents_inseparable_classes", call = NULL))
  }
  per_class <- length(estimates) %/% n_classes
  first <- seq_len(per_class)
  contrasts <- do.call(rbind, lapply(seq_len(n_classes)[-1L], function(class) {
    row <- matrix(0, per_class, length(estimates))
    row[cbind(first, first)] <- -1
    row[cbind(first, (class - 1L) * per_class + first)] <- 1
    row
  }))
  difference <- as.vector(contrasts %*% estimates)
  covariance <- contrasts %*% crossprod(clustered) %*% t(contrasts)
  statistic <- if (anyNA(covariance)) NA_real_ else tryCatch(
    as.numeric(crossprod(difference, solve(covariance, difference))),
    error = function(condition) {
      warning(warningCondition(paste(
        "The covariance of the class contrasts is singular, so the overall test",
        "has no statistic."), class = "latents_singular_information", call = NULL))
      NA_real_
    })
  df <- length(difference)
  result <- data.frame(level = level, method = method, statistic = statistic, df = df,
                       p_value = stats::pchisq(statistic, df, lower.tail = FALSE),
                       row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "vcov_type") <- vcov_type
  column <- .latents_column
  .latents_table(result, NULL,
                 sprintf("Distal outcome: equal across classes (%s, Wald)", method),
                 list(column("Wald", "statistic"), column("df", "df", "integer"),
                      column("p", "p_value", "p")))
}
