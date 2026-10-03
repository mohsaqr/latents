# Ordinal outcomes in regression mixtures: mixture_regression(family =
# "ordinal"), phase E7 of validation/ENGINE_DESIGN.md. The GLM block of the
# regression mixtures with a cumulative-logit density, under the same
# structures (rows, groups, two levels) and the same EM, scores and tables.
#
# Within class k the outcome Y in 1..C follows a proportional-odds regression
# (McCullagh 1980):
#
#   P(Y <= c | x, k) = F(t_kc - eta_ik),   F = plogis,  c = 1..C-1,
#   eta_ik = x_i' beta_k + z_i' gamma + offset_i,
#
# with ordered thresholds t_k1 < ... < t_k,C-1 that take the intercept's
# place. Write a = t_k,y - eta and b = t_k,y-1 - eta (t_k0 = -Inf,
# t_kC = Inf) and g = a - b, the gap between the two thresholds. Then
#
#   P(Y = y) = F(a) - F(b) = F(a) (1 - F(b)) (1 - exp(-g)),
#
# so log P(Y = y) is a sum of three terms each computed in its own stable
# form, with no cancellation in either tail and no special case for the
# first or last category (an infinite bound contributes zero).
#
# The joint log likelihood of (thresholds, slopes) is concave (Pratt 1981),
# so the class M-step is Newton's method on the natural thresholds and
# slopes, with steps halved until the posterior-weighted objective does not
# fall and the thresholds stay ordered. For packing and inference the
# thresholds are the unconstrained coordinates of MASS::polr(): the first
# threshold and the logs of the gaps between consecutive thresholds, so
# every coordinate vector is a valid, ordered model.
#
# References: McCullagh, P. (1980). Regression models for ordinal data.
# Journal of the Royal Statistical Society B, 42, 109-142. Pratt, J. W.
# (1981). Concavity of the log likelihood. Journal of the American Statistical
# Association, 76, 103-106.

#' Resolve an ordinal outcome to category codes
#'
#' An ordered factor or a factor is read in its level order; whole numbers
#' are read as categories in increasing order. Categories are those observed
#' (the model frame drops unused factor levels), because a category no row
#' takes has no finite threshold. With `levels`, the outcome of new data is
#' coded against the fitted categories instead.
#'
#' @param response The model response.
#' @param levels `NULL`, or the fitted category labels.
#' @return A list of numeric `y` (codes 1..C), `trials` (ones), `levels`
#'   (labels), `values` (the categories as given: factor levels or numbers)
#'   and `kind` (`"ordered"`, `"factor"` or `"numeric"`); raises
#'   `latents_bad_data`.
#' @noRd
.ordinal_response <- function(response, levels = NULL) {
  bad_data <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  whole <- function(v) all(abs(v - round(v)) < sqrt(.Machine$double.eps))
  if (is.numeric(response) && !is.matrix(response) && whole(response)) {
    response <- round(as.vector(response))
  }
  if (!is.null(levels)) {
    # New data: any type whose values name the fitted categories.
    codes <- match(as.character(response), levels)
    if (anyNA(codes)) {
      bad_data(sprintf(
        "The outcome has categories the model was not fitted with: %s.",
        paste(sprintf("`%s`", unique(as.character(response)[is.na(codes)])),
              collapse = ", ")))
    }
    return(list(y = as.numeric(codes), trials = rep(1, length(codes))))
  }
  usable <- is.factor(response) ||
    (is.numeric(response) && !is.matrix(response) && whole(response))
  if (!usable) {
    bad_data(paste("An ordinal outcome must be an ordered factor, a factor",
                   "(its levels taken in order) or whole-number categories."))
  }
  kind <- if (is.ordered(response)) "ordered" else if (is.factor(response))
    "factor" else "numeric"
  values <- if (is.factor(response)) levels(response) else sort(unique(response))
  labels <- as.character(values)
  if (length(values) < 2L) {
    bad_data("An ordinal outcome needs at least two observed categories.")
  }
  list(y = as.numeric(match(as.character(response), labels)),
       trials = rep(1, length(response)), levels = labels, values = values,
       kind = kind)
}

#' Refuse an ordinal mixture that its data cannot identify
#'
#' With one outcome per class assignment (rows as the classified units), a
#' mixture of two-category outcomes is a mixture of Bernoullis, itself a
#' Bernoulli (Follmann and Lambert 1991), so the same rule as for
#' `family = "binomial"` applies; and with no predictors at all a mixture of
#' categorical distributions is again one categorical distribution, whatever
#' the number of categories. Otherwise the classes are identified by how the
#' predictors shift the category probabilities, which needs continuous
#' predictors; `class_level = "group"` identifies them through each group's
#' repeated outcomes.
#'
#' @param n_categories Number of categories C.
#' @param n_classes Number of classes.
#' @param nesting The nesting.
#' @param n_predictors Number of design columns (class-specific and shared).
#' @return `NULL`, invisibly; raises `latents_not_identified`.
#' @noRd
.ordinal_check_identified <- function(n_categories, n_classes, nesting,
                                      n_predictors) {
  if (n_classes < 2L || identical(nesting, "group")) return(invisible(NULL))
  if (n_categories == 2L) {
    stop(errorCondition(paste(
      "A two-category ordinal outcome is binary, and a mixture of binary",
      "outcomes with one outcome per class assignment is not identified (a",
      "mixture of Bernoulli distributions is itself a Bernoulli; Follmann and",
      "Lambert 1991). With two categories the ordinal model is the logistic",
      "mixture of `family = \"binomial\"`; give each class assignment several",
      "outcomes with `id` and `class_level = \"group\"`."),
      class = "latents_not_identified", call = NULL))
  }
  if (n_predictors == 0L) {
    stop(errorCondition(paste(
      "An ordinal mixture without predictors and with one outcome per class",
      "assignment is not identified: a mixture of categorical distributions",
      "is itself a categorical distribution. Add predictors, or give each",
      "class assignment several outcomes with `id` and",
      "`class_level = \"group\"`."),
      class = "latents_not_identified", call = NULL))
  }
  invisible(NULL)
}

#' Threshold labels, `"lower|upper"` as MASS::polr() names them
#' @param spec The model specification.
#' @return A character vector of length C - 1.
#' @noRd
.ordinal_threshold_labels <- function(spec) {
  levels <- spec$category_levels
  paste(levels[-length(levels)], levels[-1L], sep = "|")
}

#' Starting thresholds: the logits of the pooled cumulative proportions
#'
#' The proportional-odds fit without predictors, the same for every class;
#' each category is observed, so every proportion is inside (0, 1).
#' @param spec The model specification.
#' @return A `(C - 1) x K` matrix.
#' @noRd
.ordinal_start_thresholds <- function(spec) {
  cumulative <- cumsum(tabulate(spec$y, spec$n_categories)) / spec$n
  start <- stats::qlogis(cumulative[-spec$n_categories])
  matrix(start, length(start), spec$n_classes,
         dimnames = list(.ordinal_threshold_labels(spec),
                         paste0("class_", seq_len(spec$n_classes))))
}

#' log(1 - exp(-gap)) for gap > 0, accurate at both ends
#'
#' `-expm1(-gap)` near zero and `log1p(-exp(-gap))` beyond log 2 (Maechler
#' 2012, "Accurately computing log(1 - exp(-|a|))"); an infinite gap gives
#' zero, and a gap that is not positive (thresholds out of order) gives
#' `-Inf`, an impossible model.
#' @param gap Numeric vector or matrix.
#' @return The same shape as `gap`.
#' @noRd
.ordinal_log1mexp <- function(gap) {
  out <- gap
  out[] <- -Inf
  small <- gap > 0 & gap <= log(2)
  large <- gap > log(2)
  out[small] <- log(-expm1(-gap[small]))
  out[large] <- log1p(-exp(-gap[large]))
  out
}

#' The bounds of every row's category under every class
#'
#' @param spec The model specification.
#' @param params The parameter list (`thresholds`).
#' @param eta Linear predictor matrix, `n x K`.
#' @return A list of `n x K` matrices: `upper` (a = t_y - eta), `lower`
#'   (b = t_y-1 - eta) and `gap` (a - b).
#' @noRd
.ordinal_bounds <- function(spec, params, eta) {
  extended <- rbind(-Inf, params$thresholds, Inf)
  upper <- extended[spec$y + 1L, , drop = FALSE]
  lower <- extended[spec$y, , drop = FALSE]
  list(upper = upper - eta, lower = lower - eta, gap = upper - lower)
}

#' Row-by-class log probabilities of the observed categories
#' @param spec The model specification.
#' @param params The parameter list.
#' @param eta Linear predictor matrix, `n x K`.
#' @return An `n x K` matrix.
#' @noRd
.ordinal_log_density <- function(spec, params, eta) {
  bounds <- .ordinal_bounds(spec, params, eta)
  stats::plogis(bounds$upper, log.p = TRUE) +
    stats::plogis(bounds$lower, lower.tail = FALSE, log.p = TRUE) +
    .ordinal_log1mexp(bounds$gap)
}

#' First and second derivatives of the log probability in its two bounds
#'
#' With r = 1 / (exp(g) - 1) and s = exp(g) / (exp(g) - 1)^2:
#' d/da = 1 - F(a) + r, d/db = -F(b) - r, d2/da2 = -F(a)(1 - F(a)) - s,
#' d2/db2 = -F(b)(1 - F(b)) - s, d2/da db = s. An infinite bound gives zero
#' in every term that involves it. The derivative in eta is
#' -(d/da + d/db) = F(a) + F(b) - 1.
#'
#' @inheritParams .ordinal_log_density
#' @return A list of `n x K` matrices `upper`, `lower`, `upper_upper`,
#'   `lower_lower`, `upper_lower`.
#' @noRd
.ordinal_derivatives <- function(spec, params, eta) {
  bounds <- .ordinal_bounds(spec, params, eta)
  f_upper <- stats::plogis(bounds$upper)
  s_upper <- stats::plogis(bounds$upper, lower.tail = FALSE)
  f_lower <- stats::plogis(bounds$lower)
  s_lower <- stats::plogis(bounds$lower, lower.tail = FALSE)
  decay <- exp(-bounds$gap)
  shortfall <- -expm1(-bounds$gap)
  ratio <- decay / shortfall
  curvature <- decay / shortfall^2
  list(upper = s_upper + ratio, lower = -f_lower - ratio,
       upper_upper = -f_upper * s_upper - curvature,
       lower_lower = -f_lower * s_lower - curvature,
       upper_lower = curvature)
}

#' Which threshold bounds each row's category from above and from below
#' @param y Category codes.
#' @param n_categories Number of categories C.
#' @return A list of `n x (C - 1)` indicator matrices `upper` (threshold y)
#'   and `lower` (threshold y - 1).
#' @noRd
.ordinal_indicators <- function(y, n_categories) {
  thresholds <- seq_len(n_categories - 1L)
  list(upper = outer(y, thresholds, `==`) * 1,
       lower = outer(y - 1, thresholds, `==`) * 1)
}

#' Whether every class's thresholds are strictly increasing
#' @noRd
.ordinal_ordered <- function(thresholds) {
  all(is.finite(thresholds)) && all(diff(thresholds) > 0)
}

#' M-step for an ordinal regression mixture
#'
#' Newton-Raphson on the posterior-weighted complete-data log likelihood,
#' `sum_ik tau_ik log P_k(y_i)`, jointly over every class's thresholds and
#' slopes and the slopes shared by every class (`common`), which couple the
#' classes. The objective is concave, so the Newton direction ascends; a step
#' is halved until the objective does not fall and the thresholds stay
#' ordered, so the update can only increase the EM likelihood.
#'
#' @param spec The model specification.
#' @param params The current parameter list.
#' @param tau Row posteriors, `n x K`.
#' @param max_newton Maximum Newton iterations.
#' @return The parameter list with `thresholds`, `beta` and `common` updated
#'   and `separation` set when a threshold or slope diverged.
#' @noRd
.ordinal_update <- function(spec, params, tau, max_newton = 50L) {
  n_classes <- spec$n_classes
  m <- spec$n_categories - 1L
  p <- ncol(spec$x)
  q <- ncol(spec$z)
  own <- seq_len(m + p)
  shared <- m + p + seq_len(q)
  indicators <- .ordinal_indicators(spec$y, spec$n_categories)
  # A row's log probability depends on its two bounds a = t_y - eta and
  # b = t_y-1 - eta; these are their derivatives in (thresholds, beta, gamma).
  upper_design <- cbind(indicators$upper, -spec$x, -spec$z)
  lower_design <- cbind(indicators$lower, -spec$x, -spec$z)
  pack <- function(pars) c(as.vector(rbind(pars$thresholds, pars$beta)), pars$common)
  unpack <- function(theta, pars) {
    block <- matrix(theta[seq_len((m + p) * n_classes)], m + p, n_classes)
    pars$thresholds[] <- block[seq_len(m), , drop = FALSE]
    pars$beta[] <- block[m + seq_len(p), , drop = FALSE]
    if (q > 0L) pars$common[] <- theta[(m + p) * n_classes + seq_len(q)]
    pars
  }
  objective <- function(pars) .mixture_regression_objective(spec, pars, tau)
  newton_step <- function(pars) {
    eta <- .mixture_linear_predictors(spec, pars)
    d <- .ordinal_derivatives(spec, pars, eta)
    blocks <- lapply(seq_len(n_classes), function(k) {
      w <- tau[, k]
      cross <- crossprod(upper_design, lower_design * (w * d$upper_lower[, k]))
      list(gradient = drop(crossprod(upper_design, w * d$upper[, k]) +
                             crossprod(lower_design, w * d$lower[, k])),
           information = -(crossprod(upper_design, upper_design * (w * d$upper_upper[, k])) +
                             crossprod(lower_design, lower_design * (w * d$lower_lower[, k])) +
                             cross + t(cross)))
    })
    gradient <- c(unlist(lapply(blocks, function(b) b$gradient[own])),
                  if (q > 0L) Reduce(`+`, lapply(blocks, function(b) b$gradient[shared])))
    information <- .mixture_block_diagonal(lapply(blocks, function(b) {
      b$information[own, own, drop = FALSE]
    }))
    if (q > 0L) {
      cross <- do.call(rbind, lapply(blocks, function(b) {
        b$information[own, shared, drop = FALSE]
      }))
      common_block <- Reduce(`+`, lapply(blocks, function(b) {
        b$information[shared, shared, drop = FALSE]
      }))
      information <- rbind(cbind(information, cross), cbind(t(cross), common_block))
    }
    # As for the GLM families: a ridge far below the data scale keeps an
    # empty class solvable without moving a well-determined solution.
    ridge <- 1e-10 * (1 + max(abs(diag(information))))
    .mixture_solve(information + diag(ridge, nrow(information)), gradient)
  }
  current <- objective(params)
  iteration <- 0L
  # Newton iterations are inherently sequential; each needs the previous point.
  repeat {
    iteration <- iteration + 1L
    direction <- newton_step(params)
    theta <- pack(params)
    step <- 1
    accepted <- FALSE
    while (step > 1e-8 && !accepted) {
      candidate <- unpack(theta + step * direction, params)
      value <- if (.ordinal_ordered(candidate$thresholds)) objective(candidate) else -Inf
      accepted <- is.finite(value) && value >= current - 1e-12 * (1 + abs(current))
      if (!accepted) step <- step / 2
    }
    if (!accepted) break
    gain <- value - current
    params <- candidate
    current <- value
    if (gain <= 1e-12 * (1 + abs(current)) || iteration >= max_newton) break
  }
  # A threshold beyond 30 on the logit scale is a category the class almost
  # never (or almost always) reaches: separation, as for a GLM coefficient.
  params$separation <- any(abs(params$beta) > 30) ||
    (q > 0L && any(abs(params$common) > 30)) || any(abs(params$thresholds) > 30)
  params
}

#' Jacobian of one class's natural thresholds in their packed coordinates
#'
#' The coordinates are the first threshold and the logs of the gaps, so
#' d t_j / d first = 1 and d t_j / d log gap_l = gap_l for l <= j.
#' @param thresholds One class's ordered thresholds.
#' @return A `(C - 1) x (C - 1)` matrix, rows thresholds, columns coordinates.
#' @noRd
.ordinal_threshold_jacobian <- function(thresholds) {
  m <- length(thresholds)
  scale <- c(1, diff(thresholds))
  (row(diag(m)) >= col(diag(m))) * matrix(scale, m, m, byrow = TRUE)
}

#' Pack every class's thresholds into named unconstrained coordinates
#' @param spec The model specification.
#' @param params The parameter list.
#' @return A named numeric vector: per class, `threshold.class_k.<a|b>` and
#'   `log_threshold_gap.class_k.<b|c>` for the gaps.
#' @noRd
.ordinal_pack <- function(spec, params) {
  m <- spec$n_categories - 1L
  class_names <- paste0("class_", seq_len(spec$n_classes))
  coordinates <- rbind(params$thresholds[1L, , drop = FALSE],
                       log(diff(params$thresholds)))
  prefix <- c("threshold.", rep("log_threshold_gap.", m - 1L))
  stats::setNames(as.vector(coordinates),
                  paste0(rep(prefix, spec$n_classes), rep(class_names, each = m), ".",
                         rep(.ordinal_threshold_labels(spec), spec$n_classes)))
}

#' Natural thresholds from packed coordinates, one class per column
#' @param values Coordinates, `(C - 1) * K` of them in packing order.
#' @param m Number of thresholds per class.
#' @return An `m x K` matrix of ordered thresholds.
#' @noRd
.ordinal_natural <- function(values, m) {
  coordinates <- matrix(values, m)
  increments <- rbind(coordinates[1L, , drop = FALSE],
                      exp(coordinates[-1L, , drop = FALSE]))
  # A lower-triangular matrix of ones accumulates each column.
  (row(diag(m)) >= col(diag(m))) %*% increments
}

#' Per-row scores of the threshold coordinates
#'
#' Fisher's identity: each row's score is its posterior-weighted
#' complete-data score, here in the natural thresholds, carried to the packed
#' coordinates by the Jacobian.
#' @param spec The model specification.
#' @param params The parameter list.
#' @param tau Row posteriors.
#' @param derivatives From `.ordinal_derivatives()`.
#' @return An `n x ((C - 1) K)` matrix in packing order.
#' @noRd
.ordinal_threshold_scores <- function(spec, params, tau, derivatives) {
  indicators <- .ordinal_indicators(spec$y, spec$n_categories)
  do.call(cbind, lapply(seq_len(spec$n_classes), function(k) {
    natural <- indicators$upper * (tau[, k] * derivatives$upper[, k]) +
      indicators$lower * (tau[, k] * derivatives$lower[, k])
    natural %*% .ordinal_threshold_jacobian(params$thresholds[, k])
  }))
}

#' Category probabilities of one class
#' @param eta Linear predictors (a vector).
#' @param thresholds That class's ordered thresholds.
#' @return A `length(eta) x C` matrix whose rows sum to one.
#' @noRd
.ordinal_probabilities <- function(eta, thresholds) {
  extended <- c(-Inf, thresholds, Inf)
  gaps <- diff(extended)
  n_categories <- length(gaps)
  matrix(vapply(seq_len(n_categories), function(category) {
    exp(stats::plogis(extended[category + 1L] - eta, log.p = TRUE) +
          stats::plogis(extended[category] - eta, lower.tail = FALSE, log.p = TRUE) +
          .ordinal_log1mexp(gaps[category]))
  }, numeric(length(eta))), length(eta), n_categories)
}

#' Expected category score of every row under every class
#'
#' The categories are scored 1..C in their order, so the class "mean" of an
#' ordinal outcome is sum_c c P_k(Y = c | x), on the scale of the category
#' codes.
#' @param eta Linear predictor matrix, `n x K`.
#' @param thresholds A `(C - 1) x K` threshold matrix.
#' @return An `n x K` matrix.
#' @noRd
.ordinal_expected_scores <- function(eta, thresholds) {
  n_categories <- nrow(thresholds) + 1L
  matrix(vapply(seq_len(ncol(thresholds)), function(k) {
    drop(.ordinal_probabilities(eta[, k], thresholds[, k]) %*% seq_len(n_categories))
  }, numeric(nrow(eta))), nrow(eta), ncol(thresholds))
}

#' Category values from codes, in the type the outcome was given
#'
#' A factor outcome comes back as a factor with the fitted levels (ordered
#' if it was), numbers as the fitted category values, so a simulated outcome
#' can replace the observed column and be refitted.
#' @param spec The model specification.
#' @param codes Category codes.
#' @return A factor or numeric vector.
#' @noRd
.ordinal_categories <- function(spec, codes) {
  if (identical(spec$category_kind, "numeric")) return(spec$category_values[codes])
  factor(spec$category_levels[codes], levels = spec$category_levels,
         ordered = identical(spec$category_kind, "ordered"))
}

#' Draw one outcome per row from its class's category probabilities
#' @param spec The model specification.
#' @param params The parameter list.
#' @param classes Each row's class.
#' @return Category values, as `.ordinal_categories()`.
#' @noRd
.ordinal_draw <- function(spec, params, classes) {
  eta <- .mixture_linear_predictors(spec, params)
  probabilities <- Reduce(`+`, lapply(seq_len(spec$n_classes), function(k) {
    .ordinal_probabilities(eta[, k], params$thresholds[, k]) * (classes == k)
  }))
  .ordinal_categories(spec, .latents_draw_rows(probabilities))
}

#' Category probabilities of every row under every class, as a tidy table
#' @param spec The (possibly new-data) specification.
#' @param params The parameter list.
#' @param priors Row-by-class prior probabilities.
#' @return A data frame, one row per data row, class and category: `row`,
#'   `class`, `category`, `prior`, `probability`.
#' @noRd
.ordinal_probability_table <- function(spec, params, priors) {
  eta <- .mixture_linear_predictors(spec, params)
  n_categories <- spec$n_categories
  probabilities <- vapply(seq_len(spec$n_classes), function(k) {
    .ordinal_probabilities(eta[, k], params$thresholds[, k])
  }, matrix(0, spec$n, n_categories))
  # probabilities is n x C x K; the table runs category fastest, then class.
  arranged <- aperm(array(probabilities, c(spec$n, n_categories, spec$n_classes)),
                    c(2L, 3L, 1L))
  data.frame(row = rep(spec$kept_rows, each = n_categories * spec$n_classes),
             class = rep(rep(paste0("class_", seq_len(spec$n_classes)),
                             each = n_categories), spec$n),
             category = rep(spec$category_levels, spec$n * spec$n_classes),
             prior = rep(as.vector(t(priors)), each = n_categories),
             probability = as.vector(arranged))
}

#' Threshold rows of the coefficient table
#'
#' The natural thresholds with delta-method standard errors from the packed
#' coordinates (an exact Jacobian).
#' @param spec The model specification.
#' @param inference The inference list (`theta`, `vcov`).
#' @param level Confidence level.
#' @return A data frame with `class`, `term` (`"threshold:<a|b>"`) and the
#'   Wald columns.
#' @noRd
.ordinal_threshold_table <- function(spec, inference, level) {
  m <- spec$n_categories - 1L
  labels <- .ordinal_threshold_labels(spec)
  class_names <- paste0("class_", seq_len(spec$n_classes))
  do.call(rbind, lapply(class_names, function(class) {
    coordinates <- paste0(c("threshold.", rep("log_threshold_gap.", m - 1L)), class,
                          ".", labels)
    natural <- drop(.ordinal_natural(inference$theta[coordinates], m))
    jacobian <- .ordinal_threshold_jacobian(natural)
    std_error <- .inference_delta_se(
      jacobian, inference$vcov[coordinates, coordinates, drop = FALSE], "rowsums")
    cbind(data.frame(class = class, term = paste0("threshold:", labels)),
          .inference_wald(natural, std_error, level))
  }))
}

#' Each ordinal class's expected-score trajectory with a delta-method band
#'
#' The expected category score depends on the class's thresholds as well as
#' its linear predictor, so its standard error comes from the delta method
#' over every packed estimate, and the band is formed on the score scale,
#' clipped to the range of the categories.
#' @param x A trajectory fit with `family = "ordinal"`.
#' @param design The regression design of the rows (class-specific columns,
#'   then shared ones).
#' @param inference The inference list.
#' @param z The critical value.
#' @param classes The class names.
#' @return As `.trajectory_band()`.
#' @noRd
.ordinal_trajectory_band <- function(x, design, inference, z, classes) {
  spec <- x$spec
  do.call(rbind, lapply(seq_along(classes), function(k) {
    delta <- .mixture_delta(inference$theta, inference$vcov, function(theta) {
      params <- .mixture_unpack(spec, theta, x$params)
      eta <- design %*% c(params$beta[, k], params$common)
      drop(.ordinal_expected_scores(eta, params$thresholds[, k, drop = FALSE]))
    })
    data.frame(class = classes[k], estimate = delta$estimate,
               std_error = delta$std_error,
               conf_low = pmax(delta$estimate - z * delta$std_error, 1),
               conf_high = pmin(delta$estimate + z * delta$std_error, spec$n_categories),
               stringsAsFactors = FALSE)
  }))
}
