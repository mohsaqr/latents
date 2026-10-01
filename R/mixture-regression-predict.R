# predict(), simulate() and plot() for mixture_regression() fits.

#' Rebuild a specification on new data with the fitted design
#'
#' Uses the fit's terms and factor levels, so a new data frame yields columns
#' in the fitted order. The outcome is resolved only when `with_response`.
#'
#' @param object A `latents_mixture_regression` fit.
#' @param newdata A data frame.
#' @param with_response Whether the outcome must be present.
#' @return A specification list usable by the engine.
#' @noRd
.mixture_new_spec <- function(object, newdata, with_response) {
  spec <- object$spec
  stopifnot("`newdata` must be a data frame" = is.data.frame(newdata))
  model_terms <- if (with_response) spec$terms else
    stats::delete.response(spec$terms)
  frame <- stats::model.frame(model_terms, newdata, xlev = spec$xlevels,
                              na.action = stats::na.fail)
  full_design <- stats::model.matrix(model_terms, frame,
                                     contrasts.arg = spec$contrasts)
  new <- spec
  new$n <- nrow(newdata)
  new$x <- full_design[, colnames(spec$x), drop = FALSE]
  new$z <- full_design[, colnames(spec$z), drop = FALSE]
  new$offset <- as.numeric(stats::model.offset(frame) %||% rep(0, new$n))
  new$kept_rows <- seq_len(new$n)
  if (with_response) {
    response <- .mixture_response(stats::model.response(frame), spec$family)
    new$y <- response$y
    new$trials <- response$trials
    new$log_normalizer <- switch(spec$family,
      gaussian = rep(0, new$n),
      binomial = lchoose(response$trials, response$y),
      poisson = -lgamma(response$y + 1))
  } else {
    new$y <- rep(0, new$n)
    new$trials <- rep(1, new$n)
    new$log_normalizer <- rep(0, new$n)
  }
  membership_matrix <- function(name) {
    blueprint <- spec[[paste0(name, "_design")]] %||%
      .mixture_membership_design(spec[[name]], spec$model_data)
    frame <- stats::model.frame(blueprint$terms, newdata,
                                xlev = blueprint$xlevels,
                                na.action = stats::na.fail)
    stats::model.matrix(blueprint$terms, frame,
                        contrasts.arg = blueprint$contrasts)
  }
  row_membership <- membership_matrix("membership")
  if (!is.null(spec$id) && (!identical(spec$nesting, "observation") ||
                            with_response || spec$id %in% names(newdata))) {
    if (!spec$id %in% names(newdata)) {
      stop(errorCondition(sprintf(
        "`newdata` needs the grouping column `%s` for this model.", spec$id),
        class = "latents_bad_data", call = NULL))
    }
    if (anyNA(newdata[[spec$id]])) {
      stop(errorCondition("The grouping column in `newdata` has missing values.",
                          class = "latents_bad_data", call = NULL))
    }
    new$group_levels <- sort(unique(newdata[[spec$id]]))
    new$group_index <- match(newdata[[spec$id]], new$group_levels)
    new$n_groups <- length(new$group_levels)
  }
  new$w <- switch(spec$nesting,
    observation = row_membership,
    group = .mixture_group_design(row_membership, new$group_index,
                                 "membership"),
    "two-level" = row_membership[, colnames(spec$w), drop = FALSE])
  if (identical(spec$nesting, "two-level")) {
    new$v <- .mixture_group_design(
      membership_matrix("group_membership"),
      new$group_index, "group_membership")
  }
  new
}

#' Prior class probabilities of every row, before seeing its outcome
#' @noRd
.mixture_row_priors <- function(spec, params) {
  switch(spec$nesting,
    observation = exp(.mixture_log_softmax(spec$w, params$gamma)),
    group = exp(.mixture_log_softmax(spec$w, params$gamma))[
      spec$group_index, , drop = FALSE],
    "two-level" = {
      eta <- exp(.mixture_log_softmax(spec$v, params$delta))
      Reduce(`+`, lapply(seq_len(spec$n_group_classes), function(h) {
        eta[spec$group_index, h] * exp(.mixture_log_softmax(
          .mixture_two_level_design(spec, h), params$class_logits))
      }))
    })
}

#' Predict from a mixture-of-regressions fit
#'
#' @param object A fit from [mixture_regression()].
#' @param newdata `NULL` for the fitted rows, or a data frame with the
#'   predictors (and the membership covariates and `id` the model uses).
#' @param type `"response"`: the mean outcome averaged over the classes with
#'   their prior probabilities (the covariates, but not the outcome, inform
#'   them). `"class_response"`: the mean outcome under every class, one row
#'   per row and class. `"posterior"`: the posterior class probabilities
#'   given the outcome, which `newdata` must then contain.
#' @param ... Unused.
#' @return A base `data.frame`. For `"response"`: `row` and `fitted`. For
#'   `"class_response"`: `row`, `class`, `prior` and `fitted`, one row per
#'   data row and class. For `"posterior"`: the same columns as
#'   `get_results(fit, "assignments")`.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' new_students <- data.frame(hours = c(2, 10))
#' predict(fit, new_students, type = "class_response")
#' @importFrom stats predict
#' @export
predict.latents_mixture_regression <- function(object, newdata = NULL,
                                   type = c("response", "class_response",
                                            "posterior"), ...) {
  type <- match.arg(type)
  spec <- if (is.null(newdata)) object$spec else
    .mixture_new_spec(object, newdata, identical(type, "posterior"))
  params <- object$params
  if (identical(type, "posterior")) {
    refit <- object
    refit$spec <- spec
    refit$expectation <- .mixture_expectation(spec, params, weighted = FALSE)
    return(.mixture_assignment_table(refit))
  }
  means <- .mixture_class_means(spec, params)
  priors <- .mixture_row_priors(spec, params)
  if (identical(type, "response")) {
    return(data.frame(row = spec$kept_rows, fitted = rowSums(priors * means)))
  }
  data.frame(row = rep(spec$kept_rows, spec$n_classes),
             class = rep(paste0("class_", seq_len(spec$n_classes)),
                         each = spec$n),
             prior = as.vector(priors), fitted = as.vector(means))
}

#' Simulate outcomes from a mixture-of-regressions fit
#'
#' Draws class memberships from the fitted mixing model (group classes first
#' in the two-level model, one class per group under `class_level =
#' "group"`), then outcomes from each row's class regression, at the fitted
#' rows' predictors.
#'
#' @param object A fit from [mixture_regression()].
#' @param nsim Number of simulated outcome vectors.
#' @param seed `NULL` or an integer; the caller's random state is restored.
#' @param ... Unused.
#' @return A base `data.frame` with one row per fitted row and one column per
#'   simulation, `sim_1`, `sim_2`, ...; for a binomial fit with more than one
#'   trial the columns hold success counts.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' head(simulate(fit, nsim = 2, seed = 1))
#' @importFrom stats simulate
#' @export
simulate.latents_mixture_regression <- function(object, nsim = 1, seed = NULL, ...) {
  stopifnot("`nsim` must be a single positive integer" = .mixture_is_count(nsim))
  .mixture_with_seed(seed, {
    draws <- lapply(seq_len(nsim), function(i) .mixture_draw(object))
    stats::setNames(as.data.frame(draws), paste0("sim_", seq_len(nsim)))
  })
}

#' One simulated outcome vector
#' @noRd
.mixture_draw <- function(object) {
  spec <- object$spec
  params <- object$params
  draw_category <- function(probabilities) {
    cumulative <- t(apply(probabilities, 1L, cumsum))
    cumulative <- matrix(cumulative, nrow(probabilities))
    1L + rowSums(stats::runif(nrow(probabilities)) > cumulative)
  }
  classes <- switch(spec$nesting,
    observation = draw_category(exp(.mixture_log_softmax(spec$w,
                                                        params$gamma))),
    group = draw_category(exp(.mixture_log_softmax(
      spec$w, params$gamma)))[spec$group_index],
    "two-level" = {
      group_class <- draw_category(exp(.mixture_log_softmax(
        spec$v, params$delta)))[spec$group_index]
      by_group_class <- lapply(seq_len(spec$n_group_classes), function(h) {
        exp(.mixture_log_softmax(.mixture_two_level_design(spec, h),
                                params$class_logits))
      })
      probabilities <- Reduce(`+`, lapply(seq_len(spec$n_group_classes),
                                          function(h) {
        by_group_class[[h]] * (group_class == h)
      }))
      draw_category(probabilities)
    })
  classes <- pmin(classes, spec$n_classes)
  eta <- .mixture_linear_predictors(spec, params)[cbind(seq_len(spec$n),
                                                      classes)]
  switch(spec$family,
    gaussian = stats::rnorm(spec$n, eta, sqrt(params$sigma2[classes])),
    binomial = stats::rbinom(spec$n, spec$trials, stats::plogis(eta)),
    poisson = stats::rpois(spec$n, exp(eta)))
}

#' Plot a mixture-of-regressions fit
#'
#' @param x A fit from [mixture_regression()].
#' @param what `"fitted"`: the outcome against one numeric predictor, rows
#'   marked by modal class (colour and symbol), each class's regression line
#'   drawn with the other predictors at their mean (or reference level).
#'   `"coefficients"`: every class's estimates with confidence intervals.
#'   `"posteriors"`: the distribution of each unit's largest posterior
#'   probability, by class.
#' @param predictor For `"fitted"`, the numeric predictor for the horizontal
#'   axis; defaults to the first one in the formula.
#' @param level Confidence level for `"coefficients"`.
#' @param main Plot title; `NULL` for the default.
#' @param ... Unused.
#' @return `x`, invisibly. Called for its plot.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' plot(fit)
#' plot(fit, what = "coefficients")
#' @export
plot.latents_mixture_regression <- function(x, what = c("fitted", "coefficients",
                                            "posteriors"),
                                predictor = NULL, level = 0.95, main = NULL,
                                ...) {
  what <- match.arg(what)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  switch(what,
    fitted = .mixture_plot_fitted(x, predictor, main),
    coefficients = .mixture_plot_coefficients(x, level, main),
    posteriors = .mixture_plot_posteriors(x, main))
  invisible(x)
}

#' @noRd
.mixture_plot_fitted <- function(x, predictor, main) {
  spec <- x$spec
  frame <- .mixture_model_data(x)
  predictors <- all.vars(stats::delete.response(spec$terms))
  numeric_columns <- predictors[vapply(frame[predictors], is.numeric, logical(1))]
  predictor <- predictor %||% numeric_columns[1L]
  if (is.na(predictor) || !predictor %in% numeric_columns) {
    stop(errorCondition(paste(
      "`what = \"fitted\"` needs a numeric predictor in the formula; name",
      "one with `predictor`."), class = "latents_bad_argument", call = NULL))
  }
  k <- spec$n_classes
  colours <- .multilpa_palette(k)
  symbols <- .multilpa_symbols(k)
  lines_types <- .multilpa_linetypes(k)
  class <- max.col(x$expectation$tau, ties.method = "first")
  observed <- if (identical(spec$family, "binomial") && any(spec$trials > 1)) {
    spec$y / spec$trials
  } else spec$y
  graphics::plot(frame[[predictor]], observed, pch = symbols[class],
                 col = colours[class], bg = grDevices::adjustcolor(
                   colours[class], alpha.f = 0.4),
                 xlab = predictor, ylab = spec$response_name,
                 main = main %||% sprintf("%d-class mixture regression", k))
  grid_values <- seq(min(frame[[predictor]]), max(frame[[predictor]]),
                     length.out = 100L)
  typical <- .mixture_typical_rows(frame, predictor, grid_values)
  new_spec <- .mixture_new_spec(x, typical, with_response = FALSE)
  means <- .mixture_class_means(new_spec, x$params)
  if (identical(spec$family, "binomial")) means <- stats::plogis(
    .mixture_linear_predictors(new_spec, x$params))
  invisible(lapply(seq_len(k), function(j) {
    graphics::lines(grid_values, means[, j], col = colours[j],
                    lty = lines_types[j], lwd = 2)
  }))
  graphics::legend("topleft", legend = paste("class", seq_len(k)),
                   col = colours, pch = symbols, pt.bg = colours,
                   lty = lines_types, lwd = 2, bty = "n")
}

#' The data columns the model used, in fitted row order
#' @noRd
.mixture_model_data <- function(x) {
  x$spec$model_data %||% stop(errorCondition(
    "This fit does not carry its data; refit it with this version of latents.",
    class = "latents_bad_argument", call = NULL))
}

#' Rows varying one predictor, others held at mean or reference level
#' @noRd
.mixture_typical_rows <- function(predictors, predictor, grid_values) {
  typical <- lapply(predictors, function(column) {
    if (is.numeric(column)) mean(column) else if (is.factor(column)) {
      factor(levels(column)[1L], levels = levels(column))
    } else column[1L]
  })
  rows <- as.data.frame(lapply(typical, rep, length(grid_values)),
                        stringsAsFactors = FALSE)
  names(rows) <- names(predictors)
  rows[[predictor]] <- grid_values
  rows
}

#' @noRd
.mixture_plot_coefficients <- function(x, level, main) {
  table <- .mixture_table(x, "coefficients", level)
  k <- x$spec$n_classes
  labels <- unique(table$term)
  classes <- unique(table$class)
  colours <- .multilpa_palette(length(classes))
  symbols <- .multilpa_symbols(length(classes))
  term_position <- match(table$term, labels)
  class_position <- match(table$class, classes)
  offset <- (class_position - (length(classes) + 1) / 2) * 0.15
  y <- term_position + offset
  graphics::par(mar = c(4, 9, 3, 1))
  range_x <- range(c(table$conf_low, table$conf_high), na.rm = TRUE)
  graphics::plot(table$estimate, y, xlim = range_x, yaxt = "n",
                 pch = symbols[class_position],
                 col = colours[class_position], bg = colours[class_position],
                 xlab = sprintf("Estimate with %g%% interval", 100 * level),
                 ylab = "", main = main %||% "Class regression coefficients")
  graphics::segments(table$conf_low, y, table$conf_high, y,
                     col = colours[class_position])
  graphics::abline(v = 0, lty = 3, col = "#999999")
  graphics::axis(2, at = seq_along(labels), labels = labels, las = 1)
  graphics::legend("topright", legend = classes, col = colours, pch = symbols,
                   pt.bg = colours, bty = "n")
}

#' @noRd
.mixture_plot_posteriors <- function(x, main) {
  spec <- x$spec
  posterior <- if (identical(spec$nesting, "group")) x$expectation$group_tau else
    x$expectation$tau
  class <- max.col(posterior, ties.method = "first")
  largest <- posterior[cbind(seq_len(nrow(posterior)), class)]
  k <- spec$n_classes
  colours <- .multilpa_palette(k)
  symbols <- .multilpa_symbols(k)
  graphics::stripchart(split(largest, factor(class, levels = seq_len(k))),
                       method = "jitter", vertical = TRUE, pch = symbols,
                       col = colours, bg = colours,
                       group.names = paste("class", seq_len(k)),
                       ylab = "Largest posterior probability", ylim = c(0, 1),
                       main = main %||% "Classification certainty")
  graphics::abline(h = 1 / k, lty = 3, col = "#999999")
}
