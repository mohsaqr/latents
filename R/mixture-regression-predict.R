utils::globalVariables(c("predictor_value", "outcome", "class_label",
                         "position", "certainty"))

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
  if (any(!is.finite(full_design)) || any(!is.finite(new$offset))) {
    stop(errorCondition("The evaluated regression design and offsets must be finite.",
                        class = "latents_bad_data", call = NULL))
  }
  new$kept_rows <- seq_len(new$n)
  if (with_response) {
    response <- .mixture_response(stats::model.response(frame), spec$family,
                                  spec$category_levels)
    new$y <- response$y
    new$trials <- response$trials
    new$log_normalizer <- switch(spec$family,
      gaussian = rep(0, new$n),
      binomial = lchoose(response$trials, response$y),
      poisson = ,
      negative_binomial = -lgamma(response$y + 1),
      ordinal = rep(0, new$n))
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
    design <- stats::model.matrix(blueprint$terms, frame,
                                   contrasts.arg = blueprint$contrasts)
    if (any(!is.finite(design))) {
      stop(errorCondition("The evaluated membership design must be finite.",
                          class = "latents_bad_data", call = NULL))
    }
    design
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

#' @rdname latents-predict
#' @importFrom stats predict
#' @export
predict.latents_mixture_regression <- function(object, newdata = NULL,
                                   type = c("response", "class_response",
                                            "posterior", "probabilities"), ...,
                                   trials = NULL) {
  type <- match.arg(type)
  if (identical(type, "probabilities") && !identical(object$spec$family, "ordinal")) {
    stop(errorCondition(paste(
      "`type = \"probabilities\"` gives category probabilities of an ordinal",
      "outcome; use \"response\" or \"class_response\" for this family."),
      class = "latents_bad_argument", call = NULL))
  }
  spec <- if (is.null(newdata)) object$spec else
    .mixture_prediction_data(.mixture_new_spec(
      object, newdata, identical(type, "posterior")))
  if (identical(spec$family, "binomial")) {
    resolved <- .mixture_prediction_trials(object, newdata, spec, trials,
                                            identical(type, "posterior"))
    spec$trials <- resolved
  } else if (!is.null(trials)) {
    stop(errorCondition("`trials` applies only to binomial predictions.",
                        class = "latents_bad_argument", call = NULL))
  }
  params <- object$params
  if (identical(type, "posterior")) {
    refit <- object
    refit$spec <- spec
    refit$expectation <- .mixture_expectation(spec, params, weighted = FALSE)
    return(.mixture_assignment_table(refit))
  }
  priors <- .mixture_row_priors(spec, params)
  if (identical(type, "probabilities")) {
    return(.ordinal_probability_table(spec, params, priors))
  }
  means <- .mixture_class_means(spec, params)
  if (identical(type, "response")) {
    return(data.frame(row = spec$kept_rows, fitted = rowSums(priors * means)))
  }
  data.frame(row = rep(spec$kept_rows, spec$n_classes),
             class = rep(paste0("class_", seq_len(spec$n_classes)),
                         each = spec$n),
             prior = as.vector(priors), fitted = as.vector(means))
}

#' Translate evaluated prediction-data failures into the public condition
#' @param expression A prediction specification or response expression.
#' @return Its value, or an error of class `latents_bad_data`.
#' @noRd
.mixture_prediction_data <- function(expression) {
  tryCatch(expression, error = function(condition) {
    stop(errorCondition(conditionMessage(condition), class = "latents_bad_data",
                        call = NULL))
  })
}

#' Resolve binomial exposure without changing the units of a prediction
#'
#' The training response gives the fitted trial counts. New counts can come
#' from an explicit vector/column, or the evaluated response when supplied.
#' Without either, an aggregated fit cannot assume one trial per new row.
#' @param object A fitted mixture regression.
#' @param newdata `NULL` or the new data frame.
#' @param spec The specification prepared for prediction.
#' @param trials `NULL`, a numeric vector, or a column name.
#' @param posterior Whether the outcome is being used to classify rows.
#' @return A numeric trial-count vector of length `spec$n`.
#' @noRd
.mixture_prediction_trials <- function(object, newdata, spec, trials, posterior) {
  bad_data <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  if (!is.null(trials)) {
    if (is.character(trials) && length(trials) == 1L && !is.na(trials)) {
      frame <- newdata %||% object$spec$model_data
      if (!trials %in% names(frame)) bad_data("`trials` names no column of the prediction data.")
      trials <- frame[[trials]]
    }
    if (!is.numeric(trials) || !is.null(dim(trials)) ||
        !length(trials) %in% c(1L, spec$n) || any(!is.finite(trials)) ||
        any(trials < 0) ||
        any(abs(trials - round(trials)) > sqrt(.Machine$double.eps))) {
      bad_data("`trials` must be finite, non-negative whole numbers: one value or one per prediction row.")
    }
    if (length(trials) == 1L) trials <- rep(trials, spec$n)
    trials <- as.numeric(trials)
    if (posterior && !isTRUE(all.equal(trials, as.numeric(spec$trials)))) {
      bad_data("For posterior prediction, `trials` must agree with the observed binomial response.")
    }
    return(trials)
  }
  if (is.null(newdata) || posterior) return(spec$trials)
  response <- attr(object$spec$terms, "variables")[[2L]]
  counts <- any(object$spec$trials != 1) || is.matrix(.mixture_prediction_data(
    eval(response, object$spec$model_data, environment(object$spec$terms))))
  if (!counts) return(rep(1, spec$n))
  if (all(all.vars(response) %in% names(newdata))) {
    return(.mixture_prediction_data(.mixture_response(
      eval(response, newdata, environment(object$spec$terms)), "binomial"))$trials)
  }
  bad_data(paste("Binomial count predictions need `trials` or the response",
                 "columns in `newdata`; trial counts cannot be inferred from predictors alone."))
}

#' @rdname latents-predict
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
  classes <- switch(spec$nesting,
    observation = .latents_draw_rows(exp(.mixture_log_softmax(spec$w,
                                                              params$gamma))),
    group = .latents_draw_rows(exp(.mixture_log_softmax(
      spec$w, params$gamma)))[spec$group_index],
    "two-level" = {
      group_class <- .latents_draw_rows(exp(.mixture_log_softmax(
        spec$v, params$delta)))[spec$group_index]
      by_group_class <- lapply(seq_len(spec$n_group_classes), function(h) {
        exp(.mixture_log_softmax(.mixture_two_level_design(spec, h),
                                params$class_logits))
      })
      probabilities <- Reduce(`+`, lapply(seq_len(spec$n_group_classes),
                                          function(h) {
        by_group_class[[h]] * (group_class == h)
      }))
      .latents_draw_rows(probabilities)
    })
  classes <- pmin(classes, spec$n_classes)
  if (identical(spec$family, "ordinal")) return(.ordinal_draw(spec, params, classes))
  eta <- .mixture_linear_predictors(spec, params)[cbind(seq_len(spec$n),
                                                      classes)]
  switch(spec$family,
    gaussian = stats::rnorm(spec$n, eta, sqrt(params$sigma2[classes])),
    binomial = stats::rbinom(spec$n, spec$trials, stats::plogis(eta)),
    poisson = stats::rpois(spec$n, exp(eta)),
    negative_binomial = stats::rnbinom(spec$n, size = 1 / params$sigma2[classes],
                                       mu = exp(eta)))
}

#' @rdname latents-plots
#' @export
plot.latents_mixture_regression <- function(x, what = NULL, predictor = NULL,
                                            level = 0.95, main = NULL, subtitle = NULL,
                                            time = NULL, facet = FALSE, persons = 12L,
                                            max_persons = 120L, ...) {
  stopifnot("`x` must be a `latents_mixture_regression` fit" =
              inherits(x, "latents_mixture_regression"))
  .latents_warn_unused_arguments(
    list(...), "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  trajectory <- identical(x$spec$nesting, "group")
  trajectory_views <- .mixture_trajectory_views()
  # Naming a `predictor` asks for the fitted-values view, as before.
  what <- match.arg(what %||% if (trajectory && is.null(predictor)) "trajectories" else
    "fitted",
                    union(.mixture_plot_views(x), trajectory_views))
  if (what %in% trajectory_views && !trajectory) {
    stop(errorCondition(sprintf(paste(
      "`what = \"%s\"` draws trajectories of persons: fit with `id` and",
      "`class_level = \"group\"`."), what), class = "latents_bad_argument", call = NULL))
  }
  if (trajectory && !identical(what, "fitted")) {
    .gg_require()
    key <- .growth_class_key(x)
    return(switch(what,
      trajectories = .growth_plot_trajectories(x, key, time, facet, max_persons, 0,
                                               level, main, subtitle),
      individuals = .growth_plot_individuals(x, key, time, persons, main, subtitle),
      coefficients = .growth_plot_coefficients(x, key, level, main, subtitle),
      classification = ,
      posteriors = .growth_plot_classification(x, key, main, subtitle)))
  }
  .gg_require()
  drawn <- switch(what,
    fitted = .mixture_plot_fitted(x, predictor, main, subtitle),
    coefficients = .mixture_plot_coefficients(x, level, main, subtitle),
    posteriors = .mixture_plot_posteriors(x, main, subtitle))
  .gg_draw_on_call(drawn)
}

#' The views of a mixture regression that classifies persons
#' @return A character vector of `what` values.
#' @noRd
.mixture_trajectory_views <- function() {
  c("trajectories", "individuals", "classification")
}

#' The views a mixture regression's plot() draws
#'
#' Every fit draws the fitted values, coefficients and posteriors; a fit that
#' classifies persons (`class_level = "group"`) also draws trajectories.
#' @param x A `latents_mixture_regression` fit.
#' @return A character vector of `what` values.
#' @noRd
.mixture_plot_views <- function(x) {
  c("fitted", "coefficients", "posteriors",
    if (identical(x$spec$nesting, "group")) .mixture_trajectory_views())
}

#' @noRd
.mixture_plot_fitted <- function(x, predictor, main, subtitle) {
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
  classes <- paste("class", seq_len(k))
  modal <- max.col(x$expectation$tau, ties.method = "first")
  observed <- if (identical(spec$family, "binomial") && any(spec$trials > 1)) {
    spec$y / spec$trials
  } else spec$y
  grid_values <- seq(min(frame[[predictor]]), max(frame[[predictor]]),
                     length.out = 100L)
  typical <- .mixture_typical_rows(frame, predictor, grid_values)
  new_spec <- .mixture_new_spec(x, typical, with_response = FALSE)
  means <- .mixture_class_means(new_spec, x$params)
  if (identical(spec$family, "binomial")) means <- stats::plogis(
    .mixture_linear_predictors(new_spec, x$params))
  points <- data.frame(predictor_value = frame[[predictor]], outcome = observed,
                       class_label = factor(classes[modal], levels = classes))
  curves <- data.frame(predictor_value = rep(grid_values, times = k),
                       outcome = as.vector(means),
                       class_label = factor(rep(classes, each = length(grid_values)),
                                            levels = classes))
  ggplot2::ggplot(points, ggplot2::aes(predictor_value, outcome, colour = class_label,
                                       fill = class_label, shape = class_label)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.45) +
    ggplot2::geom_line(data = curves, ggplot2::aes(linetype = class_label),
                       linewidth = 1) +
    .mixture_class_scales(classes) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = predictor, y = spec$response_name) +
    .gg_titles(main, subtitle, sprintf("%d-class mixture regression", k),
               "Rows marked by their most likely class; a line per class")
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
  # A model without other predictors (`y ~ 1`) still needs one row per grid
  # value.
  rows <- if (length(typical) == 0L) data.frame(row.names = seq_along(grid_values)) else
    stats::setNames(as.data.frame(lapply(typical, rep, length(grid_values)),
                                  stringsAsFactors = FALSE), names(predictors))
  rows[[predictor]] <- grid_values
  rows
}

#' @noRd
.mixture_plot_coefficients <- function(x, level, main, subtitle) {
  table <- .mixture_table(x, "coefficients", level)
  terms <- unique(table$term)
  classes <- unique(table$class)
  # Classes side by side within each term, the first term at the top.
  offset <- (match(table$class, classes) - (length(classes) + 1) / 2) * 0.18
  table$position <- length(terms) + 1L - match(table$term, terms) + offset
  table$class_label <- factor(sub("_", " ", table$class),
                              levels = sub("_", " ", classes))
  ggplot2::ggplot(table, ggplot2::aes(estimate, position, colour = class_label,
                                      fill = class_label, shape = class_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "#999999") +
    ggplot2::geom_segment(ggplot2::aes(x = conf_low, xend = conf_high,
                                       yend = position), linewidth = 0.8,
                         show.legend = FALSE) +
    ggplot2::geom_point(size = 2.8, colour = "grey20") +
    .mixture_class_scales(levels(table$class_label), lines = FALSE) +
    ggplot2::scale_y_continuous(breaks = rev(seq_along(terms)), labels = terms) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = sprintf("Estimate with %g%% interval", 100 * level), y = NULL) +
    .gg_titles(main, subtitle, "Class regression coefficients", NULL)
}

#' @noRd
.mixture_plot_posteriors <- function(x, main, subtitle) {
  spec <- x$spec
  posterior <- if (identical(spec$nesting, "group")) x$expectation$group_tau else
    x$expectation$tau
  k <- spec$n_classes
  classes <- paste("class", seq_len(k))
  modal <- max.col(posterior, ties.method = "first")
  certainty <- posterior[cbind(seq_len(nrow(posterior)), modal)]
  # A deterministic spread (golden-ratio offsets of each unit's rank within its
  # class), as the fits' own case strips use: no random jitter, no RNG drawn.
  within <- stats::ave(certainty, modal, FUN = \(v) rank(v, ties.method = "first"))
  cases <- data.frame(class_label = factor(classes[modal], levels = classes),
                      position = modal + 0.5 * ((within * 0.618) %% 1 - 0.5),
                      certainty = certainty)
  ggplot2::ggplot(cases, ggplot2::aes(position, certainty, colour = class_label,
                                      fill = class_label, shape = class_label)) +
    ggplot2::geom_hline(yintercept = 1 / k, linetype = "dotted", colour = "#999999") +
    ggplot2::geom_point(size = 1.6, alpha = 0.6) +
    .mixture_class_scales(classes, lines = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_len(k), labels = classes) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = NULL, y = "Largest posterior probability") +
    .gg_titles(main, subtitle, "Classification certainty",
               "Each unit's largest posterior probability; dotted line: chance")
}

#' Colour, fill, shape and line-type scales for mixture classes
#' @param classes Class labels, in order.
#' @param lines Whether to add the line-type scale.
#' @return A list of ggplot2 scales.
#' @noRd
.mixture_class_scales <- function(classes, lines = TRUE) {
  colours <- stats::setNames(rep_len(.gg_okabe_ito, length(classes)), classes)
  scales <- list(
    ggplot2::scale_colour_manual(values = colours),
    ggplot2::scale_fill_manual(values = colours),
    ggplot2::scale_shape_manual(values = stats::setNames(
      rep_len(.gg_shapes, length(classes)), classes)))
  if (!lines) return(scales)
  c(scales, list(ggplot2::scale_linetype_manual(values = stats::setNames(
    rep_len(c("solid", "dashed", "dotdash", "dotted", "longdash"),
            length(classes)), classes))))
}
