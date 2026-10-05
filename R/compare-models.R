# Comparing fitted models on the same data, in one call.

utils::globalVariables(c("model", "criterion", "value"))

#' Compare fitted models on the same data
#'
#' One row per model with its log likelihood, number of parameters,
#' information criteria, classification quality and convergence, so
#' alternative specifications (for example trajectory classes alone against
#' a growth mixture model) are compared without assembling tables by hand.
#' `delta_bic` is each model's BIC above the smallest, and `bic_weight` the
#' Schwarz weight, exp(-delta/2) normalized over the models compared (Wagenmakers
#' and Farrell 2004): the approximate posterior probability of each model if
#' one of them is true, with equal prior probabilities.
#'
#' Information criteria are comparable only on the same data: every model
#' must be fitted to the same outcome, persons (or groups) and observations,
#' or the comparison is refused.
#'
#' @param ... Two or more fits from [mixture_regression()] (with or without
#'   `random`). Name them to label the rows; unnamed fits are labelled by
#'   their expression.
#' @return A base `data.frame` of class `latents_comparison` with one row per
#'   model: `model`, `n_classes`, `random`, `n_parameters`, `log_likelihood`,
#'   `aic`, `bic`, `sabic`, `icl`, `entropy`, `smallest_share`, `delta_bic`,
#'   `bic_weight`, `converged`. Raises `latents_bad_argument` for fits it
#'   cannot compare and `latents_incomparable_models` for fits to different
#'   data.
#' @references
#' Wagenmakers, E.-J., & Farrell, S. (2004). AIC model selection using Akaike
#' weights. *Psychonomic Bulletin & Review*, 11, 192--196.
#' @examples
#' \donttest{
#' trajectories <- mixture_regression(attendance ~ sequence, course_engagement,
#'                                    n_classes = 2, id = "student",
#'                                    class_level = "group", seed = 1)
#' growth <- mixture_regression(attendance ~ sequence, course_engagement,
#'                              n_classes = 2, id = "student",
#'                              class_level = "group", random = "intercept",
#'                              seed = 1)
#' compare_models(trajectories = trajectories, growth = growth)
#' }
#' @export
compare_models <- function(...) {
  fits <- list(...)
  expressions <- vapply(as.list(substitute(list(...)))[-1L], deparse1, character(1))
  labels <- names(fits) %||% rep("", length(fits))
  labels[!nzchar(labels)] <- expressions[!nzchar(labels)]
  supported <- vapply(fits, function(fit) {
    inherits(fit, "latents_mixture_regression") || inherits(fit, "latents_growth_mixture")
  }, logical(1))
  if (length(fits) < 2L || !all(supported)) {
    stop(errorCondition(paste(
      "`compare_models()` takes two or more fits from mixture_regression(),",
      "with or without `random`."), class = "latents_bad_argument", call = NULL))
  }
  if (anyDuplicated(labels)) {
    stop(errorCondition("Every model needs a distinct name.",
                        class = "latents_bad_argument", call = NULL))
  }
  signature <- function(fit) {
    spec <- fit$spec
    list(y = spec$y, trials = spec$trials, family = spec$family,
         categories = spec$category_levels, groups = spec$group_index,
         group_levels = spec$group_levels, clusters = spec$cluster_index,
         cluster_levels = spec$cluster_levels,
         weights = spec$sampling_weights, rows = spec$kept_rows,
         response = spec$response_name)
  }
  reference <- signature(fits[[1L]])
  same <- vapply(fits, function(fit) isTRUE(all.equal(signature(fit), reference)),
                 logical(1))
  if (!all(same)) {
    stop(errorCondition(sprintf(paste(
      "%s %s fitted to different data from %s: information criteria compare",
      "models only on the same outcome and observations."),
      paste(labels[!same], collapse = ", "), if (sum(!same) > 1L) "are" else "is",
      labels[1L]), class = "latents_incomparable_models", call = NULL))
  }
  rows <- lapply(seq_along(fits), function(i) {
    fit_row <- get_results(fits[[i]], "fit")
    data.frame(model = labels[i], n_classes = as.integer(fit_row$n_classes),
               random = .growth_random_description(fits[[i]]$spec),
               n_parameters = as.integer(fit_row$n_parameters),
               log_likelihood = fit_row$log_likelihood, aic = fit_row$aic,
               bic = fit_row$bic, sabic = fit_row$sabic, icl = fit_row$icl,
               entropy = fit_row$entropy, smallest_share = fit_row$smallest_share,
               converged = isTRUE(fit_row$converged), stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  table$delta_bic <- table$bic - min(table$bic)
  relative <- exp(-table$delta_bic / 2)
  table$bic_weight <- relative / sum(relative)
  table <- table[c("model", "n_classes", "random", "n_parameters", "log_likelihood",
                   "aic", "bic", "sabic", "icl", "entropy", "smallest_share",
                   "delta_bic", "bic_weight", "converged")]
  column <- .latents_column
  table <- .latents_table(table, title = "Model comparison", display = list(
    column("Model", "model", "text"), column("Classes", "n_classes", "integer"),
    column("Random", "random", "text"), column("k", "n_parameters", "integer"),
    column("LogLik", "log_likelihood"), column("BIC", "bic"),
    column(if (isTRUE(l10n_info()$`UTF-8`)) "\u0394BIC" else "dBIC", "delta_bic"), column("Weight", "bic_weight", "share"),
    column("Entropy", "entropy")))
  attr(table, "marks") <- ifelse(table$delta_bic == 0, "  <- best", "")
  class(table) <- c("latents_comparison", class(table))
  table
}

#' Plot a model comparison
#'
#' Each model's AIC, BIC and sample-size adjusted BIC, lowest best; the
#' model with the smallest BIC is marked.
#'
#' @param x A `latents_comparison` from [compare_models()].
#' @param main,subtitle Title and subtitle.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @export
plot.latents_comparison <- function(x, main = NULL, subtitle = NULL, ...) {
  .gg_require()
  table <- as.data.frame(x)
  criteria <- c(aic = "AIC", bic = "BIC", sabic = "SABIC")
  long <- do.call(rbind, lapply(names(criteria), function(column) {
    data.frame(model = factor(table$model, levels = table$model),
               criterion = factor(criteria[[column]], levels = criteria),
               value = table[[column]], best = table[[column]] == min(table[[column]]),
               stringsAsFactors = FALSE)
  }))
  ggplot2::ggplot(long, ggplot2::aes(model, value, group = criterion,
                                     colour = criterion, shape = criterion)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_point(data = long[long$best, , drop = FALSE], size = 5.2, shape = 21,
                        fill = NA, colour = "black", stroke = 0.6) +
    ggplot2::scale_colour_manual(values = .gg_okabe_ito[c(1L, 4L, 3L)]) +
    ggplot2::scale_shape_manual(values = c(16L, 17L, 15L)) +
    ggplot2::labs(x = NULL, y = "Information criterion (lower is better)") +
    .gg_theme() + .gg_top_legend() +
    .gg_titles(main, subtitle, "Model comparison",
               "Each criterion's best model is circled")
}
