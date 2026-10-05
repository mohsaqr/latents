# Plots of the growth mixture model: class trajectories with the spread of
# persons around them, individual fits, coefficients, predicted random
# effects and classification.

utils::globalVariables(c("trajectory_time", "observed", "predicted", "class_mean",
                         "class_label", "person_label", "effect_x", "effect_y",
                         "estimate", "conf_low", "conf_high", "label_y", "probability",
                         "person", "spread_low", "spread_high", "observed_mean",
                         "observed_low", "observed_high", "value_label",
                         "position", "posterior_class", "note"))

#' One row per class: share, display rank, label, colour, line type and shape
#' @noRd
.growth_class_key <- function(x) {
  posterior <- .trajectory_view(x)$posterior
  weights <- x$spec$sampling_weights %||% rep(1, nrow(posterior))
  share <- colSums(posterior * weights) / sum(weights)
  n_classes <- length(share)
  ordering <- order(-share, seq_len(n_classes))
  rank <- integer(n_classes)
  rank[ordering] <- seq_len(n_classes)
  label <- sprintf("Class %d (%s)", seq_len(n_classes), .gg_percent(share))
  data.frame(class = .growth_class_names(x), rank = rank, share = share,
             label = factor(label, levels = label[ordering]),
             colour = rep_len(.gg_okabe_ito, n_classes)[rank],
             linetype = rep_len(c("solid", "dashed", "dotdash", "longdash", "twodash",
                                  "dotted"), n_classes)[rank],
             shape = rep_len(.gg_shapes, n_classes)[rank],
             stringsAsFactors = FALSE)
}

#' Manual colour, fill, line-type and shape scales for the class key
#' @noRd
.growth_scales <- function(key, legend = FALSE) {
  levels <- levels(key$label)
  at <- match(levels, as.character(key$label))
  guide <- if (isTRUE(legend)) "legend" else "none"
  values <- function(column) stats::setNames(key[[column]][at], levels)
  list(ggplot2::scale_colour_manual(values = values("colour"), guide = guide),
       ggplot2::scale_fill_manual(values = values("colour"), guide = guide),
       ggplot2::scale_linetype_manual(values = values("linetype"), guide = guide),
       ggplot2::scale_shape_manual(values = values("shape"), guide = guide))
}

#' Plot a growth mixture model
#'
#' @param x A `latents_growth_mixture` fit.
#' @param what `"trajectories"` (each class's mean trajectory with its
#'   confidence band, the spread of its persons' own trajectories, and the
#'   observed class means at each time), `"individuals"` (a sample of persons
#'   across each class's range of classification certainty: observations,
#'   their class trajectory dashed, their own predicted curve solid),
#'   `"coefficients"` (estimates with confidence intervals and p-values, by
#'   term), `"random"` (predicted random effects with the class's
#'   model-implied 95% ellipse; degenerate covariances are flagged), or
#'   `"classification"` (one bar per person, split by their posterior class
#'   probabilities). `"posteriors"` is a synonym of `"classification"`.
#' @param time The numeric variable to draw along; by default the first
#'   numeric variable of `random`.
#' @param facet For `"trajectories"`, `TRUE` draws one panel per class.
#' @param persons For `"individuals"`, how many persons to show, or a vector
#'   of their identifiers.
#' @param max_persons For `"trajectories"`, the most persons drawn behind the
#'   class trajectories (an even sample by class); `0` draws none.
#' @param spread For `"trajectories"`, the share of persons' class-implied
#'   trajectories the light band holds (from the random effects); `0` hides
#'   it.
#' @param level Confidence level of the bands and intervals.
#' @param main,subtitle Title and subtitle; defaults describe the view.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @examples
#' few_students <- subset(growth_scores, student <= 60)
#' fit <- mixture_regression(score ~ wave, few_students, n_classes = 2,
#'                           id = "student", class_level = "group",
#'                           random = "intercept", seed = 1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
#' @export
plot.latents_growth_mixture <- function(x, what = c("trajectories", "individuals",
                                                    "coefficients", "random",
                                                    "classification", "posteriors"),
                                        time = NULL, facet = FALSE, persons = 12L,
                                        max_persons = 120L, spread = 0.8, level = 0.95,
                                        main = NULL, subtitle = NULL, ...) {
  .gg_require()
  what <- match.arg(what)
  stopifnot("`spread` must be a single number in [0, 1)" =
              is.numeric(spread) && length(spread) == 1L && spread >= 0 && spread < 1,
            "`facet` must be TRUE or FALSE" = isTRUE(facet) || isFALSE(facet))
  key <- .growth_class_key(x)
  switch(what,
    trajectories = .growth_plot_trajectories(x, key, time, facet, max_persons, spread,
                                             level, main, subtitle),
    individuals = .growth_plot_individuals(x, key, time, persons, main, subtitle),
    coefficients = .growth_plot_coefficients(x, key, level, main, subtitle),
    random = .growth_plot_random(x, key, main, subtitle),
    classification = ,
    posteriors = .growth_plot_classification(x, key, main, subtitle))
}

#' Posterior-weighted observed means of every class at each observed time
#'
#' Drawn only when there are at most 30 distinct times; a mean needs at least
#' one person's worth of posterior weight.
#' @noRd
.growth_observed_means <- function(x, time, level) {
  view <- .trajectory_view(x)
  spec <- view$spec
  posterior <- view$posterior[spec$group_index, , drop = FALSE]
  outcome <- .trajectory_observed(view)
  times <- spec$model_data[[time]]
  distinct <- sort(unique(times))
  if (length(distinct) > 30L) return(NULL)
  z <- stats::qnorm(1 - (1 - level) / 2)
  do.call(rbind, lapply(seq_len(ncol(posterior)), function(k) {
    do.call(rbind, lapply(distinct, function(t) {
      at <- times == t
      w <- posterior[at, k] *
        (spec$sampling_weights %||% rep(1, spec$n_groups))[spec$group_index[at]]
      total <- sum(w)
      if (total < 1) return(NULL)
      mean <- sum(w * outcome[at]) / total
      effective <- total^2 / sum(w^2)
      sd <- sqrt(sum(w * (outcome[at] - mean)^2) / total)
      data.frame(class = .growth_class_names(x)[k], trajectory_time = t,
                 observed_mean = mean, observed_low = mean - z * sd / sqrt(effective),
                 observed_high = mean + z * sd / sqrt(effective), stringsAsFactors = FALSE)
    }))
  }))
}

#' Where `spread` of each class's persons' trajectories lie
#'
#' From the random effects alone: mean +- z sqrt(z(t)' G z(t)), with other
#' random-effect variables at their typical values.
#' @noRd
.growth_spread <- function(x, bands, time, spread) {
  if (spread <= 0) return(NULL)
  grid <- unique(bands$trajectory_time)
  variables <- setdiff(all.vars(x$spec$random), x$spec$id)
  rows <- if (length(variables) == 0L) data.frame(row = seq_along(grid)) else
    .mixture_typical_rows(x$spec$model_data[unique(c(time, variables))], time, grid)
  random_design <- stats::model.matrix(x$spec$random, rows)
  covariances <- .growth_covariances(x$spec, x$params)
  z <- stats::qnorm(1 - (1 - spread) / 2)
  do.call(rbind, lapply(seq_along(covariances), function(k) {
    variance <- rowSums((random_design %*% covariances[[k]]) * random_design)
    mean <- bands$estimate[bands$class == .growth_class_names(x)[k]]
    data.frame(class = .growth_class_names(x)[k], trajectory_time = grid,
               spread_low = mean - z * sqrt(variance), spread_high = mean + z * sqrt(variance),
               stringsAsFactors = FALSE)
  }))
}

#' @noRd
.growth_plot_trajectories <- function(x, key, time, facet, max_persons, spread, level,
                                      main, subtitle) {
  time <- .growth_time_variable(x, time)
  bands <- as.data.frame(.growth_trajectory_table(x, level, .trajectory_inference(x),
                                                  time))
  bands$trajectory_time <- bands[[time]]
  label_of <- function(class) key$label[match(class, key$class)]
  bands$class_label <- label_of(bands$class)
  spread_band <- if (inherits(x, "latents_growth_mixture")) {
    .growth_spread(x, bands, time, spread)
  }
  means <- .growth_observed_means(x, time, level)
  observed <- as.data.frame(.growth_individual_table(x, time))
  observed$trajectory_time <- observed[[time]]
  observed$person <- observed[[x$spec$id]]
  people <- unique(observed[c("person", "class")])
  per_class <- floor(max_persons / nrow(key))
  shown <- if (per_class < 1L) NULL else
    unlist(lapply(split(people$person, people$class), function(ids) {
      ids[unique(round(seq(1, length(ids), length.out = min(length(ids), per_class))))]
    }), use.names = FALSE)
  observed <- observed[observed$person %in% shown, , drop = FALSE]
  observed <- observed[order(observed$person, observed$trajectory_time), , drop = FALSE]
  observed$class_label <- label_of(observed$class)
  layers <- list()
  if (nrow(observed) > 0L) {
    layers$persons <- ggplot2::geom_line(
      data = observed, ggplot2::aes(trajectory_time, observed, group = person,
                                    colour = class_label),
      alpha = 0.10, linewidth = 0.3)
  }
  if (!is.null(spread_band)) {
    spread_band$class_label <- label_of(spread_band$class)
    layers$spread <- ggplot2::geom_ribbon(
      data = spread_band, ggplot2::aes(trajectory_time, ymin = spread_low,
                                       ymax = spread_high, fill = class_label),
      alpha = 0.10)
  }
  layers$band <- ggplot2::geom_ribbon(
    data = bands, ggplot2::aes(trajectory_time, ymin = conf_low, ymax = conf_high,
                               fill = class_label), alpha = 0.30)
  layers$mean <- ggplot2::geom_line(
    data = bands, ggplot2::aes(trajectory_time, estimate, colour = class_label,
                               linetype = class_label), linewidth = 1.3)
  if (!is.null(means)) {
    means$class_label <- label_of(means$class)
    layers$observed_errors <- ggplot2::geom_linerange(
      data = means, ggplot2::aes(trajectory_time, ymin = observed_low, ymax = observed_high,
                                 colour = class_label), linewidth = 0.5)
    layers$observed_means <- ggplot2::geom_point(
      data = means, ggplot2::aes(trajectory_time, observed_mean, fill = class_label,
                                 shape = class_label), colour = "white", size = 2.4,
      stroke = 0.5)
  }
  caption <- paste(c(
    sprintf("Line and dark band: class mean trajectory and its %s%% confidence band.",
            format(100 * level)),
    if (!is.null(spread_band))
      sprintf("Light band: where %s%% of the class's persons' trajectories lie.",
              format(100 * spread)),
    if (!is.null(means)) "Points: observed class means (posterior-weighted) with intervals.",
    if (nrow(observed) > 0L) "Thin lines: persons, coloured by their most likely class."),
    collapse = "\n")
  plot <- ggplot2::ggplot() + layers + .growth_scales(key) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(x = time, y = x$spec$response_name) + .gg_theme()
  description <- if (inherits(x, "latents_growth_mixture")) {
    sprintf("%d persons, %d classes; random %s (%s covariance)", length(x$stats$n),
            x$spec$n_classes, .growth_random_description(x$spec), x$spec$random_covariance)
  } else {
    sprintf("%d persons, %d classes; no random effects (latent class growth)",
            x$spec$n_groups, x$spec$n_classes)
  }
  if (isTRUE(facet)) {
    return(plot + ggplot2::facet_wrap(~ class_label) +
             .gg_titles(main, subtitle, "Class trajectories", description,
                        caption = caption))
  }
  ends <- bands[bands$trajectory_time == max(bands$trajectory_time), , drop = FALSE]
  span <- diff(range(c(bands$conf_low, bands$conf_high,
                       spread_band$spread_low, spread_band$spread_high)))
  ends$label_y <- .gg_spread(ends$estimate, 0.07 * span)
  right <- max(bands$trajectory_time)
  width <- diff(range(bands$trajectory_time))
  plot +
    ggplot2::geom_text(data = ends,
                       ggplot2::aes(right + 0.02 * width, label_y, label = class_label,
                                    colour = class_label),
                       hjust = 0, size = 3.7, fontface = "bold") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 115, 10, 12)) +
    .gg_titles(main, subtitle, "Class trajectories", description, caption = caption)
}

#' @noRd
.growth_plot_individuals <- function(x, key, time, persons, main, subtitle) {
  time <- .growth_time_variable(x, time)
  rows <- as.data.frame(.growth_individual_table(x, time))
  rows$trajectory_time <- rows[[time]]
  rows$person <- rows[[x$spec$id]]
  assignments <- as.data.frame(.growth_assignment_table(x))
  if (length(persons) == 1L && is.numeric(persons)) {
    # Persons spread over each class's range of certainty, from its most to
    # its least clearly classified, so uncertain members are shown too.
    by_class <- split(assignments, assignments$class)
    per_class <- ceiling(persons / length(by_class))
    chosen <- unlist(lapply(by_class, function(members) {
      ranked <- members[[x$spec$id]][order(-members$probability)]
      ranked[unique(round(seq(1, length(ranked), length.out = min(per_class,
                                                                   length(ranked)))))]
    }), use.names = FALSE)
    chosen <- utils::head(chosen, persons)
  } else chosen <- persons
  rows <- rows[rows$person %in% chosen, , drop = FALSE]
  if (nrow(rows) == 0L) {
    stop(errorCondition("None of `persons` is in the fitted data.",
                        class = "latents_bad_argument", call = NULL))
  }
  rows$class_label <- key$label[match(rows$class, key$class)]
  probability <- assignments$probability[match(rows$person, assignments[[x$spec$id]])]
  rows$person_label <- factor(sprintf("%s %s \u00b7 %s, p = %.2f", x$spec$id, rows$person,
                                      sub("class_", "class ", rows$class), probability))
  rows <- rows[order(rows$person, rows$trajectory_time), , drop = FALSE]
  ggplot2::ggplot(rows, ggplot2::aes(trajectory_time, colour = class_label)) +
    ggplot2::geom_line(ggplot2::aes(y = class_mean), linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = predicted), linewidth = 0.95) +
    ggplot2::geom_point(ggplot2::aes(y = observed, fill = class_label, shape = class_label),
                        colour = "white", size = 2.2, stroke = 0.4) +
    ggplot2::facet_wrap(~ person_label) +
    .growth_scales(key, legend = TRUE) +
    ggplot2::labs(x = time, y = x$spec$response_name) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::guides(linetype = "none", fill = "none", shape = "none") +
    .gg_titles(main, subtitle, "Individual trajectories",
               paste("Points are observations; the dashed line is the class trajectory,",
                     "the solid line the person's own predicted curve"))
}

#' @noRd
.growth_plot_coefficients <- function(x, key, level, main, subtitle) {
  table <- as.data.frame(get_results(x, "coefficients", level = level))
  class_rows <- table$class != "common"
  shared <- "Shared by all classes"
  table$class_label <- ifelse(class_rows,
                              as.character(key$label[match(table$class, key$class)]),
                              shared)
  levels <- c(levels(key$label), if (any(!class_rows)) shared)
  table$class_label <- factor(table$class_label, levels = rev(levels))
  # One number of decimals per panel (term), as in the tables, at most three.
  decimals <- vapply(split(table$estimate, table$term), function(values) {
    min(.latents_column_decimals(values), 3L)
  }, integer(1))[table$term]
  number <- function(values) sprintf("%.*f", unname(decimals), values)
  table$value_label <- sprintf("%s [%s, %s]  %s", number(table$estimate),
                               number(table$conf_low), number(table$conf_high),
                               ifelse(table$p_value < 0.001, "p < .001",
                                      sprintf("p = %.3f", table$p_value)))
  colours <- c(stats::setNames(as.character(key$colour), as.character(key$label)),
               stats::setNames("grey35", shared))
  shapes <- c(stats::setNames(key$shape, as.character(key$label)),
              stats::setNames(21L, shared))
  ggplot2::ggplot(table, ggplot2::aes(estimate, class_label, colour = class_label)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.4,
                        linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = conf_low, xmax = conf_high), width = 0.16,
                           linewidth = 0.7, orientation = "y") +
    ggplot2::geom_point(ggplot2::aes(fill = class_label, shape = class_label),
                        colour = "white", size = 3.2, stroke = 0.5) +
    # The estimate and interval read at each row's right edge, inside the panel.
    ggplot2::geom_text(ggplot2::aes(label = value_label), x = Inf, hjust = 1.02,
                       vjust = -1.25, size = 3.05, show.legend = FALSE) +
    ggplot2::facet_wrap(~ term, scales = "free_x") +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::scale_shape_manual(values = shapes, guide = "none") +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5, 0.95))) +
    ggplot2::labs(x = "Estimate", y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.3),
                   panel.grid.major.y = ggplot2::element_blank()) +
    .gg_titles(main, subtitle, "Trajectory coefficients",
               sprintf("Estimates with %s%% confidence intervals and Wald p-values",
                       format(100 * level)))
}

#' @noRd
.growth_plot_random <- function(x, key, main, subtitle) {
  effects <- as.data.frame(.growth_random_effects_table(x))
  terms <- colnames(x$spec$random_design)
  effects$class_label <- key$label[match(effects$class, key$class)]
  covariances <- .growth_covariances(x$spec, x$params)
  flagged <- key$class[vapply(seq_len(nrow(key)), function(k) {
    .growth_class_is_degenerate(x, k)
  }, logical(1))]
  notes <- if (length(flagged) > 0L) {
    data.frame(class_label = key$label[match(flagged, key$class)],
               note = "Degenerate covariance: not separable from zero or each other",
               stringsAsFactors = FALSE)
  }
  if (length(terms) == 1L) {
    effects$effect_x <- effects[[terms[1L]]]
    plot <- ggplot2::ggplot(effects, ggplot2::aes(effect_x, class_label,
                                                  colour = class_label)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
      ggplot2::geom_jitter(ggplot2::aes(fill = class_label, shape = class_label),
                           colour = "white", height = 0.18, width = 0, alpha = 0.85,
                           size = 2, stroke = 0.3) +
      .growth_scales(key) +
      ggplot2::labs(x = sprintf("Predicted random effect: %s", terms[1L]), y = NULL) +
      .gg_theme() +
      .gg_titles(main, subtitle, "Predicted random effects",
                 "Each point is a person's deviation from their class trajectory")
    if (!is.null(notes)) {
      plot <- plot + ggplot2::geom_text(data = notes, ggplot2::aes(label = note),
                                        x = Inf, hjust = 1.02, vjust = -1.2, size = 3.1,
                                        colour = "#D55E00", inherit.aes = FALSE)
    }
    return(plot)
  }
  effects$effect_x <- effects[[terms[1L]]]
  effects$effect_y <- effects[[terms[2L]]]
  ellipses <- do.call(rbind, lapply(seq_len(nrow(key)), function(k) {
    shape <- .gg_ellipse(c(0, 0), covariances[[k]][1:2, 1:2])
    data.frame(effect_x = shape$x, effect_y = shape$y, class_label = key$label[k],
               stringsAsFactors = FALSE)
  }))
  plot <- ggplot2::ggplot(effects, ggplot2::aes(effect_x, effect_y, colour = class_label)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_path(data = ellipses, ggplot2::aes(linetype = class_label),
                       linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(fill = class_label, shape = class_label),
                        colour = "white", alpha = 0.85, size = 1.9, stroke = 0.3) +
    ggplot2::facet_wrap(~ class_label) +
    .growth_scales(key) +
    ggplot2::labs(x = sprintf("Random effect: %s", terms[1L]),
                  y = sprintf("Random effect: %s", terms[2L])) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.3)) +
    .gg_titles(main, subtitle, "Predicted random effects",
               paste("Each point is a person's deviation from their class trajectory;",
                     "the ellipse holds 95% of the class's random effects"))
  if (!is.null(notes)) {
    plot <- plot + ggplot2::geom_text(data = notes, ggplot2::aes(label = note),
                                      x = -Inf, y = Inf, hjust = -0.03, vjust = 1.6,
                                      size = 3.1, colour = "#D55E00", inherit.aes = FALSE)
  }
  plot
}

#' @noRd
.growth_plot_classification <- function(x, key, main, subtitle) {
  posterior <- .trajectory_view(x)$posterior
  classes <- .growth_class_names(x)
  modal <- max.col(posterior, ties.method = "first")
  certainty <- posterior[cbind(seq_along(modal), modal)]
  # Persons ordered by their class (largest first) and, within it, certainty.
  ordering <- order(key$rank[modal], -certainty)
  position <- integer(length(modal))
  position[ordering] <- seq_along(ordering)
  long <- data.frame(
    position = rep(position, times = length(classes)),
    class_label = rep(key$label[modal], times = length(classes)),
    posterior_class = rep(key$label, each = length(modal)),
    probability = as.vector(posterior), stringsAsFactors = FALSE)
  average <- vapply(seq_along(classes), function(k) mean(certainty[modal == k]), numeric(1))
  subtitle_text <- paste0("One bar per person, split by their class probabilities. ",
                          "Average probability of the assigned class: ",
                          paste(sprintf("%s %.2f", sub(" \\(.*", "", key$label), average),
                                collapse = ", "))
  ggplot2::ggplot(long, ggplot2::aes(position, probability, fill = posterior_class)) +
    ggplot2::geom_col(width = 1, colour = NA) +
    ggplot2::facet_grid(~ class_label, scales = "free_x", space = "free_x") +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    .growth_scales(key, legend = TRUE)[2L] +
    ggplot2::labs(x = "Persons, by assigned class and certainty",
                  y = "Posterior probability") +
    .gg_theme() + .gg_top_legend() +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(6, "pt"),
                   panel.grid.major.y = ggplot2::element_blank()) +
    .gg_titles(main, subtitle, "Classification", subtitle_text)
}
