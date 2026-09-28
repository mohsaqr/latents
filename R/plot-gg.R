# The drawing engine behind every plot() method of the fitted families.
#
# Each view returns a ggplot object. Views read the same internal helpers the
# tables use -- .multilpa_standardization(), .multilpa_mean_error_matrix(),
# .multilpa_plot_observations(), the fit's own posteriors and effective counts
# -- so a picture can never disagree with the numbers a user prints.
# Shared rules, applied by every view:
#   * one profile share: the posterior (effective) proportion;
#   * one display order: largest profile first, profile number breaking ties;
#   * one colour and one shape per profile, in that order (Okabe-Ito);
#   * direct labels where a legend would make the reader match colours.

utils::globalVariables(c(
  "x", "y", "ymin", "ymax", "xmin", "xmax", "xend", "yend", "lower", "upper",
  "mean", "value", "label", "profile_label", "indicator", "fill_value",
  "text_colour", "position", "group", "height", "criterion", "series",
  "n_profiles", "assigned_label", "posterior_label", "average_posterior",
  "q1", "q3", "median", "count", "bar_label", "case", "certainty", "panel",
  "uncertainty", "xvar", "yvar", "cases", "share", "fill", "y0", "row_label",
  "column_label", "class_label", "time_position", "cell_text", "column",
  "probability", "failed_y", "title", "block", "category_label",
  "classes_label"))

.gg_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                   "#CC79A7", "#999999", "#000000", "#F0E442")
.gg_shapes <- c(21L, 22L, 24L, 23L, 25L)

#' Stop unless ggplot2 is installed
#' @noRd
.gg_require <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(errorCondition(
      "Plots need the ggplot2 package: install.packages(\"ggplot2\").",
      class = "latents_missing_package", call = NULL))
  }
  invisible(TRUE)
}

#' Whether plots can be drawn, telling the caller once when they cannot
#'
#' diagnostics() and report() draw by default; without ggplot2 they still
#' print everything else and say why nothing was drawn.
#' @noRd
.gg_available <- function() {
  if (requireNamespace("ggplot2", quietly = TRUE)) return(TRUE)
  .multilpa_notice(
    "Plots need the ggplot2 package; install it, or pass `plots = FALSE`.",
    class = "latents_missing_package")
  FALSE
}

#' The house theme: white field, horizontal guides only, left-aligned titles
#' @noRd
.gg_theme <- function(base_size = 11) {
  el <- ggplot2::element_text
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = el(face = "bold", size = ggplot2::rel(1.15),
                      margin = ggplot2::margin(b = 3)),
      plot.subtitle = el(colour = "grey35", size = ggplot2::rel(0.88),
                         margin = ggplot2::margin(b = 10)),
      plot.caption = el(colour = "grey45", size = ggplot2::rel(0.78),
                        hjust = 0, margin = ggplot2::margin(t = 8)),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90",
                                                 linewidth = 0.35),
      axis.title = el(colour = "grey30", size = ggplot2::rel(0.88)),
      axis.text = el(colour = "grey25"),
      strip.text = el(face = "bold", hjust = 0, size = ggplot2::rel(0.95),
                      margin = ggplot2::margin(b = 4)),
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 16, 10, 12)
    )
}

#' Legend across the top, left-aligned, without a title
#' @noRd
.gg_top_legend <- function() {
  ggplot2::theme(legend.position = "top", legend.justification = "left",
                 legend.title = ggplot2::element_blank(),
                 legend.margin = ggplot2::margin(0, 0, 0, 0),
                 legend.key.size = ggplot2::unit(0.9, "lines"))
}

#' A percentage for a label, never rounded down to a misleading zero
#' @noRd
.gg_percent <- function(proportion) {
  percent <- 100 * proportion
  ifelse(percent < 1, "<1%", sprintf("%.0f%%", percent))
}

#' One row per profile, indexed by profile number: share, display rank,
#' label, colour and shape
#'
#' The share is the profile's effective count over all effective counts, as
#' every table reports it. A noise component's mass stays in the denominator,
#' so the shares of the real profiles are what they are in the fit.
#' @noRd
.gg_profile_key <- function(x) {
  n_profiles <- x$n_profiles
  counts <- x$effective_profile_counts
  share <- (counts / sum(counts))[seq_len(n_profiles)]
  ordering <- order(-share, seq_len(n_profiles))
  rank <- integer(n_profiles)
  rank[ordering] <- seq_len(n_profiles)
  label <- sprintf("Profile %d (%s)", seq_len(n_profiles), .gg_percent(share))
  data.frame(
    profile = seq_len(n_profiles),
    rank = rank,
    share = share,
    count = counts[seq_len(n_profiles)],
    label = factor(label, levels = label[ordering]),
    colour = rep_len(.gg_okabe_ito, n_profiles)[rank],
    shape = rep_len(.gg_shapes, n_profiles)[rank],
    stringsAsFactors = FALSE
  )
}

#' Manual profile scales for the aesthetics a view maps
#' @noRd
.gg_scales <- function(key, aesthetics = c("colour", "fill", "shape"),
                       legend = FALSE) {
  levels <- levels(key$label)
  at <- match(levels, as.character(key$label))
  colours <- stats::setNames(key$colour[at], levels)
  shapes <- stats::setNames(key$shape[at], levels)
  guide <- if (isTRUE(legend)) "legend" else "none"
  scales <- list(
    colour = ggplot2::scale_colour_manual(values = colours, guide = guide),
    fill = ggplot2::scale_fill_manual(values = colours, guide = guide),
    shape = ggplot2::scale_shape_manual(values = shapes, guide = guide)
  )
  scales[aesthetics]
}

#' Title and subtitle: the caller's when given, the view's default otherwise
#' @noRd
.gg_titles <- function(main, subtitle, default_main, default_subtitle,
                       caption = NULL) {
  ggplot2::labs(title = main %||% default_main,
                subtitle = subtitle %||% default_subtitle,
                caption = caption)
}

#' Push labels apart so no two are closer than `gap`, keeping their order
#'
#' Subtracting `gap * rank` turns "at least `gap` apart" into "non-decreasing",
#' which a cumulative maximum enforces; adding it back and re-centring on the
#' original mean keeps the block where the data are.
#' @noRd
.gg_spread <- function(y, gap) {
  if (length(y) < 2L) return(y)
  ordering <- order(y)
  steps <- gap * (seq_along(y) - 1)
  spread <- cummax(y[ordering] - steps) + steps
  spread <- spread - mean(spread) + mean(y)
  out <- numeric(length(y))
  out[ordering] <- spread
  out
}

#' Ink for text on a fill: dark on light fills, white on dark ones
#' @noRd
.gg_ink <- function(fills) {
  if (length(fills) == 0L) return(character())
  rgb <- grDevices::col2rgb(fills)
  luminance <- 0.299 * rgb[1L, ] + 0.587 * rgb[2L, ] + 0.114 * rgb[3L, ]
  ifelse(luminance > 150, "grey10", "white")
}

#' The white-to-blue ramp for a probability, as a colour string
#' @noRd
.gg_probability_fill <- function(values) {
  ramp <- grDevices::colorRamp(c("#FFFFFF", "#0072B2"))
  inside <- pmin(pmax(values, 0), 1)
  inside[is.na(inside)] <- 0
  grDevices::rgb(ramp(inside), maxColorValue = 255)
}

#' A shared probability fill scale from 0 to 1, with its key
#' @noRd
.gg_probability_scale <- function(name) {
  ggplot2::scale_fill_gradient(
    low = "#FFFFFF", high = "#0072B2", limits = c(0, 1), name = name,
    na.value = "#FFFFFF",
    guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(7, "lines"),
                                     barwidth = ggplot2::unit(0.7, "lines")))
}

#' Theme adjustments shared by the tile views
#' @noRd
.gg_tile_theme <- function() {
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 legend.position = "right",
                 legend.title = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                      colour = "grey30"),
                 axis.text.y = ggplot2::element_text(hjust = 1))
}

#' "2 group classes" when there are several, nothing when there is one
#' @noRd
.gg_group_class_text <- function(x) {
  if (isTRUE(x$n_group_classes > 1L)) {
    sprintf(", %d group classes", x$n_group_classes)
  } else ""
}

# -- measurement --------------------------------------------------------------

#' Profile means (and 95% intervals) in long form, on the requested scale
#'
#' @param errors A profiles-by-indicators matrix of standard errors, a string
#'   saying why there are none, or `NULL` when intervals were not asked for.
#' @return A list: `frame` (one row per profile and indicator), `vars`,
#'   `note` (the subtitle fragment about intervals) and `has_errors`.
#' @noRd
.gg_mean_frame <- function(x, key, scale, errors) {
  vars <- .multilpa_continuous_names(x)
  if (length(vars) == 0L) {
    stop(errorCondition(
      "This model has no continuous indicators; plot `what = \"responses\"` instead.",
      class = "latents_no_continuous", call = NULL))
  }
  reason <- if (is.character(errors)) errors else NULL
  if (is.character(errors)) errors <- NULL
  means <- x$means
  if (identical(scale, "standardized")) {
    basis <- .multilpa_standardization(x, vars)
    means <- sweep(sweep(means, 2L, basis$centre, "-"), 2L, basis$spread, "/")
    if (!is.null(errors)) errors <- sweep(errors, 2L, basis$spread, "/")
  }
  n_profiles <- x$n_profiles
  cells <- expand.grid(position = seq_along(vars), profile = seq_len(n_profiles))
  z <- stats::qnorm(0.975)
  frame <- data.frame(
    profile = cells$profile,
    position = cells$position,
    indicator = factor(vars[cells$position], levels = vars),
    mean = means[cbind(cells$profile, cells$position)])
  spread <- if (is.null(errors)) NA_real_ else
    errors[cbind(cells$profile, cells$position)]
  frame$lower <- frame$mean - z * spread
  frame$upper <- frame$mean + z * spread
  frame$profile_label <- key$label[frame$profile]
  frame$rank <- key$rank[frame$profile]
  note <- if (!is.null(errors) && any(is.finite(spread))) {
    "; whiskers are 95% intervals"
  } else if (!is.null(reason)) paste0("; ", reason) else ""
  list(frame = frame, vars = vars, note = note,
       has_errors = !is.null(errors) && any(is.finite(spread)))
}

#' Lines across indicators with direct labels at the right edge
#'
#' Shared by the profile means and the categorical response curves: both are
#' one series per profile across indicators, optionally with intervals.
#' @noRd
.gg_series_plot <- function(frame, key, labels_on, has_errors, y_label,
                            x_labels, ylim = NULL, reference = NULL) {
  n_x <- length(x_labels)
  k <- nrow(key)
  width <- if (has_errors) min(0.3, 0.08 * k) else 0
  frame$x <- frame$position +
    (frame$rank - (k + 1) / 2) * width / max(k - 1, 1)
  span <- range(c(frame$mean, frame$lower, frame$upper), na.rm = TRUE)
  last <- frame[frame$position == n_x, , drop = FALSE]
  last$y <- .gg_spread(last$mean, gap = 0.075 * max(diff(span), 1e-8))
  last$x <- n_x + 0.2 + width / 2
  plot <- ggplot2::ggplot(frame, ggplot2::aes(x, mean, colour = profile_label,
                                              group = profile_label))
  if (!is.null(reference)) {
    plot <- plot + ggplot2::geom_hline(yintercept = reference,
                                       colour = "grey55", linewidth = 0.4,
                                       linetype = "22")
  }
  if (has_errors) {
    plot <- plot + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper), width = 0.06,
      linewidth = 0.45, na.rm = TRUE)
  }
  plot <- plot +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = profile_label,
                                     fill = profile_label),
                        size = 2.6, colour = "white", stroke = 0.6)
  if (isTRUE(labels_on)) {
    plot <- plot + ggplot2::geom_text(
      data = last, ggplot2::aes(x = x, y = y, label = profile_label),
      hjust = 0, size = 3.5, fontface = "bold")
  }
  plot +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n_x), labels = x_labels,
      limits = c(0.7, n_x + if (isTRUE(labels_on)) 1.35 else 0.3),
      expand = ggplot2::expansion(0)) +
    .gg_scales(key, legend = !isTRUE(labels_on)) +
    ggplot2::coord_cartesian(ylim = ylim, clip = "off") +
    ggplot2::labs(x = NULL, y = y_label) +
    .gg_theme() +
    (if (!isTRUE(labels_on)) .gg_top_legend())
}

#' @noRd
.gg_view_profiles <- function(x, scale, labels, errors, main, subtitle) {
  key <- .gg_profile_key(x)
  means <- .gg_mean_frame(x, key, scale, errors)
  standardized <- identical(scale, "standardized")
  .gg_series_plot(means$frame, key, labels, means$has_errors,
                  y_label = if (standardized) "Standardized mean" else
                    "Estimated mean",
                  x_labels = means$vars,
                  reference = if (standardized) 0) +
    .gg_titles(main, subtitle,
               sprintf("Profile means across %d indicator%s",
                       length(means$vars),
                       if (length(means$vars) == 1L) "" else "s"),
               sprintf("%d profiles%s; %s scale%s", x$n_profiles,
                       .gg_group_class_text(x),
                       if (standardized) "standardized" else "input",
                       means$note))
}

#' @noRd
.gg_view_bars <- function(x, scale, errors, main, subtitle) {
  key <- .gg_profile_key(x)
  means <- .gg_mean_frame(x, key, scale, errors)
  dodge <- ggplot2::position_dodge(width = 0.82)
  plot <- ggplot2::ggplot(means$frame,
                          ggplot2::aes(indicator, mean, fill = profile_label)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_col(position = dodge, width = 0.78)
  if (means$has_errors) {
    plot <- plot + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper, group = profile_label),
      position = dodge, width = 0.22, linewidth = 0.4, colour = "grey20",
      na.rm = TRUE)
  }
  plot +
    .gg_scales(key, "fill", legend = TRUE) +
    .gg_titles(main, subtitle, "Profile means",
               sprintf("%d profiles; %s scale; bars start at zero%s",
                       x$n_profiles,
                       if (identical(scale, "standardized")) "standardized" else
                         "input", means$note)) +
    ggplot2::labs(x = NULL, y = if (identical(scale, "standardized"))
      "Standardized mean" else "Estimated mean") +
    .gg_theme() +
    .gg_top_legend()
}

#' Profile means in observed standard deviations, with a colour key
#'
#' The observed centre and spread of each indicator come from
#' .multilpa_standardization(), the same map as
#' `as.data.frame(scale = "standardized")`, so the tiles and the table agree.
#' @noRd
.gg_view_heatmap <- function(x, main, subtitle) {
  if (length(.multilpa_continuous_names(x)) == 0L) {
    return(.gg_view_response_heatmap(x, main, subtitle))
  }
  key <- .gg_profile_key(x)
  means <- .gg_mean_frame(x, key, "standardized", NULL)$frame
  limit <- max(abs(means$mean), 1e-8)
  means$text_colour <- ifelse(abs(means$mean) > 0.55 * limit, "white",
                              "grey15")
  means$profile_label <- factor(means$profile_label,
                                levels = rev(levels(key$label)))
  ggplot2::ggplot(means, ggplot2::aes(indicator, profile_label, fill = mean)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean),
                                    colour = text_colour), size = 3.4) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient2(
      low = "#D33F6A", mid = "white", high = "#4A6FE3", midpoint = 0,
      limits = c(-limit, limit), name = "SDs from\nobserved mean",
      guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(7, "lines"),
                                       barwidth = ggplot2::unit(0.7, "lines"))) +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    .gg_titles(main, subtitle, "Profile means in standard deviations",
               paste("Each profile mean minus the indicator's observed mean,",
                     "in observed standard deviations")) +
    ggplot2::labs(x = NULL, y = NULL) +
    .gg_theme() +
    .gg_tile_theme()
}

#' Every response probability as tiles, one block per indicator
#' @noRd
.gg_view_response_heatmap <- function(x, main, subtitle) {
  blocks <- x$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) {
    stop(errorCondition("This model has no indicators to draw.",
                        class = "latents_nothing_to_plot", call = NULL))
  }
  key <- .gg_profile_key(x)
  binary <- all(vapply(blocks, ncol, integer(1)) == 2L)
  # A binary indicator shows its last category, since the other is its
  # complement; a polytomous one shows every category, so a profile's whole
  # response distribution reads across its row.
  cells <- do.call(rbind, lapply(names(blocks), \(indicator) {
    block <- blocks[[indicator]]
    kept <- if (ncol(block) == 2L) ncol(block) else seq_len(ncol(block))
    do.call(rbind, lapply(kept, \(index) {
      data.frame(indicator = indicator, category_label = colnames(block)[index],
                 profile = seq_len(nrow(block)),
                 probability = block[, index], stringsAsFactors = FALSE)
    }))
  }))
  cells$indicator <- factor(cells$indicator, levels = names(blocks))
  cells$profile_label <- factor(as.character(key$label[cells$profile]),
                                levels = rev(levels(key$label)))
  cells$text_colour <- .gg_ink(.gg_probability_fill(cells$probability))
  shown <- if (binary) sprintf(" (category %s shown)",
                               paste(unique(cells$category_label),
                                     collapse = "/")) else ""
  # A binary indicator is one column named by the indicator; a polytomous one
  # is a block of columns named by category, under the indicator's strip.
  cells$column_label <- if (binary) cells$indicator else
    factor(cells$category_label, levels = unique(cells$category_label))
  plot <- ggplot2::ggplot(cells, ggplot2::aes(column_label, profile_label,
                                              fill = probability)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", probability),
                                    colour = text_colour), size = 3.1) +
    ggplot2::scale_colour_identity() +
    .gg_probability_scale("Probability") +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    .gg_titles(main, subtitle, "Response probabilities",
               sprintf("Probability of each response in each profile%s",
                       shown)) +
    ggplot2::labs(x = NULL, y = NULL) +
    .gg_theme() +
    .gg_tile_theme()
  if (binary) {
    return(plot + ggplot2::scale_x_discrete(expand = c(0, 0)) +
             ggplot2::theme(axis.text.x = ggplot2::element_text(
               angle = 30, hjust = 1)))
  }
  plot +
    ggplot2::facet_grid(cols = ggplot2::vars(indicator), scales = "free_x",
                        space = "free_x") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   panel.spacing = ggplot2::unit(0.6, "lines"),
                   strip.text = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5))
}

#' Each indicator's distribution per profile: density, quartiles, cases
#' @noRd
.gg_view_raincloud <- function(x, scale, main, subtitle) {
  vars <- .multilpa_continuous_names(x)
  if (length(vars) == 0L) {
    stop(errorCondition(
      "This model has no continuous indicators; plot `what = \"heatmap\"` instead.",
      class = "latents_no_continuous", call = NULL))
  }
  key <- .gg_profile_key(x)
  observed <- .multilpa_plot_observations(x, vars)
  values <- observed$values
  if (identical(scale, "standardized")) {
    basis <- .multilpa_standardization(x, vars)
    values <- sweep(sweep(values, 2L, basis$centre, "-"), 2L, basis$spread,
                    "/")
  }
  k <- nrow(key)
  y0 <- (k + 1L - key$rank)
  pieces <- lapply(seq_along(vars), \(index) {
    by_profile <- lapply(seq_len(k), \(profile) {
      own <- values[observed$profile == profile, index]
      own <- own[is.finite(own)]
      if (length(own) == 0L) return(NULL)
      label <- as.character(key$label[profile])
      density <- if (length(unique(own)) > 1L) {
        d <- stats::density(own, n = 256L, cut = 3)
        data.frame(x = d$x, height = d$y, y0 = y0[profile],
                   profile_label = label, indicator = vars[index])
      }
      quartiles <- stats::quantile(own, c(0.25, 0.5, 0.75), names = FALSE)
      ordering <- rank(own, ties.method = "first")
      list(
        density = density,
        box = data.frame(q1 = quartiles[1L], median = quartiles[2L],
                         q3 = quartiles[3L], y0 = y0[profile],
                         profile_label = label, indicator = vars[index]),
        points = data.frame(
          x = own,
          y = y0[profile] - 0.2 + 0.07 * ((ordering * 0.618) %% 1 - 0.5),
          profile_label = label, indicator = vars[index]))
    })
    by_profile <- Filter(Negate(is.null), by_profile)
    density <- do.call(rbind, lapply(by_profile, `[[`, "density"))
    # One height scale per indicator, so the ridges within a panel compare.
    if (!is.null(density)) {
      density$height <- 0.62 * density$height / max(density$height)
    }
    list(density = density,
         box = do.call(rbind, lapply(by_profile, `[[`, "box")),
         points = do.call(rbind, lapply(by_profile, `[[`, "points")))
  })
  bind <- \(part) {
    out <- do.call(rbind, lapply(pieces, `[[`, part))
    if (is.null(out)) return(NULL)
    out$indicator <- factor(out$indicator, levels = vars)
    out$profile_label <- factor(out$profile_label, levels = levels(key$label))
    out
  }
  density <- bind("density")
  box <- bind("box")
  points <- bind("points")
  plot <- ggplot2::ggplot()
  if (!is.null(density)) {
    plot <- plot +
      ggplot2::geom_ribbon(data = density,
                           ggplot2::aes(x = x, ymin = y0, ymax = y0 + height,
                                        fill = profile_label,
                                        group = profile_label),
                           alpha = 0.55, colour = NA) +
      ggplot2::geom_line(data = density,
                         ggplot2::aes(x = x, y = y0 + height,
                                      colour = profile_label,
                                      group = profile_label),
                         linewidth = 0.45)
  }
  plot +
    ggplot2::geom_segment(data = box,
                          ggplot2::aes(x = q1, xend = q3, y = y0 - 0.07,
                                       yend = y0 - 0.07,
                                       colour = profile_label),
                          linewidth = 1.6) +
    ggplot2::geom_point(data = box, ggplot2::aes(x = median, y = y0 - 0.07),
                        shape = 21, fill = "white", colour = "grey15",
                        size = 1.6, stroke = 0.5) +
    ggplot2::geom_point(data = points,
                        ggplot2::aes(x = x, y = y, colour = profile_label),
                        size = 0.7, alpha = 0.55) +
    ggplot2::facet_wrap(ggplot2::vars(indicator), scales = "free_x") +
    ggplot2::scale_y_continuous(breaks = y0, labels = as.character(key$label),
                                expand = ggplot2::expansion(add = c(0.35, 0.1))) +
    .gg_scales(key, c("colour", "fill")) +
    .gg_titles(main, subtitle, "Indicator distributions by profile",
               paste("Density, interquartile range with median, and each",
                     "case, grouped by its most likely profile")) +
    ggplot2::labs(x = if (identical(scale, "standardized"))
      "Standardized value" else NULL, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.3),
                   panel.spacing = ggplot2::unit(1.2, "lines"))
}

#' Categorical response probabilities, one line per profile
#' @noRd
.gg_view_responses <- function(x, category, labels, errors, main, subtitle) {
  blocks <- x$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) {
    stop(errorCondition("This model has no categorical indicators.",
                        class = "latents_no_categorical", call = NULL))
  }
  key <- .gg_profile_key(x)
  vars <- names(blocks)
  n_profiles <- x$n_profiles
  chosen <- vapply(blocks, \(block) {
    index <- if (identical(category, "last")) ncol(block) else
      if (identical(category, "first")) 1L else
        if (is.numeric(category)) as.integer(category) else
          match(as.character(category), colnames(block))
    if (is.na(index) || index < 1L || index > ncol(block)) {
      stop(errorCondition(sprintf(
        "`category` does not match a category of an indicator with categories %s.",
        paste(colnames(block), collapse = ", ")),
        class = "latents_unknown_category", call = NULL))
    }
    index
  }, integer(1))
  category_labels <- vapply(seq_along(blocks), \(index) {
    colnames(blocks[[index]])[chosen[[index]]]
  }, character(1))
  reason <- if (is.character(errors)) errors else NULL
  if (is.character(errors)) errors <- NULL
  cells <- expand.grid(position = seq_along(vars), profile = seq_len(n_profiles))
  frame <- data.frame(
    profile = cells$profile, position = cells$position,
    mean = vapply(seq_len(nrow(cells)), \(row) {
      blocks[[cells$position[row]]][cells$profile[row],
                                    chosen[[cells$position[row]]]]
    }, numeric(1)))
  spread <- if (is.null(errors)) rep(NA_real_, nrow(frame)) else
    .multilpa_match_error(errors, "response", frame$profile,
                          sprintf("%s:%s", vars[frame$position],
                                  category_labels[frame$position]))
  # A probability's interval is clipped to [0, 1], where it lives.
  frame$lower <- pmax(0, frame$mean - stats::qnorm(0.975) * spread)
  frame$upper <- pmin(1, frame$mean + stats::qnorm(0.975) * spread)
  frame$profile_label <- key$label[frame$profile]
  frame$rank <- key$rank[frame$profile]
  has_errors <- any(is.finite(spread))
  held <- isTRUE(attr(errors, "held"))
  note <- if (has_errors) {
    if (held) "; whiskers: 95% intervals where estimable" else
      "; whiskers: 95% intervals"
  } else if (!is.null(reason)) paste0("; ", reason) else ""
  # The axis names each indicator with its category when they differ, so a
  # mixed set cannot be misread as one shared category.
  shared <- length(unique(category_labels)) == 1L
  axis_labels <- if (shared) vars else sprintf("%s=%s", vars, category_labels)
  .gg_series_plot(frame, key, labels, has_errors,
                  y_label = if (shared)
                    sprintf("P(response = %s)", category_labels[[1L]]) else
                      "Response probability",
                  x_labels = axis_labels, ylim = c(0, 1)) +
    .gg_titles(main, subtitle, "Response probabilities by profile",
               sprintf("%d profiles%s, %d categorical indicator%s%s",
                       n_profiles, .gg_group_class_text(x), length(vars),
                       if (length(vars) == 1L) "" else "s", note))
}

# -- structure ----------------------------------------------------------------

#' Profile prevalence within each group class, one line per class
#' @noRd
.gg_view_probabilities <- function(x, labels, main, subtitle) {
  probabilities <- x$profile_probabilities
  if (is.null(probabilities)) {
    stop(errorCondition("This model has no profile prevalence to draw.",
                        class = "latents_nothing_to_plot", call = NULL))
  }
  key <- .gg_profile_key(x)
  n_classes <- nrow(probabilities)
  profile_order <- order(key$rank)
  cells <- expand.grid(position = seq_len(x$n_profiles),
                       class = seq_len(n_classes))
  frame <- data.frame(
    position = cells$position,
    mean = probabilities[cbind(cells$class, profile_order[cells$position])],
    class = cells$class)
  class_labels <- sprintf("Class %d (%s)", seq_len(n_classes),
                          .gg_percent(x$group_probabilities))
  frame$profile_label <- factor(class_labels[frame$class],
                                levels = class_labels)
  frame$rank <- frame$class
  class_key <- data.frame(
    label = factor(class_labels, levels = class_labels),
    colour = rep_len(.gg_okabe_ito, n_classes),
    shape = rep_len(.gg_shapes, n_classes), rank = seq_len(n_classes))
  frame$lower <- NA_real_
  frame$upper <- NA_real_
  .gg_series_plot(frame, class_key, labels, FALSE,
                  y_label = "Probability within group class",
                  x_labels = as.character(key$label[profile_order]),
                  ylim = c(0, 1)) +
    .gg_titles(main, subtitle, "Profile prevalence by group class",
               sprintf("%d group class%s over %d observed groups", n_classes,
                       if (n_classes == 1L) "" else "es", x$n_groups))
}

#' Each group's assigned profile at each position, blocked by group class
#' @noRd
.gg_view_sequences <- function(x, labels, cell_labels, main, subtitle) {
  layout <- .multilpa_sequences(x, format = "wide")
  ordering <- order(layout$group_class, layout$group)
  codes <- .multilpa_sequence_matrix(x)[ordering, , drop = FALSE]
  classes <- layout$group_class[ordering]
  key <- .gg_profile_key(x)
  time_values <- sort(unique(x$time_values))
  n_groups <- nrow(codes)
  cells <- expand.grid(row = seq_len(n_groups), column = seq_len(ncol(codes)))
  cells$profile <- codes[cbind(cells$row, cells$column)]
  cells <- cells[!is.na(cells$profile), , drop = FALSE]
  cells$class_label <- factor(sprintf("Class %d", classes[cells$row]),
                              levels = sprintf("Class %d", sort(unique(classes))))
  cells$profile_label <- factor(as.character(key$label[cells$profile]),
                                levels = levels(key$label))
  cells$cell_text <- as.character(cells$profile)
  cells$text_colour <- .gg_ink(key$colour[cells$profile])
  # A digit fits a cell only while the grid stays sparse enough to hold one.
  stamp <- isTRUE(cell_labels) && n_groups <= 60L && ncol(codes) <= 40L
  plot <- ggplot2::ggplot(cells, ggplot2::aes(column, row,
                                              fill = profile_label)) +
    # Cell borders only while rows are tall enough for a border to separate
    # them; on a dense grid they would cover the colour.
    ggplot2::geom_tile(colour = if (n_groups <= 40L) "white" else NA,
                       linewidth = 0.3)
  if (stamp) {
    plot <- plot + ggplot2::geom_text(
      ggplot2::aes(label = cell_text, colour = text_colour), size = 2.6) +
      ggplot2::scale_colour_identity()
  }
  if (length(unique(classes)) > 1L && isTRUE(labels)) {
    plot <- plot + ggplot2::facet_grid(rows = ggplot2::vars(class_label),
                                       scales = "free_y", space = "free_y",
                                       switch = "y")
  }
  plot +
    ggplot2::scale_x_continuous(breaks = seq_len(ncol(codes)),
                                labels = as.character(time_values),
                                expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    .gg_scales(key, "fill", legend = TRUE) +
    .gg_titles(main, subtitle, "Assigned profile in sequence order",
               sprintf("%d groups in %d classes; %d profiles across %d positions",
                       n_groups, x$n_group_classes, x$n_profiles,
                       ncol(codes))) +
    ggplot2::labs(x = x$time, y = sprintf("Groups (%d)", n_groups)) +
    .gg_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0,
                                                             face = "bold"),
                   panel.spacing = ggplot2::unit(0.4, "lines"),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
}

#' The estimated transition matrix, one panel per group class
#'
#' A row with no data support is labelled in parentheses: it is uniform by
#' construction rather than estimated.
#' @noRd
.gg_view_transitions <- function(x, main, subtitle) {
  probabilities <- x$transition_probabilities
  n_profiles <- x$n_profiles
  n_classes <- x$n_group_classes
  empty <- x$empty_transition_rows
  cells <- expand.grid(to = seq_len(n_profiles), from = seq_len(n_profiles),
                       class = seq_len(n_classes))
  cells$probability <- probabilities[cbind(cells$from, cells$to, cells$class)]
  unestimated <- if (is.null(empty)) rep(FALSE, nrow(cells)) else
    as.logical(empty[cbind(cells$from, cells$class)])
  cells$row_label <- ifelse(unestimated,
                            sprintf("(Profile %d)", cells$from),
                            sprintf("Profile %d", cells$from))
  row_levels <- unique(cells$row_label[order(-cells$from)])
  cells$row_label <- factor(cells$row_label, levels = row_levels)
  cells$column_label <- factor(sprintf("%d", cells$to),
                               levels = sprintf("%d", seq_len(n_profiles)))
  cells$class_label <- factor(
    sprintf("Group class %d (%s of groups)", cells$class,
            .gg_percent(x$group_probabilities[cells$class])),
    levels = sprintf("Group class %d (%s of groups)", seq_len(n_classes),
                     .gg_percent(x$group_probabilities)))
  cells$text_colour <- .gg_ink(.gg_probability_fill(cells$probability))
  ggplot2::ggplot(cells, ggplot2::aes(column_label, row_label,
                                      fill = probability)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", probability),
                                    colour = text_colour), size = 3.4) +
    ggplot2::scale_colour_identity() +
    .gg_probability_scale("Probability") +
    ggplot2::facet_wrap(ggplot2::vars(class_label), scales = "free_y") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    .gg_titles(main, subtitle, "Estimated transitions",
               paste("Row: current profile; column: next profile;",
                     "each row sums to one")) +
    ggplot2::labs(x = "Next profile", y = NULL) +
    .gg_theme() +
    .gg_tile_theme()
}

# -- classification -----------------------------------------------------------

#' Effective profile sizes, largest first, labelled with count and share
#' @noRd
.gg_view_sizes <- function(x, main, subtitle) {
  if (is.null(x$effective_profile_counts) ||
      length(x$effective_profile_counts) == 0L) {
    stop(errorCondition(
      "This model carries no effective profile counts to draw.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  key <- .gg_profile_key(x)
  key$profile_label <- factor(as.character(key$label),
                              levels = rev(levels(key$label)))
  key$bar_label <- sprintf("%.1f  (%.1f%%)", key$count, 100 * key$share)
  total <- sum(x$effective_profile_counts)
  ggplot2::ggplot(key, ggplot2::aes(count, profile_label, fill = label)) +
    ggplot2::geom_col(width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = bar_label), hjust = -0.08,
                       size = 3.4, colour = "grey20") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.24))) +
    .gg_scales(key, "fill") +
    .gg_titles(main, subtitle, "Profile sizes",
               sprintf(paste("Effective number of cases (sum of posterior",
                             "probabilities); %s in total"),
                       format(round(total), big.mark = ","))) +
    ggplot2::labs(x = "Effective cases", y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey90",
                                                              linewidth = 0.35))
}

#' Per-case certainty and relative entropy from the fit's own posteriors
#' @noRd
.gg_case_values <- function(x) {
  posteriors <- x$subject_posteriors
  k <- ncol(posteriors)
  assigned <- max.col(posteriors, ties.method = "first")
  # -sum p log p per case, relative to its maximum log(K): 0 is certain,
  # 1 is a flat posterior.
  entropy <- -rowSums(posteriors * log(pmax(posteriors, .Machine$double.xmin)))
  data.frame(profile = assigned,
             certainty = posteriors[cbind(seq_len(nrow(posteriors)), assigned)],
             entropy = if (k > 1L) entropy / log(k) else 0)
}

#' Strip plot of one per-case value by profile, with each profile's mean
#' @noRd
.gg_case_strip <- function(x, value, limits, reference, reference_label,
                           default_main, default_subtitle, x_label, main,
                           subtitle) {
  key <- .gg_profile_key(x)
  cases <- .gg_case_values(x)
  k <- nrow(key)
  y0 <- k + 1L - key$rank
  cases$value <- cases[[value]]
  cases$y0 <- y0[cases$profile]
  cases$profile_label <- key$label[cases$profile]
  within <- stats::ave(cases$value, cases$profile,
                       FUN = \(v) rank(v, ties.method = "first"))
  cases$y <- cases$y0 + 0.5 * ((within * 0.618) %% 1 - 0.5)
  means <- stats::aggregate(value ~ profile + y0, data = cases, FUN = mean)
  means$label <- sprintf("mean %.2f", means$value)
  at_right <- isTRUE(all.equal(reference, limits[2L]))
  ggplot2::ggplot(cases, ggplot2::aes(value, y, colour = profile_label)) +
    ggplot2::geom_vline(xintercept = reference, linetype = "22",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::annotate("text", x = reference, y = k + 0.45,
                      label = reference_label, vjust = 1,
                      hjust = if (at_right) 1.08 else -0.08,
                      size = 3, colour = "grey45") +
    ggplot2::geom_point(size = 1.1, alpha = 0.55) +
    ggplot2::geom_segment(data = means,
                          ggplot2::aes(x = value, xend = value,
                                       y = y0 - 0.32, yend = y0 + 0.32),
                          colour = "grey10", linewidth = 0.8,
                          inherit.aes = FALSE) +
    ggplot2::geom_text(data = means,
                       ggplot2::aes(x = value, y = y0 + 0.34, label = label),
                       colour = "grey10", size = 3, vjust = 0,
                       hjust = if (identical(value, "certainty")) 1.05 else
                         -0.05,
                       inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(limits = limits,
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = y0, labels = as.character(key$label),
                                limits = c(0.5, k + 0.5),
                                expand = ggplot2::expansion(0)) +
    .gg_scales(key, "colour") +
    ggplot2::coord_cartesian(clip = "off") +
    .gg_titles(main, subtitle, default_main, default_subtitle) +
    ggplot2::labs(x = x_label, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.3))
}

#' @noRd
.gg_view_posteriors <- function(x, main, subtitle) {
  .multilpa_require_case_posteriors(x, "posteriors")
  k <- ncol(x$subject_posteriors)
  .gg_case_strip(
    x, value = "certainty", limits = c(1 / k, 1),
    reference = 1 / k, reference_label = "chance",
    default_main = "Classification certainty",
    default_subtitle = paste("Posterior probability of each case's assigned",
                             "profile; the profile mean is its average",
                             "posterior probability"),
    x_label = "Posterior probability of the assigned profile",
    main = main, subtitle = subtitle)
}

#' @noRd
.gg_view_entropy <- function(x, main, subtitle) {
  .multilpa_require_case_posteriors(x, "entropy")
  cases <- .gg_case_values(x)
  .gg_case_strip(
    x, value = "entropy", limits = c(0, 1),
    reference = 1, reference_label = "flat",
    default_main = "Classification uncertainty",
    default_subtitle = sprintf(paste(
      "Each case's entropy relative to a flat posterior (0 = certain,",
      "1 = flat); relative entropy = 1 - mean = %.3f"),
      1 - mean(cases$entropy)),
    x_label = "Relative case entropy", main = main, subtitle = subtitle)
}

#' Average posterior probability matrix, assigned by posterior profile
#' @noRd
.gg_view_avepp <- function(x, main, subtitle) {
  posteriors <- x$subject_posteriors
  if (is.null(posteriors) || ncol(posteriors) < 1L) {
    stop(errorCondition(
      "This model carries no individual posteriors to average.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  key <- .gg_profile_key(x)
  averages <- .multilpa_average_posterior_matrix(posteriors)
  k <- ncol(averages)
  assigned <- tabulate(max.col(posteriors, ties.method = "first"), k)
  display <- order(key$rank)
  cells <- expand.grid(column = seq_len(k), row = seq_len(k))
  cells$average_posterior <- averages[cbind(cells$row, cells$column)]
  assigned_levels <- sprintf("Assigned %d (n = %d)", seq_len(k), assigned)
  cells$assigned_label <- factor(assigned_levels[cells$row],
                                 levels = rev(assigned_levels[display]))
  cells$posterior_label <- factor(sprintf("Profile %d", cells$column),
                                  levels = sprintf("Profile %d", display))
  cells$text_colour <- .gg_ink(.gg_probability_fill(cells$average_posterior))
  cells$label <- ifelse(is.na(cells$average_posterior), "",
                        sprintf("%.2f", cells$average_posterior))
  diagonal <- cells[cells$row == cells$column, , drop = FALSE]
  ggplot2::ggplot(cells, ggplot2::aes(posterior_label, assigned_label,
                                      fill = average_posterior)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_tile(data = diagonal, fill = NA, colour = "grey15",
                       linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = text_colour),
                       size = 3.4) +
    ggplot2::scale_colour_identity() +
    .gg_probability_scale("Average\nposterior") +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    .gg_titles(main, subtitle, "Average posterior probability",
               paste("Mean posterior each assigned group places on each",
                     "profile; the outlined diagonal is classification",
                     "accuracy")) +
    ggplot2::labs(x = NULL, y = NULL) +
    .gg_theme() +
    .gg_tile_theme()
}

# -- cases --------------------------------------------------------------------

#' The observations with their profile and assignment certainty
#' @noRd
.gg_observed_cases <- function(x, vars) {
  observed <- .multilpa_plot_observations(x, vars)
  kept <- x$subject_profiles >= 1L & x$subject_profiles <= x$n_profiles
  posteriors <- x$subject_posteriors[kept, , drop = FALSE]
  list(values = observed$values, profile = observed$profile,
       certainty = posteriors[cbind(seq_len(nrow(posteriors)),
                                    observed$profile)])
}

#' Every case as a line across indicators, one panel per profile
#'
#' Indicators are put on one scale, the observed z-score of each, so lines
#' across them mean something whatever the units.
#' @noRd
.gg_view_parallel <- function(x, main, subtitle) {
  vars <- .multilpa_continuous_names(x)
  if (length(vars) < 2L) {
    stop(errorCondition(
      "Parallel coordinates need at least two continuous indicators.",
      class = "latents_no_continuous", call = NULL))
  }
  key <- .gg_profile_key(x)
  cases <- .gg_observed_cases(x, vars)
  basis <- .multilpa_standardization(x, vars)
  z <- sweep(sweep(cases$values, 2L, basis$centre, "-"), 2L, basis$spread, "/")
  long <- data.frame(
    case = rep(seq_len(nrow(z)), times = length(vars)),
    position = rep(seq_along(vars), each = nrow(z)),
    value = as.vector(z),
    profile = rep(cases$profile, times = length(vars)),
    certainty = rep(cases$certainty, times = length(vars)))
  long <- long[is.finite(long$value), , drop = FALSE]
  long$panel <- key$label[long$profile]
  context <- do.call(rbind, lapply(levels(key$label), \(panel) {
    cbind(long[c("case", "position", "value")], panel = panel)
  }))
  context$panel <- factor(context$panel, levels = levels(key$label))
  means <- .gg_mean_frame(x, key, "standardized", NULL)$frame
  means$panel <- means$profile_label
  ggplot2::ggplot() +
    ggplot2::geom_line(data = context,
                       ggplot2::aes(position, value, group = case),
                       colour = "grey88", linewidth = 0.25) +
    ggplot2::geom_line(data = long,
                       ggplot2::aes(position, value, group = case,
                                    colour = panel, alpha = certainty),
                       linewidth = 0.35) +
    ggplot2::geom_line(data = means,
                       ggplot2::aes(position, mean, group = panel),
                       colour = "grey10", linewidth = 1.1) +
    ggplot2::geom_point(data = means,
                        ggplot2::aes(position, mean, fill = panel,
                                     shape = panel),
                        colour = "grey10", size = 2.4, stroke = 0.6) +
    ggplot2::facet_wrap(ggplot2::vars(panel)) +
    ggplot2::scale_x_continuous(breaks = seq_along(vars), labels = vars,
                                expand = ggplot2::expansion(add = 0.15)) +
    ggplot2::scale_alpha_continuous(range = c(0.12, 0.7), limits = c(0, 1),
                                    guide = "none") +
    .gg_scales(key) +
    .gg_titles(main, subtitle, "Cases across indicators",
               paste("Each line is a case, faded by the certainty of its",
                     "assignment; the dark line is the profile mean;",
                     "other cases in grey")) +
    ggplot2::labs(x = NULL, y = "Standardized value (z)") +
    .gg_theme() +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.2, "lines"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Points of a 95% ellipse of a bivariate normal
#' @noRd
.gg_ellipse <- function(centre, covariance, level = 0.95, n = 90L) {
  decomposition <- eigen(covariance, symmetric = TRUE)
  radius <- sqrt(stats::qchisq(level, df = 2L))
  angle <- seq(0, 2 * pi, length.out = n)
  axes <- decomposition$vectors %*%
    diag(sqrt(pmax(decomposition$values, 0)), 2L)
  circle <- rbind(cos(angle), sin(angle)) * radius
  points <- t(axes %*% circle) + rep(centre, each = n)
  data.frame(x = points[, 1L], y = points[, 2L])
}

#' One profile's covariance for a pair of indicators
#'
#' A diagonal structure stores variances only, so the pair's covariance is
#' zero there.
#' @noRd
.gg_pair_covariance <- function(x, profile, i, j) {
  if (!is.null(x$covariances)) {
    return(x$covariances[c(i, j), c(i, j), profile])
  }
  diag(x$variances[profile, c(i, j)], 2L)
}

#' Scatter-plot matrix: cases by profile with each profile's 95% ellipse
#' @noRd
.gg_view_pairs <- function(x, main, subtitle) {
  vars <- .multilpa_continuous_names(x)
  d <- length(vars)
  if (d < 2L) {
    stop(errorCondition(
      "A scatter-plot matrix needs at least two continuous indicators.",
      class = "latents_no_continuous", call = NULL))
  }
  key <- .gg_profile_key(x)
  cases <- .gg_observed_cases(x, vars)
  labels <- as.character(key$label[cases$profile])
  pairs <- expand.grid(i = seq_len(d), j = seq_len(d))
  pairs <- pairs[pairs$i < pairs$j, , drop = FALSE]
  scatter <- do.call(rbind, Map(\(i, j) {
    data.frame(x = cases$values[, i], y = cases$values[, j],
               uncertainty = 1 - cases$certainty, profile_label = labels,
               xvar = vars[i], yvar = vars[j])
  }, pairs$i, pairs$j))
  scatter <- scatter[is.finite(scatter$x) & is.finite(scatter$y), ,
                     drop = FALSE]
  ellipses <- do.call(rbind, Map(\(i, j) {
    do.call(rbind, lapply(seq_len(x$n_profiles), \(profile) {
      ellipse <- .gg_ellipse(x$means[profile, c(i, j)],
                             .gg_pair_covariance(x, profile, i, j))
      cbind(ellipse, profile_label = as.character(key$label[profile]),
            xvar = vars[i], yvar = vars[j])
    }))
  }, pairs$i, pairs$j))
  diagonal <- do.call(rbind, lapply(seq_len(d), \(index) {
    values <- cases$values[, index]
    range_v <- range(values, na.rm = TRUE)
    curves <- do.call(rbind, lapply(seq_len(x$n_profiles), \(profile) {
      own <- values[cases$profile == profile & is.finite(values)]
      if (length(own) < 3L || !(stats::sd(own) > 0)) return(NULL)
      dens <- stats::density(own, n = 200L, from = range_v[1L],
                             to = range_v[2L])
      data.frame(x = dens$x, height = dens$y * key$share[profile],
                 profile_label = as.character(key$label[profile]))
    }))
    if (is.null(curves)) return(NULL)
    # Size-weighted densities, drawn into the indicator's own range.
    curves$y <- range_v[1L] + 0.92 * diff(range_v) *
      curves$height / max(curves$height)
    curves$ymin <- range_v[1L]
    cbind(curves, xvar = vars[index], yvar = vars[index])
  }))
  as_levels <- \(frame) {
    if (is.null(frame)) return(NULL)
    frame$xvar <- factor(frame$xvar, levels = vars)
    frame$yvar <- factor(frame$yvar, levels = vars)
    if ("profile_label" %in% names(frame)) {
      frame$profile_label <- factor(frame$profile_label,
                                    levels = levels(key$label))
    }
    frame
  }
  scatter <- as_levels(scatter)
  ellipses <- as_levels(ellipses)
  diagonal <- as_levels(diagonal)
  panels <- as_levels(unique(rbind(
    scatter[c("xvar", "yvar")],
    data.frame(xvar = vars, yvar = vars))))
  plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = panels, xmin = -Inf, xmax = Inf, ymin = -Inf,
                       ymax = Inf, fill = "grey97", colour = NA)
  if (!is.null(diagonal)) {
    plot <- plot +
      ggplot2::geom_ribbon(data = diagonal,
                           ggplot2::aes(x = x, ymin = ymin, ymax = y,
                                        fill = profile_label,
                                        group = profile_label),
                           alpha = 0.35, colour = NA) +
      ggplot2::geom_line(data = diagonal,
                         ggplot2::aes(x = x, y = y, colour = profile_label,
                                      group = profile_label),
                         linewidth = 0.5)
  }
  plot +
    ggplot2::geom_point(data = scatter,
                        ggplot2::aes(x = x, y = y, colour = profile_label,
                                     size = uncertainty),
                        alpha = 0.7, stroke = 0) +
    ggplot2::geom_path(data = ellipses,
                       ggplot2::aes(x = x, y = y, colour = profile_label,
                                    group = profile_label),
                       linewidth = 0.7) +
    ggplot2::facet_grid(rows = ggplot2::vars(yvar), cols = ggplot2::vars(xvar),
                        scales = "free", switch = "both") +
    ggplot2::scale_size_continuous(range = c(0.9, 3.4), limits = c(0, 1),
                                   guide = "none") +
    .gg_scales(key, "colour", legend = TRUE) +
    .gg_scales(key, "fill") +
    .gg_titles(main, subtitle, "Profiles in each pair of indicators",
               paste0("Cases coloured by profile, larger when less certain; ",
                      "95% ellipses of the fitted covariances\n",
                      "Diagonal: each profile's density, weighted by its size")) +
    ggplot2::labs(x = NULL, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text = ggplot2::element_text(face = "bold",
                                                      size = ggplot2::rel(0.85)),
                   strip.text.y.left = ggplot2::element_text(angle = 90),
                   panel.spacing = ggplot2::unit(0.5, "lines"),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
}

# -- candidate grids ----------------------------------------------------------

#' One information criterion's series across the grid, as a long frame
#'
#' A series is one covariance model at one number of group classes; its label
#' names only what tells the series apart.
#' @noRd
.gg_enumeration_frame <- function(grid, criteria) {
  model <- grid$model
  known <- unique(stats::na.omit(model))
  if (length(known) == 1L) model[is.na(model)] <- known
  model[is.na(model)] <- "-"
  classes_text <- sprintf("%d group class%s", grid$n_group_classes,
                          ifelse(grid$n_group_classes == 1L, "", "es"))
  several_models <- length(unique(model)) > 1L
  several_groups <- length(unique(grid$n_group_classes)) > 1L
  series <- if (several_models && several_groups) {
    paste(model, classes_text, sep = ", ")
  } else if (several_models) model else if (several_groups) classes_text else
    rep("", nrow(grid))
  do.call(rbind, Map(\(column, title) {
    values <- grid[[column]]
    eligible <- grid$converged %in% TRUE & is.finite(values)
    data.frame(title = title, series = series, n_profiles = grid$n_profiles,
               value = values, eligible = eligible)
  }, criteria, names(criteria)))
}

#' Information criteria against the number of profiles
#' @noRd
.gg_view_enumeration <- function(x, criteria, labels, mark_minimum, main,
                                 subtitle) {
  grid <- as.data.frame(x)
  long <- .gg_enumeration_frame(grid, criteria)
  long$title <- factor(long$title, levels = names(criteria))
  series_levels <- unique(long$series[order(grid$n_group_classes, grid$model)])
  long$series <- factor(long$series, levels = series_levels)
  ok <- long[long$eligible, , drop = FALSE]
  # A failed candidate is marked on its panel's floor, so a gap in a line reads
  # as a failure and not as a candidate that was never fitted.
  floors <- tapply(ok$value, ok$title, \(v) min(v) - 0.08 * max(diff(range(v)),
                                                                1e-8))
  failed <- long[!long$eligible, , drop = FALSE]
  failed$failed_y <- floors[as.character(failed$title)]
  minimum <- do.call(rbind, lapply(split(ok, ok$title, drop = TRUE),
                                   \(panel) panel[which.min(panel$value), ]))
  ends <- do.call(rbind, lapply(split(ok, list(ok$title, ok$series),
                                      drop = TRUE),
                                \(s) s[which.max(s$n_profiles), ]))
  ends <- do.call(rbind, lapply(split(ends, ends$title, drop = TRUE), \(p) {
    span <- diff(range(ok$value[ok$title == p$title[1L]]))
    p$y <- .gg_spread(p$value, gap = 0.07 * max(span, 1e-8))
    p
  }))
  ends$x <- max(ok$n_profiles) + 0.28
  n_series <- nlevels(long$series)
  colours <- stats::setNames(rep_len(.gg_okabe_ito, n_series), series_levels)
  shapes <- stats::setNames(rep_len(.gg_shapes, n_series), series_levels)
  # Direct labels while they fit beside the lines; past five series they
  # collide, and a legend keeps every series identifiable.
  labelled <- isTRUE(labels) && n_series > 1L && n_series <= 5L
  legend <- n_series > 1L && !labelled
  guide <- if (legend) "legend" else "none"
  failures <- sum(!grid$converged)
  plot <- ggplot2::ggplot(ok, ggplot2::aes(n_profiles, value, colour = series,
                                           group = series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = series, fill = series),
                        size = 2.3, colour = "white", stroke = 0.6)
  if (nrow(failed) > 0L) {
    plot <- plot + ggplot2::geom_point(
      data = failed, ggplot2::aes(x = n_profiles, y = failed_y,
                                  colour = series),
      shape = 4, size = 2.6, stroke = 1.1, inherit.aes = FALSE)
  }
  if (isTRUE(mark_minimum)) {
    plot <- plot + ggplot2::geom_point(
      data = minimum, ggplot2::aes(x = n_profiles, y = value), shape = 21,
      size = 5.2, stroke = 0.8, colour = "grey10", fill = NA,
      inherit.aes = FALSE)
  }
  if (labelled) {
    plot <- plot + ggplot2::geom_text(
      data = ends, ggplot2::aes(x = x, y = y, label = series, colour = series),
      hjust = 0, size = 3.2, fontface = "bold", inherit.aes = FALSE)
  }
  several <- length(criteria) > 1L
  plot +
    ggplot2::facet_wrap(ggplot2::vars(title), scales = "free_y",
                        nrow = if (length(criteria) <= 3L) 1L else NULL,
                        ncol = if (length(criteria) > 3L) 2L else NULL) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(grid$n_profiles)),
      expand = ggplot2::expansion(add = c(0.3, if (labelled) 1.1 else 0.3))) +
    ggplot2::scale_colour_manual(values = colours, guide = guide) +
    ggplot2::scale_fill_manual(values = colours, guide = guide) +
    ggplot2::scale_shape_manual(values = shapes, guide = guide) +
    ggplot2::coord_cartesian(clip = "off") +
    .gg_titles(main, subtitle,
               if (several) "Information criteria by candidate model" else
                 sprintf("%s by candidate model", names(criteria)),
               sprintf(paste("%d candidates%s; lower is preferred%s;",
                             "no candidate is selected automatically"),
                       nrow(grid),
                       if (failures == 0L) "" else
                         sprintf(", %d did not converge (x)", failures),
                       if (isTRUE(mark_minimum)) ", ring = lowest" else "")) +
    ggplot2::labs(x = "Number of profiles",
                  y = if (several) NULL else names(criteria)) +
    .gg_theme() +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.6, "lines"),
                   strip.text = if (several) ggplot2::element_text(
                     face = "bold", hjust = 0, size = ggplot2::rel(0.95)) else
                     ggplot2::element_blank(),
                   legend.position = if (legend) "bottom" else "none",
                   legend.title = ggplot2::element_blank())
}

#' Band between two x-intervals at two heights, eased so it leaves and
#' arrives vertically
#' @noRd
.gg_flow_band <- function(top_left, top_right, bottom_left, bottom_right,
                          y_top, y_bottom, n = 32L) {
  t <- seq(0, 1, length.out = n)
  ease <- t * t * (3 - 2 * t)
  y <- y_top + (y_bottom - y_top) * t
  left <- top_left + (bottom_left - top_left) * ease
  right <- top_right + (bottom_right - top_right) * ease
  data.frame(x = c(left, rev(right)), y = c(y, rev(y)))
}

#' How cases move between profiles as the number of profiles grows
#'
#' An icicle: each row is one solution, a bar split into its profiles by
#' posterior share. Bands between rows carry the posterior mass shared by a
#' profile and one in the next solution, the sum over cases of p(a) * p(b), so
#' each profile's bands add up exactly to its share. Colour is lineage: a
#' profile that is its dominant parent's main continuation keeps the parent's
#' colour, any other profile takes a new one. The layout follows clustering
#' trees (Zappia and Oshlack, 2018). One panel per series (covariance model
#' and number of group classes).
#' @noRd
.gg_view_tree <- function(x, main, subtitle) {
  grid <- as.data.frame(x)
  grid <- grid[grid$converged %in% TRUE, , drop = FALSE]
  model <- grid$model
  model[is.na(model)] <- ""
  keys <- unique(data.frame(model = model,
                            n_group_classes = grid$n_group_classes))
  several_groups <- length(unique(keys$n_group_classes)) > 1L
  half <- 0.16
  gap <- 0.03
  segments_of <- \(share, ordering) {
    widths <- share[ordering]
    k <- length(widths)
    right <- cumsum(widths) + gap * (seq_len(k) - 1) - gap * (k - 1) / 2
    data.frame(profile = ordering, xmin = right - widths, xmax = right,
               share = widths)
  }
  per_series <- lapply(seq_len(nrow(keys)), \(index) {
    series_model <- keys$model[index]
    classes <- keys$n_group_classes[index]
    rows_in <- model == series_model & grid$n_group_classes == classes
    counts <- sort(unique(grid$n_profiles[rows_in]))
    if (length(counts) < 2L) return(NULL)
    posteriors <- lapply(counts, \(k) {
      candidate_fit(x, n_profiles = k, n_group_classes = classes,
                    model = if (nzchar(series_model)) series_model)$subject_posteriors
    })
    n_cases <- nrow(posteriors[[1L]])
    flows <- lapply(seq_along(counts)[-1L], \(level) {
      crossprod(posteriors[[level - 1L]], posteriors[[level]]) / n_cases
    })
    first_share <- colSums(posteriors[[1L]]) / n_cases
    first_order <- order(-first_share, seq_along(first_share))
    # Order each row's profiles by the flow-weighted centre of their parents,
    # so bands cross only where cases really are reshuffled.
    rows <- Reduce(\(previous, level) {
      parent <- previous[[length(previous)]]
      flow <- flows[[level - 1L]]
      centres <- (parent$xmin + parent$xmax)[order(parent$profile)] / 2
      child_centre <- colSums(flow * centres) / colSums(flow)
      c(previous, list(segments_of(colSums(flow),
                                   order(child_centre,
                                         seq_along(child_centre)))))
    }, seq_along(counts)[-1L], init = list(segments_of(first_share,
                                                       first_order)))
    first_colours <- character(counts[1L])
    first_colours[first_order] <- rep_len(.gg_okabe_ito, counts[1L])
    lineage <- Reduce(\(state, level) {
      flow <- flows[[level - 1L]]
      dominant_parent <- apply(flow, 2L, which.max)
      main_child <- apply(flow, 1L, which.max)
      continues <- main_child[dominant_parent] == seq_len(ncol(flow))
      fresh <- rep_len(.gg_okabe_ito,
                       state$used + sum(!continues))[-seq_len(state$used)]
      colours <- character(ncol(flow))
      colours[continues] <- state$colours[[level - 1L]][
        dominant_parent[continues]]
      colours[!continues] <- fresh
      list(colours = c(state$colours, list(colours)),
           used = state$used + sum(!continues))
    }, seq_along(counts)[-1L],
    init = list(colours = list(first_colours), used = counts[1L]))$colours
    panel <- if (nzchar(series_model)) series_model else "Candidates"
    classes_label <- sprintf("%d group class%s", classes,
                             if (classes == 1L) "" else "es")
    bars <- do.call(rbind, Map(\(k, segments, colours) {
      cbind(segments, panel = panel, classes_label = classes_label,
            n_profiles = k, fill = colours[segments$profile])
    }, counts, rows, lineage))
    bands <- do.call(rbind, lapply(seq_along(counts)[-1L], \(level) {
      parent <- rows[[level - 1L]]
      child <- rows[[level]]
      flow <- flows[[level - 1L]]
      links <- expand.grid(from = parent$profile, to = child$profile)
      links$mass <- flow[cbind(links$from, links$to)]
      links <- links[links$mass > 1e-3, , drop = FALSE]
      # Outflows stack inside a parent in child order and inflows inside a
      # child in parent order, so bands tile both bars without overlap.
      links$from_rank <- match(links$from, parent$profile)
      links$to_rank <- match(links$to, child$profile)
      links <- links[order(links$from_rank, links$to_rank), , drop = FALSE]
      links$top_right <- parent$xmin[links$from_rank] +
        stats::ave(links$mass, links$from, FUN = cumsum)
      links$top_left <- links$top_right - links$mass
      links <- links[order(links$to_rank, links$from_rank), , drop = FALSE]
      links$bottom_right <- child$xmin[links$to_rank] +
        stats::ave(links$mass, links$to, FUN = cumsum)
      links$bottom_left <- links$bottom_right - links$mass
      do.call(rbind, lapply(seq_len(nrow(links)), \(i) {
        band <- .gg_flow_band(links$top_left[i], links$top_right[i],
                              links$bottom_left[i], links$bottom_right[i],
                              y_top = counts[level - 1L] + half,
                              y_bottom = counts[level] - half)
        cbind(band, panel = panel, classes_label = classes_label,
              group = sprintf("%s-%s-%d-%d-%d", panel, classes_label,
                              counts[level], links$from[i], links$to[i]),
              fill = lineage[[level - 1L]][links$from[i]])
      }))
    }))
    list(bars = bars, bands = bands)
  })
  per_series <- Filter(Negate(is.null), per_series)
  if (length(per_series) == 0L) {
    stop(errorCondition(
      "A tree needs at least two numbers of profiles for one model.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  bars <- do.call(rbind, lapply(per_series, `[[`, "bars"))
  bands <- do.call(rbind, lapply(per_series, `[[`, "bands"))
  levels <- unique(bars$panel)
  class_levels <- unique(bars$classes_label[order(bars$classes_label)])
  # A share is written only where its segment is wide enough to hold it; the
  # more structures side by side, the narrower every segment.
  min_share <- 0.02 + 0.04 * length(levels)
  bars$label <- ifelse(bars$share < min_share, "", .gg_percent(bars$share))
  bars$text_colour <- .gg_ink(bars$fill)
  bars$panel <- factor(bars$panel, levels = levels)
  bands$panel <- factor(bands$panel, levels = levels)
  bars$classes_label <- factor(bars$classes_label, levels = class_levels)
  bands$classes_label <- factor(bands$classes_label, levels = class_levels)
  counts <- sort(unique(bars$n_profiles))
  facets <- if (several_groups) {
    ggplot2::facet_grid(rows = ggplot2::vars(classes_label),
                        cols = ggplot2::vars(panel))
  } else ggplot2::facet_wrap(ggplot2::vars(panel), nrow = 1L)
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = bands,
                          ggplot2::aes(x = x, y = y, group = group,
                                       fill = fill),
                          alpha = 0.4, colour = NA) +
    ggplot2::geom_rect(data = bars,
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = n_profiles - half,
                                    ymax = n_profiles + half, fill = fill),
                       colour = "white", linewidth = 0.9) +
    ggplot2::geom_text(data = bars,
                       ggplot2::aes(x = (xmin + xmax) / 2, y = n_profiles,
                                    label = label, colour = text_colour),
                       size = 3, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    facets +
    ggplot2::scale_y_reverse(breaks = counts,
                             expand = ggplot2::expansion(add = 0.25)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    .gg_titles(main, subtitle, "How profiles split as more are added",
               paste0("Each row is one solution, divided into its profiles ",
                      "by share of cases\n",
                      "Colour follows a profile across solutions; a new ",
                      "colour marks a profile that split off"),
               caption = paste0("A band that forks: a profile split in two ",
                                "(nested, stable solutions). Bands that merge ",
                                "or cross: cases reshuffled between profiles.")) +
    ggplot2::labs(x = NULL, y = "Number of profiles") +
    .gg_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(1.5, "lines"))
}

# -- several plots ------------------------------------------------------------

#' A list of ggplot objects that draws every one when printed
#' @noRd
.gg_plots <- function(plots) {
  structure(plots, class = "latents_plots")
}

#' Print several plots
#'
#' The value of `plot(x, what = "all")`, of `plot()` on a diagnostics result
#' and of an enumeration plotted with `combine = FALSE`: a list of ggplot
#' objects, named by view. Printing draws each in turn.
#'
#' @param x A `latents_plots` list.
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the side effect of drawing.
#' @export
print.latents_plots <- function(x, ...) {
  invisible(lapply(x, print))
  invisible(x)
}
