# ggplot2 views of a fit (prototype, stages 1-2 of the plotting rebuild).
#
# Each builder returns a ggplot object and reads only the public tables of
# get_results(), so a view can never disagree with the numbers a user prints.
# Shared rules, applied by every view:
#   * one profile share: the posterior (effective) proportion from "counts";
#   * one display order: largest profile first, profile id breaking ties;
#   * one colour and one shape per profile, in that order (Okabe-Ito);
#   * direct labels where a legend would make the reader match colours.

utils::globalVariables(c(
  "x", "y", "ymin", "ymax", "xmin", "xmax", "lower", "upper", "mean",
  "label", "profile_label", "indicator", "value", "fill_value", "text_colour",
  "position", "group", "height", "criterion", "model", "n_profiles",
  "assigned_label", "posterior_label", "average_posterior", "q1", "q3",
  "median", "offset", "effective_count", "bar_label", "is_minimum", "case",
  "certainty", "panel", "uncertainty", "xvar", "yvar", "xend", "yend",
  "cases", "share", "fill"))

.gg_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                   "#CC79A7", "#999999", "#000000", "#F0E442")
.gg_shapes <- c(21L, 22L, 24L, 23L, 25L)

#' Stop unless ggplot2 is installed
#' @noRd
.gg_require <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(errorCondition(
      "These views need the ggplot2 package: install.packages(\"ggplot2\").",
      class = "latents_missing_package", call = NULL))
  }
  invisible(TRUE)
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

#' One row per profile: share, display order, label, colour and shape
#' @return A data.frame ordered for display (largest profile first).
#' @noRd
.gg_profile_key <- function(fit) {
  counts <- get_results(fit, what = "counts")
  counts <- counts[counts$level == "individuals", , drop = FALSE]
  counts <- counts[order(-counts$effective_proportion, counts$class), ,
                   drop = FALSE]
  n <- nrow(counts)
  share <- 100 * counts$effective_proportion
  share_text <- ifelse(share < 1, "<1%", sprintf("%.0f%%", share))
  label <- sprintf("Profile %d (%s)", counts$class, share_text)
  data.frame(
    profile = counts$class,
    rank = seq_len(n),
    share = counts$effective_proportion,
    effective_count = counts$effective_count,
    label = factor(label, levels = label),
    colour = rep_len(.gg_okabe_ito, n),
    shape = rep_len(.gg_shapes, n),
    stringsAsFactors = FALSE
  )
}

#' Manual profile scales for the aesthetics a view maps
#' @noRd
.gg_scales <- function(key, aesthetics = c("colour", "fill", "shape")) {
  colours <- stats::setNames(key$colour, levels(key$label))
  shapes <- stats::setNames(key$shape, levels(key$label))
  scales <- list(
    colour = ggplot2::scale_colour_manual(values = colours),
    fill = ggplot2::scale_fill_manual(values = colours),
    shape = ggplot2::scale_shape_manual(values = shapes)
  )
  scales[aesthetics]
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

#' Continuous indicator means with 95% Wald intervals, joined to the key
#' @noRd
.gg_profile_means <- function(fit, key) {
  means <- get_results(fit, what = "profiles")
  if (nrow(means) == 0L) {
    stop(errorCondition(
      "This view needs continuous indicators; the fit has none.",
      class = "latents_bad_argument", call = NULL))
  }
  indicators <- intersect(fit$vars, unique(means$indicator))
  z <- stats::qnorm(0.975)
  means$lower <- means$mean - z * means$mean_standard_error
  means$upper <- means$mean + z * means$mean_standard_error
  means$position <- match(means$indicator, indicators)
  means$indicator <- factor(means$indicator, levels = indicators)
  means$profile_label <- key$label[match(means$profile, key$profile)]
  means$rank <- key$rank[match(means$profile, key$profile)]
  means
}

# -- views --------------------------------------------------------------------

#' Profile means across indicators, labelled at the right edge
#' @noRd
.gg_view_profiles <- function(fit) {
  key <- .gg_profile_key(fit)
  means <- .gg_profile_means(fit, key)
  n_ind <- nlevels(means$indicator)
  k <- nrow(key)
  # Dodge whiskers apart so overlapping intervals stay readable.
  width <- min(0.3, 0.08 * k)
  means$x <- means$position + (means$rank - (k + 1) / 2) * width / max(k - 1, 1)

  last <- means[means$position == n_ind, , drop = FALSE]
  span <- diff(range(c(means$lower, means$upper), na.rm = TRUE))
  last$y <- .gg_spread(last$mean, gap = 0.075 * span)
  last$x <- n_ind + 0.22

  ggplot2::ggplot(means, ggplot2::aes(x, mean, colour = profile_label,
                                      group = profile_label)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           width = 0.06, linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = profile_label,
                                     fill = profile_label),
                        size = 2.6, colour = "white", stroke = 0.6) +
    ggplot2::geom_text(data = last,
                       ggplot2::aes(x = x, y = y, label = profile_label),
                       hjust = 0, size = 3.5, fontface = "bold") +
    ggplot2::scale_x_continuous(breaks = seq_len(n_ind),
                                labels = levels(means$indicator),
                                limits = c(0.7, n_ind + 1.35),
                                expand = ggplot2::expansion(0)) +
    .gg_scales(key) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Profile means",
      subtitle = sprintf(
        "%d profiles across %d indicators; whiskers are 95%% intervals",
        k, n_ind),
      x = NULL, y = "Mean") +
    .gg_theme()
}

#' Profile means as bars that start at zero
#' @noRd
.gg_view_bars <- function(fit) {
  key <- .gg_profile_key(fit)
  means <- .gg_profile_means(fit, key)
  dodge <- ggplot2::position_dodge(width = 0.82)
  ggplot2::ggplot(means, ggplot2::aes(indicator, mean, fill = profile_label)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_col(position = dodge, width = 0.78) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper,
                                        group = profile_label),
                           position = dodge, width = 0.22, linewidth = 0.4,
                           colour = "grey20", na.rm = TRUE) +
    .gg_scales(key, "fill") +
    ggplot2::labs(
      title = "Profile means",
      subtitle = "Bars start at zero; whiskers are 95% intervals",
      x = NULL, y = "Mean", fill = NULL) +
    .gg_theme() +
    ggplot2::theme(legend.position = "top",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.key.size = ggplot2::unit(0.9, "lines"))
}

#' Profile means in overall standard deviations, with a colour key
#'
#' The overall mean and SD are the mixture's own: sum(pi * mu) and
#' sqrt(sum(pi * (sigma^2 + mu^2)) - mean^2), so the view needs no raw data.
#' @noRd
.gg_view_heatmap <- function(fit) {
  key <- .gg_profile_key(fit)
  means <- .gg_profile_means(fit, key)
  means$share <- key$share[match(means$profile, key$profile)]
  moments <- lapply(split(means, means$indicator, drop = TRUE), \(rows) {
    overall <- sum(rows$share * rows$mean)
    spread <- sqrt(sum(rows$share * (rows$variance + rows$mean^2)) - overall^2)
    data.frame(indicator = rows$indicator[1L], overall = overall,
               spread = spread)
  })
  moments <- do.call(rbind, moments)
  at <- match(means$indicator, moments$indicator)
  means$fill_value <- (means$mean - moments$overall[at]) / moments$spread[at]
  limit <- max(abs(means$fill_value))
  means$text_colour <- ifelse(abs(means$fill_value) > 0.55 * limit,
                              "white", "grey15")
  means$profile_label <- factor(means$profile_label,
                                levels = rev(levels(key$label)))

  ggplot2::ggplot(means, ggplot2::aes(indicator, profile_label,
                                      fill = fill_value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", fill_value),
                                    colour = text_colour), size = 3.4) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient2(
      low = "#D33F6A", mid = "white", high = "#4A6FE3", midpoint = 0,
      limits = c(-limit, limit), name = "SDs from\noverall mean",
      guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(7, "lines"),
                                       barwidth = ggplot2::unit(0.7, "lines"))) +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      title = "How each profile departs from the overall mean",
      subtitle = "Profile mean minus the overall mean, in overall standard deviations",
      x = NULL, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.title = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                        colour = "grey30"),
                   axis.text.y = ggplot2::element_text(hjust = 1))
}

#' Each indicator's distribution per profile: density, quartiles, cases
#' @noRd
.gg_view_raincloud <- function(fit) {
  key <- .gg_profile_key(fit)
  cases <- get_results(fit, what = "assignments")
  indicators <- intersect(fit$vars,
                          unique(get_results(fit, what = "profiles")$indicator))
  k <- nrow(key)
  # Rows from the top: the largest profile sits highest.
  key$y0 <- rev(seq_len(k))[key$rank]

  pieces <- lapply(indicators, \(variable) {
    by_profile <- lapply(seq_len(k), \(row) {
      values <- cases[[variable]][cases$profile == key$profile[row]]
      values <- values[is.finite(values)]
      if (length(values) == 0L) return(NULL)
      y0 <- key$y0[row]
      label <- key$label[row]
      density <- if (length(values) >= 3L && stats::sd(values) > 0) {
        d <- stats::density(values, n = 256L, cut = 3)
        data.frame(x = d$x, height = d$y)
      }
      quartiles <- stats::quantile(values, c(0.25, 0.5, 0.75), names = FALSE)
      ordering <- rank(values, ties.method = "first")
      list(
        density = if (!is.null(density)) cbind(density, y0 = y0,
                                                profile_label = label,
                                                indicator = variable),
        box = data.frame(q1 = quartiles[1L], median = quartiles[2L],
                         q3 = quartiles[3L], y0 = y0, profile_label = label,
                         indicator = variable),
        points = data.frame(x = values,
                            y = y0 - 0.2 + 0.07 * ((ordering * 0.618) %% 1 - 0.5),
                            profile_label = label, indicator = variable)
      )
    })
    by_profile <- Filter(Negate(is.null), by_profile)
    density <- do.call(rbind, lapply(by_profile, `[[`, "density"))
    # One height scale per indicator, so ridges within a panel compare.
    if (!is.null(density)) {
      density$height <- 0.62 * density$height / max(density$height)
    }
    list(density = density,
         box = do.call(rbind, lapply(by_profile, `[[`, "box")),
         points = do.call(rbind, lapply(by_profile, `[[`, "points")))
  })
  bind <- \(part) {
    out <- do.call(rbind, lapply(pieces, `[[`, part))
    out$indicator <- factor(out$indicator, levels = indicators)
    out$profile_label <- factor(out$profile_label, levels = levels(key$label))
    out
  }
  density <- bind("density")
  box <- bind("box")
  points <- bind("points")

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = density,
                         ggplot2::aes(x = x, ymin = y0, ymax = y0 + height,
                                      fill = profile_label,
                                      group = profile_label),
                         alpha = 0.55, colour = NA) +
    ggplot2::geom_line(data = density,
                       ggplot2::aes(x = x, y = y0 + height,
                                    colour = profile_label,
                                    group = profile_label),
                       linewidth = 0.45) +
    ggplot2::geom_segment(data = box,
                          ggplot2::aes(x = q1, xend = q3, y = y0 - 0.07,
                                       yend = y0 - 0.07,
                                       colour = profile_label),
                          linewidth = 1.6) +
    ggplot2::geom_point(data = box,
                        ggplot2::aes(x = median, y = y0 - 0.07),
                        shape = 21, fill = "white", colour = "grey15",
                        size = 1.6, stroke = 0.5) +
    ggplot2::geom_point(data = points,
                        ggplot2::aes(x = x, y = y, colour = profile_label),
                        size = 0.7, alpha = 0.55) +
    ggplot2::facet_wrap(ggplot2::vars(indicator), scales = "free_x") +
    ggplot2::scale_y_continuous(breaks = key$y0, labels = key$label,
                                expand = ggplot2::expansion(add = c(0.35, 0.1))) +
    .gg_scales(key, c("colour", "fill")) +
    ggplot2::labs(
      title = "Indicator distributions by profile",
      subtitle = paste("Density, interquartile range with median, and each case,",
                       "grouped by its most likely profile"),
      x = NULL, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.3),
                   panel.spacing = ggplot2::unit(1.2, "lines"))
}

#' Effective profile sizes, largest first, labelled with count and share
#' @noRd
.gg_view_sizes <- function(fit) {
  key <- .gg_profile_key(fit)
  key$profile_label <- factor(key$label, levels = rev(levels(key$label)))
  key$bar_label <- sprintf("%.1f  (%.1f%%)", key$effective_count,
                           100 * key$share)
  total <- sum(key$effective_count)
  ggplot2::ggplot(key, ggplot2::aes(effective_count, profile_label,
                                    fill = label)) +
    ggplot2::geom_col(width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = bar_label), hjust = -0.08,
                       size = 3.4, colour = "grey20") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.22))) +
    .gg_scales(key, "fill") +
    ggplot2::labs(
      title = "Profile sizes",
      subtitle = sprintf(
        "Effective number of cases (sum of posterior probabilities); %s in total",
        format(round(total), big.mark = ",")),
      x = "Effective cases", y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey90",
                                                              linewidth = 0.35))
}

#' Per-case posterior summaries, one strip per assigned profile
#' @noRd
.gg_case_values <- function(fit) {
  cases <- get_results(fit, what = "assignments")
  posterior_columns <- grep("^posterior_profile_", names(cases), value = TRUE)
  posteriors <- as.matrix(cases[posterior_columns])
  k <- ncol(posteriors)
  entropy <- -rowSums(ifelse(posteriors > 0, posteriors * log(posteriors), 0))
  data.frame(
    profile = cases$profile,
    certainty = posteriors[cbind(seq_len(nrow(posteriors)),
                                 match(paste0("posterior_profile_", cases$profile),
                                       posterior_columns))],
    entropy = if (k > 1L) entropy / log(k) else 0
  )
}

#' Strip plot of one per-case value by profile, with each profile's mean
#' @noRd
.gg_case_strip <- function(fit, value, limits, reference, reference_label,
                           title, subtitle, x_label) {
  key <- .gg_profile_key(fit)
  cases <- .gg_case_values(fit)
  k <- nrow(key)
  key$y0 <- rev(seq_len(k))[key$rank]
  cases$value <- cases[[value]]
  cases$y0 <- key$y0[match(cases$profile, key$profile)]
  cases$profile_label <- key$label[match(cases$profile, key$profile)]
  order_within <- stats::ave(cases$value, cases$profile,
                             FUN = \(v) rank(v, ties.method = "first"))
  cases$y <- cases$y0 + 0.5 * ((order_within * 0.618) %% 1 - 0.5)
  means <- stats::aggregate(value ~ profile + y0, data = cases, FUN = mean)
  means$profile_label <- key$label[match(means$profile, key$profile)]
  means$label <- sprintf("mean %.2f", means$value)

  ggplot2::ggplot(cases, ggplot2::aes(value, y, colour = profile_label)) +
    ggplot2::geom_vline(xintercept = reference, linetype = "22",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::annotate("text", x = reference, y = k + 0.45,
                      label = reference_label, vjust = 1,
                      hjust = if (isTRUE(all.equal(reference, limits[2L])))
                        1.08 else -0.08,
                      size = 3, colour = "grey45") +
    ggplot2::geom_point(size = 1.1, alpha = 0.55) +
    ggplot2::geom_segment(data = means,
                          ggplot2::aes(x = value, xend = value,
                                       y = y0 - 0.32, yend = y0 + 0.32),
                          colour = "grey10", linewidth = 0.8) +
    ggplot2::geom_text(data = means,
                       ggplot2::aes(x = value, y = y0 + 0.34, label = label),
                       colour = "grey10", size = 3, vjust = 0,
                       hjust = if (value == "certainty") 1.05 else -0.05) +
    ggplot2::scale_x_continuous(limits = limits,
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = key$y0, labels = key$label,
                                limits = c(0.5, k + 0.5),
                                expand = ggplot2::expansion(0)) +
    .gg_scales(key, "colour") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey92",
                                                              linewidth = 0.3))
}

#' @noRd
.gg_view_posteriors <- function(fit) {
  k <- fit$n_profiles
  .gg_case_strip(
    fit, value = "certainty", limits = c(1 / k, 1),
    reference = 1 / k, reference_label = "chance",
    title = "Classification certainty",
    subtitle = paste("Posterior probability of each case's assigned profile;",
                     "the profile mean is its average posterior probability"),
    x_label = "Posterior probability of the assigned profile")
}

#' @noRd
.gg_view_entropy <- function(fit) {
  cases <- .gg_case_values(fit)
  .gg_case_strip(
    fit, value = "entropy", limits = c(0, 1),
    reference = 1, reference_label = "flat",
    title = "Classification uncertainty",
    subtitle = sprintf(paste(
      "Each case's entropy relative to a flat posterior (0 = certain, 1 = flat);",
      "relative entropy = 1 - mean = %.3f"), 1 - mean(cases$entropy)),
    x_label = "Relative case entropy")
}

#' Average posterior probability matrix, assigned by posterior profile
#' @noRd
.gg_view_avepp <- function(fit) {
  key <- .gg_profile_key(fit)
  table <- get_results(fit, what = "average_posteriors")
  table <- table[table$level == "individuals", , drop = FALSE]
  assigned_n <- table$n_assigned[match(key$profile, table$assigned_class)]
  assigned_levels <- sprintf("Assigned %d (n = %d)", key$profile, assigned_n)
  posterior_levels <- sprintf("Profile %d", key$profile)
  table$assigned_label <- factor(
    assigned_levels[match(table$assigned_class, key$profile)],
    levels = rev(assigned_levels))
  table$posterior_label <- factor(
    posterior_levels[match(table$class, key$profile)],
    levels = posterior_levels)
  table$text_colour <- ifelse(table$average_posterior > 0.55, "white", "grey15")
  diagonal <- table[table$assigned_class == table$class, , drop = FALSE]

  ggplot2::ggplot(table, ggplot2::aes(posterior_label, assigned_label,
                                      fill = average_posterior)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_tile(data = diagonal, fill = NA, colour = "grey15",
                       linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", average_posterior),
                                    colour = text_colour), size = 3.4) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient(
      low = "white", high = "#0072B2", limits = c(0, 1),
      name = "Average\nposterior",
      guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(7, "lines"),
                                       barwidth = ggplot2::unit(0.7, "lines"))) +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      title = "Average posterior probability",
      subtitle = paste("Mean posterior each assigned group places on each",
                       "profile; the outlined diagonal is classification accuracy"),
      x = NULL, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.title = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                        colour = "grey30"))
}

#' Information criteria across the candidate grid, one panel per criterion
#' @noRd
.gg_view_enumeration <- function(x) {
  grid <- as.data.frame(x)
  grid <- grid[grid$converged, , drop = FALSE]
  criteria <- c(BIC = "bic_individual", AIC = "aic", ICL = "icl_individual")
  long <- do.call(rbind, Map(\(name, column) {
    data.frame(criterion = name, model = grid$model,
               n_profiles = grid$n_profiles, value = grid[[column]])
  }, names(criteria), criteria))
  long <- long[is.finite(long$value), , drop = FALSE]
  long$criterion <- factor(long$criterion, levels = names(criteria))
  long$is_minimum <- stats::ave(long$value, long$criterion,
                                FUN = \(v) seq_along(v) == which.min(v)) == 1
  models <- unique(long$model)
  colours <- stats::setNames(rep_len(.gg_okabe_ito, length(models)), models)
  shapes <- stats::setNames(rep_len(.gg_shapes, length(models)), models)

  ends <- long[long$n_profiles == max(long$n_profiles), , drop = FALSE]
  ends <- do.call(rbind, lapply(split(ends, ends$criterion, drop = TRUE), \(p) {
    p$y <- .gg_spread(p$value, gap = 0.07 * diff(range(
      long$value[long$criterion == p$criterion[1L]])))
    p
  }))
  ends$x <- max(long$n_profiles) + 0.28
  minimum <- long[long$is_minimum, , drop = FALSE]
  profile_breaks <- sort(unique(long$n_profiles))

  ggplot2::ggplot(long, ggplot2::aes(n_profiles, value, colour = model,
                                     group = model)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = model, fill = model), size = 2.3,
                        colour = "white", stroke = 0.6) +
    ggplot2::geom_point(data = minimum, shape = 21, size = 5.2, stroke = 0.8,
                        colour = "grey10", fill = NA) +
    ggplot2::geom_text(data = if (length(models) > 1L) ends else ends[0, ],
                       ggplot2::aes(x = x, y = y, label = model),
                       hjust = 0, size = 3.2, fontface = "bold") +
    ggplot2::facet_wrap(ggplot2::vars(criterion), nrow = 1L,
                        scales = "free_y") +
    ggplot2::scale_x_continuous(
      breaks = profile_breaks,
      expand = ggplot2::expansion(add = c(0.25, if (length(models) > 1L) 0.8 else 0.25))) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Model comparison",
      subtitle = "Lower is better; the ring marks each criterion's minimum",
      x = "Number of profiles", y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.6, "lines"))
}

#' Every case as a line across indicators, one panel per profile
#'
#' Indicators are put on one scale (z-scores of each column) so that lines
#' across them mean something whatever the units. Each panel shows its
#' profile's cases in colour, faded by how certainly they belong, over all
#' other cases in grey, with the profile mean on top.
#' @noRd
.gg_view_parallel <- function(fit) {
  key <- .gg_profile_key(fit)
  cases <- get_results(fit, what = "assignments")
  indicators <- intersect(fit$vars,
                          unique(get_results(fit, what = "profiles")$indicator))
  centres <- vapply(cases[indicators], \(v) mean(v, na.rm = TRUE), numeric(1L))
  scales <- vapply(cases[indicators], \(v) stats::sd(v, na.rm = TRUE),
                   numeric(1L))
  scales[!(scales > 0)] <- 1
  certainty <- 1 - cases$uncertainty

  long <- do.call(rbind, lapply(seq_along(indicators), \(j) {
    data.frame(case = seq_len(nrow(cases)), profile = cases$profile,
               certainty = certainty, position = j,
               value = (cases[[indicators[j]]] - centres[j]) / scales[j])
  }))
  long <- long[is.finite(long$value), , drop = FALSE]
  long$profile_label <- key$label[match(long$profile, key$profile)]

  # Grey context: every case, repeated in every panel.
  context <- do.call(rbind, lapply(levels(key$label), \(panel) {
    cbind(long[c("case", "position", "value")], panel = panel)
  }))
  context$panel <- factor(context$panel, levels = levels(key$label))
  long$panel <- long$profile_label

  means <- .gg_profile_means(fit, key)
  at <- match(as.character(means$indicator), indicators)
  means$value <- (means$mean - centres[at]) / scales[at]
  means$panel <- means$profile_label

  ggplot2::ggplot() +
    ggplot2::geom_line(data = context,
                       ggplot2::aes(position, value, group = case),
                       colour = "grey88", linewidth = 0.25) +
    ggplot2::geom_line(data = long,
                       ggplot2::aes(position, value, group = case,
                                    colour = profile_label, alpha = certainty),
                       linewidth = 0.35) +
    ggplot2::geom_line(data = means,
                       ggplot2::aes(position, value, group = panel),
                       colour = "grey10", linewidth = 1.1) +
    ggplot2::geom_point(data = means,
                        ggplot2::aes(position, value, fill = profile_label,
                                     shape = profile_label),
                        colour = "grey10", size = 2.4, stroke = 0.6) +
    ggplot2::facet_wrap(ggplot2::vars(panel)) +
    ggplot2::scale_x_continuous(breaks = seq_along(indicators),
                                labels = indicators,
                                expand = ggplot2::expansion(add = 0.15)) +
    ggplot2::scale_alpha_continuous(range = c(0.12, 0.7), limits = c(0, 1)) +
    .gg_scales(key) +
    ggplot2::labs(
      title = "Cases across indicators",
      subtitle = paste("Each line is a case, faded by the certainty of its",
                       "assignment; the dark line is the profile mean;",
                       "other cases in grey"),
      x = NULL, y = "Standardized value (z)") +
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

#' Scatter-plot matrix: cases by profile with each profile's 95% ellipse
#'
#' Lower triangle only. Larger points are less certainly assigned. The
#' diagonal shows each profile's density, drawn into the variable's range.
#' @noRd
.gg_view_pairs <- function(fit) {
  key <- .gg_profile_key(fit)
  cases <- get_results(fit, what = "assignments")
  means <- get_results(fit, what = "profiles")
  covariances <- get_results(fit, what = "covariances")
  indicators <- intersect(fit$vars, unique(means$indicator))
  d <- length(indicators)
  if (d < 2L) {
    stop(errorCondition(
      "`what = \"pairs\"` needs at least two continuous indicators.",
      class = "latents_bad_argument", call = NULL))
  }
  cases$profile_label <- key$label[match(cases$profile, key$profile)]
  pairs <- expand.grid(col = seq_len(d), row = seq_len(d))
  pairs <- pairs[pairs$col < pairs$row, , drop = FALSE]

  scatter <- do.call(rbind, Map(\(i, j) {
    data.frame(x = cases[[indicators[i]]], y = cases[[indicators[j]]],
               uncertainty = cases$uncertainty,
               profile_label = cases$profile_label,
               xvar = indicators[i], yvar = indicators[j])
  }, pairs$col, pairs$row))
  scatter <- scatter[is.finite(scatter$x) & is.finite(scatter$y), ,
                     drop = FALSE]

  covariance_of <- \(profile, a, b) {
    rows <- covariances$profile == profile &
      covariances$indicator == a & covariances$indicator_2 == b
    covariances$covariance[rows][1L]
  }
  mean_of <- \(profile, a) {
    means$mean[means$profile == profile & means$indicator == a]
  }
  ellipses <- do.call(rbind, Map(\(i, j) {
    a <- indicators[i]
    b <- indicators[j]
    do.call(rbind, lapply(seq_len(nrow(key)), \(row) {
      profile <- key$profile[row]
      sigma <- matrix(c(covariance_of(profile, a, a),
                        covariance_of(profile, a, b),
                        covariance_of(profile, b, a),
                        covariance_of(profile, b, b)), 2L)
      ellipse <- .gg_ellipse(c(mean_of(profile, a), mean_of(profile, b)),
                             sigma)
      cbind(ellipse, profile_label = key$label[row], xvar = a, yvar = b)
    }))
  }, pairs$col, pairs$row))

  diagonal <- do.call(rbind, lapply(indicators, \(variable) {
    values <- cases[[variable]]
    range_v <- range(values, na.rm = TRUE)
    curves <- lapply(seq_len(nrow(key)), \(row) {
      mine <- values[cases$profile == key$profile[row] & is.finite(values)]
      if (length(mine) < 3L || !(stats::sd(mine) > 0)) return(NULL)
      dens <- stats::density(mine, n = 200L, from = range_v[1L],
                             to = range_v[2L])
      data.frame(x = dens$x, height = dens$y * key$share[row],
                 profile_label = key$label[row])
    })
    curves <- do.call(rbind, curves)
    # Mixture-weighted densities, drawn into the variable's own range.
    curves$y <- range_v[1L] + 0.92 * diff(range_v) *
      curves$height / max(curves$height)
    curves$ymin <- range_v[1L]
    cbind(curves, xvar = variable, yvar = variable)
  }))

  panels <- rbind(unique(scatter[c("xvar", "yvar")]),
                  unique(diagonal[c("xvar", "yvar")]))
  as_levels <- \(frame) {
    frame$xvar <- factor(frame$xvar, levels = indicators)
    frame$yvar <- factor(frame$yvar, levels = indicators)
    if ("profile_label" %in% names(frame)) {
      frame$profile_label <- factor(frame$profile_label,
                                    levels = levels(key$label))
    }
    frame
  }
  scatter <- as_levels(scatter)
  ellipses <- as_levels(ellipses)
  diagonal <- as_levels(diagonal)
  panels <- as_levels(panels)

  ggplot2::ggplot() +
    ggplot2::geom_rect(data = panels, xmin = -Inf, xmax = Inf, ymin = -Inf,
                       ymax = Inf, fill = "grey97", colour = NA) +
    ggplot2::geom_ribbon(data = diagonal,
                         ggplot2::aes(x = x, ymin = ymin, ymax = y,
                                      fill = profile_label,
                                      group = profile_label),
                         alpha = 0.35, colour = NA) +
    ggplot2::geom_line(data = diagonal,
                       ggplot2::aes(x = x, y = y, colour = profile_label,
                                    group = profile_label),
                       linewidth = 0.5) +
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
    .gg_scales(key, c("colour", "fill")) +
    ggplot2::labs(
      title = "Profiles in each pair of indicators",
      subtitle = paste0("Cases coloured by profile, larger when less certain; ",
                        "95% ellipses of the fitted covariances\n",
                        "Diagonal: each profile's density, weighted by its size"),
      x = NULL, y = NULL) +
    .gg_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text = ggplot2::element_text(face = "bold",
                                                      size = ggplot2::rel(0.85)),
                   strip.text.y.left = ggplot2::element_text(angle = 90),
                   panel.spacing = ggplot2::unit(0.5, "lines"),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_blank()) +
    ggplot2::guides(fill = "none")
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
#' profile and one in the next solution, sum over cases of p(a) * p(b), so
#' each profile's bands add up exactly to its share. A band fanning out is a
#' split; bands crossing are cases reshuffled. The layout follows clustering
#' trees (Zappia and Oshlack, 2018). One panel per covariance structure.
#' @noRd
.gg_view_tree <- function(x) {
  grid <- as.data.frame(x)
  grid <- grid[grid$converged, , drop = FALSE]
  half <- 0.16
  per_model <- lapply(unique(grid$model), \(model) {
    counts <- sort(unique(grid$n_profiles[grid$model == model]))
    if (length(counts) < 2L) return(NULL)
    posteriors <- lapply(counts, \(k) {
      candidate_fit(x, n_profiles = k, model = model)$subject_posteriors
    })
    n_cases <- nrow(posteriors[[1L]])

    # Order each row's profiles by the flow-weighted centre of their
    # parents, so bands cross only where cases really are reshuffled.
    first_share <- colSums(posteriors[[1L]]) / n_cases
    first_order <- order(-first_share, seq_along(first_share))
    # Segments are separated by a fixed gap and each row is centred, so a
    # split shows as one band forking into two.
    gap <- 0.03
    segments_of <- \(share, ordering) {
      widths <- share[ordering]
      k <- length(widths)
      right <- cumsum(widths) + gap * (seq_len(k) - 1) - gap * (k - 1) / 2
      data.frame(profile = ordering, xmin = right - widths, xmax = right,
                 share = widths)
    }
    rows <- Reduce(\(previous, level) {
      parent <- previous[[length(previous)]]
      flow <- crossprod(posteriors[[level - 1L]], posteriors[[level]]) / n_cases
      centres <- (parent$xmin + parent$xmax)[order(parent$profile)] / 2
      child_centre <- colSums(flow * centres) / colSums(flow)
      share <- colSums(flow)
      ordering <- order(child_centre, seq_along(child_centre))
      c(previous, list(segments_of(share, ordering)))
    }, seq_along(counts)[-1L],
    init = list(segments_of(first_share, first_order)))

    # Colour is lineage: a profile that is its dominant parent's main
    # continuation keeps the parent's colour; any other profile is new and
    # takes the next colour. Colours are indexed by profile id.
    first_colours <- character(counts[1L])
    first_colours[first_order] <- rep_len(.gg_okabe_ito, counts[1L])
    lineage <- Reduce(\(state, level) {
      flow <- crossprod(posteriors[[level - 1L]], posteriors[[level]])
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

    bars <- do.call(rbind, Map(\(k, segments, colours) {
      cbind(segments, model = model, n_profiles = k,
            fill = colours[segments$profile])
    }, counts, rows, lineage))

    bands <- do.call(rbind, lapply(seq_along(counts)[-1L], \(level) {
      parent <- rows[[level - 1L]]
      child <- rows[[level]]
      flow <- crossprod(posteriors[[level - 1L]], posteriors[[level]]) / n_cases
      pairs <- expand.grid(from = parent$profile, to = child$profile)
      pairs$mass <- flow[cbind(pairs$from, pairs$to)]
      pairs <- pairs[pairs$mass > 1e-3, , drop = FALSE]
      # Stack outflows inside a parent in child order, inflows inside a
      # child in parent order, so bands tile both bars without overlap.
      pairs$from_rank <- match(pairs$from, parent$profile)
      pairs$to_rank <- match(pairs$to, child$profile)
      out <- pairs[order(pairs$from_rank, pairs$to_rank), , drop = FALSE]
      out$top_right <- parent$xmin[out$from_rank] +
        stats::ave(out$mass, out$from, FUN = cumsum)
      out$top_left <- out$top_right - out$mass
      into <- out[order(out$to_rank, out$from_rank), , drop = FALSE]
      into$bottom_right <- child$xmin[into$to_rank] +
        stats::ave(into$mass, into$to, FUN = cumsum)
      into$bottom_left <- into$bottom_right - into$mass
      do.call(rbind, lapply(seq_len(nrow(into)), \(i) {
        band <- .gg_flow_band(into$top_left[i], into$top_right[i],
                              into$bottom_left[i], into$bottom_right[i],
                              y_top = counts[level - 1L] + half,
                              y_bottom = counts[level] - half)
        cbind(band, model = model,
              group = sprintf("%d-%d-%d", counts[level], into$from[i],
                              into$to[i]),
              fill = lineage[[level - 1L]][into$from[i]])
      }))
    }))
    list(bars = bars, bands = bands)
  })
  per_model <- Filter(Negate(is.null), per_model)
  if (length(per_model) == 0L) {
    stop(errorCondition(
      "This view needs at least two numbers of profiles for one structure.",
      class = "latents_bad_argument", call = NULL))
  }
  bars <- do.call(rbind, lapply(per_model, `[[`, "bars"))
  bands <- do.call(rbind, lapply(per_model, `[[`, "bands"))
  percent <- 100 * bars$share
  bars$label <- ifelse(percent < 1, "<1%", sprintf("%.0f%%", percent))
  # Only write a share where the segment is wide enough to hold it.
  bars$label[bars$share < 0.05] <- ""
  # Dark text on light fills, white on dark ones.
  rgb <- grDevices::col2rgb(bars$fill)
  luminance <- 0.299 * rgb[1L, ] + 0.587 * rgb[2L, ] + 0.114 * rgb[3L, ]
  bars$text_colour <- ifelse(luminance > 150, "grey10", "white")
  counts <- sort(unique(bars$n_profiles))

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
    ggplot2::facet_wrap(ggplot2::vars(model), nrow = 1L) +
    ggplot2::scale_y_reverse(breaks = counts,
                             expand = ggplot2::expansion(add = 0.25)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    ggplot2::labs(
      title = "How profiles split as more are added",
      subtitle = paste0("Each row is one solution, divided into its profiles ",
                        "by share of cases\n",
                        "Colour follows a profile across solutions; a new ",
                        "colour marks a profile that split off"),
      caption = paste0("A band that forks: a profile split in two ",
                       "(nested, stable solutions). Bands that merge or ",
                       "cross: cases reshuffled between profiles."),
      x = NULL, y = "Number of profiles") +
    .gg_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(1.5, "lines"))
}

#' Build one view as a ggplot object
#'
#' @param x A `multilpa` fit, or a `multilpa_enumeration` for
#'   `what = "enumeration"` and `"tree"`.
#' @param what The view.
#' @return A ggplot object.
#' @noRd
.gg_plot <- function(x, what = c("profiles", "bars", "heatmap", "raincloud",
                                 "parallel", "pairs", "sizes", "entropy",
                                 "posteriors", "avepp", "enumeration",
                                 "tree")) {
  .gg_require()
  what <- match.arg(what)
  if (what %in% c("enumeration", "tree")) {
    stopifnot("`x` must be an enumeration for this view" =
                inherits(x, "multilpa_enumeration"))
    return(if (identical(what, "tree")) .gg_view_tree(x) else
      .gg_view_enumeration(x))
  }
  stopifnot("`x` must be a multilpa fit" = inherits(x, "multilpa"))
  switch(what,
    profiles = .gg_view_profiles(x),
    bars = .gg_view_bars(x),
    heatmap = .gg_view_heatmap(x),
    raincloud = .gg_view_raincloud(x),
    parallel = .gg_view_parallel(x),
    pairs = .gg_view_pairs(x),
    sizes = .gg_view_sizes(x),
    entropy = .gg_view_entropy(x),
    posteriors = .gg_view_posteriors(x),
    avepp = .gg_view_avepp(x))
}
