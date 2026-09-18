#' Plot a fitted multilevel latent profile model
#'
#' Draws the profile means across indicators, or the profile prevalence within
#' each latent group class. Profiles are distinguished by colour, point symbol
#' and line type together, and labelled directly, so the plot stays readable in
#' greyscale and needs no legend.
#'
#' @param x A fitted `multilpa` model.
#' @param what `"profiles"` plots the Gaussian measurement model, one line per
#'   profile across the continuous indicators. `"responses"` plots the
#'   categorical measurement model, one line per profile across the categorical
#'   indicators, showing the probability of a chosen category. `"probabilities"`
#'   plots profile prevalence within each group class, one line per group class.
#' @param scale For `what = "profiles"`, `"raw"` plots the estimated means in
#'   input units, and `"standardized"` divides each indicator's deviation from
#'   its grand mean by that indicator's observed standard deviation. Use
#'   `"standardized"` when indicators are on different scales, where raw means
#'   make the largest-scale indicator dominate the shape.
#' @param category For `what = "responses"`, which category's probability to
#'   plot. `"last"` uses each indicator's highest category, which is the usual
#'   choice for binary indicators, `"first"` uses the lowest, or give a single
#'   category label or index used for every indicator.
#' @param labels `TRUE` prints a direct label at the right end of each series.
#' @param main,subtitle Panel title and secondary line. `NULL` for none,
#'   `waiver` is not used; pass `""` to reserve the space without text.
#' @param palette,symbols,linetypes Vectors of colours, plotting characters and
#'   line types, recycled to the number of series. Defaults are the Okabe-Ito
#'   palette and matched symbol and line-type sequences.
#' @param style A list of visual constants, as built by `.multilpa_style()`;
#'   pass named elements to override individual constants.
#' @param ... Further named visual constants, merged into `style`.
#' @return The fitted model, invisibly. Called for the side effect of drawing.
#' @details The plot shows point estimates only. It carries no standard errors,
#'   and profile order is arbitrary, so two fits must have their labels aligned
#'   before their plots are compared. With `scale = "standardized"` the
#'   standard deviations are the observed indicator standard deviations, not
#'   the within-profile residual standard deviations, so the plotted values are
#'   comparable across indicators but are not effect sizes.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' plot(fit)
#' plot(fit, scale = "standardized")
#' @export
plot.multilpa <- function(x, what = c("profiles", "responses", "probabilities"),
                        scale = c("raw", "standardized"), category = "last",
                        labels = TRUE, main = NULL, subtitle = NULL,
                        palette = NULL, symbols = NULL, linetypes = NULL,
                        style = .multilpa_style(), ...) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`category` must be a single label or index" = length(category) == 1L)
  what <- match.arg(what)
  scale <- match.arg(scale)
  style <- utils::modifyList(style, list(...))
  previous <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA)
  switch(what,
    profiles = .multilpa_plot_profiles(x, scale, labels, main, subtitle, palette,
                                     symbols, linetypes, style),
    responses = .multilpa_plot_responses(x, category, labels, main, subtitle,
                                       palette, symbols, linetypes, style),
    probabilities = .multilpa_plot_probabilities(x, labels, main, subtitle,
                                               palette, symbols, linetypes, style))
  invisible(x)
}

#' Draw categorical response probabilities across indicators
#' @param x A fitted `multilpa` model with categorical indicators.
#' @param category Which category to plot, as `"last"`, `"first"`, a label, or
#'   an index.
#' @param labels Whether to draw direct series labels.
#' @param main,subtitle Panel title and secondary line.
#' @param palette,symbols,linetypes Series aesthetics, or `NULL` for defaults.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_responses <- function(x, category, labels, main, subtitle, palette,
                                   symbols, linetypes, style) {
  blocks <- x$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) {
    stop(errorCondition("This model has no categorical indicators.",
                        class = "multilpa_no_categorical", call = NULL))
  }
  indicators <- names(blocks)
  n_profiles <- x$n_profiles
  chosen <- vapply(blocks, function(block) {
    index <- if (identical(category, "last")) ncol(block) else
      if (identical(category, "first")) 1L else
        if (is.numeric(category)) as.integer(category) else
          match(as.character(category), colnames(block))
    if (is.na(index) || index < 1L || index > ncol(block)) {
      stop(errorCondition(sprintf("`category` does not match a category of an indicator with categories %s.",
                                  paste(colnames(block), collapse = ", ")),
                          class = "multilpa_unknown_category", call = NULL))
    }
    index
  }, integer(1))
  values <- t(vapply(seq_along(blocks), function(index) {
    blocks[[index]][, chosen[[index]]]
  }, numeric(n_profiles)))
  category_labels <- vapply(seq_along(blocks), function(index) {
    colnames(blocks[[index]])[chosen[[index]]]
  }, character(1))
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  points <- if (is.null(symbols)) .multilpa_symbols(n_profiles) else
    rep(symbols, length.out = n_profiles)
  lines <- if (is.null(linetypes)) .multilpa_linetypes(n_profiles) else
    rep(linetypes, length.out = n_profiles)
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  label_text <- sprintf("Profile %d (%.0f%%)", seq_len(n_profiles), 100 * share)
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  positions <- seq_along(indicators)
  # The axis names each indicator with the category being plotted, so a mixed
  # set of category labels cannot be misread as one shared category.
  shared <- length(unique(category_labels)) == 1L
  axis_labels <- if (shared) indicators else
    sprintf("%s=%s", indicators, category_labels)
  .multilpa_panel(xlim = c(1 - 0.35, length(indicators) + 0.35), ylim = c(0, 1.02),
    xlab = "Categorical indicator",
    ylab = if (shared) sprintf("P(response = %s)", category_labels[[1L]]) else
      "Response probability",
    main = if (is.null(main)) "Response probabilities by profile" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d profiles, %d group class%s, %d categorical indicator%s", n_profiles,
      x$n_group_classes, if (x$n_group_classes == 1L) "" else "es",
      length(indicators), if (length(indicators) == 1L) "" else "s") else subtitle,
    x_at = positions, x_labels = axis_labels, y_at = seq(0, 1, by = 0.25),
    style = style)
  invisible(lapply(seq_len(n_profiles), function(profile) {
    graphics::lines(positions, values[, profile], col = colours[profile],
                    lwd = style$line_width, lty = lines[profile])
    graphics::points(positions, values[, profile], pch = points[profile],
                     bg = colours[profile], col = style$panel_fill,
                     cex = style$point_size, lwd = 1.4)
  }))
  if (isTRUE(labels)) {
    label_y <- .multilpa_spread_labels(values[length(indicators), ], 0.055)
    graphics::text(length(indicators) + 0.45, label_y, label_text,
                   adj = c(0, 0.5), col = colours,
                   cex = style$label_text_size, font = 2L)
  }
  invisible(NULL)
}

#' Draw profile means across indicators
#' @param x A fitted `multilpa` model.
#' @param scale Raw or standardized means.
#' @param labels Whether to draw direct series labels.
#' @param main,subtitle Panel title and secondary line.
#' @param palette,symbols,linetypes Series aesthetics, or `NULL` for defaults.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_profiles <- function(x, scale, labels, main, subtitle, palette,
                                  symbols, linetypes, style) {
  indicators <- .multilpa_continuous_names(x)
  n_profiles <- x$n_profiles
  means <- x$means
  if (length(indicators) == 0L) {
    stop(errorCondition("This model has no continuous indicators; plot `what = \"responses\"` instead.",
                        class = "multilpa_no_continuous", call = NULL))
  }
  if (identical(scale, "standardized")) {
    observed <- x$indicator_data
    if (is.null(observed)) {
      stop(errorCondition("This fit did not retain indicator data, so it cannot be standardized.",
                          class = "multilpa_no_indicator_data", call = NULL))
    }
    centre <- colMeans(observed, na.rm = TRUE)
    spread <- vapply(seq_along(indicators), function(index) {
      stats::sd(observed[, index], na.rm = TRUE)
    }, numeric(1))
    if (any(!is.finite(spread)) || any(spread <= 0)) {
      stop(errorCondition("An indicator has zero or undefined standard deviation.",
                          class = "multilpa_bad_scale", call = NULL))
    }
    means <- sweep(sweep(means, 2L, centre, "-"), 2L, spread, "/")
  }
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  points <- if (is.null(symbols)) .multilpa_symbols(n_profiles) else
    rep(symbols, length.out = n_profiles)
  lines <- if (is.null(linetypes)) .multilpa_linetypes(n_profiles) else
    rep(linetypes, length.out = n_profiles)
  positions <- seq_along(indicators)
  span <- range(means)
  padding <- 0.12 * max(diff(span), .Machine$double.eps)
  ylim <- c(span[1L] - padding, span[2L] + padding)
  xlim <- c(1 - 0.35, length(indicators) + 0.35)
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  label_text <- sprintf("Profile %d (%.0f%%)", seq_len(n_profiles), 100 * share)
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  default_main <- sprintf("Profile means across %d indicator%s",
                          length(indicators),
                          if (length(indicators) == 1L) "" else "s")
  default_subtitle <- sprintf("%d profiles, %d group class%s; %s scale",
    n_profiles, x$n_group_classes,
    if (x$n_group_classes == 1L) "" else "es",
    if (identical(scale, "standardized")) "standardized" else "input")
  .multilpa_panel(xlim = xlim, ylim = ylim, xlab = "Indicator",
                ylab = if (identical(scale, "standardized"))
                  "Standardized mean" else "Estimated mean",
                main = if (is.null(main)) default_main else main,
                subtitle = if (is.null(subtitle)) default_subtitle else subtitle,
                x_at = positions, x_labels = indicators, style = style)
  if (identical(scale, "standardized")) {
    graphics::abline(h = 0, col = style$muted_colour, lwd = 1, lty = 3L)
  }
  invisible(lapply(seq_len(n_profiles), function(profile) {
    values <- means[profile, ]
    graphics::lines(positions, values, col = colours[profile],
                    lwd = style$line_width, lty = lines[profile])
    graphics::points(positions, values, pch = points[profile],
                     bg = colours[profile], col = style$panel_fill,
                     cex = style$point_size, lwd = 1.4)
  }))
  if (isTRUE(labels)) {
    label_y <- .multilpa_spread_labels(means[, length(indicators)],
                                     0.055 * diff(ylim))
    graphics::text(length(indicators) + 0.45, label_y, label_text, adj = c(0, 0.5),
                   col = colours, cex = style$label_text_size, font = 2L)
  }
  invisible(NULL)
}

#' Draw profile prevalence within each group class
#' @param x A fitted `multilpa` model.
#' @param labels Whether to draw direct series labels.
#' @param main,subtitle Panel title and secondary line.
#' @param palette,symbols,linetypes Series aesthetics, or `NULL` for defaults.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_probabilities <- function(x, labels, main, subtitle, palette,
                                       symbols, linetypes, style) {
  n_types <- x$n_group_classes
  n_profiles <- x$n_profiles
  probabilities <- x$profile_probabilities
  colours <- if (is.null(palette)) .multilpa_palette(n_types) else
    rep(palette, length.out = n_types)
  points <- if (is.null(symbols)) .multilpa_symbols(n_types) else
    rep(symbols, length.out = n_types)
  lines <- if (is.null(linetypes)) .multilpa_linetypes(n_types) else
    rep(linetypes, length.out = n_types)
  positions <- seq_len(n_profiles)
  xlim <- c(1 - 0.35, n_profiles + 0.35)
  label_text <- sprintf("Class %d (%.0f%%)", seq_len(n_types),
                        100 * x$group_probabilities)
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  .multilpa_panel(xlim = xlim, ylim = c(0, 1.02), xlab = "Profile",
                ylab = "Probability within group class",
                main = if (is.null(main)) "Profile prevalence by group class" else main,
                subtitle = if (is.null(subtitle))
                  sprintf("%d group class%s over %d observed groups", n_types,
                          if (n_types == 1L) "" else "es", x$n_groups) else subtitle,
                x_at = positions, x_labels = sprintf("Profile %d", positions),
                y_at = seq(0, 1, by = 0.25), style = style)
  invisible(lapply(seq_len(n_types), function(group_type) {
    values <- probabilities[group_type, ]
    graphics::lines(positions, values, col = colours[group_type],
                    lwd = style$line_width, lty = lines[group_type])
    graphics::points(positions, values, pch = points[group_type],
                     bg = colours[group_type], col = style$panel_fill,
                     cex = style$point_size, lwd = 1.4)
  }))
  if (isTRUE(labels)) {
    label_y <- .multilpa_spread_labels(probabilities[, n_profiles], 0.055)
    graphics::text(n_profiles + 0.45, label_y, label_text, adj = c(0, 0.5),
                   col = colours, cex = style$label_text_size, font = 2L)
  }
  invisible(NULL)
}

#' Plot a class-enumeration grid
#'
#' Draws one information criterion against the number of profiles, with one
#' line per number of group classes. Candidates that failed to converge are
#' marked rather than dropped, so a gap in a line is visible as a failure and
#' not mistaken for a missing candidate.
#'
#' @param x An `multilpa_enumeration` result from [enumerate_classes()].
#' @param criterion Name of the column to plot, as it appears in
#'   `as.data.frame(x)`, for example `"bic_individual"` or `"sabic_groups"`.
#' @param labels `TRUE` prints a direct label at the right end of each series.
#' @param mark_minimum `TRUE` rings the lowest value of the criterion. This
#'   marks an extremum, it does not select a model.
#' @param main,subtitle Panel title and secondary line.
#' @param palette,symbols,linetypes Series aesthetics, recycled over the number
#'   of group-class counts.
#' @param style A list of visual constants, as built by `.multilpa_style()`.
#' @param ... Further named visual constants, merged into `style`.
#' @return The enumeration result, invisibly. Called for its drawing side effect.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' candidates <- enumerate_classes(example_data, c("score_a", "score_b"), "school",
#'                                profiles = 1:3, group_classes = 1, n_starts = 2,
#'                                seed = 1)
#' plot(candidates)
#' @export
plot.multilpa_enumeration <- function(x, criterion = "bic_individual",
                                    labels = TRUE, mark_minimum = TRUE,
                                    main = NULL, subtitle = NULL, palette = NULL,
                                    symbols = NULL, linetypes = NULL,
                                    style = .multilpa_style(), ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"),
            "`criterion` must be a single column name" =
              is.character(criterion) && length(criterion) == 1L,
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`mark_minimum` must be TRUE or FALSE" =
              isTRUE(mark_minimum) || isFALSE(mark_minimum))
  grid <- as.data.frame(x)
  if (!criterion %in% names(grid)) {
    stop(errorCondition(sprintf("`%s` is not a column of the enumeration grid.",
                                criterion),
                        class = "multilpa_unknown_criterion", call = NULL))
  }
  values <- grid[[criterion]]
  if (all(is.na(values))) {
    stop(errorCondition(sprintf("Every candidate has a missing `%s`; nothing to plot.",
                                criterion),
                        class = "multilpa_nothing_to_plot", call = NULL))
  }
  style <- utils::modifyList(style, list(...))
  previous <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA)
  class_counts <- sort(unique(grid$n_group_classes))
  profile_counts <- sort(unique(grid$n_profiles))
  colours <- if (is.null(palette)) .multilpa_palette(length(class_counts)) else
    rep(palette, length.out = length(class_counts))
  points <- if (is.null(symbols)) .multilpa_symbols(length(class_counts)) else
    rep(symbols, length.out = length(class_counts))
  lines <- if (is.null(linetypes)) .multilpa_linetypes(length(class_counts)) else
    rep(linetypes, length.out = length(class_counts))
  span <- range(values, na.rm = TRUE)
  padding <- 0.12 * max(diff(span), .Machine$double.eps)
  xlim <- c(min(profile_counts) - 0.35, max(profile_counts) + 0.35)
  failures <- sum(!grid$converged)
  label_text <- sprintf("%d group class%s", class_counts,
                        ifelse(class_counts == 1L, "", "es"))
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  .multilpa_panel(xlim = xlim, ylim = c(span[1L] - padding, span[2L] + padding),
    xlab = "Number of profiles", ylab = criterion,
    main = if (is.null(main)) sprintf("%s by candidate model", criterion) else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d candidates%s; lower is preferred, no candidate is selected automatically",
      nrow(grid),
      if (failures == 0L) "" else sprintf(", %d did not converge", failures)) else
      subtitle,
    x_at = profile_counts, x_labels = profile_counts, style = style)
  ends <- vapply(seq_along(class_counts), function(index) {
    rows <- grid$n_group_classes == class_counts[index]
    series <- data.frame(profiles = grid$n_profiles[rows], value = values[rows])
    series <- series[order(series$profiles), , drop = FALSE]
    finite <- !is.na(series$value)
    graphics::lines(series$profiles[finite], series$value[finite],
                    col = colours[index], lwd = style$line_width,
                    lty = lines[index])
    graphics::points(series$profiles[finite], series$value[finite],
                     pch = points[index], bg = colours[index],
                     col = style$panel_fill, cex = style$point_size, lwd = 1.4)
    # A failed candidate is drawn as a hollow cross on the axis baseline so the
    # gap in the line is explained rather than silently missing.
    if (any(!finite)) {
      # Offset by series so that two candidates failing at the same profile
      # count remain visibly distinct rather than drawing on top of each other.
      offset <- (index - (length(class_counts) + 1) / 2) * 0.09
      graphics::points(series$profiles[!finite] + offset,
                       rep(span[1L] - padding * 0.55, sum(!finite)),
                       pch = 4L, col = colours[index], cex = style$point_size,
                       lwd = 2)
    }
    if (any(finite)) series$value[finite][sum(finite)] else NA_real_
  }, numeric(1))
  if (isTRUE(mark_minimum)) {
    best <- which.min(values)
    graphics::points(grid$n_profiles[best], values[best], pch = 1L,
                     col = style$title_colour, cex = style$point_size * 2.4,
                     lwd = 1.8)
  }
  if (isTRUE(labels)) {
    keep <- !is.na(ends)
    if (any(keep)) {
      graphics::text(max(profile_counts) + 0.45,
                     .multilpa_spread_labels(ends[keep], 0.055 * diff(span)),
                     label_text[keep],
                     adj = c(0, 0.5), col = colours[keep],
                     cex = style$label_text_size, font = 2L)
    }
  }
  invisible(x)
}
