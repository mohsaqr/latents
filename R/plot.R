#' Plot a fitted multilevel latent profile model
#'
#' Draws the profile means across indicators, or the profile prevalence within
#' each latent group class. Series are distinguished by colour, point symbol and
#' line type together, and labelled directly, so a line plot stays readable in
#' greyscale and needs no legend; the sequence grid carries the profile number
#' in each cell for the same reason.
#'
#' @param x A fitted `multilpa` model.
#' @param data Optional. The data frame the model was fitted to, used by
#'   `what = "bars"` to put a 95% interval on every bar and ignored by every
#'   other view. A fit carries the columns it was built from, so the intervals
#'   are drawn without this being supplied; pass it only to draw them from a
#'   different frame. A model family whose standard errors are not implemented
#'   gets bars without whiskers rather than an error.
#' @param what Which view to draw. [multilpa_plot_types()] lists every value
#'   with its group and a one-line description; the `"enumeration"` row it also
#'   lists belongs to [plot.multilpa_enumeration()], not to this method.
#'
#'   The measurement model, three ways: `"profiles"` draws one line per profile
#'   across the continuous indicators, with each profile's marker area
#'   proportional to its prevalence. `"bars"` draws the same means as grouped
#'   bars, one bar per profile within each indicator, with a 95% interval on
#'   every bar when `data` is supplied. `"heatmap"` draws them as a diverging
#'   grid of standard deviations from each indicator's grand mean, which is the
#'   quickest read when there are many indicators or many profiles.
#'   `"responses"` is the categorical counterpart of `"profiles"`, one line per
#'   profile across the categorical indicators, showing the probability of a
#'   chosen category.
#'
#'   The two-level structure: `"probabilities"` plots profile prevalence within
#'   each group class, one line per group class -- the quantity the second level
#'   exists to estimate. `"sequences"` draws one row per group, one column per
#'   position, coloured and numbered by the assigned profile, with the groups
#'   grouped by their latent class; it needs a fit made with `time =`.
#'
#'   Classification quality: `"entropy"` draws each case's entropy contribution
#'   as one ridge per profile, and `"posteriors"` draws the posterior
#'   probability of each case's assigned profile the same way. Both ridges are
#'   scaled to their own maximum, following the usual ridgeline convention, so
#'   ridge height compares shapes and not profile sizes; prevalence is printed
#'   in each profile's label instead. Both read the individual posteriors alone,
#'   so every family of this package can draw them; a one-profile fit refuses
#'   them with an error of class `multilpa_nothing_to_plot`, because every case
#'   then belongs to the single profile with probability one.
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
#' @param cell_labels For `what = "sequences"`, `TRUE` prints the profile number
#'   inside each cell, so the profile is never carried by colour alone. The
#'   numbers are drawn only where the cell is wide and tall enough to hold one.
#' @param main,subtitle Panel title and secondary line. `NULL` for none; pass
#'   `""` to reserve the space without text.
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
#' plot(fit, what = "bars")
#' plot(fit, what = "heatmap")
#' plot(fit, what = "entropy")
#' multilpa_plot_types()
#' @seealso [multilpa_plot_types()] for the catalogue of views.
#' @export
plot.multilpa <- function(x, what = c("profiles", "bars", "heatmap", "responses",
                                      "probabilities", "sequences",
                                      "entropy", "posteriors"),
                        data = NULL,
                        scale = c("raw", "standardized"), category = "last",
                        labels = TRUE, main = NULL, subtitle = NULL,
                        palette = NULL, symbols = NULL, linetypes = NULL,
                        style = .multilpa_style(), cell_labels = TRUE, ...) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`cell_labels` must be TRUE or FALSE" =
              isTRUE(cell_labels) || isFALSE(cell_labels),
            "`category` must be a single label or index" = length(category) == 1L)
  what <- match.arg(what)
  scale <- match.arg(scale)
  style <- utils::modifyList(style, list(...))
  previous <- graphics::par(no.readonly = TRUE)
  # `mfg` is dropped before the state is restored: setting it switches `new` on
  # as a documented side effect, so restoring it on a device nothing has been
  # drawn to yet both warns and leaves `new = TRUE` behind. That is exactly the
  # state after a refused view, where the restore runs before anything is drawn.
  previous$mfg <- NULL
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA)
  switch(what,
    profiles = .multilpa_plot_profiles(x, scale, labels, main, subtitle, palette,
                                     symbols, linetypes, style),
    responses = .multilpa_plot_responses(x, category, labels, main, subtitle,
                                       palette, symbols, linetypes, style),
    probabilities = .multilpa_plot_probabilities(x, labels, main, subtitle,
                                               palette, symbols, linetypes, style),
    sequences = .multilpa_plot_sequences(x, labels, main, subtitle, palette,
                                        style, cell_labels),
    bars = .multilpa_plot_bars(x, scale, .multilpa_mean_error_matrix(x, data),
                               main, subtitle, palette, style),
    heatmap = .multilpa_plot_heatmap(x, main, subtitle, style),
    entropy = ,
    posteriors = .multilpa_plot_case_diagnostic(x, what, main, subtitle,
                                                palette, style))
  invisible(x)
}

#' Draw one case-level classification diagnostic, for any fitted family
#'
#' The entropy and assignment-probability panels need only the individual
#' posteriors and the effective profile counts, which every family of this
#' package carries. They are drawn from here rather than from each family's
#' `plot()` method so that `plot(diagnostics(fit))` can draw them for a fit
#' whose own method offers a different set of `what` values.
#'
#' @param x A fitted model carrying `subject_posteriors`.
#' @param what `"entropy"` or `"posteriors"`.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Fill colours, or `NULL` for the package palette.
#' @param style Visual constants.
#' @param ... Further style constants, merged into `style`.
#' @return `NULL`, invisibly. Called for the side effect of drawing.
#' @noRd
.multilpa_plot_case_diagnostic <- function(x, what = c("entropy", "posteriors"),
                                           main = NULL, subtitle = NULL,
                                           palette = NULL,
                                           style = .multilpa_style(), ...) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  what <- match.arg(what)
  .multilpa_require_case_posteriors(x, what)
  style <- utils::modifyList(style, list(...))
  previous <- graphics::par(no.readonly = TRUE)
  # See the note in plot.multilpa(): `mfg` switches `new` on when it is set, so
  # restoring it on a device nothing has been drawn to both warns and leaves
  # `new = TRUE` behind.
  previous$mfg <- NULL
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA)
  switch(what,
    entropy = .multilpa_plot_entropy(x, main, subtitle, palette, style),
    posteriors = .multilpa_plot_posteriors(x, main, subtitle, palette, style))
  invisible(NULL)
}

#' Refuse a case-level diagnostic panel that has nothing to separate
#'
#' Both panels lay the cases of each profile out side by side. With one profile
#' every case is assigned to it with probability one, so the panel would be a
#' single spike over a zero-width axis; the drawing code fails on the degenerate
#' break sequence, which reads as an internal error rather than as a refusal.
#'
#' @param x A fitted model of this package.
#' @param what The view that was asked for, for the message.
#' @return `NULL`, invisibly, when the view can be drawn.
#' @noRd
.multilpa_require_case_posteriors <- function(x, what) {
  posteriors <- x$subject_posteriors
  if (is.null(posteriors) || !is.matrix(posteriors)) {
    stop(errorCondition(sprintf(
      "This fit carries no individual posteriors, so `what = \"%s\"` has nothing to draw. Views available for this fit: %s.",
      what, .multilpa_available_views_text(x)),
      class = "multilpa_nothing_to_plot", call = NULL))
  }
  if (ncol(posteriors) < 2L) {
    stop(errorCondition(sprintf(
      paste("A single-profile fit assigns every case to that profile with probability one,",
            "so `what = \"%s\"` has nothing to separate. Views available for this fit: %s."),
      what, .multilpa_available_views_text(x)),
      class = "multilpa_nothing_to_plot", call = NULL))
  }
  invisible(NULL)
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
  vars <- names(blocks)
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
  values <- t(matrix(vapply(seq_along(blocks), function(index) {
    blocks[[index]][, chosen[[index]]]
  }, numeric(n_profiles)), n_profiles, length(blocks)))
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
  positions <- seq_along(vars)
  # The axis names each indicator with the category being plotted, so a mixed
  # set of category labels cannot be misread as one shared category.
  shared <- length(unique(category_labels)) == 1L
  axis_labels <- if (shared) vars else
    sprintf("%s=%s", vars, category_labels)
  .multilpa_panel(xlim = c(1 - 0.35, length(vars) + 0.35), ylim = c(0, 1.02),
    xlab = "Categorical indicator",
    ylab = if (shared) sprintf("P(response = %s)", category_labels[[1L]]) else
      "Response probability",
    main = if (is.null(main)) "Response probabilities by profile" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d profiles, %d group class%s, %d categorical indicator%s", n_profiles,
      x$n_group_classes, if (x$n_group_classes == 1L) "" else "es",
      length(vars), if (length(vars) == 1L) "" else "s") else subtitle,
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
    label_y <- .multilpa_spread_labels(values[length(vars), ], 0.055)
    graphics::text(length(vars) + .multilpa_label_offset(style), label_y,
                   label_text,
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
  vars <- .multilpa_continuous_names(x)
  n_profiles <- x$n_profiles
  means <- x$means
  if (length(vars) == 0L) {
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
    spread <- vapply(seq_along(vars), function(index) {
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
  positions <- seq_along(vars)
  span <- range(means)
  padding <- 0.12 * max(diff(span), .Machine$double.eps)
  ylim <- c(span[1L] - padding, span[2L] + padding)
  xlim <- c(1 - 0.35, length(vars) + 0.35)
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  label_text <- sprintf("Profile %d (%.0f%%)", seq_len(n_profiles), 100 * share)
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  default_main <- sprintf("Profile means across %d indicator%s",
                          length(vars),
                          if (length(vars) == 1L) "" else "s")
  default_subtitle <- sprintf("%d profiles, %d group class%s; %s scale",
    n_profiles, x$n_group_classes,
    if (x$n_group_classes == 1L) "" else "es",
    if (identical(scale, "standardized")) "standardized" else "input")
  .multilpa_panel(xlim = xlim, ylim = ylim, xlab = "Indicator",
                ylab = if (identical(scale, "standardized"))
                  "Standardized mean" else "Estimated mean",
                main = if (is.null(main)) default_main else main,
                subtitle = if (is.null(subtitle)) default_subtitle else subtitle,
                x_at = positions, x_labels = vars, style = style)
  if (identical(scale, "standardized")) {
    graphics::abline(h = 0, col = style$muted_colour, lwd = 1, lty = 3L)
  }
  # Point area carries the profile's share of the sample, so a reader sees how
  # much of the data each line speaks for without consulting a second table.
  sizes <- .multilpa_point_sizes(share, style)
  invisible(lapply(seq_len(n_profiles), function(profile) {
    values <- means[profile, ]
    graphics::lines(positions, values, col = colours[profile],
                    lwd = style$line_width, lty = lines[profile])
    graphics::points(positions, values, pch = points[profile],
                     bg = colours[profile], col = style$panel_fill,
                     cex = sizes[profile], lwd = 1.4)
  }))
  if (isTRUE(labels)) {
    label_y <- .multilpa_spread_labels(means[, length(vars)],
                                     0.055 * diff(ylim))
    graphics::text(length(vars) + .multilpa_label_offset(style), label_y,
                   label_text, adj = c(0, 0.5),
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
    graphics::text(n_profiles + .multilpa_label_offset(style), label_y,
                   label_text, adj = c(0, 0.5),
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
#'                                n_profiles = 1:3, n_group_classes = 1, n_starts = 2,
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
  # `mfg` is dropped before the state is restored: setting it switches `new` on
  # as a documented side effect, so restoring it on a device nothing has been
  # drawn to yet both warns and leaves `new = TRUE` behind. That is exactly the
  # state after a refused view, where the restore runs before anything is drawn.
  previous$mfg <- NULL
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
      graphics::text(max(profile_counts) + .multilpa_label_offset(style),
                     .multilpa_spread_labels(ends[keep], 0.055 * diff(span)),
                     label_text[keep],
                     adj = c(0, 0.5), col = colours[keep],
                     cex = style$label_text_size, font = 2L)
    }
  }
  invisible(x)
}

#' Draw the assignments in sequence order
#'
#' One row per group and one column per position, filled with the assigned
#' profile. Rows are blocked by latent group class and divided by a rule, so the
#' picture answers whether a class's groups look alike over time rather than
#' only how often each profile occurs in them. Each cell carries its profile
#' number as well as its colour, so the profile survives greyscale printing and
#' colour-blind vision; the numbers are dropped only where the grid is too dense
#' to hold them.
#'
#' @param x A fitted `multilpa` model made with `time =`.
#' @param labels Whether to label each class block at the right edge.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Colours, one per profile; `NULL` uses the Okabe-Ito palette.
#' @param style Visual constants, as built by `.multilpa_style()`.
#' @param cell_labels Whether to print the profile number inside each cell, so
#'   that colour is not the only channel carrying it.
#' @return `NULL`, invisibly. Called for the side effect of drawing.
#' @noRd
.multilpa_plot_sequences <- function(x, labels, main, subtitle, palette, style,
                                     cell_labels = TRUE) {
  ## sequences() owns the reshape, and the same helper supplies the rectangular
  ## layout here, so the picture cannot drift from the tidy verb behind it. Both
  ## are in increasing group order, so one ordering serves both.
  layout <- sequences(x, format = "wide")
  ordering <- order(layout$group_class, layout$group)
  codes <- .multilpa_sequence_matrix(x)[ordering, , drop = FALSE]
  classes <- layout$group_class[ordering]
  colours <- if (is.null(palette)) .multilpa_palette(x$n_profiles) else
    rep(palette, length.out = x$n_profiles)
  time_positions <- sort(unique(x$time_values))
  positions <- if (is.numeric(time_positions)) as.numeric(time_positions) else
    seq_along(time_positions)
  n_groups <- nrow(codes)
  class_labels <- sprintf("Class %d", sort(unique(classes)))
  ## The right side is sized for the class labels by the shared helper; the
  ## bottom gains room the line panels do not need, for the profile key.
  margins <- .multilpa_margins(style, if (isTRUE(labels)) class_labels else
                                 character(), style$label_text_size)
  if (isTRUE(labels)) margins[1L] <- margins[1L] + 3.2
  graphics::par(mar = margins)

  .multilpa_panel(
    xlim = c(min(positions) - 0.5, max(positions) + 0.5),
    ylim = c(n_groups + 0.5, 0.5),
    xlab = x$time, ylab = "Group",
    main = if (is.null(main)) "Assigned profile in sequence order" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d groups in %d classes; %d profiles across %d positions",
      n_groups, x$n_group_classes, x$n_profiles, ncol(codes)) else subtitle,
    x_at = positions, x_labels = as.character(time_positions), style = style)
  # image() indexes z positionally and ignores its dimnames, so the group and
  # occasion names carried by .multilpa_sequence_matrix() are dropped here rather
  # than handed to a callee that would silently discard them.
  graphics::image(x = positions, y = seq_len(n_groups), z = unname(t(codes)),
                  add = TRUE, col = colours, zlim = c(0.5, x$n_profiles + 0.5))
  .multilpa_cell_labels(codes, positions, colours, style, cell_labels)
  boundaries <- which(diff(classes) != 0) + 0.5
  if (length(boundaries)) {
    graphics::abline(h = boundaries, col = style$panel_fill, lwd = 3)
  }
  if (isTRUE(labels)) {
    graphics::text(max(positions) + 0.6,
                   vapply(split(seq_len(n_groups), classes), mean, numeric(1)),
                   class_labels, adj = c(0, 0.5),
                   col = style$text_colour, cex = style$label_text_size, font = 2L)
    graphics::legend("bottom", horiz = TRUE, bty = "n", inset = c(0, -0.30),
                     legend = sprintf("Profile %d", seq_len(x$n_profiles)),
                     fill = colours, border = NA, xpd = NA,
                     cex = style$label_text_size)
  }
  invisible(NULL)
}

#' Plot a covariate model
#'
#' Draws the parts of a covariate fit that are still fixed quantities. The
#' measurement model is one set of profile means as usual, and the assignments
#' can be laid out in sequence order. Profile prevalence cannot be drawn this
#' way: a covariate model has no single prevalence vector, because prevalence is
#' a function of each unit's covariates, so asking for it is refused rather than
#' answered with an average that no unit has.
#'
#' @param x A fitted `multilpa_covariates` model.
#' @param what `"profiles"` (the default) draws the measurement model;
#'   `"sequences"` draws the assignments in course order and needs a fit made
#'   with `time =`. `"entropy"` and `"posteriors"` draw the two case-level
#'   classification diagnostics exactly as [plot.multilpa()] draws them: they
#'   read the individual posteriors, which a covariate fit has, and say nothing
#'   about prevalence, which it does not.
#' @param scale,labels,cell_labels,main,subtitle,palette,symbols,linetypes,style,...
#'   Passed through as in [plot.multilpa()].
#' @return The fitted model, invisibly. Called for the side effect of drawing.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- fit_covariates(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                       n_group_classes = 2, profile_covariates = "x",
#'                       n_starts = 2, seed = 1)
#' plot(fit)
#' @export
plot.multilpa_covariates <- function(x, what = c("profiles", "sequences",
                                                 "entropy", "posteriors"),
                                     scale = c("raw", "standardized"),
                                     labels = TRUE, main = NULL, subtitle = NULL,
                                     palette = NULL, symbols = NULL,
                                     linetypes = NULL, style = .multilpa_style(),
                                     cell_labels = TRUE, ...) {
  stopifnot("`x` must be a fitted `multilpa_covariates` model" =
              inherits(x, "multilpa_covariates"),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`cell_labels` must be TRUE or FALSE" =
              isTRUE(cell_labels) || isFALSE(cell_labels))
  if (identical(what, "probabilities")) {
    stop(errorCondition(
      paste("A covariate model has no single profile prevalence: it varies with",
            "each unit's covariates. Use parameter_inference() for the logits."),
      class = "multilpa_nothing_to_plot", call = NULL))
  }
  what <- match.arg(what)
  scale <- match.arg(scale)
  style <- utils::modifyList(style, list(...))
  previous <- graphics::par(no.readonly = TRUE)
  # `mfg` is dropped before the state is restored: setting it switches `new` on
  # as a documented side effect, so restoring it on a device nothing has been
  # drawn to yet both warns and leaves `new = TRUE` behind. That is exactly the
  # state after a refused view, where the restore runs before anything is drawn.
  previous$mfg <- NULL
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  graphics::par(xpd = NA)
  switch(what,
    profiles = .multilpa_plot_profiles(x, scale, labels, main, subtitle, palette,
                                       symbols, linetypes, style),
    sequences = .multilpa_plot_sequences(x, labels, main, subtitle, palette,
                                         style, cell_labels),
    entropy = ,
    posteriors = .multilpa_plot_case_diagnostic(x, what, main, subtitle,
                                                palette, style))
  invisible(x)
}

#' Draw the measurement model as a heatmap of standardized profile means
#'
#' One tile per profile and indicator, filled by how far that profile's mean
#' sits from the indicator's grand mean in standard deviations, so indicators on
#' different scales are comparable within one panel. The value is printed in
#' each tile in ink chosen for contrast, which keeps the number readable without
#' a colour key and keeps the distinction from resting on colour alone.
#'
#' @param x A fitted model carrying `means` and `indicator_data`.
#' @param main,subtitle Panel title and secondary line.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_heatmap <- function(x, main, subtitle, style) {
  vars <- .multilpa_continuous_names(x)
  if (length(vars) == 0L) {
    stop(errorCondition(
      "This model has no continuous indicators; plot `what = \"responses\"` instead.",
      class = "multilpa_no_continuous", call = NULL))
  }
  # One standardization map for the whole package, so the heatmap, the bars and
  # as.data.frame(scale = "standardized") cannot disagree about what z means.
  standardization <- .multilpa_standardization(x, vars)
  z <- sweep(sweep(x$means, 2L, standardization$centre, "-"), 2L,
             standardization$spread, "/")
  n_profiles <- nrow(z)
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  rows <- rev(seq_len(n_profiles))
  graphics::par(mar = style$margins + c(0, 2.2, 0, 0))
  .multilpa_panel(
    xlim = c(0.5, length(vars) + 0.5), ylim = c(0.5, n_profiles + 0.5),
    xlab = "Indicator", ylab = "",
    main = if (is.null(main)) "Profile means in standard deviations" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "distance from each indicator's grand mean; %d profile%s", n_profiles,
      if (n_profiles == 1L) "" else "s") else subtitle,
    x_at = seq_along(vars), x_labels = vars,
    y_at = rows, y_labels = sprintf("Profile %d\n(%.0f%%)", seq_len(n_profiles),
                                    100 * share),
    style = style)
  cells <- expand.grid(indicator = seq_along(vars),
                       profile = seq_len(n_profiles))
  values <- z[cbind(cells$profile, cells$indicator)]
  fills <- .multilpa_diverging(values)
  graphics::rect(cells$indicator - 0.5, rows[cells$profile] - 0.5,
                 cells$indicator + 0.5, rows[cells$profile] + 0.5,
                 col = fills, border = style$panel_fill, lwd = 1.5)
  graphics::text(cells$indicator, rows[cells$profile], sprintf("%.2f", values),
                 col = .multilpa_ink(fills), cex = style$label_text_size)
  invisible(NULL)
}

#' Draw a ridge plot of a per-case quantity, one ridge per profile
#'
#' Shared by the entropy and posterior-probability diagnostics, which differ
#' only in what they measure and where their reference line sits.
#'
#' Each ridge is scaled to its own maximum, as `ggridges` does, so the shape of
#' a small profile stays readable. The cost is that height carries no
#' information about size, which is why every ridge is labelled with its share
#' and the subtitle says so.
#'
#' @param values Per-case numeric vector.
#' @param profile Integer profile index per case.
#' @param n_profiles Number of profiles.
#' @param share Profile shares, used to label each ridge.
#' @param limits Range of the horizontal axis.
#' @param reference Vertical reference line, or `NULL`.
#' @param reference_label Text for the reference line.
#' @param xlab,main,subtitle Axis label, title and secondary line.
#' @param palette Fill colours, or `NULL` for the package palette.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_case_ridges <- function(values, profile, n_profiles, share,
                                       limits, reference, reference_label,
                                       xlab, main, subtitle, palette, style) {
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  breaks <- seq(limits[1L], limits[2L], length.out = 22L)
  rows <- rev(seq_len(n_profiles))
  graphics::par(mar = style$margins + c(0, 2.2, 0, 0))
  .multilpa_panel(
    xlim = limits, ylim = c(0.5, n_profiles + 0.75), xlab = xlab, ylab = "",
    main = main, subtitle = subtitle,
    y_at = rows, y_labels = sprintf("Profile %d\n(%.0f%%)", seq_len(n_profiles),
                                    100 * share),
    style = style)
  invisible(lapply(seq_len(n_profiles), function(index) {
    .multilpa_ridge(values[profile == index], breaks = breaks,
                    base = rows[index] - 0.42, height = 0.88,
                    fill = colours[index], border = colours[index])
  }))
  if (!is.null(reference)) {
    # segments(), not abline(): plot.multilpa() sets xpd = NA so that direct
    # labels can sit outside the panel, and an abline() then runs the full
    # height of the device and straight through the title.
    top <- n_profiles + 0.62
    graphics::segments(reference, 0.5, reference, top, lty = 2L, lwd = 1.3,
                       col = style$muted_colour)
    graphics::text(reference, top, reference_label, pos = 4L, offset = 0.25,
                   cex = style$label_text_size * 0.9, col = style$muted_colour)
  }
  invisible(NULL)
}

#' Draw the per-case entropy contribution within each profile
#' @param x A fitted model carrying `subject_posteriors`.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Fill colours, or `NULL`.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_entropy <- function(x, main, subtitle, palette, style) {
  posteriors <- x$subject_posteriors
  n_profiles <- ncol(posteriors)
  # -sum p log p per case: zero when a case is assigned with certainty, and
  # log(K) when the posterior is flat.
  contribution <- -rowSums(posteriors * log(pmax(posteriors, .Machine$double.xmin)))
  assigned <- max.col(posteriors, ties.method = "first")
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  .multilpa_plot_case_ridges(
    values = contribution, profile = assigned, n_profiles = n_profiles,
    share = share, limits = c(0, log(n_profiles) * 1.02),
    reference = mean(contribution), reference_label = "mean",
    xlab = "Case-specific entropy contribution",
    main = if (is.null(main)) "How confidently is each case assigned?" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      paste("0 is certain, %.2f is a flat posterior over %d profiles;",
            "each ridge is scaled to its own maximum"),
      log(n_profiles), n_profiles) else subtitle,
    palette = palette, style = style)
}

#' Draw the posterior probability of the assigned profile, within each profile
#' @param x A fitted model carrying `subject_posteriors`.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Fill colours, or `NULL`.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_posteriors <- function(x, main, subtitle, palette, style) {
  posteriors <- x$subject_posteriors
  n_profiles <- ncol(posteriors)
  assigned <- max.col(posteriors, ties.method = "first")
  probability <- posteriors[cbind(seq_len(nrow(posteriors)), assigned)]
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  .multilpa_plot_case_ridges(
    values = probability, profile = assigned, n_profiles = n_profiles,
    share = share, limits = c(1 / n_profiles, 1),
    reference = NULL, reference_label = NULL,
    xlab = "Posterior probability of the assigned profile",
    main = if (is.null(main)) "Modal assignment probabilities" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      paste("mass piled at 1 is a clean separation; %.2f is a coin toss;",
            "each ridge is scaled to its own maximum"),
      1 / n_profiles) else subtitle,
    palette = palette, style = style)
}

#' The plots this package can draw
#'
#' Returns the catalogue of `what =` values accepted by the `plot()` methods,
#' with the group each belongs to, so the available views can be listed rather
#' than recalled from a help page.
#'
#' @return A base `data.frame` with one row per plot type and the columns
#'   `type` (the value to pass as `what`), `group` (`"measurement"`,
#'   `"structure"`, `"diagnostics"` or `"selection"`), and `description`.
#' @examples
#' multilpa_plot_types()
#' @export
multilpa_plot_types <- function() {
  data.frame(
    type = c("profiles", "bars", "heatmap", "responses", "probabilities",
             "sequences", "random_intercepts", "entropy", "posteriors",
             "enumeration"),
    group = c(rep("measurement", 4L), rep("structure", 3L),
              rep("diagnostics", 2L), "selection"),
    description = c(
      "Profile means across indicators, point size showing profile prevalence",
      "Profile means as grouped bars, with 95% intervals when `data` is given",
      "Profile means as standard deviations from each indicator's grand mean",
      "Categorical response probabilities, one line per profile",
      "Profile prevalence within each group class, the two-level quantity",
      "Each group's profile at each occasion, one row per group",
      "One interval per group: its posterior mean random intercept, plus and minus one posterior standard deviation",
      "Per-case entropy contribution within each profile, as ridges",
      "Posterior probability of the assigned profile, as ridges",
      "Information criteria across a candidate grid (plot an enumeration)"
    ),
    stringsAsFactors = FALSE
  )
}

#' Draw profile means as grouped bars
#'
#' The layout most readers of a profile analysis expect: indicators along the
#' axis, one bar per profile within each indicator. Position within a group
#' repeats what colour says -- the first bar is always the first profile -- and
#' each bar carries its own value, so neither the comparison nor the number
#' depends on reading a colour.
#'
#' On the standardized scale the bars run from zero in both directions, which is
#' the usual way a profile is read: above or below the average for that
#' indicator.
#'
#' @param x A fitted model carrying `means`.
#' @param scale `"raw"` or `"standardized"`.
#' @param errors Optional matrix of standard errors, drawn as whiskers.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Bar fills, or `NULL` for the package palette.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_bars <- function(x, scale, errors, main, subtitle, palette, style) {
  vars <- .multilpa_continuous_names(x)
  if (length(vars) == 0L) {
    stop(errorCondition(
      "This model has no continuous indicators; plot `what = \"responses\"` instead.",
      class = "multilpa_no_continuous", call = NULL))
  }
  n_profiles <- x$n_profiles
  means <- x$means
  if (identical(scale, "standardized")) {
    standardization <- .multilpa_standardization(x, vars)
    means <- sweep(sweep(means, 2L, standardization$centre, "-"), 2L,
                   standardization$spread, "/")
    if (!is.null(errors)) errors <- sweep(errors, 2L, standardization$spread, "/")
  }
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  width <- 0.8 / n_profiles
  offsets <- (seq_len(n_profiles) - (n_profiles + 1) / 2) * width
  centres <- seq_along(vars)
  top <- if (is.null(errors)) means else means + 1.96 * errors
  bottom <- if (is.null(errors)) means else means - 1.96 * errors
  baseline <- if (identical(scale, "standardized")) 0 else
    min(0, min(bottom) * 1.05)
  span <- range(c(baseline, top, bottom))
  padding <- 0.16 * max(diff(span), .Machine$double.eps)
  # Room under the axis title for the legend strip.
  graphics::par(mar = style$margins + c(2.6, 0, 0, 0))
  .multilpa_panel(
    xlim = c(0.4, length(vars) + 0.6),
    ylim = c(span[1L] - padding * 0.4, span[2L] + padding),
    xlab = "Indicator",
    ylab = if (identical(scale, "standardized")) "Standardized mean" else
      "Estimated mean",
    main = if (is.null(main)) "Profile means" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d profiles; %s scale%s", n_profiles,
      if (identical(scale, "standardized")) "standardized" else "input",
      if (is.null(errors)) "" else "; whiskers are 95% intervals") else subtitle,
    x_at = centres, x_labels = vars, style = style)
  graphics::abline(h = baseline, col = style$muted_colour, lwd = 1)
  invisible(lapply(seq_len(n_profiles), function(profile) {
    left <- centres + offsets[profile] - width * 0.44
    right <- centres + offsets[profile] + width * 0.44
    value <- means[profile, ]
    graphics::rect(left, baseline, right, value, col = colours[profile],
                   border = style$panel_fill, lwd = 1.2)
    if (!is.null(errors)) {
      middle <- (left + right) / 2
      graphics::segments(middle, bottom[profile, ], middle, top[profile, ],
                         col = style$text_colour, lwd = 1.2)
      graphics::segments(left + width * 0.12, top[profile, ],
                         right - width * 0.12, top[profile, ],
                         col = style$text_colour, lwd = 1.2)
      graphics::segments(left + width * 0.12, bottom[profile, ],
                         right - width * 0.12, bottom[profile, ],
                         col = style$text_colour, lwd = 1.2)
    }
    # Above a positive bar, below a negative one, so the label never sits on it.
    graphics::text((left + right) / 2,
                   ifelse(value >= baseline, top[profile, ], bottom[profile, ]),
                   sprintf("%.2f", value),
                   pos = ifelse(value >= baseline, 3L, 1L), offset = 0.3,
                   cex = style$label_text_size * 0.85, col = style$text_colour)
  }))
  .multilpa_legend_strip(sprintf("Profile %d (%.0f%%)", seq_len(n_profiles),
                                 100 * share), colours, style)
  invisible(NULL)
}

#' Measurement standard errors shaped like the means matrix
#'
#' @param x A fitted model.
#' @param data The data frame the model was fitted to, or `NULL` to use the
#'   columns the fit carries.
#' @return A profiles-by-indicators matrix of standard errors, or `NULL`.
#' @noRd
.multilpa_mean_error_matrix <- function(x, data) {
  # A fit carries the columns it was built from, so the intervals cost the
  # caller nothing to ask for and are drawn by default. Families whose standard
  # errors are not implemented refuse by condition class; the bars are still
  # worth drawing without whiskers, so the refusal is caught and reported as
  # "no errors" rather than propagated out of a plot call.
  errors <- tryCatch(
    .multilpa_measurement_errors(x, .multilpa_resolve_data(x, data)),
    multilpa_unsupported_inference = function(condition) NULL,
    multilpa_no_inference = function(condition) NULL,
    multilpa_singular_information = function(condition) NULL,
    multilpa_incomplete_fit = function(condition) NULL,
    multilpa_bad_inference_data = function(condition) NULL)
  if (is.null(errors)) return(NULL)
  vars <- .multilpa_continuous_names(x)
  cells <- expand.grid(indicator = vars, profile = seq_len(x$n_profiles),
                       stringsAsFactors = FALSE)
  values <- .multilpa_match_error(errors, "mean", cells$profile, cells$indicator)
  if (all(is.na(values))) return(NULL)
  matrix(values, x$n_profiles, length(vars), byrow = TRUE)
}
