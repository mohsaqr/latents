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
#' @param what Which view to draw. [plot_views()] lists every value
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
#'   `"raincloud"` shows what the means summarize: one panel per indicator,
#'   and in it, for the cases assigned to each profile, a density of the
#'   indicator, a box of its quartiles and median, and the observations
#'   themselves. It shows the spread within each profile and how far the
#'   profiles overlap. For a fit
#'   whose indicators are all categorical, `"heatmap"` draws the response
#'   probabilities instead: one row per profile and one column per category
#'   (a binary indicator shows its last category), each cell printing its
#'   probability.
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
#'   How big each profile is: `"sizes"` draws one bar per profile holding its
#'   effective count -- the posterior mass it carries, not the number of cases
#'   that won a tie -- with the count and the share printed on the bar.
#'
#'   Classification quality: `"entropy"` draws each case's entropy contribution
#'   as one ridge per profile, and `"posteriors"` draws the posterior
#'   probability of each case's assigned profile the same way. Both ridges are
#'   scaled to their own maximum, following the usual ridgeline convention, so
#'   ridge height compares shapes and not profile sizes; prevalence is printed
#'   in each profile's label instead. Both read the individual posteriors alone,
#'   so every family of this package can draw them; a one-profile fit refuses
#'   them with an error of class `latents_nothing_to_plot`, because every case
#'   then belongs to the single profile with probability one. `"avepp"` draws
#'   the average posterior probability matrix: one row per assigned profile, one
#'   column per profile, each cell the mean posterior that group puts on that
#'   profile. The diagonal is the avePP usually reported, and the off-diagonal
#'   says which profiles a group is confused with, which the diagonal alone
#'   cannot show.
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
#' @param intervals For `what = "profiles"` and `"responses"`, `TRUE` draws a
#'   95% interval on every mean or probability when the fit has standard
#'   errors (a probability's interval is clipped to `[0, 1]`); the profiles are then dodged
#'   apart so the whiskers stay readable. A fit without them (a bound-active or
#'   unconverged fit, or a family without inference) is drawn without whiskers
#'   and its subtitle says why.
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
#' @details `"profiles"` and `"bars"` draw 95% Wald intervals when the fit has
#'   standard errors; the other views show point estimates. Profile order is
#'   arbitrary, so two fits must have their labels aligned before their plots
#'   are compared. With `scale = "standardized"` the
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
#' plot(fit, what = "sizes")
#' plot(fit, what = "avepp")
#' plot_views()
#' @seealso [plot_views()] for the catalogue of views.
#' @export
plot.multilpa <- function(x, what = c("profiles", "bars", "heatmap", "raincloud",
                                      "responses",
                                      "probabilities", "sequences", "sizes",
                                      "entropy", "posteriors", "avepp", "all"),
                        data = NULL,
                        scale = c("raw", "standardized"), category = "last",
                        labels = TRUE, main = NULL, subtitle = NULL,
                        palette = NULL, symbols = NULL, linetypes = NULL,
                        style = .multilpa_style(), cell_labels = TRUE,
                        intervals = TRUE, ...) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"),
            "`intervals` must be TRUE or FALSE" = isTRUE(intervals) || isFALSE(intervals),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`cell_labels` must be TRUE or FALSE" =
              isTRUE(cell_labels) || isFALSE(cell_labels),
            "`category` must be a single label or index" = length(category) == 1L)
  what <- match.arg(what)
  if (identical(what, "all")) {
    return(.multilpa_plot_every_view(x, match.call(), parent.frame()))
  }
  if (what %in% c("entropy", "posteriors", "avepp")) {
    .multilpa_refuse_noise(x, sprintf("plot(what = \"%s\")", what),
                           class = "latents_nothing_to_plot")
  }
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
                                     symbols, linetypes, style,
                                     errors = if (isTRUE(intervals))
                                       .multilpa_mean_error_matrix(x, data)),
    responses = .multilpa_plot_responses(x, category, labels, main, subtitle,
                                       palette, symbols, linetypes, style,
                                       errors = if (isTRUE(intervals))
                                         .multilpa_response_error_table(x, data)),
    probabilities = .multilpa_plot_probabilities(x, labels, main, subtitle,
                                               palette, symbols, linetypes, style),
    sequences = .multilpa_plot_sequences(x, labels, main, subtitle, palette,
                                        style, cell_labels),
    bars = .multilpa_plot_bars(x, scale, .multilpa_mean_error_matrix(x, data),
                               main, subtitle, palette, style),
    heatmap = if (length(.multilpa_continuous_names(x)) == 0L) {
      .multilpa_plot_response_heatmap(x, main, subtitle, style)
    } else .multilpa_plot_heatmap(x, main, subtitle, style),
    raincloud = .multilpa_plot_raincloud(
      x, scale, main = main %||% "Indicator distributions by profile",
      subtitle = subtitle %||% paste(
        "density, quartiles and observations of each indicator among the",
        "cases assigned to each profile"),
      palette, symbols, style),
    sizes = .multilpa_plot_sizes(x, main, subtitle, palette, style),
    avepp = .multilpa_plot_avepp(x, main, subtitle, style),
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
      class = "latents_nothing_to_plot", call = NULL))
  }
  if (ncol(posteriors) < 2L) {
    stop(errorCondition(sprintf(
      paste("A single-profile fit assigns every case to that profile with probability one,",
            "so `what = \"%s\"` has nothing to separate. Views available for this fit: %s."),
      what, .multilpa_available_views_text(x)),
      class = "latents_nothing_to_plot", call = NULL))
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
                                   symbols, linetypes, style, errors = NULL) {
  missing_reason <- if (is.character(errors)) errors else NULL
  if (is.character(errors)) errors <- NULL
  blocks <- x$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) {
    stop(errorCondition("This model has no categorical indicators.",
                        class = "latents_no_categorical", call = NULL))
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
                          class = "latents_unknown_category", call = NULL))
    }
    index
  }, integer(1))
  values <- t(matrix(vapply(seq_along(blocks), function(index) {
    blocks[[index]][, chosen[[index]]]
  }, numeric(n_profiles)), n_profiles, length(blocks)))
  category_labels <- vapply(seq_along(blocks), function(index) {
    colnames(blocks[[index]])[chosen[[index]]]
  }, character(1))
  # One standard error per plotted probability, profiles in columns like
  # `values`; NA where the fit reports none.
  spread <- if (is.null(errors)) NULL else {
    t(matrix(vapply(seq_along(blocks), function(index) {
      .multilpa_match_error(errors, "response", seq_len(n_profiles),
                            rep(sprintf("%s:%s", vars[[index]],
                                        category_labels[[index]]), n_profiles))
    }, numeric(n_profiles)), n_profiles, length(blocks)))
  }
  held_any <- isTRUE(attr(errors, "held"))
  if (!is.null(spread) && all(is.na(spread))) spread <- NULL
  dodge <- if (is.null(spread)) numeric(n_profiles) else
    (seq_len(n_profiles) - (n_profiles + 1) / 2) * min(0.08, 0.3 / n_profiles)
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
    subtitle = if (is.null(subtitle)) paste0(sprintf(
      "%d profiles%s, %d categorical indicator%s", n_profiles,
      if (x$n_group_classes == 1L) "" else sprintf(", %d group classes",
                                                    x$n_group_classes),
      length(vars), if (length(vars) == 1L) "" else "s"),
      if (!is.null(spread)) {
        if (held_any) "; whiskers: 95% intervals where estimable" else
          "; whiskers: 95% intervals"
      } else
        if (!is.null(missing_reason)) paste0("; ", missing_reason) else "") else
      subtitle,
    x_at = positions, x_labels = rep("", length(positions)),
    y_at = seq(0, 1, by = 0.25), style = style)
  .multilpa_block_names(axis_labels, positions, 1, style,
                        line = style$axis_line + 0.4)
  invisible(lapply(seq_len(n_profiles), function(profile) {
    at <- positions + dodge[profile]
    graphics::lines(at, values[, profile], col = colours[profile],
                    lwd = style$line_width, lty = lines[profile])
    if (!is.null(spread)) {
      # A probability's interval is clipped to [0, 1], where it lives.
      low <- pmax(0, values[, profile] - 1.96 * spread[, profile])
      high <- pmin(1, values[, profile] + 1.96 * spread[, profile])
      shown <- !is.na(low)
      graphics::segments(at[shown], low[shown], at[shown], high[shown],
                         col = colours[profile], lwd = 1.6)
    }
    graphics::points(at, values[, profile], pch = points[profile],
                     bg = colours[profile], col = style$panel_fill,
                     cex = style$point_size, lwd = 1.4)
  }))
  if (isTRUE(labels)) {
    label_y <- .multilpa_spread_labels(values[length(vars), ],
                                       max(0.055, .multilpa_label_gap(style)))
    graphics::text(length(vars) + max(dodge) + .multilpa_label_offset(style), label_y,
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
                                  symbols, linetypes, style, errors = NULL) {
  vars <- .multilpa_continuous_names(x)
  n_profiles <- x$n_profiles
  means <- x$means
  if (length(vars) == 0L) {
    stop(errorCondition("This model has no continuous indicators; plot `what = \"responses\"` instead.",
                        class = "latents_no_continuous", call = NULL))
  }
  missing_reason <- if (is.character(errors)) errors else NULL
  if (is.character(errors)) errors <- NULL
  if (identical(scale, "standardized")) {
    standardization <- .multilpa_standardization(x, vars)
    rescale <- function(values) {
      sweep(sweep(values, 2L, standardization$centre, "-"), 2L,
            standardization$spread, "/")
    }
    means <- rescale(means)
    if (!is.null(errors)) errors <- sweep(errors, 2L, standardization$spread, "/")
  }
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  points <- if (is.null(symbols)) .multilpa_symbols(n_profiles) else
    rep(symbols, length.out = n_profiles)
  lines <- if (is.null(linetypes)) .multilpa_linetypes(n_profiles) else
    rep(linetypes, length.out = n_profiles)
  # With intervals the profiles are dodged apart, so overlapping whiskers stay
  # readable; without them the lines meet at the indicator itself.
  dodge <- if (is.null(errors)) numeric(n_profiles) else
    (seq_len(n_profiles) - (n_profiles + 1) / 2) * min(0.08, 0.3 / n_profiles)
  positions <- seq_along(vars)
  top <- if (is.null(errors)) means else means + 1.96 * errors
  bottom <- if (is.null(errors)) means else means - 1.96 * errors
  span <- range(c(top, bottom), na.rm = TRUE)
  padding <- 0.12 * max(diff(span), .Machine$double.eps)
  ylim <- c(span[1L] - padding, span[2L] + padding)
  xlim <- c(1 - 0.35, length(vars) + 0.35)
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  share <- share[seq_len(n_profiles)]
  label_text <- sprintf("Profile %d (%.0f%%)", seq_len(n_profiles), 100 * share)
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  default_main <- sprintf("Profile means across %d indicator%s",
                          length(vars),
                          if (length(vars) == 1L) "" else "s")
  default_subtitle <- paste0(sprintf("%d profiles%s; %s scale",
    n_profiles,
    if (x$n_group_classes == 1L) "" else sprintf(", %d group classes",
                                                  x$n_group_classes),
    if (identical(scale, "standardized")) "standardized" else "input"),
    if (!is.null(errors)) "; whiskers are 95% intervals" else
      if (!is.null(missing_reason)) paste0("; ", missing_reason) else "")
  .multilpa_panel(xlim = xlim, ylim = ylim, xlab = "Indicator",
                ylab = if (identical(scale, "standardized"))
                  "Standardized mean" else "Estimated mean",
                main = if (is.null(main)) default_main else main,
                subtitle = if (is.null(subtitle)) default_subtitle else subtitle,
                x_at = positions, x_labels = rep("", length(vars)), style = style)
  .multilpa_block_names(vars, positions, 1, style, line = style$axis_line + 0.4)
  if (identical(scale, "standardized")) {
    .multilpa_hline(0, col = style$muted_colour, lwd = 1, lty = 3L)
  }
  # Point area carries the profile's share of the sample, so a reader sees how
  # much of the data each line speaks for without consulting a second table.
  sizes <- .multilpa_point_sizes(share, style)
  invisible(lapply(seq_len(n_profiles), function(profile) {
    at <- positions + dodge[profile]
    values <- means[profile, ]
    graphics::lines(at, values, col = colours[profile],
                    lwd = style$line_width, lty = lines[profile])
    if (!is.null(errors)) {
      graphics::segments(at, bottom[profile, ], at, top[profile, ],
                         col = colours[profile], lwd = 1.6)
      graphics::segments(at - 0.035, top[profile, ], at + 0.035, top[profile, ],
                         col = colours[profile], lwd = 1.6)
      graphics::segments(at - 0.035, bottom[profile, ], at + 0.035,
                         bottom[profile, ], col = colours[profile], lwd = 1.6)
    }
    graphics::points(at, values, pch = points[profile],
                     bg = colours[profile], col = style$panel_fill,
                     cex = sizes[profile], lwd = 1.4)
  }))
  if (isTRUE(labels)) {
    label_y <- .multilpa_spread_labels(
      means[, length(vars)],
      max(0.055 * diff(ylim), .multilpa_label_gap(style)))
    graphics::text(length(vars) + max(dodge) + .multilpa_label_offset(style),
                   label_y, label_text, adj = c(0, 0.5),
                   col = colours, cex = style$label_text_size, font = 2L)
  }
  invisible(NULL)
}

#' The observations behind a profile plot, with each case's assigned profile
#'
#' The indicators are taken as the fit stored them, which are the values it was
#' fitted to. A centred fit reports its means on the centred scale, where the
#' stored columns would not sit, so it is refused rather than drawn wrongly.
#' @param x A fitted model.
#' @param vars The continuous indicators.
#' @return A list with `values` (cases by indicators) and `profile`.
#' @noRd
.multilpa_plot_observations <- function(x, vars) {
  if (!identical(x$centering %||% "none", "none")) {
    stop(errorCondition(paste(
      "The raincloud view draws the indicators as supplied, and this fit was",
      "centred, so its means are on a different scale."),
      class = "latents_bad_argument", call = NULL))
  }
  values <- x$indicator_data
  if (is.null(values) || is.null(x$subject_profiles)) {
    stop(errorCondition("This fit did not retain its indicators, so they cannot be drawn.",
                        class = "latents_no_indicator_data", call = NULL))
  }
  values <- matrix(as.numeric(values), nrow(values), ncol(values),
                   dimnames = list(NULL, colnames(values)))
  # A noise component's cases belong to no profile, so they are not drawn.
  kept <- x$subject_profiles >= 1L & x$subject_profiles <= x$n_profiles
  list(values = values[kept, vars, drop = FALSE],
       profile = x$subject_profiles[kept])
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
#' Draws information criteria against the number of profiles as line plots,
#' one panel per criterion, with one line per number of group classes and
#' covariance model. Candidates that failed to converge are marked on the
#' baseline, so a gap in a line reads as a failure and not as a missing
#' candidate.
#'
#' @param x An `multilpa_enumeration` result from [enumerate_classes()].
#' @param criterion One or more criterion columns of `as.data.frame(x)`, such
#'   as `"bic_individual"` or `"sabic_groups"`. Several names draw one panel
#'   each; the default draws AIC, BIC under both sample-size conventions, and
#'   ICL counted over individuals.
#' @param combine `TRUE`, the default, draws several criteria as panels of one
#'   figure. `FALSE` draws each criterion as its own full-size figure, one after
#'   another, so a report shows them as separate images.
#' @param labels `TRUE` prints a direct label at the right end of each series.
#' @param mark_minimum `TRUE` rings the lowest value among converged candidates.
#'   This marks an extremum, it does not select a model.
#' @param main,subtitle Panel title and secondary line.
#' @param palette,symbols,linetypes Series aesthetics, recycled over combinations
#'   of group-class count and covariance model.
#' @param style A list of visual constants, as built by `.multilpa_style()`.
#' @param ... Further named visual constants, merged into `style`.
#' @return The enumeration result, invisibly. Called for its drawing side effect.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' candidates <- enumerate_classes(
#'   example_data, c("score_a", "score_b"), "school", n_profiles = 1:3,
#'   n_group_classes = 1, n_starts = 2, seed = 1
#' )
#' plot(candidates)
#' @export
plot.multilpa_enumeration <- function(x, criterion = c("aic", "bic_groups",
                                                      "bic_individual",
                                                      "icl_individual"),
                                    combine = TRUE,
                                    labels = TRUE, mark_minimum = TRUE,
                                    main = NULL, subtitle = NULL, palette = NULL,
                                    symbols = NULL, linetypes = NULL,
                                    style = .multilpa_style(), ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"),
            "`criterion` must name one or more columns" =
              is.character(criterion) && length(criterion) >= 1L &&
              !anyNA(criterion),
            "`combine` must be TRUE or FALSE" = isTRUE(combine) || isFALSE(combine),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`mark_minimum` must be TRUE or FALSE" =
              isTRUE(mark_minimum) || isFALSE(mark_minimum))
  grid <- as.data.frame(x)
  unknown <- setdiff(criterion, .multilpa_enumeration_criteria())
  if (length(unknown) > 0L) {
    stop(errorCondition(sprintf(
      "%s is not an information criterion in the enumeration grid.",
      paste(sprintf("`%s`", unknown), collapse = ", ")),
      class = "latents_unknown_criterion", call = NULL))
  }
  empty <- criterion[!vapply(criterion, function(name) {
    any(grid$converged %in% TRUE & is.finite(grid[[name]]))
  }, logical(1))]
  if (length(empty) > 0L) {
    stop(errorCondition(sprintf(
      "No converged candidate has a finite %s; nothing to plot.",
      paste(sprintf("`%s`", empty), collapse = ", ")),
      class = "latents_nothing_to_plot", call = NULL))
  }
  style <- utils::modifyList(style, list(...))
  criterion <- .multilpa_distinct_criteria(grid, criterion)
  previous <- graphics::par(no.readonly = TRUE)
  # `mfg` is dropped before the state is restored: setting it switches `new` on
  # as a documented side effect, so restoring it on a device nothing has been
  # drawn to yet both warns and leaves `new = TRUE` behind.
  previous$mfg <- NULL
  on.exit(graphics::par(previous), add = TRUE, after = FALSE)
  several <- length(criterion) > 1L
  if (several && combine) {
    # One panel per criterion, laid out once here and restored once on exit,
    # so the panels do not reset each other's layout.
    columns <- if (length(criterion) <= 3L) length(criterion) else 2L
    graphics::par(mfrow = c(ceiling(length(criterion) / columns), columns))
  }
  invisible(lapply(seq_along(criterion), function(index) {
    name <- criterion[[index]]
    title <- names(criterion)[[index]]
    .multilpa_enumeration_panel(
      grid, name, labels = labels, mark_minimum = mark_minimum,
      main = if (is.null(main)) (if (several && combine) title else NULL) else main,
      subtitle = if (is.null(subtitle) && several && combine) "" else subtitle,
      palette = palette, symbols = symbols, linetypes = linetypes,
      style = style, title = title, axis_title = !(several && combine))
  }))
  invisible(x)
}

#' The criteria worth a panel, named for display
#'
#' A single-level grid has one unit per observation, so every criterion's
#' group-level and individual-level versions are the same numbers; drawing both
#' repeats a panel. Such a pair is kept once under the criterion's plain name.
#' @param grid The enumeration table.
#' @param criterion Requested criterion columns.
#' @return The criteria to draw, named by their display titles.
#' @noRd
.multilpa_distinct_criteria <- function(grid, criterion) {
  base <- sub("_(groups|individual)$", "", criterion)
  level <- ifelse(grepl("_groups$", criterion), "groups",
                  ifelse(grepl("_individual$", criterion), "individuals", ""))
  same_at_both <- vapply(base, function(stem) {
    both <- paste0(stem, c("_groups", "_individual"))
    all(both %in% names(grid)) &&
      isTRUE(all.equal(grid[[both[1L]]], grid[[both[2L]]]))
  }, logical(1))
  keep <- !duplicated(ifelse(same_at_both, base, criterion))
  titles <- ifelse(same_at_both | !nzchar(level), toupper(base),
                   sprintf("%s (%s)", toupper(base), level))
  stats::setNames(criterion[keep], titles[keep])
}

#' Draw one information criterion across an enumeration grid
#'
#' The panel [plot.multilpa_enumeration()] draws for each criterion. It sets
#' margins but leaves the graphics state to its caller, which saves it once
#' before the first panel and restores it once after the last.
#' @param grid The enumeration table.
#' @param criterion One criterion column name, already validated.
#' @param labels,mark_minimum,main,subtitle,palette,symbols,linetypes,style As
#'   for [plot.multilpa_enumeration()].
#' @return `NULL`, invisibly. Called for its drawing side effect.
#' @noRd
.multilpa_enumeration_panel <- function(grid, criterion, labels, mark_minimum,
                                        main, subtitle, palette, symbols,
                                        linetypes, style, title = criterion,
                                        axis_title = TRUE) {
  values <- grid[[criterion]]
  eligible <- grid$converged %in% TRUE & is.finite(values)
  if (!any(eligible)) {
    stop(errorCondition(sprintf("No converged candidate has a finite `%s`; nothing to plot.",
                                criterion),
                        class = "latents_nothing_to_plot", call = NULL))
  }
  graphics::par(xpd = NA)
  series_model <- grid$model
  known_models <- unique(stats::na.omit(series_model))
  if (length(known_models) == 1L) {
    series_model[is.na(series_model)] <- known_models
  }
  series_keys <- unique(data.frame(n_group_classes = grid$n_group_classes,
                                   model = series_model))
  series_keys <- series_keys[order(series_keys$n_group_classes,
                                   series_keys$model), , drop = FALSE]
  n_series <- nrow(series_keys)
  profile_counts <- sort(unique(grid$n_profiles))
  colours <- if (is.null(palette)) .multilpa_palette(n_series) else
    rep(palette, length.out = n_series)
  points <- if (is.null(symbols)) .multilpa_symbols(n_series) else
    rep(symbols, length.out = n_series)
  lines <- if (is.null(linetypes)) .multilpa_linetypes(n_series) else
    rep(linetypes, length.out = n_series)
  span <- range(values[eligible])
  padding <- 0.12 * max(diff(span), .Machine$double.eps)
  xlim <- c(min(profile_counts) - 0.35, max(profile_counts) + 0.35)
  failures <- sum(!grid$converged)
  # A label names only what tells the series apart: the covariance model, the
  # group-class count, or both; "1 group class" on every line says nothing.
  group_text <- sprintf("%d group class%s", series_keys$n_group_classes,
                        ifelse(series_keys$n_group_classes == 1L, "", "es"))
  several_models <- length(unique(series_keys$model)) > 1L
  several_groups <- length(unique(series_keys$n_group_classes)) > 1L
  label_text <- if (several_models && several_groups) {
    paste(series_keys$model, group_text, sep = ", ")
  } else if (several_models) {
    as.character(series_keys$model)
  } else if (several_groups) {
    group_text
  } else ""  # one series needs no label
  graphics::par(mar = .multilpa_margins(style, if (isTRUE(labels)) label_text else
    character(), style$label_text_size))
  .multilpa_panel(xlim = xlim, ylim = c(span[1L] - padding, span[2L] + padding),
    xlab = "Number of profiles", ylab = if (isTRUE(axis_title)) title else "",
    main = if (is.null(main)) sprintf("%s by candidate model", title) else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "%d candidates%s; lower is preferred, no candidate is selected automatically",
      nrow(grid),
      if (failures == 0L) "" else sprintf(", %d did not converge", failures)) else
      subtitle,
    x_at = profile_counts, x_labels = profile_counts, style = style)
  ends <- vapply(seq_len(n_series), function(index) {
    same_model <- if (is.na(series_keys$model[index]))
      is.na(series_model) else series_model == series_keys$model[index]
    rows <- grid$n_group_classes == series_keys$n_group_classes[index] &
      same_model
    series <- data.frame(profiles = grid$n_profiles[rows],
                         value = values[rows], eligible = eligible[rows])
    series <- series[order(series$profiles), , drop = FALSE]
    finite <- series$eligible
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
      offset <- (index - (n_series + 1) / 2) * 0.09
      graphics::points(series$profiles[!finite] + offset,
                       rep(span[1L] - padding * 0.55, sum(!finite)),
                       pch = 4L, col = colours[index], cex = style$point_size,
                       lwd = 2)
    }
    if (any(finite)) series$value[finite][sum(finite)] else NA_real_
  }, numeric(1))
  if (isTRUE(mark_minimum)) {
    best <- which.min(replace(values, !eligible, Inf))
    graphics::points(grid$n_profiles[best], values[best], pch = 1L,
                     col = style$title_colour, cex = style$point_size * 2.4,
                     lwd = 1.8)
  }
  if (isTRUE(labels)) {
    keep <- !is.na(ends)
    if (any(keep)) {
      graphics::text(max(profile_counts) + .multilpa_label_offset(style),
                     .multilpa_spread_labels(ends[keep],
                       max(0.055 * diff(span), .multilpa_label_gap(style))),
                     label_text[keep],
                     adj = c(0, 0.5), col = colours[keep],
                     cex = style$label_text_size, font = 2L)
    }
  }
  invisible(NULL)
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
  layout <- .multilpa_sequences(x, format = "wide")
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
    .multilpa_hline(boundaries, col = style$panel_fill, lwd = 3)
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
#' @param what `"profiles"` (the default) draws the measurement model, and
#'   `"bars"` and `"heatmap"` draw it as [plot.multilpa()] does (for an
#'   all-categorical fit the heatmap shows the response probabilities);
#'   `"sizes"` and `"avepp"` draw the effective profile sizes and the average
#'   posterior matrix. `"sequences"` draws the assignments in course order and
#'   needs a fit made with `time =`. `"entropy"` and `"posteriors"` draw the two case-level
#'   classification diagnostics exactly as [plot.multilpa()] draws them: they
#'   read the individual posteriors, which a covariate fit has, and say nothing
#'   about prevalence, which it does not. `"responses"` draws the categorical
#'   response probabilities, one line per profile, as [plot.multilpa()] does;
#'   the measurement model does not depend on the covariates.
#' @param data,scale,labels,cell_labels,main,subtitle,palette,symbols,linetypes,style,category,intervals,...
#'   As in [plot.multilpa()].
#' @return The fitted model, invisibly. Called for the side effect of drawing.
#' @examples
#' set.seed(5)
#' school <- rep(seq_len(16), each = 8)
#' high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
#' x <- rnorm(128)
#' profile <- ifelse(
#'   runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
#' )
#' example_data <- data.frame(
#'   school = school, x = x,
#'   y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
#'   y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
#' )
#' fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
#'                 n_group_classes = 2, profile_covariates = "x",
#'                 n_starts = 2, seed = 1)
#' plot(fit)
#' @export
plot.multilpa_covariates <- function(x, what = c("profiles", "bars", "heatmap",
                                                 "raincloud",
                                                 "responses", "sequences",
                                                 "sizes", "entropy",
                                                 "posteriors", "avepp", "all"),
                                     data = NULL,
                                     scale = c("raw", "standardized"),
                                     category = "last",
                                     labels = TRUE, main = NULL, subtitle = NULL,
                                     palette = NULL, symbols = NULL,
                                     linetypes = NULL, style = .multilpa_style(),
                                     cell_labels = TRUE, intervals = TRUE,
                                     ...) {
  stopifnot("`x` must be a fitted `multilpa_covariates` model" =
              inherits(x, "multilpa_covariates"),
            "`intervals` must be TRUE or FALSE" = isTRUE(intervals) || isFALSE(intervals),
            "`labels` must be TRUE or FALSE" = isTRUE(labels) || isFALSE(labels),
            "`cell_labels` must be TRUE or FALSE" =
              isTRUE(cell_labels) || isFALSE(cell_labels))
  if (identical(what, "probabilities")) {
    stop(errorCondition(
      paste("A covariate model has no single profile prevalence: it varies with",
            "each unit's covariates. Use parameter_inference() for the logits."),
      class = "latents_nothing_to_plot", call = NULL))
  }
  what <- match.arg(what)
  if (identical(what, "all")) {
    return(.multilpa_plot_every_view(x, match.call(), parent.frame()))
  }
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
  # The measurement and classification views are the ones plot.multilpa()
  # draws: the measurement model does not depend on the covariates, and the
  # classification views read only the posteriors.
  switch(what,
    profiles = .multilpa_plot_profiles(x, scale, labels, main, subtitle, palette,
                                       symbols, linetypes, style,
                                       errors = if (isTRUE(intervals))
                                         .multilpa_mean_error_matrix(x, data)),
    bars = .multilpa_plot_bars(x, scale, .multilpa_mean_error_matrix(x, data),
                               main, subtitle, palette, style),
    heatmap = if (length(.multilpa_continuous_names(x)) == 0L) {
      .multilpa_plot_response_heatmap(x, main, subtitle, style)
    } else .multilpa_plot_heatmap(x, main, subtitle, style),
    raincloud = .multilpa_plot_raincloud(
      x, scale, main = main %||% "Indicator distributions by profile",
      subtitle = subtitle %||% paste(
        "density, quartiles and observations of each indicator among the",
        "cases assigned to each profile"),
      palette, symbols, style),
    sequences = .multilpa_plot_sequences(x, labels, main, subtitle, palette,
                                         style, cell_labels),
    sizes = .multilpa_plot_sizes(x, main, subtitle, palette, style),
    avepp = .multilpa_plot_avepp(x, main, subtitle, style),
    entropy = ,
    posteriors = .multilpa_plot_case_diagnostic(x, what, main, subtitle,
                                                palette, style),
    responses = .multilpa_plot_responses(x, category, labels, main, subtitle,
                                         palette, symbols, linetypes, style,
                                         errors = if (isTRUE(intervals))
                                           .multilpa_response_error_table(x, data)))
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
      class = "latents_no_continuous", call = NULL))
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

#' Draw the average posterior probability matrix
#'
#' One tile per pair of profiles: the row is the profile a case was assigned to,
#' the column is a profile, and the cell is the mean posterior that the cases
#' assigned to that row put on that column. The diagonal is the classification
#' certainty usually reported as avePP, and the off-diagonal says where a
#' profile leaks -- which of the two is the problem is not visible from the
#' diagonal alone.
#'
#' The fill is the shared white-to-blue ramp taken over the whole probability
#' range, so darker is always a higher posterior wherever it sits. Every cell
#' also prints its number, so nothing rests on reading a colour. A profile no
#' case was assigned to has no row of means, and prints empty rather than `NA`.
#'
#' @param x A fitted model carrying `subject_posteriors`.
#' @param main,subtitle Panel title and secondary line.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_avepp <- function(x, main, subtitle, style) {
  probabilities <- x$subject_posteriors
  if (is.null(probabilities) || ncol(probabilities) < 1L) {
    stop(errorCondition(
      "This model carries no individual posteriors to average.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  averages <- .multilpa_average_posterior_matrix(probabilities)
  n_classes <- ncol(averages)
  assigned <- tabulate(max.col(probabilities, ties.method = "first"), n_classes)
  rows <- rev(seq_len(n_classes))
  graphics::par(mar = style$margins + c(0, 2.2, 0, 0))
  .multilpa_panel(
    xlim = c(0.5, n_classes + 0.5), ylim = c(0.5, n_classes + 0.5),
    xlab = "Posterior profile", ylab = "",
    main = if (is.null(main)) "Average posterior probability" else main,
    subtitle = if (is.null(subtitle)) paste(
      "mean posterior each assigned group puts on each profile;",
      "the diagonal is classification certainty") else subtitle,
    x_at = seq_len(n_classes),
    x_labels = sprintf("Profile %d", seq_len(n_classes)),
    y_at = rows,
    y_labels = sprintf("Assigned %d\n(n = %d)", seq_len(n_classes), assigned),
    style = style)
  cells <- expand.grid(column = seq_len(n_classes), row = seq_len(n_classes))
  values <- averages[cbind(cells$row, cells$column)]
  # A probability needs no midpoint: taken over [0, 1] the shared ramp uses its
  # upper half only, which is a plain white-to-blue sequential scale.
  fills <- .multilpa_diverging(values, limit = 1)
  graphics::rect(cells$column - 0.5, rows[cells$row] - 0.5,
                 cells$column + 0.5, rows[cells$row] + 0.5,
                 col = fills, border = style$panel_fill, lwd = 1.5)
  labels <- rep("", length(values))
  labels[!is.na(values)] <- sprintf("%.2f", values[!is.na(values)])
  graphics::text(cells$column, rows[cells$row], labels,
                 col = .multilpa_ink(fills), cex = style$label_text_size)
  invisible(NULL)
}

#' Draw how many cases each profile holds
#'
#' Effective counts rather than modal ones, because a profile's size in a
#' mixture is the posterior mass it carries and not the number of cases that
#' happened to win a tie. Each bar prints its own count and share, so the
#' comparison does not depend on reading a colour or a gridline.
#'
#' @param x A fitted model carrying `effective_profile_counts`.
#' @param main,subtitle Panel title and secondary line.
#' @param palette Bar fills, or `NULL` for the package palette.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_sizes <- function(x, main, subtitle, palette, style) {
  counts <- x$effective_profile_counts
  if (is.null(counts) || length(counts) == 0L) {
    stop(errorCondition(
      "This model carries no effective profile counts to draw.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  n_profiles <- length(counts)
  share <- counts / sum(counts)
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  graphics::par(mar = style$margins)
  .multilpa_panel(
    xlim = c(0.5, n_profiles + 0.5), ylim = c(0, max(counts) * 1.18),
    xlab = "Profile", ylab = "Effective count",
    main = if (is.null(main)) "Profile sizes" else main,
    subtitle = if (is.null(subtitle)) sprintf(
      "posterior mass per profile; %d case%s in total", round(sum(counts)),
      if (round(sum(counts)) == 1L) "" else "s") else subtitle,
    x_at = seq_len(n_profiles),
    x_labels = sprintf("Profile %d", seq_len(n_profiles)),
    style = style)
  graphics::rect(seq_len(n_profiles) - 0.36, 0,
                 seq_len(n_profiles) + 0.36, counts,
                 col = colours, border = style$panel_fill, lwd = 1.5)
  graphics::text(seq_len(n_profiles), counts,
                 sprintf("%.0f\n(%.1f%%)", counts, 100 * share),
                 pos = 3L, offset = 0.4, cex = style$label_text_size)
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
#' plot_views()
#' @export
plot_views <- function() {
  data.frame(
    type = c("profiles", "bars", "heatmap", "raincloud", "responses",
             "probabilities",
             "sequences", "transitions", "sizes", "entropy", "posteriors",
             "avepp", "enumeration", "all"),
    group = c(rep("measurement", 5L), rep("structure", 4L),
              rep("diagnostics", 3L), "selection", "every"),
    description = c(
      "Profile means across indicators, point size showing profile prevalence",
      "Profile means as grouped bars, with 95% intervals when `data` is given",
      "Profile means as standard deviations from each indicator's grand mean",
      "Each indicator's distribution by assigned profile: density, box, observations",
      "Categorical response probabilities, one line per profile",
      "Profile prevalence within each group class, the two-level quantity",
      "Each group's profile at each occasion, one row per group",
      "Estimated transition matrix, one panel per group class (a transition fit)",
      "Effective number of cases in each profile, with its share",
      "Per-case entropy contribution within each profile, as ridges",
      "Posterior probability of the assigned profile, as ridges",
      "Average posterior probability: assigned profile by posterior profile",
      "Information criteria across a candidate grid (plot an enumeration)",
      "Every view above that this fit has the ingredients for, in one call"
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
      class = "latents_no_continuous", call = NULL))
  }
  n_profiles <- x$n_profiles
  means <- x$means
  missing_reason <- if (is.character(errors)) errors else NULL
  if (is.character(errors)) errors <- NULL
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
      if (!is.null(errors)) "; whiskers are 95% intervals" else
        if (!is.null(missing_reason)) paste0("; ", missing_reason) else "") else
      subtitle,
    x_at = centres, x_labels = rep("", length(vars)), style = style)
  .multilpa_block_names(vars, centres, 1, style, line = style$axis_line + 0.4)
  .multilpa_hline(baseline, col = style$muted_colour, lwd = 1)
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
#' @return A profiles-by-indicators matrix of standard errors; or a single
#'   string saying why there are none, for the plot's subtitle.
#' @noRd
.multilpa_mean_error_matrix <- function(x, data) {
  # A fit carries the columns it was built from, so the intervals cost the
  # caller nothing to ask for and are drawn by default. Families whose standard
  # errors are not implemented refuse by condition class; the bars are still
  # worth drawing without whiskers, so the refusal is caught and reported as
  # "no errors" rather than propagated out of a plot call.
  # The reason travels to the subtitle, so a plot without whiskers says why
  # instead of looking like a plot whose intervals were forgotten.
  because <- function(reason) function(condition) reason
  errors <- tryCatch(
    .multilpa_measurement_errors(x, .multilpa_resolve_data(x, data)),
    latents_unsupported_inference = because("no intervals for this model"),
    latents_unsupported_noise = because("no intervals with a noise component"),
    latents_singular_information = because("no intervals: singular information"),
    latents_incomplete_fit = because("no intervals: incomplete fit"),
    latents_bad_inference_data = because("no intervals: data do not reproduce the fit"),
    latents_boundary_fit = because("no intervals: an estimate is at its bound"),
    latents_no_converge = because("no intervals: the fit did not converge"))
  if (is.character(errors)) return(errors)
  vars <- .multilpa_continuous_names(x)
  cells <- expand.grid(indicator = vars, profile = seq_len(x$n_profiles),
                       stringsAsFactors = FALSE)
  values <- .multilpa_match_error(errors, "mean", cells$profile, cells$indicator)
  if (all(is.na(values))) return("no intervals for this model")
  matrix(values, x$n_profiles, length(vars), byrow = TRUE)
}

#' Response-probability standard errors for the responses plot
#' @param x A fitted model.
#' @param data The fitting data, or `NULL` for the columns the fit carries.
#' @return The measurement rows of the inference table; or a string saying
#'   why there are none.
#' @noRd
.multilpa_response_error_table <- function(x, data) {
  because <- function(reason) function(condition) reason
  # Probabilities on their bound are held there, so every other one keeps its
  # interval instead of one bound-active probability removing them all.
  tryCatch(
    {
      inference <- parameter_inference(x, .multilpa_resolve_data(x, data),
                                       boundary = "fix")
      inference <- inference[inference$level == "measurement", , drop = FALSE]
      attr(inference, "held") <- length(attr(inference, "fixed_at_bound") %||%
                                          character()) > 0L
      inference
    },
    latents_unsupported_inference = because("no intervals for this model"),
    latents_unsupported_noise = because("no intervals with a noise component"),
    latents_singular_information = because("no intervals: singular information"),
    latents_incomplete_fit = because("no intervals: incomplete fit"),
    latents_bad_inference_data = because("no intervals: data do not reproduce the fit"),
    latents_boundary_fit = because("no intervals: a probability is at its bound"),
    latents_no_converge = because("no intervals: the fit did not converge"))
}

#' Draw every response probability as a heatmap
#'
#' The categorical counterpart of the profile-means heatmap, for a fit whose
#' indicators are all categorical: one row per profile and one column per
#' category of every indicator, each cell the probability of that response in
#' that profile. A binary indicator gets one column, its last category, since
#' the other is its complement; a polytomous one gets a column per category,
#' so a profile's whole response distribution reads across the row. The fill
#' is the white-to-blue probability ramp shared with the average posterior
#' matrix, and every cell prints its value.
#'
#' @param x A fitted model with categorical indicators only.
#' @param main,subtitle Panel title and secondary line.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_response_heatmap <- function(x, main, subtitle, style) {
  blocks <- x$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) {
    stop(errorCondition("This model has no indicators to draw.",
                        class = "latents_nothing_to_plot", call = NULL))
  }
  n_profiles <- x$n_profiles
  columns <- do.call(rbind, lapply(names(blocks), function(indicator) {
    block <- blocks[[indicator]]
    kept <- if (ncol(block) == 2L) ncol(block) else seq_len(ncol(block))
    data.frame(indicator = indicator, category = colnames(block)[kept],
               index = kept, stringsAsFactors = FALSE)
  }))
  columns$position <- seq_len(nrow(columns))
  binary <- all(vapply(blocks, ncol, integer(1)) == 2L)
  # A polytomous indicator's columns are labelled by category position, and
  # the indicator is named once under its block; full category labels are
  # usually too long to sit under a column. When every indicator shares one
  # set of categories the subtitle keys the positions to them.
  column_labels <- if (binary) columns$indicator else columns$index
  category_sets <- unique(lapply(blocks, colnames))
  key <- if (!binary && length(category_sets) == 1L) {
    categories <- category_sets[[1L]]
    # Labels that already begin with their position ("1 Extremely well") are
    # their own key.
    numbered <- all(startsWith(categories, as.character(seq_along(categories))))
    paste(if (numbered) categories else
      sprintf("%d = %s", seq_along(categories), categories), collapse = ", ")
  } else NULL
  share <- x$effective_profile_counts / sum(x$effective_profile_counts)
  rows <- rev(seq_len(n_profiles))
  graphics::par(mar = style$margins + c(if (binary) 1 else 1.6, 2.2,
                                       if (is.null(key)) 0 else 1, 0))
  .multilpa_panel(
    xlim = c(0.5, nrow(columns) + 0.5), ylim = c(0.5, n_profiles + 0.5),
    xlab = if (binary) "Indicator" else "", ylab = "",
    main = if (is.null(main)) "Response probabilities" else main,
    subtitle = if (is.null(subtitle)) paste0(sprintf(
      "probability of each response in each profile%s; %d profile%s",
      if (binary) sprintf(" (category %s shown)",
                          paste(unique(columns$category), collapse = "/")) else "",
      n_profiles, if (n_profiles == 1L) "" else "s"),
      if (is.null(key)) "" else paste0("\ncategories: ", key)) else subtitle,
    # A binary indicator's name goes under its column through the same
    # placement as a polytomous block's, which shrinks and staggers names that
    # do not fit rather than letting the axis drop every other one.
    x_at = columns$position,
    x_labels = if (binary) rep("", nrow(columns)) else column_labels,
    y_at = rows, y_labels = sprintf("Profile %d\n(%.0f%%)", seq_len(n_profiles),
                                    100 * share[seq_len(n_profiles)]),
    style = style)
  cells <- expand.grid(column = columns$position, profile = seq_len(n_profiles))
  values <- vapply(seq_len(nrow(cells)), function(cell) {
    column <- columns[cells$column[cell], ]
    blocks[[column$indicator]][cells$profile[cell], column$index]
  }, numeric(1))
  fills <- .multilpa_diverging(values, limit = 1)
  graphics::rect(cells$column - 0.5, rows[cells$profile] - 0.5,
                 cells$column + 0.5, rows[cells$profile] + 0.5,
                 col = fills, border = style$panel_fill, lwd = 1.5)
  # Values are printed only where they fit inside a cell; the fill still
  # carries them when the grid is too dense.
  value_size <- style$label_text_size * 0.9
  if (graphics::strwidth("0.00", cex = value_size) < 0.9) {
    graphics::text(cells$column, rows[cells$profile], sprintf("%.2f", values),
                   col = .multilpa_ink(fills), cex = value_size)
  }
  # A rule between indicators, so a polytomous item's categories read as one
  # block, and the indicator's name centred under that block.
  block <- match(columns$indicator, unique(columns$indicator))
  if (!binary) {
    edges <- which(diff(block) != 0L) + 0.5
    graphics::segments(edges, 0.5, edges, n_profiles + 0.5,
                       col = style$text_colour, lwd = 1.2)
  }
  .multilpa_block_names(unique(columns$indicator),
                        vapply(split(columns$position, block), mean, numeric(1)),
                        min(tabulate(block)), style,
                        line = style$axis_line + if (binary) 0.4 else 2.3,
                        font = if (binary) 1L else 2L)
  invisible(NULL)
}

#' Names under blocks of columns, never running together
#'
#' Base graphics silently drops axis labels that would overlap. Names that do
#' not fit the width they label are shrunk, and if still too wide alternate
#' between two lines, so every name is drawn and none touches its neighbour.
#' @param names The names, one per block.
#' @param centres Horizontal centre of each block, in user units.
#' @param width Width of the narrowest block, in user units.
#' @param style Visual constants.
#' @param line Margin line of the first row of names.
#' @param font Font face.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_block_names <- function(names, centres, width, style, line, font = 1L) {
  size <- style$axis_size
  widest <- max(graphics::strwidth(names, cex = size, font = font))
  if (widest > 0.95 * width) size <- max(0.6, size * 0.95 * width / widest)
  staggered <- max(graphics::strwidth(names, cex = size, font = font)) > 0.95 * width
  lines_at <- line + if (staggered) (seq_along(names) + 1L) %% 2L * 0.9 else 0
  graphics::mtext(names, side = 1L, at = centres, line = lines_at,
                  col = style$text_colour, cex = size, font = font)
  invisible(NULL)
}

#' Draw every view a fit supports, in one call
#'
#' `what = "all"` re-issues the caller's own `plot()` call once per supported
#' view, so every style argument the caller gave is carried into each panel
#' rather than silently reset to the defaults. A view can still refuse for a
#' reason only its own panel knows -- nothing in it to draw -- so a refusal
#' moves on to the next view and is named at the end instead of stopping the
#' sequence halfway through.
#'
#' @param x A fitted model of this package.
#' @param call The caller's `match.call()`, re-evaluated once per view.
#' @param env The caller's `parent.frame()`, where that call is evaluated.
#' @return `x`, invisibly. Called for the side effect of drawing.
#' @noRd
.multilpa_plot_every_view <- function(x, call, env) {
  # `match.call()` inside an S3 method names the method, and the methods of
  # this package are registered rather than exported, so re-evaluating that
  # call in the caller's environment fails for everyone who has not attached
  # the namespace. Dispatch through the generic instead.
  call[[1L]] <- quote(plot)
  views <- .multilpa_supported_views(x)
  if (length(views) == 0L || anyNA(views)) {
    stop(errorCondition(
      "This model family has no named plot views to draw.",
      class = "latents_no_plot", call = NULL))
  }
  drawn <- vapply(views, function(view) {
    call$what <- view
    tryCatch({
      eval(call, env)
      TRUE
    },
    latents_no_plot = function(condition) FALSE,
    latents_nothing_to_plot = function(condition) FALSE,
    latents_no_time = function(condition) FALSE,
    latents_no_categorical = function(condition) FALSE,
    latents_no_continuous = function(condition) FALSE,
    latents_no_indicator_data = function(condition) FALSE)
  }, logical(1))
  if (!all(drawn)) {
    message(sprintf("Not drawn for this model: %s.",
                    paste(views[!drawn], collapse = ", ")))
  }
  invisible(x)
}

#' Draw each indicator's distribution by profile as rainclouds
#'
#' One panel per continuous indicator. Within a panel each profile has a row,
#' and the row shows the indicator's distribution among the observations
#' assigned to that profile three ways: a density (the cloud), a box of the
#' quartiles with the median, and the observations themselves (the rain).
#' Where the profile plot shows the estimated means, this view shows the
#' spread and shape the means summarize, and how far the profiles overlap.
#' Profiles are told apart by row, colour and symbol together. The rain is
#' spread by a fixed function of the row rather than a random jitter, so
#' drawing leaves the random-number stream untouched.
#'
#' @param x A fitted model with continuous indicators.
#' @param scale `"raw"` or `"standardized"`.
#' @param main,subtitle Title and secondary line of the whole figure.
#' @param palette,symbols Per-profile colours and symbols, or `NULL`.
#' @param style Visual constants.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_plot_raincloud <- function(x, scale, main, subtitle, palette, symbols,
                                     style) {
  vars <- .multilpa_continuous_names(x)
  if (length(vars) == 0L) {
    stop(errorCondition(
      "This model has no continuous indicators; plot `what = \"heatmap\"` instead.",
      class = "latents_no_continuous", call = NULL))
  }
  observed <- .multilpa_plot_observations(x, vars)
  values <- observed$values
  if (identical(scale, "standardized")) {
    standardization <- .multilpa_standardization(x, vars)
    values <- sweep(sweep(values, 2L, standardization$centre, "-"), 2L,
                    standardization$spread, "/")
  }
  n_profiles <- x$n_profiles
  colours <- if (is.null(palette)) .multilpa_palette(n_profiles) else
    rep(palette, length.out = n_profiles)
  points <- if (is.null(symbols)) .multilpa_symbols(n_profiles) else
    rep(symbols, length.out = n_profiles)
  share <- tabulate(observed$profile, nbins = n_profiles) / length(observed$profile)
  row_labels <- sprintf("Profile %d\n(%.0f%%)", seq_len(n_profiles), 100 * share)
  columns <- if (length(vars) <= 3L) length(vars) else
    if (length(vars) == 4L) 2L else 3L
  graphics::par(mfrow = c(ceiling(length(vars) / columns), columns),
                oma = c(0.5, 0, if (is.null(subtitle) && is.null(main)) 0 else 3.2, 0),
                mar = c(3.2, 5.2, 2.2, 0.8))
  n_rain <- nrow(values)
  alpha <- max(0.1, min(0.6, 12 / sqrt(n_rain)))
  offset <- ((seq_len(n_rain) * 0.6180339887) %% 1)
  invisible(lapply(seq_along(vars), function(index) {
    column <- values[, index]
    xlim <- range(column, na.rm = TRUE)
    xlim <- xlim + c(-1, 1) * 0.04 * max(diff(xlim), .Machine$double.eps)
    rows <- rev(seq_len(n_profiles))
    .multilpa_panel(xlim = xlim, ylim = c(0.45, n_profiles + 0.55),
                    xlab = if (identical(scale, "standardized")) "Standardized value" else "",
                    ylab = "", main = vars[[index]], subtitle = NULL,
                    y_at = rows, y_labels = row_labels, style = style)
    invisible(lapply(seq_len(n_profiles), function(profile) {
      base <- rows[profile]
      own <- column[observed$profile == profile]
      own <- own[!is.na(own)]
      if (length(own) == 0L) return(NULL)
      # Cloud: a density above the row, scaled to its own maximum.
      if (length(unique(own)) > 1L) {
        density <- stats::density(own, from = min(own), to = max(own), n = 256L)
        height <- 0.42 * density$y / max(density$y)
        graphics::polygon(c(density$x, rev(density$x)),
                          c(base + 0.02 + height, rep(base + 0.02, length(height))),
                          col = grDevices::adjustcolor(colours[profile], alpha.f = 0.55),
                          border = colours[profile], lwd = 1)
      }
      # Box: quartiles, median and whiskers to the most extreme observation
      # within 1.5 interquartile ranges.
      quartiles <- stats::quantile(own, c(0.25, 0.5, 0.75), names = FALSE)
      reach <- 1.5 * diff(quartiles[c(1L, 3L)])
      whiskers <- c(min(own[own >= quartiles[1L] - reach]),
                    max(own[own <= quartiles[3L] + reach]))
      graphics::segments(whiskers[1L], base - 0.06, whiskers[2L], base - 0.06,
                         col = style$text_colour, lwd = 1)
      graphics::rect(quartiles[1L], base - 0.11, quartiles[3L], base - 0.01,
                     col = style$panel_fill, border = style$text_colour, lwd = 1)
      graphics::segments(quartiles[2L], base - 0.11, quartiles[2L], base - 0.01,
                         col = style$title_colour, lwd = 2)
      # Rain: the observations below the box.
      band <- offset[observed$profile == profile][seq_along(own)]
      graphics::points(own, base - 0.16 - 0.26 * band, pch = points[profile],
                       col = grDevices::adjustcolor(colours[profile], alpha.f = alpha),
                       bg = grDevices::adjustcolor(colours[profile], alpha.f = alpha * 0.6),
                       cex = 0.5)
    }))
  }))
  if (!is.null(main) || !is.null(subtitle)) {
    graphics::mtext(main %||% "", side = 3L, line = 1.7, outer = TRUE, adj = 0.02,
                    col = style$title_colour, cex = style$title_size, font = 2L)
    graphics::mtext(subtitle %||% "", side = 3L, line = 0.5, outer = TRUE,
                    adj = 0.02, col = style$muted_colour, cex = style$subtitle_size)
  }
  invisible(NULL)
}
