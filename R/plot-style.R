#' Colour-blind safe qualitative palette
#'
#' The Okabe-Ito palette, recycled when a model has more classes than colours.
#' Colour is never the only channel carrying a distinction in this package's
#' plots; point symbol and line type repeat it.
#'
#' @param n Number of colours required.
#' @return Character vector of `n` hexadecimal colours.
#' @noRd
.multilpa_palette <- function(n) {
  stopifnot("`n` must be a single positive integer" =
              is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1)
  colours <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
  rep(colours, length.out = n)
}

#' Point symbols distinguishing classes without relying on colour
#' @param n Number of symbols required.
#' @return Integer vector of `n` plotting characters.
#' @noRd
.multilpa_symbols <- function(n) {
  stopifnot("`n` must be a single positive integer" =
              is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1)
  rep(c(21L, 22L, 24L, 23L, 25L, 21L, 22L, 24L, 23L), length.out = n)
}

#' Line types distinguishing classes without relying on colour
#' @param n Number of line types required.
#' @return Integer vector of `n` line types.
#' @noRd
.multilpa_linetypes <- function(n) {
  stopifnot("`n` must be a single positive integer" =
              is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1)
  rep(c(1L, 2L, 4L, 3L, 5L, 6L), length.out = n)
}

#' Ink that stays legible on a given background
#'
#' Chooses black or white text for each background colour by its relative
#' luminance, so a label drawn inside a filled cell is readable whether the fill
#' is the palette's yellow or its black.
#'
#' @param background Character vector of background colours.
#' @return Character vector of the same length, each `"#000000"` or `"#FFFFFF"`.
#' @noRd
.multilpa_ink <- function(background) {
  stopifnot("`background` must be a character vector of colours" =
              is.character(background))
  if (length(background) == 0L) return(character())
  channels <- grDevices::col2rgb(background) / 255
  # Relative luminance, ITU-R BT.709 coefficients; the linearisation is skipped
  # because only the side of the threshold matters here, not the exact value.
  luminance <- as.vector(c(0.2126, 0.7152, 0.0722) %*% channels)
  ifelse(luminance > 0.55, "#000000", "#FFFFFF")
}

#' Print the code carried by each cell of a filled grid
#'
#' A grid of coloured cells would otherwise encode its series by colour alone.
#' This stamps each cell with the integer it stands for, in ink chosen to stay
#' legible on that cell, but only when the cell is large enough to hold the
#' text: the fit is measured rather than guessed, because a grid of overlapping
#' digits reads worse than the colours on their own.
#'
#' @param codes Integer matrix of codes, one row per grid row and one column per
#'   grid column, `NA` where a cell is empty.
#' @param positions Numeric positions of the grid columns on the horizontal axis.
#' @param colours The fill colour of each code, in code order.
#' @param style Visual constants, supplying the text size.
#' @param draw `FALSE` skips the labels entirely.
#' @return `TRUE` if the labels were drawn and `FALSE` if they were not,
#'   invisibly.
#' @noRd
.multilpa_cell_labels <- function(codes, positions, colours, style,
                                  draw = TRUE) {
  stopifnot("`codes` must be a matrix of codes" = is.matrix(codes),
            "`draw` must be TRUE or FALSE" = isTRUE(draw) || isFALSE(draw))
  if (!isTRUE(draw)) return(invisible(FALSE))
  filled <- !is.na(codes)
  if (!any(filled)) return(invisible(FALSE))
  text <- as.character(codes[filled])
  size <- style$label_text_size
  cell_width <- if (length(positions) > 1L) min(diff(positions)) else 1
  widest <- max(graphics::strwidth(text, units = "user", cex = size, font = 2L))
  tallest <- max(abs(graphics::strheight(text, units = "user", cex = size,
                                         font = 2L)))
  if (widest > 0.7 * cell_width || tallest > 0.7) return(invisible(FALSE))
  # which() and matrix indexing both walk the matrix by column, so the cells and
  # their codes stay aligned without a join.
  cells <- which(filled, arr.ind = TRUE)
  graphics::text(positions[cells[, 2L]], cells[, 1L], text,
                 col = .multilpa_ink(colours[codes[filled]]), cex = size,
                 font = 2L, adj = c(0.5, 0.5))
  invisible(TRUE)
}

#' Open a styled plotting panel
#'
#' Draws the background, horizontal reference grid and axes shared by every
#' plot in the package, so that no plot inherits the default base graphics
#' frame. The caller restores graphical parameters.
#'
#' @param xlim,ylim Numeric ranges of the data region.
#' @param xlab,ylab Axis titles.
#' @param main Panel title, or `NULL` for none.
#' @param subtitle Secondary line under the title, or `NULL` for none.
#' @param x_at,x_labels Positions and labels for the horizontal axis. When
#'   `x_at` is `NULL` a numeric axis is drawn.
#' @param y_at Positions for the vertical axis, or `NULL` for a numeric axis.
#' @param style A list of visual constants, as built by `.multilpa_style()`.
#' @return `NULL`, invisibly. Called for the side effect of drawing.
#' @noRd
.multilpa_panel <- function(xlim, ylim, xlab, ylab, main = NULL, subtitle = NULL,
                          x_at = NULL, x_labels = NULL, y_at = NULL,
                          y_labels = NULL, style = .multilpa_style()) {
  stopifnot("`style` must be a list of visual constants" = is.list(style))
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  usr <- graphics::par("usr")
  graphics::rect(usr[1L], usr[3L], usr[2L], usr[4L], col = style$panel_fill,
                 border = NA)
  grid_at <- if (is.null(y_at)) pretty(ylim) else y_at
  # Keep the grid and its labels inside the panel; pretty() overshoots the range.
  grid_at <- grid_at[grid_at >= ylim[1L] & grid_at <= ylim[2L]]
  # A categorical y axis names its own rows, so the reference grid would be
  # drawing lines through the middle of them.
  if (is.null(y_labels)) {
    graphics::abline(h = grid_at, col = style$grid_colour, lwd = style$grid_width)
  }
  axis_positions <- if (is.null(x_at)) pretty(xlim) else x_at
  axis_labels <- if (is.null(x_labels)) TRUE else x_labels
  graphics::axis(1L, at = axis_positions, labels = axis_labels, tick = FALSE,
                 line = style$axis_line, col.axis = style$text_colour,
                 cex.axis = style$axis_size, gap.axis = 0.25)
  graphics::axis(2L, at = if (is.null(y_labels)) grid_at else y_at,
                 labels = if (is.null(y_labels)) TRUE else y_labels,
                 tick = FALSE, las = 1L,
                 line = style$axis_line, col.axis = style$text_colour,
                 cex.axis = style$axis_size)
  graphics::mtext(xlab, side = 1L, line = style$label_line, col = style$text_colour,
                  cex = style$label_size)
  graphics::mtext(ylab, side = 2L, line = style$label_line + 0.8,
                  col = style$text_colour, cex = style$label_size)
  if (!is.null(main)) {
    graphics::mtext(main, side = 3L, line = if (is.null(subtitle)) 0.9 else 1.7,
                    adj = 0, col = style$title_colour, cex = style$title_size,
                    font = 2L)
  }
  if (!is.null(subtitle)) {
    graphics::mtext(subtitle, side = 3L, line = 0.5, adj = 0,
                    col = style$muted_colour, cex = style$subtitle_size)
  }
  invisible(NULL)
}

#' Visual constants shared by the package's plots
#'
#' Every constant the plots use is collected here and can be overridden, so no
#' colour, size or spacing is hard-coded inside a plotting method.
#'
#' @param panel_fill Panel background colour.
#' @param grid_colour,grid_width Reference grid colour and line width.
#' @param text_colour,title_colour,muted_colour Text colours for axis labels,
#'   titles, and secondary annotation.
#' @param axis_size,label_size,title_size,subtitle_size,label_text_size
#'   Character expansion for axis text, axis titles, titles, subtitles, and
#'   direct series labels.
#' @param axis_line,label_line Axis text and axis title offsets, in lines.
#' @param point_size,line_width Series point expansion and line width.
#' @param margins Plot margins in lines, as `c(bottom, left, top, right)`.
#' @param ... Further named constants, merged into the result.
#' @return A named list of visual constants.
#' @noRd
.multilpa_style <- function(panel_fill = "#FBFBFA", grid_colour = "#E4E4E1",
                          grid_width = 0.8, text_colour = "#3B3B3B",
                          title_colour = "#1A1A1A", muted_colour = "#6E6E6E",
                          axis_size = 0.82, label_size = 0.88, title_size = 1.0,
                          subtitle_size = 0.78, label_text_size = 0.78,
                          axis_line = -0.4, label_line = 2.1,
                          point_size = 1.25, line_width = 2.1,
                          margins = c(4.1, 4.6, 3.4, 5.6), ...) {
  c(list(panel_fill = panel_fill, grid_colour = grid_colour,
         grid_width = grid_width, text_colour = text_colour,
         title_colour = title_colour, muted_colour = muted_colour,
         axis_size = axis_size, label_size = label_size,
         title_size = title_size, subtitle_size = subtitle_size,
         label_text_size = label_text_size, axis_line = axis_line,
         label_line = label_line, point_size = point_size,
         line_width = line_width, margins = margins), list(...))
}

#' Place direct series labels without overlap
#'
#' Shifts labels apart along the vertical axis by a minimum separation, keeping
#' their original order. Direct labels remove the need for a legend, so the
#' reader never has to match a colour swatch to a series.
#'
#' @param y Numeric vector of preferred label positions.
#' @param minimum_gap Smallest permitted separation between adjacent labels.
#' @return Numeric vector of adjusted positions, in the input order.
#' @noRd
.multilpa_spread_labels <- function(y, minimum_gap) {
  stopifnot("`y` must be numeric" = is.numeric(y),
            "`minimum_gap` must be a single nonnegative number" =
              is.numeric(minimum_gap) && length(minimum_gap) == 1L &&
              is.finite(minimum_gap) && minimum_gap >= 0)
  if (length(y) < 2L) return(y)
  order_index <- order(y)
  sorted <- y[order_index]
  # One upward sweep is enough: each label is pushed just clear of the one
  # below it, so the result stays ordered and moves as little as possible.
  adjusted <- Reduce(function(previous, current) {
    max(current, previous + minimum_gap)
  }, sorted[-1L], accumulate = TRUE, init = sorted[1L])
  centre_shift <- (mean(sorted) - mean(adjusted))
  result <- numeric(length(y))
  result[order_index] <- adjusted + centre_shift
  result
}

#' Right margin needed for direct series labels
#'
#' Measures the widest label and converts it to margin lines, so direct labels
#' are never clipped by a fixed margin guess.
#'
#' @param labels Character vector of label text, possibly empty.
#' @param cex Character expansion the labels will be drawn at.
#' @param padding Extra margin lines added beyond the measured width.
#' @return A single number of margin lines.
#' @noRd
.multilpa_label_offset <- function(style) {
  # One "m" of the label's own type, expressed in user units. A fixed user-unit
  # gap such as 0.45 is a different physical distance on every panel, so the
  # margin reserved in inches fitted some plots and clipped others.
  graphics::strwidth("m", cex = style$label_text_size)
}

#' Right margin wide enough for the direct labels, in lines
#' @param labels Character vector of direct labels.
#' @param cex Character expansion the labels are drawn at.
#' @param padding Extra lines, covering the gap between panel and label.
#' @return A single number of margin lines.
#' @noRd
.multilpa_label_margin <- function(labels, cex, padding = 2.2) {
  stopifnot("`labels` must be character" = is.character(labels),
            "`cex` must be a single positive number" =
              is.numeric(cex) && length(cex) == 1L && is.finite(cex) && cex > 0)
  if (length(labels) == 0L) return(padding)
  line_height <- graphics::par("csi")
  ## Direct labels are drawn bold, so they must be measured bold; measuring in
  ## the regular weight under-reserves the margin and clips the last character.
  width <- max(graphics::strwidth(labels, units = "inches", cex = cex, font = 2L))
  width / line_height + padding
}

#' Plot margins with a right side sized for direct labels
#' @param style Visual constants, supplying the base margins.
#' @param labels Character vector of direct labels, possibly empty.
#' @param cex Character expansion the labels will be drawn at.
#' @return A numeric vector of four margin widths, in lines.
#' @noRd
.multilpa_margins <- function(style, labels, cex) {
  stopifnot("`style` must carry four margins" =
              is.list(style) && length(style$margins) == 4L)
  margins <- style$margins
  margins[4L] <- max(margins[4L], .multilpa_label_margin(labels, cex))
  margins
}

#' Point sizes that encode a share without becoming illegible
#'
#' The mclust-book profile plot encodes each profile's mixing proportion in the
#' point size, so a reader sees how much of the sample a line speaks for without
#' consulting a second table. Areas, not radii, are made proportional to the
#' share, because area is what the eye reads; the range is bounded so a rare
#' profile stays visible and a dominant one does not swamp the panel.
#'
#' @param share Numeric vector of shares, summing to one.
#' @param style Visual constants.
#' @return Numeric vector of `cex` values, one per share.
#' @noRd
.multilpa_point_sizes <- function(share, style) {
  stopifnot("`share` must be finite and non-negative" =
              is.numeric(share) && length(share) >= 1L &&
              all(is.finite(share)) && all(share >= 0))
  # Absolute, not min-max normalised within the plot. Normalising to the observed
  # range maps whatever spread happens to exist onto the full size range, so two
  # profiles at 52% and 48% are drawn as far apart as 80% and 20% -- the plot
  # manufactures a difference out of rounding. Here an evenly split profile is
  # always the default size, area is proportional to share, and the result is
  # bounded so a rare profile stays visible.
  reference <- share * length(share)
  pmin(2.0, pmax(0.55, sqrt(reference))) * style$point_size
}

#' Diverging fill colours around a white midpoint
#'
#' The house diverging pair, interpolated symmetrically so that zero is always
#' white and the two directions are distinguishable to a colour-blind reader.
#'
#' @param values Numeric vector to colour.
#' @param limit Value mapped to full saturation; defaults to the largest
#'   absolute value, so the scale is symmetric about zero.
#' @return Character vector of hexadecimal colours.
#' @noRd
.multilpa_diverging <- function(values, limit = max(abs(values), na.rm = TRUE)) {
  stopifnot("`values` must be numeric" = is.numeric(values),
            "`limit` must be a single finite number" =
              is.numeric(limit) && length(limit) == 1L && is.finite(limit))
  if (limit <= 0) return(rep("#FFFFFF", length(values)))
  intensity <- pmax(-1, pmin(1, values / limit))
  ends <- grDevices::col2rgb(c("#D33F6A", "#4A6FE3")) / 255
  white <- c(1, 1, 1)
  mixed <- vapply(intensity, function(t) {
    if (is.na(t)) return(c(1, 1, 1))
    end <- if (t < 0) ends[, 1L] else ends[, 2L]
    white + abs(t) * (end - white)
  }, numeric(3))
  grDevices::rgb(mixed[1L, ], mixed[2L, ], mixed[3L, ])
}

#' Draw one ridge of a binned-line ridge plot
#'
#' A base-graphics equivalent of `ggridges::geom_density_ridges(stat =
#' "binline")`: the binned counts are drawn as a step outline with a filled
#' area, sitting on a baseline at `base` and never rising above `base + height`.
#' Counts rather than a smoothed density, because these quantities are bounded
#' (a probability cannot exceed one) and a kernel would put mass outside the
#' support and invent a tail the data do not have.
#'
#' @param values Numeric vector for this ridge.
#' @param breaks Bin edges shared across every ridge, so heights are comparable.
#' @param base Vertical baseline for this ridge.
#' @param height Maximum height above the baseline.
#' @param fill,border Colours for the filled area and its outline.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_ridge <- function(values, breaks, base, height, fill, border) {
  stopifnot("`breaks` must be increasing" = length(breaks) >= 2L &&
              all(diff(breaks) > 0),
            "`height` must be a single positive number" =
              is.numeric(height) && length(height) == 1L && height > 0)
  counts <- graphics::hist(values, breaks = breaks, plot = FALSE)$counts
  if (max(counts) <= 0) return(invisible(NULL))
  scaled <- base + height * counts / max(counts)
  # Each bin contributes two x values so the outline steps rather than slopes,
  # which is what makes a bounded quantity's edges readable.
  xs <- rep(breaks, each = 2L)
  ys <- c(base, rep(scaled, each = 2L), base)
  graphics::polygon(xs, ys, col = grDevices::adjustcolor(fill, alpha.f = 0.55),
                    border = border, lwd = 1.2)
  invisible(NULL)
}

#' Draw a compact legend strip above the panel
#'
#' Grouped bars cannot carry direct end-labels the way lines can, so the key
#' sits above the panel. Each entry pairs its swatch with text, and the text
#' names the profile, so the legend itself never depends on colour.
#'
#' @param labels Character vector of entry labels.
#' @param colours Fill colour per entry.
#' @param style Visual constants.
#' @param lines Distance below the panel, in margin lines, counted the way
#'   `mtext(side = 1, line = )` counts.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_legend_strip <- function(labels, colours, style, lines = 4.4) {
  stopifnot("`labels` and `colours` must match" =
              length(labels) == length(colours))
  usr <- graphics::par("usr")
  widths <- graphics::strwidth(labels, cex = style$label_text_size)
  swatch <- graphics::strwidth("m", cex = style$label_text_size)
  gap <- swatch * 1.6
  entries <- widths + swatch * 2
  total <- sum(entries) + gap * (length(labels) - 1L)
  left <- usr[1L] + (diff(usr[1:2]) - total) / 2
  height <- graphics::strheight("M", cex = style$label_text_size)
  # `lines` counts downward from the panel edge, the way mtext() does, so the
  # caller can place the strip below the axis title without knowing the margin.
  # grconvertY() measures the other way -- upward from the outer margin edge --
  # so the two are reconciled here rather than at every call site.
  y <- graphics::grconvertY(graphics::par("mar")[1L] - lines,
                            from = "lines", to = "user")
  starts <- left + c(0, cumsum(entries + gap)[-length(entries)])
  graphics::rect(starts, y - height * 0.55, starts + swatch, y + height * 0.55,
                 col = colours, border = style$panel_fill, lwd = 1.2)
  graphics::text(starts + swatch * 1.5, y, labels, adj = c(0, 0.5),
                 cex = style$label_text_size, col = style$text_colour)
  invisible(NULL)
}
