# Model-selection evidence for an enumeration grid. Every table and plot here
# reports quantities a reader can check against the grid (its candidates
# table, or for the smallest class, each candidate's estimated sizes): how
# much a criterion improves with each added profile or group class (deltas),
# where maximum perpendicular chord distance places a knee (knees), the
# criterion over the whole grid (heatmap), the smallest class, and whether
# each estimate can be trusted. None of them selects a model.

utils::globalVariables(c("delta", "status", "tile_label", "smallest_share",
                         "n_group_label", "line", "key"))

#' The model-selection views of an enumeration grid
#' @return The view names `plot(grid, what = )` takes for model selection.
#' @noRd
.latents_selection_views <- function() {
  c("deltas", "heatmap", "smallest", "trust")
}

#' One row per candidate with the criterion the views read
#'
#' @param x A `multilpa_enumeration`.
#' @param using The criterion column.
#' @return A data frame: `candidate`, `n_profiles`, `n_group_classes`,
#'   `model`, `n_parameters`, `value`, and `usable` (converged with a finite
#'   value).
#' @noRd
.latents_selection_frame <- function(x, using) {
  criteria <- unname(.latents_grid_views(x, "The selection view")$criteria)
  if (!is.character(using) || length(using) != 1L || !using %in% criteria) {
    stop(errorCondition(sprintf(
      "`using` must name one information criterion of the grid: %s.",
      paste(sprintf("\"%s\"", criteria), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  # The criteria and diagnostics are stored; the AvePP columns, which need
  # every fit's posteriors, are not needed here.
  grid <- .latents_grid_table(x)
  value <- grid[[using]]
  data.frame(candidate = grid$candidate, n_profiles = grid$n_profiles,
             n_group_classes = grid$n_group_classes, model = grid$model,
             n_parameters = grid$n_parameters, value = value,
             usable = grid$converged %in% TRUE & is.finite(value),
             stringsAsFactors = FALSE)
}

#' The value of the neighbour one step smaller in one direction
#'
#' @param frame The output of `.latents_selection_frame()`.
#' @param direction `"profiles"` or `"group_classes"`.
#' @return For every row, the criterion of the usable candidate with one
#'   profile (or one group class) fewer and everything else equal, and that
#'   candidate's name; `NA` where there is none.
#' @noRd
.latents_smaller_neighbour <- function(frame, direction) {
  own <- if (identical(direction, "profiles")) "n_profiles" else "n_group_classes"
  other <- setdiff(c("n_profiles", "n_group_classes"), own)
  key <- function(own_value, rows) {
    paste(own_value, rows[[other]], ifelse(is.na(rows$model), "", rows$model))
  }
  usable <- frame[frame$usable, , drop = FALSE]
  at <- match(key(frame[[own]] - 1L, frame), key(usable[[own]], usable))
  list(value = usable$value[at], candidate = usable$candidate[at])
}

#' What each added profile and each added group class buys
#'
#' @param x A `multilpa_enumeration`.
#' @param using The criterion.
#' @return One row per candidate: `candidate`, `n_profiles`,
#'   `n_group_classes`, `model`, `criterion`, `value`, `delta_profile` (the
#'   criterion's decrease from the candidate with one profile fewer, the same
#'   group classes and structure; positive means this candidate is lower),
#'   `from_profile`, `delta_group_class` and `from_group_class` (the same for
#'   one group class fewer). `NA` where there is no such neighbour or this
#'   candidate did not converge.
#' @noRd
.latents_deltas_table <- function(x, using = "bic_groups") {
  frame <- .latents_selection_frame(x, using)
  profile <- .latents_smaller_neighbour(frame, "profiles")
  group <- .latents_smaller_neighbour(frame, "group_classes")
  own <- ifelse(frame$usable, frame$value, NA_real_)
  data.frame(candidate = frame$candidate, n_profiles = frame$n_profiles,
             n_group_classes = frame$n_group_classes, model = frame$model,
             criterion = using, value = own,
             delta_profile = profile$value - own,
             from_profile = ifelse(is.na(own), NA_character_, profile$candidate),
             delta_group_class = group$value - own,
             from_group_class = ifelse(is.na(own), NA_character_, group$candidate),
             stringsAsFactors = FALSE)
}

#' The maximum perpendicular chord-distance knee of one curve
#'
#' Both axes are rescaled to 0-1, and the knee is the point furthest below
#' the straight line from the first point to the last,
#' measured perpendicular to that line. This simplified Kneedle-inspired
#' statistic omits smoothing, local-maxima thresholds and sensitivity; a tiny
#' positive distance may therefore still name a knee.
#'
#' @param x,y The curve, ordered by `x`.
#' The curve is a criterion that should fall as the model grows; a line
#' whose last value is not below its first has no knee, because a bend in a
#' criterion that only rises marks no diminishing return.
#'
#' @return A list: `at` (the index of the knee, `NA` without one), `distance`
#'   (each point's perpendicular distance below the line, in the rescaled
#'   units) and `rising` (`TRUE` when the criterion does not fall along the
#'   line).
#' @noRd
.latents_kneedle <- function(x, y) {
  n <- length(x)
  none <- list(at = NA_integer_, distance = rep(NA_real_, n), rising = FALSE)
  if (n < 3L || diff(range(x)) <= 0 || diff(range(y)) <= 0) return(none)
  if (y[n] >= y[1L]) {
    none$rising <- TRUE
    return(none)
  }
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  slope <- y_scaled[n] - y_scaled[1L]
  below <- (y_scaled[1L] + slope * x_scaled) - y_scaled
  distance <- below / sqrt(1 + slope^2)
  at <- which.max(distance)
  # The end points lie on the line; a curve that never bends below it has no
  # knee.
  list(at = if (distance[at] > sqrt(.Machine$double.eps)) at else NA_integer_,
       distance = distance, rising = FALSE)
}

#' The knee of every line of the grid, in each direction
#'
#' @param x A `multilpa_enumeration`.
#' @param using The criterion.
#' @return One row per line: `direction` (`"profiles"`: along the number of
#'   profiles at fixed group classes; `"group_classes"`: along the group
#'   classes at fixed profiles; `"parameters"`: every candidate, by number of
#'   parameters), `line` (what is held fixed), `criterion`, `from`, `to`,
#'   `n_points`, `knee_candidate`, `distance`, `knee_without_largest`,
#'   `stable` and `note`.
#' @noRd
.latents_knees_table <- function(x, using = "bic_groups") {
  frame <- .latents_selection_frame(x, using)
  usable <- frame[frame$usable, , drop = FALSE]
  structure_text <- function(rows) {
    if (all(is.na(rows$model)) || length(unique(stats::na.omit(frame$model))) < 2L) "" else
      paste0(", ", rows$model[1L])
  }
  one_line <- function(rows, along, direction, line) {
    rows <- rows[order(rows[[along]]), , drop = FALSE]
    n <- nrow(rows)
    knee <- .latents_kneedle(rows[[along]], rows$value)
    shorter <- if (n >= 4L) .latents_kneedle(rows[[along]][-n], rows$value[-n])$at else
      NA_integer_
    found <- !is.na(knee$at)
    without <- if (is.na(shorter)) NA_character_ else rows$candidate[shorter]
    knee_name <- if (found) rows$candidate[knee$at] else NA_character_
    note <- if (n < 3L) "fewer than three candidates" else
      if (isTRUE(knee$rising)) "the criterion does not fall along this line" else
        if (!found) "the curve does not bend below the line" else
        if (n < 4L) "three candidates: stability not assessed" else
          if (is.na(without)) "no knee when the largest is left out" else
            if (identical(without, knee_name)) "" else
              "moves when the largest is left out"
    data.frame(direction = direction, line = line, criterion = using,
               from = rows$candidate[1L], to = rows$candidate[n], n_points = n,
               knee_candidate = knee_name,
               distance = if (found) knee$distance[knee$at] else NA_real_,
               knee_without_largest = without,
               stable = if (found && n >= 4L) identical(without, knee_name) else NA,
               note = note, stringsAsFactors = FALSE)
  }
  model_key <- ifelse(is.na(usable$model), "", usable$model)
  along_profiles <- lapply(split(usable, list(usable$n_group_classes, model_key),
                                 drop = TRUE), function(rows) {
    one_line(rows, "n_profiles", "profiles",
             sprintf("%d group class%s%s", rows$n_group_classes[1L],
                     if (rows$n_group_classes[1L] == 1L) "" else "es",
                     structure_text(rows)))
  })
  several_groups <- length(unique(usable$n_group_classes)) > 1L
  along_groups <- if (several_groups) {
    lapply(split(usable, list(usable$n_profiles, model_key), drop = TRUE),
           function(rows) {
             one_line(rows, "n_group_classes", "group_classes",
                      sprintf("%d profile%s%s", rows$n_profiles[1L],
                              if (rows$n_profiles[1L] == 1L) "" else "s",
                              structure_text(rows)))
           })
  } else list()
  by_parameters <- lapply(split(usable, model_key, drop = TRUE), function(rows) {
    one_line(rows, "n_parameters", "parameters",
             paste0("all candidates", structure_text(rows)))
  })
  result <- do.call(rbind, c(along_profiles, along_groups, by_parameters))
  if (is.null(result)) {
    return(data.frame(direction = character(), line = character(),
                      criterion = character(), from = character(),
                      to = character(), n_points = integer(),
                      knee_candidate = character(), distance = numeric(),
                      knee_without_largest = character(), stable = logical(),
                      note = character()))
  }
  result <- result[order(match(result$direction,
                               c("profiles", "group_classes", "parameters"))), ,
                   drop = FALSE]
  row.names(result) <- NULL
  result
}

#' Shared look of the selection views
#' @noRd
.latents_selection_look <- function(plot, main, subtitle, default_main,
                                    default_subtitle, x_label, legend = FALSE) {
  plot +
    .gg_titles(main, subtitle, default_main, default_subtitle) +
    ggplot2::labs(x = x_label) +
    .gg_theme() +
    ggplot2::theme(panel.spacing = ggplot2::unit(1.6, "lines"),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = if (legend) "bottom" else "none",
                   legend.title = ggplot2::element_blank())
}

#' "profiles" or "classes", as the grid's indicators make them
#' @noRd
.latents_class_word <- function(x) {
  fitted <- Filter(Negate(is.null), .subset2(x, "fits"))
  if (length(fitted) > 0L && length(.multilpa_continuous_names(fitted[[1L]])) == 0L) {
    "classes"
  } else "profiles"
}

#' Deltas view: what each added profile and group class buys
#'
#' Two panels. In each, a line holds the other count fixed, and every point
#' is the criterion's decrease from the neighbour one step smaller. Each line
#' is named at its end.
#' @noRd
.gg_selection_deltas <- function(x, using, main, subtitle) {
  deltas <- .latents_deltas_table(x, using)
  word <- .latents_grid_views(x)$word
  singular <- .latents_singular(word)
  several_models <- length(unique(deltas$model[!is.na(deltas$model)])) > 1L
  structure_suffix <- if (several_models) paste0(", ", deltas$model) else ""
  panel_profile <- sprintf("Adding a %s", singular)
  panel_group <- "Adding a group class"
  long <- rbind(
    data.frame(panel = panel_profile, candidate = deltas$candidate,
               x = deltas$n_profiles, delta = deltas$delta_profile,
               line = sprintf("%d group class%s%s", deltas$n_group_classes,
                              ifelse(deltas$n_group_classes == 1L, "", "es"),
                              structure_suffix),
               stringsAsFactors = FALSE),
    data.frame(panel = panel_group, candidate = deltas$candidate,
               x = deltas$n_group_classes, delta = deltas$delta_group_class,
               line = sprintf("%d %s%s", deltas$n_profiles,
                              ifelse(deltas$n_profiles == 1L, singular, word),
                              structure_suffix),
               stringsAsFactors = FALSE))
  long <- long[is.finite(long$delta), , drop = FALSE]
  if (nrow(long) == 0L) {
    stop(errorCondition(
      "No converged candidate has a converged neighbour one step smaller.",
      class = "latents_nothing_to_plot", call = NULL))
  }
  long$panel <- factor(long$panel, levels = c(panel_profile, panel_group))
  long$key <- paste(long$panel, long$line)
  keys <- unique(long$key)
  colours <- stats::setNames(rep_len(.gg_okabe_ito, length(keys)), keys)
  # Each line is labelled at its right end, so no legend has to be matched
  # by colour.
  ends <- do.call(rbind, lapply(split(long, long$key), \(s) s[which.max(s$x), ]))
  ends <- do.call(rbind, lapply(split(ends, ends$panel, drop = TRUE), function(p) {
    span <- diff(range(long$delta[long$panel == p$panel[1L]]))
    p$y <- .gg_spread(p$delta, gap = 0.06 * max(span, 1e-8))
    p
  }))
  plot <- ggplot2::ggplot(long, ggplot2::aes(x, delta, colour = key, group = key)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    # A line with a single point has no segment to draw.
    ggplot2::geom_line(data = long[stats::ave(long$x, long$key, FUN = length) > 1L, ,
                                   drop = FALSE], linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_text(data = ends, ggplot2::aes(x = x + 0.15, y = y, label = line),
                       hjust = 0, size = 3, fontface = "bold") +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_x_continuous(breaks = function(limits) {
      seq(ceiling(limits[1L]), floor(limits[2L]))
    }, expand = ggplot2::expansion(add = c(0.3, 1.6))) +
    ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(y = sprintf("Decrease in %s from one step smaller", using))
  .latents_selection_look(
    plot, main, subtitle, "What each added class buys",
    sprintf("Each point: how much %s fell from the candidate one step smaller",
            using),
    sprintf("Number of %s (left) or group classes (right) after the step", word))
}

#' Heatmap view: the criterion over profiles and group classes
#' @noRd
.gg_selection_heatmap <- function(x, using, main, subtitle) {
  frame <- .latents_selection_frame(x, using)
  if (!any(frame$usable)) {
    stop(errorCondition(sprintf("No converged candidate has a finite %s.", using),
                        class = "latents_nothing_to_plot", call = NULL))
  }
  lowest <- min(frame$value[frame$usable])
  frame$delta <- ifelse(frame$usable, frame$value - lowest, NA_real_)
  frame$tile_label <- ifelse(!frame$usable, "x", sprintf("+%.0f", frame$delta))
  frame$tile_label[which.min(frame$delta)] <- "lowest"
  frame$n_group_label <- factor(frame$n_group_classes,
                                levels = sort(unique(frame$n_group_classes)))
  several_models <- length(unique(frame$model[!is.na(frame$model)])) > 1L
  plot <- ggplot2::ggplot(frame, ggplot2::aes(n_profiles, n_group_label, fill = delta)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = tile_label), size = 3.2) +
    ggplot2::scale_fill_gradient(low = "#4A6FE3", high = "white", na.value = "grey85",
                                 name = sprintf("%s above the lowest", using)) +
    ggplot2::scale_x_continuous(breaks = sort(unique(frame$n_profiles))) +
    ggplot2::labs(y = "Number of group classes")
  if (several_models) plot <- plot + ggplot2::facet_wrap(ggplot2::vars(model))
  .latents_selection_look(
    plot, main, subtitle, sprintf("%s over the whole grid", using),
    "Each tile: the criterion above the lowest candidate; x = not fitted or not converged",
    sprintf("Number of %s", .latents_grid_views(x)$word)) +
    ggplot2::theme(legend.position = "right", panel.grid = ggplot2::element_blank())
}

#' Smallest-class view: the share of the smallest class of every candidate
#' @noRd
.gg_selection_smallest <- function(x, main, subtitle) {
  views <- .latents_grid_views(x)
  frame <- .latents_selection_frame(x, views$using)
  shares <- .latents_grid_smallest_shares(x)
  word <- views$word
  singular <- .latents_singular(word)
  long <- do.call(rbind, lapply(names(shares), function(level) {
    values <- shares[[level]]
    data.frame(frame[c("candidate", "n_profiles", "n_group_classes", "model")],
               panel = if (identical(level, "individuals"))
                 sprintf("Smallest %s", singular) else "Smallest group class",
               smallest_share = values,
               keep = frame$usable & is.finite(values), stringsAsFactors = FALSE)
  }))
  long <- long[long$keep, , drop = FALSE]
  if (nrow(long) == 0L) {
    stop(errorCondition("No converged candidate has two or more classes.",
                        class = "latents_nothing_to_plot", call = NULL))
  }
  long$panel <- factor(long$panel, levels = c(sprintf("Smallest %s", singular),
                                              "Smallest group class"))
  several_models <- length(unique(stats::na.omit(long$model))) > 1L
  long$line <- sprintf("%d group class%s%s", long$n_group_classes,
                       ifelse(long$n_group_classes == 1L, "", "es"),
                       if (several_models) paste0(", ", long$model) else "")
  lines <- unique(long$line[order(long$n_group_classes)])
  colours <- stats::setNames(rep_len(.gg_okabe_ito, length(lines)), lines)
  shapes <- stats::setNames(rep_len(.gg_shapes, length(lines)), lines)
  plot <- ggplot2::ggplot(long, ggplot2::aes(n_profiles, smallest_share, colour = line,
                                             group = line)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey35",
                        linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = line, fill = line), size = 2.3,
                        colour = "white", stroke = 0.6) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::scale_y_continuous(labels = function(v) sprintf("%.0f%%", 100 * v),
                                limits = c(0, NA)) +
    ggplot2::scale_x_continuous(breaks = sort(unique(long$n_profiles))) +
    ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free_y") +
    ggplot2::labs(y = "Estimated share of the smallest class")
  .latents_selection_look(
    plot, main, subtitle, "The smallest class of every candidate",
    "Dashed line: 5%, a common floor below which a class is hard to interpret",
    sprintf("Number of %s", word), legend = length(lines) > 1L)
}

#' Trust view: whether each candidate's estimate can be believed
#' @noRd
.gg_selection_trust <- function(x, main, subtitle) {
  grid <- .latents_grid_table(x)
  status <- ifelse(is.na(grid$log_likelihood), "not fitted",
                   ifelse(!(grid$converged %in% TRUE), "not converged",
                          ifelse(grid$boundary %in% TRUE, "on a boundary",
                                 ifelse(is.na(grid$n_best_replicated),
                                        "replication unavailable",
                                        ifelse(grid$n_best_replicated %in% 1L,
                                               "best found once", "best replicated")))))
  frame <- data.frame(candidate = grid$candidate, n_profiles = grid$n_profiles,
                      n_group_label = factor(grid$n_group_classes,
                                             levels = sort(unique(grid$n_group_classes))),
                      model = grid$model,
                      status = factor(status, levels = c("best replicated", "best found once",
                                                         "replication unavailable",
                                                         "on a boundary", "not converged",
                                                         "not fitted")),
                      tile_label = ifelse(is.na(grid$log_likelihood), "x",
                                          ifelse(is.na(grid$n_best_replicated), "",
                                                 sprintf("%dx", grid$n_best_replicated))),
                      stringsAsFactors = FALSE)
  several_models <- length(unique(frame$model[!is.na(frame$model)])) > 1L
  plot <- ggplot2::ggplot(frame, ggplot2::aes(n_profiles, n_group_label, fill = status)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = tile_label), size = 3.2) +
    ggplot2::scale_fill_manual(values = c("best replicated" = "#009E73",
                                          "best found once" = "#F0E442",
                                          "replication unavailable" = "#56B4E9",
                                          "on a boundary" = "#E69F00",
                                          "not converged" = "#D55E00",
                                          "not fitted" = "#999999"),
                               drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(frame$n_profiles))) +
    ggplot2::labs(y = "Number of group classes")
  if (several_models) plot <- plot + ggplot2::facet_wrap(ggplot2::vars(model))
  .latents_selection_look(
    plot, main, subtitle, "Convergence and repeated maxima",
    paste("Label: how many random starts reached the best log likelihood;",
          "blank = unavailable; x = not fitted. Replication does not prove a global maximum."),
    sprintf("Number of %s", .latents_grid_views(x)$word), legend = TRUE) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Draw one model-selection view
#' @param x A `multilpa_enumeration`.
#' @param view One of `.latents_selection_views()`.
#' @param using The criterion.
#' @param main,subtitle Titles.
#' @return A ggplot.
#' @noRd
.gg_selection_view <- function(x, view, using, main, subtitle) {
  switch(view,
         deltas = .gg_selection_deltas(x, using, main, subtitle),
         heatmap = .gg_selection_heatmap(x, using, main, subtitle),
         smallest = .gg_selection_smallest(x, main, subtitle),
         trust = .gg_selection_trust(x, main, subtitle))
}
