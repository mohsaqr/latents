# Cluster bootstrap inference, for every covariance structure.
#
# The Wald path differentiates one log variance per profile and indicator. A
# structure that constrains the volume, the shape or the orientation has fewer
# free parameters than that chart has coordinates, so its information matrix is
# singular by construction and the Wald path refuses it. The constrained model
# is still a model, and its natural parameters -- the means, the variances or
# covariances, the mixing probabilities -- are still estimable. The bootstrap
# reaches them without needing a chart at all: resample, refit inside the same
# family, and read the spread of the estimates.
#
# Groups are the resampling unit, not rows. Observations inside a group are not
# independent, and resampling rows would treat a repeated measure as fresh
# information and return intervals that are too narrow (Field & Welsh, 2007,
# Journal of the Royal Statistical Society B, 69, 369-390; Davison & Hinkley,
# 1997, Bootstrap Methods and their Application, chapter 3).

utils::globalVariables(c("estimate", "panel", "at", "series", "conf_low", "conf_high"))

#' The `multilpa()` arguments that refit a model as it was fitted
#'
#' A refit that reads `variance_model` and `covariance_model` alone silently
#' widens any structure that constrains the volume, the shape or the
#' orientation: a `VEI` fit comes back as `VVI`, a different model with two more
#' parameters. The three pieces of the structure carry all of it, and naming
#' them alongside `covariance_model` is an error, so a fit that records a
#' structure is refitted from the structure and a fit that does not is refitted
#' from the two switches it does record.
#'
#' @param x A fitted `multilpa` model.
#' @return A named list of arguments for [multilpa()].
#' @noRd
.multilpa_refit_arguments <- function(x) {
  stopifnot(inherits(x, "multilpa"))
  shared <- list(vars = x$vars, id = x$id, n_profiles = x$n_profiles,
                 n_group_classes = x$n_group_classes,
                 categorical = x$categorical %||% character(),
                 min_variance = x$min_variance,
                 min_probability = x$min_probability %||% 1e-10,
                 missing = x$missing %||% "error",
                 time = x$time,
                 centering = x$centering %||% "none",
                 scaling = x$scaling %||% "none",
                 drop_empty = x$drop_empty %||% TRUE,
                 prior = x$prior,
                 noise = isTRUE(x$noise),
                 acceleration = x$acceleration %||% "none",
                 weights = x$weights, ordinal = x$ordinal %||% character(),
                 count = x$count %||% character(),
                 count_model = x$extra_data$count_model %||% "poisson",
                 count_dispersion = x$extra_data$count_dispersion %||% "varying",
                 # The search the fit was made with; a fit from before these
                 # were recorded was made by the kmeans search and the
                 # default selection.
                 search = x$search %||% "kmeans",
                 select_start = x$select_start %||% "likelihood",
                 boundary = x$boundary_policy %||% "allow")
  c(shared, .multilpa_refit_measurement(x))
}

#' The arguments that refit a fit's covariance and its measurement extensions
#'
#' One place for every verb that refits a model (bootstrap inference, seed
#' sensitivity, the bootstrap likelihood ratio test), so a measurement option
#' cannot be carried by one refit and dropped by another. A fit that records
#' a covariance structure is refitted from the structure, one that does not
#' from its two switches; covariances for chosen pairs are refitted as the
#' same pattern (the pairs, read from the fit's pattern, and the variance
#' model, which a full structure argument would replace by a full
#' covariance); direct effects are refitted from their terms.
#'
#' @param x A fitted `multilpa` model.
#' @return A named list of [multilpa()] arguments.
#' @noRd
.multilpa_refit_measurement <- function(x) {
  extensions <- list(direct_effects = if (is.null(x$extra_data$direct)) character() else
    .latents_direct_terms(x$extra_data$direct$effects))
  if (!is.null(x$covariance_pattern)) {
    return(c(list(variance_model = x$variance_model,
                  local_dependence = .latents_pattern_pairs(x$covariance_pattern)),
             extensions))
  }
  structure <- x$covariance_structure
  family <- if (is.null(structure) || is.na(structure) ||
                !structure %in% .multilpa_structures()) {
    list(variance_model = x$variance_model,
         covariance_model = x$covariance_model %||% "diagonal")
  } else .multilpa_structure_arguments(structure)
  c(family, extensions)
}

#' The `fit_staged()` arguments that refit a staged model as it was fitted
#'
#' A staged fit is refitted by staging again, so that each resample estimates
#' its own measurement before its own group-class structure. That is what
#' carries the first stage's sampling variability into the second stage's
#' replicates; refitting with the measurement held would carry none of it.
#'
#' @param x A `multilpa` fit from [fit_staged()].
#' @return A named list of arguments for [fit_staged()], without `data` and the
#'   start and convergence controls.
#' @noRd
.multilpa_staged_refit_arguments <- function(x) {
  stopifnot(inherits(x, "multilpa"), isTRUE(x$staged))
  list(vars = x$vars, id = x$id, n_profiles = x$n_profiles,
       n_group_classes = x$n_group_classes,
       variance_model = x$variance_model,
       covariance_model = x$covariance_model %||% "diagonal",
       categorical = x$categorical %||% character(),
       min_variance = x$min_variance,
       min_probability = x$min_probability %||% 1e-10,
       missing = x$missing %||% "error", time = x$time,
       search = x$search %||% "kmeans")
}

#' Bootstrap replicates of the natural coefficients
#'
#' A staged fit is refitted through [fit_staged()], both stages on every
#' resample; any other fit through [multilpa()].
#'
#' @param x A fitted `multilpa` model.
#' @param data The fitting data.
#' @param iter How many resamples to draw.
#' @param n_starts,max_iter,tol Passed to each refit.
#' @return A list with `estimates`, an `iter` by parameters matrix whose failed
#'   rows are `NA`, and `messages`, the reason each failure gave.
#' @noRd
.multilpa_bootstrap_replicates <- function(x, data, iter, n_starts, max_iter, tol) {
  staged <- isTRUE(x$staged)
  estimator <- if (staged) fit_staged else multilpa
  arguments <- if (staged) .multilpa_staged_refit_arguments(x) else
    .multilpa_refit_arguments(x)
  controls <- .latents_refit_controls(x, n_starts = n_starts,
                                     max_iter = max_iter, tol = tol)
  if (staged) {
    # The supplied-start second stage can record no screening even when the
    # measurement searched an explicit pool. Replay that first stage's own
    # controls, including a measurement supplied separately to fit_staged().
    measurement_arguments <- .multilpa_refit_arguments(x$stage_one)
    # A supplied measurement can use other identifier and time column names.
    # One group class has unconditional measurement; replay it using the
    # resample's columns, recreating any single-level observation identifier.
    measurement_arguments["id"] <- list(
      if (isTRUE(x$stage_one$single_level)) NULL else x$id)
    measurement_arguments["time"] <- list(x$time)
    measurement_controls <- .latents_refit_controls(
      x$stage_one, n_starts = n_starts, max_iter = max_iter, tol = tol)
    measurement_arguments <- c(
      measurement_arguments[!names(measurement_arguments) %in%
                              names(measurement_controls)],
      measurement_controls)
  }
  .latents_cluster_bootstrap(
    data, x$group_index, x$n_groups, x$id, iter,
    refit = function(resampled) {
      measurement <- if (staged) list(measurement = do.call(
        multilpa, c(list(data = resampled), measurement_arguments))) else list()
      do.call(estimator, c(list(data = resampled), arguments, measurement,
                           list(n_starts = n_starts, max_iter = max_iter, tol = tol,
                                n_screen = controls$n_screen)))
    },
    extract = function(fit) {
      # A resample that misses a rare category refits that indicator with one
      # category fewer; its parameters cannot be set beside the original's,
      # so the replicate is dropped (and counted) rather than aligned.
      missed <- names(x$categorical_levels)[!vapply(names(x$categorical_levels), function(item) {
        identical(fit$categorical_levels[[item]], x$categorical_levels[[item]])
      }, logical(1))]
      if (length(missed) > 0L) {
        return(sprintf("the resample does not contain every category of %s",
                       paste(sprintf("\"%s\"", missed), collapse = ", ")))
      }
      unname(.multilpa_coefficients(.multilpa_align_labels(fit, x), "natural"))
    },
    reference_names = names(.multilpa_coefficients(x, "natural")))
}

#' Percentile bootstrap inference for a fitted model
#'
#' Called by [parameter_inference()] when `method = "bootstrap"`. The interval
#' is the percentile interval of the replicates and the standard error is their
#' standard deviation, so both are read off the same replicates and neither
#' assumes the estimate is normal around its truth.
#'
#' A [fit_staged()] result is bootstrapped by staging every resample again, so
#' the replicates vary the measurement as well as the group-class structure.
#' Its table therefore has a row for every parameter, the held measurement
#' included, and the second-stage errors carry the first stage's uncertainty.
#'
#' @param x A fitted `multilpa` model.
#' @param data The fitting data.
#' @param level Interval level.
#' @param iter How many resamples to draw.
#' @param n_starts,max_iter,tol Passed to each refit.
#' @param adjust Kept for a frame of the same shape as the Wald path's.
#' @return A `data.frame`, one row per natural coefficient.
#' @noRd
.multilpa_bootstrap_inference <- function(x, data, level, iter, n_starts,
                                          max_iter, tol, adjust) {
  staged <- isTRUE(x$staged)
  if (length(x$fixed %||% character()) > 0L && !staged) {
    stop(errorCondition(paste(
      "A fit that holds a block fixed cannot be bootstrapped here: the held",
      "values came from another fit, and resampling these data does not",
      "resample them. Bootstrap the fit the measurement came from, or fit",
      "both stages with fit_staged(), whose bootstrap refits each stage."),
      class = "latents_unsupported_inference", call = NULL))
  }
  if (!isTRUE(x$converged)) {
    stop(errorCondition(
      "The original fit did not converge; bootstrap it only once it has.",
      class = "latents_no_converge", call = NULL))
  }
  ## Every coefficient of a staged fit is resampled, the held measurement
  ## included, so every one of them gets a row.
  free_natural <- if (staged) seq_along(.multilpa_coefficients(x, "natural")) else
    .multilpa_free_index(x, "natural")
  drawn <- .multilpa_bootstrap_replicates(x, data, iter, n_starts, max_iter, tol)
  estimates <- drawn$estimates[, free_natural, drop = FALSE]
  valid <- stats::complete.cases(estimates)
  if (sum(valid) < 2L) {
    reasons <- drawn$messages[!is.na(drawn$messages)]
    stop(errorCondition(sprintf(
      "Only %d of %d resamples produced a usable fit; the first reason was: %s",
      sum(valid), iter,
      if (length(reasons) > 0L) reasons[1L] else "no reason was recorded"),
      class = "latents_bootstrap_failed", call = NULL))
  }
  kept <- estimates[valid, , drop = FALSE]
  if (sum(valid) < iter) {
    warning(warningCondition(sprintf(
      "%d of %d resamples did not produce a usable fit and were dropped.",
      iter - sum(valid), iter), class = "latents_bootstrap_dropped", call = NULL))
  }
  bounds <- c((1 - level) / 2, (1 + level) / 2)
  quantiles <- apply(kept, 2L, stats::quantile, probs = bounds, names = FALSE)
  point <- .multilpa_coefficients(x, "natural")[free_natural]
  result <- data.frame(
    .multilpa_coefficient_labels(x, "natural")[free_natural, , drop = FALSE],
    estimate = unname(point),
    standard_error = apply(kept, 2L, stats::sd),
    ## A Wald statistic would put a normal approximation back on top of the
    ## replicates the interval was read from, which is the assumption this path
    ## exists to avoid. The interval is the inference.
    statistic = NA_real_, p_value = NA_real_,
    conf_low = unname(quantiles[1L, ]), conf_high = unname(quantiles[2L, ]),
    row.names = NULL, stringsAsFactors = FALSE)
  result <- .multilpa_adjust_p(result, adjust)
  covariance <- stats::cov(kept)
  dimnames(covariance) <- list(names(point), names(point))
  attributes(result) <- c(attributes(result), list(
    covariance = covariance, covariance_unconstrained = NULL,
    level = level, method = "bootstrap", iter = iter, n_valid = sum(valid),
    replicates = kept, messages = drawn$messages,
    structure = x$covariance_structure,
    fixed = x$fixed %||% character(),
    stages_resampled = if (staged) 2L else 1L))
  .latents_bootstrap_table(result, valid)
}

#' Mark a bootstrap table, so that plot() and as.data.frame() reach its replicates
#'
#' Shared by the profile, transition and growth bootstraps. The table stays a
#' data frame (one row per parameter, its label columns before `estimate`);
#' the `replicates` attribute holds the kept replicates, one column per row,
#' `replicate_number` the resample each kept replicate came from, and
#' `replicate_rows` the row names the replicate columns belong to: `[` keeps
#' attributes when it only reorders or repeats rows, so the row names are what
#' tells a reordered copy from the table itself.
#' @param result The bootstrap table, with its `replicates` attribute.
#' @param valid One logical per resample: whether it produced a usable fit.
#' @return `result` with class `c("latents_bootstrap", "data.frame")`.
#' @noRd
.latents_bootstrap_table <- function(result, valid) {
  attr(result, "replicate_number") <- which(valid)
  attr(result, "replicate_rows") <- row.names(result)
  class(result) <- c("latents_bootstrap", "data.frame")
  result
}

#' The label columns of a bootstrap table and the column naming the kind
#' @return A list: `labels`, the columns before `estimate`, with a second-order
#'   transition table's `previous` (which follows the statistics) before
#'   `from`; `kind`, the column that names the kind of parameter
#'   (`parameter`, or `block` for a transition or growth model), or `NULL`.
#' @noRd
.latents_bootstrap_columns <- function(x) {
  labels <- names(x)[seq_len(match("estimate", names(x)) - 1L)]
  if ("previous" %in% names(x)) {
    labels <- append(labels, "previous", after = match("from", labels, nomatch = 0L) - 1L)
  }
  kind <- intersect(c("parameter", "block"), labels)
  list(labels = labels, kind = if (length(kind) > 0L) kind[[1L]] else NULL)
}

#' The kept replicates of a bootstrap table, one column per row of the table
#' @return A replicates-by-parameters matrix, or `latents_bad_argument` when a
#'   subset, reordered or modified copy of the table no longer matches them.
#' @noRd
.latents_bootstrap_matrix <- function(x) {
  kept <- attr(x, "replicates")
  if (is.null(kept) || !identical(ncol(kept), nrow(x)) ||
      !identical(attr(x, "replicate_rows"), row.names(x))) {
    stop(errorCondition(paste(
      "This table no longer matches its bootstrap replicates (a subset or a",
      "reordered or modified copy loses them). Use the table",
      "parameter_inference(method = \"bootstrap\") returned."),
      class = "latents_bad_argument", call = NULL))
  }
  kept
}

#' The replicates of a bootstrap table, one row per replicate and parameter
#' @return A base data frame: `replicate` (the resample's number), the
#'   table's label columns, and `estimate`, the parameter in that replicate.
#' @noRd
.latents_bootstrap_replicates <- function(x) {
  kept <- .latents_bootstrap_matrix(x)
  labels <- .latents_bootstrap_columns(x)$labels
  table <- as.data.frame(x)
  number <- attr(x, "replicate_number") %||% seq_len(nrow(kept))
  index <- expand.grid(parameter = seq_len(ncol(kept)), replicate = seq_len(nrow(kept)))
  data.frame(replicate = number[index$replicate],
             table[index$parameter, labels, drop = FALSE],
             estimate = kept[cbind(index$replicate, index$parameter)],
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @noRd
#' @export
as.data.frame.latents_bootstrap <- function(x, row.names = NULL, optional = FALSE,
                                            what = c("estimates", "replicates"), ...) {
  what <- .latents_match_formal(what)
  if (identical(what, "replicates")) return(.latents_bootstrap_replicates(x))
  table <- x
  attributes(table) <- list(names = names(x), row.names = attr(x, "row.names"),
                            class = "data.frame")
  if (!is.null(row.names)) row.names(table) <- row.names
  table
}

#' @noRd
#' @export
plot.latents_bootstrap <- function(x, what = NULL, main = NULL, subtitle = NULL, ...) {
  .gg_require()
  kept <- .latents_bootstrap_matrix(x)
  columns <- .latents_bootstrap_columns(x)
  table <- as.data.frame(x)
  kinds <- if (is.null(columns$kind)) character() else unique(table[[columns$kind]])
  if (!is.null(what) && !(is.character(what) && length(what) > 0L && all(what %in% kinds))) {
    stop(errorCondition(sprintf(
      "`what` must name kinds of parameter in this table: %s.",
      if (length(kinds) == 0L) "this table has none" else
        paste(sprintf("\"%s\"", kinds), collapse = ", ")),
      class = "latents_bad_argument", call = NULL))
  }
  chosen <- if (is.null(what)) rep(TRUE, nrow(table)) else table[[columns$kind]] %in% what
  # A parameter held fixed or determined by the others (the probability of a
  # single group class) has no interval to draw; the caption names it instead.
  spread <- apply(kept, 2L, function(values) diff(range(values)) > 1e-10 * max(1, abs(values)))
  shown <- which(chosen & spread)
  if (length(shown) == 0L) {
    stop(errorCondition(
      "No chosen parameter varies across the replicates, so there is nothing to plot.",
      class = c("latents_nothing_to_plot", "latents_bad_argument"), call = NULL))
  }
  frame <- .latents_bootstrap_layout(table, columns, shown)
  level <- attr(x, "level") %||% 0.95
  kind_text <- if (is.null(columns$kind)) "parameters" else
    paste(unique(table[[columns$kind]][shown]), collapse = ", ")
  main <- main %||% sprintf("Bootstrap intervals: %s", kind_text)
  subtitle <- subtitle %||% sprintf(paste(
    "Estimates with %s%% percentile intervals from %d refitted resamples (of %d)."),
    format(100 * level), nrow(kept), attr(x, "iter") %||% nrow(kept))
  flat <- which(chosen & !spread)
  flat_names <- vapply(flat, function(row) {
    parts <- vapply(columns$labels, function(name) as.character(table[[name]][row]),
                    character(1))
    paste(parts[!is.na(parts)], collapse = " · ")
  }, character(1))
  # Series are told apart by colour and shape; a panel with a single series
  # (the group classes' probabilities) is drawn in black without a key.
  dodge <- ggplot2::position_dodge(width = 0.35)
  several <- frame[!is.na(frame$series), , drop = FALSE]
  single <- frame[is.na(frame$series), , drop = FALSE]
  joined <- several[several$series_points > 1L, , drop = FALSE]
  series <- unique(several$series)
  palette <- rep_len(.gg_okabe_ito[-7L], length(series))
  shapes <- rep_len(c(16L, 17L, 15L, 18L, 1L, 2L, 0L, 5L), length(series))
  plot <- ggplot2::ggplot(frame, ggplot2::aes(at, estimate)) +
    ggplot2::geom_line(data = joined, ggplot2::aes(colour = series, group = series),
                       position = dodge, linewidth = 0.6) +
    ggplot2::geom_errorbar(data = several, ggplot2::aes(ymin = conf_low, ymax = conf_high,
                                                        colour = series, group = series),
                           position = dodge, width = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(data = several, ggplot2::aes(colour = series, shape = series,
                                                     group = series),
                        position = dodge, size = 2.4) +
    ggplot2::geom_errorbar(data = single, ggplot2::aes(ymin = conf_low, ymax = conf_high),
                           colour = "#000000", width = 0.15, linewidth = 0.6) +
    ggplot2::geom_point(data = single, colour = "#000000", size = 2.4) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_colour_manual(values = stats::setNames(palette, series), breaks = series,
                                 name = NULL) +
    ggplot2::scale_shape_manual(values = stats::setNames(shapes, series), breaks = series,
                                name = NULL) +
    ggplot2::labs(title = main, subtitle = subtitle, x = NULL, y = "Estimate",
                  caption = if (length(flat) > 0L) sprintf(
                    "Not drawn (no spread across replicates): %s.",
                    paste(flat_names, collapse = "; "))) +
    .gg_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "bottom",
                   panel.spacing = ggplot2::unit(1, "lines"),
                   panel.grid.minor = ggplot2::element_blank())
  plot
}

#' Where each drawn parameter of a bootstrap table sits on the plot
#'
#' One panel per kind of parameter, split by `level` when a kind spans
#' several (the probabilities of profiles and of group classes), which is
#' then named with its level. Within a
#' panel, the last label column that varies is the horizontal axis (the
#' indicator of a mean, the group class of a profile's probability, the
#' occasion of a transition) and the others that vary name the series that
#' a line joins (the profile, the transition); a panel where only one column
#' varies has a single series.
#' @param table The bootstrap table as a base data frame.
#' @param columns From .latents_bootstrap_columns().
#' @param shown The rows drawn.
#' @return A base data frame, one row per drawn parameter: `panel` and `at`
#'   (factors in table order), `series` (`NA` for a panel's single series),
#'   `series_points` (how many points the series has in its panel),
#'   `estimate`, `conf_low` and `conf_high`.
#' @noRd
.latents_bootstrap_layout <- function(table, columns, shown) {
  rows <- table[shown, , drop = FALSE]
  text <- function(name) {
    values <- as.character(rows[[name]])
    values[is.na(values)] <- "(none)"
    values
  }
  kind <- if (is.null(columns$kind)) rep("estimate", nrow(rows)) else text(columns$kind)
  levels_of_kind <- if ("level" %in% columns$labels) {
    tapply(text("level"), kind, function(values) length(unique(values)))[kind]
  } else rep(1L, nrow(rows))
  panel <- ifelse(levels_of_kind > 1L, paste(text("level"), kind), kind)
  others <- setdiff(columns$labels, c(columns$kind, "level"))
  pieces <- lapply(split(seq_len(nrow(rows)), factor(panel, levels = unique(panel))),
                   function(index) {
    varying <- Filter(function(name) length(unique(text(name)[index])) > 1L, others)
    at <- if (length(varying) == 0L) panel[index] else text(varying[length(varying)])[index]
    series <- if (length(varying) < 2L) rep(NA_character_, length(index)) else {
      do.call(paste, c(lapply(varying[-length(varying)], function(name) text(name)[index]),
                       sep = " · "))
    }
    data.frame(row = index, at = at, series = series, stringsAsFactors = FALSE)
  })
  placed <- do.call(rbind, pieces)
  placed <- placed[order(placed$row), , drop = FALSE]
  key <- paste(panel, placed$series)
  data.frame(panel = factor(panel, levels = unique(panel)),
             at = factor(placed$at, levels = unique(placed$at)),
             series = placed$series,
             series_points = as.vector(table(key)[key]),
             estimate = rows$estimate, conf_low = rows$conf_low, conf_high = rows$conf_high,
             row.names = NULL, stringsAsFactors = FALSE)
}
