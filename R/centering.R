# Centering the continuous indicators before the measurement model sees them.
#
# A within-person analysis asks which enrolments were unusual *for that
# student*, not which were unusual overall. The transform that asks it is
# subtracting each unit's own mean from its observations, and until now a
# caller had to do that to the data frame before calling the fitter -- which
# left the fit holding numbers that no longer matched the frame it was given,
# so every verb that checks row alignment was comparing centred values against
# raw ones.
#
# Centring here instead keeps the offsets on the fit, so the raw data can be
# reconstructed exactly and the alignment checks still work.

#' The centring schemes `multilpa()` accepts
#' @return A character vector, the first being the default.
#' @noRd
.multilpa_centering_schemes <- function() c("none", "person", "grand")

#' The scaling schemes `multilpa()` accepts
#'
#' The same levels as centring: `"person"` divides each group's values by that
#' group's own standard deviation, `"grand"` by one standard deviation per
#' indicator. With the matching centring they give the classical within-person
#' and grand-mean z-scores.
#' @return A character vector, the first being the default.
#' @noRd
.multilpa_scaling_schemes <- function() c("none", "person", "grand")

#' Centre the continuous indicators
#'
#' `"person"` subtracts each group's own mean from its rows, so a value reads
#' as a deviation from that unit's average and the profiles become profiles of
#' *change* (Quintana, 2021; Voelkle et al., 2014). `"grand"` subtracts one
#' mean per indicator, which shifts the origin without touching the
#' within-group structure. `"none"` returns the matrix unchanged.
#'
#' Missing values are ignored when the offsets are formed and stay missing
#' afterwards, so an incomplete row is centred on what its unit was observed
#' doing rather than dropped.
#'
#' @param x The numeric indicator matrix, possibly containing `NA`.
#' @param group_index One group index per row.
#' @param n_groups How many groups there are.
#' @param centering One of `.multilpa_centering_schemes()`.
#' @param scaling One of `.multilpa_scaling_schemes()`.
#' @return A list with the transformed matrix `x`, `offsets`, the matrix that
#'   was subtracted, and `scales`, the matrix that the result was divided by:
#'   each one row per group under `"person"`, one row overall under `"grand"`,
#'   and `NULL` under `"none"`.
#' @noRd
.multilpa_center_indicators <- function(x, group_index, n_groups, centering,
                                        scaling = "none") {
  stopifnot(
    "`x` must be a matrix" = is.matrix(x),
    "`centering` must be one scheme" =
      length(centering) == 1L && centering %in% .multilpa_centering_schemes(),
    "`scaling` must be one scheme" =
      length(scaling) == 1L && scaling %in% .multilpa_scaling_schemes())
  # An all-categorical fit has no Gaussian block to centre, and `as.matrix()`
  # of a zero-column frame is logical rather than numeric, so the width is
  # checked before the storage mode.
  if (ncol(x) == 0L ||
      (identical(centering, "none") && identical(scaling, "none"))) {
    return(list(x = x, offsets = NULL, scales = NULL))
  }
  stopifnot("`x` must be numeric to be centred" = is.numeric(x))
  # A standard deviation does not move with a shift, so the scales are read
  # off the raw values whatever the centring.
  scales <- .multilpa_indicator_scales(x, group_index, n_groups, scaling)
  centred <- .multilpa_center_only(x, group_index, n_groups, centering)
  scaled <- if (is.null(scales)) centred$x else if (nrow(scales) == 1L)
    sweep(centred$x, 2L, scales[1L, ], "/") else
      centred$x / scales[group_index, , drop = FALSE]
  list(x = scaled, offsets = centred$offsets, scales = scales)
}

#' The standard deviation each indicator is divided by
#'
#' The sample standard deviation (denominator `n - 1`) of the observed values:
#' one per indicator under `"grand"`, one per group and indicator under
#' `"person"`. A scale that is missing or zero -- a group observed once on an
#' indicator, or always giving the same value -- cannot be divided by, and a
#' stand-in would put that group on a different scale from the rest, so it is
#' refused with `latents_bad_data` naming the indicators and how many groups.
#'
#' @param x The numeric indicator matrix, possibly containing `NA`.
#' @param group_index,n_groups The group of each row, and how many there are.
#' @param scaling One of `.multilpa_scaling_schemes()`.
#' @return A matrix of standard deviations, or `NULL` under `"none"`.
#' @noRd
.multilpa_indicator_scales <- function(x, group_index, n_groups, scaling) {
  if (identical(scaling, "none")) return(NULL)
  observed <- !is.na(x)
  units <- if (identical(scaling, "person")) group_index else rep(1L, nrow(x))
  n_units <- if (identical(scaling, "person")) n_groups else 1L
  # Sums of values and of squares by unit give every standard deviation in
  # one pass; `reorder = TRUE` sorts them by unit index.
  filled <- ifelse(observed, x, 0)
  counts <- rowsum(observed * 1, units, reorder = TRUE)
  totals <- rowsum(filled, units, reorder = TRUE)
  present <- sort(unique(units))
  means <- totals / pmax(counts, 1)
  deviations <- (filled - means[match(units, present), , drop = FALSE]) * observed
  squares <- rowsum(deviations^2, units, reorder = TRUE)
  standard <- matrix(NA_real_, n_units, ncol(x),
                     dimnames = list(NULL, colnames(x)))
  # A unit with no observed value of an indicator has nothing to divide, so
  # its scale is one and touches nothing; only a unit with values whose
  # standard deviation is missing or zero is a problem.
  standard[present, ] <- ifelse(counts == 0, 1,
                                ifelse(counts > 1,
                                       sqrt(squares / pmax(counts - 1, 1)),
                                       NA_real_))
  # A column that does not vary is judged relative to its own magnitude, so
  # the verdict does not depend on the units it was recorded in.
  magnitude <- matrix(apply(abs(x), 2L, max, na.rm = TRUE), n_units, ncol(x),
                      byrow = TRUE)
  degenerate <- is.na(standard) |
    standard <= sqrt(.Machine$double.eps) * pmax(magnitude, 1e-300)
  if (any(degenerate)) {
    columns <- colnames(x)[colSums(degenerate) > 0L]
    stop(errorCondition(sprintf(paste(
      "`scaling = \"%s\"` divides by a standard deviation that is zero or",
      "undefined for %s%s. %s"),
      scaling, paste(sprintf("`%s`", columns), collapse = ", "),
      if (identical(scaling, "person"))
        sprintf(" in %d group(s)", sum(rowSums(degenerate) > 0L)) else "",
      if (identical(scaling, "person"))
        "Every group needs two or more different values of each indicator; drop the groups that do not have them, or scale with `scaling = \"grand\"`." else
        "A constant indicator identifies no profile; drop it."),
      class = "latents_bad_data", call = NULL))
  }
  standard
}

#' Subtract the centring offsets, without scaling
#' @noRd
.multilpa_center_only <- function(x, group_index, n_groups, centering) {
  if (identical(centering, "none")) return(list(x = x, offsets = NULL))
  if (identical(centering, "grand")) {
    offsets <- matrix(colMeans(x, na.rm = TRUE), nrow = 1L, ncol = ncol(x),
                      dimnames = list(NULL, colnames(x)))
    return(list(x = sweep(x, 2L, offsets[1L, ], "-"), offsets = offsets))
  }
  # One mean per group and indicator. A group with no observed value for an
  # indicator has no mean to centre on, and subtracting a global stand-in would
  # silently mix the two scales, so that cell's offset is zero and the values
  # stay on the raw scale -- recorded, because `offsets` is returned.
  # `reorder = TRUE` sorts the sums by group index, which is the order of
  # `present`; first-occurrence order would pair them only when the indices
  # happen to appear in increasing order.
  observed <- !is.na(x)
  totals <- rowsum(ifelse(observed, x, 0), group_index, reorder = TRUE)
  counts <- rowsum(observed * 1, group_index, reorder = TRUE)
  present <- sort(unique(group_index))
  offsets <- matrix(0, nrow = n_groups, ncol = ncol(x),
                    dimnames = list(NULL, colnames(x)))
  offsets[present, ] <- ifelse(counts > 0, totals / pmax(counts, 1), 0)
  list(x = x - offsets[group_index, , drop = FALSE], offsets = offsets)
}

#' Undo the centring a fit applied
#'
#' The fit stores the indicators as the measurement model saw them. Every verb
#' that checks a supplied frame against the fit needs the values the caller
#' supplied instead, which is the centred matrix with its offsets added back.
#'
#' @param x A fitted model of this package.
#' @return The indicator matrix on the scale of the data that was handed in, or
#'   `NULL` when the fit kept no indicators.
#' @noRd
.multilpa_uncentered_indicators <- function(x) {
  stored <- x$indicator_data
  if (is.null(stored)) return(NULL)
  # Undo in the reverse order of the transform: multiply back, then add back.
  scales <- x$scaling_factors
  if (!is.null(scales)) {
    stored <- if (nrow(scales) == 1L) sweep(stored, 2L, scales[1L, ], "*") else
      stored * scales[x$group_index, , drop = FALSE]
  }
  offsets <- x$centering_offsets
  if (is.null(offsets)) return(stored)
  if (nrow(offsets) == 1L) return(sweep(stored, 2L, offsets[1L, ], "+"))
  stored + offsets[x$group_index, , drop = FALSE]
}

#' Put a supplied frame on the scale the model was fitted on
#'
#' A diagnostic that compares observed associations with model-implied ones
#' has to work on the scale the model was estimated on. The caller supplies raw
#' columns, because that is what they have; the fit knows what it subtracted.
#'
#' @param x A fitted model of this package.
#' @param observed A numeric matrix of indicator columns, on the input scale.
#' @return The same matrix, centred as the fit centred its own.
#' @noRd
.multilpa_center_like <- function(x, observed) {
  if (ncol(observed) == 0L) return(observed)
  offsets <- x$centering_offsets
  if (!is.null(offsets)) {
    aligned <- offsets[, colnames(observed), drop = FALSE]
    observed <- if (nrow(aligned) == 1L) sweep(observed, 2L, aligned[1L, ], "-") else
      observed - aligned[x$group_index, , drop = FALSE]
  }
  scales <- x$scaling_factors
  if (is.null(scales)) return(observed)
  aligned <- scales[, colnames(observed), drop = FALSE]
  if (nrow(aligned) == 1L) return(sweep(observed, 2L, aligned[1L, ], "/"))
  observed / aligned[x$group_index, , drop = FALSE]
}

#' Was a fit made on transformed indicators?
#'
#' A centred or scaled fit reports its means on the transformed scale, so the
#' views and stages that put estimates beside raw values refuse it alike.
#' @param x A fitted model.
#' @return A single logical.
#' @noRd
.multilpa_is_transformed <- function(x) {
  !identical(x$centering %||% "none", "none") ||
    !identical(x$scaling %||% "none", "none")
}

#' Refuse a centring that leaves an indicator with nothing to model
#'
#' Subtracting each unit's own mean removes exactly the between-unit variation.
#' An indicator observed once per unit has none of anything left afterwards,
#' and a constant column identifies no profile, so it is refused here rather
#' than surfacing as a degenerate variance later.
#'
#' A person-centred column is constant exactly when every unit's observed
#' values are all equal (its deviations sum to zero, so a constant deviation
#' is zero), which is read off the raw values without a tolerance. Judging the
#' centred values against an absolute threshold instead made the verdict
#' depend on the units: indicators recorded in units of `1e-9` were refused
#' under either scheme while `centering = "none"` fitted them.
#'
#' @param x The centred indicator matrix.
#' @param centering The scheme that produced it.
#' @param raw,group_index The indicator matrix before centring and the unit
#'   of each row. Without them the centred values are compared with an
#'   absolute `sqrt(.Machine$double.eps)`, which suits values of order one
#'   only.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_check_centered <- function(x, centering, raw = NULL,
                                     group_index = NULL) {
  if (identical(centering, "none") || ncol(x) == 0L) return(invisible(NULL))
  exact <- !is.null(raw) && !is.null(group_index)
  stopifnot("`raw` must match `x`" = !exact || identical(dim(raw), dim(x)),
            "`group_index` must have one unit per row" =
              !exact || length(group_index) == nrow(x))
  flat <- vapply(seq_len(ncol(x)), function(column) {
    if (exact) {
      observed <- !is.na(raw[, column])
      units <- if (identical(centering, "person")) group_index[observed] else
        rep(1L, sum(observed))
      within <- tapply(raw[observed, column], units,
                       function(values) length(unique(values)) < 2L)
      return(sum(observed) == 0L || all(within))
    }
    values <- x[, column]
    values <- values[!is.na(values)]
    length(values) == 0L || diff(range(values)) < sqrt(.Machine$double.eps)
  }, logical(1))
  if (!any(flat)) return(invisible(NULL))
  stop(errorCondition(sprintf(paste(
    "Centring left %s constant, so it can identify no profile. With",
    "`centering = \"%s\"` every unit needs more than one observation of an",
    "indicator for anything of it to survive."),
    paste(sprintf("`%s`", colnames(x)[flat]), collapse = ", "), centering),
    class = "latents_bad_data", call = NULL))
}
