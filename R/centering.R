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
#' @param centering One of [.multilpa_centering_schemes()].
#' @return A list with the centred matrix `x`, and `offsets`, the matrix that
#'   was subtracted: one row per group under `"person"`, one row overall under
#'   `"grand"`, and `NULL` under `"none"`.
#' @noRd
.multilpa_center_indicators <- function(x, group_index, n_groups, centering) {
  stopifnot(
    "`x` must be a matrix" = is.matrix(x),
    "`centering` must be one scheme" =
      length(centering) == 1L && centering %in% .multilpa_centering_schemes())
  # An all-categorical fit has no Gaussian block to centre, and `as.matrix()`
  # of a zero-column frame is logical rather than numeric, so the width is
  # checked before the storage mode.
  if (ncol(x) == 0L || identical(centering, "none")) {
    return(list(x = x, offsets = NULL))
  }
  stopifnot("`x` must be numeric to be centred" = is.numeric(x))
  if (identical(centering, "grand")) {
    offsets <- matrix(colMeans(x, na.rm = TRUE), nrow = 1L, ncol = ncol(x),
                      dimnames = list(NULL, colnames(x)))
    return(list(x = sweep(x, 2L, offsets[1L, ], "-"), offsets = offsets))
  }
  # One mean per group and indicator. A group with no observed value for an
  # indicator has no mean to centre on, and subtracting a global stand-in would
  # silently mix the two scales, so that cell's offset is zero and the values
  # stay on the raw scale -- recorded, because `offsets` is returned.
  observed <- !is.na(x)
  totals <- rowsum(ifelse(observed, x, 0), group_index, reorder = FALSE)
  counts <- rowsum(observed * 1, group_index, reorder = FALSE)
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
  offsets <- x$centering_offsets
  if (is.null(offsets) || ncol(observed) == 0L) return(observed)
  aligned <- offsets[, colnames(observed), drop = FALSE]
  if (nrow(aligned) == 1L) return(sweep(observed, 2L, aligned[1L, ], "-"))
  observed - aligned[x$group_index, , drop = FALSE]
}

#' Refuse a centring that leaves an indicator with nothing to model
#'
#' Subtracting each unit's own mean removes exactly the between-unit variation.
#' An indicator observed once per unit has none of anything left afterwards,
#' and a constant column identifies no profile, so it is refused here rather
#' than surfacing as a degenerate variance later.
#'
#' @param x The centred indicator matrix.
#' @param centering The scheme that produced it.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_check_centered <- function(x, centering) {
  if (identical(centering, "none") || ncol(x) == 0L) return(invisible(NULL))
  flat <- vapply(seq_len(ncol(x)), function(column) {
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
