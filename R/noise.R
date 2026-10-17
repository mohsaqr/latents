# A uniform noise component, as in mclust (`Mclust(..., initialization =
# list(noise = ))`, `me(..., Vinv = )`). Banfield, J. D., & Raftery, A. E.
# (1993). Model-based Gaussian and non-Gaussian clustering. Biometrics, 49,
# 803-821.
#
# The mixture gains one component with constant density `1 / V` over the
# hypervolume `V` of the data, so observations no Gaussian profile explains
# well --- outliers, scatter --- are absorbed by it rather than distorting a
# profile. It has no measurement parameters; its only parameter is its mixing
# proportion. During estimation it is carried as one extra, last, column of
# the profile probabilities and the posteriors; the fitted object splits it
# out again (`noise_probability`, `noise_posteriors`), so every field that
# describes a Gaussian profile keeps its shape.

#' The logarithm of mclust's default data hypervolume
#'
#' `mclust::hypvol()`: the smaller of the volume of the axis-aligned bounding
#' box of the data and the volume of the bounding box aligned with its
#' principal components. In logs, since a product of ranges overflows.
#'
#' @param x A complete numeric matrix.
#' @return A single number, `log(V)`.
#' @noRd
.multilpa_log_hypervolume <- function(x) {
  stopifnot("`x` must be a complete numeric matrix" =
              is.matrix(x) && is.numeric(x) && !anyNA(x))
  log_ranges <- function(values) {
    sum(log(apply(values, 2L, function(column) diff(range(column)))))
  }
  if (ncol(x) == 1L) return(log_ranges(x))
  centred <- sweep(x, 2L, colMeans(x), "-")
  # princomp's scores: the centred data on the eigenvectors of its covariance.
  # The divisor does not change the eigenvectors, and a range does not care
  # about a sign.
  rotation <- eigen(crossprod(centred) / nrow(x), symmetric = TRUE)$vectors
  min(log_ranges(x), log_ranges(centred %*% rotation))
}

#' Refuse a noise component this package does not define
#'
#' @param noise The [multilpa()] argument.
#' @param n_group_classes,categorical,fixed The corresponding arguments.
#' @param n_covariates How many membership covariates were named.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_check_noise <- function(noise, n_group_classes, categorical, fixed,
                                  n_covariates) {
  if (!(isTRUE(noise) || isFALSE(noise))) {
    stop(errorCondition("`noise` must be TRUE or FALSE.",
                        class = "latents_bad_argument", call = NULL))
  }
  if (isFALSE(noise)) return(invisible(NULL))
  refuse <- function(message) {
    stop(errorCondition(message, class = "latents_unsupported_noise", call = NULL))
  }
  if (!isTRUE(as.integer(n_group_classes) == 1L)) {
    refuse(paste(
      "A noise component is available for one group class (or `id = NULL`)",
      "only: with several, its proportion would have to differ by group class,",
      "which is not implemented. Set `n_group_classes = 1`."))
  }
  if (length(categorical) > 0L) {
    refuse("A noise component is a uniform density over continuous indicators; it cannot be combined with `categorical` ones.")
  }
  if (length(fixed) > 0L) {
    refuse("A noise component cannot be combined with `fixed`.")
  }
  if (n_covariates > 0L) {
    refuse("A noise component cannot be combined with `profile_covariates` or `group_covariates`.")
  }
  invisible(NULL)
}

#' Take the noise proportion out of a supplied start
#'
#' A start for a noise fit carries `noise_probability` beside the profile
#' probabilities, which together sum to one. The profile probabilities are
#' returned renormalized over the Gaussian profiles, which is the shape every
#' other start has, so the ordinary validation applies to them unchanged.
#'
#' @param start The caller's start list.
#' @return A list with `start` and `noise_probability`.
#' @noRd
.multilpa_noise_start <- function(start) {
  probability <- start$noise_probability
  if (!is.numeric(probability) || length(probability) != 1L ||
      !is.finite(probability) || probability <= 0 || probability >= 1) {
    stop(errorCondition(paste(
      "A start for a fit with `noise = TRUE` needs `noise_probability`, a single",
      "number strictly between zero and one, beside `profile_probabilities`;",
      "the two sum to one."), class = "latents_bad_start", call = NULL))
  }
  profile <- start$profile_probabilities
  if (is.numeric(profile) && is.matrix(profile) &&
      any(abs(rowSums(profile) + probability - 1) > 1e-8)) {
    stop(errorCondition(
      "`profile_probabilities` and `noise_probability` must sum to one.",
      class = "latents_bad_start", call = NULL))
  }
  start$noise_probability <- NULL
  if (is.numeric(profile) && is.matrix(profile)) {
    start$profile_probabilities <- profile / (1 - probability)
  }
  list(start = start, noise_probability = probability)
}

#' Add the noise component to a set of Gaussian parameters
#' @param parameters Parameters whose profile probabilities sum to one.
#' @param probability The noise proportion.
#' @param log_density The noise component's log density, `-log(V)`.
#' @return The parameters with one more, last, profile-probability column and
#'   `noise_log_density` set.
#' @noRd
.multilpa_add_noise <- function(parameters, probability, log_density) {
  parameters$profile_probabilities <- cbind(
    parameters$profile_probabilities * (1 - probability), probability)
  parameters$noise_log_density <- log_density
  parameters
}

#' Split a fitted noise component back out of the parameters and posteriors
#' @param parameters The fitted parameters, noise last.
#' @param posteriors The observation posteriors, noise last.
#' @return A list with the Gaussian `parameters` and `posteriors`, and the
#'   `noise_probability` and `noise_posteriors`.
#' @noRd
.multilpa_split_noise <- function(parameters, posteriors) {
  last <- ncol(posteriors)
  gaussian <- seq_len(last - 1L)
  noise_probability <- parameters$profile_probabilities[1L, last]
  parameters$profile_probabilities <-
    parameters$profile_probabilities[, gaussian, drop = FALSE]
  parameters$noise_log_density <- NULL
  list(parameters = parameters,
       posteriors = posteriors[, gaussian, drop = FALSE],
       noise_probability = unname(noise_probability),
       noise_posteriors = unname(posteriors[, last]),
       classes = {
         modal <- max.col(posteriors, ties.method = "first")
         modal[modal == last] <- 0L
         modal
       })
}

#' Every observation's class posteriors, the noise component included
#'
#' @param x A fitted model of this package.
#' @return The observation-by-profile posterior matrix, with a final `noise`
#'   column when the fit has a noise component. Its rows sum to one either way.
#' @noRd
.multilpa_class_posteriors <- function(x) {
  if (!isTRUE(x$noise)) return(x$subject_posteriors)
  cbind(x$subject_posteriors, noise = x$noise_posteriors)
}

#' Which observations have at least one observed indicator
#'
#' A row with every indicator missing contributes nothing to the likelihood,
#' so its posterior is not a classification of that row but the prior its
#' group gives it. The sample size of the individual-level criteria already
#' leaves such rows out; this is the same set of rows for every statistic of
#' classification quality.
#'
#' @param x A fitted model of this package.
#' @return A logical vector, one per observation; all `TRUE` for a fit that
#'   did not record it.
#' @noRd
.multilpa_informative_rows <- function(x) {
  n_rows <- NROW(x$subject_posteriors)
  # `drop_empty = FALSE` keeps the convention of 0.9.37 and earlier, in which
  # every row was counted.
  if (isFALSE(x$drop_empty)) return(rep(TRUE, n_rows))
  x$informative %||% rep(TRUE, n_rows)
}

#' The group sample size of the group-convention criteria
#' @param x A fitted model of this package.
#' @return The number of groups with an observed row, or every group under
#'   `drop_empty = FALSE` and for fits that did not record it.
#' @noRd
.multilpa_group_sample_size <- function(x) {
  if (isFALSE(x$drop_empty)) return(x$n_groups)
  x$n_informative_groups %||% x$n_groups
}

#' Validate the `drop_empty` switch
#' @noRd
.multilpa_check_drop_empty <- function(drop_empty) {
  if (!(isTRUE(drop_empty) || isFALSE(drop_empty))) {
    stop(errorCondition("`drop_empty` must be TRUE or FALSE.",
                        class = "latents_bad_argument", call = NULL))
  }
  invisible(NULL)
}

#' Which groups have at least one row with an observed indicator
#'
#' A group whose every row observes nothing is unclassified for the same
#' reason as such a row: its posterior is the group-class prior. Mplus leaves
#' such a group out of its clusters; so do the group sample size and the
#' group-level classification statistics here.
#'
#' @param x A fitted model of this package.
#' @return A logical vector, one per group (row of `group_posteriors`), or
#'   `NULL` for a fit without group posteriors.
#' @noRd
.multilpa_informative_groups <- function(x) {
  if (is.null(x$group_posteriors)) return(NULL)
  n_groups <- nrow(x$group_posteriors)
  if (isFALSE(x$drop_empty)) return(rep(TRUE, n_groups))
  informative <- .multilpa_informative_rows(x)
  if (all(informative) || is.null(x$group_index)) return(rep(TRUE, n_groups))
  observed <- rowsum(informative * 1, x$group_index, reorder = TRUE)
  present <- sort(unique(x$group_index))
  result <- rep(FALSE, n_groups)
  result[present] <- observed[, 1L] > 0
  result
}

#' The group posteriors of the groups that were classified
#' @param x A fitted model of this package.
#' @return `group_posteriors` restricted to informative groups, or `NULL`.
#' @noRd
.multilpa_classified_group_posteriors <- function(x) {
  if (is.null(x$group_posteriors)) return(NULL)
  x$group_posteriors[.multilpa_informative_groups(x), , drop = FALSE]
}

#' The class posteriors of the observations that were classified
#'
#' Entropy, average posterior probability and the entropy charge of ICL, AWE
#' and CLC describe how well observed values separate the classes, so they are
#' taken over the rows with at least one observed indicator. Class sizes are
#' estimates of the population and keep every row.
#'
#' @param x A fitted model of this package.
#' @return [.multilpa_class_posteriors()] restricted to informative rows.
#' @noRd
.multilpa_classified_posteriors <- function(x) {
  posteriors <- .multilpa_class_posteriors(x)
  if (is.null(posteriors)) return(NULL)
  posteriors[.multilpa_informative_rows(x), , drop = FALSE]
}

#' Refuse a verb that does not yet account for a noise component
#' @param x A fitted model of this package.
#' @param verb What was asked for, for the message.
#' @param class Extra condition classes to carry.
#' @return `NULL`, invisibly, for a fit without one.
#' @noRd
.multilpa_refuse_noise <- function(x, verb, class = character()) {
  if (!isTRUE(x$noise)) return(invisible(NULL))
  stop(errorCondition(sprintf(paste(
    "%s does not account for a noise component yet: its Gaussian profiles'",
    "posteriors do not sum to one, and the result would be for a different",
    "model. Refit with `noise = FALSE` for it."), verb),
    class = c("latents_unsupported_noise", class), call = NULL))
}
