# What trajectory tables and plots read from a fit, for both trajectory
# models: mixture_regression() with persons as the classified units
# (`id`, `class_level = "group"`; latent class growth analysis) and its
# growth mixture extension (`random =`).

#' The pieces of a trajectory fit the tables and plots read
#'
#' @param x A `latents_growth_mixture` fit, or a `latents_mixture_regression`
#'   fit that classifies persons.
#' @return A list: `growth` (whether there are random effects), `spec`,
#'   `posterior` (persons by classes, rows summing to one), `modal`,
#'   `classes`, `family`, `inverse_link`, `weights` (per person, or `NULL`).
#' @noRd
.trajectory_view <- function(x) {
  growth <- inherits(x, "latents_growth_mixture")
  if (!growth && !(inherits(x, "latents_mixture_regression") &&
                   identical(x$spec$nesting, "group"))) {
    stop(errorCondition(paste(
      "Trajectory views need persons as the classified units: fit with `id`",
      "and `class_level = \"group\"`."), class = "latents_bad_argument", call = NULL))
  }
  spec <- x$spec
  posterior <- if (growth) x$expectation$posterior else x$expectation$group_tau
  family <- spec$family %||% "gaussian"
  list(growth = growth, spec = spec, posterior = posterior,
       modal = max.col(posterior, ties.method = "first"),
       classes = paste0("class_", seq_len(spec$n_classes)), family = family,
       inverse_link = switch(family, gaussian = identity, binomial = stats::plogis,
                             poisson = exp),
       weights = spec$sampling_weights)
}

#' The estimates and their covariance, for either kind of trajectory fit
#' @noRd
.trajectory_inference <- function(x, vcov_type = NULL) {
  if (inherits(x, "latents_growth_mixture")) .growth_resolve_inference(x, vcov_type) else
    .mixture_resolve_inference(x, vcov_type)
}

#' The outcome on the response scale (a proportion for binomial counts)
#' @noRd
.trajectory_observed <- function(view) {
  if (identical(view$family, "binomial")) view$spec$y / view$spec$trials else view$spec$y
}

#' Each class's mean trajectory at given rows, with a Wald band
#'
#' The class's linear predictor is linear in its coefficients, so its
#' standard error is exact given their covariance; the band is formed on the
#' link scale and carried to the response scale.
#'
#' @param x A trajectory fit.
#' @param rows A data frame of predictor values.
#' @param inference From `.trajectory_inference()`.
#' @param level Confidence level.
#' @return A data frame: `class`, `estimate`, `std_error` (link scale),
#'   `conf_low`, `conf_high`, one row per class and row of `rows`.
#' @noRd
.trajectory_band <- function(x, rows, inference, level) {
  view <- .trajectory_view(x)
  spec <- view$spec
  design <- .growth_design(x, rows)
  z <- stats::qnorm(1 - (1 - level) / 2)
  do.call(rbind, lapply(view$classes, function(class) {
    names <- c(sprintf("coefficient.%s.%s", class, colnames(spec$x)),
               if (ncol(spec$z) > 0L) sprintf("coefficient.common.%s", colnames(spec$z)))
    estimate <- as.vector(design %*% inference$theta[names])
    covariance <- inference$vcov[names, names, drop = FALSE]
    std_error <- sqrt(pmax(rowSums((design %*% covariance) * design), 0))
    data.frame(class = class, estimate = view$inverse_link(estimate),
               std_error = std_error,
               conf_low = view$inverse_link(estimate - z * std_error),
               conf_high = view$inverse_link(estimate + z * std_error),
               stringsAsFactors = FALSE)
  }))
}

#' Each class's coefficients (class-specific then common)
#' @noRd
.trajectory_coefficients <- function(x, k) c(x$params$beta[, k], x$params$common)

#' Whether a fit is a trajectory model (persons classified)
#' @noRd
.trajectory_is_fit <- function(x) {
  inherits(x, "latents_growth_mixture") ||
    (inherits(x, "latents_mixture_regression") && identical(x$spec$nesting, "group"))
}

#' Three-step inputs of a trajectory model: person posteriors and one
#' distal outcome value per person
#'
#' The outcome must be constant within persons (it is measured once, or is a
#' person-level quantity), must not be the trajectory outcome, and `data`
#' must be the data frame the model was fitted to. A fit with membership
#' covariates is refused, as for profile models: the correction has no
#' classification-error adjustment conditional on those covariates.
#'
#' @return A list with `pieces` (as `.multilpa_level_assignments()` builds)
#'   and `values`.
#' @noRd
.trajectory_three_step_inputs <- function(x, data, outcome) {
  view <- .trajectory_view(x)
  spec <- view$spec
  if (!spec$intercept_only) {
    stop(errorCondition(paste(
      "Three-step correction is unavailable for a fit that already uses",
      "membership covariates: it has no classification-error adjustment",
      "conditional on them. Fit the trajectory model without `membership`."),
      class = "latents_unsupported_three_step", call = NULL))
  }
  if (identical(outcome, spec$response_name) ||
      outcome %in% all.vars(stats::delete.response(spec$terms))) {
    stop(errorCondition(sprintf(
      "`%s` is part of the trajectory model, so it is not a distal outcome.",
      outcome), class = c("multilpa_indicator_reused", "latents_bad_outcome"),
      call = NULL))
  }
  if (nrow(data) < max(spec$kept_rows) ||
      !.multilpa_same_values(data[[spec$id]][spec$kept_rows], spec$model_data[[spec$id]])) {
    stop(errorCondition(
      "`data` must be the data frame the trajectory model was fitted to.",
      class = "latents_bad_inference_data", call = NULL))
  }
  by_person <- split(data[[outcome]][spec$kept_rows], spec$group_index)
  constant <- vapply(by_person, function(v) length(unique(v)) == 1L, logical(1))
  if (!all(constant)) {
    stop(errorCondition(sprintf(paste(
      "`%s` varies within persons; a distal outcome of a trajectory model is",
      "one value per person."), outcome),
      class = "latents_bad_outcome", call = NULL))
  }
  list(pieces = list(posteriors = view$posterior, modal = view$modal,
                     units = spec$group_levels,
                     group_index = seq_len(spec$n_groups),
                     n_classes = spec$n_classes),
       values = vapply(by_person, function(v) v[[1L]], numeric(1), USE.NAMES = FALSE))
}
