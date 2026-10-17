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
#'   `classes`, `family`, `inverse_link`, `response_mean` (a function of
#'   linear predictors and their classes giving the outcome's mean: the
#'   inverse link, or for the ordinal family the expected category score
#'   under that class's thresholds), `weights` (per person, or `NULL`).
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
  inverse_link <- switch(family, gaussian = identity, binomial = stats::plogis,
                         poisson = exp, negative_binomial = exp)
  response_mean <- if (identical(family, "ordinal")) {
    function(linear, classes) {
      scores <- .ordinal_expected_scores(
        matrix(linear, length(linear), spec$n_classes), x$params$thresholds)
      scores[cbind(seq_along(linear), classes)]
    }
  } else function(linear, classes) inverse_link(linear)
  list(growth = growth, spec = spec, posterior = posterior,
       modal = max.col(posterior, ties.method = "first"),
       classes = paste0("class_", seq_len(spec$n_classes)), family = family,
       inverse_link = inverse_link, response_mean = response_mean,
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
#' link scale and carried to the response scale. An ordinal class's mean (the
#' expected category score) also depends on its thresholds, so its band is
#' the delta method's on the score scale and `std_error` is on that scale.
#'
#' @param x A trajectory fit.
#' @param rows A data frame of predictor values.
#' @param inference From `.trajectory_inference()`.
#' @param level Confidence level.
#' @return A data frame: `class`, `estimate`, `std_error` (response scale),
#'   `conf_low`, `conf_high`, one row per class and row of `rows`.
#' @noRd
.trajectory_band <- function(x, rows, inference, level) {
  view <- .trajectory_view(x)
  spec <- view$spec
  design <- .growth_design(x, rows)
  z <- .inference_critical(level, "one_minus")
  if (identical(view$family, "ordinal")) {
    return(.ordinal_trajectory_band(x, design, inference, z, view$classes))
  }
  do.call(rbind, lapply(view$classes, function(class) {
    names <- c(sprintf("coefficient.%s.%s", class, colnames(spec$x)),
               if (ncol(spec$z) > 0L) sprintf("coefficient.common.%s", colnames(spec$z)))
    estimate <- as.vector(design %*% inference$theta[names]) + attr(design, "offset")
    covariance <- inference$vcov[names, names, drop = FALSE]
    std_error <- .inference_delta_se(design, covariance, "rowsums")
    response <- view$inverse_link(estimate)
    response_derivative <- switch(view$family, gaussian = rep(1, length(estimate)),
                                  binomial = response * (1 - response),
                                  poisson = , negative_binomial = response)
    data.frame(class = class, estimate = response,
               std_error = std_error * response_derivative,
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

#' Check a trajectory fit can take a three-step correction
#'
#' A multilevel growth mixture is refused (its persons are dependent within
#' clusters), as is a fit with membership covariates (the correction has no
#' classification-error adjustment conditional on them), and `data` must be
#' the data frame the model was fitted to.
#'
#' @param x A trajectory fit.
#' @param data The data frame passed to the three-step verb.
#' @return The fit's `.trajectory_view()`.
#' @noRd
.trajectory_check_three_step <- function(x, data) {
  view <- .trajectory_view(x)
  spec <- view$spec
  if (!is.null(spec$cluster)) {
    stop(errorCondition(paste(
      "Three-step correction is unavailable for a multilevel growth mixture:",
      "its persons are dependent within clusters, and the correction treats",
      "them as independent."), class = "latents_unsupported_three_step", call = NULL))
  }
  if (!spec$intercept_only) {
    stop(errorCondition(paste(
      "Three-step correction is unavailable for a fit that already uses",
      "membership covariates: it has no classification-error adjustment",
      "conditional on them. Fit the trajectory model without `membership`."),
      class = "latents_unsupported_three_step", call = NULL))
  }
  .trajectory_check_data(x, data)
  view
}

#' Require the data frame a trajectory model was fitted to
#' @param x A trajectory fit.
#' @param data A data frame.
#' @return `NULL`, invisibly, or `latents_bad_inference_data`.
#' @noRd
.trajectory_check_data <- function(x, data) {
  spec <- x$spec
  if (!is.data.frame(data) || !spec$id %in% names(data) ||
      nrow(data) < max(spec$kept_rows) ||
      !.multilpa_same_values(data[[spec$id]][spec$kept_rows], spec$model_data[[spec$id]])) {
    stop(errorCondition(
      "`data` must be the data frame the trajectory model was fitted to.",
      class = "latents_bad_inference_data", call = NULL))
  }
  invisible(NULL)
}

#' The person-level assignments of a trajectory fit, as three-step reads them
#' @return A list as `.multilpa_level_assignments()` builds.
#' @noRd
.trajectory_pieces <- function(view) {
  list(posteriors = view$posterior, modal = view$modal,
       units = view$spec$group_levels, group_index = seq_len(view$spec$n_groups),
       n_classes = view$spec$n_classes)
}

#' Whether a variable is part of a trajectory model's own regression
#' @noRd
.trajectory_uses <- function(spec, variable) {
  identical(variable, spec$response_name) ||
    variable %in% all.vars(stats::delete.response(spec$terms))
}

#' Three-step inputs of a trajectory model: person posteriors and one
#' distal outcome value per person
#'
#' The outcome must be constant within persons (it is measured once, or is a
#' person-level quantity) and must not be the trajectory outcome; the fit and
#' `data` are checked by `.trajectory_check_three_step()`.
#'
#' @return A list with `pieces` (as `.multilpa_level_assignments()` builds)
#'   and `values`.
#' @noRd
.trajectory_three_step_inputs <- function(x, data, outcome) {
  view <- .trajectory_check_three_step(x, data)
  spec <- view$spec
  if (.trajectory_uses(spec, outcome)) {
    stop(errorCondition(sprintf(
      "`%s` is part of the trajectory model, so it is not a distal outcome.",
      outcome), class = c("multilpa_indicator_reused", "latents_bad_outcome"),
      call = NULL))
  }
  by_person <- split(data[[outcome]][spec$kept_rows], spec$group_index)
  constant <- vapply(by_person, function(v) length(unique(v)) == 1L, logical(1))
  if (!all(constant)) {
    stop(errorCondition(sprintf(paste(
      "`%s` varies within persons; a distal outcome of a trajectory model is",
      "one value per person."), outcome),
      class = "latents_bad_outcome", call = NULL))
  }
  list(pieces = .trajectory_pieces(view),
       values = vapply(by_person, function(v) v[[1L]], numeric(1), USE.NAMES = FALSE))
}

#' R3STEP inputs of a trajectory model: person posteriors and covariates
#'
#' Each covariate must be constant within persons, because it predicts a
#' person's class, and must not be part of the trajectory model.
#'
#' @return A list with `pieces` and `design` (persons by intercept and
#'   covariates).
#' @noRd
.trajectory_r3step_inputs <- function(x, data, covariates) {
  view <- .trajectory_check_three_step(x, data)
  spec <- view$spec
  reused <- covariates[vapply(covariates, .trajectory_uses, logical(1), spec = spec)]
  if (length(reused) > 0L) {
    stop(errorCondition(sprintf(
      "%s already part of the trajectory model, so not an external predictor.",
      paste(sprintf("`%s` is", reused), collapse = ", ")),
      class = c("multilpa_indicator_reused", "latents_bad_covariate"), call = NULL))
  }
  values <- vapply(covariates, function(covariate) {
    by_person <- split(data[[covariate]][spec$kept_rows], spec$group_index)
    if (!all(vapply(by_person, function(v) length(unique(v)) == 1L, logical(1)))) {
      stop(errorCondition(sprintf(paste(
        "`%s` varies within persons; a covariate of trajectory class",
        "membership is one value per person."), covariate),
        class = "latents_bad_covariate", call = NULL))
    }
    vapply(by_person, function(v) v[[1L]], numeric(1), USE.NAMES = FALSE)
  }, numeric(spec$n_groups)) |> matrix(spec$n_groups, length(covariates))
  colnames(values) <- covariates
  list(pieces = .trajectory_pieces(view),
       design = .multilpa_r3step_check_rank(cbind(`(Intercept)` = 1, values)))
}
