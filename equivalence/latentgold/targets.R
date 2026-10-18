# This package's targets for one Latent GOLD case: what compare-functions.R
# compares Latent GOLD's output against. Sourced by make-kit.R, which writes
# them to targets/, and by test-latentgold.R, which builds them from the
# current package so that a regression is caught.

#' Posteriors of one level as one row per unit and one column per class
#'
#' @param long A long posterior table from `get_results()`.
#' @param unit The column naming the unit, `"row"` or `"group"`.
#' @param class The column naming the class, `"profile"` or `"group_class"`.
#' @return A data frame with `unit` and `class_1`, `class_2`, ...
.wide_posteriors <- function(long, unit, class) {
  classes <- sort(unique(long[[class]]))
  units <- unique(long[[unit]])
  columns <- lapply(classes, function(k) {
    rows <- long[[class]] == k
    long$posterior[rows][match(units, long[[unit]][rows])]
  })
  stopifnot("posterior table is not complete" = !anyNA(unlist(columns)))
  data.frame(unit = units, stats::setNames(columns, sprintf("class_%d", classes)))
}

#' Everything compare.R needs from this package for one case
#'
#' @param name The case's name.
#' @param case One element of [lg_cases()].
#' @param fit This package's fit to it.
#' @return A list of targets.
.targets <- function(name, case, fit) {
  # A target must be the optimum the search keeps finding, not one start's
  # luck; otherwise a disagreement with Latent GOLD is uninterpretable.
  starts <- get_results(fit, "starts")
  replicated <- sum(abs(starts$log_likelihood - fit$log_likelihood) < 1e-6,
                    na.rm = TRUE)
  if (replicated < 2L) {
    stop(errorCondition(sprintf(
      "Case `%s`: the best log likelihood was reached by %d start(s); a target needs at least 2.",
      name, replicated), class = "multilpa_unreplicated_target", call = NULL))
  }
  individual <- get_results(fit, "posteriors", format = "wide")
  profile_columns <- grep("^posterior_profile_", names(individual), value = TRUE)
  posteriors <- data.frame(id = case$data$id[individual$row],
                           stats::setNames(individual[profile_columns],
                                           sub("^posterior_profile_", "class_", profile_columns)))
  groups <- if (fit$n_group_classes > 1L) {
    wide <- .wide_posteriors(get_results(fit, "group_posteriors"), "group", "group_class")
    names(wide)[names(wide) == "unit"] <- "g"
    wide
  }
  inference <- if (name %in% c("c01_continuous_varying", "c02_continuous_equal",
                               "c07_covariates")) {
    parameter_inference(fit)
  }
  steps <- if (identical(name, "c09_three_step")) {
    list(classification_errors = get_results(fit, "classification_errors",
                                             level = "individuals"),
         distal_bch = three_step(fit, case$data, "distal", method = "bch"),
         distal_modal = three_step(fit, case$data, "distal", method = "modal"),
         distal_proportional = three_step(fit, case$data, "distal",
                                          method = "proportional"),
         # The retained listing reports model-based information SEs. The
         # package's default is now robust, a different covariance estimator.
         covariate_ml = r3step(fit, case$data, "x", vcov_type = "observed"))
  }
  list(case = name, description = case$description, anchors = case$anchors,
       lgs = names(case$files), data_file = case$data_file,
       continuous = fit$continuous, variance_model = fit$variance_model,
       covariance_model = fit$covariance_model %||% "diagonal",
       missing = fit$missing %||% "error",
       log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters,
       n_best_replicated = replicated,
       criteria = get_results(fit, "information_criteria"),
       n_observations = fit$n_observations, n_groups = fit$n_groups,
       posteriors = posteriors, group_posteriors = groups,
       profiles = get_results(fit, "profiles"),
       group_class_probabilities = if (fit$n_group_classes > 1L) fit$group_probabilities,
       parameters = inference,
       residuals = if (identical(name, "c08_bivariate_residuals"))
         get_results(fit, "residuals"),
       steps = steps,
       data = case$data,
       multilpa_version = as.character(utils::packageVersion("latents")))
}

