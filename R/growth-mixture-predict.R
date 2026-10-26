# predict() for growth mixture models: class trajectories of new or fitted
# persons, their class posteriors given their outcomes, and their own curves
# (empirical Bayes random effects), as lcmm::predictY, predictClass and
# predictRE define them (Proust-Lima, Philipps and Liquet 2017).

#' Rebuild a growth specification on new persons
#'
#' The fitted terms, factor levels and membership codings evaluate `newdata`,
#' so its design columns follow the fitted order. Every person in `newdata` is
#' a new person, and in a multilevel growth mixture every cluster is a new
#' cluster, classified from its own persons.
#'
#' @param object A `latents_growth_mixture` fit.
#' @param newdata A data frame with the predictors, `id` (and `cluster`).
#' @param with_response Whether the outcome must be present.
#' @return A growth specification for the engine, without sampling weights.
#' @noRd
.growth_new_spec <- function(object, newdata, with_response) {
  spec <- object$spec
  new <- .mixture_prediction_data(.mixture_new_spec(object, newdata, with_response))
  new$random_design <- .growth_random_design(spec, newdata)
  new$sampling_weights <- NULL
  new$row_weights <- NULL
  if (.growth_is_correlated(spec)) {
    times <- newdata[[spec$residual_time]]
    if (is.null(times) || anyNA(times) || any(!is.finite(times))) {
      stop(errorCondition(sprintf(
        "`newdata` needs the residual time column `%s`, finite on every row.",
        spec$residual_time), class = "latents_bad_data", call = NULL))
    }
    if (identical(spec$residual_correlation, "brownian") && any(times < 0)) {
      stop(errorCondition(
        "Brownian residual times in `newdata` must be nonnegative.",
        class = "latents_bad_data", call = NULL))
    }
    new$residual_times <- times
  }
  if (!is.null(spec$cluster)) {
    if (!spec$cluster %in% names(newdata) || anyNA(newdata[[spec$cluster]])) {
      stop(errorCondition(sprintf(
        "`newdata` needs the cluster column `%s`, without missing values.",
        spec$cluster), class = "latents_bad_data", call = NULL))
    }
    blueprint <- spec$group_membership_design
    frame <- .mixture_prediction_data(stats::model.frame(
      blueprint$terms, newdata, xlev = blueprint$xlevels, na.action = stats::na.fail))
    row_group <- stats::model.matrix(blueprint$terms, frame,
                                     contrasts.arg = blueprint$contrasts)
    row_membership <- new$w[new$group_index, , drop = FALSE]
    new <- utils::modifyList(new, .growth_cluster_spec(
      newdata[[spec$cluster]], new$group_index, row_membership, row_group,
      spec$n_cluster_classes, check_count = FALSE))
  }
  new
}

#' Random-effects design of new rows in the fitted coding
#' @param spec The fitted growth specification.
#' @param newdata Rows on which to evaluate the random formula.
#' @return A finite matrix with the fitted columns and contrasts.
#' @noRd
.growth_random_design <- function(spec, newdata) {
  random_terms <- spec$random_terms %||%
    stats::terms(stats::model.frame(spec$random, spec$model_data))
  xlevels <- spec$random_xlevels %||% stats::.getXlevels(
    random_terms, stats::model.frame(random_terms, spec$model_data))
  random_frame <- .mixture_prediction_data(
    stats::model.frame(random_terms, newdata, xlev = xlevels,
                       na.action = stats::na.fail))
  design <- stats::model.matrix(random_terms, random_frame,
    contrasts.arg = attr(spec$random_design, "contrasts"))
  if (!identical(colnames(design), colnames(spec$random_design)) ||
      any(!is.finite(design))) {
    stop(errorCondition(paste(
      "The random-effects design of `newdata` must be finite and have the",
      "fitted columns."), class = "latents_bad_data", call = NULL))
  }
  design
}

#' Each row's mean under every class (persons' random effects at zero)
#' @return A rows-by-classes matrix.
#' @noRd
.growth_class_means <- function(spec, params) {
  design <- cbind(spec$x, spec$z)
  vapply(seq_len(spec$n_classes), function(k) {
    as.vector(design %*% .growth_coefficients(params, k))
  }, numeric(spec$n)) |> matrix(spec$n, spec$n_classes) + spec$offset
}

#' @rdname latents-predict
#' @export
predict.latents_growth_mixture <- function(object, newdata = NULL,
                                           type = c("response", "class_response",
                                                    "posterior", "conditional",
                                                    "class_conditional",
                                                    "random_effects"), ...) {
  type <- match.arg(type)
  .latents_warn_unused_arguments(list(...), "predict()",
    "A growth mixture prediction takes `newdata` and `type`.")
  with_response <- type %in% c("posterior", "conditional", "class_conditional",
                               "random_effects")
  spec <- if (is.null(newdata)) object$spec else
    .growth_new_spec(object, newdata, with_response)
  params <- object$params
  classes <- .growth_class_names(object)
  means <- .growth_class_means(spec, params)
  person <- spec$group_index
  if (type %in% c("response", "class_response")) {
    priors <- .growth_class_prior(spec, params)[person, , drop = FALSE]
    if (identical(type, "response")) {
      return(data.frame(row = spec$kept_rows, fitted = rowSums(priors * means)))
    }
    return(data.frame(row = rep(spec$kept_rows, spec$n_classes),
                      class = rep(classes, each = spec$n),
                      prior = as.vector(priors), fitted = as.vector(means),
                      stringsAsFactors = FALSE))
  }
  # New persons are classified in the fit's standard coordinates and carried
  # back, as the fitted persons were.
  expectation <- if (is.null(newdata)) object$expectation else
    .growth_fitted_expectation(object, spec)
  if (identical(type, "posterior")) {
    refit <- object
    refit$spec <- spec
    refit$expectation <- expectation
    return(.growth_assignment_table(refit))
  }
  posterior <- expectation$posterior
  if (identical(type, "random_effects")) {
    terms <- colnames(spec$random_design)
    effects <- do.call(rbind, expectation$means)
    out <- data.frame(id = rep(spec$group_levels, spec$n_classes),
                      class = rep(classes, each = spec$n_groups),
                      posterior = as.vector(posterior),
                      matrix(effects, ncol = length(terms)),
                      stringsAsFactors = FALSE)
    names(out) <- c(spec$id, "class", "posterior", terms)
    return(out)
  }
  # A person's own curve under class k adds E(b | y, class k), the class's
  # posterior mean of their random effects, to the class trajectory.
  own <- vapply(seq_len(spec$n_classes), function(k) {
    rowSums(spec$random_design * expectation$means[[k]][person, , drop = FALSE])
  }, numeric(spec$n)) |> matrix(spec$n, spec$n_classes) + means
  row_posterior <- posterior[person, , drop = FALSE]
  if (identical(type, "conditional")) {
    return(data.frame(row = spec$kept_rows, fitted = rowSums(row_posterior * own)))
  }
  data.frame(row = rep(spec$kept_rows, spec$n_classes),
             class = rep(classes, each = spec$n),
             posterior = as.vector(row_posterior), fitted = as.vector(own),
             stringsAsFactors = FALSE)
}
