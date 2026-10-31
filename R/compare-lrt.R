# The likelihood-ratio test of compare_models(test = "lrt"): when one fit is
# a restriction of the next, and the chi-square reference for that pair.

#' Add likelihood-ratio tests to a comparison of profile and class models
#'
#' Each row is tested against the row before it, as [stats::anova()] does: the
#' earlier model must be nested in the later one, with the same numbers of
#' profiles and group classes, so the restriction holds parameters at interior
#' values and twice the log-likelihood difference has a chi-square reference
#' with the difference in free parameters as degrees of freedom (Wilks, 1938).
#'
#' @param table A `latents_comparison` with its `models` attribute.
#' @return `table` with `lr_statistic`, `lr_df` and `lr_p_value` (`NA` in the
#'   first row) and three more displayed columns.
#' @noRd
.latents_comparison_lrt <- function(table) {
  models <- attr(table, "models")
  labels <- names(models)
  lapply(models, function(fit) {
    if (.latents_is_weighted(fit)) {
      .latents_refuse_weights(paste(
        "`compare_models(test = \"lrt\")`: under pseudo maximum likelihood",
        "twice the log-likelihood difference is not chi-square distributed.",
        "Compare weighted fits on their information criteria"))
    }
    .multilpa_refuse_noise(fit, "`compare_models(test = \"lrt\")`")
    if (!is.null(fit$prior)) {
      stop(errorCondition(paste(
        "`compare_models(test = \"lrt\")` needs maximum-likelihood fits; a fit",
        "with `prior` is estimated at a posterior mode. Refit without `prior`."),
        class = "latents_unsupported_prior", call = NULL))
    }
  })
  converged <- vapply(models, function(fit) isTRUE(fit$converged), logical(1))
  if (!all(converged)) {
    stop(errorCondition(sprintf(paste(
      "%s did not converge, so the likelihood-ratio test would compare",
      "unfinished fits."), paste(labels[!converged], collapse = ", ")),
      class = "latents_no_converge", call = NULL))
  }
  tests <- lapply(seq_along(models)[-1L], function(i) {
    .latents_lrt_pair(models[[i - 1L]], models[[i]], labels[i - 1L], labels[i])
  })
  table$lr_statistic <- c(NA_real_, vapply(tests, `[[`, numeric(1), "statistic"))
  table$lr_df <- c(NA_integer_, vapply(tests, `[[`, integer(1), "df"))
  table$lr_p_value <- stats::pchisq(table$lr_statistic, table$lr_df, lower.tail = FALSE)
  bound <- labels[vapply(models, function(fit) isTRUE(fit$boundary), logical(1))]
  if (length(bound) > 0L) {
    warning(warningCondition(sprintf(paste(
      "%s %s a variance at `min_variance`: a parameter on its bound can make",
      "the chi-square reference of the likelihood-ratio test inaccurate."),
      paste(bound, collapse = ", "), if (length(bound) > 1L) "have" else "has"),
      class = "latents_boundary", call = NULL))
  }
  column <- .latents_column
  attr(table, "display") <- c(attr(table, "display"), list(
    column("LR", "lr_statistic"), column("df", "lr_df", "integer"),
    column("p", "lr_p_value", "p")))
  table
}

#' The likelihood-ratio test of one model against the next
#'
#' @param null,alternative Consecutive fits of the comparison.
#' @param null_label,alternative_label Their row labels.
#' @return A list with `statistic` (twice the log-likelihood difference, zero
#'   when the alternative is lower only by optimization tolerance) and `df`.
#'   Raises `latents_bad_nesting` for a pair that is not nested, differs in
#'   class counts or does not add parameters, and
#'   `latents_reversed_likelihood` when the alternative's likelihood is lower
#'   beyond tolerance.
#' @noRd
.latents_lrt_pair <- function(null, alternative, null_label, alternative_label) {
  refuse <- function(reason) {
    stop(errorCondition(sprintf(paste(
      "%s is not nested in %s for a likelihood-ratio test: %s. List the models",
      "from the most restricted to the most general."),
      null_label, alternative_label, reason),
      class = "latents_bad_nesting", call = NULL))
  }
  transition <- c("multilpa_transitions", "multilpa_lta")
  counts <- function(fit) {
    if (inherits(fit, transition)) {
      arguments <- .lta_arguments(fit)
      c(arguments$n_profiles, arguments$n_group_classes)
    } else c(fit$n_profiles, fit$n_group_classes %||% 1L)
  }
  if (!identical(as.integer(counts(null)), as.integer(counts(alternative)))) {
    refuse(paste(
      "they differ in the number of profiles or group classes, where the",
      "chi-square reference does not hold; use bootstrap_lrt() or",
      "get_results(grid, \"lmr\")"))
  }
  reason <- if (inherits(null, transition) && inherits(alternative, transition)) {
    if (!isTRUE(.lta_arguments(null)$mover_stayer) &&
        isTRUE(.lta_arguments(alternative)$mover_stayer)) {
      "a stayer class adds a mixing proportion whose null value is on its bound; use bootstrap_lrt()"
    } else if (!.lta_nested(null, alternative)) {
      "every option of the first must be at most as free as in the second"
    }
  } else if (!inherits(null, transition) && !inherits(alternative, transition)) {
    .multilpa_lrt_nesting(null, alternative)
  } else {
    "a transition model and a profile model are different models of the data"
  }
  if (!is.null(reason)) refuse(reason)
  df <- as.integer(alternative$n_parameters - null$n_parameters)
  if (df <= 0L) {
    refuse(sprintf("the second estimates %d free parameters against %d",
                   alternative$n_parameters, null$n_parameters))
  }
  tolerance <- max(.latents_refit_controls(null)$tol,
                   .latents_refit_controls(alternative)$tol)
  window <- 100 * tolerance * (1 + abs(null$log_likelihood))
  statistic <- 2 * (alternative$log_likelihood - null$log_likelihood)
  if (statistic < -window) {
    stop(errorCondition(sprintf(paste(
      "%s has a lower likelihood than %s, in which it is nested, by %.3g:",
      "refit %s with more starts first."), alternative_label, null_label,
      -statistic / 2, alternative_label),
      class = "latents_reversed_likelihood", call = NULL))
  }
  list(statistic = max(0, statistic), df = df)
}

#' Why one profile or class model is not a restriction of another
#'
#' The indicators, their coding and the estimation settings must agree. What
#' the null may restrict: the covariance structure, letter by letter
#' (volume, shape, orientation; `I` inside `E` inside `V`); covariances of
#' chosen pairs (a subset of the alternative's pairs, or all of them inside a
#' full covariance that extends the same diagonal); joint pairs of discrete
#' items; direct effects; membership covariates; profile slopes shared across
#' group classes; and count dispersions shared across profiles.
#'
#' @param null,alternative Fits of class `multilpa` or `multilpa_covariates`.
#' @return `NULL` when `null` is nested in `alternative`, otherwise the
#'   reason as one phrase.
#' @noRd
.multilpa_lrt_nesting <- function(null, alternative) {
  # A pair modelled jointly is one categorical indicator named "a:b"; as
  # data it is still the two nominal items.
  nominal <- function(fit) {
    coded <- fit$categorical %||% character()
    pairs <- fit$local_dependence
    joint <- if (is.null(pairs)) rep(NA_integer_, length(coded)) else
      match(coded, pairs$name[pairs$kind == "categorical"])
    joint_pairs <- if (is.null(pairs)) pairs else pairs[pairs$kind == "categorical", ]
    sort(c(coded[is.na(joint)],
           unlist(joint_pairs[joint[!is.na(joint)], c("item_1", "item_2")],
                  use.names = FALSE)))
  }
  setting <- function(fit) {
    list(vars = fit$original_vars %||% fit$vars, id = fit$id,
         single_level = isTRUE(fit$single_level),
         categorical = nominal(fit),
         ordinal = sort(fit$ordinal %||% character()),
         count = sort(fit[["count"]] %||% character()),
         count_model = fit$extra_data$count_model %||% "poisson",
         missing = fit$missing %||% "error", centering = fit$centering %||% "none",
         scaling = fit$scaling %||% "none", min_variance = fit$min_variance,
         min_probability = fit$min_probability, drop_empty = fit$drop_empty %||% TRUE)
  }
  if (!identical(setting(null), setting(alternative))) {
    return(paste(
      "they differ in their indicators, how the indicators are coded, or an",
      "estimation setting (missing data, centering, scaling, bounds)"))
  }
  held <- function(fit) fit$fixed %||% character()
  if (length(held(null)) > 0L || length(held(alternative)) > 0L) {
    return("a fit with `fixed` measurement blocks has a conditional likelihood")
  }
  if (!.multilpa_lrt_covariance_nested(null, alternative)) {
    return("the covariance structure of the first is not a restriction of the second's")
  }
  pairs <- function(fit, kind) {
    pairs <- fit$local_dependence
    if (is.null(pairs) || nrow(pairs) == 0L) return(character())
    pairs$name[pairs$kind == kind]
  }
  if (!all(pairs(null, "categorical") %in% pairs(alternative, "categorical"))) {
    return("the first models a pair of items jointly that the second does not")
  }
  included <- function(field) {
    all((null[[field]] %||% character()) %in% (alternative[[field]] %||% character()))
  }
  if (!all(names(null$direct_effects) %in% names(alternative$direct_effects))) {
    return("the first has a direct effect the second does not")
  }
  if (!included("profile_covariates") || !included("group_covariates")) {
    return("the first has a membership covariate the second does not")
  }
  rank <- function(value, order) match(value, order)
  if (rank(null$profile_slopes %||% "shared", c("shared", "varying")) >
      rank(alternative$profile_slopes %||% "shared", c("shared", "varying"))) {
    return("the first lets profile slopes vary by group class and the second does not")
  }
  if (rank(null$extra_data$count_dispersion %||% "varying", c("equal", "varying")) >
      rank(alternative$extra_data$count_dispersion %||% "varying", c("equal", "varying"))) {
    return("the first lets count dispersions vary by profile and the second does not")
  }
  NULL
}

#' Is one fit's Gaussian covariance a restriction of another's
#'
#' Letter by letter for the structures, with `I` inside `E` inside `V`. A fit
#' with covariances for chosen pairs extends its diagonal (EEI or VVI): it
#' contains that diagonal and every subset of its pairs, and is contained in a
#' larger set of pairs on a diagonal at least as free, and in the full
#' structures that contain its covariance (EEE and up for one shared matrix,
#' VVV for one per profile).
#' @param null,alternative Fits of class `multilpa` or `multilpa_covariates`.
#' @return One logical value.
#' @noRd
.multilpa_lrt_covariance_nested <- function(null, alternative) {
  continuous <- function(fit) length(.multilpa_continuous_names(fit)) > 0L
  if (!continuous(null) && !continuous(alternative)) return(TRUE)
  letters_within <- function(inner, outer) {
    all(vapply(seq_len(3L), function(position) {
      rank <- function(code) match(substr(code, position, position), c("I", "E", "V"))
      rank(inner) <= rank(outer)
    }, logical(1)))
  }
  diagonal <- function(fit) if (identical(fit$variance_model, "equal")) "EEI" else "VVI"
  pairs <- function(fit) {
    if (is.null(fit$covariance_pattern)) character() else
      .latents_pattern_pairs(fit$covariance_pattern)
  }
  structure <- function(fit) {
    if (is.null(fit$covariance_pattern)) .latents_comparison_structure(fit) else diagonal(fit)
  }
  if (length(pairs(null)) == 0L) {
    return(letters_within(structure(null), structure(alternative)))
  }
  if (length(pairs(alternative)) > 0L) {
    return(letters_within(diagonal(null), diagonal(alternative)) &&
             all(pairs(null) %in% pairs(alternative)))
  }
  full <- if (identical(diagonal(null), "EEI")) "EEE" else "VVV"
  letters_within(full, structure(alternative))
}
