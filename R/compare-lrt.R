# The likelihood-ratio test of compare_models(test = "lrt"): when one fit is
# a restriction of the next, and the chi-square reference for that pair.

#' Add likelihood-ratio tests to a comparison of profile and class models
#'
#' Each row is tested against the row before it, as [stats::anova()] does: the
#' earlier model must be nested in the later one, with the same numbers of
#' profiles and group classes, and the restriction must be regular (see
#' `.latents_lrt_singular_covariance()`), so twice the log-likelihood
#' difference has a chi-square reference with the difference in free
#' parameters as degrees of freedom.
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
  bound <- labels[vapply(models, .latents_boundary_active, logical(1))]
  if (length(bound) > 0L) {
    warning(warningCondition(sprintf(paste(
      "%s %s a parameter at a numerical bound: a boundary can make",
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
#'   class counts or does not add parameters;
#'   `latents_singular_restriction` (also a `latents_bad_nesting`) for a
#'   nested covariance restriction that is singular at the null; and
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
    } else if (!.lta_lrt_designs_nested(null, alternative)) {
      paste("a transition or initial-state covariate of the first has values",
            "the second's covariates do not reproduce")
    }
  } else if (!inherits(null, transition) && !inherits(alternative, transition)) {
    .multilpa_lrt_nesting(null, alternative)
  } else {
    "a transition model and a profile model are different models of the data"
  }
  if (!is.null(reason)) refuse(reason)
  if (.latents_lrt_singular_covariance(null, alternative)) {
    stop(errorCondition(sprintf(paste(
      "%s (%s) is nested in %s (%s), but the restriction is singular: when",
      "every profile is spherical, the orientation %s adds has no value, so",
      "twice the log-likelihood difference has no chi-square reference.",
      "Compare these two structures on their information criteria."),
      null_label, .latents_lrt_structure(null), alternative_label,
      .latents_lrt_structure(alternative), alternative_label),
      class = c("latents_singular_restriction", "latents_bad_nesting"), call = NULL))
  }
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
#' group classes; and count dispersions shared across profiles. Ordinal
#' indicators must have the same categories in the same order. A predictor
#' both fits name must enter both with the same values: each design of the
#' first (profile membership, group-class membership, the direct effects on
#' each indicator) must lie in the column space of the second's, so that
#' setting the added coefficients to zero gives the first model back.
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
  # Ordinal categories in their fitted order: the order sets the scores the
  # ordinal likelihood gives them, so the same labels in another order are a
  # different measurement model.
  ordinal_levels <- function(fit) {
    levels <- fit$extra_data$ordinal_levels
    if (length(levels) == 0L) NULL else levels[sort(names(levels))]
  }
  setting <- function(fit) {
    list(vars = fit$original_vars %||% fit$vars, id = fit$id,
         single_level = isTRUE(fit$single_level),
         categorical = nominal(fit),
         ordinal = sort(fit$ordinal %||% character()),
         ordinal_levels = ordinal_levels(fit),
         count = sort(fit[["count"]] %||% character()),
         count_model = fit$extra_data$count_model %||% "poisson",
         missing = fit$missing %||% "error", centering = fit$centering %||% "none",
         scaling = fit$scaling %||% "none", min_variance = fit$min_variance,
         min_probability = fit$min_probability, drop_empty = fit$drop_empty %||% TRUE)
  }
  if (!identical(setting(null), setting(alternative))) {
    return(paste(
      "they differ in their indicators, how the indicators are coded (including",
      "the order of ordinal categories), or an estimation setting (missing",
      "data, centering, scaling, bounds)"))
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
  # The options as the fitter stores them, least free first. A value outside
  # the vocabulary is refused by name rather than ranked as missing.
  ranks <- c(slopes = .latents_lrt_option_rank(null$profile_slopes %||% "shared",
                                               alternative$profile_slopes %||% "shared",
                                               c("shared", "group_class")),
             dispersion = .latents_lrt_option_rank(
               null$extra_data$count_dispersion %||% "varying",
               alternative$extra_data$count_dispersion %||% "varying",
               c("equal", "varying")))
  if (is.na(ranks[["slopes"]])) {
    return("a fit records a `profile_slopes` value other than \"shared\" or \"group_class\"")
  }
  if (is.na(ranks[["dispersion"]])) {
    return("a fit records a `count_dispersion` value other than \"equal\" or \"varying\"")
  }
  if (ranks[["slopes"]] > 0L) {
    return("the first lets profile slopes vary by group class and the second does not")
  }
  if (ranks[["dispersion"]] > 0L) {
    return("the first lets count dispersions vary by profile and the second does not")
  }
  .multilpa_lrt_predictors_nested(null, alternative)
}

#' Compare one option of two fits on its ordered vocabulary
#' @param null_value,alternative_value The two fits' stored values.
#' @param order The vocabulary, least free first.
#' @return The null's rank minus the alternative's (positive when the null is
#'   freer), or `NA_integer_` when either value is not one string of `order`.
#' @noRd
.latents_lrt_option_rank <- function(null_value, alternative_value, order) {
  stopifnot("`order` must be a character vector" = is.character(order))
  position <- function(value) {
    if (!is.character(value) || length(value) != 1L) NA_integer_ else match(value, order)
  }
  position(null_value) - position(alternative_value)
}

#' Does every column of one design lie in the column space of another
#'
#' Both designs have the same rows: the comparison has already established
#' that the two fits share their observations and groups. A column lies in
#' the space when its least-squares residual on the other design is zero up to
#' rounding, relative to the column's own length.
#' @param inner,outer Numeric matrices, or `NULL` for no design.
#' @return One logical value; `TRUE` when `inner` has no columns.
#' @noRd
.latents_lrt_spans <- function(inner, outer) {
  if (is.null(inner) || ncol(inner) == 0L) return(TRUE)
  if (is.null(outer) || ncol(outer) == 0L || nrow(inner) != nrow(outer) ||
      !all(is.finite(inner)) || !all(is.finite(outer))) {
    return(FALSE)
  }
  residual <- qr.resid(qr(outer), inner)
  length_of <- function(columns) sqrt(colSums(columns^2))
  # Relative to the column alone, so the verdict does not depend on the units
  # a predictor was recorded in; a zero column has a zero residual.
  all(length_of(residual) <= sqrt(.Machine$double.eps) * length_of(inner))
}

#' Why the predictors of one profile model do not restrict another's
#'
#' The membership models are multinomial logits on the stacked profile design
#' (group-class intercepts and profile covariates, once per group class) and
#' on the group design (intercept and group covariates), and a direct effect
#' shifts an indicator's mean by the same linear predictor in every profile.
#' Setting the second fit's added coefficients to zero gives the first fit
#' back only when each design of the first lies in the column space of the
#' second's matching design; a shared predictor with other values does not.
#' Names alone cannot show this.
#' @param null,alternative Fits of class `multilpa` or `multilpa_covariates`
#'   whose predictor names already nest.
#' @return `NULL` when nested, otherwise the reason as one phrase.
#' @noRd
.multilpa_lrt_predictors_nested <- function(null, alternative) {
  stacked <- function(fit) {
    if (is.null(fit$profile_design)) NULL else do.call(rbind, fit$profile_design)
  }
  # A fit without membership covariates has group-class intercepts only,
  # which every alternative's design carries.
  if (length(null$profile_covariates) > 0L &&
      !.latents_lrt_spans(stacked(null), stacked(alternative))) {
    return(paste("a profile membership covariate of the first has values the",
                 "second's covariates do not reproduce"))
  }
  if (length(null$group_covariates) > 0L &&
      !.latents_lrt_spans(null$group_design, alternative$group_design)) {
    return(paste("a group-class membership covariate of the first has values the",
                 "second's covariates do not reproduce"))
  }
  direct <- function(fit, indicator) {
    carrier <- fit$extra_data$direct
    rows <- carrier$effects$indicator == indicator
    carrier$covariates[, carrier$effects$covariate_column[rows], drop = FALSE]
  }
  targets <- unique(null$extra_data$direct$effects$indicator %||% character())
  # The profile means absorb a constant, so each indicator's predictors are
  # judged together with an intercept.
  reproduced <- vapply(targets, function(indicator) {
    outer <- direct(alternative, indicator)
    .latents_lrt_spans(direct(null, indicator), cbind(1, outer))
  }, logical(1))
  if (!all(reproduced)) {
    return(paste("a direct-effect covariate of the first has values the second's",
                 "covariates do not reproduce"))
  }
  NULL
}

#' Do the covariate designs of one transition model restrict another's
#'
#' The initial-state logits are fitted on one design and the transition
#' logits on the designs of each move, stacked; with second-order
#' transitions, the move into occasion 2 and the later moves have
#' coefficients of their own, so those two blocks are judged apart. Each block
#' of the first must lie in the column space of the second's.
#' @param null,alternative Transition fits already nested by `.lta_nested()`.
#' @return One logical value; `TRUE` when either fit has no stored designs
#'   (a homogeneous fit without covariates, whose intercepts every
#'   alternative carries).
#' @noRd
.lta_lrt_designs_nested <- function(null, alternative) {
  if (is.null(null$designs) || is.null(alternative$designs)) return(TRUE)
  second_order <- identical(as.integer(.lta_arguments(alternative)$order), 2L)
  blocks <- function(designs) {
    moves <- designs$transition
    parts <- if (second_order) list(moves[1L], moves[-1L]) else list(moves)
    c(list(designs$initial), lapply(parts, function(part) do.call(rbind, part)))
  }
  inner <- blocks(null$designs)
  outer <- blocks(alternative$designs)
  all(vapply(seq_along(inner), function(block) {
    .latents_lrt_spans(inner[[block]], outer[[block]])
  }, logical(1)))
}

#' A fit's Gaussian covariance structure, as the comparison names it
#' @param fit A profile, class or transition fit.
#' @return One string: one of the fourteen structures, the label of a pattern
#'   of chosen covariances, or `NA` without continuous indicators.
#' @noRd
.latents_lrt_structure <- function(fit) {
  if (inherits(fit, c("multilpa_transitions", "multilpa_lta"))) {
    arguments <- .lta_arguments(fit)
    return(arguments$model %||% .multilpa_resolve_structure(
      arguments$variance_model, arguments$covariance_model %||% "diagonal"))
  }
  if (length(.multilpa_continuous_names(fit)) == 0L) return(NA_character_)
  .latents_comparison_structure(fit)
}

#' Is a nested covariance restriction singular at the null
#'
#' The covariance of profile k is `lambda_k D_k A_k t(D_k)` (volume, shape,
#' orientation). EVE and VVE share one orientation across profiles with
#' shapes that vary; EEV and VEV share one shape with orientations that vary.
#' At a spherical null (EII or VII), every `A_k` is the identity, so those
#' orientations have no value: near the null these alternatives are a union
#' of linear spaces, one per orientation, not one smooth extension, and
#' Wilks' chi-square reference does not hold. The review study
#' of EII within EVE rejected 8.0% of 500 null data sets at the nominal 5%
#' (validation/review-2026-11-01/lrt-orientation.log). Every other nested
#' pair of the fourteen structures is regular at a generic null: the full
#' structures (EEE, VEE, EVV, VVV) and the diagonal ones are smooth sets, and
#' a non-spherical null has distinct variances, which fix the orientation.
#' With one profile, or one continuous indicator, there is no orientation to
#' lose.
#' @param null,alternative Fits already found nested.
#' @return One logical value.
#' @noRd
.latents_lrt_singular_covariance <- function(null, alternative) {
  counted <- function(fit) {
    if (inherits(fit, c("multilpa_transitions", "multilpa_lta"))) {
      arguments <- .lta_arguments(fit)
      continuous <- fit$continuous %||%
        setdiff(arguments$vars, c(arguments$categorical, arguments$ordinal, arguments$count))
      c(profiles = arguments$n_profiles, continuous = length(continuous))
    } else c(profiles = fit$n_profiles, continuous = length(.multilpa_continuous_names(fit)))
  }
  size <- counted(alternative)
  isTRUE(.latents_lrt_structure(null) %in% c("EII", "VII") &&
           .latents_lrt_structure(alternative) %in% c("EVE", "VVE", "EEV", "VEV") &&
           size[["profiles"]] >= 2L && size[["continuous"]] >= 2L)
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
