# parameter_inference() for growth mixture models: one tidy table of the
# natural parameters, by Wald (delta method from the stored covariance) or by
# a person bootstrap whose replicates are matched to the fitted classes.

#' The natural parameters of a growth fit, each as a function of the parameters
#'
#' Coefficients and membership logits are coordinates of the estimation
#' vector; residual standard deviations, random-effect variances, standard
#' deviations and correlations, and class shares are transformations of it.
#' Each row carries the scale its Wald interval is formed on, and whether a
#' Wald test of zero is meaningful (not for a variance, whose zero is a bound).
#'
#' @param x A `latents_growth_mixture` fit.
#' @return A list with `labels` (a data frame: `block`, `class`, `term`),
#'   `scales`, `tested` and `values`, a function of a parameter list and its
#'   specification returning every natural parameter in row order.
#' @noRd
.growth_natural_parameters <- function(x) {
  spec <- x$spec
  names_theta <- .growth_names(spec)
  classes <- .growth_class_names(x)
  n_classes <- spec$n_classes
  terms <- colnames(spec$random_design)
  n_random <- length(terms)
  split_name <- function(names, prefix) {
    parts <- strsplit(sub(prefix, "", names), ".", fixed = TRUE)
    data.frame(class = vapply(parts, `[`, character(1), 1L),
               term = vapply(parts, function(p) paste(p[-1L], collapse = "."),
                             character(1)), stringsAsFactors = FALSE)
  }
  coefficient_names <- names_theta[startsWith(names_theta, "coefficient.")]
  membership_names <- names_theta[startsWith(names_theta, "membership.") |
                                    startsWith(names_theta, "group_membership.")]
  sd_classes <- if (identical(spec$variance, "equal")) 1L else seq_len(n_classes)
  random_classes <- if (identical(spec$random_covariance, "equal")) 1L else
    seq_len(n_classes)
  pairs <- if (isTRUE(spec$random_diagonal)) matrix(integer(), 0L, 2L) else
    which(upper.tri(diag(n_random)), arr.ind = TRUE)
  class_label <- function(k, shared) if (shared) "all" else classes[k]
  random_labels <- do.call(rbind, lapply(random_classes, function(k) {
    shared <- identical(spec$random_covariance, "equal")
    rbind(
      data.frame(block = "random_variance", class = class_label(k, shared),
                 term = terms, stringsAsFactors = FALSE),
      data.frame(block = "random_sd", class = class_label(k, shared),
                 term = terms, stringsAsFactors = FALSE),
      data.frame(block = rep("random_correlation", nrow(pairs)),
                 class = rep(class_label(k, shared), nrow(pairs)),
                 term = paste(terms[pairs[, 1L]], terms[pairs[, 2L]], sep = ":"),
                 stringsAsFactors = FALSE))
  }))
  labels <- rbind(
    cbind(block = "coefficient", split_name(coefficient_names, "^coefficient\\.")),
    if (length(membership_names) > 0L)
      cbind(block = ifelse(startsWith(membership_names, "group_"),
                           "group_membership", "membership"),
            split_name(membership_names, "^(group_)?membership\\.")),
    data.frame(block = "residual_sd",
               class = vapply(sd_classes, class_label, character(1),
                              shared = identical(spec$variance, "equal")),
               term = "residual", stringsAsFactors = FALSE),
    random_labels,
    data.frame(block = "share", class = classes, term = "share",
               stringsAsFactors = FALSE),
    if (.growth_is_correlated(spec))
      data.frame(block = "residual_process", class = "all",
                 term = names(x$params$correlation), stringsAsFactors = FALSE))
  rownames(labels) <- NULL
  scales <- c(coefficient = "identity", membership = "identity",
              group_membership = "identity", residual_sd = "log",
              random_variance = "log", random_sd = "log",
              random_correlation = "atanh", share = "logit",
              residual_process = "log")[labels$block]
  values <- function(params, spec) {
    theta <- .growth_pack(spec, params)
    covariances <- .growth_covariances(spec, params)
    random <- unlist(lapply(random_classes, function(k) {
      covariance <- covariances[[k]]
      variance <- diag(covariance)
      c(variance, sqrt(variance),
        covariance[pairs] / sqrt(variance[pairs[, 1L]] * variance[pairs[, 2L]]))
    }))
    weights <- spec$sampling_weights %||% rep(1, spec$n_groups)
    prior <- .growth_class_prior(spec, params)
    c(theta[coefficient_names], theta[membership_names],
      sqrt(params$sigma2[sd_classes]), random,
      colSums(prior * weights) / sum(weights),
      if (.growth_is_correlated(spec)) params$correlation) |> unname()
  }
  # A coordinate the correlated fit left on its floor or bound has no Wald or
  # bootstrap distribution worth reporting, as a degenerate covariance has not.
  on_bound <- x$residual_boundary %||% character()
  # Only the floored class's own residual-SD row (GR2-7).
  sigma_name <- ifelse(labels$class == "all", "log_sigma.shared",
                       paste0("log_sigma.", labels$class))
  residual_held <- labels$block == "residual_sd" & sigma_name %in% on_bound |
    labels$block == "residual_process" &
    sprintf("residual_correlation.log_%s", labels$term) %in% on_bound
  list(labels = labels, scales = unname(scales),
       tested = labels$block %in% c("coefficient", "membership", "group_membership"),
       degenerate = residual_held | labels$block %in% c("random_variance", "random_sd",
                                                        "random_correlation") &
         vapply(labels$class, function(label) {
           if (identical(label, "all")) length(x$degenerate_classes %||% character()) > 0L
           else .growth_class_is_degenerate(x, match(label, classes))
         }, logical(1)),
       values = values)
}

#' Reorder a growth fit's classes
#'
#' Class-specific coefficients, residual variances, random-effect covariances
#' (or proportional scales) and membership logits follow the classes; the
#' membership logits are rebased on the new first class.
#'
#' @param spec The growth specification.
#' @param params A parameter list.
#' @param order An integer vector: new class `k` is old class `order[k]`.
#' @return The reordered parameter list.
#' @noRd
.growth_permute_params <- function(spec, params, order) {
  params$beta <- params$beta[, order, drop = FALSE]
  params$sigma2 <- params$sigma2[order]
  if (identical(spec$random_covariance, "varying")) {
    params$random_covariance <- params$random_covariance[order]
  }
  if (identical(spec$random_covariance, "proportional")) {
    params$random_scale <- params$random_scale[order]
  }
  params$gamma <- params$gamma[, order, drop = FALSE] -
    params$gamma[, order[1L]]
  params
}

#' Match a refit's classes to the original fit's
#'
#' The refit classifies the original persons, and its classes are matched to
#' the original's by their posterior memberships of those persons, which
#' identifies a class by whom it contains rather than by its trajectory alone.
#'
#' @param x The original fit.
#' @param refit A fit to a resample.
#' @return The refit's parameters in the original class order.
#' @noRd
.growth_align_refit <- function(x, refit) {
  order <- .growth_match_classes(x, refit)$order
  .growth_permute_params(refit$spec, refit$params, order)
}

#' The order that matches a refit's classes to the original fit's
#' @param x The original fit.
#' @param refit Another fit of the same model.
#' @return A list with `order` (new class `k` is the refit's class
#'   `order[k]`) and `posterior`, the refit's posteriors of the original
#'   persons in its own labels.
#' @noRd
.growth_match_classes <- function(x, refit) {
  persons <- .growth_new_spec(refit, x$spec$model_data, TRUE)
  posterior <- .growth_fitted_expectation(refit, persons)$posterior
  list(order = .multilpa_match_order(t(x$expectation$posterior), t(posterior)),
       posterior = posterior)
}

#' The mixture_regression() arguments that refit a growth fit as it was fitted
#' @param x A `latents_growth_mixture` fit.
#' @param weights The weight column of the resampling data, or `NULL`.
#' @return A named list, without `data` and the start and convergence controls.
#' @noRd
.growth_refit_arguments <- function(x, weights) {
  spec <- x$spec
  list(formula = spec$formula, n_classes = spec$n_classes, family = "gaussian",
       id = spec$id, class_level = "group", common = spec$common,
       membership = spec$membership, variance = spec$variance,
       min_variance = x$settings$min_variance, random = spec$random,
       random_covariance = spec$random_covariance,
       random_diagonal = isTRUE(spec$random_diagonal),
       select_start = x$settings$select_start, vcov_type = "none", weights = weights,
       cluster = spec$cluster, n_group_classes = spec$n_cluster_classes %||% 1L,
       group_membership = spec$group_membership %||% ~1,
       residual_correlation = spec$residual_correlation %||% "none",
       residual_time = spec$residual_time)
}

#' Person bootstrap inference for a growth fit
#'
#' Persons are drawn with replacement (a person drawn twice is two persons),
#' each resample is refitted from fresh starts, and its classes are matched to
#' the original's before the natural parameters are read. A refit's start
#' failures, degenerate starts and boundary variances are expected and
#' muffled; one that does not converge is dropped with its reason.
#'
#' @return The bootstrap table and its attributes.
#' @noRd
.growth_bootstrap_inference <- function(x, natural, level, iter, n_starts, max_iter,
                                        tol, adjust) {
  if (!is.null(x$spec$cluster)) {
    stop(errorCondition(paste(
      "The bootstrap is not implemented for a multilevel growth mixture: its",
      "replicates would also need their group classes matched. Use",
      "`method = \"wald\"`."), class = "latents_unsupported_inference", call = NULL))
  }
  spec <- x$spec
  data <- spec$model_data
  weight_column <- NULL
  if (.latents_is_weighted(x)) {
    weight_column <- ".latents_bootstrap_weight"
    data[[weight_column]] <- spec$sampling_weights[spec$group_index]
  }
  arguments <- .growth_refit_arguments(x, weight_column)
  muffled <- c("latents_failed_starts", "latents_degenerate_start",
               "latents_random_boundary", "latents_unconverged")
  drawn <- .latents_cluster_bootstrap(
    data, spec$group_index, spec$n_groups, spec$id, iter,
    refit = function(resampled) {
      withCallingHandlers(
        do.call(mixture_regression, c(list(data = resampled), arguments,
                                      list(n_starts = n_starts, max_iter = max_iter,
                                           tol = tol))),
        warning = function(w) {
          if (inherits(w, muffled)) invokeRestart("muffleWarning")
        })
    },
    extract = function(refit) natural$values(.growth_align_refit(x, refit), refit$spec),
    reference_names = paste(natural$labels$block, natural$labels$class,
                            natural$labels$term, sep = "."))
  valid <- stats::complete.cases(drawn$estimates)
  if (sum(valid) < 2L) {
    reasons <- drawn$messages[!is.na(drawn$messages)]
    stop(errorCondition(sprintf(
      "Only %d of %d resamples produced a usable fit; the first reason was: %s",
      sum(valid), iter,
      if (length(reasons) > 0L) reasons[1L] else "no reason was recorded"),
      class = "latents_bootstrap_failed", call = NULL))
  }
  if (sum(valid) < iter) {
    warning(warningCondition(sprintf(
      "%d of %d resamples did not produce a usable fit and were dropped.",
      iter - sum(valid), iter), class = "latents_bootstrap_dropped", call = NULL))
  }
  kept <- drawn$estimates[valid, , drop = FALSE]
  bounds <- c((1 - level) / 2, (1 + level) / 2)
  quantiles <- apply(kept, 2L, stats::quantile, probs = bounds, names = FALSE)
  # A Wald statistic would put a normal approximation back on top of the
  # replicates the interval was read from: the interval is the inference.
  result <- data.frame(natural$labels,
                       estimate = natural$values(x$params, spec),
                       standard_error = apply(kept, 2L, stats::sd),
                       statistic = NA_real_, p_value = NA_real_, p_adjusted = NA_real_,
                       conf_low = quantiles[1L, ], conf_high = quantiles[2L, ],
                       row.names = NULL, stringsAsFactors = FALSE)
  covariance <- stats::cov(kept)
  structure(result, covariance = covariance, level = level, method = "bootstrap",
            iter = iter, n_valid = sum(valid), replicates = kept,
            messages = drawn$messages, adjust = adjust)
}

#' @details For a growth mixture ([mixture_regression()] with `random`), the
#'   table has one row per natural parameter: `block` (`"coefficient"`,
#'   `"membership"` and, for a multilevel growth mixture,
#'   `"group_membership"`, `"residual_sd"`, `"random_variance"`,
#'   `"random_sd"`, `"random_correlation"`, `"share"`), `class` (`"all"` for
#'   a quantity shared by every class), `term`, `estimate`, `standard_error`,
#'   `statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`. Wald
#'   rows use the fit's covariance (`vcov_type`) with the delta method;
#'   intervals of standard deviations and variances are formed on the log
#'   scale, of correlations on the Fisher scale and of shares on the logit
#'   scale, and only coefficients and logits are tested against zero.
#'   `method = "bootstrap"` resamples persons, refits each resample from
#'   `n_starts` starts and matches its classes to the original by their
#'   posterior memberships of the original persons; it is refused for a
#'   multilevel growth mixture (`latents_unsupported_inference`). In a
#'   simulation of 1000 datasets per design, Wald intervals covered
#'   0.936-0.965 with well-separated classes and 0.885-0.943 with weakly
#'   separated ones, whose standard errors were 6-30% too small
#'   (validation/GROWTH_SIMULATION.md in the source repository); the
#'   bootstrap is preferable for weakly separated classes. A class
#'   whose random-effect covariance is degenerate, or (with correlated
#'   residuals) a residual standard deviation or residual process held on its
#'   bound, raises `latents_boundary_fit` under `boundary = "error"`, naming
#'   which; under `"fix"` those rows have no standard error or interval. `data`, when
#'   given, must be the fitted data; `step` is not used (the information is
#'   the derivative of the analytic score).
#' @rdname parameter_inference
#' @export
parameter_inference.latents_growth_mixture <- function(
    x, data = NULL, level = 0.95, step = 1e-4,
    vcov_type = c("observed", "robust", "opg"),
    adjust = .multilpa_p_adjust_methods, method = c("wald", "bootstrap"),
    iter = 199L, n_starts = 10L, max_iter = 1000L, tol = 1e-8, seed = NULL,
    boundary = c("error", "fix")) {
  .latents_require(
    "`level` must be a single finite number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1,
    "`iter` must be a single whole number of at least 2" =
      is.numeric(iter) && length(iter) == 1L && is.finite(iter) && iter >= 2 &&
      iter == floor(iter))
  adjust <- .latents_match_formal(adjust)
  method <- .latents_match_formal(method)
  boundary <- .latents_match_formal(boundary)
  .multilpa_check_seed(seed)
  if (!is.null(data)) .trajectory_check_data(x, data)
  if (!isTRUE(x$converged)) {
    stop(errorCondition("The fit did not converge; infer from it only once it has.",
                        class = "latents_no_converge", call = NULL))
  }
  natural <- .growth_natural_parameters(x)
  if (any(natural$degenerate) && identical(boundary, "error")) {
    # Name what is on the bound: a random-effect covariance, or (with
    # correlated residuals) a residual SD or the residual process.
    blocks <- unique(natural$labels$block[natural$degenerate])
    causes <- c(
      if (any(startsWith(blocks, "random_"))) "a random-effect covariance is degenerate",
      if ("residual_sd" %in% blocks) "a residual standard deviation is at its floor",
      if ("residual_process" %in% blocks) "the residual process is on its bound")
    stop(errorCondition(sprintf(paste(
      "In this fit %s, so those parameters sit on a bound and have no Wald or",
      "bootstrap distribution worth reporting. Use `boundary = \"fix\"` to",
      "report the other parameters."), paste(causes, collapse = " and ")),
      class = "latents_boundary_fit", call = NULL))
  }
  if (identical(method, "bootstrap")) {
    result <- .mixture_with_seed(seed, .growth_bootstrap_inference(
      x, natural, level, as.integer(iter), n_starts, max_iter, tol, adjust))
  } else {
    vcov_type <- if (missing(vcov_type)) NULL else .latents_match_formal(vcov_type)
    inference <- .growth_resolve_inference(x, vcov_type)
    spec <- x$spec
    rows <- lapply(seq_len(nrow(natural$labels)), function(row) {
      .growth_delta_table(inference, function(v) {
        natural$values(.growth_unpack(spec, v, x$params), spec)[row]
      }, natural$scales[row], level)
    })
    wald <- do.call(rbind, rows)
    statistic <- ifelse(natural$tested & wald$std_error > 0,
                        wald$estimate / wald$std_error, NA_real_)
    result <- data.frame(natural$labels, estimate = wald$estimate,
                         standard_error = wald$std_error, statistic = statistic,
                         p_value = 2 * stats::pnorm(-abs(statistic)),
                         p_adjusted = NA_real_,
                         conf_low = wald$conf_low, conf_high = wald$conf_high,
                         row.names = NULL, stringsAsFactors = FALSE)
    attributes(result) <- c(attributes(result), list(
      covariance_unconstrained = inference$vcov, level = level, method = "wald",
      vcov_type = inference$vcov_type, adjust = adjust))
  }
  result[natural$degenerate, c("standard_error", "statistic", "p_value",
                               "conf_low", "conf_high")] <- NA_real_
  result$p_adjusted <- stats::p.adjust(result$p_value, method = adjust)
  result
}
