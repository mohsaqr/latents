# Tables, S3 methods and plots for multilpa(family = "additive") fits.

utils::globalVariables(c("group_class_label", "estimate", "conf_low", "parameter",
                         "conf_high", "component", "intercept_mean",
                         "offset_x", "panel"))

.additive_tables <- function() {
  c("parameters", "group_classes", "groups", "intercepts",
    "conditional_intercepts", "assignments", "classification", "fit", "starts")
}

#' Parameter table with standard errors, or with a classed notice and `NA`
#' standard errors when inference is refused for this fit
#' @noRd
.additive_parameters_or_notice <- function(x, level, vcov_type, adjust) {
  refusals <- c("latents_boundary_fit", "latents_no_converge",
                "latents_singular_information", "latents_too_few_groups")
  inference <- tryCatch(.additive_inference(x, vcov_type),
    error = function(condition) {
      if (!inherits(condition, refusals)) stop(condition)
      .multilpa_notice(sprintf("No standard errors in this table: %s",
                               conditionMessage(condition)),
                       class = "latents_no_standard_errors")
      NULL
    })
  .additive_parameter_table(x, inference, level, adjust)
}

#' The covariance a fit's tables use by default: observed, or the sandwich
#' for a weighted (pseudo likelihood) fit, whose observed errors are refused
#' @noRd
.additive_default_vcov <- function(x) {
  .latents_weighted_vcov(.latents_is_weighted(x), "observed", TRUE)
}

#' Effective groups per class, or `NA` where inference does not apply
#'
#' Computed without warning: the parameter table carries the warning, and a
#' refused fit has no intervals whose calibration could be in question. The
#' covariance is the fit's default, so the count agrees with the
#' `latents_weak_class` warning of the default parameter table.
#' @noRd
.additive_effective_or_na <- function(x) {
  refusals <- c("latents_boundary_fit", "latents_no_converge",
                "latents_singular_information", "latents_too_few_groups")
  inference <- tryCatch(
    withCallingHandlers(.additive_inference(x, .additive_default_vcov(x),
                                            warn_weak = FALSE),
                        latents_unconverged = function(w) {
                          invokeRestart("muffleWarning")
                        }),
    error = function(condition) {
      if (!inherits(condition, refusals)) stop(condition)
      NULL
    })
  if (is.null(inference)) return(rep(NA_real_, length(x$group_probabilities)))
  .additive_effective_groups(x, inference)$effective_groups
}

#' Relative entropy of the group classification (1 = perfectly separated)
#' @noRd
.additive_entropy <- function(posterior, weights = NULL) {
  .multilpa_relative_entropy(posterior, weights)
}

#' Cross-tabulate modal group classes against a known group-level variable
#' @noRd
.additive_recovery <- function(x, data, truth) {
  if (is.null(data) || is.null(truth)) {
    stop(errorCondition("The `recovery` table needs `data` and `truth`.",
                        class = "latents_bad_argument", call = NULL))
  }
  if (!is.character(truth) || length(truth) != 1L || !truth %in% names(data)) {
    stop(errorCondition("`truth` must name one column of `data`.",
                        class = "latents_bad_argument", call = NULL))
  }
  supplied <- .additive_prepare(data, x$vars, x$id)
  if (!identical(supplied$ids, x$sufficient_statistics$ids) ||
      !identical(supplied$sizes, x$sufficient_statistics$sizes)) {
    stop(errorCondition("`data` does not have the groups this fit was made from.",
                        class = "latents_bad_inference_data", call = NULL))
  }
  values <- data[[truth]]
  per_group <- tapply(as.character(values), supplied$index,
                      function(v) length(unique(v)))
  if (anyNA(values) || any(per_group != 1L)) {
    stop(errorCondition(sprintf(
      "`%s` must be complete and constant within each group: group classes are group-level.",
      truth), class = "latents_bad_data", call = NULL))
  }
  group_truth <- as.character(values)[match(seq_len(supplied$n_groups),
                                            supplied$index)]
  # Groups of weight zero are not in the sample (ADD2-2).
  used <- .additive_used_groups(x)
  assigned <- names(x$group_probabilities)[x$group_assignments][used]
  group_truth <- group_truth[used]
  counts <- as.data.frame(table(assigned = assigned, truth = group_truth),
                          stringsAsFactors = FALSE, responseName = "n")
  counts$share <- counts$n / as.vector(table(assigned)[counts$assigned])
  counts <- counts[order(counts$assigned, counts$truth), , drop = FALSE]
  rownames(counts) <- NULL
  counts
}

#' Build one named table of an additive fit
#' @noRd
.additive_table <- function(x, what, level = 0.95, vcov_type = "observed",
                            adjust = "none", data = NULL, truth = NULL) {
  class_names <- names(x$group_probabilities)
  posterior <- x$group_posteriors
  ids <- rownames(posterior)
  assigned <- x$group_assignments
  modal_posterior <- posterior[cbind(seq_along(assigned), assigned)]
  # Groups of weight zero still have a posterior and an assignment (the
  # per-group tables keep them), but they are in no count, entropy or
  # classification table (review 2026-10-26 round 2, ADD2-2).
  used <- .additive_used_groups(x)
  used_posterior <- posterior[used, , drop = FALSE]
  used_assigned <- assigned[used]
  used_weights <- (x$sampling_weights %||% rep(1, nrow(posterior)))[used]
  id_frame <- stats::setNames(data.frame(ids, stringsAsFactors = FALSE), x$id)
  switch(what,
    parameters = .additive_parameters_or_notice(x, level, vcov_type, adjust),
    group_classes = data.frame(
      group_class = class_names,
      weight = unname(x$group_probabilities),
      # Effective counts are weighted, like the class weights they sum to.
      count = unname(colSums(posterior * (x$sampling_weights %||% 1))),
      n_assigned = tabulate(used_assigned, length(class_names)),
      mean_posterior = vapply(seq_along(class_names), function(h) {
        if (any(used_assigned == h)) {
          stats::weighted.mean(used_posterior[used_assigned == h, h],
                                 used_weights[used_assigned == h])
        } else NA_real_
      }, numeric(1)),
      effective_groups = .additive_effective_or_na(x),
      stringsAsFactors = FALSE),
    groups = cbind(id_frame, data.frame(
      n = unname(x$group_sizes), group_class = class_names[assigned],
      posterior = modal_posterior,
      stats::setNames(as.data.frame(unname(posterior)),
                      paste0("probability_", class_names)),
      stringsAsFactors = FALSE)),
    intercepts = {
      stats <- x$sufficient_statistics
      cells <- expand.grid(group = seq_along(ids), indicator = seq_along(x$vars))
      cbind(stats::setNames(data.frame(ids[cells$group],
                                       stringsAsFactors = FALSE), x$id),
            data.frame(indicator = x$vars[cells$indicator],
                       n = unname(x$group_sizes)[cells$group],
                       average = as.vector(stats$averages),
                       intercept_mean = as.vector(x$intercept_means),
                       intercept_variance = as.vector(x$intercept_variances),
                       stringsAsFactors = FALSE))
    },
    conditional_intercepts = {
      cells <- expand.grid(group = seq_along(ids), indicator = seq_along(x$vars),
                           group_class = seq_along(class_names))
      cbind(stats::setNames(data.frame(ids[cells$group],
                                       stringsAsFactors = FALSE), x$id),
            data.frame(indicator = x$vars[cells$indicator],
                       group_class = class_names[cells$group_class],
                       posterior = posterior[cbind(cells$group, cells$group_class)],
                       intercept_mean = as.vector(x$conditional_intercept_means),
                       intercept_variance =
                         as.vector(x$conditional_intercept_variances),
                       stringsAsFactors = FALSE))
    },
    assignments = {
      index <- x$sufficient_statistics$index
      cbind(data.frame(row = seq_along(index)),
            stats::setNames(data.frame(ids[index], stringsAsFactors = FALSE),
                            x$id),
            data.frame(group_class = class_names[assigned[index]],
                       posterior = modal_posterior[index],
                       assignment_level = "group", stringsAsFactors = FALSE))
    },
    classification = {
      cells <- expand.grid(assigned = seq_along(class_names),
                           group_class = seq_along(class_names))
      data.frame(
        assigned = class_names[cells$assigned],
        group_class = class_names[cells$group_class],
        mean_posterior = mapply(function(a, h) {
          if (any(used_assigned == a)) {
            stats::weighted.mean(used_posterior[used_assigned == a, h],
                                 used_weights[used_assigned == a])
          } else NA_real_
        }, cells$assigned, cells$group_class),
        n_assigned = tabulate(used_assigned, length(class_names))[cells$assigned],
        stringsAsFactors = FALSE)
    },
    fit = data.frame(
      family = x$family, between_variance = x$between_variance,
      n_group_classes = length(class_names), log_likelihood = x$log_likelihood,
      n_parameters = x$n_parameters, aic = x$aic, bic = x$bic,
      bic_individual = x$bic_individual,
      n_groups = x$n_used_groups %||% x$n_groups,
      n_obs = x$sufficient_statistics$n_used_obs %||% x$n_obs,
      entropy = .additive_entropy(used_posterior, used_weights),
      smallest_class = min(x$group_probabilities),
      min_effective_groups = min(.additive_effective_or_na(x)),
      converged = x$converged,
      boundary = any(x$between_zero) || any(x$within_floor),
      kkt = if (is.na(x$kkt)) NA else isTRUE(x$kkt),
      n_starts = x$n_starts %||% nrow(x$starts), n_best_replicated = x$n_best_replicated,
      stringsAsFactors = FALSE),
    starts = x$starts,
    recovery = .additive_recovery(x, data, truth))
}

#' @noRd
#' @export
get_results.multilpa_additive <- function(x, what = "parameters", level = 0.95,
                                          vcov_type = c("observed", "robust", "opg"),
                                          adjust = .multilpa_p_adjust_methods,
                                          data = NULL, truth = NULL, ...) {
  .multilpa_reject_extra_arguments(list(...), "get_results()", paste(
    "An additive fit's tables take `level`, `vcov_type` and `adjust`, and",
    "`what = \"recovery\"` also `data` and `truth`."))
  what <- .multilpa_match_table(what, c(.additive_tables(), "recovery"))
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), .latents_match_formal(vcov_type),
                                      missing(vcov_type))
  adjust <- .latents_match_formal(adjust)
  .mixture_require_level(level)
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.additive_tables(), function(name) {
      .additive_table(x, name, level, vcov_type, adjust)
    }), .additive_tables()))
  }
  .additive_table(x, what, level, vcov_type, adjust, data, truth)
}

#' @noRd
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.multilpa_additive <- function(x, row.names = NULL, optional = FALSE,
                                            what = "parameters", ...) {
  get_results.multilpa_additive(x, what = what, ...)
}

#' @noRd
#' @param object A `multilpa_additive` fit.
#' @param scale `"natural"` (means, variances, class weights) or
#'   `"unconstrained"` (means, log variances, class-1-baseline logits). Names of
#'   `coef()` and the rows of `vcov()` agree on either scale.
#' @export
coef.multilpa_additive <- function(object, scale = c("natural", "unconstrained"),
                                   ...) {
  scale <- .latents_match_formal(scale)
  class_names <- names(object$group_probabilities)
  theta <- .additive_pack(
    list(means = unname(object$means), between = unname(object$between_variances),
         within = unname(object$within_variances),
         weights = unname(object$group_probabilities)),
    object$structure, class_names, object$vars)
  if (identical(scale, "unconstrained")) return(theta)
  .additive_natural(theta, length(class_names), object$vars,
                    object$structure, class_names)$estimate
}

#' @noRd
#' @param type Covariance type, as `vcov_type`.
#' @export
vcov.multilpa_additive <- function(object, type = c("observed", "robust", "opg"),
                                   scale = c("natural", "unconstrained"), ...) {
  # `vcov_type = "robust"`, the name parameter_inference() uses, used to fall
  # into `...` and return the observed-information covariance without a word
  # (INF-3).
  .latents_warn_unused_arguments(list(...), "vcov()",
                                 "Choose the covariance with `type =`.")
  scale <- .latents_match_arg(scale, c("natural", "unconstrained"), "scale")
  inference <- .additive_inference(object, .latents_weighted_vcov(
    .latents_is_weighted(object),
    .latents_match_arg(type, c("observed", "robust", "opg"), "type"), missing(type)))
  if (identical(scale, "unconstrained")) return(inference$vcov)
  class_names <- names(object$group_probabilities)
  jacobian <- .additive_natural(inference$theta, length(class_names), object$vars,
                                object$structure, class_names)$jacobian
  .inference_delta_covariance(jacobian, inference$vcov)
}

#' @noRd
#' @param parm Natural parameter names or positions; all by default.
#' @export
confint.multilpa_additive <- function(object, parm, level = 0.95,
                                      type = c("observed", "robust", "opg"), ...) {
  .multilpa_reject_extra_arguments(list(...), "confint()",
                                   "Choose the covariance with `type =`.")
  .mixture_require_level(level)
  table <- .additive_parameter_table(object, .additive_inference(object,
    .latents_weighted_vcov(.latents_is_weighted(object),
                           .latents_match_arg(type, c("observed", "robust", "opg"), "type"),
                           missing(type))),
                                     level)
  bounds <- cbind(table$conf_low, table$conf_high)
  percent <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                          trim = TRUE, scientific = FALSE, digits = 3), "%")
  dimnames(bounds) <- list(names(coef(object)), percent)
  # A name or position that is not an estimate is refused, as for the other
  # families, rather than returned as an NA row or a subscript error (INF-5).
  if (!missing(parm)) {
    bounds <- bounds[.mixture_confint_rows(parm, rownames(bounds)), , drop = FALSE]
  }
  bounds
}

#' @noRd
#' @export
logLik.multilpa_additive <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = nobs(object), class = "logLik")
}

#' @noRd
#' @export
# Groups of weight zero are not observations, as `nobs.glm` counts only
# non-zero weights; `n_used_groups` is absent from fits made before 0.97.2.
nobs.multilpa_additive <- function(object, ...) object$n_used_groups %||% object$n_groups

#' @rdname parameter_inference
#' @param ... Further arguments; unused by most methods.
#' @export
parameter_inference.multilpa_additive <- function(
    x, data = NULL, level = 0.95, step = 1e-4,
    vcov_type = c("observed", "robust", "opg"),
    adjust = .multilpa_p_adjust_methods, method = c("wald", "bootstrap"),
    ...) {
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), .latents_match_formal(vcov_type),
                                      missing(vcov_type))
  adjust <- .latents_match_formal(adjust)
  method <- .latents_match_formal(method)
  if (identical(method, "bootstrap")) {
    stop(errorCondition(
      "Bootstrap inference is not yet available for `family = \"additive\"`.",
      class = "latents_unsupported_inference", call = NULL))
  }
  if (!is.null(data)) {
    supplied <- .additive_prepare(data, x$vars, x$id, x$weights)
    stored <- x$sufficient_statistics
    # The sufficient statistics do not depend on the order of a group's rows,
    # but their floating-point sums do, by a few ulps; an exact comparison
    # refused the same data with rows reordered within groups (ADD-A2,
    # 2026-10-26). A changed or swapped value moves a mean by far more.
    if (!identical(supplied$ids, stored$ids) ||
        !identical(supplied$sizes, stored$sizes) ||
        !isTRUE(all.equal(supplied$averages, stored$averages,
                          tolerance = .additive_identity_tolerance)) ||
        !isTRUE(all.equal(supplied$scatter, stored$scatter,
                          tolerance = .additive_identity_tolerance)) ||
        !isTRUE(all.equal(supplied$sampling_weights, stored$sampling_weights,
                          tolerance = 1e-12))) {
      stop(errorCondition("`data` does not reproduce the data this fit was made from.",
                          class = "latents_bad_inference_data", call = NULL))
    }
  }
  inference <- .additive_inference(x, vcov_type, step)
  table <- .additive_parameter_table(x, inference, level, adjust)
  class_names <- names(x$group_probabilities)
  jacobian <- .additive_natural(inference$theta, length(class_names), x$vars,
                                x$structure, class_names)$jacobian
  structure(table,
            covariance = .inference_delta_covariance(jacobian, inference$vcov),
            covariance_unconstrained = inference$vcov,
            hessian = -inference$information, gradient = inference$gradient,
            scaled_score = inference$scaled_score,
            condition_ratio = inference$condition_ratio, level = level,
            step = step, method = "wald", vcov_type = vcov_type,
            adjust = adjust,
            group_scores = if (vcov_type == "observed") NULL else
              inference$group_scores)
}

#' @noRd
#' @export
summary.multilpa_additive <- function(object, level = 0.95,
                                      vcov_type = c("observed", "robust", "opg"),
                                      ...) {
  .multilpa_reject_extra_arguments(list(...), "summary()",
                                   "It takes `level` and `vcov_type`.")
  .mixture_require_level(level)
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(object),
                                      .latents_match_formal(vcov_type), missing(vcov_type))
  tables <- c("fit", "parameters", "group_classes")
  structure(stats::setNames(lapply(tables, function(name) {
    .additive_table(object, name, level, vcov_type)
  }), tables), class = "summary_multilpa_additive")
}

#' @noRd
#' @export
print.summary_multilpa_additive <- function(x, digits = 4L, ...) {
  headings <- c(fit = "Model fit", parameters = "Parameters",
                group_classes = "Group classes")
  invisible(lapply(names(x), function(name) {
    cat(headings[[name]], "\n", sep = "")
    print(x[[name]], digits = digits, row.names = FALSE)
    cat("\n")
  }))
  cat("Every table: get_results(x, what = ), e.g. \"groups\", \"intercepts\", \"starts\".\n")
  invisible(x)
}

#' @noRd
#' @export
as.data.frame.summary_multilpa_additive <- function(x, row.names = NULL,
                                                    optional = FALSE, ...) {
  .mixture_summary_results(x, "parameters", list(...), "as.data.frame()")
}

#' @noRd
#' @export
get_results.summary_multilpa_additive <- function(x, what = "parameters", ...) {
  .mixture_summary_results(x, what, list(...))
}

#' Class key: label with share, colour and shape, largest class first
#' @noRd
.additive_key <- function(x) {
  class_names <- names(x$group_probabilities)
  share <- unname(x$group_probabilities)
  label <- sprintf("Group class %d (%s)", seq_along(class_names),
                   .gg_percent(share))
  # A block shared by every class is drawn once, in grey, with its own shape;
  # the classes never take that grey or that shape.
  label <- c(label, "All group classes (shared)")
  data.frame(group_class = c(class_names, "shared"),
             label = factor(label, levels = label),
             colour = c(rep_len(.gg_class_palette, length(class_names)),
                        .gg_reserved_grey),
             shape = c(rep_len(setdiff(.gg_shapes, 23L), length(class_names)), 23L),
             stringsAsFactors = FALSE)
}

#' @rdname latents-plots
#' @export
plot.multilpa_additive <- function(x, what = c("means", "variances", "intercepts"),
                                   level = 0.95, main = NULL, subtitle = NULL,
                                   ...) {
  .latents_warn_unused_arguments(
    list(...), "plot()",
    "Style the returned plot with ggplot2, for example `+ ggplot2::theme()`.")
  what <- .latents_match_formal(what)
  .gg_require()
  key <- .additive_key(x)
  # The refusal notice is replaced by the caption below, so only that one
  # classed message is muffled.
  table <- withCallingHandlers(
    .additive_table(x, "parameters", level, .additive_default_vcov(x)),
    latents_no_standard_errors = function(condition) {
      invokeRestart("muffleMessage")
    })
  has_errors <- any(is.finite(table$standard_error))
  caption <- if (has_errors) {
    sprintf("Bars: %.0f%% Wald intervals.", 100 * level)
  } else "No intervals: inference is not available for this fit."
  switch(what,
    means = {
      frame <- subset(table, parameter == "mean")
      frame$label <- key$label[match(frame$group_class, key$group_class)]
      frame$label <- droplevels(frame$label)
      key <- subset(key, label %in% levels(frame$label))
      key$label <- droplevels(key$label)
      frame$indicator <- factor(frame$indicator, levels = x$vars)
      plot <- ggplot2::ggplot(frame, ggplot2::aes(
        x = indicator, y = estimate, colour = label, fill = label, shape = label,
        group = label)) +
        ggplot2::geom_line(linewidth = 0.6, alpha = 0.7) +
        ggplot2::geom_point(size = 3, colour = "grey20")
      if (has_errors) {
        plot <- plot + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.12)
      }
      plot + .gg_scales(data.frame(label = key$label, colour = key$colour,
                                   shape = key$shape), legend = TRUE) +
        .gg_theme() + .gg_top_legend() +
        ggplot2::labs(x = NULL, y = "Group-class mean") +
        .gg_titles(main, subtitle, "Group-class means",
                   "Means of the group intercepts in each class", caption)
    },
    variances = {
      # One panel per level; within a panel each class (or the shared block)
      # has the class key's colour and shape, as in the means view.
      frame <- subset(table, parameter == "variance")
      frame$label <- key$label[match(frame$group_class, key$group_class)]
      frame$label <- droplevels(frame$label)
      key <- subset(key, label %in% levels(frame$label))
      key$label <- droplevels(key$label)
      frame$panel <- factor(ifelse(frame$level == "within", "Within groups",
                                   "Between groups"),
                            levels = c("Within groups", "Between groups"))
      frame$indicator <- factor(frame$indicator, levels = x$vars)
      dodge <- ggplot2::position_dodge(width = 0.5)
      plot <- ggplot2::ggplot(frame, ggplot2::aes(
        x = indicator, y = estimate, colour = label, fill = label,
        shape = label, group = label)) +
        ggplot2::geom_point(size = 3, colour = "grey20", position = dodge)
      if (has_errors) {
        plot <- plot + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15,
          position = dodge)
      }
      # Points are filled and shaped by class; colour maps only the interval
      # bars, so its scale is added only when there are bars to colour.
      plot + ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free_y") +
        .gg_scales(data.frame(label = key$label, colour = key$colour,
                              shape = key$shape),
                   aesthetics = c(if (has_errors) "colour", "fill", "shape"),
                   legend = TRUE) +
        ggplot2::expand_limits(y = 0) +
        .gg_theme() + .gg_top_legend() +
        ggplot2::labs(x = NULL, y = "Variance") +
        .gg_titles(main, subtitle, "Within- and between-group variances",
                   sprintf("%s family", .additive_family_label(x$family)),
                   caption)
    },
    intercepts = {
      frame <- .additive_table(x, "intercepts")
      groups <- .additive_table(x, "groups")
      frame$group_class <- groups$group_class[match(frame[[x$id]], groups[[x$id]])]
      frame$label <- key$label[match(frame$group_class, key$group_class)]
      frame$indicator <- factor(frame$indicator, levels = x$vars)
      # A deterministic horizontal offset per class keeps the classes apart
      # without jitter, so the figure is reproducible.
      n_classes <- length(x$group_probabilities)
      key <- subset(key, group_class != "shared")
      key$label <- droplevels(key$label)
      frame$offset_x <- as.integer(frame$indicator) +
        (match(frame$group_class, key$group_class) - (n_classes + 1) / 2) * 0.18
      means <- subset(table, parameter == "mean")
      if (identical(x$structure$means, "equal")) {
        # A shared mean is marked at every class's position.
        means <- do.call(rbind, lapply(names(x$group_probabilities), function(h) {
          transform(means, group_class = h)
        }))
      }
      means$label <- key$label[match(means$group_class, key$group_class)]
      means$offset_x <- match(means$indicator, x$vars) +
        (match(means$group_class, key$group_class) - (n_classes + 1) / 2) * 0.18
      ggplot2::ggplot(frame, ggplot2::aes(x = offset_x, y = intercept_mean,
                                          colour = label, fill = label,
                                          shape = label)) +
        ggplot2::geom_point(size = 1.6, alpha = 0.55, stroke = 0.3) +
        ggplot2::geom_point(data = means, ggplot2::aes(y = estimate), size = 4,
                            colour = "grey10", stroke = 0.8) +
        ggplot2::scale_x_continuous(breaks = seq_along(x$vars), labels = x$vars) +
        .gg_scales(data.frame(label = key$label, colour = key$colour,
                              shape = key$shape), legend = TRUE) +
        .gg_theme() + .gg_top_legend() +
        ggplot2::labs(x = NULL, y = "Posterior group intercept") +
        .gg_titles(main, subtitle, "Group intercepts by modal group class",
                   "Small marks: groups; large marks: class means", NULL)
    })
}
