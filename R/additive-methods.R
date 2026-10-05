# Tables, S3 methods and plots for multilpa(family = "additive") fits.

utils::globalVariables(c("group_class_label", "estimate", "conf_low", "parameter",
                         "conf_high", "component", "intercept_mean",
                         "offset_x"))

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

#' Effective groups per class, or `NA` where inference does not apply
#'
#' Computed without warning: the parameter table carries the warning, and a
#' refused fit has no intervals whose calibration could be in question.
#' @noRd
.additive_effective_or_na <- function(x) {
  refusals <- c("latents_boundary_fit", "latents_no_converge",
                "latents_singular_information", "latents_too_few_groups")
  inference <- tryCatch(
    withCallingHandlers(.additive_inference(x, "observed", warn_weak = FALSE),
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
.additive_entropy <- function(posterior) {
  if (ncol(posterior) < 2L) return(NA_real_)
  terms <- ifelse(posterior > 0, posterior * log(posterior), 0)
  1 + sum(terms) / (nrow(posterior) * log(ncol(posterior)))
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
  assigned <- names(x$group_probabilities)[x$group_assignments]
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
  id_frame <- stats::setNames(data.frame(ids, stringsAsFactors = FALSE), x$id)
  switch(what,
    parameters = .additive_parameters_or_notice(x, level, vcov_type, adjust),
    group_classes = data.frame(
      group_class = class_names,
      weight = unname(x$group_probabilities),
      count = unname(colSums(posterior)),
      n_assigned = tabulate(assigned, length(class_names)),
      mean_posterior = vapply(seq_along(class_names), function(h) {
        if (any(assigned == h)) mean(posterior[assigned == h, h]) else NA_real_
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
          if (any(assigned == a)) mean(posterior[assigned == a, h]) else NA_real_
        }, cells$assigned, cells$group_class),
        n_assigned = tabulate(assigned, length(class_names))[cells$assigned],
        stringsAsFactors = FALSE)
    },
    fit = data.frame(
      family = x$family, between_variance = x$between_variance,
      n_group_classes = length(class_names), log_likelihood = x$log_likelihood,
      n_parameters = x$n_parameters, aic = x$aic, bic = x$bic,
      bic_individual = x$bic_individual, n_groups = x$n_groups,
      n_obs = x$n_obs, entropy = .additive_entropy(posterior),
      smallest_class = min(x$group_probabilities),
      min_effective_groups = min(.additive_effective_or_na(x)),
      converged = x$converged,
      boundary = any(x$between_zero) || any(x$within_floor),
      kkt = if (is.na(x$kkt)) NA else isTRUE(x$kkt),
      n_starts = nrow(x$starts), n_best_replicated = x$n_best_replicated,
      stringsAsFactors = FALSE),
    starts = x$starts,
    recovery = .additive_recovery(x, data, truth))
}

#' Tables of an additive group-class fit
#'
#' Tidy tables of a fit from `multilpa(family = "additive")`. Every table is a
#' base `data.frame`; the fit itself never needs to be reached into.
#'
#' @param x A `multilpa_additive` fit.
#' @param what Which table; `"all"` returns a named list of every table.
#' @param level Confidence level of the intervals in `"parameters"`.
#' @param vcov_type `"observed"` (inverse observed information), `"robust"`
#'   (CR0 sandwich over groups) or `"opg"` (inverse outer product of the group
#'   scores).
#' @param adjust p-value adjustment across the class means, one of
#'   [stats::p.adjust.methods]; `"none"` by default, as elsewhere in latents.
#' @param data,truth For `what = "recovery"` only: the data the fit was made
#'   from and the name of a column, constant within groups, holding a known
#'   group classification.
#' @param ... Unused.
#' @section Tables:
#' \describe{
#'   \item{`parameters`}{One row per natural parameter: `level` (`"between"`
#'     for class means and between-group variances, `"within"` for the shared
#'     within-group variances, `"group_class"` for class weights),
#'     `group_class` (`"shared"` for a parameter common to all classes),
#'     `indicator` (`NA` for weights), `parameter` (`"mean"`, `"variance"`,
#'     `"weight"`), `estimate`, `standard_error`, `statistic`, `p_value`,
#'     `p_adjusted`, `conf_low`, `conf_high`. Only means carry a Wald test.
#'     Intervals for variances are formed on the log scale and for weights on
#'     the logit scale. When inference is refused (boundary, unconverged,
#'     singular information) the standard-error columns are `NA` and a
#'     `latents_no_standard_errors` message says why.}
#'   \item{`group_classes`}{One row per group class: `weight`, `count`
#'     (summed posterior), `n_assigned` (modal), `mean_posterior` among the
#'     groups assigned to it, and `effective_groups`: the number of perfectly
#'     classified groups that would estimate the class weight as precisely as
#'     this fit does, `w^2 (1 - w) / Var(w)`. It equals the class's expected
#'     group count when classes are perfectly separated and falls when a
#'     class is small or poorly separated. Below 50, intervals for class
#'     parameters are flagged with a `latents_weak_class` warning (threshold
#'     and its simulation evaluation: `validation/ADDITIVE_SIMULATION.md` in
#'     the source repository). `NA` for one class or when inference is refused.}
#'   \item{`groups`}{One row per group (named by the `id` column): size `n`,
#'     modal `group_class`, its `posterior`, and `probability_group_class_h`.}
#'   \item{`intercepts`}{One row per group and indicator: `n`, the observed
#'     group `average`, and the posterior mean and variance of the group's
#'     intercept, mixed over its class posterior (so the variance includes
#'     uncertainty about class membership). Conditional on the estimated
#'     parameters; not a statement about their sampling uncertainty.}
#'   \item{`conditional_intercepts`}{As `intercepts`, one row per group,
#'     indicator and group class, conditional on that class, with the class
#'     `posterior`.}
#'   \item{`assignments`}{One row per data row: `row`, the `id` column, the
#'     group's modal `group_class` and `posterior`, and
#'     `assignment_level = "group"`: every row inherits its group's class.}
#'   \item{`classification`}{One row per assigned and group class: the
#'     `mean_posterior` of `group_class` among groups assigned to `assigned`.}
#'   \item{`fit`}{One row: log likelihood, `n_parameters`, `aic`, `bic`
#'     (penalized by the number of groups), `bic_individual` (by rows),
#'     `n_groups`, `n_obs`, relative `entropy` over groups, the smallest class
#'     weight, `min_effective_groups`, `converged`, `boundary`, `kkt` and start
#'     replication.}
#'   \item{`recovery`}{Needs `data` and `truth`: one row per modal group
#'     class and known value, with the number of groups `n` and their `share`
#'     of the assigned class. Not part of `what = "all"`.}
#'   \item{`starts`}{One row per start: likelihood, convergence, iterations,
#'     number of between variances held at zero, error message of a failed
#'     start, and which start was selected.}
#' }
#' @return A base `data.frame`, or a named list of them for `what = "all"`.
#' @examples
#' set.seed(1)
#' ratings <- data.frame(
#'   team = rep(seq_len(40), each = 6),
#'   climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'     rnorm(240),
#'   support = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'     rnorm(240))
#' fit <- multilpa(ratings, c("climate", "support"), "team",
#'                 n_group_classes = 2, family = "additive", n_starts = 3,
#'                 seed = 1)
#' get_results(fit)
#' get_results(fit, "group_classes")
#' get_results(fit, "parameters", vcov_type = "robust")
#' @export
get_results.multilpa_additive <- function(x, what = "parameters", level = 0.95,
                                          vcov_type = c("observed", "robust", "opg"),
                                          adjust = .multilpa_p_adjust_methods,
                                          data = NULL, truth = NULL, ...) {
  what <- match.arg(what, c(.additive_tables(), "recovery", "all"))
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), match.arg(vcov_type),
                                      missing(vcov_type))
  adjust <- match.arg(adjust)
  stopifnot("`level` must be a single number in (0, 1)" =
              is.numeric(level) && length(level) == 1L && level > 0 &&
              level < 1)
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.additive_tables(), function(name) {
      .additive_table(x, name, level, vcov_type, adjust)
    }), .additive_tables()))
  }
  .additive_table(x, what, level, vcov_type, adjust, data, truth)
}

#' @rdname get_results.multilpa_additive
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.multilpa_additive <- function(x, row.names = NULL, optional = FALSE,
                                            what = "parameters", ...) {
  get_results.multilpa_additive(x, what = what, ...)
}

#' @rdname get_results.multilpa_additive
#' @param object A `multilpa_additive` fit.
#' @param scale `"natural"` (means, variances, class weights) or
#'   `"unconstrained"` (means, log variances, class-1-baseline logits). Names of
#'   `coef()` and the rows of `vcov()` agree on either scale.
#' @export
coef.multilpa_additive <- function(object, scale = c("natural", "unconstrained"),
                                   ...) {
  scale <- match.arg(scale)
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

#' @rdname get_results.multilpa_additive
#' @param type Covariance type, as `vcov_type`.
#' @export
vcov.multilpa_additive <- function(object, type = c("observed", "robust", "opg"),
                                   scale = c("natural", "unconstrained"), ...) {
  scale <- match.arg(scale)
  inference <- .additive_inference(object, .latents_weighted_vcov(
    .latents_is_weighted(object), match.arg(type), missing(type)))
  if (identical(scale, "unconstrained")) return(inference$vcov)
  class_names <- names(object$group_probabilities)
  jacobian <- .additive_natural(inference$theta, length(class_names), object$vars,
                                object$structure, class_names)$jacobian
  .inference_delta_covariance(jacobian, inference$vcov)
}

#' @rdname get_results.multilpa_additive
#' @param parm Natural parameter names or positions; all by default.
#' @export
confint.multilpa_additive <- function(object, parm, level = 0.95,
                                      type = c("observed", "robust", "opg"), ...) {
  table <- .additive_parameter_table(object, .additive_inference(object,
    .latents_weighted_vcov(.latents_is_weighted(object), match.arg(type), missing(type))),
                                     level)
  bounds <- cbind(table$conf_low, table$conf_high)
  percent <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                          trim = TRUE, scientific = FALSE, digits = 3), "%")
  dimnames(bounds) <- list(names(coef(object)), percent)
  if (!missing(parm)) bounds <- bounds[parm, , drop = FALSE]
  bounds
}

#' @rdname get_results.multilpa_additive
#' @export
logLik.multilpa_additive <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' @rdname get_results.multilpa_additive
#' @export
nobs.multilpa_additive <- function(object, ...) object$n_groups

#' @rdname parameter_inference
#' @param ... Unused by most methods; accepted so every method shares the
#'   generic's signature.
#' @export
parameter_inference.multilpa_additive <- function(
    x, data = NULL, level = 0.95, step = 1e-4,
    vcov_type = c("observed", "robust", "opg"),
    adjust = .multilpa_p_adjust_methods, method = c("wald", "bootstrap"),
    ...) {
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x), match.arg(vcov_type),
                                      missing(vcov_type))
  adjust <- match.arg(adjust)
  method <- match.arg(method)
  if (identical(method, "bootstrap")) {
    stop(errorCondition(
      "Bootstrap inference is not yet available for `family = \"additive\"`.",
      class = "latents_unsupported_inference", call = NULL))
  }
  if (!is.null(data)) {
    supplied <- .additive_prepare(data, x$vars, x$id)
    stored <- x$sufficient_statistics
    if (!identical(supplied$ids, stored$ids) ||
        !identical(supplied$sizes, stored$sizes) ||
        !isTRUE(all.equal(supplied$averages, stored$averages)) ||
        !isTRUE(all.equal(supplied$scatter, stored$scatter))) {
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

#' @rdname latents-summary
#' @export
summary.multilpa_additive <- function(object, level = 0.95,
                                      vcov_type = c("observed", "robust", "opg"),
                                      ...) {
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(object),
                                      match.arg(vcov_type), missing(vcov_type))
  tables <- c("fit", "parameters", "group_classes")
  structure(stats::setNames(lapply(tables, function(name) {
    .additive_table(object, name, level, vcov_type)
  }), tables), class = "summary_multilpa_additive")
}

#' @rdname latents-print
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

#' @rdname latents-as-data-frame
#' @export
as.data.frame.summary_multilpa_additive <- function(x, row.names = NULL,
                                                    optional = FALSE, ...) {
  x$parameters
}

#' Class key: label with share, colour and shape, largest class first
#' @noRd
.additive_key <- function(x) {
  class_names <- names(x$group_probabilities)
  share <- unname(x$group_probabilities)
  label <- sprintf("Group class %d (%s)", seq_along(class_names),
                   .gg_percent(share))
  # A block shared by every class is drawn once, in grey, with its own shape.
  label <- c(label, "All group classes (shared)")
  data.frame(group_class = c(class_names, "shared"),
             label = factor(label, levels = label),
             colour = c(rep_len(.gg_okabe_ito, length(class_names)), "#999999"),
             shape = c(rep_len(.gg_shapes, length(class_names)), 23L),
             stringsAsFactors = FALSE)
}

#' Plot an additive group-class fit
#'
#' `"means"` shows each group class's mean per indicator with its Wald
#' interval (when the fit supports inference). `"variances"` sets the shared
#' within-group variance beside each between-group variance. `"intercepts"`
#' shows every group's posterior intercept mean, marked by its modal class,
#' with the class means over them. Classes are distinguished by colour and by
#' point shape.
#'
#' @param x A `multilpa_additive` fit.
#' @param what Which view.
#' @param level Confidence level of the intervals.
#' @param main,subtitle Optional title and subtitle.
#' @param ... Unused.
#' @return A ggplot object. Raises `latents_missing_package` without ggplot2.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   ratings <- data.frame(
#'     team = rep(seq_len(40), each = 6),
#'     climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'       rnorm(240),
#'     support = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
#'       rnorm(240))
#'   fit <- multilpa(ratings, c("climate", "support"), "team",
#'                   n_group_classes = 2, family = "additive", n_starts = 3,
#'                   seed = 1)
#'   plot(fit)
#'   plot(fit, what = "intercepts")
#' }
#' @export
plot.multilpa_additive <- function(x, what = c("means", "variances", "intercepts"),
                                   level = 0.95, main = NULL, subtitle = NULL,
                                   ...) {
  what <- match.arg(what)
  .gg_require()
  key <- .additive_key(x)
  # The refusal notice is replaced by the caption below, so only that one
  # classed message is muffled.
  table <- withCallingHandlers(
    .additive_table(x, "parameters", level),
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
      frame <- subset(table, parameter == "variance")
      frame$component <- paste(
        ifelse(frame$level == "within", "Within groups", "Between groups"),
        ifelse(frame$group_class == "shared", "(shared)",
               paste0("(", sub("group_class_", "class ", frame$group_class), ")")))
      components <- unique(frame$component)
      frame$component <- factor(frame$component, levels = components)
      frame$indicator <- factor(frame$indicator, levels = x$vars)
      # Shared components grey; class-specific ones by the class palette.
      shared <- grepl("(shared)", components, fixed = TRUE)
      colours <- stats::setNames(ifelse(shared, "#999999", rep_len(
        .gg_okabe_ito, length(components))), components)
      shapes <- stats::setNames(rep_len(c(23L, .gg_shapes), length(components)),
                                components)
      dodge <- ggplot2::position_dodge(width = 0.5)
      plot <- ggplot2::ggplot(frame, ggplot2::aes(
        x = indicator, y = estimate, colour = component, fill = component,
        shape = component)) +
        ggplot2::geom_point(size = 3, colour = "grey20", position = dodge)
      if (has_errors) {
        plot <- plot + ggplot2::geom_errorbar(
          ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15,
          position = dodge)
      }
      plot + ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::scale_shape_manual(values = shapes) +
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
