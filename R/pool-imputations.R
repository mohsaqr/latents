utils::globalVariables(c("estimate", "label", "source", "conf_low", "conf_high"))

# Fitting a model to multiply imputed data and pooling it with Rubin's rules.
#
# Missing indicators are integrated out by `missing = "fiml"`; missing
# covariates cannot be, because the model conditions on them and has no
# distribution to integrate over. Multiple imputation is the standard answer:
# impute the covariates (with the indicators in the imputation model), fit each
# completed data set, and pool. Pooling a mixture has one step a regression
# does not: the imputations label their profiles and group classes
# independently, so they are aligned to the first before anything is averaged.

#' Fit a model to multiply imputed data and pool it with Rubin's rules
#'
#' Fits [multilpa()] to every completed data set, aligns each fit's profile and
#' group-class labels to those of the first, and pools the estimates with
#' Rubin's rules, returning a `latents_pooled` object. Use it for missing
#' covariates, which `missing = "fiml"` cannot handle because the model
#' conditions on covariates. Impute them first (for example with the mice
#' package, with the indicators and the group identifier in the imputation
#' model) and pass the completed data.
#'
#' Pooling is on the reported (natural) scale: the pooled estimate is the mean
#' over imputations, and its variance is the mean within-imputation variance
#' plus `(1 + 1/m)` times the between-imputation variance. Degrees of freedom
#' follow Rubin (1987) with a large complete-data sample. Pooled intervals are
#' symmetric t intervals on this scale, so for probabilities and variances
#' they can extend outside the parameter space, unlike the transformed Wald
#' intervals of [parameter_inference()]. A structural constant, such as the
#' probability of one for a single class, has zero within- and
#' between-imputation variance and zero missing information. When the
#' within-imputation variance is zero and the between-imputation variance is
#' positive, `riv` is infinite and `fmi` is one.
#'
#' @section Choosing the imputation model:
#' The pooled results are only as valid as the imputations. The covariate and
#' the indicators are related through the latent classes, which an imputation
#' model linear in the indicators (mice's defaults) does not capture. In
#' `validation/pool-imputations-coverage.R` (a membership slope with 30% of
#' its covariate missing at random given the indicators, 200 replications,
#' ten imputations), 95% intervals covered 0.925 under predictive mean
#' matching and 0.935 under Bayesian linear regression, and 0.995 when the
#' covariate was drawn approximately from its conditional distribution given
#' the known generating parameters. That reference arm omits uncertainty in
#' the imputation parameters, so it does not validate proper multiple
#' imputation, and the study covers one design only. The imputation model
#' should account for the latent structure and the clustering, include the
#' indicators and relevant nonlinearities, and propagate parameter
#' uncertainty; imputing separately within assigned classes ignores
#' classification uncertainty. Complete-case analysis covered 0.955 in that
#' design, but selecting on the indicators changes their within-class
#' distributions, and the case-control argument for unbiased logistic slopes
#' does not carry over when the outcome is latent and its measurement model
#' is estimated.
#'
#' @param imputed The completed data sets: a list of at least two data frames
#'   with the same columns and rows in the same order, or a `mids` object from
#'   the mice package, which is completed here.
#' @param vars,id,n_profiles As in [multilpa()].
#' @param ... Further arguments for [multilpa()], such as `n_group_classes`,
#'   `profile_covariates`, `group_covariates`, `categorical`, `n_starts` and
#'   `seed`, used for every imputation.
#' @param level Confidence level of the pooled intervals (default 0.95).
#' @param vcov_type `"observed"` (default), `"robust"` or `"opg"`, as in
#'   [parameter_inference()], applied to every imputation's fit. For fits
#'   with sampling weights the default is `"robust"`, and the other types are
#'   refused with `latents_unsupported_weights`.
#' @param adjust Multiplicity correction applied to the pooled p-values to
#'   give `p_adjusted`, as in [parameter_inference()]: `"none"` (default) or
#'   a method of [stats::p.adjust()]. The per-imputation tables are computed
#'   without correction.
#' @return An object of class `latents_pooled`. [as.data.frame()] and
#'   `get_results(x, "estimates")` give one row per parameter, with the columns
#'   of [parameter_inference()] (`level`, `outcome`, `term`, `parameter`,
#'   `estimate`, `standard_error`, `statistic`, `p_value`, `p_adjusted`,
#'   `conf_low`, `conf_high`) plus `df` (Rubin's degrees of freedom),
#'   `within` and `between` (the two variance components), `riv` (the
#'   relative increase in variance due to missingness) and `fmi` (the fraction
#'   of missing information). `statistic` and `p_value` are `NA` where
#'   [parameter_inference()] reports no test (variances and probabilities).
#'   `get_results(x, "imputations")` gives each imputation's aligned estimates,
#'   one row per imputation and parameter, and `get_results(x, "fits")` one row
#'   per imputation with its log likelihood, convergence and the label
#'   permutation that aligned it.
#' @section Conditions:
#'   `latents_bad_data` when `imputed` is not a list of at least two data frames
#'   with identical columns and row counts; `latents_pooling_failed` when any
#'   imputation's fit or inference fails (the message names the imputation and
#'   gives the original reason); `latents_missing_package` for a `mids` object
#'   without mice installed. Warnings from the individual fits are passed on
#'   as they are raised.
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
#' Wiley.
#'
#' Meng, X.-L. (1994). Multiple-imputation inferences with uncongenial
#' sources of input. *Statistical Science*, 9(4), 538--558.
#'
#' van Buuren, S. (2018). *Flexible Imputation of Missing Data* (2nd ed.),
#' section 5.2. Chapman & Hall/CRC.
#' @examples
#' set.seed(1)
#' # Two stand-in completed data sets. With real missing covariates, pass the
#' # `mids` object mice::mice() returns instead.
#' completed <- lapply(seq_len(2), function(i) {
#'   within(subset(course_engagement, student <= 40),
#'          previous_grade <- previous_grade + rnorm(length(previous_grade), sd = 0.1))
#' })
#' pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"),
#'                            "student", n_profiles = 2, n_group_classes = 1,
#'                            profile_covariates = "previous_grade",
#'                            n_starts = 2, seed = 1)
#' pooled
#' get_results(pooled, "fits")
#' @export
pool_imputations <- function(imputed, vars, id, n_profiles, ..., level = 0.95,
                             vcov_type = c("observed", "robust", "opg"),
                             adjust = .multilpa_p_adjust_methods) {
  stopifnot(
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1)
  vcov_defaulted <- missing(vcov_type)
  vcov_type <- match.arg(vcov_type)
  adjust <- match.arg(adjust)
  completed <- .multilpa_completed_data(imputed)
  call <- match.call()
  fits <- lapply(seq_along(completed), function(index) {
    .multilpa_pool_step(index, "fitting", {
      multilpa(completed[[index]], vars, id, n_profiles, ...)
    })
  })
  reference <- fits[[1L]]
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(reference), vcov_type,
                                      vcov_defaulted)
  alignments <- lapply(fits, .multilpa_pool_align, reference = reference)
  aligned <- lapply(alignments, `[[`, "fit")
  tables <- lapply(seq_along(aligned), function(index) {
    .multilpa_pool_step(index, "inference", {
      parameter_inference(aligned[[index]], vcov_type = vcov_type,
                          adjust = "none")
    })
  })
  keys <- lapply(tables, .multilpa_pool_keys)
  if (!all(vapply(keys, identical, logical(1), keys[[1L]]))) {
    stop(errorCondition(paste(
      "The imputations' fits do not report the same parameters, so they cannot",
      "be pooled. This happens when a categorical indicator has a category in",
      "one completed data set that another lacks."),
      class = "latents_pooling_failed", call = NULL))
  }
  pooled <- .multilpa_rubin(tables, level, adjust)
  structure(list(
    estimates = pooled,
    imputations = .multilpa_pool_long(tables),
    fits = .multilpa_pool_fit_frame(fits, alignments),
    covariance = .multilpa_rubin_covariance(tables, aligned, vcov_type),
    m = length(completed), level = level, vcov_type = vcov_type,
    adjust = adjust, model = class(reference)[1L],
    n_profiles = reference$n_profiles,
    n_group_classes = reference$n_group_classes, call = call),
    class = "latents_pooled")
}

#' The completed data sets, as a checked list of data frames
#' @param imputed A list of data frames or a `mids` object.
#' @return A list of at least two data frames.
#' @noRd
.multilpa_completed_data <- function(imputed) {
  bad <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  if (inherits(imputed, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop(errorCondition(paste(
        "A `mids` object needs the mice package to be completed.",
        "Install it with install.packages(\"mice\")."),
        class = "latents_missing_package", call = NULL))
    }
    imputed <- mice::complete(imputed, action = "all")
  }
  if (!is.list(imputed) || is.data.frame(imputed) ||
      !all(vapply(imputed, is.data.frame, logical(1)))) {
    bad("`imputed` must be a list of completed data frames, or a `mids` object.")
  }
  if (length(imputed) < 2L) {
    bad("`imputed` must hold at least two completed data sets; one is not pooling.")
  }
  first <- imputed[[1L]]
  same_shape <- vapply(imputed, function(frame) {
    identical(names(frame), names(first)) && nrow(frame) == nrow(first)
  }, logical(1))
  if (!all(same_shape)) {
    bad(sprintf(paste(
      "Every completed data set must have the same columns and rows; set%s %s",
      "differ%s from the first."),
      if (sum(!same_shape) > 1L) "s" else "",
      paste(which(!same_shape), collapse = ", "),
      if (sum(!same_shape) > 1L) "" else "s"))
  }
  unname(imputed)
}

#' Run one imputation's step, naming the imputation if it fails
#'
#' Rubin's rules need every imputation, so a failure is not dropped: it is
#' raised again as `latents_pooling_failed`, with the imputation's number and
#' the original condition's class and message.
#' @param index The imputation's number.
#' @param stage `"fitting"` or `"inference"`.
#' @param expr The step.
#' @return The step's value.
#' @noRd
.multilpa_pool_step <- function(index, stage, expr) {
  tryCatch(expr, error = function(error) {
    stop(errorCondition(sprintf(
      "Imputation %d failed at %s (%s): %s", index, stage,
      class(error)[1L], conditionMessage(error)),
      class = "latents_pooling_failed", call = NULL))
  })
}

#' The identifying columns of an inference table, as one key per row
#' @param table A [parameter_inference()] table.
#' @return A character vector.
#' @noRd
.multilpa_pool_keys <- function(table) {
  paste(table$level, table$outcome, table$term, table$parameter, sep = "\r")
}

#' Pool inference tables with Rubin's rules
#'
#' Rubin (1987): the pooled estimate is the mean `Q`, the total variance
#' `T = U + (1 + 1/m) B` with `U` the mean within-imputation variance and `B`
#' the between-imputation variance, `r=(1 + 1/m) * B / U`, degrees of freedom
#' `(m - 1)(1 + 1/r)^2` (infinite when `B` is zero), and the fraction of
#' missing information `(r + 2 / (df + 3)) / (r + 1)`.
#'
#' @param tables Aligned [parameter_inference()] tables, one per imputation.
#' @param level Interval level.
#' @param adjust A [stats::p.adjust()] method for `p_adjusted`.
#' @return A data frame, one row per parameter.
#' @noRd
.multilpa_rubin <- function(tables, level, adjust = "none") {
  m <- length(tables)
  estimates <- vapply(tables, `[[`, numeric(nrow(tables[[1L]])), "estimate")
  variances <- vapply(tables, function(table) table$standard_error^2,
                      numeric(nrow(tables[[1L]])))
  estimates <- matrix(estimates, ncol = m)
  variances <- matrix(variances, ncol = m)
  pooled <- rowMeans(estimates)
  within <- rowMeans(variances)
  between <- apply(estimates, 1L, stats::var)
  total <- within + (1 + 1 / m) * between
  riv <- (1 + 1 / m) * between / within
  riv[!is.na(within) & !is.na(between) & within == 0 & between == 0] <- 0
  df <- ifelse(between > 0, (m - 1) * (1 + 1 / riv)^2, Inf)
  fmi <- ifelse(is.infinite(riv), 1, (riv + 2 / (df + 3)) / (riv + 1))
  standard_error <- sqrt(total)
  ## A test is reported only where the imputations' tables report one.
  tested <- !is.na(tables[[1L]]$statistic)
  statistic <- ifelse(tested, pooled / standard_error, NA_real_)
  p_value <- ifelse(tested, 2 * stats::pt(-abs(statistic), df), NA_real_)
  critical <- stats::qt((1 + level) / 2, df)
  data.frame(
    tables[[1L]][c("level", "outcome", "term", "parameter")],
    estimate = pooled, standard_error = standard_error,
    statistic = statistic, p_value = p_value,
    p_adjusted = stats::p.adjust(p_value, method = adjust),
    conf_low = pooled - critical * standard_error,
    conf_high = pooled + critical * standard_error,
    df = df, within = within, between = between, riv = riv, fmi = fmi,
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Pool the imputations' covariance matrices
#'
#' A covariate-free inference table carries its covariance; a covariate fit's
#' comes from [vcov()] with the same `vcov_type`. Either way its diagonal must
#' reproduce the table's squared standard errors, which is what ties its rows
#' to the table's rows.
#' @param tables Aligned [parameter_inference()] tables.
#' @param fits The aligned fits the tables came from.
#' @param vcov_type As in [parameter_inference()].
#' @return The total covariance `U + (1 + 1/m) B`, named by parameter.
#' @noRd
.multilpa_rubin_covariance <- function(tables, fits, vcov_type) {
  m <- length(tables)
  covariances <- lapply(seq_len(m), function(index) {
    attr(tables[[index]], "covariance") %||%
      vcov(fits[[index]], vcov_type = vcov_type)
  })
  invisible(lapply(seq_len(m), function(index) {
    stopifnot("each covariance must match its inference table row for row" =
                isTRUE(all.equal(sqrt(diag(covariances[[index]])),
                                 tables[[index]]$standard_error,
                                 check.attributes = FALSE, tolerance = 1e-8)))
  }))
  within <- Reduce(`+`, covariances) / m
  estimates <- do.call(rbind, lapply(tables, `[[`, "estimate"))
  total <- within + (1 + 1 / m) * stats::cov(estimates)
  dimnames(total) <- dimnames(covariances[[1L]])
  total
}

#' Every imputation's aligned estimates, one row per imputation and parameter
#' @param tables Aligned [parameter_inference()] tables.
#' @return A data frame.
#' @noRd
.multilpa_pool_long <- function(tables) {
  do.call(rbind, lapply(seq_along(tables), function(index) {
    data.frame(imputation = index,
               tables[[index]][c("level", "outcome", "term", "parameter",
                                 "estimate", "standard_error")],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

#' One row per imputation: its fit and the relabelling that aligned it
#' @param fits The unaligned fits.
#' @param alignments The alignments from `.multilpa_pool_align()`.
#' @return A data frame.
#' @noRd
.multilpa_pool_fit_frame <- function(fits, alignments) {
  data.frame(
    imputation = seq_along(fits),
    log_likelihood = vapply(fits, `[[`, numeric(1), "log_likelihood"),
    n_parameters = vapply(fits, function(fit) as.numeric(fit$n_parameters),
                          numeric(1)),
    converged = vapply(fits, function(fit) isTRUE(fit$converged), logical(1)),
    profile_order = vapply(alignments, function(alignment) {
      paste(alignment$profile_order, collapse = ",")
    }, character(1)),
    group_order = vapply(alignments, function(alignment) {
      paste(alignment$group_order, collapse = ",")
    }, character(1)),
    stringsAsFactors = FALSE)
}

#' @noRd
#' @export
get_results.latents_pooled <- function(x, what = c("estimates", "imputations", "fits"),
                                       ...) {
  what <- match.arg(what)
  x[[what]]
}

#' @noRd
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_pooled <- function(x, row.names = NULL, optional = FALSE,
                                         what = "estimates", ...) {
  get_results.latents_pooled(x, what = what)
}

#' @noRd
#' @export
print.latents_pooled <- function(x, digits = 4L, rows = 20L, ...) {
  cat(sprintf(paste(
    "Pooled over %d imputations (Rubin's rules): %s, %d profiles,",
    "%d group classes\n"), x$m, x$model, x$n_profiles, x$n_group_classes))
  cat(sprintf("Largest fraction of missing information: %.3f\n",
              max(x$estimates$fmi)))
  shown <- utils::head(x$estimates[c("level", "outcome", "term", "parameter",
                                     "estimate", "standard_error", "conf_low",
                                     "conf_high", "fmi")], rows)
  print(shown, digits = digits, row.names = FALSE)
  if (nrow(x$estimates) > rows) {
    cat(sprintf("... %d more rows: as.data.frame(x)\n",
                nrow(x$estimates) - rows))
  }
  invisible(x)
}

#' @noRd
#' @export
summary.latents_pooled <- function(object, ...) {
  fits <- object$fits
  identity_profiles <- paste(seq_len(object$n_profiles), collapse = ",")
  identity_groups <- paste(seq_len(object$n_group_classes), collapse = ",")
  data.frame(
    m = object$m, model = object$model, n_profiles = object$n_profiles,
    n_group_classes = object$n_group_classes,
    all_converged = all(fits$converged),
    log_likelihood_min = min(fits$log_likelihood),
    log_likelihood_max = max(fits$log_likelihood),
    relabelled = sum(fits$profile_order != identity_profiles |
                       fits$group_order != identity_groups),
    fmi_max = max(object$estimates$fmi),
    fmi_median = stats::median(object$estimates$fmi),
    stringsAsFactors = FALSE)
}

#' @rdname latents-model-methods
#' @export
vcov.latents_pooled <- function(object, ...) {
  object$covariance
}

#' @rdname latents-plots
#' @export
plot.latents_pooled <- function(x, what = c("profile", "group", "measurement"),
                                main = NULL, subtitle = NULL, level = NULL,
                                ...) {
  stopifnot("`x` must be a `latents_pooled` result" = inherits(x, "latents_pooled"))
  .latents_warn_unused_arguments(
    list(...), "plot()",
    "Choose the parameters with `what`, and style the plot with ggplot2.")
  what <- .latents_former_argument(what, !missing(what), level,
                                   c("what", "level"))
  what <- match.arg(what, c("profile", "group", "measurement"))
  estimates <- x$estimates[x$estimates$level %in% what, , drop = FALSE]
  imputations <- x$imputations[x$imputations$level %in% what, , drop = FALSE]
  if (nrow(estimates) == 0L) {
    stop(errorCondition(sprintf("This fit has no `%s` parameters to plot.", what),
                        class = "latents_bad_argument", call = NULL))
  }
  .gg_require()
  label <- function(table) paste(table$outcome, table$term, table$parameter)
  sources <- c("Pooled", "One imputation")
  rows <- rbind(
    data.frame(label = label(imputations), estimate = imputations$estimate,
               conf_low = NA_real_, conf_high = NA_real_, source = sources[2L]),
    data.frame(label = label(estimates), estimate = estimates$estimate,
               conf_low = estimates$conf_low, conf_high = estimates$conf_high,
               source = sources[1L]))
  # The first parameter at the top, as the tables list them.
  rows$label <- factor(rows$label, levels = rev(unique(label(estimates))))
  rows$source <- factor(rows$source, levels = sources)
  pooled <- rows[rows$source == sources[1L], , drop = FALSE]
  colours <- stats::setNames(.gg_okabe_ito[c(1L, 2L)], sources)
  drawn <- ggplot2::ggplot(rows, ggplot2::aes(estimate, label, colour = source,
                                     shape = source)) +
    ggplot2::geom_segment(data = pooled,
                          ggplot2::aes(x = conf_low, xend = conf_high,
                                       y = label, yend = label),
                          linewidth = 0.9, show.legend = FALSE) +
    ggplot2::geom_point(size = 2.6, stroke = 0.8) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16L, 1L), sources)) +
    .gg_theme() + .gg_top_legend() +
    ggplot2::labs(x = sprintf("Estimate with pooled %g%% interval", 100 * x$level),
                  y = NULL) +
    .gg_titles(main, subtitle, sprintf("Pooled over %d imputations", x$m),
               sprintf("%s%s parameters; open circles are the single imputations",
                       toupper(substr(what, 1L, 1L)), substring(what, 2L)))
  .gg_draw_on_call(drawn)
}
