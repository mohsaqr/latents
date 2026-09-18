#' Classification entropy sum
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @return Sum of `-p log p` over every cell, the classification entropy.
#' @noRd
.multilpa_entropy_sum <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  positive <- probabilities > 0
  -sum(probabilities[positive] * log(probabilities[positive]))
}

#' Relative classification entropy
#' @param probabilities Posterior probability matrix, rows summing to one.
#' @return Entropy on the zero-to-one scale, or `NA_real_` for a single class.
#' @noRd
.multilpa_relative_entropy <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  if (ncol(probabilities) == 1L) return(NA_real_)
  1 - .multilpa_entropy_sum(probabilities) /
    (nrow(probabilities) * log(ncol(probabilities)))
}

#' Information criteria for a fitted multilevel latent profile model
#'
#' Returns every information criterion the package computes, in one tidy table,
#' with the sample-size convention stated on each row. Multilevel mixtures admit
#' two defensible sample sizes: the number of independent groups and the number
#' of individuals. Criteria that depend on a sample size are therefore reported
#' once per convention rather than silently fixing one.
#'
#'   The verb is named for what it returns. Classification sharpness and class
#'   separation are reported by [classification_table()] and [entropy_table()]
#'   rather than folded into this table, so `information_criteria()` stays a
#'   pure likelihood-penalty summary.
#'
#' @param object A fitted `multilpa` model.
#' @seealso [classification_table()] and [entropy_table()] for the
#'   classification diagnostics that accompany these criteria.
#' @return A base `data.frame` with one row per criterion and convention, and
#'   columns `criterion`, `convention`, `n`, `value`, `penalty`, and
#'   `definition`. `log_likelihood` and `aic` do not depend on a sample size and
#'   carry convention `"none"` with `n` equal to `NA_integer_`. Lower values
#'   indicate a preferred model for every criterion except `log_likelihood`.
#' @details Let `q` be the number of free parameters, `n` the chosen sample
#'   size, and `EN` the classification entropy of the level matching that
#'   convention. The criteria are `aic = -2L + 2q`, `bic = -2L + q log(n)`,
#'   `sabic = -2L + q log((n + 2) / 24)`, `caic = -2L + q (log(n) + 1)`,
#'   `awe = -2(L - EN) + 2q (1.5 + log(n))`, and `icl = -2L + q log(n) + 2 EN`.
#'   The entropy-based criteria use group-level posteriors under the `"groups"`
#'   convention and individual-level posteriors under the `"individuals"`
#'   convention. This is a stated per-level choice, not a unique multilevel
#'   definition; use one convention consistently across compared candidates.
#' @references Schwarz, G. (1978). Estimating the dimension of a model.
#'   Annals of Statistics, 6, 461--464. Sclove, S. L. (1987). Application of
#'   model-selection criteria to some problems in multivariate analysis.
#'   Psychometrika, 52, 333--343. Bozdogan, H. (1987). Model selection and
#'   Akaike's information criterion. Psychometrika, 52, 345--370.
#'   Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
#'   non-Gaussian clustering. Biometrics, 49, 803--821. Biernacki, C., Celeux,
#'   G., & Govaert, G. (2000). Assessing a mixture model for clustering with the
#'   integrated completed likelihood. IEEE Transactions on Pattern Analysis and
#'   Machine Intelligence, 22, 719--725.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' information_criteria(fit)
#' @export
information_criteria <- function(object) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object))
  q <- object$n_parameters
  log_likelihood <- object$log_likelihood
  conventions <- data.frame(
    convention = c("groups", "individuals"),
    n = c(object$n_groups, object$n_observations))
  entropies <- c(groups = .multilpa_entropy_sum(object$group_posteriors),
                 individuals = .multilpa_entropy_sum(object$subject_posteriors))
  scale_free <- data.frame(
    criterion = c("log_likelihood", "aic"),
    convention = "none", n = NA_integer_,
    value = c(log_likelihood, -2 * log_likelihood + 2 * q),
    penalty = c(NA_real_, 2 * q),
    definition = c("maximized observed-data log likelihood", "-2L + 2q"))
  scaled <- do.call(rbind, lapply(seq_len(nrow(conventions)), function(row) {
    convention <- conventions$convention[row]
    n <- conventions$n[row]
    entropy_sum <- entropies[[convention]]
    penalties <- c(bic = q * log(n),
                   sabic = q * log((n + 2) / 24),
                   caic = q * (log(n) + 1),
                   awe = 2 * q * (1.5 + log(n)) + 2 * entropy_sum,
                   icl = q * log(n) + 2 * entropy_sum)
    data.frame(criterion = names(penalties), convention = convention, n = n,
               value = -2 * log_likelihood + penalties,
               penalty = unname(penalties),
               definition = c("-2L + q log(n)", "-2L + q log((n + 2) / 24)",
                              "-2L + q (log(n) + 1)",
                              "-2(L - EN) + 2q (1.5 + log(n))",
                              "-2L + q log(n) + 2 EN"))
  }))
  result <- rbind(scale_free, scaled)
  row.names(result) <- NULL
  result
}

#' Classification quality for a fitted multilevel latent profile model
#'
#' Summarizes how sharply the posterior probabilities separate classes, at the
#' individual level, the group level, or both. Reports the diagnostics normally
#' expected alongside a mixture solution: modal counts, model-estimated class
#' sizes, average posterior probability in the assigned class, and the odds of
#' correct classification.
#'
#' @param object A fitted `multilpa` model.
#' @param level `"individuals"` for latent profiles, `"groups"` for latent group
#'   classes, or `"both"` to stack them in one table.
#' @param detail `FALSE` returns one row per class. `TRUE` returns the full
#'   cross-tabulation, one row per assigned class and class, which is the
#'   average-posterior matrix reported by mixture software.
#' @return A base `data.frame`. With `detail = FALSE` there is one row per class
#'   and level, with columns `level`, `class`, `n_modal`, `proportion_modal`,
#'   `estimated_n`, `estimated_proportion`, `average_posterior`, and
#'   `odds_correct_classification`. With `detail = TRUE` there is one row per
#'   `level`, `assigned_class`, and `class`, with columns `n_assigned` and
#'   `average_posterior`. `average_posterior` is the mean posterior probability
#'   of `class` among individuals or groups whose modal class is
#'   `assigned_class`; the diagonal of that matrix is the summary column of the
#'   same name. Rows for classes with no modal members carry `NA_real_`
#'   averages, and `odds_correct_classification` is `NA_real_` when it is not
#'   defined, that is when the average posterior or the estimated proportion is
#'   zero or one.
#' @details The odds of correct classification for class `k` is
#'   `(p / (1 - p)) / (r / (1 - r))`, where `p` is the average posterior in the
#'   assigned class and `r` is the model-estimated class proportion. Values near
#'   one indicate classification no better than the class proportion alone;
#'   values of five or more are conventionally read as adequate separation.
#'   Modal assignment discards classification uncertainty, so `n_modal` and
#'   `estimated_n` differ whenever entropy is below one.
#' @references Nagin, D. S. (2005). Group-Based Modeling of Development.
#'   Harvard University Press.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' classification_table(fit)
#' @export
classification_table <- function(object, level = c("individuals", "groups", "both"),
                                 detail = FALSE) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object),
            "`detail` must be TRUE or FALSE" =
              isTRUE(detail) || isFALSE(detail))
  level <- match.arg(level)
  levels_wanted <- if (level == "both") c("individuals", "groups") else level
  posteriors <- list(individuals = object$subject_posteriors,
                     groups = object$group_posteriors)
  result <- do.call(rbind, lapply(levels_wanted, function(which_level) {
    probabilities <- posteriors[[which_level]]
    .multilpa_classification_rows(probabilities, which_level, detail)
  }))
  row.names(result) <- NULL
  result
}

#' Assemble classification diagnostics for one level
#' @param probabilities Posterior probability matrix.
#' @param level Label recorded in the `level` column.
#' @param detail Whether to return the full average-posterior cross-tabulation.
#' @return A tidy `data.frame` of classification diagnostics.
#' @noRd
.multilpa_classification_rows <- function(probabilities, level, detail) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  n_classes <- ncol(probabilities)
  classes <- seq_len(n_classes)
  assigned <- max.col(probabilities, ties.method = "first")
  n_modal <- tabulate(assigned, nbins = n_classes)
  # Mean posterior for every class within each modal assignment group.
  average <- t(vapply(classes, function(class) {
    rows <- assigned == class
    if (!any(rows)) rep(NA_real_, n_classes) else
      colMeans(probabilities[rows, , drop = FALSE])
  }, numeric(n_classes)))
  if (isTRUE(detail)) {
    return(data.frame(
      level = level,
      assigned_class = rep(classes, each = n_classes),
      class = rep(classes, times = n_classes),
      n_assigned = rep(n_modal, each = n_classes),
      average_posterior = as.vector(t(average))))
  }
  estimated_n <- colSums(probabilities)
  estimated_proportion <- estimated_n / nrow(probabilities)
  diagonal <- average[cbind(classes, classes)]
  odds <- .multilpa_odds_correct(diagonal, estimated_proportion)
  data.frame(level = level, class = classes, n_modal = n_modal,
             proportion_modal = n_modal / nrow(probabilities),
             estimated_n = estimated_n,
             estimated_proportion = estimated_proportion,
             average_posterior = diagonal,
             odds_correct_classification = odds)
}

#' Odds of correct classification
#' @param average Average posterior probability in the assigned class.
#' @param proportion Model-estimated class proportion.
#' @return Odds ratio, or `NA_real_` where the ratio is undefined.
#' @noRd
.multilpa_odds_correct <- function(average, proportion) {
  stopifnot("`average` and `proportion` must have equal length" =
              length(average) == length(proportion))
  defined <- !is.na(average) & average > 0 & average < 1 &
    proportion > 0 & proportion < 1
  odds <- rep(NA_real_, length(average))
  odds[defined] <- (average[defined] / (1 - average[defined])) /
    (proportion[defined] / (1 - proportion[defined]))
  odds
}

#' Relative entropy of a fitted multilevel latent profile model
#'
#' Reports the zero-to-one relative entropy at each level, alongside the raw
#' classification entropy the entropy-penalized information criteria use.
#'
#' @param object A fitted `multilpa` model.
#' @return A base `data.frame` with one row per level and columns `level`,
#'   `n_classes`, `n_units`, `entropy_sum`, and `relative_entropy`.
#'   `relative_entropy` is `1 - EN / (n log K)` and is `NA_real_` when a level
#'   has a single class, where it is undefined rather than perfect.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' entropy_table(fit)
#' @export
entropy_table <- function(object) {
  stopifnot("`object` must be a fitted model of this package" = .multilpa_any_fit(object))
  posteriors <- list(individuals = object$subject_posteriors,
                     groups = object$group_posteriors)
  result <- do.call(rbind, lapply(names(posteriors), function(level) {
    probabilities <- posteriors[[level]]
    data.frame(level = level, n_classes = ncol(probabilities),
               n_units = nrow(probabilities),
               entropy_sum = .multilpa_entropy_sum(probabilities),
               relative_entropy = .multilpa_relative_entropy(probabilities))
  }))
  row.names(result) <- NULL
  result
}
