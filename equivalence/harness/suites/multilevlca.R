# Two-level latent class analysis against multilevLCA (Lyrvall, Di Mari, Bakk,
# Oser & Kuha 2025, Multivariate Behavioral Research 60(4):731-747).
#
# multilevLCA fits the same nonparametric two-level mixture and is the only
# other CRAN package to publish peer-reviewed worked numbers on bundled public
# data. Its parameter blocks map onto multilpa's as
#
#   vOmega <-> group_probabilities      P(group class)
#   mPi    <-> profile_probabilities    P(profile | group class), transposed
#   mPhi   <-> response_probabilities   P(item = 1 | profile), binary items
#
# and its two information criteria correspond to multilpa's two sample-size
# conventions: BIClow to bic_individual, BIChigh to bic_groups.
#
# `fixedpars = 0` is required. multilevLCA defaults to `fixedpars = 1`, a
# two-step estimator whose measurement model is held fixed from step one, which
# is deliberately not the one-step maximum likelihood multilpa computes.

suite_multilevlca <- function() {
  require_suite_packages(
    "multilevLCA",
    reason = "the third independent two-level implementation and the dataTOY and dataIEA data")
  rbind(
    .multilevlca_case(3L, 2L), .multilevlca_case(2L, 2L), .multilevlca_case(3L, 3L),
    .multilevlca_published()
  )
}

#' Compare one dataTOY configuration against multilevLCA
#' @param profiles Number of lower-level latent classes.
#' @param group_classes Number of higher-level latent classes.
#' @return A `data.frame` of compared quantities.
.multilevlca_case <- function(profiles, group_classes) {
  observed <- get(utils::data(list = "dataTOY", package = "multilevLCA",
                              envir = environment()), envir = environment())
  items <- paste0("Y_", seq_len(10))
  set.seed(2023L)
  reference <- multilevLCA::multiLCA(observed, Y = items, iT = profiles,
                                     id_high = "id_high", iM = group_classes,
                                     fixedpars = 0, verbose = FALSE)
  fit <- multilpa(observed, vars = items, id = "id_high",
                  n_profiles = profiles, n_group_classes = group_classes,
                  categorical = items, n_starts = 30L, seed = 1L,
                  tol = 1e-12, max_iter = 20000L)

  # Binary items, so the second column is P(item = 1 | profile).
  positives <- flatten_blocks(lapply(fit$response_probabilities,
                                     \(block) unname(as.matrix(block))[, 2L, drop = FALSE]))
  profile_order <- align_classes(t(unname(reference$mPhi)), positives)
  prevalence <- unname(fit$profile_probabilities)[, profile_order, drop = FALSE]
  group_order <- align_classes(t(unname(reference$mPi)), prevalence)

  label <- sprintf("dataTOY, %d profiles x %d group classes", profiles, group_classes)
  rbind(
    compare_values(sprintf("%s: maximised log-likelihood", label),
                   utils::tail(reference$LLKSeries, 1L), as.numeric(logLik(fit)),
                   tolerance = 1e-5),
    compare_values(sprintf("%s: AIC", label), reference$AIC, fit$aic, tolerance = 1e-4),
    compare_values(sprintf("%s: BIC, individual sample size", label),
                   reference$BIClow, fit$bic_individual, tolerance = 1e-4),
    compare_values(sprintf("%s: BIC, group sample size", label),
                   reference$BIChigh, fit$bic_groups, tolerance = 1e-4),
    compare_values(sprintf("%s: group-class probability %d", label, seq_len(group_classes)),
                   sort(as.numeric(reference$vOmega)),
                   sort(as.numeric(fit$group_probabilities)), tolerance = 1e-4),
    compare_values(sprintf("%s: P(profile %d | group class %d)", label,
                           rep(seq_len(profiles), each = group_classes),
                           rep(seq_len(group_classes), profiles)),
                   as.numeric(t(unname(reference$mPi))),
                   as.numeric(prevalence[group_order, , drop = FALSE]), tolerance = 1e-4),
    compare_values(sprintf("%s: P(%s = 1 | profile %d)", label,
                           rep(items, each = profiles), rep(seq_len(profiles), length(items))),
                   as.numeric(t(unname(reference$mPhi))),
                   as.numeric(positives[profile_order, , drop = FALSE]), tolerance = 1e-4)
  ) |> transform(source = sprintf("multilevLCA %s",
                                  utils::packageVersion("multilevLCA")),
                 dataset = "dataTOY", precision = "machine")
}

#' Published two-level results for the ICCS 2016 citizenship data
#'
#' Lyrvall et al. (2025) report a five-class, two-group-class model for twelve
#' binary citizenship-norm items measured on students within countries.
#' multilevLCA deletes incomplete rows, leaving 83,060 of 90,221.
#'
#' multilpa is scored at multilevLCA's estimates rather than re-optimized. That
#' is both faster and a stronger claim: it tests whether the two packages write
#' down the same likelihood and count the same free parameters, instead of
#' testing whether two different multi-start schemes reach the same mode on a
#' five-by-two landscape.
#'
#' @return A `data.frame` of compared quantities.
.multilevlca_published <- function() {
  observed <- get(utils::data(list = "dataIEA", package = "multilevLCA",
                              envir = environment()), envir = environment())
  items <- c("obey", "rights", "local", "work", "envir", "vote",
             "history", "respect", "news", "protest", "discuss", "party")
  complete <- observed[stats::complete.cases(observed[, c(items, "COUNTRY")]), , drop = FALSE]
  set.seed(2023L)
  reference <- multilevLCA::multiLCA(observed, Y = items, iT = 5L,
                                     id_high = "COUNTRY", iM = 2L,
                                     fixedpars = 0, verbose = FALSE)
  # mPhi holds P(item = 1 | class); multilpa carries both categories in sorted
  # level order, so category one is the zero response.
  responses <- lapply(seq_along(items), function(item) {
    positive <- as.numeric(reference$mPhi[item, ])
    # `cbind()` names the second column after the expression it was given, and
    # multilpa refuses a start whose response labels do not match the data's
    # own category levels. The matrix carries no labels.
    unname(cbind(1 - positive, positive))
  })
  evaluated <- multilpa(complete, vars = items, id = "COUNTRY",
                        n_profiles = 5L, n_group_classes = 2L, categorical = items,
                        n_starts = 1L, max_iter = 0L,
                        start = starting_values(list(
                          profile_probabilities = t(unname(reference$mPi)),
                          group_probabilities = as.numeric(reference$vOmega),
                          response_probabilities = responses)))
  against_package <- rbind(
    compare_values("dataIEA 5 profiles x 2 group classes: analysis sample size",
                   83060, evaluated$n_informative, tolerance = 0.5),
    compare_values("dataIEA 5 profiles x 2 group classes: log-likelihood at multilevLCA's estimates",
                   utils::tail(reference$LLKSeries, 1L), as.numeric(logLik(evaluated)),
                   tolerance = 1e-4)
  ) |> transform(source = sprintf("multilevLCA %s",
                                  utils::packageVersion("multilevLCA")),
                 dataset = "dataIEA", precision = "machine")
  against_paper <- compare_printed(
    sprintf("dataIEA 5 profiles x 2 group classes: published %s",
            c("BIC, individual sample size", "BIC, group sample size", "AIC")),
    c(869122.92, 868554.62, 868479.34),
    c(evaluated$bic_individual, evaluated$bic_groups, evaluated$aic),
    digits = 2) |>
    transform(source = "Lyrvall et al. (2025) Multivar Behav Res 60(4):731-747",
              dataset = "dataIEA", precision = "printed")
  rbind(against_package, against_paper)
}
