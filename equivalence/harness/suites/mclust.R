# Gaussian single-level limit against mclust (Scrucca, Fraley, Murphy & Raftery
# 2023, Model-Based Clustering, Classification, and Density Estimation Using
# mclust in R. Chapman & Hall/CRC).
#
# With one group class and continuous indicators multilpa is an ordinary
# Gaussian mixture, so mclust's parameterizations map onto its two covariance
# switches:
#
#   mclust EEI <-> covariance_model "diagonal", variance_model "equal"
#   mclust VVI <-> covariance_model "diagonal", variance_model "varying"
#   mclust EEE <-> covariance_model "full",     variance_model "equal"
#   mclust VVV <-> covariance_model "full",     variance_model "varying"
#
# mclust reports BIC on the opposite sign convention, as 2 * logLik - df *
# log(n), so its value is the negative of multilpa's `bic_individual`.
#
# Published targets are the book's Example 3.2 (thyroid) and the mclust
# vignette's density-estimation and clustering sections (acidity, faithful).
# The `diabetes` data are deliberately excluded: mclust's own help page flags
# them as flawed, and three mutually inconsistent published analyses exist.

#' Map an mclust model name onto multilpa's covariance switches
#' @param model An mclust model name.
#' @return A list with `covariance_model` and `variance_model`.
.mclust_covariance <- function(model) {
  switch(model,
         EEI = list(covariance_model = "diagonal", variance_model = "equal"),
         VVI = list(covariance_model = "diagonal", variance_model = "varying"),
         EEE = list(covariance_model = "full", variance_model = "equal"),
         VVV = list(covariance_model = "full", variance_model = "varying"),
         E   = list(covariance_model = "diagonal", variance_model = "equal"),
         V   = list(covariance_model = "diagonal", variance_model = "varying"),
         stop(sprintf("no multilpa equivalent for mclust model '%s'", model)))
}

suite_mclust <- function() {
  require_suite_packages(
    "mclust",
    reason = "the Gaussian-mixture reference and the acidity and thyroid data")
  # Mclust() resolves its model-search function by name, so the namespace has to
  # be attached rather than only loaded.
  attached <- "package:mclust" %in% search()
  if (!attached) {
    suppressMessages(library(mclust))
    on.exit(detach("package:mclust", unload = FALSE), add = TRUE)
  }
  cases <- list(
    list(dataset = "thyroid", variables = c("RT3U", "T4", "T3", "TSH", "DTSH"),
         model = "VVI", g = 3L,
         published = c(bic = -4777.9069, log_likelihood = -2303.0232, df = 32),
         published_digits = c(4, 4, 0),
         citation = "Scrucca et al. (2023) mclust book, Example 3.2"),
    list(dataset = "acidity", variables = NULL, model = "E", g = 2L,
         published = c(bic = -392.0723108, log_likelihood = -185.9493052, df = 4),
         published_digits = c(7, 7, 0),
         citation = "mclust vignette, density estimation"),
    list(dataset = "faithful", variables = c("eruptions", "waiting"),
         model = "EEE", g = 3L,
         published = c(bic = -2314.316, log_likelihood = -1126.326, df = 11),
         published_digits = c(3, 3, 0),
         citation = "mclust vignette, clustering"),
    list(dataset = "thyroid", variables = c("RT3U", "T4", "T3", "TSH", "DTSH"),
         model = "VVV", g = 3L, published = NULL, published_digits = NULL,
         citation = NULL),
    list(dataset = "thyroid", variables = c("RT3U", "T4", "T3", "TSH", "DTSH"),
         model = "EEI", g = 2L, published = NULL, published_digits = NULL,
         citation = NULL),
    list(dataset = "faithful", variables = c("eruptions", "waiting"),
         model = "VVV", g = 2L, published = NULL, published_digits = NULL,
         citation = NULL)
  )
  do.call(rbind, lapply(cases, function(case) do.call(.mclust_case, case)))
}

#' Compare one mclust model against multilpa
#' @param dataset Dataset name, from mclust or datasets.
#' @param variables Columns to use, or `NULL` for a univariate dataset.
#' @param model mclust model name.
#' @param g Number of mixture components.
#' @param published Named numeric vector of published values, or `NULL`.
#' @param published_digits Decimals each published value was printed to.
#' @param citation Source of the published values, or `NULL`.
#' @return A `data.frame` of compared quantities.
.mclust_case <- function(dataset, variables, model, g, published,
                         published_digits, citation) {
  package <- if (dataset %in% c("faithful")) "datasets" else "mclust"
  raw <- get(utils::data(list = dataset, package = package,
                         envir = environment()), envir = environment())
  frame <- if (is.null(variables)) data.frame(y = as.numeric(raw)) else
    as.data.frame(raw)[, variables, drop = FALSE]
  indicators <- names(frame)
  frame$unit <- factor(seq_len(nrow(frame)))

  # Mclust() names its columns from deparse(substitute(data)), so the argument
  # has to be a plain symbol: an inline expression deparses to several lines and
  # `colnames<-` then fails on univariate data.
  model_data <- if (is.null(variables)) as.numeric(raw) else
    as.matrix(frame[, indicators, drop = FALSE])
  # Two mclust fits, because two different claims are being made. The converged
  # fit is what multilpa is compared against; the default-tolerance fit is what
  # the published values were produced from. mclust's default EM tolerance is
  # 1e-5 on the log-likelihood and stops measurably short of the maximum: on
  # faithful it gives -1126.326236 where the converged value is -1126.315928.
  # Comparing multilpa's converged estimate against a published under-converged
  # one would report that 0.010 gap as if it were a disagreement between
  # implementations.
  reference <- mclust::Mclust(model_data, G = g, modelNames = model,
                              verbose = FALSE,
                              control = mclust::emControl(tol = c(1e-13, 1e-13),
                                                          itmax = c(1e5, 1e5)))
  as_published <- mclust::Mclust(model_data, G = g, modelNames = model,
                                 verbose = FALSE)
  covariance <- .mclust_covariance(model)
  fit <- multilpa(frame, vars = indicators, id = "unit", n_profiles = g,
                  n_group_classes = 1L, n_starts = 60L, seed = 1L,
                  tol = 1e-13, max_iter = 20000L,
                  covariance_model = covariance$covariance_model,
                  variance_model = covariance$variance_model)

  # mclust stores means as indicators-by-component for multivariate data and as
  # a plain vector for univariate data; both become components-by-indicators.
  reference_means <- matrix(as.numeric(reference$parameters$mean),
                            nrow = g, byrow = TRUE)
  profile_order <- align_classes(reference_means, unname(fit$means))
  label <- sprintf("%s, %s, %d components", dataset, model, g)

  # mclust initializes from model-based hierarchical clustering and takes no
  # random restarts, so the two programs can land on different modes. Where they
  # do, parameters are not comparable; what is comparable, and what is checked
  # instead, is that multilpa's maximum is no lower.
  shortfall <- reference$loglik - as.numeric(logLik(fit))
  same_mode <- abs(shortfall) <= 1e-4

  rows <- rbind(
    compare_values(sprintf("%s: free parameters", label),
                   reference$df, fit$n_parameters, tolerance = 0.5),
    compare_values(sprintf("%s: log-likelihood shortfall against mclust", label),
                   0, max(shortfall, 0), tolerance = 1e-4))
  if (same_mode) {
    rows <- rbind(rows,
      compare_values(sprintf("%s: maximised log-likelihood", label),
                     reference$loglik, as.numeric(logLik(fit)), tolerance = 1e-8),
      compare_values(sprintf("%s: BIC on mclust's sign convention", label),
                     reference$bic, -fit$bic_individual, tolerance = 1e-7),
      compare_values(sprintf("%s: mixing proportion %d", label, seq_len(g)),
                     sort(as.numeric(reference$parameters$pro)),
                     sort(as.numeric(fit$profile_probabilities)), tolerance = 1e-6),
      compare_values(sprintf("%s: mean of %s in component %d", label,
                             rep(indicators, each = g), rep(seq_len(g), length(indicators))),
                     as.numeric(reference_means),
                     as.numeric(unname(fit$means)[profile_order, , drop = FALSE]),
                     tolerance = 1e-6, scale = "relative"))
  }
  rows <- transform(rows, source = sprintf("mclust %s", utils::packageVersion("mclust")),
                    dataset = dataset, precision = "machine")

  if (is.null(published)) return(rows)
  # The published value is checked against the mclust run that produced it, at
  # the precision it was printed to. This establishes where the number came
  # from; it is not a second equivalence claim, and stating it against the
  # converged fit would charge mclust's stopping rule to multilpa.
  rbind(rows, compare_printed(
    sprintf("%s: published %s reproduced by mclust at its default tolerance",
            label, c("BIC", "log-likelihood", "free parameters")),
    as.numeric(published[c("bic", "log_likelihood", "df")]),
    c(as_published$bic, as_published$loglik, as_published$df),
    digits = published_digits) |>
      transform(source = citation, dataset = dataset, precision = "printed"))
}
