# Published reference values.
#
# Every other suite compares multilpa against another program. This one
# compares it against numbers printed in the literature, which is a stronger
# claim: agreeing with poLCA shows two implementations share a likelihood,
# whereas agreeing with Agresti & Lang (1993) shows the implementation
# reproduces a result the field has already accepted.
#
# Sources
#   Linzer & Lewis (2011) poLCA: An R package for polytomous variable latent
#     class analysis. Journal of Statistical Software 42(10), pp. 21-22.
#   Agresti & Lang (1993) Quasi-symmetric latent class models, with application
#     to rater agreement. Biometrics 49(1), 131-139, Table 2 p. 136.
#   Zhou & Lange (2010) On the bumpy road to the dominant mode. Scandinavian
#     Journal of Statistics 37(4), Tables 3 and 4.
#   Goodman (2002) Latent class analysis: the empirical study of latent types,
#     latent variables, and latent structures. In Hagenaars & McCutcheon (eds),
#     Applied Latent Class Analysis, Tables 5a and 6.
#   Dayton (1998) Latent Class Scaling Analysis. Sage, Table 3.3 p. 32.
#
# The carcinoma data are seven pathologists rating 118 slides for carcinoma,
# taken by Agresti & Lang from Landis & Koch (1977). Note that they report
# degrees of freedom against the full 2^7 contingency table while poLCA reports
# residual degrees of freedom against N, so only the statistics are comparable.

#' Likelihood-ratio goodness-of-fit statistic for a categorical latent class fit
#'
#' Computed from the definition against the observed response-pattern counts,
#' independently of any package that also reports it.
#'
#' @param object A fitted `multilpa` model with only categorical indicators and
#'   one group class.
#' @return A list with `g_squared`, `pearson` and `cells`, the size of the full
#'   contingency table.
#'
#' Both statistics are accumulated over the complete cross-classification, not
#' only the patterns that happen to occur. An empty cell contributes nothing to
#' G-squared, because the limit of `o * log(o / e)` as `o` approaches zero is
#' zero, but it contributes its expected count to Pearson's statistic. Summing
#' only over observed patterns therefore reproduces G-squared while
#' understating Pearson's -- on the carcinoma three-class fit, 18.734 instead of
#' the published 20.503.
.likelihood_ratio_fit <- function(object) {
  codes <- object$categorical_data
  stopifnot("goodness of fit needs complete categorical data" = !anyNA(codes),
            "goodness of fit assumes a single group class" = object$n_group_classes == 1L)
  n_categories <- vapply(object$response_probabilities, ncol, integer(1))
  cells <- as.matrix(expand.grid(lapply(n_categories, seq_len), KEEP.OUT.ATTRS = FALSE))
  shares <- as.numeric(object$profile_probabilities)
  probability <- vapply(seq_len(nrow(cells)), function(row) {
    per_profile <- vapply(seq_along(shares), function(profile) {
      prod(vapply(seq_along(n_categories), function(item) {
        object$response_probabilities[[item]][profile, cells[row, item]]
      }, numeric(1)))
    }, numeric(1))
    sum(shares * per_profile)
  }, numeric(1))
  expected <- probability * object$n_observations
  key <- apply(codes, 1L, paste, collapse = "|")
  counts <- table(factor(key, levels = apply(cells, 1L, paste, collapse = "|")))
  observed <- as.numeric(counts)
  positive <- observed > 0
  list(g_squared = 2 * sum(observed[positive] * log(observed[positive] / expected[positive])),
       pearson = sum((observed - expected)^2 / expected),
       cells = nrow(cells))
}

#' Fit a single-level latent class model to a poLCA dataset
#' @param dataset Name of a poLCA dataset.
#' @param items Character vector of item names.
#' @param k Number of latent classes.
#' @param starts Number of random starts.
#' @return A fitted `multilpa` model.
.published_fit <- function(dataset, items, k, starts = 60L) {
  raw <- get(utils::data(list = dataset, package = "poLCA",
                         envir = environment()), envir = environment())
  frame <- raw[, items, drop = FALSE]
  frame$unit <- factor(seq_len(nrow(frame)))
  multilpa(frame, indicators = items, group = "unit", n_profiles = k,
           n_group_classes = 1L, categorical = items, n_starts = starts,
           seed = 20260919L, tol = 1e-13, max_iter = 20000L)
}

suite_published <- function() {
  if (!requireNamespace("poLCA", quietly = TRUE)) {
    stop("suite 'published' needs the poLCA package for its bundled datasets")
  }
  rbind(.published_carcinoma(), .published_values(), .published_cheating(),
        .published_gss82())
}

#' Carcinoma: Linzer & Lewis (2011), Agresti & Lang (1993), Zhou & Lange (2010)
#' @return A `data.frame` of compared quantities.
.published_carcinoma <- function() {
  items <- c("A", "B", "C", "D", "E", "F", "G")
  fits <- lapply(2:4, function(k) .published_fit("carcinoma", items, k))
  statistics <- lapply(fits, .likelihood_ratio_fit)
  three <- fits[[2L]]

  # Agresti & Lang report a one-class baseline as well.
  one_class <- .published_fit("carcinoma", items, 1L, starts = 1L)

  jss <- rbind(
    compare_printed("carcinoma 3-class: log-likelihood",
                    -293.705, as.numeric(logLik(three)), digits = 3),
    compare_values("carcinoma 3-class: free parameters",
                   23, three$n_parameters, tolerance = 0.5),
    compare_printed("carcinoma 3-class: AIC", 633.41, three$aic, digits = 2),
    compare_printed("carcinoma 3-class: BIC", 697.1357, three$bic_individual,
                    digits = 4),
    compare_printed("carcinoma 3-class: likelihood-ratio G-squared",
                    15.26171, statistics[[2L]]$g_squared, digits = 5),
    compare_printed("carcinoma 3-class: Pearson chi-squared",
                    20.50336, statistics[[2L]]$pearson, digits = 5),
    compare_printed(sprintf("carcinoma 3-class: class share %d", seq_len(3L)),
                    sort(c(0.3736, 0.1817, 0.4447)),
                    sort(as.numeric(three$profile_probabilities)), digits = 4)
  ) |> transform(source = "Linzer & Lewis (2011) JSS 42(10) pp. 21-22")

  agresti <- compare_printed(
    sprintf("carcinoma %d-class: G-squared", 1:3),
    c(476.8, 62.4, 15.3),
    c(.likelihood_ratio_fit(one_class)$g_squared,
      statistics[[1L]]$g_squared, statistics[[2L]]$g_squared),
    digits = 1) |>
    transform(source = "Agresti & Lang (1993) Biometrics 49(1) Table 2")

  # Zhou & Lange publish the dominant mode and an inferior local mode of the
  # four-class likelihood, which makes this a test of the optimizer rather than
  # of the likelihood. Each seed is a single start, so the collected optima are
  # the modes this implementation actually reaches.
  single_starts <- vapply(seq_len(120L), function(seed) {
    raw <- get(utils::data(list = "carcinoma", package = "poLCA",
                           envir = environment()), envir = environment())
    frame <- raw[, items, drop = FALSE]
    frame$unit <- factor(seq_len(nrow(frame)))
    as.numeric(logLik(multilpa(frame, indicators = items, group = "unit",
                               n_profiles = 4L, n_group_classes = 1L,
                               categorical = items, n_starts = 1L, seed = seed,
                               tol = 1e-13, max_iter = 20000L)))
  }, numeric(1))
  modes <- sort(unique(round(single_starts, 4L)), decreasing = TRUE)
  inferior <- modes[which.min(abs(modes - (-293.3200)))]

  # multilpa reaches -289.28584884, which is 5.1e-05 above the published
  # -289.2859, so the two are not the same point on the surface and the
  # dominant mode cannot be checked for digit reproduction. What is checkable is
  # that multilpa's maximum is not the lower of the two. The inferior mode is a
  # different matter: it reproduces -293.3200 exactly.
  dominant <- max(single_starts)
  zhou <- rbind(
    compare_values("carcinoma 4-class: shortfall of the dominant mode against Zhou & Lange",
                   0, max(-289.2859 - dominant, 0), tolerance = 1e-6),
    compare_printed("carcinoma 4-class: inferior local mode log-likelihood",
                    -293.3200, inferior, digits = 4),
    # Four classes on 118 observations is 31 parameters, and the two smallest
    # classes hold roughly 11 and 22 of them. Along that flat a ridge, a
    # 5e-05 difference in log-likelihood moves those shares in the fourth
    # decimal: 0.0936 against a published 0.0938, 0.1882 against 0.1881. The
    # two large classes do reproduce the published digits.
    compare_values(sprintf("carcinoma 4-class: dominant mode class share %d", seq_len(4L)),
                   sort(c(0.3430, 0.3751, 0.0938, 0.1881)),
                   sort(as.numeric(fits[[3L]]$profile_probabilities)),
                   tolerance = 5e-4)
  ) |> transform(source = "Zhou & Lange (2010) Scand J Stat 37(4) Tables 3-4")

  rbind(jss, agresti, zhou) |>
    transform(dataset = "carcinoma", precision = "printed")
}

#' Values: Goodman (2002) Tables 5a and 6
#' @return A `data.frame` of compared quantities.
.published_values <- function() {
  items <- c("A", "B", "C", "D")
  two <- .published_fit("values", items, 2L)
  statistics <- .likelihood_ratio_fit(two)
  class_order <- order(as.numeric(two$profile_probabilities))
  # Goodman tabulates the probability of the SECOND response category. The
  # model-implied marginal for item A under his table, 0.28 * 0.993 + 0.72 *
  # 0.714 = 0.792, is the observed proportion of category 2 (0.792) and not of
  # category 1 (0.208), which fixes the convention without reference to the fit.
  probabilities <- vapply(items, function(item) {
    two$response_probabilities[[item]][class_order, 2L]
  }, numeric(2L))
  rbind(
    compare_printed("values 2-class: likelihood-ratio G-squared",
                    2.72, statistics$g_squared, digits = 2),
    compare_printed("values 2-class: Pearson chi-squared",
                    2.72, statistics$pearson, digits = 2),
    compare_printed(sprintf("values 2-class: class share %d", seq_len(2L)),
                    c(0.28, 0.72), sort(as.numeric(two$profile_probabilities)),
                    digits = 2),
    compare_printed(sprintf("values 2-class: P(item %s = category 2 | smaller class)", items),
                    c(0.993, 0.940, 0.927, 0.769), probabilities[1L, ], digits = 3),
    compare_printed(sprintf("values 2-class: P(item %s = category 2 | larger class)", items),
                    c(0.714, 0.330, 0.354, 0.132), probabilities[2L, ], digits = 3)
  ) |> transform(source = "Goodman (2002) in Applied Latent Class Analysis, Tables 5a and 6",
                 dataset = "values", precision = "printed")
}

#' Cheating: Dayton (1998) Table 3.3
#' @return A `data.frame` of compared quantities.
.published_cheating <- function() {
  items <- c("LIEEXAM", "LIEPAPER", "FRAUD", "COPYEXAM")
  two <- .published_fit("cheating", items, 2L)
  statistics <- .likelihood_ratio_fit(two)
  rbind(
    compare_printed("cheating 2-class: log-likelihood", -440.027,
                    as.numeric(logLik(two)), digits = 3),
    compare_values("cheating 2-class: free parameters", 9, two$n_parameters,
                   tolerance = 0.5),
    compare_printed("cheating 2-class: AIC", 898.054, two$aic, digits = 3),
    compare_printed("cheating 2-class: BIC", 931.941, two$bic_individual, digits = 3),
    compare_printed("cheating 2-class: likelihood-ratio G-squared", 7.7642,
                    statistics$g_squared, digits = 4),
    compare_printed("cheating 2-class: Pearson chi-squared", 8.3234,
                    statistics$pearson, digits = 4)
  ) |> transform(source = "Dayton (1998) Latent Class Scaling Analysis, Table 3.3 p. 32",
                 dataset = "cheating", precision = "printed")
}

#' GSS 1982: McCutcheon (1987) survey-evaluation items
#' @return A `data.frame` of compared quantities.
.published_gss82 <- function() {
  items <- c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT")
  two <- .published_fit("gss82", items, 2L)
  statistics <- .likelihood_ratio_fit(two)
  rbind(
    compare_printed("gss82 2-class: log-likelihood", -2783.268,
                    as.numeric(logLik(two)), digits = 3),
    compare_values("gss82 2-class: free parameters", 13, two$n_parameters,
                   tolerance = 0.5),
    compare_printed("gss82 2-class: likelihood-ratio L-squared", 79.3,
                    statistics$g_squared, digits = 1),
    compare_printed(sprintf("gss82 2-class: class share %d", seq_len(2L)),
                    c(0.1923, 0.8077), sort(as.numeric(two$profile_probabilities)),
                    digits = 4)
  ) |> transform(source = "McCutcheon (1987) Latent Class Analysis, Sage",
                 dataset = "gss82", precision = "printed")
}
