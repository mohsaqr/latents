# Functions for compare.R: read a Latent GOLD 6.1 listing of any of this kit's
# models, express each latents fit in Latent GOLD's parameterization with
# delta-method standard errors, and compare estimates, errors, fit statistics
# and every row's posteriors.

# The listing reader, delta method and class matching of the wider kit. This
# file is sourced from the project root (compare.R) or from equivalence/ (the
# testthat file), so the wider kit is found from either.
source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

#' Parse one value of a Latent GOLD fit summary
#' @noRd
.lg_value <- function(lines, label) {
  hit <- grep(sprintf("^%s\t", label), lines, value = TRUE)[1L]
  as.numeric(trimws(strsplit(hit, "\t", fixed = TRUE)[[1L]][2L]))
}

#' Read the fit summary and every parameter table of a Latent GOLD listing
#'
#' Latent GOLD prints one to five term columns, depending on whether a table
#' holds conditioned terms (`| Cluster`), plain regressions, or unconditioned
#' variances. A row holds its term columns and then as many statistic columns
#' as the table's header, so each row's width is its tab count less the
#' header's, plus one. (The intercept `1` is a term, not a number, so the
#' width cannot be found from where the numbers start.) The listing repeats
#' the latent covariances in a last table, which is skipped. Reference
#' categories and fixed loadings (printed with "." as their error) are
#' dropped.
#'
#' @param path A `.lst` file (Latin-1).
#' @return A list: `log_likelihood`, `n_parameters`, `bic`, `n_cases`,
#'   `standard_errors` and `parameters` (`term`, `estimate`, `standard_error`).
read_listing_any <- function(path) {
  lines <- sub("\r$", "", iconv(readLines(path, warn = FALSE), "latin1", "UTF-8"))
  start <- grep("^Regression Parameters", lines)[1L]
  stops <- grep("^(Paired Comparisons|Profile|ProbMeans|Bivariate Residuals|Estimated Values|Iteration Detail)",
                lines)
  end <- c(stops[stops > start], length(lines) + 1L)[[1L]] - 1L
  block <- lines[seq.int(start, end)]
  tabs <- function(x) lengths(regmatches(x, gregexpr("\t", x, fixed = TRUE)))
  headers <- grep("^term\t", block)
  titles <- trimws(sub("\t.*$", "", block[headers - 1L]))
  bounds <- c(headers, length(block) + 2L)
  tables <- lapply(seq_along(headers), function(i) {
    if (identical(titles[[i]], "Variances / Covariances continuous latent")) return(NULL)
    rows <- block[seq.int(headers[[i]] + 1L, bounds[[i + 1L]] - 2L)]
    rows <- rows[nzchar(trimws(sub("\t.*$", "", rows)))]
    if (!length(rows)) return(NULL)
    width <- tabs(rows) - tabs(block[[headers[[i]]]]) + 1L
    fields <- strsplit(rows, "\t", fixed = TRUE)
    value <- function(f, at) {
      x <- if (length(f) >= at) trimws(f[[at]]) else ""
      if (grepl("^-?([0-9]+[.]?[0-9]*|[.][0-9]+)(e[-+]?[0-9]+)?$", x)) as.numeric(x) else NA_real_
    }
    data.frame(
      term = vapply(seq_along(fields), function(r) paste(gsub(" ", "", fields[[r]][seq_len(width[[r]])],
                                                              fixed = TRUE), collapse = ""), character(1)),
      estimate = vapply(seq_along(fields), function(r) value(fields[[r]], width[[r]] + 1L), numeric(1)),
      standard_error = vapply(seq_along(fields), function(r) value(fields[[r]], width[[r]] + 2L),
                              numeric(1)))
  })
  parameters <- do.call(rbind, tables)
  stopifnot("every parameter row has a numeric estimate" = !anyNA(parameters$estimate))
  parameters <- parameters[!is.na(parameters$standard_error) &
                             !duplicated(parameters$term), , drop = FALSE]
  row.names(parameters) <- NULL
  se_line <- grep("^standard errors\t", lines, value = TRUE)[1L]
  list(log_likelihood = .lg_value(lines, "Log-likelihood \\(LL\\)"),
       n_parameters = as.integer(.lg_value(lines, "Number of parameters \\(Npar\\)")),
       bic = .lg_value(lines, "BIC \\(based on LL\\)"),
       n_cases = .lg_value(lines, "Number of cases"),
       standard_errors = trimws(strsplit(se_line, "\t", fixed = TRUE)[[1L]][2L]),
       parameters = parameters)
}

#' Membership logits against Latent GOLD's first class
#' @param gamma Function of a latents class returning its logit (0 for the
#'   reference), for the term `term`.
#' @noRd
.class_logit_terms <- function(perm, gamma, latent, term = "1") {
  later <- seq_along(perm)[-1L]
  stats::setNames(vapply(later, function(j) gamma(perm[j]) - gamma(perm[1L]), numeric(1)),
                  sprintf("%s(%d)<-%s", latent, later, term))
}

#' Latent GOLD terms of a categorical latent class model
#'
#' Item terms are log odds of each category against the first, per class.
#' @param levels Named list: each item's category labels in Latent GOLD order.
#' @param covariate `NULL`, or the name of a membership covariate.
#' @param membership Whether to include the class logits (not for the
#'   measurement part of a transition model).
lca_map <- function(perm, levels, covariate = NULL, latent = "Cluster", membership = TRUE) {
  function(theta) {
    probability <- function(k, v, c) theta[[sprintf("measurement.response.profile_%d.%s:%s", k, v, c)]]
    item_terms <- unlist(lapply(names(levels), function(v) {
      categories <- levels[[v]]
      unlist(lapply(seq_along(perm), function(j) {
        stats::setNames(vapply(categories[-1L], function(c)
          log(probability(perm[j], v, c) / probability(perm[j], v, categories[[1L]])), numeric(1)),
          sprintf("%s(%s)<-1|%s(%d)", v, categories[-1L], latent, j))
      }))
    }))
    logits <- if (!membership) NULL else if (is.null(covariate)) {
      share <- function(k) theta[[sprintf("profile.probability.profile_%d.group_class_1", k)]]
      .class_logit_terms(perm, function(k) log(share(k)), latent)
    } else {
      coefficient <- function(k, term) {
        name <- sprintf("profile.coefficient.profile_%d.%s", k, term)
        if (name %in% names(theta)) theta[[name]] else 0
      }
      c(.class_logit_terms(perm, function(k) coefficient(k, "(Intercept)"), latent),
        .class_logit_terms(perm, function(k) coefficient(k, covariate), latent, covariate))
    }
    c(logits, item_terms)
  }
}

#' Latent GOLD terms of a single-level profile model under one structure
#' @param model `"VVI"`, `"EEI"` or `"VVV"`.
lpa_map <- function(perm, indicators, model) {
  function(theta) {
    mean <- function(k, v) theta[[sprintf("measurement.mean.profile_%d.%s", k, v)]]
    share <- function(k) theta[[sprintf("profile.probability.profile_%d.group_class_1", k)]]
    later <- seq_along(perm)[-1L]
    means <- unlist(lapply(indicators, function(v) c(
      stats::setNames(mean(perm[1L], v), sprintf("%s<-1", v)),
      stats::setNames(vapply(later, function(j) mean(perm[j], v) - mean(perm[1L], v), numeric(1)),
                      sprintf("%s<-Cluster(%d)", v, later)))))
    variances <- switch(model,
      EEI = stats::setNames(vapply(indicators, function(v)
        theta[[sprintf("measurement.variance.shared.%s", v)]], numeric(1)), indicators),
      VVI = unlist(lapply(indicators, function(v) stats::setNames(
        vapply(seq_along(perm), function(j)
          theta[[sprintf("measurement.variance.profile_%d.%s", perm[j], v)]], numeric(1)),
        sprintf("%s|Cluster(%d)", v, seq_along(perm))))),
      VVV = {
        covariance <- function(k, a, b) theta[[sprintf("measurement.covariance.profile_%d.%s:%s", k, a, b)]]
        c(unlist(lapply(indicators, function(v) stats::setNames(
            vapply(seq_along(perm), function(j) covariance(perm[j], v, v), numeric(1)),
            sprintf("%s|Cluster(%d)", v, seq_along(perm))))),
          stats::setNames(vapply(seq_along(perm), function(j)
            covariance(perm[j], indicators[[2L]], indicators[[1L]]), numeric(1)),
            sprintf("%s<->%s|Cluster(%d)", indicators[[1L]], indicators[[2L]], seq_along(perm))))
      })
    c(.class_logit_terms(perm, function(k) log(share(k)), "Cluster"), means, variances)
  }
}

#' Latent GOLD terms of a latent class growth model (long format)
lcga_map <- function(perm) {
  function(theta) {
    gamma <- function(k) if (k == 1L) 0 else theta[[sprintf("membership.class_%d.(Intercept)", k)]]
    classes <- seq_along(perm)
    c(.class_logit_terms(perm, gamma, "Class"),
      stats::setNames(vapply(classes, function(j)
        theta[[sprintf("coefficient.class_%d.(Intercept)", perm[j])]], numeric(1)),
        sprintf("y<-1|Class(%d)", classes)),
      stats::setNames(vapply(classes, function(j)
        theta[[sprintf("coefficient.class_%d.time", perm[j])]], numeric(1)),
        sprintf("y<-time|Class(%d)", classes)),
      stats::setNames(vapply(classes, function(j)
        exp(2 * theta[[sprintf("log_sigma.class_%d", perm[j])]]), numeric(1)),
        sprintf("y|Class(%d)", classes)))
  }
}

#' Latent GOLD terms of a growth mixture written as a factor model (wide)
#'
#' F1 and F2 are the random intercept and slope: their class means are the
#' latents coefficients, their (co)variances the Cholesky products, and every
#' wave's residual variance the class's exp(2 log_sigma).
#' @param random `"slope"` (class-specific 2 x 2 covariance) or `"intercept"`
#'   (one random-intercept variance shared by the classes).
#' @param waves Number of waves.
gmm_map <- function(perm, random, waves) {
  function(theta) {
    gamma <- function(k) if (k == 1L) 0 else theta[[sprintf("membership.class_%d.(Intercept)", k)]]
    classes <- seq_along(perm)
    by_class <- function(f, name) stats::setNames(vapply(classes, function(j) f(perm[j]), numeric(1)),
                                                  sprintf(name, classes))
    covariance <- if (identical(random, "slope")) {
      a <- function(k) theta[[sprintf("random.log_cholesky.class_%d.(Intercept):(Intercept)", k)]]
      b <- function(k) theta[[sprintf("random.log_cholesky.class_%d.time:time", k)]]
      lower <- function(k) theta[[sprintf("random.cholesky.class_%d.time:(Intercept)", k)]]
      c(by_class(function(k) exp(2 * a(k)), "F1|Class(%d)"),
        by_class(function(k) lower(k)^2 + exp(2 * b(k)), "F2|Class(%d)"),
        by_class(function(k) exp(a(k)) * lower(k), "F1<->F2|Class(%d)"))
    } else {
      c(F1 = exp(2 * theta[["random.log_cholesky.shared.(Intercept):(Intercept)"]]))
    }
    residuals <- unlist(lapply(seq_len(waves), function(t)
      by_class(function(k) exp(2 * theta[[sprintf("log_sigma.class_%d", k)]]),
               paste0("y", t, "|Class(%d)"))))
    c(.class_logit_terms(perm, gamma, "Class"),
      by_class(function(k) theta[[sprintf("coefficient.class_%d.(Intercept)", k)]], "F1<-1|Class(%d)"),
      by_class(function(k) theta[[sprintf("coefficient.class_%d.time", k)]], "F2<-1|Class(%d)"),
      covariance, residuals)
  }
}

#' Latent GOLD terms of a two-state latent transition model without covariates
#'
#' Initial and transition terms are log odds of state 2 against state 1 (in
#' Latent GOLD's labels), given the previous state. latents reports them as
#' probabilities, named `transition_probability.<to>.group_class_1:<from>`.
#' @param categorical Whether the indicators are binary (log-odds terms) or
#'   Gaussian (intercept, contrast and state variances).
lta_map <- function(perm, indicators, categorical) {
  function(theta) {
    initial <- function(k) theta[[sprintf("profile.initial_probability.profile_%d.group_class_1", k)]]
    move <- function(from, to) {
      theta[[sprintf("profile.transition_probability.profile_%d.group_class_1:profile_%d", to, from)]]
    }
    dynamics <- c(
      `State[=0](2)<-1` = log(initial(perm[2L]) / initial(perm[1L])),
      `State(2)<-1|State[-1](1)` = log(move(perm[1L], perm[2L]) / move(perm[1L], perm[1L])),
      `State(2)<-1|State[-1](2)` = log(move(perm[2L], perm[2L]) / move(perm[2L], perm[1L])))
    measurement <- if (categorical) {
      lca_map(perm, stats::setNames(rep(list(c("0", "1")), length(indicators)), indicators),
              latent = "State", membership = FALSE)(theta)
    } else {
      mean <- function(k, v) theta[[sprintf("measurement.mean.profile_%d.%s", k, v)]]
      variance <- function(k, v) theta[[sprintf("measurement.variance.profile_%d.%s", k, v)]]
      measurement_terms(mean, variance, perm, indicators, "State")
    }
    c(dynamics, measurement)
  }
}

#' Latent GOLD terms of a two-state transition model with a transition covariate
lta_covariate_map <- function(perm, indicators) {
  function(theta) {
    move <- function(from, to, term) {
      if (from == to) 0 else
        theta[[sprintf("transition.group_class_1.profile_%d->profile_%d.%s", from, to, term)]]
    }
    initial <- c(theta[["initial.group_class_1.profile_1.(Intercept)"]], 0)
    logit <- function(previous, term) move(perm[previous], perm[2L], term) -
      move(perm[previous], perm[1L], term)
    mean <- function(k, v) theta[[sprintf("mean.all.profile_%d.%s", k, v)]]
    variance <- function(k, v) exp(theta[[sprintf("log_variance.all.profile_%d.%s", k, v)]])
    c(`State[=0](2)<-1` = initial[perm[2L]] - initial[perm[1L]],
      `State(2)<-1|State[-1](1)` = logit(1L, "(Intercept)"),
      `State(2)<-1|State[-1](2)` = logit(2L, "(Intercept)"),
      `State(2)<-z|State[-1](1)` = logit(1L, "z"),
      `State(2)<-z|State[-1](2)` = logit(2L, "z"),
      measurement_terms(mean, variance, perm, indicators, "State"))
  }
}

#' Latent GOLD terms of a two-class mixture regression of a count (NB2) or an
#' ordinal outcome (cumulative logit) on x
#'
#' Latent GOLD's `poisson overdispersed` reports the NB2 dispersion (variance
#' mu + d mu^2), latents' exp(log_dispersion). Its `cumlogit` intercepts model
#' "category c or higher", so each is minus latents' threshold below c; the
#' thresholds are the first plus the exponentiated gaps.
#' @param family `"negative_binomial"` or `"ordinal"`.
#' @param categories Number of ordinal categories.
family_regression_map <- function(perm, family, categories = 4L) {
  function(theta) {
    gamma <- function(k) if (k == 1L) 0 else theta[[sprintf("membership.class_%d.(Intercept)", k)]]
    classes <- seq_along(perm)
    slopes <- stats::setNames(vapply(classes, function(j)
      theta[[sprintf("coefficient.class_%d.x", perm[j])]], numeric(1)), sprintf("y<-x|Cluster(%d)", classes))
    if (identical(family, "negative_binomial")) {
      return(c(.class_logit_terms(perm, gamma, "Cluster"),
               stats::setNames(vapply(classes, function(j)
                 theta[[sprintf("coefficient.class_%d.(Intercept)", perm[j])]], numeric(1)),
                 sprintf("y<-1|Cluster(%d)", classes)),
               slopes,
               stats::setNames(vapply(classes, function(j)
                 exp(theta[[sprintf("log_dispersion.class_%d", perm[j])]]), numeric(1)),
                 sprintf("y|Cluster(%d)", classes))))
    }
    thresholds <- function(k) {
      first <- theta[[sprintf("threshold.class_%d.1|2", k)]]
      gaps <- vapply(seq_len(categories - 2L) + 1L, function(c)
        exp(theta[[sprintf("log_threshold_gap.class_%d.%d|%d", k, c, c + 1L)]]), numeric(1))
      cumsum(c(first, gaps))
    }
    intercepts <- unlist(lapply(classes, function(j) stats::setNames(
      -thresholds(perm[j]), sprintf("y(%d)<-1|Cluster(%d)", seq_len(categories - 1L) + 1L, j))))
    c(.class_logit_terms(perm, gamma, "Cluster"), intercepts, slopes)
  }
}

# Tight convergence and many starts on both sides, so gaps measure Latent
# GOLD's printed rounding, not where two searches stopped.
fit_case <- function(case) {
  data <- case$data
  starts <- list(n_starts = 50L, seed = 1L, tol = 1e-12)
  switch(case$kind,
    lcga = do.call(mixture_regression, c(list(y ~ time, data, 2L, id = "id",
                                              class_level = "group", max_iter = 100000L), starts)),
    gmm_slope = do.call(mixture_regression, c(list(y ~ time, data, 2L, id = "id",
                                                   class_level = "group", random = ~time,
                                                   random_covariance = "varying",
                                                   max_iter = 100000L), starts)),
    gmm_intercept = do.call(mixture_regression, c(list(y ~ time, data, 2L, id = "id",
                                                       class_level = "group", random = ~1,
                                                       random_covariance = "equal",
                                                       max_iter = 100000L), starts)),
    lca = , lca_enumeration = do.call(lca, c(list(data, case$items, n_classes = case$classes,
                                                  max_iter = 20000L), starts)),
    # The membership logits are flat enough that 1e-12 stops 1e-5 short in one
    # item logit (checked: 1e-15 moves it within Latent GOLD's rounding).
    lca_covariate = do.call(lca, c(list(data, case$items, n_classes = case$classes,
                                        profile_covariates = "x", max_iter = 200000L),
                                   utils::modifyList(starts, list(tol = 1e-15)))),
    lpa = , lpa_enumeration = do.call(lpa, c(list(data, setdiff(names(data), "id"),
                                                  n_profiles = case$classes, model = case$model,
                                                  max_iter = 20000L), starts)),
    lta = do.call(lta, c(list(data, case$items, "id", 2L, time = "time"),
                         if (case$categorical) list(categorical = case$items), starts)),
    lta_covariate = do.call(lta, c(list(data, case$items, "id", 2L, time = "time",
                                        transition_covariates = "z"), starts)),
    nb_regression = do.call(mixture_regression, c(list(y ~ x, data, 2L, family = "negative_binomial",
                                                       max_iter = 100000L),
                                                  utils::modifyList(starts, list(tol = 1e-14)))),
    # EM moves slowly along the ordinal thresholds: at 1e-14 one threshold sat
    # 9e-6 short of the maximum (outside Latent GOLD's rounding); 1e-16 reaches it.
    ordinal_regression = do.call(mixture_regression, c(list(y ~ x, data, 2L, family = "ordinal",
                                                            max_iter = 1000000L),
                                                       utils::modifyList(starts, list(tol = 1e-16)))))
}

.case_levels <- function(case) {
  stats::setNames(lapply(case$items, function(v) as.character(sort(unique(case$data[[v]])))),
                  case$items)
}

maps_for <- function(case) {
  k <- case$classes %||% 2L
  perms <- permutations(k)
  switch(case$kind,
    lcga = lapply(perms, lcga_map),
    gmm_slope = lapply(perms, gmm_map, random = "slope", waves = 5L),
    gmm_intercept = lapply(perms, gmm_map, random = "intercept", waves = 5L),
    lca = lapply(perms, lca_map, levels = .case_levels(case)),
    lca_covariate = lapply(perms, lca_map, levels = .case_levels(case), covariate = "x"),
    lpa = lapply(perms, lpa_map, indicators = c("y1", "y2", "y3")[seq_len(
      ncol(case$data) - 1L)], model = case$model),
    lta = lapply(perms, lta_map, indicators = case$items, categorical = case$categorical),
    lta_covariate = lapply(perms, lta_covariate_map, indicators = case$items),
    nb_regression = lapply(perms, family_regression_map, family = "negative_binomial"),
    ordinal_regression = lapply(perms, family_regression_map, family = "ordinal"))
}

#' Posterior class probabilities of a latents fit, one row per Latent GOLD row
#' @noRd
.latents_posteriors <- function(fit, case) {
  if (case$kind %in% c("lcga", "gmm_slope", "gmm_intercept", "nb_regression", "ordinal_regression")) {
    return(fit$expectation$tau)
  }
  fit$subject_posteriors
}

#' Compare every row's posteriors under the class matching chosen for the terms
#' @noRd
compare_posteriors <- function(fit, case, name, returned, perm) {
  path <- file.path(returned, paste0(name, "_posteriors.txt"))
  lg <- utils::read.delim(path, check.names = FALSE)
  columns <- grep("#[0-9]+$", names(lg), value = TRUE)
  ours <- .latents_posteriors(fit, case)
  stopifnot("one latents posterior row per Latent GOLD row" = nrow(ours) == nrow(lg),
            "the same number of classes" = ncol(ours) == length(columns))
  gap <- max(abs(as.matrix(lg[columns]) - ours[, perm, drop = FALSE]))
  data.frame(case = name, quantity = sprintf("posteriors (%d rows x %d classes)", nrow(lg),
                                             length(columns)),
             latent_gold = NA_real_, latents = NA_real_, difference = gap, tolerance = 2e-5,
             agrees = gap <= 2e-5, n_values = nrow(lg) * length(columns))
}

#' Fit statistics of an enumeration step: likelihood, parameters, BIC
#' @noRd
compare_fit_statistics <- function(fit, listing, name) {
  ours <- c(as.numeric(stats::logLik(fit)), fit$n_parameters, stats::BIC(fit))
  rows <- data.frame(case = name, quantity = c("log likelihood", "number of parameters", "BIC"),
                     latent_gold = c(listing$log_likelihood, listing$n_parameters, listing$bic),
                     latents = ours, tolerance = c(6e-5, 0, 1.2e-4))
  rows$difference <- abs(rows$latents - rows$latent_gold)
  rows$agrees <- rows$difference <= rows$tolerance
  rows$n_values <- 1L
  rows[c("case", "quantity", "latent_gold", "latents", "difference", "tolerance", "agrees",
         "n_values")]
}

#' A latents start at Latent GOLD's printed solution of an enumeration step
#'
#' Beyond the true number of classes a mixture likelihood has many local
#' maxima (and, for VVI, unbounded spikes at the variance floor), so the
#' enumeration is compared at Latent GOLD's solution: latents is started there
#' and must reach the same stationary point, likelihood, parameters and BIC.
#' @noRd
.start_from_listing <- function(listing, case) {
  value <- function(term, default = NA_real_) {
    at <- match(term, listing$parameters$term)
    if (is.na(at)) default else listing$parameters$estimate[[at]]
  }
  k <- case$classes
  logits <- vapply(seq_len(k), function(j) if (j == 1L) 0 else value(sprintf("Cluster(%d)<-1", j)),
                   numeric(1))
  shares <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
  given <- function(j) if (k == 1L) "" else sprintf("|Cluster(%d)", j)
  if (identical(case$kind, "lpa_enumeration")) {
    vars <- c("y1", "y2", "y3")
    means <- t(vapply(seq_len(k), function(j) vapply(vars, function(v)
      value(sprintf("%s<-1", v)) + if (j == 1L) 0 else value(sprintf("%s<-Cluster(%d)", v, j)),
      numeric(1)), numeric(length(vars))))
    variances <- t(vapply(seq_len(k), function(j) vapply(vars, function(v)
      value(sprintf("%s|Cluster(%d)", v, j)), numeric(1)), numeric(length(vars))))
    return(list(means = matrix(means, k), variances = matrix(variances, k),
                profile_probabilities = matrix(shares, 1L), group_probabilities = 1))
  }
  responses <- lapply(stats::setNames(case$items, case$items), function(v) {
    odds <- vapply(seq_len(k), function(j) exp(value(sprintf("%s(1)<-1%s", v, given(j)))),
                   numeric(1))
    matrix(c(rep(1, k), odds) / rep(1 + odds, 2L), k, 2L, dimnames = list(NULL, c("0", "1")))
  })
  list(means = matrix(numeric(), k, 0L), variances = matrix(numeric(), k, 0L),
       profile_probabilities = matrix(shares, 1L), group_probabilities = 1,
       response_probabilities = responses)
}

#' Refit an enumeration step from Latent GOLD's solution
#' @noRd
fit_at_reference <- function(case, listing) {
  start <- .start_from_listing(listing, case)
  arguments <- list(n_starts = 1L, seed = 1L, tol = 1e-12, max_iter = 20000L, start = start)
  if (identical(case$kind, "lpa_enumeration")) {
    do.call(lpa, c(list(case$data, c("y1", "y2", "y3"), n_profiles = case$classes,
                        model = case$model), arguments))
  } else {
    do.call(lca, c(list(case$data, case$items, n_classes = case$classes), arguments))
  }
}

.case_function <- c(lcga = "growth", gmm_slope = "growth", gmm_intercept = "growth",
                    lca = "lca", lca_covariate = "lca", lca_enumeration = "enumeration",
                    lpa = "lpa", lpa_enumeration = "enumeration", lta = "lta",
                    lta_covariate = "lta", nb_regression = "mixture_regression",
                    ordinal_regression = "mixture_regression")

#' Fit and compare every case
#' @return One row per compared quantity, with the function it checks.
compare_all <- function(cases, returned) {
  search <- data.frame(case = character(), latents_search = numeric(), latent_gold = numeric())
  rows <- lapply(names(cases), function(name) {
    case <- cases[[name]]
    fit <- fit_case(case)
    listing <- read_listing_any(file.path(returned, paste0(name, ".lst")))
    out <- if (case$kind %in% c("lca_enumeration", "lpa_enumeration")) {
      search <<- rbind(search, data.frame(case = name, latents_search = as.numeric(stats::logLik(fit)),
                                          latent_gold = listing$log_likelihood))
      compare_fit_statistics(fit_at_reference(case, listing), listing, name)
    } else {
      theta <- stats::coef(fit)
      covariance <- stats::vcov(fit)
      stopifnot("coef() and vcov() name the same estimates" =
                  identical(names(theta), rownames(covariance)))
      maps <- maps_for(case)
      terms <- compare_case(listing, maps, theta, covariance,
                            as.numeric(stats::logLik(fit)), fit$n_parameters, name)
      perm_index <- which.min(vapply(maps, function(map) {
        value <- map(theta)
        at <- match(listing$parameters$term, names(value))
        if (anyNA(at)) Inf else max(abs(value[at] - listing$parameters$estimate))
      }, numeric(1)))
      rbind(transform(terms, n_values = 1L),
            compare_posteriors(fit, case, name, returned,
                               permutations(case$classes %||% 2L)[[perm_index]]))
    }
    out$fn <- .case_function[[case$kind]]
    out
  })
  out <- do.call(rbind, rows)
  attr(out, "search") <- transform(search, reached = latents_search >= latent_gold - 1e-4)
  out
}
