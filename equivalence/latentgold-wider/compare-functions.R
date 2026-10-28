# Functions for compare.R: read a Latent GOLD 6.1 listing, and express a
# latents fit in Latent GOLD's parameterization with delta-method standard
# errors, so both are compared term by term.

#' Read the fit summary and parameter tables of a Latent GOLD listing
#'
#' @param path A `.lst` file (Latin-1).
#' @return A list: `log_likelihood`, `n_parameters`, `standard_errors` (the
#'   type the listing records) and `parameters`, a data frame with one row per
#'   estimated (not reference-fixed) term: `term` (the term's fields joined,
#'   e.g. `"y<-x|Cluster(2)"`), `estimate` and `standard_error`.
read_listing <- function(path) {
  lines <- iconv(readLines(path, warn = FALSE), "latin1", "UTF-8")
  value_after <- function(label) {
    hit <- grep(sprintf("^%s\t", label), lines, value = TRUE)[1L]
    trimws(strsplit(hit, "\t", fixed = TRUE)[[1L]][2L])
  }
  start <- grep("^Regression Parameters", lines)[1L]
  end <- grep("^Paired Comparisons", lines)[1L]
  block <- lines[seq.int(start + 1L, end - 1L)]
  fields <- strsplit(block, "\t", fixed = TRUE)
  rows <- Filter(function(f) length(f) >= 7L && !f[[1L]] %in% c("term", "Variances"), fields)
  parameters <- data.frame(
    term = vapply(rows, function(f) paste(gsub(" ", "", f[1:5], fixed = TRUE), collapse = ""),
                  character(1)),
    estimate = vapply(rows, function(f) as.numeric(f[[6L]]), numeric(1)),
    standard_error = vapply(rows, function(f) suppressWarnings(as.numeric(f[[7L]])), numeric(1)))
  # A reference category is printed as 0.0000 with "." for its error.
  parameters <- subset(parameters, !is.na(standard_error))
  rownames(parameters) <- NULL
  list(log_likelihood = as.numeric(value_after("Log-likelihood \\(LL\\)")),
       n_parameters = as.integer(value_after("Number of parameters \\(Npar\\)")),
       standard_errors = value_after("standard errors"),
       parameters = parameters)
}

#' Delta-method standard errors of a function of the estimates
#'
#' @param map Function of the named estimate vector returning a named vector.
#' @param theta,vcov Estimates and their covariance (matching names).
#' @param step Central-difference step.
#' @return A data frame: `term`, `estimate`, `standard_error`.
delta_terms <- function(map, theta, vcov, step = 1e-6) {
  value <- map(theta)
  jacobian <- vapply(seq_along(theta), function(i) {
    up <- theta
    down <- theta
    up[i] <- up[i] + step
    down[i] <- down[i] - step
    (map(up) - map(down)) / (2 * step)
  }, numeric(length(value)))
  jacobian <- matrix(jacobian, nrow = length(value))
  variance <- rowSums((jacobian %*% vcov) * jacobian)
  data.frame(term = names(value), estimate = unname(value),
             standard_error = sqrt(pmax(variance, 0)))
}

#' Latent GOLD terms of a mixture regression
#' @param perm `perm[j]` is the latents class shown as Latent GOLD's class j.
#' @param family `"gaussian"`, `"poisson"` or `"binomial"`.
#' @return A function of the estimate vector.
regression_map <- function(perm, family) {
  outcome <- if (identical(family, "binomial")) "y(1)" else "y"
  function(theta) {
    beta <- function(k, term) theta[[sprintf("coefficient.class_%d.%s", k, term)]]
    gamma <- function(k, term) {
      if (k == 1L) 0 else theta[[sprintf("membership.class_%d.%s", k, term)]]
    }
    classes <- seq_along(perm)
    later <- classes[-1L]
    c(stats::setNames(vapply(later, function(j) gamma(perm[j], "(Intercept)") -
                               gamma(perm[1L], "(Intercept)"), numeric(1)),
                      sprintf("Cluster(%d)<-1", later)),
      stats::setNames(vapply(later, function(j) gamma(perm[j], "z") - gamma(perm[1L], "z"),
                             numeric(1)), sprintf("Cluster(%d)<-z", later)),
      stats::setNames(vapply(classes, function(j) beta(perm[j], "(Intercept)"), numeric(1)),
                      sprintf("%s<-1|Cluster(%d)", outcome, classes)),
      stats::setNames(vapply(classes, function(j) beta(perm[j], "x"), numeric(1)),
                      sprintf("%s<-x|Cluster(%d)", outcome, classes)),
      if (identical(family, "gaussian")) {
        stats::setNames(vapply(classes, function(j)
          exp(2 * theta[[sprintf("log_sigma.class_%d", perm[j])]]), numeric(1)),
          sprintf("y|Cluster(%d)", classes))
      })
  }
}

#' Measurement terms shared by the profile, two-level and transition maps
#' @param mean,variance Functions `(class, indicator)` of the latents values.
#' @param perm Class order, as for [regression_map()].
#' @param latent Latent GOLD's name of the variable (`"Cluster"`, `"State"`).
measurement_terms <- function(mean, variance, perm, indicators, latent) {
  classes <- seq_along(perm)
  later <- classes[-1L]
  unlist(lapply(indicators, function(v) c(
    stats::setNames(mean(perm[1L], v), sprintf("%s<-1", v)),
    stats::setNames(vapply(later, function(j) mean(perm[j], v) - mean(perm[1L], v),
                           numeric(1)), sprintf("%s<-%s(%d)", v, latent, later)),
    stats::setNames(vapply(classes, function(j) variance(perm[j], v), numeric(1)),
                    sprintf("%s|%s(%d)", v, latent, classes)))))
}

#' Latent GOLD terms of a single-level profile model (natural-scale estimates)
profile_map <- function(perm, indicators) {
  function(theta) {
    mean <- function(k, v) theta[[sprintf("measurement.mean.profile_%d.%s", k, v)]]
    variance <- function(k, v) theta[[sprintf("measurement.variance.profile_%d.%s", k, v)]]
    share <- function(k) theta[[sprintf("profile.probability.profile_%d.group_class_1", k)]]
    later <- seq_along(perm)[-1L]
    c(stats::setNames(vapply(later, function(j) log(share(perm[j]) / share(perm[1L])),
                             numeric(1)), sprintf("Cluster(%d)<-1", later)),
      measurement_terms(mean, variance, perm, indicators, "Cluster"))
  }
}

#' Latent GOLD terms of a two-level profile model with two group classes
#' @param group_perm Group-class order, as `perm` is for the profiles.
two_level_map <- function(perm, group_perm, indicators) {
  function(theta) {
    mean <- function(k, v) theta[[sprintf("measurement.mean.profile_%d.%s", k, v)]]
    variance <- function(k, v) theta[[sprintf("measurement.variance.profile_%d.%s", k, v)]]
    within <- function(k, m) theta[[sprintf("profile.probability.profile_%d.group_class_%d", k, m)]]
    group <- function(m) theta[[sprintf("group.probability.group_class_%d", m)]]
    logit <- function(m) log(within(perm[2L], m) / within(perm[1L], m))
    c(`GClass(2)<-1` = log(group(group_perm[2L]) / group(group_perm[1L])),
      `Cluster(2)<-1` = logit(group_perm[1L]),
      `Cluster(2)<-GClass(2)` = logit(group_perm[2L]) - logit(group_perm[1L]),
      measurement_terms(mean, variance, perm, indicators, "Cluster"))
  }
}

#' Latent GOLD terms of a two-state transition model (lta() estimates)
transition_map <- function(perm, indicators) {
  function(theta) {
    mean <- function(k, v) theta[[sprintf("mean.all.profile_%d.%s", k, v)]]
    variance <- function(k, v) exp(theta[[sprintf("log_variance.all.profile_%d.%s", k, v)]])
    initial <- c(theta[["initial.group_class_1.profile_1.(Intercept)"]], 0)
    move <- function(from, to) {
      if (from == to) 0 else
        theta[[sprintf("transition.group_class_1.profile_%d->profile_%d.(Intercept)", from, to)]]
    }
    c(`State[=0](2)<-1` = initial[perm[2L]] - initial[perm[1L]],
      `State(2)<-1|State[-1](1)` = move(perm[1L], perm[2L]) - move(perm[1L], perm[1L]),
      `State(2)<-1|State[-1](2)` = move(perm[2L], perm[2L]) - move(perm[2L], perm[1L]),
      measurement_terms(mean, variance, perm, indicators, "State"))
  }
}

#' Compare a latents fit with a listing under the best class matching
#'
#' @param listing From [read_listing()].
#' @param maps A list of candidate maps, one per class matching.
#' @param theta,vcov The latents estimates and covariance.
#' @param log_likelihood,n_parameters The latents fit's.
#' @param case Case label.
#' @param tolerance Half a unit of Latent GOLD's last printed digit plus the
#'   two programs' convergence slack.
#' @return One row per compared quantity: `case`, `quantity`, `latent_gold`,
#'   `latents`, `difference`, `tolerance`, `agrees`.
compare_case <- function(listing, maps, theta, vcov, log_likelihood, n_parameters, case,
                         tolerance = 5e-5 + 2e-6) {
  candidates <- lapply(maps, delta_terms, theta = theta, vcov = vcov)
  gap <- vapply(candidates, function(terms) {
    at <- match(listing$parameters$term, terms$term)
    if (anyNA(at)) Inf else max(abs(terms$estimate[at] - listing$parameters$estimate))
  }, numeric(1))
  terms <- candidates[[which.min(gap)]]
  at <- match(listing$parameters$term, terms$term)
  if (anyNA(at)) {
    stop(sprintf("%s: latents has no term for %s", case,
                 paste(listing$parameters$term[is.na(at)], collapse = ", ")))
  }
  rows <- rbind(
    data.frame(quantity = c("log likelihood", "number of parameters"),
               latent_gold = c(listing$log_likelihood, listing$n_parameters),
               latents = c(log_likelihood, n_parameters),
               tolerance = c(5e-5 + 1e-5, 0)),
    data.frame(quantity = paste("estimate", listing$parameters$term),
               latent_gold = listing$parameters$estimate, latents = terms$estimate[at],
               tolerance = tolerance),
    data.frame(quantity = paste(listing$standard_errors, "s.e.", listing$parameters$term),
               latent_gold = listing$parameters$standard_error,
               latents = terms$standard_error[at], tolerance = tolerance))
  rows$difference <- abs(rows$latents - rows$latent_gold)
  rows$agrees <- rows$difference <= rows$tolerance
  cbind(case = case, rows[c("quantity", "latent_gold", "latents", "difference", "tolerance",
                             "agrees")])
}

#' All orderings of `k` classes
permutations <- function(k) {
  if (k == 1L) return(list(1L))
  unlist(lapply(seq_len(k), function(first) lapply(permutations(k - 1L), function(rest)
    c(first, setdiff(seq_len(k), first)[rest]))), recursive = FALSE)
}

# Tight convergence on both sides, so the gaps measure Latent GOLD's printed
# rounding, not where two searches stopped. Mixture regressions run to
# tol = 1e-15: at 1e-12 the flat membership intercept of r03 (s.e. 0.69) was
# still 1.8e-4 short of its maximum (log likelihood 4e-8 lower).
fit_case <- function(case) {
  data <- case$data
  switch(case$kind,
    regression = {
      binomial <- identical(case$family, "binomial")
      mixture_regression(y ~ x, data, 2L, family = case$family, membership = ~z,
                         id = if (binomial) "id", class_level = if (binomial) "group" else "observation",
                         weights = if (identical(case$se, "weighted")) "w",
                         vcov_type = if (identical(case$se, "standard")) "observed" else "robust",
                         n_starts = 30L, seed = 1L, tol = 1e-15, max_iter = 100000L)
    },
    profiles = lpa(data, c("y1", "y2", "y3"), n_profiles = 3L, weights = "w",
                   n_starts = 30L, seed = 1L, tol = 1e-12, max_iter = 20000L),
    transitions = lta(data, c("y1", "y2"), id = "id", n_profiles = 2L, time = "time",
                      weights = "w", n_starts = 30L, seed = 1L, tol = 1e-12),
    two_level = multilpa(data, c("y1", "y2"), id = "g", n_profiles = 2L,
                         n_group_classes = 2L, weights = "w", n_starts = 30L, seed = 1L,
                         tol = 1e-12, max_iter = 20000L))
}

maps_for <- function(case) {
  switch(case$kind,
    regression = lapply(permutations(2L), regression_map, family = case$family),
    profiles = lapply(permutations(3L), profile_map, indicators = c("y1", "y2", "y3")),
    transitions = lapply(permutations(2L), transition_map, indicators = c("y1", "y2")),
    two_level = unlist(lapply(permutations(2L), function(perm)
      lapply(permutations(2L), function(group_perm)
        two_level_map(perm, group_perm, c("y1", "y2")))), recursive = FALSE))
}

#' Fit and compare every case
#' @param cases The kit's `cases.rds`.
#' @param returned Folder of the retained listings.
#' @return The rows of [compare_case()] for every case.
compare_all <- function(cases, returned) {
  do.call(rbind, lapply(names(cases), function(name) {
    case <- cases[[name]]
    fit <- fit_case(case)
    listing <- read_listing(file.path(returned, paste0(name, ".lst")))
    theta <- stats::coef(fit)
    covariance <- stats::vcov(fit)
    stopifnot("coef() and vcov() name the same estimates" =
                identical(names(theta), rownames(covariance)))
    compare_case(listing, maps_for(case), theta, covariance,
                 as.numeric(stats::logLik(fit)), fit$n_parameters, name)
  }))
}
