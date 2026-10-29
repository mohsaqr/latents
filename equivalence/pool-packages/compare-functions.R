# Functions for compare.R: pool_imputations() against two implementations of
# Rubin's rules that do not share code with mice (which the harness uses):
# mitml::testEstimates() and mitools::MIcombine(). Each is given latents' own
# per-imputation estimates and variances (get_results(x, "imputations")), so
# only the pooling arithmetic is compared: the pooled estimate, its standard
# error, Rubin's degrees of freedom (complete-data df infinite), the relative
# increase in variance and the fraction of missing information.
#
# Tolerances, fixed before the comparison: 1e-10 of max(1, |value|) for the
# estimate, standard error, relative increase and fraction; 1e-8 of the
# degrees of freedom (which reach 1e28 for parameters the imputations do not
# move).

.pool_case <- function(name, data, vars, classes, profile_covariates = NULL) {
  list(name = name, data = data, vars = vars, classes = classes,
       profile_covariates = profile_covariates)
}

pool_cases <- function() {
  iris_missing <- withr::with_seed(1L, {
    d <- datasets::iris[c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")]
    d$Sepal.Width[sample.int(150L, 15L)] <- NA
    d$Petal.Length[sample.int(150L, 15L)] <- NA
    d
  })
  srl_missing <- withr::with_seed(1L, {
    d <- latents::srl[c("cognitive_strategies", "intrinsic_value", "self_regulation",
                        "self_efficacy", "test_anxiety")]
    invisible(lapply(c("intrinsic_value", "self_regulation", "test_anxiety"), function(v) {
      d[[v]][sample.int(nrow(d), 30L)] <<- NA
    }))
    d
  })
  list(
    .pool_case("iris_k2", iris_missing, names(iris_missing), 2L),
    .pool_case("iris_k3", iris_missing, names(iris_missing), 3L),
    .pool_case("srl_k2", srl_missing, names(srl_missing), 2L),
    .pool_case("srl_covariate_k2", srl_missing, setdiff(names(srl_missing), "test_anxiety"), 2L,
               profile_covariates = "test_anxiety"))
}

fit_pool_case <- function(case) {
  imputed <- mice::mice(case$data, m = 5L, method = "pmm", seed = 1L, printFlag = FALSE)
  arguments <- list(imputed, case$vars, id = NULL, n_profiles = case$classes, n_starts = 10L,
                    seed = 1L, tol = 1e-12)
  if (!is.null(case$profile_covariates)) arguments$profile_covariates <- case$profile_covariates
  do.call(pool_imputations, arguments)
}

#' The per-imputation estimates as an m x p matrix, and their variances
#' @noRd
.per_imputation <- function(pooled) {
  rows <- get_results(pooled, "imputations")
  key <- paste(rows$level, rows$outcome, rows$term, rows$parameter, sep = "|")
  m <- max(rows$imputation)
  keys <- unique(key)
  estimate <- t(vapply(seq_len(m), function(i)
    rows$estimate[rows$imputation == i][match(keys, key[rows$imputation == i])], numeric(length(keys))))
  variance <- t(vapply(seq_len(m), function(i)
    rows$standard_error[rows$imputation == i][match(keys, key[rows$imputation == i])]^2,
    numeric(length(keys))))
  colnames(estimate) <- colnames(variance) <- keys
  list(estimate = estimate, variance = variance)
}

#' One row per quantity the reference reports, compared with latents'
#' @noRd
.quantity_rows <- function(case, reference, ours, theirs) {
  labels <- c(estimate = "pooled estimates", standard_error = "pooled standard errors",
              df = "degrees of freedom", riv = "relative increase in variance",
              fmi = "fraction of missing information")
  tolerance <- c(estimate = 1e-10, standard_error = 1e-10, df = 1e-8, riv = 1e-10, fmi = 1e-10)
  reported <- intersect(names(labels), names(theirs))
  differences <- vapply(reported, function(q)
    max(abs(theirs[[q]] - ours[[q]]) / pmax(1, abs(theirs[[q]]))), numeric(1))
  data.frame(case = case$name, reference = reference,
             quantity = sprintf("%s (%d values)", labels[reported], length(ours$estimate)),
             reference_value = NA, latents = NA, difference = unname(differences),
             tolerance = unname(tolerance[reported]), n_values = length(ours$estimate),
             agrees = unname(differences <= tolerance[reported]))
}

compare_pool_case <- function(case) {
  pooled <- fit_pool_case(case)
  per <- .per_imputation(pooled)
  estimates <- get_results(pooled, "estimates")
  key <- paste(estimates$level, estimates$outcome, estimates$term, estimates$parameter, sep = "|")
  at <- match(colnames(per$estimate), key)
  stopifnot("every imputed quantity is pooled" = !anyNA(at))
  # Parameters every imputation gives the same value (and the same error) have
  # no between-imputation variance; their relative increase is 0/0 in the
  # references' arithmetic, so they are compared where it is defined.
  moved <- apply(per$estimate, 2L, stats::var) > 0
  ours <- list(estimate = estimates$estimate[at], standard_error = estimates$standard_error[at],
               df = estimates$df[at], riv = estimates$riv[at], fmi = estimates$fmi[at])
  ours <- lapply(ours, function(v) v[moved])
  m <- nrow(per$estimate)
  mitml_result <- mitml::testEstimates(qhat = lapply(seq_len(m), function(i) per$estimate[i, moved]),
                                       uhat = lapply(seq_len(m), function(i) per$variance[i, moved]))
  table <- mitml_result$estimates
  mitml_values <- list(estimate = table[, "Estimate"], standard_error = table[, "Std.Error"],
                       df = table[, "df"], riv = table[, "RIV"], fmi = table[, "FMI"])
  mi <- mitools::MIcombine(results = lapply(seq_len(m), function(i) per$estimate[i, moved]),
                           variances = lapply(seq_len(m), function(i)
                             diag(per$variance[i, moved], sum(moved))), df.complete = Inf)
  # mitools reports no relative increase in variance; only what it reports is compared.
  mitools_values <- list(estimate = unname(mi$coefficients),
                         standard_error = unname(sqrt(diag(mi$variance))), df = unname(mi$df),
                         fmi = unname(mi$missinfo))
  rbind(.quantity_rows(case, "mitml", ours, lapply(mitml_values, unname)),
        .quantity_rows(case, "mitools", ours, mitools_values))
}

compare_all_pool <- function(cases = pool_cases()) {
  out <- do.call(rbind, lapply(cases, compare_pool_case))
  out$fn <- "pool_imputations"
  out
}
