# Multiple-imputation pooling, pool_imputations(), against mice's Rubin's rules
# (mice::pool.scalar(), van Buuren & Groothuis-Oudshoorn 2011; Rubin 1987)
# and against a likelihood written from scratch for each imputation's fit.
#
# Data: srl (300 students, five continuous scales) with 30 values (10%) set
# missing completely at random in each of intrinsic_value, self_regulation
# and test_anxiety under base-R seed 1, imputed by mice(m = 5, method = "pmm",
# seed = 1). Two models, each pooled over the five completed data sets:
#
#   A. five indicators, two profiles, constant class probabilities;
#      level = 0.95, adjust = "none" (the documented defaults);
#   B. four indicators, two profiles, test_anxiety as a membership covariate
#      (the use case the help page names: a covariate imputed, then pooled);
#      level = 0.90, adjust = "BH", so the level and the correction are
#      exercised away from their defaults.
#
# Every table the result reports is compared column by column, through the
# verb ledger (compare_verb()), so a value nobody checked fails the run:
#
#   table        keys                                       reference
#   estimates    level, outcome, term, parameter            mice::pool.scalar() per row
#   covariance   row, column (the names vcov() gives)       Ubar + (1 + 1/m) B, base R
#   coef         parameter                                  pool.scalar's qbar
#   confint      parameter (at a level other than `level`)  qbar -/+ t_df sqrt(T)
#   summary      one row                                    see below
#   imputations  imputation, level, outcome, term, param.   from-scratch maximum
#   fits         imputation                                 from-scratch maximum
#
# Two questions, kept apart (they have different sources and tolerances):
#
# 1. Is the pooling arithmetic right? As the review's design asks, mice's
#    pool.scalar(Q, U, n = Inf, rule = "rubin1987") is applied to the
#    package's own per-imputation estimates Q and squared standard errors U
#    (the `imputations` table), so this isolates the arithmetic at machine
#    precision. Correspondence: qbar <-> estimate; ubar <-> within; b <->
#    between; t <-> standard_error^2; r <-> riv; df <-> df (n = Inf gives
#    Rubin's (m - 1)(1 + 1/r)^2, the package's "large complete-data sample");
#    fmi <-> fmi, (r + 2/(df + 3)) / (r + 1) in both. statistic = qbar /
#    sqrt(t) and p = 2 P(t_df > |statistic|) only where ?pool_imputations says
#    a test is reported (not for variances and probabilities); p_adjusted is
#    stats::p.adjust() over those p-values; intervals are qbar -/+
#    qt((1 + level)/2, df) sqrt(t). Where pool.scalar's arithmetic is 0/0 --
#    the single group class's probability, the constant one, with ubar = b = 0
#    -- the reference takes the limit the help page documents (no variance, no
#    missing information: riv 0, df Inf, fmi 0); ubar = 0 < b would take riv
#    Inf, df m - 1, fmi 1, and is asserted absent here. The pooled covariance
#    is mean_i(V_i) + (1 + 1/m) cov(Q), with V_i each imputation's vcov() from
#    a direct multilpa() call, relabelled by the map below.
#
# 2. Is each imputation's fit, as pooled, the maximum of the right likelihood
#    on the right completed data set, relabelled consistently? The reference
#    (pool-imputations-reference.R) writes the Gaussian mixture likelihood and
#    its analytic gradient in base R, carries a start to the maximum by Newton
#    steps (numDeriv Hessian of the analytic gradient), and takes standard
#    errors from the inverse observed information on the natural scale. The
#    start is a direct multilpa() call's estimate; the fixed point is the
#    reference's own. Profiles are aligned by the permutation minimising the
#    squared distance between profile means and imputation 1's, and logits are
#    re-expressed against the new last profile. `fits$profile_order` is
#    compared against that permutation as text (with two profiles a
#    permutation is its own inverse, so this cannot tell the two conventions
#    apart). The direct multilpa() calls also give a same-computation check
#    (1e-10): the pooled fits are the fits a user gets by completing each data
#    set with mice::complete() and calling multilpa() with the same arguments.
#
# Summary row: m and the design (n_profiles 2, n_group_classes 1) from the
# call; log_likelihood_min/max from the from-scratch maxima; relabelled from
# the reference permutations; fmi_max/median from pool.scalar; model and
# all_converged from the direct calls (bookkeeping, same computation).
#
# The package's fits are run with tol = 1e-15, as the glca suite tightens its,
# so a comparison against a fixed point measures the implementation and not
# where the EM's stopping rule fired: at 1e-13 the covariate fits stopped
# 1.2e-6 to 1.9e-6 (up to 8e-6 standard errors) short of the maximum in their
# intercept, with the likelihood 5e-11 to 9e-11 below it; at 1e-15 the gap is
# 2.6e-7 at most, in at most 36 EM sweeps.

source(file.path("equivalence", "harness", "pool-imputations-reference.R"), local = environment())

suite_pool_imputations <- function() {
  require_suite_packages(
    "mice", "numDeriv",
    reason = "mice's Rubin's rules (pool.scalar) and the imputations, and the reference observed information")
  imputed <- .pool_imputed_srl()
  scales <- names(latents::srl)
  cases <- list(
    list(label = "srl, 2 profiles", vars = scales, covariates = character(),
         level = 0.95, adjust = "none", interval_level = 0.8),
    list(label = "srl, 2 profiles, test_anxiety on membership",
         vars = setdiff(scales, "test_anxiety"), covariates = "test_anxiety",
         level = 0.9, adjust = "BH", interval_level = 0.99))
  stack_cases(lapply(cases, function(case) do.call(.pool_case, c(list(imputed = imputed), case))))
}

#' Evaluate code under a seed, restoring the caller's random state
#' @param seed Integer seed.
#' @param code Code to evaluate.
#' @return The value of `code`.
.pool_with_seed <- function(seed, code) {
  saved <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(saved)) {
    if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", saved, globalenv()), add = TRUE)
  set.seed(seed)
  code
}

#' srl with 10% MCAR in three scales, imputed five times by mice
#' @return A `mids` object.
.pool_imputed_srl <- function() {
  observed <- latents::srl
  holed <- .pool_with_seed(1L, Reduce(function(frame, scale) {
    frame[[scale]][sample(nrow(frame), 30L)] <- NA
    frame
  }, c("intrinsic_value", "self_regulation", "test_anxiety"), observed))
  stopifnot("30 values must be missing in each of three scales" =
              identical(unname(colSums(is.na(holed))), c(0, 30, 0, 30, 30)))
  # mice() calls set.seed(seed) itself; the caller's stream is restored.
  imputed <- .pool_with_seed(1L, mice::mice(holed, m = 5L, method = "pmm", seed = 1L,
                                            printFlag = FALSE))
  stopifnot("every completed data set must be complete" =
              all(vapply(seq_len(imputed$m), function(i) {
                !anyNA(mice::complete(imputed, i))
              }, logical(1))))
  imputed
}

#' The text key of a parameter row
#' @param table A `data.frame` with `level`, `outcome`, `term`, `parameter`.
#' @return A character vector.
.pool_key <- function(table) {
  paste(table$level, table$outcome, table$term, table$parameter, sep = "\r")
}

#' Give the missing term of the group row a text key
#'
#' The single group class's probability has no term (`NA`). compare_tables()
#' joins NA keys but counts a row with a missing key as not found, so the key
#' is written as text on both sides; a row the fit reports with a term would
#' then fail to join.
#' @param table A parameter table.
#' @return The table with `term` NA replaced by `"(none)"`.
.pool_text_term <- function(table) {
  table$term <- ifelse(is.na(table$term), "(none)", table$term)
  table
}

#' One pooled model against both references
#' @param imputed The `mids` object.
#' @param label Label of the case.
#' @param vars,covariates Indicators and membership covariates.
#' @param level,adjust Passed to pool_imputations().
#' @param interval_level The level confint() is asked for.
#' @return A `harness_case`.
.pool_case <- function(imputed, label, vars, covariates, level, adjust, interval_level) {
  m <- imputed$m
  # The tightest stopping rule (see the header): the covariate EM settled in
  # at most 36 sweeps, and the comparison is then with the maximum, not the stop.
  pooled <- pool_imputations(imputed, vars, id = NULL, n_profiles = 2L,
                             profile_covariates = covariates, n_starts = 3L, seed = 1L,
                             tol = 1e-15, max_iter = 20000L, level = level, adjust = adjust)
  direct <- lapply(seq_len(m), function(i) {
    multilpa(mice::complete(imputed, i), vars, id = NULL, n_profiles = 2L,
             profile_covariates = covariates, n_starts = 3L, seed = 1L,
             tol = 1e-15, max_iter = 20000L)
  })
  layout <- pool_reference_layout(vars, 2L, covariates)
  reported <- layout$reported
  keys <- reported[c("level", "outcome", "term", "parameter")]
  names_all <- pool_reference_names(keys)

  per_imputation <- lapply(seq_len(m), function(i) {
    .pool_imputation_reference(direct[[i]], mice::complete(imputed, i), layout)
  })
  target <- per_imputation[[1L]]$means
  perms <- lapply(per_imputation, function(one) pool_reference_permutation(one$means, target))
  maps <- lapply(perms, pool_reference_relabel, layout = layout)

  # Question 2: each imputation's fit, against the from-scratch maximum.
  expected_imputations <- do.call(rbind, lapply(seq_len(m), function(i) {
    one <- per_imputation[[i]]
    covariance <- maps[[i]] %*% one$reference$covariance %*% t(maps[[i]])
    data.frame(imputation = i, keys,
               estimate = as.vector(maps[[i]] %*% one$reference$estimate),
               standard_error = sqrt(pmax(diag(covariance), 0)), stringsAsFactors = FALSE)
  }))
  expected_fits <- data.frame(
    imputation = seq_len(m),
    log_likelihood = vapply(per_imputation, function(one) one$log_likelihood, numeric(1)),
    n_parameters = sum(reported$free),
    profile_order = vapply(perms, paste, character(1), collapse = ","),
    # One group class has nothing to match.
    group_order = "1", stringsAsFactors = FALSE)
  imputation_keys <- c("imputation", "level", "outcome", "term", "parameter")
  obtained_fit_tables <- list(
    imputations = keyed_table(.pool_text_term(get_results(pooled, "imputations")),
                              imputation_keys),
    fits = keyed_table(get_results(pooled, "fits"), "imputation"))

  # The same computation, called directly: aligned by the same maps.
  direct_imputations <- do.call(rbind, lapply(seq_len(m), function(i) {
    one <- per_imputation[[i]]
    covariance <- maps[[i]] %*% one$direct_covariance %*% t(maps[[i]])
    data.frame(imputation = i, keys,
               estimate = as.vector(maps[[i]] %*% one$direct_estimate),
               standard_error = sqrt(pmax(diag(covariance), 0)), stringsAsFactors = FALSE)
  }))
  direct_fits <- data.frame(
    imputation = seq_len(m),
    log_likelihood = vapply(direct, function(fit) fit$log_likelihood, numeric(1)),
    n_parameters = vapply(direct, function(fit) as.numeric(fit$n_parameters), numeric(1)),
    converged = vapply(direct, function(fit) isTRUE(fit$converged), logical(1)))
  consistency <- compare_tables(
    sprintf("%s, against direct multilpa() calls", label), obtained_fit_tables,
    list(imputations = .pool_text_term(direct_imputations), fits = direct_fits),
    # The same arithmetic on the same data, apart from the order of a
    # relabelling: machine precision with a margin for the numerical
    # derivatives behind the standard errors.
    tolerance = c(default = 1e-10))

  fit_case <- compare_verb(
    sprintf("%s, per-imputation fits", label), obtained_fit_tables,
    list(imputations = .pool_text_term(expected_imputations), fits = expected_fits),
    tolerance = c(
      # Two maxima of the same likelihood, the package's stopped by its EM rule
      # at tol = 1e-15 and the reference's by Newton steps below 1e-11. The
      # package default: two converged optimisers agree to about the square
      # root of their stopping tolerance.
      default = 1e-6,
      # Inverse observed information at maxima 1e-6 apart: one computed by the
      # package, one by numDeriv's Richardson Jacobian of the analytic
      # gradient. 1e-6 absolute is 1e-4 of the smallest standard error
      # (0.0077), a gap the curvature cannot move by a parameter gap of 1e-6.
      `imputations$standard_error` = 1e-6,
      # A likelihood at two points 1e-6 apart near its maximum differs by the
      # square of that gap times the curvature, far below 1e-10 relative.
      `fits$log_likelihood` = 1e-10),
    extra = consistency)
  fit_case$comparisons$source <- ifelse(
    grepl("against direct multilpa", fit_case$comparisons$quantity, fixed = TRUE),
    sprintf("latents %s multilpa() called directly (same computation)",
            utils::packageVersion("latents")),
    sprintf("base-R Gaussian mixture likelihood, numDeriv %s", utils::packageVersion("numDeriv")))

  # Question 1: the pooling arithmetic, on the package's own Q and U.
  imputations <- get_results(pooled, "imputations")
  row_of <- match(.pool_key(imputations), .pool_key(keys))
  stopifnot("every imputation row must be a reference parameter" = !anyNA(row_of),
            "every imputation must report every parameter once" =
              nrow(imputations) == m * nrow(keys))
  cell <- cbind(row_of, imputations$imputation)
  q <- matrix(NA_real_, nrow(keys), m)
  u <- q
  q[cell] <- imputations$estimate
  u[cell] <- imputations$standard_error^2
  rubin <- .pool_rubin(q, u)
  tested <- !(keys$parameter %in% c("variance", "probability"))
  statistic <- ifelse(tested, rubin$qbar / sqrt(rubin$t), NA_real_)
  p_value <- ifelse(tested, 2 * stats::pt(-abs(statistic), rubin$df), NA_real_)
  p_adjusted <- p_value
  p_adjusted[tested] <- stats::p.adjust(p_value[tested], method = adjust)
  half_width <- function(at) stats::qt((1 + at) / 2, rubin$df) * sqrt(rubin$t)
  expected_estimates <- data.frame(
    keys, estimate = rubin$qbar, standard_error = sqrt(rubin$t), statistic = statistic,
    p_value = p_value, p_adjusted = p_adjusted, conf_low = rubin$qbar - half_width(level),
    conf_high = rubin$qbar + half_width(level), df = rubin$df, within = rubin$ubar,
    between = rubin$b, riv = rubin$r, fmi = rubin$fmi, stringsAsFactors = FALSE)
  within <- Reduce(`+`, lapply(seq_len(m), function(i) {
    maps[[i]] %*% per_imputation[[i]]$direct_covariance %*% t(maps[[i]])
  })) / m
  total <- within + (1 + 1 / m) * stats::cov(t(q))
  dimnames(total) <- list(names_all, names_all)
  intervals <- confint(pooled, level = interval_level)
  obtained_pool_tables <- list(
    estimates = keyed_table(.pool_text_term(get_results(pooled, "estimates")),
                            c("level", "outcome", "term", "parameter")),
    covariance = keyed_table(.pool_long_matrix(vcov(pooled)), c("row", "column")),
    coef = keyed_table(data.frame(parameter = names(coef(pooled)),
                                  estimate = unname(coef(pooled))), "parameter"),
    confint = keyed_table(data.frame(parameter = rownames(intervals),
                                     lower = intervals[, 1L], upper = intervals[, 2L]),
                          "parameter"),
    summary = keyed_table(data.frame(.row = 1L, summary(pooled)), ".row"))
  expected_pool_tables <- list(
    estimates = .pool_text_term(expected_estimates),
    covariance = .pool_long_matrix(total),
    coef = data.frame(parameter = names_all, estimate = rubin$qbar),
    confint = data.frame(parameter = names_all,
                         lower = rubin$qbar - half_width(interval_level),
                         upper = rubin$qbar + half_width(interval_level)),
    summary = data.frame(
      .row = 1L, m = m, model = class(direct[[1L]])[1L], n_profiles = 2L,
      n_group_classes = 1L,
      all_converged = all(vapply(direct, function(fit) isTRUE(fit$converged), logical(1))),
      log_likelihood_min = min(expected_fits$log_likelihood),
      log_likelihood_max = max(expected_fits$log_likelihood),
      relabelled = sum(vapply(perms, function(perm) any(perm != seq_along(perm)), logical(1))),
      fmi_max = max(rubin$fmi), fmi_median = stats::median(rubin$fmi),
      stringsAsFactors = FALSE))
  pool_case <- compare_verb(
    sprintf("%s, pooled", label), obtained_pool_tables, expected_pool_tables,
    tolerance = c(
      # Rubin's rules on identical inputs: the same arithmetic in a different
      # order (mice reaches df through lambda = (1 + 1/m) b / t, the package
      # through r), so agreement is to rounding.
      default = 1e-10,
      # df = (m - 1) / lambda^2 doubles lambda's relative rounding, and lambda
      # rests on b, a variance of five nearly equal numbers; the gap is
      # judged relative to df (which runs to 1e5), as the review's design sets.
      `estimates$df` = 1e-8,
      # The two log-likelihood columns are the fits' own (question 2's
      # tolerance and reason).
      `summary$log_likelihood_min` = 1e-10, `summary$log_likelihood_max` = 1e-10))
  pool_case$comparisons$source <- sprintf(
    "mice %s pool.scalar (Rubin 1987)", utils::packageVersion("mice"))

  case <- stack_cases(list(pool_case, fit_case))
  case$comparisons <- transform(case$comparisons, dataset = "srl (10% MCAR, mice pmm m = 5)",
                                precision = "machine")
  case
}

#' One imputation: the from-scratch maximum, and the direct call's tables
#' @param fit The direct multilpa() fit on this completed data set.
#' @param completed The completed data set.
#' @param layout From [pool_reference_layout()].
#' @return A list: `reference` (estimate and covariance over the reported
#'   parameters), `log_likelihood`, `means` (K x p, for alignment),
#'   `direct_estimate` and `direct_covariance` (the direct call's, in the
#'   layout's order).
.pool_imputation_reference <- function(fit, completed, layout) {
  reported <- layout$reported
  table <- parameter_inference(fit)
  index <- match(.pool_key(reported), .pool_key(table))
  stopifnot("the fit must report exactly the layout's parameters" =
              !anyNA(index) && nrow(table) == nrow(reported))
  names_all <- pool_reference_names(reported)
  covariance <- vcov(fit)
  stopifnot("vcov() must name every reported parameter" =
              all(names_all %in% rownames(covariance)) && nrow(covariance) == nrow(reported))
  y <- as.matrix(completed[layout$vars])
  design <- if (length(layout$covariates) == 0L) NULL else
    cbind(1, as.matrix(completed[layout$covariates]))
  maximum <- pool_reference_maximise(table$estimate[index][reported$free], layout, y, design)
  at_maximum <- pool_reference_reported(maximum$theta, maximum$covariance, layout)
  list(reference = at_maximum, log_likelihood = maximum$log_likelihood,
       means = matrix(at_maximum$estimate[reported$block == "mean"], layout$n_profiles,
                      length(layout$vars), byrow = TRUE),
       direct_estimate = table$estimate[index],
       direct_covariance = covariance[names_all, names_all])
}

#' mice's Rubin's rules, one parameter at a time
#'
#' Rows whose within- and between-imputation variances are both zero take the
#' documented limit (riv 0, df Inf, fmi 0); a row with zero within-variance
#' and positive between-variance is refused, since nothing here exercises it.
#' @param q,u Parameters-by-imputations estimates and squared standard errors.
#' @return A list of numeric vectors: qbar, ubar, b, t, r, df, fmi.
.pool_rubin <- function(q, u) {
  rows <- lapply(seq_len(nrow(q)), function(r) {
    mice::pool.scalar(q[r, ], u[r, ], n = Inf, k = 1, rule = "rubin1987")
  })
  take <- function(name) vapply(rows, function(row) as.numeric(row[[name]]), numeric(1))
  out <- lapply(stats::setNames(nm = c("qbar", "ubar", "b", "t", "r", "df", "fmi")), take)
  constant <- out$ubar == 0 & out$b == 0
  stopifnot("zero within-variance with positive between-variance is not exercised" =
              !any(out$ubar == 0 & out$b > 0),
            "pool.scalar must be finite wherever the variances are not both zero" =
              all(is.finite(c(out$r, out$df, out$fmi))[rep(!constant, 3L)]))
  out$r[constant] <- 0
  out$df[constant] <- Inf
  out$fmi[constant] <- 0
  out
}

#' A named square matrix as one row per cell
#' @param x A matrix with identical row and column names.
#' @return A `data.frame` with `row`, `column` and `covariance`.
.pool_long_matrix <- function(x) {
  data.frame(row = rep(rownames(x), times = ncol(x)),
             column = rep(colnames(x), each = nrow(x)),
             covariance = as.vector(x), stringsAsFactors = FALSE)
}
