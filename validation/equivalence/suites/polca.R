# Single-level limit against poLCA (Linzer & Lewis 2011, J Stat Soft 42(10)).
#
# With one group class the two-level likelihood factorizes over individuals, so
# multilpa must reduce exactly to ordinary latent class analysis. Three claims
# are separated here, because a converged-to-converged comparison conflates
# them:
#   1. optimum      - both programs find the same maximised log-likelihood;
#   2. likelihood   - multilpa scores poLCA's own estimates identically, which
#                     tests the likelihood function rather than the optimizer;
#   3. formula      - a closed-form latent class log-likelihood written here,
#                     independent of both programs, agrees at those estimates.

#' Closed-form latent class log-likelihood, written from the definition
#'
#' @param codes Integer matrix of category codes, observations by items.
#' @param shares Numeric vector of class shares.
#' @param probabilities List of class-by-category matrices, one per item.
#' @return The observed-data log-likelihood, a single numeric value.
.textbook_lca_log_likelihood <- function(codes, shares, probabilities) {
  per_class <- vapply(seq_along(shares), function(class) {
    contributions <- vapply(seq_len(ncol(codes)), function(item) {
      code <- codes[, item]
      value <- log(probabilities[[item]][class, ][code])
      value[is.na(code)] <- 0          # observed-data likelihood: skip missing
      value
    }, numeric(nrow(codes)))
    log(shares[[class]]) + rowSums(contributions)
  }, numeric(nrow(codes)))
  peak <- apply(per_class, 1L, max)
  sum(peak + log(rowSums(exp(per_class - peak))))
}

suite_polca <- function() {
  require_suite_packages(
    "poLCA",
    reason = "the single-level categorical reference and its five bundled datasets")
  cases <- list(
    list(dataset = "carcinoma", items = c("A", "B", "C", "D", "E", "F", "G"),
         classes = 2:4),
    list(dataset = "cheating",
         items = c("LIEEXAM", "LIEPAPER", "FRAUD", "COPYEXAM"), classes = 2:3),
    list(dataset = "gss82",
         items = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"), classes = 2:3),
    list(dataset = "values", items = c("A", "B", "C", "D"), classes = 2L),
    list(dataset = "election",
         items = c("MORALG", "CARESG", "KNOWG", "LEADG", "DISHONG", "INTELG"),
         classes = 2:3)
  )
  do.call(rbind, lapply(cases, function(case) {
    do.call(rbind, lapply(case$classes, function(k) .polca_case(case$dataset, case$items, k)))
  }))
}

#' Run one dataset-by-class-count poLCA comparison
#' @param dataset Name of a poLCA dataset.
#' @param items Character vector of manifest item names.
#' @param k Number of latent classes.
#' @return A `data.frame` of compared quantities.
.polca_case <- function(dataset, items, k) {
  raw <- get(utils::data(list = dataset, package = "poLCA",
                         envir = environment()), envir = environment())
  # Compared on complete cases only. poLCA's `na.rm = FALSE` path does not
  # report the observed-data log-likelihood: on `election` it returns
  # -10848.839695 where both multilpa and a textbook likelihood written from
  # the definition give -10834.784505 at poLCA's own estimates, and its class
  # shares sum to 0.99216 rather than one. Equivalence under missing data is
  # established against glca instead, which does agree to 6e-08.
  observed <- raw[stats::complete.cases(raw[, items, drop = FALSE]), items, drop = FALSE]
  frame <- observed
  frame$unit <- factor(seq_len(nrow(frame)))

  formula <- stats::as.formula(paste("cbind(", paste(items, collapse = ","), ") ~ 1"))
  set.seed(20260919L)
  # poLCA's convergence test is on the absolute log-likelihood change; 1e-10 is
  # the tightest it reaches in practice, and asking for less simply exhausts
  # maxiter on every replication.
  reference <- poLCA::poLCA(formula, data = observed, nclass = k, nrep = 30L,
                            maxiter = 10000L, tol = 1e-10, verbose = FALSE,
                            na.rm = TRUE, calc.se = FALSE)

  fit <- multilpa(frame, vars = items, id = "unit", n_profiles = k,
                  n_group_classes = 1L, categorical = items, n_starts = 30L,
                  seed = 20260919L, tol = 1e-12, max_iter = 20000L)

  # Score poLCA's own estimates under the multilpa likelihood.
  at_reference <- multilpa(frame, vars = items, id = "unit", n_profiles = k,
                           n_group_classes = 1L, categorical = items, n_starts = 1L,
                           max_iter = 0L,
                           # poLCA labels its response columns "Pr(1)", "Pr(2)";
                           # multilpa checks a supplied start's labels against
                           # the data's own category levels and refuses a
                           # mismatch, so the foreign labels are dropped rather
                           # than renamed to something they do not mean.
                           start = starting_values(list(
                             profile_probabilities = matrix(reference$P, 1L, k),
                             group_probabilities = 1,
                             response_probabilities = lapply(reference$probs, unname))))

  codes <- vapply(items, function(item) as.integer(factor(observed[[item]])),
                  integer(nrow(observed)))
  textbook <- .textbook_lca_log_likelihood(codes, reference$P, reference$probs)

  class_order <- align_classes(flatten_blocks(fit$response_probabilities),
                               flatten_blocks(reference$probs))
  aligned <- lapply(reference$probs, function(block) block[class_order, , drop = FALSE])
  label <- sprintf("%s, %d classes", dataset, k)
  n_categories <- vapply(fit$response_probabilities, ncol, integer(1))
  cell_names <- unlist(lapply(items, function(item) {
    sprintf("%s: P(%s = %d | profile %d)", label, item,
            rep(seq_len(ncol(aligned[[item]])), each = k), rep(seq_len(k), ncol(aligned[[item]])))
  }))

  # Parameters are only worth comparing where they are identified. The
  # contingency table has prod(categories) - 1 degrees of freedom; when that
  # barely exceeds the number of free parameters the likelihood has a flat
  # ridge, and two programs reaching the same maximum can still report quite
  # different estimates. `cheating` with three classes is such a case: four
  # binary items give 15 degrees of freedom against 14 free parameters, and
  # poLCA and multilpa agree on the log-likelihood to 3e-08 while reporting
  # class shares of 0.035/0.111/0.854 and 0.045/0.078/0.877.
  residual_df <- prod(vapply(fit$response_probabilities, ncol, integer(1))) - 1L -
    fit$n_parameters
  identified <- residual_df >= 3L

  fit_rows <- rbind(
    compare_values(sprintf("%s: maximised log-likelihood", label),
                   reference$llik, as.numeric(logLik(fit)), tolerance = 1e-6),
    compare_values(sprintf("%s: free parameters", label),
                   reference$npar, fit$n_parameters, tolerance = 0.5),
    compare_values(sprintf("%s: BIC", label),
                   reference$bic, fit$bic_individual, tolerance = 1e-5),
    compare_values(sprintf("%s: AIC", label),
                   reference$aic, fit$aic, tolerance = 1e-5),
    compare_values(sprintf("%s: log-likelihood at poLCA's estimates", label),
                   reference$llik, as.numeric(logLik(at_reference)), tolerance = 1e-6),
    compare_values(sprintf("%s: textbook log-likelihood at poLCA's estimates", label),
                   reference$llik, textbook, tolerance = 1e-6),
    # poLCA defines this as min(N, prod(categories) - 1) - npar: the response
    # pattern table bounds the degrees of freedom once it is smaller than the
    # sample.
    compare_values(sprintf("%s: residual degrees of freedom", label),
                   reference$resid.df,
                   min(fit$n_informative, prod(n_categories) - 1L) - fit$n_parameters,
                   tolerance = 0.5))
  if (!identified) {
    return(transform(fit_rows, source = "poLCA 1.6.0.1", dataset = dataset))
  }
  rbind(fit_rows,
    # poLCA's convergence test is on the absolute change in the log-likelihood
    # and reaches about 1e-10, which leaves its parameters determined to roughly
    # 1e-5. Tightening it further only exhausts maxiter.
    compare_values(sprintf("%s: class share %d", label, seq_len(k)),
                   sort(reference$P), sort(as.numeric(fit$profile_probabilities)),
                   tolerance = 1e-4),
    compare_values(cell_names,
                   unlist(lapply(aligned, as.numeric)),
                   unlist(lapply(fit$response_probabilities, as.numeric)),
                   tolerance = 1e-4)
  ) |> transform(source = "poLCA 1.6.0.1", dataset = dataset)
}
