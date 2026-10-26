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
  stack_cases(lapply(cases, function(case) do.call(.mclust_case, case)))
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
                              control = mclust::emControl(tol = c(1e-15, 1e-15),
                                                          itmax = c(1e5, 1e5)))
  # mclust stops on a relative change in the log-likelihood. Where the
  # likelihood is flat that leaves the estimates short of the fixed point: on
  # faithful EEE with three components one more EM step still moves a mean by
  # 2.7e-7, and a borderline posterior sits 1.1e-6 from its limit. The winning
  # solution is refined with mclust's own EM, without a stopping tolerance, as
  # the poLCA suite refines poLCA's, so the reference is mclust's estimate and
  # not where its stopping rule happened to fire.
  reference <- .mclust_refine(reference, model_data)
  as_published <- mclust::Mclust(model_data, G = g, modelNames = model,
                                 verbose = FALSE)
  covariance <- .mclust_covariance(model)
  fit <- multilpa(frame, vars = indicators, id = "unit", n_profiles = g,
                  n_group_classes = 1L, n_starts = 60L, seed = 1L,
                  tol = 1e-15, max_iter = 20000L,
                  covariance_model = covariance$covariance_model,
                  variance_model = covariance$variance_model)

  label <- sprintf("%s, %s, %d components", dataset, model, g)

  # mclust initializes from model-based hierarchical clustering and takes no
  # random restarts, so the two programs can land on different modes. Where they
  # do, parameters are not comparable; what is comparable, and what is checked
  # instead, is that multilpa's maximum is no lower.
  shortfall <- reference$loglik - as.numeric(logLik(fit))
  same_mode <- abs(shortfall) <= 1e-4

  if (!same_mode) {
    # Different modes: nothing but the shortfall is comparable, and the ledger
    # records that every output went uncompared for this reason.
    rows <- rbind(
      compare_values(sprintf("%s: free parameters", label),
                     reference$df, fit$n_parameters, tolerance = 0.5,
                     output = "model$n_parameters"),
      compare_values(sprintf("%s: log-likelihood shortfall against mclust", label),
                     0, max(shortfall, 0), tolerance = 1e-4))
    inventory <- fit_outputs(fit)
    case <- harness_case(rows, output_ledger(label, inventory, rows, exempt = stats::setNames(
      rep("mclust stopped at a different mode, so its estimates describe another solution",
          length(unique(inventory$table))), unique(inventory$table))))
  } else {
    case <- compare_fit(
      label, fit, .mclust_reference(reference, model, frame[indicators]),
      # The standard errors come from a numerical Hessian (Richardson
      # extrapolation) of the likelihood at mclust's estimates; the package's
      # come from its own analytic or numerical information at its own.
      tolerance = c(default = 1e-6, `profiles$mean_standard_error` = 1e-4,
                    `profiles$variance_standard_error` = 1e-4),
      extra = compare_values(sprintf("%s: log-likelihood shortfall against mclust", label),
                             0, max(shortfall, 0), tolerance = 1e-4))
  }
  source_label <- sprintf("mclust %s", utils::packageVersion("mclust"))
  case$comparisons <- transform(case$comparisons, source = source_label,
                                dataset = dataset, precision = "machine")

  if (is.null(published)) return(case)
  # The published value is checked against the mclust run that produced it, at
  # the precision it was printed to. This establishes where the number came
  # from; it is not a second equivalence claim, and stating it against the
  # converged fit would charge mclust's stopping rule to multilpa.
  case$comparisons <- rbind(case$comparisons, .with_output_columns(compare_printed(
    sprintf("%s: published %s reproduced by mclust at its default tolerance",
            label, c("BIC", "log-likelihood", "free parameters")),
    as.numeric(published[c("bic", "log_likelihood", "df")]),
    c(as_published$bic, as_published$loglik, as_published$df),
    digits = published_digits) |>
      transform(source = citation, dataset = dataset, precision = "printed")))
  case
}

#' Continue mclust's EM from its own solution to the fixed point
#'
#' @param reference A converged `Mclust` fit.
#' @param data The data it was fitted to.
#' @return `reference` with refined `parameters`, `z` and `loglik`.
.mclust_refine <- function(reference, data) {
  refined <- mclust::em(data = data, modelName = reference$modelName,
                        parameters = reference$parameters,
                        control = mclust::emControl(tol = 0, itmax = 20000L))
  stopifnot("mclust's refining EM must finish without an error code" =
              identical(as.numeric(attr(refined, "returnCode")), 0))
  log_likelihood <- sum(mclust::dens(data, modelName = reference$modelName,
                                     parameters = refined$parameters, logarithm = TRUE))
  stopifnot("refining must not lose likelihood" =
              log_likelihood >= reference$loglik - 1e-9)
  reference$parameters <- refined$parameters
  reference$z <- refined$z
  reference$loglik <- log_likelihood
  reference
}

#' mclust's solution in the harness's reference shape
#'
#' Covariances are expanded to one matrix per component whatever the model, so
#' the shared and diagonal structures are read the same way. The per-unit
#' log-likelihood contributions are mclust's own `dens()`; the standard errors
#' come from [.mclust_standard_errors()].
#'
#' @param reference The converged `Mclust` fit.
#' @param model Its model name.
#' @param data The indicator data frame.
#' @return A `harness_reference`.
.mclust_reference <- function(reference, model, data) {
  indicators <- names(data)
  g <- reference$G
  j <- length(indicators)
  means <- matrix(as.numeric(reference$parameters$mean), g, j, byrow = TRUE,
                  dimnames = list(NULL, indicators))
  variance <- reference$parameters$variance
  sigma <- if (j == 1L) {
    array(rep_len(variance$sigmasq, g), c(1L, 1L, g))
  } else array(variance$sigma, c(j, j, g))
  dimnames(sigma) <- list(indicators, indicators, NULL)
  harness_reference(
    log_likelihood = reference$loglik, n_parameters = reference$df,
    posteriors = reference$z, means = means, covariances = sigma,
    standard_errors = .mclust_standard_errors(means, sigma, reference$parameters$pro,
                                              model, as.matrix(data)),
    profile_probabilities = matrix(reference$parameters$pro, 1L),
    group_probabilities = 1,
    group_log_likelihood = stats::setNames(
      mclust::dens(as.matrix(data), modelName = reference$modelName,
                   parameters = reference$parameters, logarithm = TRUE),
      seq_len(nrow(data))),
    data = data)
}

#' Observed-information standard errors of a Gaussian mixture, from scratch
#'
#' The log-likelihood is written here in base R (Cholesky log densities and a
#' log-sum-exp), parameterised by the free natural parameters of the mclust
#' model: the means, the distinct covariance elements (shared across components
#' for an `E` structure, diagonal for an `I` one) and the first `G - 1` mixing
#' proportions. Its Hessian at mclust's estimates is taken by Richardson
#' extrapolation (numDeriv). At a maximum the inverse observed information
#' transforms with the Jacobian of any reparameterisation, so these are the
#' standard errors of the natural parameters whatever scale the package
#' optimises on.
#'
#' @return A list with `mean` and `variance`, components-by-indicators.
.mclust_standard_errors <- function(means, sigma, pro, model, data) {
  require_suite_packages("numDeriv", reason = "the reference observed information")
  g <- nrow(means)
  j <- ncol(means)
  shared <- substr(model, 1L, 1L) == "E"
  diagonal <- nchar(model) == 1L || substr(model, 3L, 3L) == "I"
  cells <- if (diagonal) cbind(seq_len(j), seq_len(j)) else
    which(lower.tri(diag(j), diag = TRUE), arr.ind = TRUE)
  blocks <- if (shared) 1L else g
  unpack <- function(theta) {
    mu <- matrix(theta[seq_len(g * j)], g, j)
    covariance_values <- matrix(theta[g * j + seq_len(nrow(cells) * blocks)],
                                nrow(cells), blocks)
    covariances <- lapply(seq_len(g), function(k) {
      m <- matrix(0, j, j)
      m[cells] <- covariance_values[, if (shared) 1L else k]
      m[cells[, 2:1, drop = FALSE]] <- covariance_values[, if (shared) 1L else k]
      m
    })
    first <- theta[g * j + nrow(cells) * blocks + seq_len(g - 1L)]
    list(mu = mu, covariances = covariances, pro = c(first, 1 - sum(first)))
  }
  log_likelihood <- function(theta) {
    parameters <- unpack(theta)
    if (any(parameters$pro <= 0)) return(-Inf)
    by_component <- vapply(seq_len(g), function(k) {
      root <- chol(parameters$covariances[[k]])
      standardised <- backsolve(root, t(data) - parameters$mu[k, ], transpose = TRUE)
      log(parameters$pro[[k]]) - 0.5 * colSums(standardised^2) -
        sum(log(diag(root))) - 0.5 * j * log(2 * pi)
    }, numeric(nrow(data)))
    by_component <- matrix(by_component, nrow(data))
    peak <- apply(by_component, 1L, max)
    sum(peak + log(rowSums(exp(by_component - peak))))
  }
  theta <- c(as.vector(means),
             as.vector(vapply(seq_len(blocks), function(k) matrix(sigma[, , k], j, j)[cells],
                              numeric(nrow(cells)))),
             pro[seq_len(g - 1L)])
  covariance <- solve(-numDeriv::hessian(log_likelihood, theta))
  errors <- sqrt(diag(covariance))
  variance_rows <- matrix(errors[g * j + seq_len(nrow(cells) * blocks)], nrow(cells), blocks)
  on_diagonal <- cells[, 1L] == cells[, 2L]
  variance_errors <- t(vapply(seq_len(g), function(k) {
    values <- variance_rows[on_diagonal, if (shared) 1L else k]
    values[order(cells[on_diagonal, 1L])]
  }, numeric(j)))
  list(mean = matrix(errors[seq_len(g * j)], g, j, dimnames = dimnames(means)),
       variance = matrix(variance_errors, g, j, dimnames = dimnames(means)))
}
