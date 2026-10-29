# Functions for compare.R: fit Gaussian latent profile models with independent
# R implementations (Rmixmod for VVI, EEI, VVV and EEE; EMCluster and mixtools
# for VVV; mixtools for EEE) and with lpa(), put every solution in one shape,
# and compare the maximized likelihood, shares, means, (co)variances and every
# row's posteriors.
#
# Tolerances, fixed before any comparison was run: the log likelihood to 1e-6
# plus 1e-10 of its size; shares, means and (co)variances to 1e-4 of
# max(1, |value|); posteriors to 1e-4. Every package runs on one core.
# flexmix's FLXMCmvnorm is not used: it estimates covariances with the
# unbiased divisor, so its fixed point is not the maximum likelihood estimate.

.lpa_case <- function(name, data, classes, model) {
  references <- switch(model, VVI = "Rmixmod", EEI = "Rmixmod",
                       VVV = c("Rmixmod", "EMCluster", "mixtools"), EEE = c("Rmixmod", "mixtools"))
  list(name = sprintf("%s_%s_k%d", name, model, classes), data = data, vars = names(data),
       classes = classes, model = model, references = references)
}

#' The datasets, structures and class counts compared
lpa_cases <- function() {
  data_of <- function(name, package) {
    get(utils::data(list = name, package = package, envir = environment()), envir = environment())
  }
  faithful <- datasets::faithful
  iris <- datasets::iris[c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")]
  diabetes <- data_of("diabetes", "mclust")[c("glucose", "insulin", "sspg")]
  banknote <- data_of("banknote", "mclust")[c("Length", "Left", "Right", "Bottom", "Top", "Diagonal")]
  crabs <- data_of("crabs", "MASS")[c("FL", "RW", "CL", "CW", "BD")]
  thyroid <- data_of("thyroid", "mclust")[c("RT3U", "T4", "T3", "TSH", "DTSH")]
  list(
    .lpa_case("faithful", faithful, 2L, "VVV"), .lpa_case("faithful", faithful, 2L, "EEE"),
    .lpa_case("faithful", faithful, 3L, "VVI"), .lpa_case("faithful", faithful, 3L, "EEI"),
    .lpa_case("iris", iris, 3L, "VVV"), .lpa_case("iris", iris, 3L, "EEE"),
    .lpa_case("iris", iris, 2L, "VVI"), .lpa_case("iris", iris, 3L, "EEI"),
    .lpa_case("diabetes", diabetes, 3L, "VVV"), .lpa_case("diabetes", diabetes, 3L, "VVI"),
    .lpa_case("banknote", banknote, 2L, "VVV"), .lpa_case("banknote", banknote, 2L, "EEE"),
    .lpa_case("crabs", crabs, 2L, "VVV"), .lpa_case("crabs", crabs, 2L, "EEE"),
    .lpa_case("thyroid", thyroid, 3L, "VVI"), .lpa_case("thyroid", thyroid, 3L, "VVV"))
}

#' A solution in one shape
#' @return `log_likelihood`, `shares` (K), `means` (K x p), `covariances`
#'   (list of K p x p matrices; one repeated for an equal structure) and
#'   `posteriors` (n x K or `NULL`).
.solution <- function(log_likelihood, shares, means, covariances, posteriors = NULL) {
  list(log_likelihood = log_likelihood, shares = as.numeric(shares),
       means = unname(as.matrix(means)), covariances = lapply(covariances, function(m) unname(m)),
       posteriors = if (is.null(posteriors)) NULL else unname(as.matrix(posteriors)))
}

#' lpa() started at a solution (one start), or by its own search
#' @param start `NULL` for the search, or a `.solution()` to start from.
#' @param search Passed to lpa() when searching.
fit_latents_lpa <- function(case, start = NULL, search = NULL) {
  arguments <- list(case$data, case$vars, n_profiles = case$classes, model = case$model,
                    seed = 1L, tol = 1e-12, max_iter = 100000L)
  arguments <- if (is.null(start)) c(arguments, list(n_starts = 50L), if (!is.null(search))
    list(search = search)) else c(arguments, list(n_starts = 1L, start = .as_start(start, case)))
  fit <- do.call(lpa, arguments)
  theta <- stats::coef(fit)
  k <- case$classes
  p <- length(case$vars)
  value <- function(name) theta[[name]]
  means <- t(vapply(seq_len(k), function(j) vapply(case$vars, function(v)
    value(sprintf("measurement.mean.profile_%d.%s", j, v)), numeric(1)), numeric(p)))
  covariance_of <- function(j) {
    owner <- if (case$model %in% c("EEI", "EEE")) "shared" else sprintf("profile_%d", j)
    out <- matrix(0, p, p)
    invisible(lapply(seq_len(p), function(a) lapply(seq_len(a), function(b) {
      if (case$model %in% c("VVI", "EEI")) {
        if (a == b) out[a, a] <<- value(sprintf("measurement.variance.%s.%s", owner, case$vars[[a]]))
      } else {
        entry <- value(sprintf("measurement.covariance.%s.%s:%s", owner, case$vars[[a]], case$vars[[b]]))
        out[a, b] <<- entry
        out[b, a] <<- entry
      }
    })))
    out
  }
  .solution(as.numeric(stats::logLik(fit)),
            vapply(seq_len(k), function(j)
              value(sprintf("profile.probability.profile_%d.group_class_1", j)), numeric(1)),
            means, lapply(seq_len(k), covariance_of), fit$subject_posteriors)
}

#' A reference solution as an lpa() start
#' @noRd
.as_start <- function(solution, case) {
  out <- list(means = solution$means, profile_probabilities = matrix(solution$shares, 1L),
              group_probabilities = 1)
  if (case$model %in% c("VVI", "EEI")) {
    out$variances <- t(vapply(solution$covariances, diag, numeric(length(case$vars))))
  } else {
    out$covariances <- simplify2array(solution$covariances)
  }
  out
}

.rmixmod_models <- c(VVI = "Gaussian_pk_Lk_Bk", EEI = "Gaussian_pk_L_B",
                     VVV = "Gaussian_pk_Lk_Ck", EEE = "Gaussian_pk_L_C")

fit_rmixmod <- function(case) {
  fit <- Rmixmod::mixmodCluster(
    case$data, nbCluster = case$classes,
    models = Rmixmod::mixmodGaussianModel(listModels = .rmixmod_models[[case$model]]),
    strategy = Rmixmod::mixmodStrategy(nbTry = 20L, epsilonInAlgo = 1e-12,
                                       nbIterationInAlgo = 100000L), seed = 1L)
  best <- fit@bestResult
  if (!is.finite(best@likelihood) || anyNA(best@parameters@mean) || best@likelihood == 0) return(NULL)
  .solution(best@likelihood, best@parameters@proportions, best@parameters@mean,
            best@parameters@variance, best@proba)
}

#' A symmetric matrix from EMCluster's row-by-row lower triangle
#' @noRd
.lt_to_matrix <- function(lower, p) {
  out <- matrix(0, p, p)
  index <- which(lower.tri(out, diag = TRUE), arr.ind = TRUE)
  out[index[order(index[, "row"], index[, "col"]), , drop = FALSE]] <- lower
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  out
}

fit_emcluster <- function(case) {
  x <- as.matrix(case$data)
  fits <- lapply(seq_len(20L), function(s) withr::with_seed(s, {
    printed <- utils::capture.output(fitted <- EMCluster::emcluster(
      x, EMCluster::init.EM(x, nclass = case$classes, EMC = EMCluster::.EMControl()),
      EMC = EMCluster::.EMControl(em.iter = 100000L, em.eps = 1e-12)))
    fitted
  }))
  p <- ncol(x)
  # EMCluster keeps a start whose covariance collapses onto a few points (an
  # infinite or implausibly high likelihood); only proper solutions compete.
  scale <- max(apply(x, 2L, stats::var))
  proper <- vapply(fits, function(f) {
    is.finite(f$llhdval) && all(vapply(seq_len(case$classes), function(j) {
      min(eigen(.lt_to_matrix(f$LTSigma[j, ], p), symmetric = TRUE, only.values = TRUE)$values) >
        1e-8 * scale
    }, logical(1)))
  }, logical(1))
  if (!any(proper)) return(NULL)
  fits <- fits[proper]
  best <- fits[[which.max(vapply(fits, function(f) f$llhdval, numeric(1)))]]
  covariances <- lapply(seq_len(case$classes), function(j) .lt_to_matrix(best$LTSigma[j, ], p))
  .solution(best$llhdval, best$pi, best$Mu, covariances, EMCluster::e.step(x, emobj = best)$Gamma)
}

fit_mixtools <- function(case) {
  x <- as.matrix(case$data)
  fits <- lapply(seq_len(20L), function(s) withr::with_seed(s, tryCatch({
    printed <- utils::capture.output(fitted <- mixtools::mvnormalmixEM(
      x, k = case$classes, arbvar = identical(case$model, "VVV"), epsilon = 1e-12,
      maxit = 100000L))
    fitted
  }, error = function(e) {
    message(sprintf("mixtools start %d failed: %s", s, conditionMessage(e)))
    NULL
  })))
  fits <- Filter(Negate(is.null), fits)
  best <- fits[[which.max(vapply(fits, function(f) f$loglik, numeric(1)))]]
  covariances <- if (is.list(best$sigma)) best$sigma else rep(list(best$sigma), case$classes)
  .solution(best$loglik, best$lambda, do.call(rbind, best$mu), covariances, best$posterior)
}

.fitters <- list(Rmixmod = fit_rmixmod, EMCluster = fit_emcluster, mixtools = fit_mixtools)

.permutations <- function(k) {
  if (k == 1L) return(list(1L))
  unlist(lapply(seq_len(k), function(first) lapply(.permutations(k - 1L), function(rest)
    c(first, setdiff(seq_len(k), first)[rest]))), recursive = FALSE)
}

.relative_gap <- function(reference, ours) max(abs(reference - ours) / pmax(1, abs(reference)))

#' Compare one package's solution with latents' under the closest class order
#' @noRd
compare_solution <- function(case, package, ours, reference) {
  perms <- .permutations(case$classes)
  gap <- vapply(perms, function(perm) .relative_gap(reference$means[perm, , drop = FALSE], ours$means),
                numeric(1))
  perm <- perms[[which.min(gap)]]
  covariance_gap <- max(vapply(seq_len(case$classes), function(j)
    .relative_gap(reference$covariances[[perm[j]]], ours$covariances[[j]]), numeric(1)))
  k <- case$classes
  p <- length(case$vars)
  n_covariance <- if (case$model %in% c("VVI", "EEI")) p else p * (p + 1L) / 2
  n_covariance <- n_covariance * if (case$model %in% c("EEI", "EEE")) 1L else k
  rows <- data.frame(
    case = case$name, reference = package,
    quantity = c("maximized log likelihood", sprintf("shares (%d values)", k),
                 sprintf("means (%d values)", k * p),
                 sprintf("%s (%d values)", if (case$model %in% c("VVI", "EEI")) "variances" else
                   "covariances", n_covariance)),
    reference_value = c(reference$log_likelihood, NA, NA, NA),
    latents = c(ours$log_likelihood, NA, NA, NA),
    difference = c(abs(reference$log_likelihood - ours$log_likelihood),
                   .relative_gap(reference$shares[perm], ours$shares), min(gap), covariance_gap),
    tolerance = c(1e-6 + 1e-10 * abs(ours$log_likelihood), 1e-4, 1e-4, 1e-4),
    n_values = c(1L, k, k * p, n_covariance))
  if (!is.null(reference$posteriors)) {
    rows <- rbind(rows, data.frame(
      case = case$name, reference = package,
      quantity = sprintf("posteriors (%d rows x %d classes)", nrow(ours$posteriors), k),
      reference_value = NA, latents = NA,
      difference = max(abs(reference$posteriors[, perm, drop = FALSE] - ours$posteriors)),
      tolerance = 1e-4, n_values = length(ours$posteriors)))
  }
  rows$agrees <- rows$difference <= rows$tolerance
  rows
}

#' Fit every case with lpa() and its references and compare
#' @return One row per compared quantity, with `fn = "lpa"`.
#' Fit every case with its references and compare lpa() at each reference's
#' solution: started there, lpa() must reach the same stationary point.
#' lpa()'s own search (default and emEM) is recorded beside the best reference
#' maximum in the "search" attribute; it is not part of the comparison.
#' A reference that returns no solution gives a failed row, for the register.
#' @return One row per compared quantity, with `fn = "lpa"`.
compare_all_lpa <- function(cases = lpa_cases()) {
  search <- list()
  out <- do.call(rbind, lapply(cases, function(case) {
    references <- lapply(stats::setNames(case$references, case$references),
                         function(package) .fitters[[package]](case))
    rows <- do.call(rbind, lapply(names(references), function(package) {
      reference <- references[[package]]
      if (is.null(reference)) {
        return(data.frame(case = case$name, reference = package,
                          quantity = "a solution from the reference", reference_value = NA,
                          latents = NA, difference = NA, tolerance = NA, n_values = 0L,
                          agrees = FALSE))
      }
      compare_solution(case, package, fit_latents_lpa(case, start = reference), reference)
    }))
    found <- Filter(Negate(is.null), references)
    search[[case$name]] <<- data.frame(
      case = case$name,
      latents_default = fit_latents_lpa(case)$log_likelihood,
      latents_emEM = fit_latents_lpa(case, search = "emEM")$log_likelihood,
      best_reference = max(vapply(found, function(r) r$log_likelihood, numeric(1))))
    rows
  }))
  out$fn <- "lpa"
  search <- do.call(rbind, search)
  search$default_reaches <- search$latents_default >= search$best_reference - 1e-6
  search$emEM_reaches <- search$latents_emEM >= search$best_reference - 1e-6
  row.names(search) <- NULL
  attr(out, "search") <- search
  out
}
