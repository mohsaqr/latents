# Functions for compare.R: class enumeration and the bootstrap likelihood-ratio
# test against independent implementations.
#
#   enumerate_lpa() against mclust::mclustBIC() (four structures, one to four
#   profiles, three datasets; mclust's EM run to 1e-12, not its default 1e-5,
#   which stops up to 0.04 short) and enumerate_lca() against poLCA (one to four
#   classes, three datasets): wherever both reach the same maximum, the log
#   likelihood, the parameter count and the BIC of every grid cell; the
#   parameter count everywhere. Cells where the two searches reach different
#   maxima are recorded in the "search" attribute (not compared).
#
#   bootstrap_lrt() against mclust::mclustBootstrapLRT() on mclust's diabetes
#   (VVV): the observed statistic, and the bootstrap p-value within Monte
#   Carlo error; mclust's vignette prints the statistics. mclust's
#   hierarchical initialisation stops at lower maxima than lpa()'s 50 starts
#   at two and four profiles (a proper solution: mclust's own EM from
#   latents' posteriors stays there), so the statistics are compared at
#   mclust's solutions (lpa() started there); latents' own maxima are recorded
#   in the "search" attribute.
#
# Tolerances, fixed before the comparison: the log likelihood to 1e-6 plus
# 1e-10 of its size and the BIC (-2 LL plus a constant) to twice that;
# parameter counts exactly; a bootstrap p-value
# within four Monte Carlo standard errors of the two estimates combined
# (sqrt(p (1 - p) (1 / B1 + 1 / B2)), with p their pooled value) plus 1 / B;
# printed statistics to half a unit of their last digit.

# mclust's EM stops at a relative change of 1e-5 by default, up to 0.04 short
# of the maximum here (faithful EEI k4); every mclust fit is run to 1e-12.
.mclust_control <- function() mclust::emControl(tol = c(1e-12, 1e-12))

#' mclust's diabetes as its vignette uses it
#'
#' The vignette (mclust 6.1.2, "data corrections (see help("diabetes"))")
#' renames the columns glucose, insulin, sspg to fpg, glucose, insulin and
#' sets `glucose[104] <- 455`, which is the installed data's insulin[104]
#' (45 there).
#' @noRd
.diabetes_vignette <- function() {
  diabetes <- get(utils::data(list = "diabetes", package = "mclust", envir = environment()),
                  envir = environment())[c("glucose", "insulin", "sspg")]
  diabetes$insulin[[104L]] <- 455
  diabetes
}

.cell_row <- function(case, reference, quantity, value, obtained, tolerance) {
  data.frame(case = case, reference = reference, quantity = quantity, reference_value = value,
             latents = obtained, difference = abs(obtained - value), tolerance = tolerance,
             n_values = 1L, kind = "R package", agrees = abs(obtained - value) <= tolerance)
}

#' enumerate_lpa() against mclustBIC() on one dataset
#' @noRd
compare_mclust_grid <- function(name, data) {
  models <- c("EEI", "VVI", "EEE", "VVV")
  ours <- as.data.frame(enumerate_lpa(data, names(data), n_profiles = 1:4, model = models,
                                      n_starts = 50L, seed = 1L, tol = 1e-12, max_iter = 100000L))
  bic <- mclust::mclustBIC(data, G = 1:4, modelNames = models, control = .mclust_control(),
                           verbose = FALSE)
  n <- nrow(data)
  cells <- expand.grid(k = 1:4, model = models, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(cells)), function(i) {
    k <- cells$k[[i]]
    model <- cells$model[[i]]
    theirs_bic <- -bic[as.character(k), model]
    npar <- mclust::nMclustParams(model, ncol(data), k)
    theirs_ll <- (-theirs_bic + npar * log(n)) / 2
    mine <- ours[ours$n_profiles == k & ours$model == model, ]
    case <- sprintf("%s_%s_k%d", name, model, k)
    count <- .cell_row(case, "mclust", "number of parameters", npar, mine$n_parameters, 0)
    same <- is.finite(theirs_ll) && abs(mine$log_likelihood - theirs_ll) <= 1e-6 + 1e-10 * abs(theirs_ll)
    list(rows = rbind(count, if (same) rbind(
      .cell_row(case, "mclust", "log likelihood", theirs_ll, mine$log_likelihood,
                1e-6 + 1e-10 * abs(theirs_ll)),
      .cell_row(case, "mclust", "BIC", theirs_bic, mine$bic_individual,
                2 * (1e-6 + 1e-10 * abs(theirs_ll))))),
      search = data.frame(case = case, reference = "mclust", latents = mine$log_likelihood,
                          reference_ll = theirs_ll, same_maximum = same))
  })
  list(rows = do.call(rbind, lapply(rows, `[[`, "rows")),
       search = do.call(rbind, lapply(rows, `[[`, "search")))
}

#' enumerate_lca() against poLCA on one dataset
#' @noRd
compare_polca_grid <- function(name, data) {
  items <- names(data)
  ours <- as.data.frame(enumerate_lca(data, items, n_classes = 1:4, n_starts = 50L, seed = 1L,
                                      tol = 1e-12, max_iter = 100000L))
  formula <- stats::as.formula(sprintf("cbind(%s) ~ 1", paste(items, collapse = ", ")))
  coded <- as.data.frame(lapply(data, function(v) as.integer(v - min(v) + 1L)))
  rows <- lapply(1:4, function(k) {
    printed <- utils::capture.output(fit <- withr::with_seed(1L, poLCA::poLCA(
      formula, coded, nclass = k, nrep = if (k == 1L) 1L else 30L, maxiter = 100000L,
      tol = 1e-12, verbose = FALSE, calc.se = FALSE)))
    mine <- ours[ours$n_profiles == k, ]
    case <- sprintf("%s_k%d", name, k)
    count <- .cell_row(case, "poLCA", "number of parameters", fit$npar, mine$n_parameters, 0)
    same <- abs(mine$log_likelihood - fit$llik) <= 1e-6 + 1e-10 * abs(fit$llik)
    list(rows = rbind(count, if (same) rbind(
      .cell_row(case, "poLCA", "log likelihood", fit$llik, mine$log_likelihood,
                1e-6 + 1e-10 * abs(fit$llik)),
      .cell_row(case, "poLCA", "BIC", fit$bic, mine$bic_individual,
                2 * (1e-6 + 1e-10 * abs(fit$llik))))),
      search = data.frame(case = case, reference = "poLCA", latents = mine$log_likelihood,
                          reference_ll = fit$llik, same_maximum = same))
  })
  list(rows = do.call(rbind, lapply(rows, `[[`, "rows")),
       search = do.call(rbind, lapply(rows, `[[`, "search")))
}

#' mclust's VVV solution with k profiles as an lpa() start
#' @noRd
.mclust_start <- function(model) {
  sigma <- model$parameters$variance$sigma
  list(means = t(unname(as.matrix(model$parameters$mean))),
       profile_probabilities = matrix(model$parameters$pro %||% 1, 1L),
       group_probabilities = 1, covariances = unname(sigma))
}

#' bootstrap_lrt() against mclustBootstrapLRT() and mclust's vignette
#' @noRd
compare_blrt <- function() {
  diabetes <- .diabetes_vignette()
  replicates <- 199L
  reference <- withr::with_seed(1L, mclust::mclustBootstrapLRT(
    diabetes, modelName = "VVV", nboot = replicates, maxG = 3L, control = .mclust_control(),
    verbose = FALSE))
  # The vignette's statistics are mclust's at its default tolerance.
  default_ll <- mclust::mclustBIC(diabetes, G = 1:4, modelNames = "VVV",
                                  verbose = FALSE)[, "VVV"]
  # mclust's BIC is 2 LL - npar log n.
  default_ll <- (default_ll + vapply(1:4, function(k) mclust::nMclustParams("VVV", ncol(diabetes), k),
                                      numeric(1)) * log(nrow(diabetes))) / 2
  settings <- list(data = diabetes, vars = names(diabetes), model = "VVV", seed = 1L,
                   tol = 1e-12, max_iter = 100000L)
  # mclust's maxima (its hierarchical initialisation, as mclustBootstrapLRT
  # fits them), and lpa() started at each; one profile has one maximum.
  # (Mclust() and summary() re-evaluate their calls unqualified, which needs
  # mclust attached; summaryMclustBIC() of the BIC table is Mclust()'s model.)
  theirs <- lapply(1:4, function(k) mclust::summaryMclustBIC(
    mclust::mclustBIC(diabetes, G = k, modelNames = "VVV", control = .mclust_control(),
                      verbose = FALSE), diabetes))
  at_mclust <- lapply(1:4, function(k) {
    arguments <- c(settings, list(n_profiles = k),
                   if (k == 1L) list(n_starts = 50L) else
                     list(n_starts = 1L, start = .mclust_start(theirs[[k]])))
    do.call(lpa, arguments)
  })
  own <- lapply(1:4, function(k) do.call(lpa, c(settings, list(n_profiles = k, n_starts = 50L))))
  vignette <- c("369.98537", "121.10652", "19.90071")
  rows <- lapply(1:3, function(k) {
    # Replicates are refitted with lpa()'s own 50-start search.
    test <- bootstrap_lrt(at_mclust[[k]], at_mclust[[k + 1L]], data = diabetes,
                          iter = replicates, seed = 1L, n_starts = 50L)
    result <- as.data.frame(test)
    statistic <- result$statistic[[1L]]
    p <- result$p_value[[1L]]
    theirs_p <- reference$p.value[[k]]
    pooled <- (p + theirs_p) / 2
    p_tolerance <- 4 * sqrt(max(pooled * (1 - pooled), 1 / replicates) * 2 / replicates) +
      1 / replicates
    case <- sprintf("diabetes_VVV_k%d_vs_k%d", k, k + 1L)
    digits <- nchar(sub("^-?[0-9]*[.]?", "", vignette[[k]]))
    boot <- get_results(test, "replicates")
    list(rows = rbind(
      .cell_row(case, "mclust", "bootstrap LRT statistic", reference$obs[[k]], statistic,
                2 * (1e-6 + 1e-10 * abs(theirs[[k + 1L]]$loglik))),
      .cell_row(case, "mclust", sprintf("bootstrap p-value (%d replicates each)", replicates),
                theirs_p, p, p_tolerance),
      .cell_row(case, "mclust vignette (Scrucca et al.), diabetes", "bootstrap LRT statistic",
                as.numeric(vignette[[k]]), statistic, 0.5 * 10^(-digits))),
      bootstrap = data.frame(case = case, reference = "mclust",
                             latents_statistic = statistic, reference_statistic = reference$obs[[k]],
                             reference_default_tolerance = 2 * (default_ll[[k + 1L]] - default_ll[[k]]),
                             vignette = as.numeric(vignette[[k]]),
                             latents_mean = mean(boot$statistic[boot$valid]),
                             reference_mean = mean(reference$boot[, k]),
                             latents_p = p, reference_p = theirs_p))
  })
  search <- do.call(rbind, lapply(1:4, function(k) {
    data.frame(case = sprintf("diabetes_VVV_k%d", k), reference = "mclust",
               latents = own[[k]]$log_likelihood, reference_ll = theirs[[k]]$loglik,
               same_maximum = abs(own[[k]]$log_likelihood - theirs[[k]]$loglik) <=
                 1e-6 + 1e-10 * abs(theirs[[k]]$loglik))
  }))
  at <- do.call(rbind, lapply(2:4, function(k) {
    .cell_row(sprintf("diabetes_VVV_k%d", k), "mclust", "log likelihood at mclust's solution",
              theirs[[k]]$loglik, at_mclust[[k]]$log_likelihood,
              1e-6 + 1e-10 * abs(theirs[[k]]$loglik))
  }))
  out <- rbind(at, do.call(rbind, lapply(rows, `[[`, "rows")))
  attr(out, "search") <- search
  attr(out, "bootstrap") <- do.call(rbind, lapply(rows, `[[`, "bootstrap"))
  out
}

compare_all_enumeration <- function() {
  data_of <- function(name, package) {
    get(utils::data(list = name, package = package, envir = environment()), envir = environment())
  }
  zero_one <- function(frame) as.data.frame(lapply(frame, function(v) { v <- as.integer(v); v - min(v) }))
  lpa_sets <- list(faithful = datasets::faithful,
                   iris = datasets::iris[c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")],
                   diabetes = .diabetes_vignette())
  lca_sets <- list(carcinoma = zero_one(data_of("carcinoma", "poLCA")),
                   values = zero_one(data_of("values", "poLCA")),
                   cheating = zero_one(data_of("cheating", "poLCA")[c("LIEEXAM", "LIEPAPER", "FRAUD", "COPYEXAM")]))
  grids <- c(Map(compare_mclust_grid, names(lpa_sets), lpa_sets),
             Map(compare_polca_grid, names(lca_sets), lca_sets))
  enumeration <- do.call(rbind, lapply(grids, `[[`, "rows"))
  enumeration$fn <- "enumeration"
  blrt <- compare_blrt()
  blrt_search <- attr(blrt, "search")
  bootstrap <- attr(blrt, "bootstrap")
  blrt$fn <- "bootstrap_lrt"
  out <- rbind(enumeration, blrt)
  row.names(out) <- NULL
  # The bootstrap's own-search fits are the grid's diabetes VVV cells again.
  search <- unique(rbind(do.call(rbind, lapply(grids, `[[`, "search")), blrt_search))
  row.names(search) <- NULL
  attr(out, "search") <- search
  attr(out, "tables") <- list(bootstrap = bootstrap)
  out
}
