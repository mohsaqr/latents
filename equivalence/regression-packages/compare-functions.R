# Functions for compare.R: fit mixtures of regressions with independent R
# implementations (mixtools regmixEM and MoEClust for Gaussian outcomes,
# flexmix for Poisson and binomial outcomes) and with mixture_regression(),
# and compare the maximized likelihood, class shares, coefficients, residual
# standard deviations and posteriors. Published results (SAS PROC FMM Example
# 40.3; Grun and Leisch 2008; Benaglia et al. 2009) are compared to their
# printed precision.
#
# Tolerances, fixed before any comparison was run: software references, the
# log likelihood to 1e-6 plus 1e-10 of its size and the rest to 1e-4 of
# max(1, |value|); published values, half a unit of the last printed digit.
# Every package runs on one core.

.case <- function(name, family, formula, data, classes, references, common = NULL,
                  membership = ~1, id = NULL) {
  list(name = name, family = family, formula = formula, data = data, classes = classes,
       references = references, common = common, membership = membership, id = id)
}

#' The Ames salmonella assay of Margolin, Kaplan and Zeiger (1981), as printed
#' in the SAS/STAT 14.2 PROC FMM documentation (Example 40.3)
#' @noRd
salmonella <- function() {
  dose <- c(0, 10, 33, 100, 333, 1000)
  counts <- rbind(c(15, 21, 29), c(16, 18, 21), c(16, 26, 33), c(27, 41, 60),
                  c(33, 38, 41), c(20, 27, 42))
  data.frame(dose = rep(dose, each = 3L), num = as.vector(t(counts)),
             logd = log(rep(dose, each = 3L) + 10))
}

regression_cases <- function() {
  data_of <- function(name, package) {
    get(utils::data(list = name, package = package, envir = environment()), envir = environment())
  }
  co2 <- data_of("CO2data", "mixtools")
  tone <- data_of("tonedata", "mixtools")
  no <- data_of("NOdata", "mixtools")
  npreg <- data_of("NPreg", "flexmix")
  betablocker <- transform(data_of("betablocker", "flexmix"), Failures = Total - Deaths)
  simulated <- withr::with_seed(20261301L, {
    z <- stats::rnorm(400L)
    x <- stats::rnorm(400L)
    class <- 1L + (stats::runif(400L) < stats::plogis(-0.2 + 1.2 * z))
    data.frame(x = x, z = z,
               y = c(1, -1)[class] + c(0.8, -0.5)[class] * x + stats::rnorm(400L, sd = c(0.5, 0.8)[class]))
  })
  gaussian <- c("mixtools", "MoEClust")
  list(
    .case("CO2data_k2", "gaussian", CO2 ~ GNP, co2, 2L, gaussian),
    .case("tonedata_k2", "gaussian", tuned ~ stretchratio, tone, 2L, gaussian),
    .case("NOdata_k2", "gaussian", Equivalence ~ NO, no, 2L, gaussian),
    .case("NPreg_gaussian_k2", "gaussian", yn ~ x + I(x^2), npreg, 2L, gaussian),
    .case("simulated_membership_k2", "gaussian", y ~ x, simulated, 2L, "MoEClust", membership = ~z),
    .case("NPreg_poisson_k2", "poisson", yp ~ x, npreg, 2L, "flexmix"),
    .case("salmonella_k2", "poisson", num ~ dose + logd, salmonella(), 2L, "flexmix",
          common = ~ dose + logd),
    .case("betablocker_k3", "binomial", cbind(Deaths, Failures) ~ Treatment, betablocker, 3L,
          "flexmix", common = ~Treatment, id = "Center"))
}

#' mixture_regression() fit and its solution in one shape
fit_latents_regression <- function(case) {
  fit <- mixture_regression(case$formula, case$data, case$classes, family = case$family,
                            common = case$common, membership = case$membership, id = case$id,
                            class_level = if (is.null(case$id)) "observation" else "group",
                            n_starts = 50L, seed = 1L, tol = 1e-14, max_iter = 200000L)
  theta <- stats::coef(fit)
  k <- case$classes
  terms <- colnames(fit$spec$x)
  common <- colnames(fit$spec$z) %||% character()
  coefficients <- t(vapply(seq_len(k), function(j) c(
    vapply(terms, function(t) theta[[sprintf("coefficient.class_%d.%s", j, t)]], numeric(1)),
    vapply(common, function(t) theta[[sprintf("coefficient.common.%s", t)]], numeric(1))),
    numeric(length(terms) + length(common))))
  sigma <- if (identical(case$family, "gaussian")) {
    vapply(seq_len(k), function(j) exp(theta[[sprintf("log_sigma.class_%d", j)]]), numeric(1))
  }
  tau <- fit$expectation$tau
  list(fit = fit, log_likelihood = as.numeric(stats::logLik(fit)),
       shares = colMeans(tau), coefficients = unname(coefficients),
       sigma = sigma, posteriors = unname(tau))
}

fit_mixtools <- function(case) {
  frame <- stats::model.frame(case$formula, case$data)
  y <- stats::model.response(frame)
  x <- stats::model.matrix(case$formula, frame)[, -1L, drop = FALSE]
  fits <- lapply(seq_len(20L), function(s) withr::with_seed(s, tryCatch({
    printed <- utils::capture.output(fitted <- mixtools::regmixEM(
      y, x, k = case$classes, epsilon = 1e-14, maxit = 200000L))
    fitted
  }, error = function(e) {
    message(sprintf("%s: regmixEM start %d failed: %s", case$name, s, conditionMessage(e)))
    NULL
  })))
  fits <- Filter(Negate(is.null), fits)
  best <- fits[[which.max(vapply(fits, function(f) f$loglik, numeric(1)))]]
  list(log_likelihood = best$loglik, shares = best$lambda, coefficients = unname(t(best$beta)),
       sigma = best$sigma, posteriors = unname(best$posterior))
}

fit_moeclust <- function(case) {
  frame <- stats::model.frame(case$formula, case$data)
  response <- stats::model.response(frame)
  network <- cbind(case$data, .y = response)
  expert <- stats::reformulate(attr(stats::terms(case$formula), "term.labels"))
  # MoEClust starts from one hierarchical clustering by default; 20 random
  # starts (seeded) give it the same chance as the other packages.
  printed <- utils::capture.output(fit <- withr::with_seed(1L, MoEClust::MoE_clust(
    response, G = case$classes, expert = expert, gating = case$membership, modelNames = "V",
    network.data = network, control = MoEClust::MoE_control(tol = c(1e-14, 1e-14, 1e-14),
                                                              itmax = c(200000L, 200000L, 200000L),
                                                              init.z = "random", nstarts = 20L,
                                                              verbose = FALSE))))
  coefficients <- t(vapply(fit$expert, stats::coef, numeric(ncol(stats::model.matrix(case$formula, frame)))))
  # MoEClust keeps the likelihood of every iteration; the last is the maximum reached.
  list(log_likelihood = utils::tail(as.numeric(fit$loglik), 1L), shares = colMeans(fit$z),
       coefficients = unname(coefficients),
       sigma = sqrt(as.numeric(fit$parameters$variance$sigmasq)), posteriors = unname(fit$z))
}

fit_flexmix <- function(case) {
  common <- case$common
  formula <- if (is.null(common)) case$formula else stats::update(case$formula, . ~ 1)
  if (!is.null(case$id)) {
    formula <- stats::as.formula(paste(deparse(formula), "|", case$id))
  }
  model <- if (is.null(common)) flexmix::FLXMRglm(family = case$family) else
    flexmix::FLXMRglmfix(family = case$family, fixed = common)
  fits <- lapply(seq_len(20L), function(s) withr::with_seed(s, tryCatch(
    flexmix::flexmix(formula, data = case$data, k = case$classes, model = model,
                     control = list(tolerance = 1e-14, iter.max = 200000L, minprior = 0)),
    error = function(e) {
      message(sprintf("%s: flexmix start %d failed: %s", case$name, s, conditionMessage(e)))
      NULL
    })))
  fits <- Filter(function(f) !is.null(f) && f@k == case$classes, fits)
  best <- fits[[which.max(vapply(fits, function(f) f@logLik, numeric(1)))]]
  parameters <- flexmix::parameters(best)
  rows <- rownames(parameters)
  intercepts <- parameters[grep("Intercept", rows), , drop = FALSE]
  slopes <- parameters[setdiff(rows, grep("Intercept", rows, value = TRUE)), , drop = FALSE]
  coefficients <- t(rbind(intercepts, slopes))
  posteriors <- flexmix::posterior(best)
  if (!is.null(case$id)) posteriors <- posteriors[!duplicated(case$data[[case$id]]), , drop = FALSE]
  list(log_likelihood = best@logLik, shares = colMeans(posteriors),
       coefficients = unname(coefficients), sigma = NULL, posteriors = unname(posteriors))
}

.fitters <- list(mixtools = fit_mixtools, MoEClust = fit_moeclust, flexmix = fit_flexmix)

.permutations <- function(k) {
  if (k == 1L) return(list(1L))
  unlist(lapply(seq_len(k), function(first) lapply(.permutations(k - 1L), function(rest)
    c(first, setdiff(seq_len(k), first)[rest]))), recursive = FALSE)
}

.relative_gap <- function(reference, ours) max(abs(reference - ours) / pmax(1, abs(reference)))

compare_solution <- function(case, package, ours, reference) {
  perms <- .permutations(case$classes)
  gap <- vapply(perms, function(perm)
    .relative_gap(reference$coefficients[perm, , drop = FALSE], ours$coefficients), numeric(1))
  perm <- perms[[which.min(gap)]]
  k <- case$classes
  rows <- data.frame(
    case = case$name, reference = package,
    quantity = c("maximized log likelihood", sprintf("class shares (%d values)", k),
                 sprintf("coefficients (%d values)", length(ours$coefficients))),
    reference_value = c(reference$log_likelihood, NA, NA), latents = c(ours$log_likelihood, NA, NA),
    difference = c(abs(reference$log_likelihood - ours$log_likelihood),
                   .relative_gap(reference$shares[perm], ours$shares), min(gap)),
    tolerance = c(1e-6 + 1e-10 * abs(ours$log_likelihood), 1e-4, 1e-4),
    n_values = c(1L, k, length(ours$coefficients)))
  if (!is.null(reference$sigma)) {
    rows <- rbind(rows, data.frame(case = case$name, reference = package,
                                   quantity = sprintf("residual standard deviations (%d values)", k),
                                   reference_value = NA, latents = NA,
                                   difference = .relative_gap(reference$sigma[perm], ours$sigma),
                                   tolerance = 1e-4, n_values = k))
  }
  if (!is.null(reference$posteriors) && nrow(reference$posteriors) == nrow(ours$posteriors)) {
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

#' A printed value compared to half a unit of its last printed digit
#' @noRd
.printed <- function(case, source, quantity, printed, obtained) {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  tolerance <- 0.5 * 10^(-digits) + 1e-12
  value <- as.numeric(printed)
  data.frame(case = case, reference = source, quantity = quantity, reference_value = value,
             latents = obtained, difference = abs(obtained - value), tolerance = tolerance,
             n_values = 1L, agrees = abs(obtained - value) <= tolerance)
}

#' The published values these cases reproduce
#' @noRd
compare_published <- function(results) {
  sal <- results[["salmonella_k2"]]
  theta <- stats::coef(sal$fit)
  errors <- sqrt(diag(stats::vcov(sal$fit)))
  big <- if (theta[["coefficient.class_1.(Intercept)"]] < theta[["coefficient.class_2.(Intercept)"]]) 1L else 2L
  small <- 3L - big
  logit <- if (big == 1L) -theta[["membership.class_2.(Intercept)"]] else theta[["membership.class_2.(Intercept)"]]
  # SAS's component 1 is the larger component, with the smaller intercept.
  sas <- "SAS/STAT 14.2 PROC FMM Example 40.3 (Wang et al. 1996 model)"
  one <- mixture_regression(num ~ dose + logd, salmonella(), 1L, family = "poisson",
                            n_starts = 1L, seed = 1L, tol = 1e-14)
  rows <- rbind(
    .printed("salmonella_k2", sas, "-2 log likelihood", "121.8", -2 * sal$log_likelihood),
    .printed("salmonella_k2", sas, "component 1 intercept", "1.9097",
             theta[[sprintf("coefficient.class_%d.(Intercept)", big)]]),
    .printed("salmonella_k2", sas, "component 2 intercept", "2.4770",
             theta[[sprintf("coefficient.class_%d.(Intercept)", small)]]),
    .printed("salmonella_k2", sas, "common dose", "-0.00126", theta[["coefficient.common.dose"]]),
    .printed("salmonella_k2", sas, "common logd", "0.3639", theta[["coefficient.common.logd"]]),
    .printed("salmonella_k2", sas, "component 1 intercept s.e.", "0.2654",
             errors[[sprintf("coefficient.class_%d.(Intercept)", big)]]),
    .printed("salmonella_k2", sas, "component 2 intercept s.e.", "0.2731",
             errors[[sprintf("coefficient.class_%d.(Intercept)", small)]]),
    .printed("salmonella_k2", sas, "common dose s.e.", "0.000273", errors[["coefficient.common.dose"]]),
    .printed("salmonella_k2", sas, "common logd s.e.", "0.06602", errors[["coefficient.common.logd"]]),
    .printed("salmonella_k2", sas, "component 1 mixing probability", "0.8173",
             stats::plogis(logit)),
    .printed("salmonella_k2", sas, "component 1 logit", "1.4984", logit),
    .printed("salmonella_k2", sas, "component 1 logit s.e.", "0.6875",
             errors[["membership.class_2.(Intercept)"]]),
    .printed("salmonella_k1", sas, "one class: -2 log likelihood", "136.3",
             -2 * as.numeric(stats::logLik(one))),
    .printed("salmonella_k1", sas, "one class: logd", "0.3198",
             stats::coef(one)[["coefficient.class_1.logd"]]))
  beta <- results[["betablocker_k3"]]
  co2 <- results[["CO2data_k2"]]
  gl <- "Grun and Leisch (2008) JSS 28(4), betablocker"
  bg <- "Benaglia et al. (2009) JSS 32(6), CO2data"
  rbind(rows,
        .printed("betablocker_k3", gl, "log likelihood", "-159.3605", beta$log_likelihood),
        .printed("betablocker_k3", gl, "common Treatment", "-0.2582",
                 stats::coef(beta$fit)[["coefficient.common.TreatmentTreated"]]),
        .printed("CO2data_k2", bg, "log likelihood", "-66.93977", co2$log_likelihood),
        .printed("CO2data_k2", bg, "residual standard deviations (larger)", "2.0493", max(co2$sigma)),
        .printed("CO2data_k2", bg, "residual standard deviations (smaller)", "0.8094", min(co2$sigma)),
        .printed("CO2data_k2", bg, "larger share", "0.7549", max(co2$shares)))
}

#' Fit every case, compare with its references and with published values
#' @return One row per compared quantity, with `fn = "mixture_regression"`.
compare_all_regression <- function(cases = regression_cases()) {
  results <- lapply(stats::setNames(cases, vapply(cases, `[[`, "", "name")), fit_latents_regression)
  software <- do.call(rbind, lapply(cases, function(case) {
    ours <- results[[case$name]]
    do.call(rbind, lapply(case$references, function(package) {
      compare_solution(case, package, ours, .fitters[[package]](case))
    }))
  }))
  software$kind <- "R package"
  published <- compare_published(results)
  published$kind <- "published"
  out <- rbind(software, published)
  out$fn <- "mixture_regression"
  out
}
